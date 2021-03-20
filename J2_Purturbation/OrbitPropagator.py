import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import tools as t
import planetary_data as pd

plt.style.use('dark_background')

def null_perts():
    return {
        'J2' : False,
        'aero' : False,
        'moon_grav' : False,
        'solar_grav' : False
    }

class OrbitPropagator:
    def __init__(self,state0,tspan,dt,coes = False,deg = True,cb = pd.earth, perts = null_perts()):
        if coes:
            self.r0, self.v0, _ = t.coes2rv(state0, deg = deg, mu = cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]

        self.tspan = tspan
        self.dt = dt
        self.cb = cb

        self.n_steps = int(np.ceil(self.tspan/self.dt))

        self.ys = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))

        self.y0 = self.r0.tolist() +self.v0.tolist()
        self.ts[0] = 0
        self.ys[0] = self.y0
        self.step = 1

        self.solver = ode(self.diff_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0, 0)

        self.perts = perts

        self.propagate_orbit()

    def diff_q(self,t,y):
        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        norm_r = np.linalg.norm(r)

        a = -r*self.cb['mu']/norm_r**3

        if self.perts['J2']:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 -1)
            ty = r[1] / norm_r * (5 * z2 / r2 -1)
            tz = r[2] / norm_r * (5 * z2 / r2 -3)

            a_j2 = 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])

            a += a_j2

        return [vx,vy,vz,a[0],a[1],a[2]]

    def calculate_coes(self, degrees = False):
        print("Calculating COES")

        self.coes = np.zeros((self.n_steps, 6))

        for n in range(self.step):
            self.coes[n,:] = t.rv2coes(self.rs[n,:], self.vs[n,:], mu = self.cb['mu'],degrees = degrees)

    def plot_coes(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'COEs', figsize = (16,8)):
        print('Plotting...')
        fig,axs = plt.subplots(nrows=2,ncols=3,figsize = figsize)
        fig.suptitle(title, fontsize = 20)

        if hours:
            ts = self.ts / 3600.0
            xlabel = "time elapsed(h)"
        elif self.days:
            ts = self.ts / 3600.0 / 24.0
            xlabel = "time elapsed(d)"
        else:
            ts = self.ts
            xlabel = "time elapsed(s)"

        axs[0,0].plot(ts, self.coes[:,3])
        axs[0,0].set_title ('true anomaly - time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('angle(deg)')

        axs[1,0].plot(ts, self.coes[:,0])
        axs[1,0].set_title ('semimajor axis - time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('semimajor axis(deg)')
        axs[1,0].set_xlabel(xlabel)

        axs[0,1].plot(ts, self.coes[:,1])
        axs[0,1].set_title ('eccentricity - time')
        axs[0,1].grid(True)
        
        axs[0,2].plot(ts, self.coes[:,4])
        axs[0,2].set_title ('argument of periapse - time')
        axs[0,2].grid(True)
        
        axs[1,1].plot(ts, self.coes[:,2])
        axs[1,1].set_title ('inclination - time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('angle(deg)')
        axs[1,1].set_xlabel(xlabel)

        axs[1,2].plot(ts, self.coes[:,5])
        axs[1,2].set_title ('RAAN - time')
        axs[1,2].grid(True)
        axs[1,2].set_xlabel(xlabel)

        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title+'.png', dpi = 300)

    def propagate_orbit(self):
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] =self.solver.t
            self.ys[self.step] = self.solver.y
            self.step+=1

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]