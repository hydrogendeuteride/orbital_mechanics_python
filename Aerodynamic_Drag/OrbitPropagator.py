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
        'solar_grav' : False,
        'aero' : False,
        'Cd': 0,
        'A': 0,
        'mu':0,
        'third_bodies':[],
        'srp':False,
        'srp_custom_func' : False,
        'CR':0,
        'B':0,
        'oblateness':False,
        'J3':False,
        'J4':False,
        'J5':False,
        'J6':False,
        'J7':False,
        'relativity': False,
        'thrust':0,
        'thrust_direction':0,
        'rho':0,
        'C20':False,
        'custom_pert':False
    }

class OrbitPropagator:
    def __init__(self,state0,tspan,dt,coes = False,mass0 = 0,deg = True,cb = pd.earth, perts = null_perts()):
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

        self.mass = mass0

        self.propagate_orbit()

    def diff_q(self,t_,y):
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

        if self.perts['aero']:
            z = norm_r - self.cb['radius']
            rho = t.calc_atmospheric_density(z)
            v_rel = v - np.cross(self.cb['atm_rot_vector'], r)
            a -= v_rel * 0.5 * rho *np.linalg.norm(v_rel) * self.perts['Cd'] * self.perts['A'] / self.mass

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

    def calculate_apoapse_periapse(self):
        self.apoapses = self.coes[:,0] * (1 + self.coes[:,1])
        self.periapses = self.coes[:,0] * (1 - self.coes[:,1])

    def plot_apoapse_periapse(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'Apoaose and Periapse -time', dpi = 500):
        plt.figure(figsize=(20,10))
        
        if hours:
            ts = self.ts / 3600.0
            x_unit = 'Hours'
        elif days:
            ts = self.ts / (3600.0*24.0)
            x_unit = 'Days'
        else:
            ts = self.ts
            x_unit = "seconds"

        plt.plot(ts, self.apoapses, 'b', label = 'Apoapse')
        plt.plot(ts, self.periapses, 'r', label = 'Periapse')
        plt.grid(True)
        plt.xlabel('Time(%s)' % x_unit)
        plt.ylabel('Altitude')
        plt.title(title)

        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png', dpi = dpi)

    def plot_alts(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'Distance -time', dpi = 500, figsize = (16,8)):
        if hours:
            ts = self.ts / 3600.0
            x_unit = 'Hours'
        elif days:
            ts = self.ts / (3600.0*24.0)
            x_unit = 'Days'
        else:
            ts = self.ts
            x_unit = "seconds"

        plt.figure(figsize = figsize)
        plt.plot(ts, self.alts, 'w')
        plt.grid(True)
        plt.xlabel('Time(%s)' % x_unit)
        plt.ylabel('Altitude')
        plt.title(title)

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png', dpi = dpi)

    def propagate_orbit(self):
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] =self.solver.t
            self.ys[self.step] = self.solver.y
            self.step+=1

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]
        self.ts = self.ts[:self.step]
        self.alts = np.reshape((np.linalg.norm(self.rs, axis = 1) - self.cb['radius']), (self.step,1))