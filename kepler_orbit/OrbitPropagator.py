import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import tools as t
import planetary_data as pd

class OrbitPropagator:
    def __init__(self,state0,tspan,dt,coes = False,deg = True,cb = pd.earth):
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

        self.propagate_orbit()

    def diff_q(self,t,y):
        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        norm_r = np.linalg.norm(r)

        ax,ay,az = -r*self.cb['mu']/norm_r**3

        return [vx,vy,vz,ax,ay,az]

    def propagate_orbit(self):
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] =self.solver.t
            self.ys[self.step] = self.solver.y
            self.step+=1

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]