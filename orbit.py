import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

earth_radius = 6378.0
earth_mu = 398600.0

def plot(r):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111, projection = '3d')

    r_plot = earth_radius
    
    ax.plot(r[:,0],r[:,1],r[:,2],'w')
    ax.plot([r[0,0]],[r[0,1]],[r[0,2]],'wo')

    _u, _v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x = r_plot * np.cos(_u)*np.sin(_v)
    _y = r_plot * np.sin(_u)*np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap = 'Blues', alpha = 0.6)

    l = r_plot * 2.0
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
    u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w, color = 'w')
   
    max_val = np.max(np.abs(r))

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')

    plt.legend(['Trajectory', 'starting position'])
    plt.show()

def diff_q(t,y,mu):
    rx,ry,rz,vx,vy,vz = y
    r = np.array([rx,ry,rz])

    norm_r = np.linalg.norm(r)

    ax,ay,az = -r*mu/norm_r**3

    return [vx,vy,vz,ax,ay,az]

if __name__ == '__main__':
    r_mag = earth_radius+500.0
    v_mag = np.sqrt(earth_mu/r_mag)

    r0 = [r_mag,0,0]
    v0 = [0, v_mag,0]
    
    tspan = 100*60.0

    dt = 100.0

    n_steps = int(np.ceil(tspan/dt))

    ys = np.zeros((n_steps, 6))
    ts = np.zeros((n_steps, 1))

    y0 = r0+v0
    ys[0] = np.array(y0)
    step = 1

    solver = ode(diff_q)
    solver.set_integrator('lsoda')
    solver.set_initial_value(y0,0)
    solver.set_f_params(earth_mu)

    while solver.successful() and step < n_steps:
        solver.integrate(solver.t+dt)
        ts[step] =solver.t
        ys[step] = solver.y
        step+=1

    rs = ys[:,:3]

    plot(rs)