import numpy as np
import matplotlib.pyplot as plt
import planetary_data as pd
import mpl_toolkits.mplot3d as Axes3D

d2r = np.pi/180.0

def plot_n_orbits(rs, labels, cb = pd.earth, show_plot = False, save_plot = False, title = 'Test'):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111, projection = '3d')

    r_plot = cb['radius']

    n = 0
    for r in rs:    
        ax.plot(r[:,0],r[:,1],r[:,2],label = labels[n])
        ax.plot([r[0,0]],[r[0,1]],[r[0,2]])
        n+=1

    _u, _v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x = r_plot * np.cos(_u)*np.sin(_v)
    _y = r_plot * np.sin(_u)*np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap = 'Blues', alpha = 0.6)

    l = r_plot * 2.0
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
    u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w, color = 'black')
    
    max_val = np.max(np.abs(rs))

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')

    ax.set_title(title)

    plt.legend()
        
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+'.png', dpi = 300)