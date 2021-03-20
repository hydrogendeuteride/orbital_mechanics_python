import numpy as np
import matplotlib.pyplot as plt
import datetime
import planetary_data as pd
import mpl_toolkits.mplot3d as Axes3D
import math as m

d2r = np.pi/180.0
r2d = 1.0 / d2r

def norm(v):
    return np.linalg.norm(v)

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

def coes2rv(coes, deg = False, mu = pd.earth['mu']):
    if deg:
        a,e,i,ta,aop,raan, date = coes
        i *= d2r
        ta *= d2r
        aop *= d2r
        raan *= d2r
    else:
        a,e,i,ta,aop,raan, date = coes

    E = ecc_anomaly([ta, e], 'tae')

    r_norm = a * (1 - e **2) / (1 + e * np.cos(ta))

    r_perif = r_norm * np.array([m.cos(ta), m.sin(ta), 0])
    v_perif = np.sqrt(mu * a) / r_norm * np.array([-m.sin(E), m.cos(E) * m.sqrt(1 - e**2), 0])

    perif2eci = np.transpose(eci2perif(raan, aop, i))

    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r,v, date

def rv2coes(r, v, mu = pd.earth['mu'], degrees = False, print_results = False):
    r_norm = norm(r)

    h = np.cross(r,v)
    h_norm = norm(h)

    i = m.acos(h[2] / h_norm)
    e = ((norm(v)**2 - mu / r_norm) * r - np.dot(r,v) * v) / mu
    e_norm = norm(e)
    
    N = np.cross([0,0,1], h)
    N_norm = norm(N)

    raan = m.acos(N[0] / N_norm)

    if N[1]<0:
        raan - 2 * np.pi - raan

    aop = m.acos(np.dot(N,e) / N_norm / e_norm)

    if e[2] < 0:
        aop = 2 * np.pi - aop

    ta = m.acos(np.dot(e,r) / e_norm / r_norm)

    if np.dot(r , v) < 0:
        ta = 2 * np.pi - ta

    a = r_norm * (1 + e_norm * m.cos(ta)) / (1 - e_norm ** 2)
    if print_results:
        print('a', a)
        print('e', e_norm)
        print('i',i * r2d)
        print('RAAN', raan * r2d)
        print('AOP', aop * r2d)
        print('TA',ta * r2d)

    if degrees:
        return [a,e_norm, i * r2d, ta * r2d, aop * r2d, raan * r2d]
    else:
        return [a,e_norm, i, ta, aop, raan]

def ecc_anomaly(arr, method, tol = 1e-8):
    if method == 'newton':
        Me, e = arr
        if Me < np.pi/2.0 :
            E0 = Me + e/2.0
        else:
            E0 = Me - e
        
        for n in range(200):
            ratio = (E0 - e*m.sin(E0)-Me) / (1-e * m.cos(E0))
            if abs(ratio)<tol:
                if n==0:
                    return E0
                else:
                    return E1
            else:
                E1 = E0 - ratio
                E0 = E1
        return False

    elif method == 'tae':
        ta, e =arr
        return 2 * m.atan(np.sqrt((1-e)/(1+e)) * m.tan(ta/2.0))
    else:
        print('INVALID METHOD FOR ECCENTRIC ANOMALY')
    
def eci2perif(raan, aop, i):
    row0 = [-m.sin(raan) * m.cos(i) * m.sin(aop) + m.cos(raan) * m.cos(aop), m.cos(raan) * m.cos(i) * m.sin(aop) + m.sin(raan) * m.cos(aop), m.sin(i) * m.sin(aop)]
    row1 = [-m.sin(raan) * m.cos(i) * m.cos(aop) - m.cos(raan) * m.sin(aop), m.cos(raan) * m.cos(i) * m.cos(aop) - m.sin(raan) * m.sin(aop), m.sin(i) * m.cos(aop)]
    row2 = [m.sin(raan) * m.sin(i), -m.cos(raan) * m.sin(i), m.cos(i)]

    return np.array([row0, row1, row2])

def tle2coes(tle_filename, mu = pd.earth['mu']):
    with open(tle_filename, 'r') as f:
        lines = f.readlines()
    
    line0 = lines[0].strip()
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    epoch = line1[3]
    year, month, day, hour = calc_epoch(epoch)

    i = float(line2[2]) * d2r
    raan = float(line2[3]) * d2r
    e_string = line2[4]
    e = float('0.'+ e_string)
    aop = float(line2[5]) * d2r
    Me = float(line2[6]) * d2r
    mean_motion = float(line2[7])

    T = 1 / mean_motion * 24 * 3600

    a = (T**2 * mu / 4.0 / np.pi**2) ** (1/3.0)

    E = ecc_anomaly([Me, e], 'newton')

    ta = true_anomaly([E, e])

    r_mag = a * (1 - e * np.cos(E))

    return a,e,i,ta,aop, raan, [year, month, day, hour]

def calc_epoch(epoch):
    year = int('20' + epoch[:2])
    epoch = epoch[2:].split('.')
    day_of_year = int(epoch[0])-1
    hour = float('0.'+epoch[1]) * 24.0
    date = datetime.date(year,1,1)+datetime.timedelta(day_of_year)

    month = float(date.month)
    day = float(date.day)

    return year, month, day, hour

def true_anomaly(arr):
    E, e = arr
    return 2 * np.arctan(np.sqrt((1+e) / (1 - e)) * np.tan(E / 2.0))

def tle2rv(tle_filename):
    return coes2rv(tle2coes(tle_filename))