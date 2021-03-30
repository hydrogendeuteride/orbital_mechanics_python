import numpy as np

G_meter = 6.67408e-11
G_Kmeter = G_meter * 10**-9

sun = {
    'name' : 'Sun',
    'mass' : 1.989e30, 
    'mu' : 1.32712e11,
    'radius' : 695700.0
}

atm = np.array([[63.096, 2.059e-4], [251.189, 5.909e-11], [1000.0, 3.561e-15]])
earth = {
    'name' : 'Earth',
    'mass' : 5.972e24, 
    'mu' : 5.972e24 * G_Kmeter,
    'radius' : 6378.0,
    'J2' : -1.082635854e-3,
    'zs': atm[:,0],
    'rhos': atm[:,1] * 10**8,
    'atm_rot_vector': np.array([0.0,0.0,72.9211e-6])
}
