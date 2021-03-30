import numpy as np
from math import sqrt
from sys import path
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
import spice_tools as st
import spiceypy as spice
import lambert_tools as lt

dt = 1000.0

STEPS = 100000
FRAME = 'ECLIPJ2000'
OBSERVER = 'SUN'

date0 = '2005-12-01'
datef = '2006-03-01'

HEADER = 't,rx,ry,rz'

cb = pd.sun

if __name__ == '__main__':

    spice.furnsh('Ephemeris/spicefile/solar_system.mk')

    et0 = spice.utc2et(date0)
    etf = spice.utc2et(datef)

    transfer_time = etf - et0

    time_arr = np.linspace(et0, etf, 10000)
    states_earth = st.get_ephemeris_data('EARTH', time_arr, FRAME ,OBSERVER)
    states_venus = st.get_ephemeris_data('VENUS', time_arr, FRAME ,OBSERVER)
    
    r0 = states_earth[0,:3]

    statef_venus = st.get_ephemeris_data('VENUS', [etf], FRAME, OBSERVER)

    rf = states_venus[-1,:3]

    v0, vf = lt.lamberts_universal_variables(r0, rf, transfer_time, mu = cb['mu'])

    state0_sc = r0.tolist() + v0.tolist()

    op_sc = OP(r0, v0, transfer_time, dt, cb = cb)

    t.plot_n_orbits([states_earth[:,:3], states_venus[:,:3], op_sc.rs], labels = ['Earth', 'Venus', 'Spacecraft'], cb = cb, show_plot = True, title = 'Earth to Venus Transfer')