import numpy as np
from math import sqrt
from sys import path
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import spice_tools as st
import spiceypy as spice

tspan = 3600 * 24 * 1.0
dt = 100.0

STEPS = 100000
FRAME = 'ECLIPJ2000'
OBSERVER = 'SUN'


cb = pd.earth

if __name__ == '__main__':

    spice.furnsh('Ephemeris/spicefile/solar_system.mk')
    ids, names, tcs_sec, tcs_cal = st.get_objects('Ephemeris/spicefile/de440s.bsp', display = True)

    names = [f for f in names if 'BARYCENTER' in f]

    times = st.tc2array(tcs_sec[0], STEPS)

    rs = []

    for name in names:
        rs.append(st.get_ephemeris_data(name, times, FRAME, OBSERVER))

    t.plot_n_orbits(rs, names, show_plot = True, AU = True, cb = pd.sun, figsize = (20, 10))