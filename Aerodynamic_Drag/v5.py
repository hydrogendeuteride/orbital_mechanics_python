import numpy as np
from math import sqrt
from sys import path
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

tspan = 3600 * 24 * 2.0
dt = 100.0

cb = pd.earth

if __name__ == '__main__':

    perts = null_perts()
    perts['aero'] = True
    perts['Cd'] = 2.2
    perts['A'] = (1e-3)**2 /4.0

    mass0 = 10.0

    op0 = OP(t.tle2coes('C:\\Users\\soft_kitty\\VScode Source\\python\\kepler_orbit\\ISS.txt'), tspan, dt, coes = True, deg = False, perts = perts, mass0 = mass0)
    op0.plot_alts(show_plot = True, hours = True)
    t.plot_n_orbits([op0.rs], labels = ['ISS'], show_plot = True)
    op0.calculate_coes(degrees = True)
    op0.plot_coes(show_plot = True, hours = True)
    op0.calculate_apoapse_periapse()
    op0.plot_apoapse_periapse(show_plot = True, hours = True)