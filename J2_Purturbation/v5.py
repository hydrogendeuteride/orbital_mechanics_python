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
    perts['J2'] = True

    op0 = OP(t.tle2coes('C:\\Users\\soft_kitty\\VScode Source\\python\\kepler_orbit\\ISS.txt'), tspan, dt, coes = True, deg = False, perts = perts)

    op0.calculate_coes()
    #t.plot_n_orbits([op0.rs], labels = ['ISS'], show_plot = True)
    op0.plot_coes(show_plot = True, hours = True)