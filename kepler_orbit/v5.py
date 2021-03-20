import numpy as np
from math import sqrt
from sys import path
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

tspan = 3600 * 24 * 1.0
dt = 10.0

cb = pd.earth

if __name__ == '__main__':
    op0 = OP(t.tle2coes('C:\\Users\\soft_kitty\\VScode Source\\python\\kepler_orbit\\ISS.txt'), tspan, dt, coes = True, deg = False)
    op1 = OP(t.tle2coes('C:\\Users\\soft_kitty\\VScode Source\\python\\kepler_orbit\\HST.txt'), tspan, dt, coes = True, deg = False)

    c0 = [cb['radius'] + 35800.0, 0.0,0.0,0.0,0.0,0.0 ,[2021,3,17,00]]
    op2 = OP(c0, tspan, dt, coes = True)

    t.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels = ['ISS', 'HST', 'GEO'], show_plot = True)