import math
import numpy as np

import tools as t
import planetary_data as pd

pi = math.pi

def lamberts_universal_variables(r0, r1, deltat, tm = 1, mu = pd.earth['mu'], tol = 1e-6, max_steps = 200, psi = 0, psi_u = 4*pi**2, psi_l = -4 * pi):
    sqrt_mu = math.sqrt(mu)
    r0_norm = t.norm(r0)
    r1_norm = t.norm(r1)

    gamma = np.dot(r0, r1) / r0_norm /r1_norm

    beta = tm * math.sqrt(1 - gamma**2)
    
    A = tm * math.sqrt(r0_norm * r1_norm * (1+gamma))

    if A == 0:
        return np.array([0,0,0]), np.array([0,0,0])

    c2 = 0.5
    c3 = 1 / 6.0

    step = 0
    solved = False

    for n in range(max_steps):
        B = r0_norm + r1_norm + A *(psi * c3 - 1) / math.sqrt(c2)

        if A > 0.0 and B < 0.0:
            psi_l += pi
            B *= -1

        chi3 = math.sqrt(B / c2) ** 3
        deltat_ = (chi3 * c3 + A * np.sqrt(B)) / sqrt_mu

        if abs(deltat - deltat_) < tol:
            solved = True
            break

        if deltat_ <= deltat:
            psi_l = psi
        else:
            psi_u = psi

        psi = (psi_u + psi_l) / 2.0
        c2 = C2(psi)
        c3 = C3(psi)
        print(psi)

    if not solved:
        print('Lamberts UV variables did not converge')
        return np.array([0,0,0]), np.array([0,0,0])

    f = 1 - B / r0_norm
    g = A * math.sqrt(B / mu)
    gdot = 1 - B / r1_norm

    v0 = (r1 - f * r0) / g
    v1 = (gdot * r1 - r0) / g

    return v0, v1

def C2(psi):
        return (1 - math.cos(math.sqrt(psi))) / psi

def C3(psi):
        return (math.sqrt(psi) - math.sin(math.sqrt(psi))) / (psi * math.sqrt(psi))