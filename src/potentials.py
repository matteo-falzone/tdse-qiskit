"""
potentials.py
-------------
Collection of potential energy functions V(x) for 1D quantum simulations.
"""

import numpy as np

#Free particle
def free(x):
    return np.zeros_like(x)

#Harmonic oscillator
def harmonic(x, omega=1, m=1):
    return 0.5 * m * omega**2 * x**2

#Rectangular potential barrier
def barrier(x, V0, a):
    return np.where(np.abs(x) < a, V0, 0.0)