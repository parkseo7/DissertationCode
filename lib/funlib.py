from __future__ import division, print_function

import os
import numpy as np
from numpy.polynomial.polynomial import polyval
import math
from math import pi


# ORDER PARAMETER

def order(Y):
    '''
    Returns the one-dimensional array of order values from Y.
    '''
    
    # Take the exponential
    expY = np.exp(1j*Y)
    N = Y.shape[1]
    
    # Summation
    orderY = np.abs(np.sum(expY, 1) / N)
    
    return orderY


# 2D ANALYSIS

def Omega2D(parameters):
    '''
    Determines the fixed-point equations for Omega, Delta, given initial guesses
    Omega0, Delta0.
    '''
    
    # Parameters
    g = parameters['g']
    w0 = parameters['omega0']
    gain = parameters['gain']
    tau0 = parameters['tau0']
    
    # Here Delta = Delta_12
    Delta_fun = lambda u: np.arcsin((w0 - u)/g)
    
    # Fixed-point equation for Omega:
    Omega_fun = lambda u: u - w0 - g*np.sin(-u*(tau0 + gain*(w0 - u)/g) + Delta_fun(u))
    
    return Omega_fun


def eig2D_det(Omega, Delta, parameters):
    '''
    Returns the 2x2 determinant complex eigenvalue criterion, in the form of
    an exponential polynomial. Here, Omega, Delta are (one of the) solutions
    to the fixed-point equation given by Omega2D. We assume that Delta > 0.
    '''
    
    # Parameters
    g = parameters['g']
    w0 = parameters['omega0']
    gain = parameters['gain']
    tau0 = parameters['tau0']
    
    # Defined parameters
    k = 1 - Omega*gain*np.cos(Delta)
    tauE = tau0 + gain*np.sin(Delta)
    C1 = g*np.cos(-Omega*tauE + Delta)
    C2 = g*np.cos(Delta)
    
    # Polynomials
    P3 = 1
    P2 = 1 + C1 + C2
    P1 = C1*k + C2*(1+C1)
    P0 = C1*C2
    
    Q1 = -C1*C2
    Q0 = -C1*C2
    
    poly = lambda z: (P3*z**3 + P2*z**2 + P1*z + P0) + (Q1*z + Q0)*np.exp(-z*tauE)
    
    return poly


def eig2D_cubic(Omega, Delta, parameters):
    '''
    Returns the coefficients of the quartic polynomial P(z) + Q(z) from the
    exponential polynomial in our eigenvalue equation.
    '''
    
    # Parameters
    g = parameters['g']
    w0 = parameters['omega0']
    gain = parameters['gain']
    tau0 = parameters['tau0']
    
    # Defined parameters
    k = Omega*gain*np.cos(Delta)
    tauE = tau0 + gain*Delta
    C_12 = g*np.cos(-Omega*tauE + Delta)
    C_21 = g*np.cos(Delta)
    C = C_12 + C_21
    C_2 = C_12*C_21
    
    b_3 = 1
    b_2 = 1 + C_12 + C_21
    b_1 = C_12*(1-k) + C_21
    b_0 = 0
    
    return np.array([b_3, b_2, b_1, b_0])


def quadratic_roots(coeffs):
    '''
    Returns the two branch roots of the quadratic using the coefficients.
    '''
    
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]
    
    r1 = (-b + np.sqrt(b**2 - 4*a*c, dtype='complex64'))/(2*a)
    r2 = (-b - np.sqrt(b**2 - 4*a*c, dtype='complex64'))/(2*a)
    
    return (r1, r2)


# ND ANALYSIS

def Omega_infty(Omega, delta, parameters, L=pi, steps=100):
    '''
    Computes the right-side integral as a Riemann sum with N steps,
    at delay tau, and a Gaussian distribution of differences at mean 0 and
    variance sigma^2.
    '''
    
    w0 = parameters['omega0']
    g = parameters['g']
    gain = parameters['gain']
    tau0 = parameters['tau0']
    
    N = steps
    z0 = np.zeros(N)
    delta2 = delta**2
    Delta = -L + 2*L*np.arange(N) / N
    
    if delta2 == 0:
        return w0 + np.sin(-Omega*tau0) 
    else:
        gauss = ((np.sqrt(2*pi*delta2))**-1)*np.exp(-Delta**2 / (2*delta2))
        sin_arr = np.sin(-Omega*np.maximum(tau0 + gain*Delta, z0) + Delta)*gauss
        return w0 + g*2*L*np.sum(sin_arr) / N


# SUPPLEMENTARY

def sign_log(x):
    '''
    Given any float x, returns the sign log scale of x.
    '''
    
    return np.sign(x)*np.log(1 + np.abs(x))


if __name__ == '__main__':
    pass
    