#!/usr/bin/env python
"""
Created on 28-02-2014

Script is used for modeling the solar core via solving the four linked
differential equations.

@author Kristoffer Braekken
"""
import numpy as np

from numpy import pi

"""CONSTANTS"""

OPACITY_FILE = "../data/opacity.txt" # File with opacities

"""RHS FUNCTIONS"""

def rhs_r(r, rho):
    """
    @return dr / dm
    """
    return 1. / ( 4. * pi r*r * rho )

def rhs_P(r, m):
    """
    @return dP / dm
    """
    return - ( G * m ) / ( 4. * pi * r**4 )

def rhs_L():
    """
    @return dL / dm
    """
    return epsilon()

def rhs_T():
    """
    @return dT / dm
    """
    return ( 3. * kappa() * L() ) / ( 256. * pi*pi * SIGMA * r**4 * T**3 )

"""OTHER FUNCTIONS"""

def kappa(T, rho, opacityFile=OPACITY_FILE):
    """
    Function reads from opacity file (defined as constant by default). And uses
    temperature T and density rho to find the kappa which is closest.

    @return kappa in unit [cm^2 g^-1]
    """
    # log10(T)
    log10T = log(T)

    # log10(R)
    R = float(rho) / ( T / 1e6 )
    log10R = log(R)

    # Traverse file to find correct log10kappa
    opacities = open(opacityFile, 'r')

    """
    First, traverse header line to find correct column by finding the
    log10(T) which is closest.

    Exclude the first input ( [1:] ) since that is text.
    """
    headerline = \
            np.asarray( opacities.readline().strip().split()[1:], dtype=np.float64 )

    closest = headerline[0]
    for T in headerline[1:]:
        if abs(log10T-T) > abs(log10T-closes):
            closest = T

    # Which column does this T correspond to?

    # Returned from file is log10kappa
    return 10**log10kappa
