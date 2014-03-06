#!/usr/bin/env python
"""
Created on Thu 6 March 2014

Contains testing routines for `SolarCoreModel.py`.

@author Kristoffer Braekken
"""
import SolarCoreModel

from numpy import log10

def opacity_test(tol=1.e-10):
    """
    Function for testing that the opacity is fetched correctly.
    """
    # Test values
    T = 10**(5.) # Feth 5.00 row
    rho = 1.e-6 # Fetch -5.0 column
    rho /= 1.e3; rho *= 1./1e6 # Convert to SI units [kg m^-3]

    ans = log10(SolarCoreModel.kappa(T, rho))
    if abs(ans - (-0.068)) < tol:
        print 'Sucess.'
    else:
        print 'Fail.\n10**kappa =', ans, 'and not -0.068.'

if __name__ == '__main__':
    opacity_test()
