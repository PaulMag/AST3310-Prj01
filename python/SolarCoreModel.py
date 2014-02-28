#!/usr/bin/env python
"""
Created on 28-02-2014

Script is used for modeling the solar core via solving the four linked
differential equations.

@author Kristoffer Braekken
"""
import numpy as np

from numpy import pi, log

"""CONSTANTS"""

OPACITY_FILE = "../data/opacity.txt" # File with opacities

"""RHS FUNCTIONS"""

def rhs_r(r, rho):
    """
    @return dr / dm
    """
    return 1. / ( 4. * pi * r*r * rho )

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
    R = float(rho) / ( T / 1.e6 )
    log10R = log(R)

    # Begin process of traversing file to find correct log10kappa
    opacities = open(opacityFile, 'r')

    """
    First, traverse header line to find correct column by finding the
    log10(R) which is closest.

    Exclude the first input ( [1:] ) since that is text.
    """
    headerline = \
            np.asarray( opacities.readline().strip().split()[1:], dtype=np.float64 )

    opacities.close()

    closest = headerline[0]
    for R in headerline[1:]:
        if abs(log10R-R) < abs(log10R-closest):
            closest = R

    # Which column does this R correspond to?
    columnNo = 0
    for R in headerline:
        if R == closest:
            break
        columnNo += 1

    # Using the log10R, can now find the correct row to use
    opacities = open(opacityFile, 'r')
    opacities.readline() # Skip first line

    closest = False
    for line in opacities:
        if len(line.strip()) == 0: continue # Skip empty lines

        row = \
            np.asarray( line.strip().split(), dtype=np.float64 )

        if closest == False:
            closest = row[0]

        if abs(log10T-row[0]) < abs(log10T-closest):
            closest = row[0]

    opacities.close()

    # Which row does this log(R) correspond to?
    opacities = open(opacityFile, 'r')
    opacities.readline() # Skip first line

    rowNo = 0
    for line in opacities:
        if len(line.strip()) == 0: continue # Skip empty lines

        row = \
            np.asarray( line.strip().split(), dtype=np.float64 )

        if row[0] == closest: break

        rowNo += 1

    opacities.close()

    # Fetch correct value
    opacities = open(opacityFile, 'r')
    opacities.readline() # Skip first line

    current_row = 0
    current_column = 0
    for line in opacities:
        if len(line.strip()) == 0: continue # Skip empty lines

        row = \
            np.asarray( line.strip().split(), dtype=np.float64 )

        if current_row == rowNo:
            for k in row[1:]:
                if current_column == columnNo:
                    log10kappa = k
                    break
                current_column += 1
            break
        current_row += 1

    opacities.close()

    # Returned from file is log10kappa
    print log10kappa
    return 10**log10kappa

if __name__ == '__main__':
    T = 10**(3.75)
    rho = 5.623413252e-8
    kappa(T, rho)
