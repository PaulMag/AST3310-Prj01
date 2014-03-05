#!/usr/bin/env python
"""
Created on 28-02-2014

Script is used for modeling the solar core via solving the four linked
differential equations.

@author Kristoffer Braekken
"""
import numpy as np

from numpy import pi, log10

"""CONSTANTS"""

_OPACITY_FILE = "../data/opacity.txt" # File with opacities

_L_SUN = 3.846e26 # [W]
_R_SUN = 6.96e8 # [m]
_M_SUN = 1.989e30 # [kg]
_SIGMA = 5.67e-8 # [W m^-2 K^-4]
_K_B = 1.382e-23 # [m^2 kg s^-2 K^-1]
# TODO add constants for particle masses

"""INITIAL PARAMETERS"""

_L0 = _L_SUN
_R0 = 0.5*_R_SUN
_M0 = 0.7*_M_SUN
_RHO0 = 1.e3 # [kg m^-1]
_T0 = 1.e5 # [K]
_P0 = 1.e11 # [Pa]

_X0 = 0.7
_Y3_0 = 1.e-10
_Y0 = 0.29
_Z0 = 0.01
_Z0_7LI = 1.e-5
_Z0_7BE = 1.e-5

_H_MASS = 1.6738e-27 # [kg]
_HE3_MASS = 5.0081e-27 # [kg]
_HE4_MASS = 6.6464e-27 # [kg]
_LI7_MASS = 7.01600455 # [kg]
_BE7_MASS = 7.01692983 # [kg]

"""CLASSES"""

class Compound(object):
    """
    Describes a compound. Used in conjunction with calculating energy
    production rates so that the kronecker delta can be decided generally.

    Keeps ratios and mass of one particle in order.

    @field rho Particle density.
    """
    def __init__(self, name, mass, ratio):
        """
        @param name Identifier.
        @param mass Mass of one particle.
        @param ratio How much of the total mass is this compound.
        """
        self.name,self.m,self.r = name, mass, ratio
        self.relative_density = self.m / self.r

    def r(self, rho, T, other_compound):
        """
        @param rho The current mass density of the solar core.
        @param T The current temperature.
        @param other_compound The other compound to react with.
        ®return Rate per unit mass.
        """
        #TODO is it correct that equal compounds is 1?
        if self.name == other_compound.name:
            kronecker_delta = 1
        else:
            kronecker_delta = 0

        return ( self.r(rho) * other_compound.r(rho) ) \
                        / ( rho * (1 + kronecker_delta) ) * \
                        lambda_function(self, other_compound, T)

    def n(self, rho):
        """
        @return Number density of this compound.
        """
        return rho*self.relative_density

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

def kappa(T, rho, opacityFile=_OPACITY_FILE):
    """
    Function reads from opacity file (defined as constant by default). And uses
    temperature T and density rho to find the kappa which is closest.

    @return kappa in unit [cm^2 g^-1]
    """
    # log10(T)
    log10T = log10(T)

    # log10(R)
    R = float(rho) / ( T / 1.e6 )
    log10R = log10(R)

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
    return 10**log10kappa

def opacity_test(tol=1.e-10):
    """
    Function for testing that the opacity is fetched correctly.
    """
    # Test values
    T = 10**(5.) # Feth 5.00 row
    rho = 1.e-6 # Fetch -5.0 column

    ans = log10(kappa(T, rho))
    if abs(ans - (-0.068)) < tol:
        print 'Sucess.'
    else:
        print 'Fail.\n10**kappa =', ans, 'and not -0.068.'

def create_compounds():
    """
    @return List of compounds from starting parameters.
    """
    compounds = {}

    compounds['H'] = Compound( 'H', _H_MASS, _X0 )
    compounds['He3'] = Compound( 'He3', _HE3_MASS, _Y3_0 )
    compounds['He4'] = Compound( 'He4', _HE3_MASS, _Y0 - _Y3_0 )
    compounds['Li7'] = Compound( 'Li7', _HE3_MASS, _Z0_7BE )
    compounds['Be7'] = Compound( 'Be7', _HE3_MASS, _Z0_7LI )

    return compounds

"""MAIN INTEGRATION PROCESS"""

def integrate(dm):
    """
    Performs integration.
    """
    #TODO Write function (skeleton)
    pass

if __name__ == '__main__':
    opacity_test()
