#!/usr/bin/env python
"""
Created on 28-02-2014

Script is used for modeling the solar core via solving the four linked
differential equations.

@author Kristoffer Braekken
"""
import sys
import numpy as np

from numpy import pi, log10, exp, zeros

"""CONSTANTS"""

_OPACITY_FILE = "../data/opacity.txt" # File with opacities

"""PHYSICAL CONSTANTS"""

_L_SUN = 3.846e33 # [erg/s]
_R_SUN = 6.96e10 # [cm]
_M_SUN = 1.989e33 # [g]

_G = 6.67384e-8 # [cm^3 g^-1 s^-2]
_C = 3.e10 # [cm s^-1]
_SIGMA = 5.6704e-5 # Stefan-boltzmann [erg cm^-2 s^-1 K^-4]
_K_B = 1.3806488e-16 # [erg K^-1]
_N_A = 6.0221413e23 # Avogadro's constant [mol^-1]
_M_U = 1.66e-24 # [g]

_H_MASS = 1.6738e-24 # [g]
_HE3_MASS = 5.0081e-24 # [g]
_HE4_MASS = 6.6464e-24 # [g]
_LI7_MASS = 1.16503486e-23 # [g]
_BE7_MASS = 1.16518851e-23 # [g]
_E_MASS = 9.10938291e-28 # [g]

"""NUCLEAR ENERGY VALUES"""

# PP I
_Q_H_H = 1.88576182e-6 # [erg]
_Q_D_HE = 8.80235805e-6 # [erg]
_Q_HE3_HE3 = 2.06039906e-5 # [erg]

# PP II
_Q_HE3_ALPHA = 2.54105203e-6 # [erg]
_Q_BE7_E = 7.85066517e-8 # [erg]
_Q_LI7_H = 2.77977634e-5 # [erg]

# PP III
_Q_BE7_H = 2.77176546e-7 # [erg]
_Q_B8 = 1.34054113e-5 # [erg]
_Q_BE8 = 4.79851881e-6 # [erg]

"""INITIAL PARAMETERS"""

_L0 = _L_SUN
_R0 = 0.5*_R_SUN
_M0 = 0.7*_M_SUN
_RHO0 = 1 # [g cm^-3]
_T0 = 1.e5 # [K]
_P0 = 1.e12 # [g cm^-1 s^-2]

_X0 = 0.7
_Y3_0 = 1.e-10
_Y0 = 0.29
_Z0 = 0.01
_Z0_7LI = 1.e-5
_Z0_7BE = 1.e-5

"""IONIZATION"""

_MU0 = 1. / ( _X0 + _Y0 / 4. + _Z0 / 2. )
_E = _MU0 * ( _X0 + (1 + 2) * _Y0 / 4. )
_MU = _MU0 / (1 + _E)

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

    def rate(self, rho, T, other_compound):
        """
        @param rho The current mass density of the solar core.
        @param T The current temperature.
        @param other_compound The other compound to react with.
        @return Rate per unit mass.
        """
        #TODO is it correct that equal compounds is 1?
        if self.name == other_compound.name:
            kronecker_delta = 1
        else:
            kronecker_delta = 0

        return ( self.n(rho) * other_compound.n(rho) ) \
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
    return - ( _G * m ) / ( 4. * pi * r**4 )

def rhs_L(rho, T, compounds):
    """
    @return dL / dm
    """
    return epsilon(rho, T, compounds)

def rhs_T(T, rho, L, r):
    """
    @return dT / dm
    """
    return ( 3. * kappa(T, rho) * L ) / ( 256. * pi*pi * _SIGMA * r**4 * T**3 )

"""OTHER FUNCTIONS"""

def kappa(T, rho, opacityFile=_OPACITY_FILE):
    """
    Function reads from opacity file (defined as constant by default). And uses
    temperature T and density rho to find the kappa which is closest.

    @param rho Density in CGS units [g cm^-3]
    @param T Temperature in SI units [K]

    @return kappa in CGS units [cm^2 g^-1]
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
    return 10**log10kappa # [ cm^2 g^-1 ]

def epsilon(rho, T, compounds):
    """
    @param compunds List of compounds.
    @return Energy produced times reaction rate. Also, the dictionairy
    containing energy levels for each particular chain. Probably convenient for
    testing for incredible values.
    """
    energy_chains = {}

    # PP I, II, III
    energy_chains['all'] = _Q_H_H * compounds['H'].rate(rho,T,compounds['H'])

    # PP I
    energy_chains['PPI'] = _Q_HE3_HE3 * compounds['He3'].rate(rho,T,compounds['He3'])

    # PP II
    energy_chains['PPII'] = _Q_HE3_ALPHA * compounds['He3'].rate(rho,T,compounds['He4'])
    energy_chains['PPII'] += _Q_BE7_E * compounds['Be7'].rate(rho,T,compounds['e-'])
    energy_chains['PPII'] += _Q_LI7_H * compounds['Li7'].rate(rho,T,compounds['H'])

    # PP III
    energy_chains['PPIII'] = (_Q_BE7_H+_Q_BE8 + _Q_B8) \
                    * compounds['Be7'].rate(rho,T,compounds['H'])

    eps = sum([energy_chains[key] for key in ['PPI','PPII','PPIII']])
    return eps, energy_chains

def ideal(P, T):
    """
    Ideal equation of state.
    """
    a = 4. * _SIGMA / _C
    P_rad = ( a / 3. ) * T**4
    P_g = P - P_rad
    rho = P_g * _MU * _M_U / (_K_B * T)
    return rho

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
    compounds['e-'] = Compound( 'e-', _E_MASS, 1)

    # Special case, e relative density
    compounds['e-'].relative_density = sum(
            [c.relative_density for c in \
            [compounds[s] for s in \
            ['He3', 'He4', 'H'] ] ] )

    return compounds

def lambda_function(i, j, T):
    """
    @param i,j The two reacting compounds.
    @param T Temperature in core.
    @return Appropriate lambda based on reaction taking place.
    """
    T9 = T / 1.e9

    # PP I, II, III
    if i.name == 'H' and j.name == 'H':
        l = 4.01e-15 * T9**(-2./3) * exp( -3.380 * T9**(-1./3) ) \
                * (1 + 0.123 * T9**(1./3) + 1.09 * T9**(2./3) + \
                0.938 * T9 )

    # PP I
    elif i.name == 'He3' and j.name == 'He3':
        l = 6.04e10 * T9**(-2./3) * exp(-12.276 * T9**(-1./3)) \
                * (1 + 0.034 * T9**(1./3) - 0.522 * T9**(2./3) \
                - 0.124 * T9 + 0.353 * T9**(4./3) \
                + 0.213 * T9**(-5./3) )

    # PP II, III
    elif (i.name == 'He3' and j.name == 'He4') or \
         (i.name == 'He4' and j.name == 'He3'):
        T_star = T9 / (1. + 4.95e-2 * T9)

        l = 5.61e6 * T_star**(5./6) * T9**(-3./2) \
                * exp(-12.826 * T_star**(-1./3))

    # PP II
    elif (i.name == 'Be7' and j.name == 'e-') or \
         (i.name == 'e-' and j.name == 'Be7'):
        l = 1.34e10 * T9**(-1./2) * (1 - 0.537 * T9**(1./3) \
                + 3.86 * T9**(2./3) + 0.0027 * T9**(-1.) \
                * exp(2.515e-3 * T9**(-1.) ) )

    elif (i.name == 'Li7' and j.name == 'H') or \
         (i.name == 'H' and j.name == 'Li7'):
        T_star = T9 / ( 1 + 0.759*T9)

        l = 1.096e9 * T9**(-2./3) * exp(-8.472 * T9**(-1./3)) \
                - 4.830e8 * T_star**(5./6) * T9**(-2./3) \
                * exp(-8.472 * T_star**(-1./3)) \
                + 1.06e10 * T9**(-3./2) * exp(-30.442 * T9**(-1.))
    
    # PP III
    elif (i.name == 'Be7' and j.name == 'H') or \
         (i.name == 'H' and j.name == 'Be7'):
        l = 3.11e5 * T9**(-2./3) * exp(-10.262 * T9**(-1./3)) \
                + 2.53e3 * T9**(-3./2) * exp(-7.306 * T9**(-1.))

    return l / _N_A

"""MAIN INTEGRATION PROCESS"""

def integrate_FE(dm, tol=1e-10):
    """
    Performs integration using the ForwardEuler scheme.
    """
    if dm > 0:
        print 'dm must be negative, returning.'
        return

    compounds = create_compounds()

    N = int(abs(_M0 / float(dm)))
    m = np.arange(N)*abs(dm)
    m = m[::-1] # Reverse, start at outside

    # Variable parameters
    r = zeros(N)
    r[0] = _R0

    P = zeros(N)
    P[0] = _P0

    L = zeros(N)
    L[0] = _L0

    T = zeros(N)
    T[0] = _T0

    rho = zeros(N)
    rho[0] = _RHO0

    # HACK: For transmitting values to datawriter
    initials = {}
    initials['M0'] = _M0
    initials['R0'] = _R0
    initials['P0'] = _P0
    initials['T0'] = _T0
    initials['RHO0'] = _RHO0
    initials['L0'] = _L0

    print rho[0]
    for i in range(1,N):
        # Print progress
        # percent = (i / float(N-1)) * 100
        # sys.stdout.write('Progress: %4.2f %s\r' % (percent, '%'))
        # sys.stdout.flush()

        rho[i] = ideal(P[i-1], T[i-1])
        print rho[i]

        L_new = rhs_L(rho[i], T[i-1], compounds)
        L[i] = L[i-1] + dm*L_new[0]
        # print L_new[1]

        r[i] = r[i-1] + dm*rhs_r(r[i-1], rho[i])
        P[i] = P[i-1] + dm*rhs_P(r[i-1], m[i-1])
        T[i] = T[i-1] + dm*rhs_T(T[i-1], rho[i], L[i-1], r[i-1])

        if (abs(m[i-1]) < tol) or (abs(r[i]) < tol) or (abs(L[i]) < tol):
            print 'Integration complete before loop finished. Returning.'
            return r[:i+1], m[:i+1], P[:i+1], L[:i+1], T[:i+1], rho[:i+1], initials

    return r, m, P, L, T, rho, initials

if __name__ == '__main__':
    dm = float(sys.argv[1])
    if abs(dm) <= 1e23:
        print 'Too tiny dm, try higher than -1e23.'
        sys.exit(1)

    r, m, P, L, T, initials = integrate_FE(dm)
