#!/usr/bin/env python
"""
Created on 28-02-2014

Script is used for modeling the solar core via solving the four linked
differential equations.

@author Kristoffer Braekken
"""
from numpy import pi, log10, exp, zeros

"""CONSTANTS"""

_OPACITY_FILE = "../data/opacity.txt" # File with opacities

"""PHYSICAL CONSTANTS"""

_L_SUN = 3.846e26 # [W]
_R_SUN = 6.96e8 # [m]
_M_SUN = 1.989e30 # [kg]

_SIGMA = 5.67e-8 # [W m^-2 K^-4]
_K_B = 1.382e-23 # [m^2 kg s^-2 K^-1]
_N_A = 6.0221413e23 # Avogadro's constant

_H_MASS = 1.6738e-27 # [kg]
_HE3_MASS = 5.0081e-27 # [kg]
_HE4_MASS = 6.6464e-27 # [kg]
_LI7_MASS = 7.01600455 # [kg]
_BE7_MASS = 7.01692983 # [kg]
_E_MASS = 9.10938291 # [kg]

"""NUCLEAR ENERGY VALUES"""

# PP I
_Q_H_H = 1.177 # [MeV]
_Q_D_HE = 5.494 # [MeV]
_Q_HE3_HE3 = 12.860 # [MeV]

# PP II
_Q_HE3_ALPHA = 1.586 # [MeV]
_Q_BE7_E = 0.049 # [MeV]
_Q_LI7_H = 17.346 # [MeV]

# PP III
_Q_BE7_H = 0.137 # [MeV]
_Q_B8 = 8.367 # [MeV]
_Q_BE8 = 2.995 # [MeV]

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
    return - ( G * m ) / ( 4. * pi * r**4 )

def rhs_L(rho, T, compounds):
    """
    @return dL / dm
    """
    return epsilon(rho, T, compounds)

def rhs_T(T, rho, L, r):
    """
    @return dT / dm
    """
    return ( 3. * kappa(T, rho) * L ) / ( 256. * pi*pi * SIGMA * r**4 * T**3 )

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
    energy_chains['PPIII'] = (_Q_BE7_H+_Q_BE8_Q_B8) \
                    * compounds['Be7'].rate(rho,T,compounds['H'])

    eps = np.sum([energy_chains[key] for key in ['PPI','PPII','PPIII']])
    return eps, energy_chains

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
    
    # PP III
    elif (i.name == 'Be7' and j.name == 'H') or \
         (i.name == 'H' and j.name == 'Be7'):
        l = 3.11e5 * T9**(-2./3) * exp(-10.262 * T9**(-1./3)) \
                + 2.53e3 * T9**(-3./2) * exp(-7.306 * T9**(-1.))

    return _N_A * l

"""MAIN INTEGRATION PROCESS"""

def integrate_FE(dm, tol=1e-10):
    """
    Performs integration using the ForwardEuler scheme.
    """
    if dm > 0:
        print 'dm must be negative, returning.'
        return

    compounds = create_compounds()

    m = _M0
    rho = _RHO0
    N = int(abs(m / float(dm)))
    print N

    # Variable parameters
    r = zeros(N)
    r[0] = _R0

    P = zeros(N)
    P[0] = _P0

    L = zeros(N)
    L[0] = _L0

    T = zeros(N)
    T[0] = _T0

    for i in range(1,N):
        L[i] = L[i-1] + dm*rhs_L(rho, T[i-1], compounds)
        r[i] = r[i-1] + dm*rhs_r(r[i-1], rho)
        P[i] = P[i-1] + dm*rhs_P(r[i-1], m)
        T[i] = T[i-1] + dm*rhs_T(T[i-1], rho, L[i-1], r[i-1])

        m -= dm

        if (abs(m) < tol) or (abs(r[i]) < tol) or (abs(L[i]) < tol):
            print 'Integration complete before loop finished. Returning.'
            return r[:i+1], P[:i+1], L[:i+1], T[:i+1]

    return r, P, L, T

if __name__ == '__main__':
    import sys
    dm = float(sys.argv[1])

    integrate_FE(dm)
