#!/usr/bin/env python
"""
Created on Fri 7 March 2014

For backing up SI constants before converting them to CGS units in main
program.

@author Kristoffer Braekken
"""

"""PHYSICAL CONSTANTS"""

_L_SUN = 3.846e26 # [W]
_R_SUN = 6.96e8 # [m]
_M_SUN = 1.989e30 # [kg]

_G = 6.67384e-11 # [m^3 kg^-1 s^-2]
_C = 3.e8 # [m s^-1]
_SIGMA = 5.67e-8 # [W m^-2 K^-4]
_K_B = 1.382e-23 # [m^2 kg s^-2 K^-1]
_N_A = 6.0221413e23 # Avogadro's constant

_H_MASS = 1.6738e-27 # [kg]
_HE3_MASS = 5.0081e-27 # [kg]
_HE4_MASS = 6.6464e-27 # [kg]
_LI7_MASS = 1.16503486e-26 # [kg]
_BE7_MASS = 1.16518851e-26 # [kg]
_E_MASS = 9.10938291e-31 # [kg]

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

"""IONIZATION"""

_MU0 = 1. / ( _X0 + _Y0 / 4. + _Z0 / 2. )
_E = _MU0 * ( _X0 + (1 + 2) * _Y0 / 4. )
_MU = _MU0 / (1 + _E)
