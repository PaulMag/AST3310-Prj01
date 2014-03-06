#!/usr/bin/env python
"""
Created on Thu 6 March 2014

Contains methods for writing datafiles from simulation data from the
`SolarCoreModel.py` script.

@author Kristoffer Braekken
"""
import sys,os,SolarCoreModel

"""CONSTANTS"""

_DATAPATH = '../data/'

"""FUNCTIONS"""

def dataTable(r_arr, m_arr, P_arr, L_arr, T_arr, rho_arr, initial_params, outputFile):
    """
    Writes a neat table to the given file. Table has following shape:
    _______________________________________________________
    m / M | r / R | P [Pa] | T [K] | rho [kg / m^3] | L / L
    _______________________________________________________
    0         0        0       0      0                0
    0         0        0       0      0                0
    0         0        0       0      0                0

    @param m,r,P,T,rho,L Arrays of all the values.
    @param initial_params Dictionairy containing initial values. Necessary keys
    are M0, r0, P0, T0, rho0, L0 (values at outer edge).
    @param outputFile Filename to write to.
    """
    if os.path.exists(outputFile):
        print 'File already exists.'
        return

    ofile = open(outputFile, 'w')
    ofile.write('______________________________________________\n')
    ofile.write('m / M0 | r / R0 | P [Pa] | T [K] | rho [kg m^-3] | L / L0\n')
    ofile.write('______________________________________________\n')

    r_arr /= initial_params['R0']
    m_arr /= initial_params['M0']
    L_arr /= initial_params['L0']

    for r,m,P,L,T,rho in zip(r_arr,m_arr,P_arr,L_arr,T_arr,rho_arr):
        ofile.write('%7f.4 %7f.4 %7g %7g %7g 7f.4\n' % (m,r,P,T,rho,L))

    ofile.close()

if __name__ == '__main__':
    if not os.path.exists(_DATAPATH):
        os.mkdir(_DATAPATH)

    dm = float(sys.argv[1])

    r,m,P,L,T,rho,init = SolarCoreModel.integrate_FE(dm)
    dataTable( r,m,P,L,T,rho,init, os.path.join(_DATAPATH,'table.dat') )
