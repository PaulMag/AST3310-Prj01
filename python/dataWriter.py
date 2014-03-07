#!/usr/bin/env python
"""
Created on Thu 6 March 2014

Contains methods for writing datafiles from simulation data from the
`SolarCoreModel.py` script.

@author Kristoffer Braekken
"""
import sys,os,SolarCoreModel
from time import gmtime, strftime

"""CONSTANTS"""

_DATAPATH = '../data/'

"""FUNCTIONS"""

def dataTable(r_arr, m_arr, P_arr, L_arr, T_arr, rho_arr, initial_params, outputFile):
    """
    Writes a neat table to the given file.

    @param m,r,P,T,rho,L Arrays of all the values.
    @param initial_params Dictionairy containing initial values. Necessary keys
    are M0, r0, P0, T0, rho0, L0 (values at outer edge).
    @param outputFile Filename to write to.
    """
    if os.path.exists(outputFile):
        print 'File already exists.'
        return

    ofile = open(outputFile, 'w')

    r_arr /= initial_params['R0']
    m_arr /= initial_params['M0']
    L_arr /= initial_params['L0']

    for r,m,P,L,T,rho in zip(r_arr,m_arr,P_arr,L_arr,T_arr,rho_arr):
        ofile.write('%12.4f %12.4f %12.3e %12.3e %12.3e %12.4f\n' % (m,r,P,T,rho,L))

    ofile.close()

if __name__ == '__main__':
    if not os.path.exists(_DATAPATH):
        os.mkdir(_DATAPATH)

    dm = float(sys.argv[1])

    r,m,P,L,T,rho,init = SolarCoreModel.integrate_FE(dm)
    dataTable( r,m,P,L,T,rho,init,
            os.path.join( _DATAPATH, strftime('table_%d-%m-%Y_at_%H:%M:%S.dat',gmtime() ) ) )
