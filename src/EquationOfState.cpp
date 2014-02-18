#include "EquationOfState.hpp"

EquationOfState :: EquationOfState( double (*getRho)(double T, double P), 
                                    double (*getP)(double rho, double T),
                                    double (*getT)(double rho, double P) )
{
  getRho = getRho;
  getP = getP;
  getT = getT;
}
