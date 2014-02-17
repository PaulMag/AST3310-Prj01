#include "rhs.hpp";
#include "helper_functions.hpp";
#include "constants.cpp";

/*
 * @return dR / dm
 */
double rhs_R(double rho) {
  return 1.0 / ( 4.0 * PI * r*r * rho );
}

/*
 * @return dP / dm
 */
double rhs_P(double m, double R) {
  return - ( G * m ) / ( 4.0 * PI * pow(R,4) );
}

/*
 * @return dL / dm
 */
double rhs_L(double T, double rho) {
  return epsilon(T,rho);
}

/*
 * @return dT / dm
 */
double rhs_T(double L, double T, double R, double rho) {
  double sigma = STEFAN_BOLTZMANN;
  return - ( 3.0 * kappa(T,R,rho) * L ) / ( 256.0 * PI*PI * sigma * pow(R,4) * pow(T,3) );
}
