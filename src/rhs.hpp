/*
 * @return dR / dm
 */
double rhs_R(double rho, double R);

/*
 * @return dP / dm
 */
double rhs_P(double m, double R);

/*
 * @return dL / dm
 */
double rhs_L(double T, double rho);

/*
 * @return dT / dm
 */
double rhs_T(double L, double T, double R, double rho);
