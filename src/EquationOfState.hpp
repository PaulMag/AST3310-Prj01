class EquationOfState {
  public:
    EquationOfState( double (*getRho)(double T, double P), 
                     double (*getP)(double rho, double T),
                     double (*getT)(double rho, double P) );
  private:
    double getRho(double T, double P);
    double getP(double rho, double T);
    double getT(double rho, double P);
};
