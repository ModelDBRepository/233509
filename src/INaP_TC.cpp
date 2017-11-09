
#include "INaP_TC.h"
#include <math.h>


INaP_TC::INaP_TC(double v) {
    G_NaP = 1; //2; //////////////////////////////////////////////
    Tad = pow(3,((Cels-23.5)/10));
   // Tad = 1;
    m0 = 1.0 / (1+exp( -(v+49)/5.0));
}

double INaP_TC::Cels = 36, INaP_TC::E_NaP = 50;

void INaP_TC::calc(double m, double &fm, double v){

  iNaP = G_NaP*m*(v - E_NaP);
  m_inf = 1.0 / (1+exp( -(v+49)/5.0));

  alpha = 0.091*(v + 38) /(1 - exp(-(v+38)/5));
  beta  = -0.062*(v + 38)/(1 - exp((v + 38)/5));

  tau_m = 1/(alpha + beta) / Tad;
  fm = -(1/tau_m)*(m - m_inf);                                  

}
