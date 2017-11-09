
#include "IKs_TC.h"
#include <math.h>


IKs_TC::IKs_TC(double v) {
    G_Ks = 1; //2; //////////////////////////////////////////////

    Tad = pow(3,((Cels-23.5)/10));
   // Tad = 1;

    m0 = 1.0 / (1+exp(-(v+43)/17));
    h0 = 1.0/(1+exp((v+58)/10.6)); }

double IKs_TC::Cels = 36, IKs_TC::E_K = -80;

void IKs_TC::calc(double m, double h, double &fm, double &fh, double v){
  iKs = G_Ks*m*h*(v - E_K);

  m_inf = 1.0 / (1+exp(-(v+43)/17));
  h_inf = 1.0/(1+exp((v+58)/10.6));

  tau_m = (1.0/( exp((v-81)/25.6)+exp(-(v+132)/18) )+ 9.9) / Tad;

  tau_h = (1.0/(exp((v-1329)/200)+exp(-(v+130)/7.1))+ 120) / Tad;

  fm = -(1/tau_m)*(m - m_inf);                                  
  fh = -(1/tau_h)*(h - h_inf);
}
