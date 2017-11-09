
#include "ICaL_TC.h"
#include <math.h>

ICaL_TC::ICaL_TC(double v) {
     G_CaL = 1; //2;///////////////////////////////////

     Shift_m = 0;
     Shift_h = 0;

     Phi_m = pow(Qm,((Cels-24)/10));
     Phi_h = pow(Qh,((Cels-24)/10));

  //   Phi_m = 1;
  //   Phi_h = 1;

  //   m0 = 1 / (1+exp(-(v+10 - Shift_m)/4));
  //   h0 = 1 / (1+exp((v+25 - Shift_h)/2));

     eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
  // eca = 120;
}

double ICaL_TC::Ca_0 = 2, ICaL_TC::Cels = 36;
double ICaL_TC::Qm = 3.55, ICaL_TC::Qh = 3; //2.8;


void ICaL_TC::init(double v) {
    m0 = 1 / (1+exp(-(v+10 - Shift_m)/4));
    h0 = 1 / (1+exp((v+25 - Shift_h)/2));

}

void ICaL_TC::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.) {
   // 	printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
    }


  eca = eca0 * log(ratio);

  iCaL = G_CaL*m*m*h*(v - eca);

  m_inf = 1 / (1+exp(-(v+10 - Shift_m)/4));
  h_inf = 1 / (1+exp((v+25 - Shift_h)/2));

  tau_m = (0.4 + 0.7/(exp(-(v+5-Shift_m)/15)+exp((v+5-Shift_m)/15)) ) / Phi_m;

  tau_h = (300 + 100/(exp(-(v+40-Shift_h)/9.5) + exp((v+40-Shift_h)/9.5))) / Phi_h;


  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}


