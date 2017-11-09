
#include "IA_TC.h" 
#include <math.h>


IA_TC::IA_TC(double v) {
    G_A = 1; //2; //////////////////////////////////////////////

    Tad = pow(3,((Cels-23.5)/10));
   // Tad = 1;

    m0 = 1.0 / (1+exp(-(v+36)/8.5));
    h0 = 1.0/(1+exp((v+78)/6)); } 

double IA_TC::Cels = 36, IA_TC::E_K = -80;

void IA_TC::calc(double m, double h, double &fm, double &fh, double v, double x){
  iA = G_A*m*m*m*m*h*(v - E_K);

  m_inf = 1.0 / (1+exp(-(v+36)/8.5));
  tau_m = (1.0/( exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7) ) +0.37) / Tad;

  h_inf = 1.0/(1+exp((v+78)/6));
  tau_h = (1.0/(exp((v+46.05)/5)+exp(-(v+238.4)/37.45)) ) / Tad;
  if(v >= -63) 
      tau_h = 19.0/Tad;



  fm = -(1/tau_m)*(m - m_inf);                                  
  fh = -(1/tau_h)*(h - h_inf);
}
