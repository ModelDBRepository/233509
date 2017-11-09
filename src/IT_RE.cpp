
#include "IT_RE.h" 
#include <math.h>
#include <iostream>


IT_RE::IT_RE(double v) {
    G_Ca = 1.75;
    Qm = 5; //2.5; //3; 
    Qh = 3; //2.5; //5;
    Shift = 0;

    Phi_m = pow(Qm,((Cels-24)/10));
    Phi_h = pow(Qh,((Cels-24)/10));


   // m0 = 1/(1 + exp(-(v + 52 - Shift)/7.4));
   // h0 = 1/(1 + exp((v + 80 - Shift)/5));


    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);

    ECA = 120;
    } 

double IT_RE::Ca_0 = 2, IT_RE::Cels = 36;


void IT_RE::init(double v) {

	 m0 = 1/(1 + exp(-(v + 52 - Shift)/7.4));
	 h0 = 1/(1 + exp((v + 80 - Shift)/5));
}

void IT_RE::calc(double m, double h, double &fm, double &fh, double v, double cai, double x) {
  ratio = Ca_0/cai;
  if(ratio <= 0.) {
  	//printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  }
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca);

 //  iT = G_Ca*m*m*h*(v - ECA);

  m_inf = 1/(1 + exp(-(v + 52 - Shift)/7.4));
  tau_m = 3 + 1/(exp((v + 27 - Shift)/10) + exp(-(v + 102 - Shift)/15));

  h_inf = 1/(1 + exp((v + 80 - Shift)/5));
  tau_h = 85 + 1/(exp((v + 48 - Shift)/4) + exp(-(v + 407 - Shift)/50));
  
  fm = -Phi_m*(m - m_inf)/tau_m;
  fh = -Phi_h*(h - h_inf)/tau_h;
}
