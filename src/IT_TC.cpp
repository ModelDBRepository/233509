
#include "IT_TC.h" 
#include <math.h>

IT_TC::IT_TC(double v) {
	 K_Tau_h = 1;

	 G_Ca = 2; //2;///////////////////////////////////

     Shift_m = 0;
     Shift_h = 0;

     Phi_m = pow(Qm,((Cels-24)/10));
     Phi_h = pow(Qh,((Cels-24)/10));

  //   m0 = 1 / (1+exp(-(v+59 - Shift_m)/6.2));
  //   h0 = 1 / (1+exp((v+83 - Shift_h)/4.0));


     eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
}

double IT_TC::Ca_0 = 2, IT_TC::Cels = 36;
double IT_TC::Qm = 3.55, IT_TC::Qh = 3; //2.8;


void IT_TC::init(double v) {
	m0 = 1 / (1+exp(-(v+59 - Shift_m)/6.2));
    h0 = 1 / (1+exp((v+83 - Shift_h)/4.0));
}

void IT_TC::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.) {
   // 	printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
    }
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca); 

  m_inf = 1 / (1+exp(-(v+59 - Shift_m)/6.2));//////////////////////////////////////
  h_inf = 1 / (1+exp((v+83- Shift_h)/4.0));  // 4.0!!!;

  tau_m = (1/(exp(-(v+132 - Shift_m)/16.7)+exp((v+16.8 - Shift_m)/18.2)) + 0.612) / Phi_m;

  // tau_h = (30.8 + (211.4 + exp((v - Shift_h + 113.2)/5))/
  //         (1+exp((v - Shift_h + 84)/3.2))) / Phi_h;

  tau_h =  ( (exp((v+467 -Shift_h)/66.6)) ) / Phi_h;
  if(v >= (-80+Shift_h) )
      tau_h = ((exp(-(v+22 -Shift_h)/10.5)) +28) / Phi_h;

  if(v <= -40 )
      tau_h = K_Tau_h * tau_h;



  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}


