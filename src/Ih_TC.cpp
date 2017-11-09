
#include "Ih_TC.h"
#include <math.h>


Ih_TC::Ih_TC(double v) {
    G_H = 1; //2; //////////////////////////////////////////////
  //  Tad = pow(3,((Cels-23.5)/10));

    Tad = 1;
    Shift_m = 0;

    //   m0 = 1.0 / (1+exp((v+75)/5.5));
}

double Ih_TC::Cels = 36, Ih_TC::E_H = -43;


void Ih_TC::init(double v) {
	 m0 = 1.0 / (1+exp((v+75 - Shift_m)/5.5));
}


void Ih_TC::calc(double m, double &fm, double v){

  iH = G_H*m*(v - E_H);

  m_inf = 1.0 / (1+exp((v+75 - Shift_m)/5.5));
  tau_m = (1.0/( exp(-(14.59+0.086*(v-Shift_m))) + exp(-1.87+0.0701*(v-Shift_m)))) / Tad;

  fm = -(1/tau_m)*(m - m_inf);                                  

}
