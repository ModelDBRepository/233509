
// From Inoue and Strowbridge, 2008

#include "Ican_TC.h"
#include <math.h>


Ican_TC::Ican_TC(double v) {
    G_CAN = 1; //2; //////////////////////////////////////////////
    Tad = 1;
    Shift_m = 0;

    Kd = 0.2;  // mM

   // m0 = 1.0/(1+exp(-(v+43)/5.2));

}

double Ican_TC::Cels = 36, Ican_TC::E_CAN = 10;

void Ican_TC::init(double v) {
	m0 = 1.0/(1+exp(-(v+43)/5.2));
}

void Ican_TC::calc(double m, double &fm, double ca, double v){

  fCa = ca/(ca+Kd);

  iCAN = G_CAN*fCa*m*(v - E_CAN);

  m_inf = 1.0/(1+exp(-(v+43)/5.2));

  tau_m = (  2.7/(exp(-(v+55)/15)+exp((v+55)/15)) + 1.6  ) / Tad;

  fm = -(1/tau_m)*(m - m_inf);                                  

}
