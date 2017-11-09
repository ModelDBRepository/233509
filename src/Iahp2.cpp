
#include "Iahp2.h"
#include <math.h>
#include "../Constants.h"

Iahp2::Iahp2(double v, double ca) {
    G_AHP = 1;
    Tad = 1;
    Shift_m = 0;

}

double Iahp2::Cels = 36, Iahp2::E_K = -80;

void Iahp2::init(double v, double ca) {

    m0 = 48*ca*ca/(48*ca*ca+0.09);
}

void Iahp2::calc(double m, double &fm, double ca, double v){

  iAHP = G_AHP*m*(v - EK);

  m_inf =  48*ca*ca/(48*ca*ca + 0.09);
  tau_m = 1/(48*ca*ca + 0.09);

  fm = -(1/tau_m)*(m - m_inf);
 // fm = alpha*(1-m)-beta*m;

}
