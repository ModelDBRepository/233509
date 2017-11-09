
#include "Iahp_TC.h"
#include <math.h>
#include "../Constants.h"

Iahp_TC::Iahp_TC(double v, double ca) {
    G_AHP = 1; //2; //////////////////////////////////////////////

    Tad = pow(3,((Cels-23.5)/10));
    //Tad = 1;

    fV  = 500*exp((v-65)/27);
    fCa = (ca-0.015)/(1-exp(-(ca-0.015)/0.0013));
    alpha = fV*fCa;
    beta = 0.05;
    m0 = alpha / (alpha+beta);
    Shift_m = 0;
}

double Iahp_TC::Cels = 36, Iahp_TC::E_K = -80;

void Iahp_TC::calc(double m, double &fm, double ca, double v){

  iAHP = G_AHP*m*(v - EK);

  fV  = 500*exp((v-65)/27);
  fCa = (ca-0.015)/(1-exp(-(ca-0.015)/0.0013));

  if ((ca-0.015)<1e-5) {
  fCa = 1/(1/0.0013 + 0.5*(0.015-ca)/(0.0013*0.0013));
  }


  alpha = fV*fCa;
  beta = 0.05;

  fm = alpha*(1-m)-beta*m;

}
