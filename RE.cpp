
#include <stdio.h>
#include "Constants.h"
#include "src/INaK.h"
#include "src/Iahp2.h"
#include "src/Ican_TC.h"

//#include "src/ICa.h"
#include "RE.h"


RE::RE():IT_RE(V0), INaK(V0), Ca_RE(), Iahp2(V0, Cai0), Ican_TC(V0) {
        E_l = -70;
        G_l = 0.05;
        G_kl = 0.0;
        S_RE = 1.43e-4;
        I_Stim = 0;
        }


void RE::init(double *y) {
       INaK::init(V0);
       IT_RE::init(V0);
       Iahp2::init(V0, Cai0);
       Ican_TC::init(V0);

        y[0] = V0;
	    y[1] = Cai0;
        y[2] = IT_RE::m0;
        y[3] = IT_RE::h0;
        y[4] = INaK::m0;
        y[5] = INaK::h0;
        y[6] = INaK::n0;
        y[7] = Iahp2::m0;
        y[8] = Ican_TC::m0;
        }

double RE::V0 = -70, RE::Cai0 = 5e-5;

void RE::calc(double x, double *y, double *f){
    IT_RE::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    INaK::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], x);
    Ca_RE::calc(y[1], f[1], iT, x);
    Iahp2::calc(y[7], f[7], y[1], y[0]);
    Ican_TC::calc(y[8], f[8], y[1], y[0]);

    i_Total = S_RE*1e3*(-G_l*(y[0] - E_l) -G_kl*(y[0] - EK) -iNa -iK -iT -iAHP -iCAN) + I_Stim;
    f[0] = i_Total/(1*S_RE*1e3);
}

void RE::setStim(double I_Stim_Param){
  I_Stim = I_Stim_Param; 
}


