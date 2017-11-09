
#include <stdio.h>
#include "Constants.h"
#include "TC.h"

TC::TC():INaK(V0), Ih_TC(V0), IT0_TC(V0), IT_TC(V0), ICaL_TC(V0), Ca(),
      Iahp2(V0, Cai0), Ican_TC(V0) {
        G_nal = 0.0;
	    G_kl  = 0.01;
        E_l   = -70;
        I_Stim = 0;

        INaK::G_Na = 90;
        INaK::G_K = 10;
        S_TC = 2.9e-4;  // 2.9e-4cm^2 !!!
        }
        
void TC::init(double *y) {
       INaK::init(V0);
       IT0_TC::init(V0);
       IT_TC::init(V0);
       Ih_TC::init(V0);
       Iahp2::init(V0, Cai0);
       Ican_TC::init(V0);
       ICaL_TC::init(V0);

        y[0] = V0;
	    y[1] = Cai0;
        y[2] = IT_TC::m0;
        y[3] = IT_TC::h0;
        y[4] = Ih_TC::m0;
        y[5] = INaK::m0;
        y[6] = INaK::h0;
        y[7] = INaK::n0;
        y[8] = Iahp2::m0;
        y[9] = IT0_TC::m0;
        y[10] = IT0_TC::h0;
        y[11] = Ican_TC::m0;
        y[12] = ICaL_TC::m0;
        y[13] = ICaL_TC::h0;
        // y[8] = IA_TC::m0;
        // y[9] = IA_TC::h0;
        // y[10] = INaP_TC::m0;
        // y[11] = IKs_TC::m0;
        // y[12] = IKs_TC::h0;
        }


double TC::Cai0 = 5e-5;  // in mM or 50 nM
// double TC::G_l = 0.01;    //0.05; //0.01;/// mS/cm^2
double TC::V0 = -65;    // -70 !!!

void TC::calc(double x, double *y, double *f){
    IT_TC::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    Ih_TC::calc(y[4], f[4], y[0]);
    INaK::calc(y[5], y[6], y[7], f[5], f[6], f[7], y[0], x);
    Iahp2::calc(y[8], f[8], y[1], y[0]);
    IT0_TC::calc(y[9], y[10], f[9], f[10], y[0], y[1], x);
    Ican_TC::calc(y[11], f[11], y[1], y[0]);
    ICaL_TC::calc(y[12], y[13], f[12], f[13], y[0], y[1], x);
    Ca::calc(y[1], f[1], iT, iT0, iCaL, x);

    //   IA_TC::calc(y[8], y[9], f[8], f[9], y[0], x);
    //   INaP_TC::calc(y[10], f[10], y[0]);
    //   IKs_TC::calc(y[11], y[12], f[11],f[12], y[0]);

    i_Total = S_TC*1e3*(-G_l*(y[0] - E_l) -G_kl*(y[0] - EK) -iNa -iK -iH -iT0 -iT -iCaL
    		 -iAHP -iCAN)  + I_Stim;
    f[0] = i_Total/(S_TC*1e3);
}

void TC::setStim(double I_Stim_Param){
  I_Stim = I_Stim_Param; 
}



