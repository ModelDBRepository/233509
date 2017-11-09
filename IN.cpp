
#include "Constants.h"
#include <stdio.h>
#include "src/IT_TC.h"
#include "src/Ih_TC.h"
#include "src/INaK.h"
#include "src/Iahp2.h"
#include "src/Ican_TC.h"
//#include "src/ICaL_TC.h"
//#include "src/IT0_TC.h"
//#include "src/IA_TC.h"
//#include "src/INaP_TC.h"
//#include "src/IKs_TC.h"

#include "IN.h"

IN::IN():Ih_TC(V0), IT_TC(V0), INaK(V0), Ca_IN(), Iahp2(V0, Cai0), Ican_TC(V0) {
        G_nal = 0.01;
	    G_kl = 0.012;
        E_l = -70;
        I_Stim = 0;

        INaK::G_Na = 90;
        INaK::G_K = 10;
        S_IN = 1.7e-4;  // cm^2
        }
        
void IN::init(double *y) {
        INaK::init(V0);
        IT_TC::init(V0);
        Ih_TC::init(V0);
        Iahp2::init(V0, Cai0);
        Ican_TC::init(V0);

        y[0] = V0;
	    y[1] = Cai0;
        y[2] = IT_TC::m0;
        y[3] = IT_TC::h0;
        y[4] = Ih_TC::m0;
        y[5] = INaK::m0;
        y[6] = INaK::h0;
        y[7] = INaK::n0;
        y[8] = Iahp2::m0;
        y[9] = Ican_TC::m0;

       /*
        y[8] = IA_TC::m0;
        y[9] = IA_TC::h0;
        y[10] = INaP_TC::m0;
        y[11] = IKs_TC::m0;
        y[12] = IKs_TC::h0;
        y[14] = IT0_TC::m0;
        y[15] = IT0_TC::h0;
        y[17] = ICaL_TC::m0;
        y[18] = ICaL_TC::h0;
        */
        }


double IN::Cai0 = 5e-5;  // in mM or 50 nM
// double TC::G_l = 0.01;    //0.05; //0.01;/// mS/cm^2
double IN::V0 = -70;

void IN::calc(double x, double *y, double *f){
    IT_TC::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    Ih_TC::calc(y[4], f[4], y[0]);
    INaK::calc(y[5], y[6], y[7], f[5], f[6], f[7], y[0], x);
    Iahp2::calc(y[8], f[8], y[1], y[0]);
    Ican_TC::calc(y[9], f[9], y[1], y[0]);
    Ca_IN::calc(y[1], f[1], iT, x);

    //  IT0_TC::calc(y[14], y[15], f[14], f[15], y[0], y[1], x);
    //  ICaL_TC::calc(y[17], y[18], f[17], f[18], y[0], y[1], x);
    //  IA_TC::calc(y[8], y[9], f[8], f[9], y[0], x);
    //  INaP_TC::calc(y[10], f[10], y[0]);
    //  IKs_TC::calc(y[11], y[12], f[11],f[12], y[0]);

    i_Total = S_IN*1e3*(-G_l*(y[0] - E_l) - iT - iH - iNa - iK
    		 -iAHP -iCAN - G_kl*(y[0] - EK)) + I_Stim;
    f[0] = i_Total/(S_IN*1e3);
}

void IN::setStim(double I_Stim_Param){
  I_Stim = I_Stim_Param; 
}



