
#include "AMPA.h"
#include <math.h>
#include <stdio.h>
#include "../Constants.h"

AMPA::AMPA() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 0;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
	E = 1;
	Use = F_USE_SPIKE;
	Tr = TAU_STD;
}


double AMPA::E_AMPA = 0, AMPA::Cdur = 0.3, AMPA::Cmax = 0.5, AMPA::Deadtime = 1;
double AMPA::Cdel = 2;
double AMPA::Alpha = 0.94, AMPA::Beta = 0.18;  // alpha = 0.94; beta = 0.18


void AMPA::calc(double g_AMPA, double x, double y_post,
		   double prespiketime, int FD) {
  
    q = ((x - lastrelease) - Cdur);

    if(q > Deadtime) {
           if( (x - prespiketime)<0.002 ) {
              if( (x - lastspike) > (Cdel + Cdur) ){
                 lastspike = x;
                 s = 1; } }  //the flag that spike was but wasn't utilized yet

           if((s == 1) && ((x - lastspike) > Cdel)) {
              s = 0;         //spike was utilized
              C = Cmax;
              R0 = R;

              if (FD == 1) {
              E = 1 - (1 - E*(1-Use)) * exptable(-q/Tr);
              } else {
              E = 1;
              }

              lastrelease = x;  }
    } else if (q < 0) {
    } else if (C == Cmax) {
            R1 = R;
            C = 0.;
    }

    if (C > 0) {
       R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
    } else {
       R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
    }
   I = g_AMPA * E * R * (y_post - E_AMPA);

}



double AMPA::exptable(double z) {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return(0);
}
