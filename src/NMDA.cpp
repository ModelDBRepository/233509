
#include "NMDA.h"
#include <math.h>
#include <stdio.h>
#include "../Constants.h"

NMDA::NMDA() {
	R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
	lastspike = -100;
	s = 0;
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
	Rtau = 1 / ((Alpha * Cmax) + Beta);
	E = 1;
	Use = F_USE_SPIKE;
	Tr = TAU_STD;
    I = 0;
}

double NMDA::E_NMDA = 0, NMDA::Cdur = 0.3, NMDA::Cmax = 0.5;
double NMDA::Deadtime = 1;
double NMDA::Cdel = 2;
double NMDA::Alpha = 1, NMDA::Beta = 0.0067;   // alpha = 1(0.25); beta = 0.0067 (0.025)


void NMDA::calc(double g_NMDA, double x, double y_post,
		double prespiketime, int FD) {
  
    q = ((x - lastrelease) - Cdur);

    if(q > Deadtime) {
           if( (x - prespiketime)<0.002 ) {
              if( (x - lastspike) > (Cdel + Cdur) ){
                 lastspike = x;
                 s = 1; } }  //the flag that spike was but wasn't utilized

           if((s == 1) && ((x - lastspike) > Cdel)) {
              s = 0;         //spike was utilized
              C = Cmax;
              R0 = R;

              if (FD==1) {
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
   fn = 1/(1+exp(-(y_post - (-25))/12.5));
   I = g_NMDA * R * fn * E * (y_post - E_NMDA);
}



double NMDA::exptable(double z) {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return(0);
}


