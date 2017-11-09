#include "Ca_TC.h"
#include <math.h>

Ca::Ca() {Taur = 5; D = 1.; //0.1;
         drive0 = 10.0/(2.*96489.); }   // UNIT of current is: uA/cm^2

double Ca::Ca_inf = 5e-5; // 2.4e-4 mM  || 5e-5
double Ca::K_T = 0.0001, Ca::K_d = 0.0001;

void Ca::calc(double cai, double &fcai, double iT, double iT0, double iCaL, double x) {
  drive = -drive0 * (iT+iT0+iCaL) / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
}
