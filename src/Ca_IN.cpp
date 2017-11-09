#include "Ca_IN.h"
#include <math.h>

Ca_IN::Ca_IN() {Taur = 5; D = 1.; //0.1;
         drive0 = 10.0/(2.*96489.); }   // UNIT of current is: uA/cm^2

double Ca_IN::Ca_inf = 5e-5; // 2.4e-4 mM
double Ca_IN::K_T = 0.0001, Ca_IN::K_d = 0.0001;

void Ca_IN::calc(double cai, double &fcai, double iT, double x) {
  drive = -drive0 * (iT) / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
}
