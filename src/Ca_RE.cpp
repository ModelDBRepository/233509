#include "Ca_RE.h"
#include <math.h>

Ca_RE::Ca_RE() {Taur = 10; D = 1.; //0.1;
         drive0 = 10.0/(2.*96489.); }

double Ca_RE::Ca_inf = 5e-5; // 2.4e-4 mM || 5e-5
double Ca_RE::K_T = 0.0001, Ca_RE::K_d = 0.0001;

void Ca_RE::calc(double cai, double &fcai, double iT, double x) {
  drive = -drive0 * (iT) / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
}
