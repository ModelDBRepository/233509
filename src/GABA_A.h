
//#include "../Constants.h"

class GABA_A {

  static double Cdur, Cmax, Deadtime, Cdel;;
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike, Use, Tr;
  double q, Rinf, Rtau;
  double exptable(double z);

public:
  double I, E_GABA, Alpha, Beta;
  double E;
  GABA_A();
  void calc(double, double, double, double, int);
};
