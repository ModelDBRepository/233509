

class AMPA {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime,Cdel;
  static double Alpha, Beta;
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike, Use, Tr;
  double q, Rinf, Rtau;

  double exptable(double z);

 public:
   double I, E;
   AMPA();
   void calc(double, double, double, double, int);
};
