class NMDA {

  static double E_NMDA;
  static double Cdur, Cmax, Deadtime, Cdel;
  static double Alpha, Beta;
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike, Use, Tr;
  double q, Rinf, Rtau, fn;
  double exptable(double z);

public:
  double I, E;
  NMDA();      // Constructor
  void calc(double g_NMDA, double x, double y_post,
	    double prespiketime, int FD);
};



