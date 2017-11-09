

class Iahp_TC {
  static double E_K, Cels;
  double alpha, beta, fV, fCa, Tad;
public:
  double iAHP, G_AHP, m0, Shift_m;
  Iahp_TC(double v, double ca);
  void calc(double m, double &fm, double ca, double v);
};
