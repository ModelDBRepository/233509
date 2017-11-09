

class Iahp2 {
  static double E_K, Cels;
  double alpha, beta, m_inf, tau_m, Tad;
public:
  double iAHP, G_AHP, m0, Shift_m;
  Iahp2(double v, double ca);

  void init(double v, double ca);
  void calc(double m, double &fm, double ca, double v);
};
