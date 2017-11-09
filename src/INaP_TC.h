

class INaP_TC {
  static double E_NaP, Cels;
  double alpha, beta, m_inf, tau_m, Tad;
public:
  double iNaP, m0, G_NaP;
  INaP_TC(double v);
  void calc(double m, double &fm, double v);
};
