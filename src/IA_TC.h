

class IA_TC {
  static double E_K, Cels;     
  double m_inf, tau_m, h_inf, tau_h, Tad;                                 
public:
  double iA, m0, h0, G_A;
  IA_TC(double v); 
  void calc(double m, double h, double &fm, double &fh, double v, double x);
};
