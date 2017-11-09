

class IKs_TC {
  static double E_K, Cels;     
  double m_inf, tau_m, h_inf, tau_h, Tad;                                 
public:
  double iKs, m0, h0, G_Ks;
  IKs_TC(double v);
  void calc(double m, double h, double &fm, double &fh, double v);
};
