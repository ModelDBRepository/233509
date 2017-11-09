

class INaK {
  static double Cels;
  double Alpha1, Beta1, Alpha2, Beta2, Alpha3, Beta3, v2, v2K, Phi;
  double tau_m, m_inf, tau_h, h_inf, tau_n, n_inf;
 // double ENa;
  
public:
  static double E_Na, E_K;
  double iK, iNa, m0, h0, n0;
  double G_Na, G_K, Vtr, VtrK, k_h, k_n;
  
  INaK(double v);
  void init(double v);
  void calc(double m, double h, double n, double &fm, double &fh, double &fn, double v, double x);
};
