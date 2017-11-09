
class IT_RE {
  static double Ca_0, Cels;
  double m_inf, tau_m, h_inf, tau_h, ratio, eca, Phi_m, Phi_h, eca0;
  double ECA;

  
public:
  double iT, m0, h0, Qm, Qh;
  double G_Ca;
  double Shift;
  IT_RE(double v); 

  void init(double v);
  void calc(double m, double h, double &fm, double &fh, double v, double cai, double x);
};
