

class IT_TC {
  static double Ca_0, Cels, Qm, Qh;
  double m_inf, tau_m, h_inf, tau_h, Phi_h, Phi_m, ratio, eca, eca0; 
public:
  double iT, m0, h0, Shift_m, Shift_h;
  double G_Ca;
  double K_Tau_h;

  IT_TC(double v); 
  void init(double v);
  void calc(double m, double h, double &fm, double &fh,  
            double v, double cai, double x);
};
