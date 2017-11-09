

class Ican_TC {
  static double E_CAN, Cels;
  double m_inf, tau_m, Kd, fCa, Tad;
public:
  double iCAN, m0, G_CAN, Shift_m;
  Ican_TC(double v);

  void init(double v);
  void calc(double m, double &fm, double ca, double v);
};
