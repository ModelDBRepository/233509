

class Ih_TC {
  static double E_H, Cels;
  double m_inf, tau_m, Tad;
public:
  double iH, m0, G_H, Shift_m;
  Ih_TC(double v);

  void init(double v);
  void calc(double m, double &fm, double v);
};
