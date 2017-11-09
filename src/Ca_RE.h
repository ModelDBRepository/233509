
class Ca_RE {
  static double Ca_inf, K_T, K_d;
  double drive, drive0;                                 
public:
  double Taur, D;
  Ca_RE();
  void calc(double cai, double &fcai, double iT, double x);
};
