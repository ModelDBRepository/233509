
class Ca {
  static double Ca_inf, K_T, K_d;
  double drive, drive0;                                 
public:
  double Taur, D;
  Ca();
  void calc(double cai, double &fcai, double iT, double iT0, double iCaL, double x);
};
