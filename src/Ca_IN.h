
class Ca_IN {
  static double Ca_inf, K_T, K_d;
  double drive, drive0;                                 
public:
  double Taur, D;
  Ca_IN();
  void calc(double cai, double &fcai, double iT, double x);
};
