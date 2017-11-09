

class AMPA_ext {
  double Deadtime; 

  int spike;
  double Tau;
  // double w;
  double TR;
  double q;
  double g_spike;
  double exptable(double z);
  double random_number;
  
 public:
  double w;
  double g;
  double lastrelease;

  AMPA_ext();
  void init(unsigned int seek);
  void calc(double g_Extern_ampa, double x);
};
