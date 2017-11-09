
#include "src/IT_RE.h"
#include "src/Ca_RE.h"



class RE: public IT_RE, public INaK, public Ca_RE, public Iahp2, public Ican_TC {
  static double Cai0, V0;
public:
  double G_kl, G_l, E_l, S_RE, I_Stim, i_Total;

  RE(); 
  void init(double *y);
  void calc(double x, double *y, double *f); 
  // Set stimulation
  void setStim(double);
};  
