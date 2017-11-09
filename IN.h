

#include "src/Ca_IN.h"
//#include "src/IT_TC.h"
//#include "src/Ih_TC.h"
//#include "src/INaK.h"
//#include "src/IA_TC.h"
//#include "src/INaP_TC.h"
//#include "src/IKs_TC.h"
//#include "src/Iahp_TC.h"
//#include "src/IT0_TC.h"
//#include "src/Ican_TC.h"
//#include "src/ICaL_TC.h"


class IN: public IT_TC, public Ih_TC, public INaK, public Ca_IN,
     public Iahp2, public Ican_TC {
  static double Cai0, V0;
public:
  double G_l, G_nal, G_kl, S_IN, E_l, I_Stim, i_Total;

  IN();
  void init(double *y);
  void calc(double x, double *y, double *f); 
    // Set stimulation
  void setStim(double);
};  
