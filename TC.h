

#include "src/INaK.h"
#include "src/Ih_TC.h"
#include "src/IT0_TC.h"
#include "src/IT_TC.h"
#include "src/ICaL_TC.h"
#include "src/Ca_TC.h"
#include "src/Iahp2.h"
#include "src/Ican_TC.h"
//#include "src/IA_TC.h"
//#include "src/INaP_TC.h"
//#include "src/IKs_TC.h"


class TC: public INaK, public Ih_TC, public IT0_TC, public IT_TC, public ICaL_TC,
          public Ca, public Iahp2, public Ican_TC {
  static double Cai0, V0;
public:
  double G_l, G_nal, G_kl, S_TC, E_l, I_Stim, i_Total;

  TC();     
  void init(double *y);
  void calc(double x, double *y, double *f); 
    // Set stimulation
  void setStim(double);
};  
