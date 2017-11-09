
#include "../Constants.h"

class CONN {
 // protected:
private:
   int i, j;

public:

   int N_Pre_TC1;
   int N_Pre_TC2;
   int N_Pre_IN;
   int N_Pre_RE;
   int N_Pre_RE_GAP;

   int Pre_TC1_X[N_TC1];
   int Pre_TC1_Y[N_TC1];
   double Pre_TC1_D[N_TC1];

   int Pre_TC2_X[N_TC2];
   int Pre_TC2_Y[N_TC2];
   double Pre_TC2_D[N_TC2];

   int Pre_IN_X[N_IN];
   int Pre_IN_Y[N_IN];

   int Pre_RE_X[N_RE];
   int Pre_RE_Y[N_RE];

   int Pre_RE_X_GAP[N_RE];
   int Pre_RE_Y_GAP[N_RE];
   double Pre_RE_D_GAP[N_RE];

   int C_TC1_TC1[N_TC1_X][N_TC1_Y];
   int C_RE_RE[N_RE_X][N_RE_Y];

   int F_TC1_TC2;
   int F_TC2_TC1;
   int F_RE_RE;

   CONN();
   void INIT();

};

