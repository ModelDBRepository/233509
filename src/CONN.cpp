#include "CONN.h"
#include <math.h>
#include <stdio.h>


CONN::CONN() {
  i = 0;
  N_Pre_TC1 = 0;
  N_Pre_TC2 = 0;
  N_Pre_IN  = 0;
  N_Pre_RE  = 0;
  N_Pre_RE_GAP  = 0;
}


void CONN::INIT() {

    F_TC1_TC2=0;
    F_TC2_TC1=0;
    F_RE_RE = 0;

 for(i=0; i<N_TC1_X; i++)
   for(j=0;j<N_TC1_Y; j++)
      C_TC1_TC1[i][j]=0;

 for(i=0; i<N_RE_X; i++)
   for(j=0;j<N_RE_Y; j++)
      C_RE_RE[i][j]=0;

   for (i=0; i<N_TC1; i++) {
	 Pre_TC1_X[i] = 0;
	 Pre_TC1_Y[i] = 0;
	 Pre_TC1_D[i] = 0.0;
   }

   for (i=0; i<N_TC2; i++) {
	 Pre_TC2_X[i] = 0;
	 Pre_TC2_Y[i] = 0;
	 Pre_TC2_D[i] = 0.0;
   }

   for (i=0; i<N_IN; i++)  {
	 Pre_IN_X[i] = 0;
	 Pre_IN_Y[i] = 0;
   }

   for (i=0; i<N_RE; i++)  {
	 Pre_RE_X[i] = 0;
	 Pre_RE_Y[i] = 0;
	 Pre_RE_X_GAP[i] = 0;
	 Pre_RE_Y_GAP[i] = 0;
	 Pre_RE_D_GAP[i] = 0.0;
   }


}




