
#include "AMPA_ext.h"
#include <math.h>
#include <stdlib.h>

AMPA_ext::AMPA_ext(){
  Tau = 5;
  Deadtime = 0.2;   // 2.0 !!!
  g_spike = 0; 
  TR = 25;
  w = 0.1;  // 0.1
  q = 0;
  g = 0;
  spike = 0;
  random_number = 0;
  lastrelease = 0;
}



void AMPA_ext::init(unsigned int seek) {srand(seek);}
void AMPA_ext::calc(double g_Extern_ampa, double x){
  
  q = x - lastrelease;
  if (q > Deadtime) {
    if (q > TR) {                  
      lastrelease = x;
      random_number = 1.0*rand()/(RAND_MAX + 1.0);
      if(random_number < 0.000001) random_number = 0.000001;
      TR = -(log(random_number))/(w);
      spike = 1;
    }
  }
  if (spike == 1){
    g = g + g_Extern_ampa;
    g_spike = g;
    spike = 0;    
  } else {                              
    g = g_spike * exptable(-(x - (lastrelease))/Tau);
  }
  
  
  
}

double AMPA_ext::exptable(double z) {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return(0);
}
