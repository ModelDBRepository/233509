
#include "Stimulation.h"
#include "../Constants.h"
#include <math.h>
#include <stdlib.h>

STIM::STIM(){
	k = 0;
	I = 0;
	ON = 1;

	FLAG_VAR_F = 1;

	FLAG_F = 1;
    FLAG_I = 1;

	Stim_T0  = STIM_START_TIME;
	Stim_Dur = STIM_DURATION;
	Stim_Amp = STIM_AMP;                 // 50 pA

	F0 = STIM_F0;                        // Initial frequency
	T_PULSE_ON = PULSE_ON_DURATION;      // Pulse ON duration
	T_Step     = STIM_STEP_DURATION;      // Duration of each step frequency stimulation
	STEP_F     = STIM_F_INCREMENT;

}


void STIM::init() {
	F = F0;
	Pulse_T0  = Stim_T0;
	Pulse_T00 = Stim_T0;
	Step_T0   = Stim_T0;
}

void STIM::calc(double x){

  if (x>=Stim_T0 && x<=Stim_T0+Stim_Dur) {

     if (FLAG_VAR_F == 1) {
	   if (x>=Step_T0 + T_Step) {
          Step_T0 = Step_T0 + T_Step;
     	  Pulse_T0  = Step_T0;
    	  Pulse_T00 = Step_T0;
      	  ON = 1;

      	  if (FLAG_F == 1) {
      	     F = F + STEP_F;
      	  } else {
      	     F = F - STEP_F;
      	  }

      	  if (F > STIM_MAX_F) {
             F = F-2*STEP_F;
          	 FLAG_F = 0;
      	   }
         }
      }


   	 if (F <= 0) {
          FLAG_I = 0;
     } else {
    	  FLAG_I = 1;

    	  T_Pulse = 1/F*1000;

    	  if (x>=Pulse_T0) {
    	  	  if (ON==0) {
    	   		Pulse_T00 = Pulse_T00 + T_Pulse;
    	        ON = 1;
    	   	  }
    	    }

    	  if (x>=Pulse_T00+T_PULSE_ON) {
    	    if (ON==1) {
    	   		Pulse_T0 = Pulse_T0 + T_Pulse;
    	        ON = 0;
    	  	  }
    	    }
     }


       if (ON==1) {
	      I = Stim_Amp;
       } else {
         I = 0;
       }

      if (FLAG_I==0) {
    	  I = 0;
      }

     }  else  {
        I = 0;
    }

}


  
  



