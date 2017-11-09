
// ************************************************************************************************************************
// An unified thalamic network model to generate delta, spindle, alpha and gamma oscillations
// Developed by Guoshi Li, University of North Carolina at Chapel Hill (guoshi_li@med.unc.edu)

// Simulation results are presented in the associated paper:
// Guoshi Li, Craig S Henriquez and Flavio Fröhlich (2017) Unified Thalamic Model Generates Multiple Distinct Oscillations
// with State-dependent Entrainment by Stimulation
// PLOS Computational Biology (In Press).
// ****************************************************************************************************************************


// Standard
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Own
#include "Constants.h"
#include "TC.h"
#include "IN.h"
#include "RE.h"
#include "src/CONN.h"
#include "src/AMPA_ext.h"
#include "src/AMPA.h"
#include "src/NMDA.h"
#include "src/GABA_A.h"
#include "src/Stimulation.h"


//#define PI 3.141592654

/*******************************************************************
// Declarations External Functions
*******************************************************************/

void rk(unsigned, void (unsigned, double, double*, double*), double, double, 
                                double*, double*, double*, double*);
void fun(unsigned, double, double*, double*);

double gaussrand();


/*******************************************************************
// Declarations External Variables (File pointers)
*******************************************************************/
FILE *finput, *f_re, *f_tc1, *f_tc2, *f_in, *f_input, *f_spike, *f_ca;
FILE *f_I, *f_SI, *f_Iacs;
FILE *f_E, *f_E_TC_RE, *f_E_RE_TC, *f_W_HTC_IN;
FILE *f_m_T, *f_h_T, *f_m_H, *f_h_KS, *f_m_AHP, *f_m_A, *f_h_A;
FILE *f_connect;


/*******************************************************************
// Indices
*******************************************************************/
int i_tc1[N_TC1_X][N_TC1_Y];  // high-threshold bursting TC cell (HTC)
int i_tc2[N_TC2_X][N_TC2_Y];  // regular low-threshold   TC cell (RTC)

int i_in[N_IN_X][N_IN_Y];     // LGN interneurons
int i_re[N_RE_X][N_RE_Y];     // RE cells

/*******************************************************************
// Create cells and synapses
*******************************************************************/
TC    tc1_cell[N_TC1_X][N_TC1_Y];
TC    tc2_cell[N_TC2_X][N_TC2_Y];

IN    in_cell[N_IN_X][N_IN_Y];
RE    re_cell[N_RE_X][N_RE_Y];

STIM  Stim, Stim_Spin;

AMPA_ext ampa_ext_tc1[N_TC1_X][N_TC1_Y][N_INPUT_TC1];
AMPA_ext ampa_ext_tc2[N_TC2_X][N_TC2_Y][N_INPUT_TC2];
AMPA_ext ampa_ext_in[N_IN_X][N_IN_Y];
AMPA_ext ampa_ext_re[N_RE_X][N_RE_Y];

AMPA     ampa_TC1_IN[N_IN_X][N_IN_Y][N_TC1];
NMDA     nmda_TC1_IN[N_IN_X][N_IN_Y][N_TC1];
AMPA     ampa_TC1_RE[N_RE_X][N_RE_Y][N_TC1];
NMDA     nmda_TC1_RE[N_RE_X][N_RE_Y][N_TC1];

AMPA     ampa_TC2_RE[N_RE_X][N_RE_Y][N_TC2];
NMDA     nmda_TC2_RE[N_RE_X][N_RE_Y][N_TC2];
AMPA     ampa_TC2_TC2[N_TC2_X][N_TC2_Y][N_TC2];
NMDA     nmda_TC2_TC2[N_TC2_X][N_TC2_Y][N_TC2];

GABA_A   gaba_IN_TC2[N_TC2_X][N_TC2_Y][N_IN];

GABA_A   gaba_RE_TC1[N_TC1_X][N_TC1_Y][N_RE];
GABA_A   gaba_RE_TC2[N_TC2_X][N_TC2_Y][N_RE];
GABA_A   gaba_RE_IN[N_IN_X][N_IN_Y][N_RE];
GABA_A   gaba_RE_RE[N_RE_X][N_RE_Y][N_RE];

CONN conn_tc1[N_TC1_X][N_TC1_Y];
CONN conn_tc2[N_TC2_X][N_TC2_Y];
CONN conn_in[N_IN_X][N_IN_Y];
CONN conn_re[N_RE_X][N_RE_Y];


/*******************************************************************
// Spike times
*******************************************************************/
double TC1_SPT[N_TC1_X][N_TC1_Y];
double TC2_SPT[N_TC2_X][N_TC2_Y];
double IN_SPT[N_IN_X][N_IN_Y];
double RE_SPT[N_RE_X][N_RE_Y];

/*******************************************************************
//  External current inputs
*******************************************************************/
double input_tc1[N_TC1_X][N_TC1_Y];
double input_tc2[N_TC2_X][N_TC2_Y];
double input_in[N_IN_X][N_IN_Y];
double input_re[N_RE_X][N_RE_Y];

/*******************************************************************
//Short-term dynamics synapses
*******************************************************************/

double TC2_E[N_TC2_X][N_TC2_Y];
double RE_E[N_RE_X][N_RE_Y];

// Afferent input conductance
double g_Extern_TC1, g_Extern_TC2, g_Extern_IN, g_Extern_RE;

/*******************************************************************
// Main function
*******************************************************************/

int main(int argc,char **argv)
{
  
  srand(RAND_SEED);

  double y_ini[N_EQ], f_ini[N_EQ];
  double y1[N_EQ], y2[N_EQ];
  

  double t = 0;        // time start
  double tau = 0;      //
  double h = 0.02;     // Default: 0.02; Step width solver

  double gL_TC1, gKL_TC1, gL_TC2, gKL_TC2;
  double gL_IN, gKL_IN, gH_IN, gCAN_IN, gL_RE, gKL_RE;
  double d = 0;

  int    i,j, k, i1,j1, ii, jj;
  double mi,  mj;  // indices

  int    PRE, POST;

  int FLAG_SPINDLE_INPUT = 0;

  char filepath[] = "data/";
  char filename[50];


 // Print the beginning time of simulation
  printf("\n \n \nProgram started on ");
  time_t t_start=time(0);
  
  char * mytime= ctime(&t_start);
  printf(mytime);


 // Set afferent input conductance for each oscillation state
  if (OSC == 1) {             // Delta
	  g_Extern_TC1 = 0.0001;  // uS
	  g_Extern_TC2 = 0.0001;
	  g_Extern_IN  = 0.0001;
	  g_Extern_RE  = 0.0001;

  } else if (OSC == 2) {      // Spindle
	  g_Extern_TC1 = 0.0003;  // uS
	  g_Extern_TC2 = 0.0003;
	  g_Extern_IN  = 0.0003;
	  g_Extern_RE  = 0.0003;

      FLAG_SPINDLE_INPUT   = 1;


  } else if (OSC == 3){       // Alpha
      g_Extern_TC1 = 0.0015;  // uS
	  g_Extern_TC2 = 0.0015;
	  g_Extern_IN  = 0.0015;
	  g_Extern_RE  = 0.0015;

  }

  else {                      // Gamma
	  g_Extern_TC1 = 0.017;   // uS
	  g_Extern_TC2 = 0.017;
	  g_Extern_IN  = 0.0015;
	  g_Extern_RE  = 0.0015;

  }

  // Parameters of spindle-triggering input
  Stim_Spin.FLAG_VAR_F = 0;
  Stim_Spin.Stim_T0    = SPIN_T0;  // spindle input start time
  Stim_Spin.Stim_Amp   = SPIN_AMP; // spindle input amplitude
  Stim_Spin.F0         = SPIN_F;   // spindle input frequency
  Stim_Spin.T_PULSE_ON = SPIN_TT;  // spindle input duration


 // Constant current Injections (not used)
  for(i=0; i<N_TC1_X; i++)
	for(j=0; j<N_TC1_Y; j++)
		input_tc1[i][j] = 0.0;

  for(i=0; i<N_TC2_X; i++)
	for(j=0; j<N_TC2_Y; j++)
		input_tc2[i][j] = 0.0;

  for(i=0; i<N_IN_X; i++)
	for(j=0; j<N_IN_Y; j++)
		input_in[i][j] = 0.0;

  for(i=0; i< N_RE_X; i++)
	for(j=0; j< N_RE_Y; j++)
		input_re[i][j] = 0.0;


// Set parameters for each cell type
// High-threshold TC cells
  for(i=0; i<N_TC1_X; i++)
	for(j=0;j<N_TC1_Y; j++) {

	  gL_TC1  = GL1_TC1 + (GL2_TC1 - GL1_TC1)*(rand()/(RAND_MAX + 1.0));

	  if (OSC == 1) {
           gKL_TC1 = GKL_TC1_S1;
      } else if (OSC == 2) {
    	   gKL_TC1 = GKL_TC1_S2;
      } else  {
    	   gKL_TC1 = GKL_TC1_S3;
      }

      tc1_cell[i][j].G_l   = gL_TC1;
      tc1_cell[i][j].G_kl  = gKL_TC1;

      tc1_cell[i][j].G_Na  = 90;
      tc1_cell[i][j].G_K   = 10;

      tc1_cell[i][j].G_Ca0 = 2.1;
      tc1_cell[i][j].G_Ca  = 3.0;
      tc1_cell[i][j].G_CaL = 0.5;

      tc1_cell[i][j].G_AHP = 0.3;
      tc1_cell[i][j].G_CAN = 0.5;
      tc1_cell[i][j].G_H   = 0.01;

      tc1_cell[i][j].E_l  = -70;
      tc1_cell[i][j].Taur =  10;
      tc1_cell[i][j].D    = 0.5;

      tc1_cell[i][j].INaK::Vtr  = -30;
      tc1_cell[i][j].INaK::VtrK = -30;

      tc1_cell[i][j].INaK::k_h = 1;
      tc1_cell[i][j].INaK::k_n = 4;

      tc1_cell[i][j].IT0_TC::Shift_m = -3;
      tc1_cell[i][j].IT0_TC::Shift_h = -3;

      tc1_cell[i][j].IT_TC::Shift_m = 25;
      tc1_cell[i][j].IT_TC::Shift_h = 25;

      tc1_cell[i][j].IT_TC::K_Tau_h = 1.0;

      tc1_cell[i][j].Ih_TC::Shift_m = 0;
	}


// Low-threshold TC cells
  for(i=0; i<N_TC2_X; i++)
	for(j=0;j<N_TC2_Y; j++) {

	  gL_TC2  = GL1_TC2 + (GL2_TC2 - GL1_TC2)*(rand()/(RAND_MAX + 1.0));

	  if (OSC==1) {
	         gKL_TC2 = GKL_TC2_S1;
	   } else if (OSC==2) {
	         gKL_TC2 = GKL_TC2_S2;
	   } else  {
	         gKL_TC2 = GKL_TC2_S3;
	   }

	  tc2_cell[i][j].G_l   = gL_TC2;
      tc2_cell[i][j].G_nal = 0.0 ;
      tc2_cell[i][j].G_kl  = gKL_TC2;

      tc2_cell[i][j].G_Na  = 90;
      tc2_cell[i][j].G_K   = 10;

      tc2_cell[i][j].G_Ca0 = 2.1;
      tc2_cell[i][j].G_Ca  = 0.6;
      tc2_cell[i][j].G_CaL = 0.3;

      tc2_cell[i][j].G_AHP = 0.1;
      tc2_cell[i][j].G_CAN = 0.6;
      tc2_cell[i][j].G_H   = 0.01;

      tc2_cell[i][j].E_l  = -70;
      tc2_cell[i][j].Taur = 10;
      tc2_cell[i][j].D    = 0.5;

      tc2_cell[i][j].INaK::Vtr  = -40;
      tc2_cell[i][j].INaK::VtrK = -40;

      tc2_cell[i][j].INaK::k_h = 1;
      tc2_cell[i][j].INaK::k_n = 4;

      tc2_cell[i][j].IT0_TC::Shift_m = -3;
      tc2_cell[i][j].IT0_TC::Shift_h = -3;

      tc2_cell[i][j].IT_TC::Shift_m = 25;
      tc2_cell[i][j].IT_TC::Shift_h = 25;

      tc2_cell[i][j].Ih_TC::Shift_m = 0;
	}

//=================================================================
//  IN Cells
  for(i=0; i<N_IN_X; i++)
	for(j=0;j<N_IN_Y; j++) {

	  gL_IN  = GL1_IN + (GL2_IN - GL1_IN)*(rand()/(RAND_MAX + 1.0));

	  if (OSC==1) {
	       gKL_IN = GKL_IN_S1;

	  } else if (OSC==2) {
	   	   gKL_IN = GKL_IN_S2;

	  } else  {
	   	   gKL_IN  = GKL_IN_S3;

	  }

      in_cell[i][j].G_l   = gL_IN;
      in_cell[i][j].G_kl  = gKL_IN;

      in_cell[i][j].G_H   = 0.05;
      in_cell[i][j].G_CAN = 0.1;

      in_cell[i][j].G_Na  = 90;
      in_cell[i][j].G_K   = 10;
      in_cell[i][j].G_Ca  = 2.5;

      in_cell[i][j].G_AHP = 0.2;

      in_cell[i][j].E_l  = -60;
      in_cell[i][j].Taur =  10;
      in_cell[i][j].D    = 0.5;

      in_cell[i][j].INaK::Vtr  = -30;
      in_cell[i][j].INaK::VtrK = -30;

      in_cell[i][j].INaK::k_h = 1;
      in_cell[i][j].INaK::k_n = 4;

      in_cell[i][j].IT_TC::Shift_m = 25;
      in_cell[i][j].IT_TC::Shift_h = 25;

      in_cell[i][j].Ih_TC::Shift_m = 0;
  }



//=====================================================================
// RE Cells
  for(i=0; i<N_RE_X; i++)
	for(j=0;j<N_RE_Y; j++) {

	  if (OSC==1) {
		   gKL_RE = GKL_RE_S1;
	  } else if (OSC==2) {
	   	   gKL_RE = GKL_RE_S2;
	  } else  {
	   	   gKL_RE = GKL_RE_S3;
	  }

	  gL_RE  = GL1_RE + (GL2_RE - GL1_RE)*(rand()/(RAND_MAX + 1.0));

	  re_cell[i][j].G_l  = gL_RE;
	  re_cell[i][j].G_kl = gKL_RE;

      re_cell[i][j].G_Na  = 90;
      re_cell[i][j].G_K   = 10;

      re_cell[i][j].G_Ca  = 1.3;

      re_cell[i][j].G_AHP = 0.2;
      re_cell[i][j].G_CAN = 0.2;

      re_cell[i][j].E_l = -60;
      re_cell[i][j].Taur = 100;
      re_cell[i][j].D    = 0.5;

      re_cell[i][j].INaK::Vtr  = -40;
      re_cell[i][j].INaK::VtrK = -40;

      re_cell[i][j].INaK::k_h = 1;
      re_cell[i][j].INaK::k_n = 1;

      re_cell[i][j].IT_RE::Shift = -3;

	}


// Change GABA reversal potentials for RE cells

  for(i=0; i<N_RE_X; i++)
	for(j=0;j<N_RE_Y; j++)
		for(k=0; k<N_IN; k++) {
		gaba_RE_RE[i][j][k].E_GABA = E_GABA_RE;
	}

// Change Random input rate to TC cells
  for(i=0; i<N_TC1_X; i++)
	for(j=0;j<N_TC1_Y; j++)
	   for (k=0; k<N_INPUT_TC1;k++){
       ampa_ext_tc1[i][j][k].w = RANDOM_F_TC;
   }

  for(i=0; i<N_TC2_X; i++)
	for(j=0;j<N_TC2_Y; j++)
	  for (k=0; k<N_INPUT_TC2;k++){
       ampa_ext_tc2[i][j][k].w = RANDOM_F_TC;
   }


// =====================================================
//         Index Cells
//======================================================
// Index TC1
 for(i=0; i<N_TC1_X; i++)
   for(j=0;j<N_TC1_Y; j++)
      i_tc1[i][j] = (j + i*N_TC1_Y) * EQ_TC;

// Index TC2
 for(i=0; i<N_TC2_X; i++)
   for(j=0;j<N_TC2_Y; j++)
      i_tc2[i][j] = N_EQ_TC1 + (j + i*N_TC2_Y) * EQ_TC;

 // Indices INs
 for(i=0; i<N_IN_X; i++)
   for(j=0;j<N_IN_Y; j++)
	  i_in[i][j] = N_EQ_TC1 + N_EQ_TC2 + (j + i*N_IN_Y) * EQ_IN;

 // Indices REs
 for(i=0; i<N_RE_X; i++)
   for(j=0;j<N_RE_Y; j++)
	i_re[i][j] = N_EQ_TC1 + N_EQ_TC2 + N_EQ_IN + (j + i*N_RE_Y) * EQ_RE;


 // =====================================================
 //          Initialization
 //======================================================
 Stim.init();
 Stim_Spin.init();

 // Initialize arrays
for(i=N_EQ-1; i>=0; --i)
   y_ini[i] = 0, f_ini[i] = 0, y1[i] = 0, y2[i] = 0;


  // Initialize Cells
 for(i=0; i<N_TC1_X; i++)
   for(j=0; j<N_TC1_Y; j++){
	tc1_cell[i][j].init(y_ini + i_tc1[i][j]);
 }

 for(i=0; i<N_TC2_X; i++)
   for(j=0; j<N_TC2_Y; j++){
	tc2_cell[i][j].init(y_ini + i_tc2[i][j]);
 }

 for(i=0; i<N_IN_X; i++)
   for(j=0; j<N_IN_Y; j++){
	in_cell[i][j].init(y_ini + i_in[i][j]);
 }

 // RE
 for(i=0; i<N_RE_X; i++)
   for(j=0; j<N_RE_Y; j++){
	re_cell[i][j].init(y_ini + i_re[i][j]);
 }


// Initialize the connection arrays
 for(i=0; i<N_TC1_X; i++)
    for(j=0; j<N_TC1_Y; j++)
   	 conn_tc1[i][j].INIT();

 for(i=0; i<N_TC2_X; i++)
    for(j=0; j<N_TC2_Y; j++)
   	 conn_tc2[i][j].INIT();

 for(i=0; i<N_IN_X; i++)
    for(j=0; j<N_IN_Y; j++)
   	 conn_in[i][j].INIT();

 for(i=0; i<N_RE_X; i++)
    for(j=0; j<N_RE_Y; j++)
   	 conn_re[i][j].INIT();

 // Initialize spike times
  for(i=0; i<N_TC1_X; i++)
 	for(j=0; j<N_TC1_Y; j++)
       TC1_SPT[i][j] = -10;

  for(i=0; i<N_TC2_X; i++)
 	for(j=0; j<N_TC2_Y; j++){
       TC2_SPT[i][j] = -10;
       TC2_E[i][j] = 1;
 	}

  for(i=0; i<N_IN_X; i++)
 	for(j=0; j<N_IN_Y; j++)
       IN_SPT[i][j] = -10;

  for(i=0; i<N_RE_X; i++)
 	for(j=0; j<N_RE_Y; j++) {
       RE_SPT[i][j] = -10;
       RE_E[i][j] = 1;
 	}




//==================================================================================
//           Connectivity of the Network
//==================================================================================

// HTC cells to HTC cells connections (Gap Junctions)
    for (i=0; i<N_TC1_X; i++)
      for (j=0; j<N_TC1_Y; j++) {
        POST = i*N_TC1_Y+j;

        for (i1=0; i1<N_TC1_X; i1++)
         for (j1=0; j1<N_TC1_Y; j1++) {

           PRE = i1*N_TC1_Y+j1;

           d = sqrt((i1-i)*(i1-i)+(j1-j)*(j1-j));  //

     		 if ( (PRE!=POST) && (conn_tc1[i][j].C_TC1_TC1[i1][j1]==0) )  {
    	       if ( ((rand()/(RAND_MAX + 1.0)) < P_TC1_TC1_GAP) && (d <= D_TC1_TC1) ) {

               // for the post_junctional cell
    	    	 k = conn_tc1[i][j].N_Pre_TC1;
    	    	 conn_tc1[i][j].C_TC1_TC1[i1][j1] = 1;
                 conn_tc1[i][j].Pre_TC1_X[k] = i1;
    	         conn_tc1[i][j].Pre_TC1_Y[k] = j1;
    	         conn_tc1[i][j].Pre_TC1_D[k] = d;
    	         conn_tc1[i][j].N_Pre_TC1++;


               // for the pre_junctional cell
    	         k = conn_tc1[i1][j1].N_Pre_TC1;
      	    	 conn_tc1[i1][j1].C_TC1_TC1[i][j] = 1;
                 conn_tc1[i1][j1].Pre_TC1_X[k] = i;
      	         conn_tc1[i1][j1].Pre_TC1_Y[k] = j;
      	         conn_tc1[i1][j1].Pre_TC1_D[k] = d;
      	         conn_tc1[i1][j1].N_Pre_TC1++;

           }
   	     }
       }
   }

// RTC cells to HTC cells connections (Gap Junctions)
    // Select a small subset of RTC cells to form gap jucntion with HTC cells
     for (i=0; i<N_TC2_X; i++)
        for (j=0; j<N_TC2_Y; j++) {
    	   if( rand()/(RAND_MAX + 1.0) < P_TC2_GAP)   {
    	   conn_tc2[i][j].F_TC2_TC1 = 1;
    	   }
       }

      for (i=0; i<N_TC2_X; i++)    // TC2 cells
         for (j=0; j<N_TC2_Y; j++) {
           POST = i*N_TC2_Y+j;

            if( conn_tc2[i][j].F_TC2_TC1 == 1)   {
               for (i1=0; i1<N_TC1_X; i1++)      // TC1 cells
                for (j1=0; j1<N_TC1_Y; j1++) {
                   PRE = i1*N_TC1_Y+j1;

                   mi = (N_TC1_X-1.0)/(N_TC2_X-1.0)*i;  // Projected x index in the HTC plane
                   mj = (N_TC1_Y-1.0)/(N_TC2_Y-1.0)*j;  // Projected y index in the HTC plane

                   d = sqrt( (i1-mi)*(i1-mi)+(j1-mj)*(j1-mj) );

           		   if ( (d <= D_TC2_TC1) && ((rand()/(RAND_MAX + 1.0)) < P_TC2_TC1_GAP) ) {

                         // for the post_junctional cell
              	    	 k = conn_tc2[i][j].N_Pre_TC1;
              	    	 conn_tc2[i][j].Pre_TC1_X[k] = i1;
              	         conn_tc2[i][j].Pre_TC1_Y[k] = j1;
              	         conn_tc2[i][j].Pre_TC1_D[k] = d;
              	         conn_tc2[i][j].N_Pre_TC1++;

                         // for the pre_junctional cell
              	         k = conn_tc1[i1][j1].N_Pre_TC2;
                	     conn_tc1[i1][j1].Pre_TC2_X[k] = i;
                	     conn_tc1[i1][j1].Pre_TC2_Y[k] = j;
                	     conn_tc1[i1][j1].Pre_TC2_D[k] = d;
                	     conn_tc1[i1][j1].N_Pre_TC2++;
              		    }
              	     }

              }
         }

// RTC cells to RTC cells connection (Chemical synapses between RTC cells)
// This connection is NOT active!!!

    for (i=0; i<N_TC2_X; i++)
        for (j=0; j<N_TC2_Y; j++) {
           POST = i*N_TC2_Y+j;

          for (i1=0; i1<N_TC2_X; i1++)
           for (j1=0; j1<N_TC2_Y; j1++) {
             PRE = i1*N_TC2_Y+j1;
             d = sqrt((i1-i)*(i1-i)+(j1-j)*(j1-j));

       		 if ( PRE!=POST  )
      	       if ( ((rand()/(RAND_MAX + 1.0)) < P_TC2_TC2) && (d <= D_TC2_TC2) ) {

      	    	 k = conn_tc2[i][j].N_Pre_TC2;
      	    	// conn_tc2[i][j].C_TC1_TC1[i1][j1] = 1;
                 conn_tc2[i][j].Pre_TC2_X[k] = i1;
      	         conn_tc2[i][j].Pre_TC2_Y[k] = j1;
      	         conn_tc2[i][j].N_Pre_TC2++;

      		    }
      	      }
      	   }



  // HTC cells to IN connections
      for (i=0; i<N_IN_X; i++)
        for (j=0; j<N_IN_Y; j++) {

          POST = i*N_IN_Y+j;

          for (i1=0; i1<N_TC1_X; i1++)
           for (j1=0; j1<N_TC1_Y; j1++) {
            PRE = i1*N_TC1_Y+j1;

       	       if ( (rand()/(RAND_MAX + 1.0)) < P_TC1_IN ) {

       	    	 k = conn_in[i][j].N_Pre_TC1;
                 conn_in[i][j].Pre_TC1_X[k] = i1;
      	         conn_in[i][j].Pre_TC1_Y[k] = j1;
      	         conn_in[i][j].N_Pre_TC1++;
            }

      	   }
      	 }


 // HTC cells to RE cells connections
          for (i=0; i<N_RE_X; i++)
            for (j=0; j<N_RE_Y; j++) {
              POST = i*N_RE_Y+j;

              for (i1=0; i1<N_TC1_X; i1++)
               for (j1=0; j1<N_TC1_Y; j1++) {
                PRE = i1*N_TC1_Y+j1;

           	       if ( (rand()/(RAND_MAX + 1.0)) < P_TC1_RE ) {

           	    	 k = conn_re[i][j].N_Pre_TC1;
                     conn_re[i][j].Pre_TC1_X[k] = i1;
          	         conn_re[i][j].Pre_TC1_Y[k] = j1;
          	         conn_re[i][j].N_Pre_TC1++;
                }
          	   }
        }



 // RTC cells to RE cells connections
     for (i=0; i<N_RE_X; i++)
       for (j=0; j<N_RE_Y; j++)  {

    	 for (i1=0; i1<N_TC2_X; i1++)
   	        for (j1=0; j1<N_TC2_Y; j1++)  {

   	    	if ((rand()/(RAND_MAX + 1.0)) < P_TC2_RE) {
   	    	 k = conn_re[i][j].N_Pre_TC2;
   	    	 conn_re[i][j].Pre_TC2_X[k] = i1;
   	         conn_re[i][j].Pre_TC2_Y[k] = j1;
   	         conn_re[i][j].N_Pre_TC2++;
   	      }
        }
      }


   // INs to RTC cells connections
            for (i=0; i<N_TC2_X; i++)
              for (j=0; j<N_TC2_Y; j++) {
                POST = i*N_TC2_Y+j;

                for (i1=0; i1<N_IN_X; i1++)
                 for (j1=0; j1<N_IN_Y; j1++) {
                  PRE = i1*N_IN_Y+j1;

                    if ( (rand()/(RAND_MAX + 1.0)) < P_IN_TC2 ) {

             	    	 k = conn_tc2[i][j].N_Pre_IN;
                         conn_tc2[i][j].Pre_IN_X[k] = i1;
            	         conn_tc2[i][j].Pre_IN_Y[k] = j1;
            	         conn_tc2[i][j].N_Pre_IN++;
                     }

           	   }
           	 }


 // RE cells to HTC cells connections
         for (i=0; i<N_TC1_X; i++)
           for (j=0; j<N_TC1_Y; j++) {

             for (i1=0; i1<N_RE_X; i1++)
               for (j1=0; j1<N_RE_Y; j1++) {

        	      if ((rand()/(RAND_MAX + 1.0)) < P_RE_TC1) {
        	       k = conn_tc1[i][j].N_Pre_RE;
        	       conn_tc1[i][j].Pre_RE_X[k] = i1;
        	       conn_tc1[i][j].Pre_RE_Y[k] = j1;
        	       conn_tc1[i][j].N_Pre_RE++;
        	     }
             }
          }

// RE cells to RTC cells connections
    for (i=0; i<N_TC2_X; i++)
      for (j=0; j<N_TC2_Y; j++) {

        for (i1=0; i1<N_RE_X; i1++)
          for (j1=0; j1<N_RE_Y; j1++) {

          if ((rand()/(RAND_MAX + 1.0)) < P_RE_TC2) {
   	       k = conn_tc2[i][j].N_Pre_RE;
   	       conn_tc2[i][j].Pre_RE_X[k] = i1;
   	       conn_tc2[i][j].Pre_RE_Y[k] = j1;
   	       conn_tc2[i][j].N_Pre_RE++;
   	       }

        }
     }


 // RE cells to INs connections
        for (i=0; i<N_IN_X; i++)
          for (j=0; j<N_IN_Y; j++) {

            for (i1=0; i1<N_RE_X; i1++)
              for (j1=0; j1<N_RE_Y; j1++) {

              if ((rand()/(RAND_MAX + 1.0)) < P_RE_IN) {
       	       k = conn_in[i][j].N_Pre_RE;
       	       conn_in[i][j].Pre_RE_X[k] = i1;
       	       conn_in[i][j].Pre_RE_Y[k] = j1;
       	       conn_in[i][j].N_Pre_RE++;
       	       }

            }
         }


// RE to RE connection (Gap junctions)
   // Select a small subset of RE cells to form gap jucntions with other neighboring RE cells
     for (i=0; i<N_RE_X; i++)
       for (j=0; j<N_RE_Y; j++) {
      	   if( rand()/(RAND_MAX + 1.0) < P_RE_GAP)   {
         	   conn_re[i][j].F_RE_RE = 1;
       	   }
        }

    for (i=0; i<N_RE_X; i++)
      for (j=0; j<N_RE_Y; j++) {
          POST = i*N_RE_Y+j;

        if( conn_re[i][j].F_RE_RE == 1)   {

           for (i1=0; i1<N_RE_X; i1++)
              for (j1=0; j1<N_RE_Y; j1++) {
                PRE = i1*N_RE_Y+j1;

                d = sqrt( (i1-i)*(i1-i)+(j1-j)*(j1-j) );

       		    if ( (PRE!=POST) && (conn_re[i][j].C_RE_RE[i1][j1]==0) )  {
       	          if ( ((rand()/(RAND_MAX + 1.0)) < P_RE_RE_GAP) && (d <= D_RE_RE) ) {

                   // for the post_junctional cell
        	    	 k = conn_re[i][j].N_Pre_RE_GAP;
        	    	 conn_re[i][j].C_RE_RE[i1][j1] = 1;
                     conn_re[i][j].Pre_RE_X_GAP[k] = i1;
        	         conn_re[i][j].Pre_RE_Y_GAP[k] = j1;
        	         conn_re[i][j].Pre_RE_D_GAP[k] = d;
        	         conn_re[i][j].N_Pre_RE_GAP++;

                   // for the pre_junctional cell
        	         k = conn_re[i1][j1].N_Pre_RE_GAP;
          	    	 conn_re[i1][j1].C_RE_RE[i][j] = 1;
                     conn_re[i1][j1].Pre_RE_X_GAP[k] = i;
          	         conn_re[i1][j1].Pre_RE_Y_GAP[k] = j;
          	         conn_re[i1][j1].Pre_RE_D_GAP[k] = d;
          	         conn_re[i1][j1].N_Pre_RE_GAP++;
        		    }
         	    }
           }
        }
     }



 // RE to RE connections (GABAergic inhibitory synapses)
    for (i=0; i<N_RE_X; i++)
       for (j=0; j<N_RE_Y; j++)  {
    	 POST = i*N_RE_Y+j;

     	 for (i1=0; i1<N_RE_X; i1++)
      	    for (j1=0; j1<N_RE_Y; j1++)  {
               PRE = i1*N_RE_Y+j1;

               d = sqrt( (i1-i)*(i1-i) + (j1-j)*(j1-j) );

               if ( PRE!=POST )  {
      	       if ( (rand()/(RAND_MAX + 1.0)) < P_RE_RE ) {
      	    	 k = conn_re[i][j].N_Pre_RE;
      	    	 conn_re[i][j].Pre_RE_X[k] = i1;
      	         conn_re[i][j].Pre_RE_Y[k] = j1;
      	         conn_re[i][j].N_Pre_RE++;
      	        }
              }
      	    }
      }




// Save network connectivity to file
//HTC-HTC connectivity
  f_connect = fopen("data/Con_HTC_HTC", "w");
  fprintf(f_connect, "The HTC-->HTC Connection Is: \n");

   for(i=0; i<N_TC1_X; i++)
     for(j=0; j<N_TC1_Y; j++) {
   	      fprintf(f_connect, "\n\nHTC[%d][%d] has %d gap junctional connections with HTCs: \n", i,j, conn_tc1[i][j].N_Pre_TC1);

    	  for(i1=0; i1<conn_tc1[i][j].N_Pre_TC1; i1++) {
     	      if (i1%5==0)
    	    	fprintf(f_connect, "\n");
    	      fprintf(f_connect, "HTC(%d, %d) d=%2.1f \t", conn_tc1[i][j].Pre_TC1_X[i1], conn_tc1[i][j].Pre_TC1_Y[i1],
    	    		  conn_tc1[i][j].Pre_TC1_D[i1]);

    	  }
      }

   fclose(f_connect);

// RTC-->HTC connectivity
      f_connect = fopen("data/Con_RTC_HTC", "w");
      fprintf(f_connect, "The RTC-->HTC Connection Is: \n");

       for(i=0; i<N_TC1_X; i++)
         for(j=0; j<N_TC1_Y; j++) {
       	  fprintf(f_connect, "\n\nHTC[%d][%d] has %d gap junctional connections with RTCs: \n", i,j, conn_tc1[i][j].N_Pre_TC2);

        	  for(i1=0; i1<conn_tc1[i][j].N_Pre_TC2; i1++) {
         	      if (i1%5==0)
        	    	fprintf(f_connect, "\n");
        	      fprintf(f_connect, "RTC(%d, %d) d=%3.2f \t", conn_tc1[i][j].Pre_TC2_X[i1], conn_tc1[i][j].Pre_TC2_Y[i1],
        	    		  conn_tc1[i][j].Pre_TC2_D[i1]);

        	  }
          }

       fclose(f_connect);


// HTC-->RTC connectivity
   f_connect = fopen("data/Con_HTC_RTC", "w");
   fprintf(f_connect, "The HTC-->RTC Connection Is: \n");

    for(i=0; i<N_TC2_X; i++)
      for(j=0; j<N_TC2_Y; j++) {
    	  fprintf(f_connect, "\n\nRTC[%d][%d] has %d gap junctional connections with HTCs: \n", i,j, conn_tc2[i][j].N_Pre_TC1);

     	  for(i1=0; i1<conn_tc2[i][j].N_Pre_TC1; i1++) {
      	      if (i1%5==0)
     	    	fprintf(f_connect, "\n");
     	      fprintf(f_connect, "HTC(%d, %d) d=%3.2f \t", conn_tc2[i][j].Pre_TC1_X[i1], conn_tc2[i][j].Pre_TC1_Y[i1],
     	    		  conn_tc2[i][j].Pre_TC1_D[i1]);

     	  }
       }

    fclose(f_connect);

    //RE-RE connectivity
    f_connect = fopen("data/Con_RE_RE", "w");
    fprintf(f_connect, "The RE-->RE Connection Is: \n");

    for(i=0; i<N_RE_X; i++)
       for(j=0; j<N_RE_Y; j++) {
     	      fprintf(f_connect, "\n\nRE[%d][%d] has %d gap junctional connections with REs: \n", i,j, conn_re[i][j].N_Pre_RE_GAP);

        	  for(i1=0; i1<conn_re[i][j].N_Pre_RE_GAP; i1++) {
         	      if (i1%5==0)
        	    	fprintf(f_connect, "\n");
        	      fprintf(f_connect, "RE(%d, %d) d=%2.1f \t", conn_re[i][j].Pre_RE_X_GAP[i1], conn_re[i][j].Pre_RE_Y_GAP[i1],
        	    	conn_re[i][j].Pre_RE_D_GAP[i1]);

        	  }
          }

       fclose(f_connect);


 // Open output files

 f_tc1 = fopen("data/tc1_all", "w");
 f_tc2 = fopen("data/tc2_all", "w");
 f_in  = fopen("data/in_all", "w");
 f_re  = fopen("data/re_all", "w");
 f_E   = fopen("data/TC2_E", "w");
 f_W_HTC_IN  = fopen("data/W_HTC_IN", "w");
 f_E_TC_RE  = fopen("data/E_TC_RE", "w");
 f_E_RE_TC  = fopen("data/E_RE_TC", "w");
 f_I    = fopen("data/I_STIM", "w");
 f_SI   = fopen("data/I_SPIN", "w");


 // Stim current
 fprintf(f_I, "%f \t", Stim.I);
 fprintf(f_I, "\n");

 fprintf(f_SI, "%f \t", Stim_Spin.I);
 fprintf(f_SI, "\n");


 // W_HTC_IN
 fprintf(f_W_HTC_IN, "%f \t", ampa_TC1_IN[0][0][0].E);
 fprintf(f_W_HTC_IN, "\n");


 // RTC-RE short-term depression
 fprintf(f_E_TC_RE, "%f \t", ampa_TC2_RE[0][0][0].E);
 fprintf(f_E_TC_RE, "\n");

 // RE-RTC short-term depression
 fprintf(f_E_RE_TC, "%f \t", gaba_RE_TC2[0][0][0].E);
 fprintf(f_E_RE_TC, "\n");


// HTC cells
  fprintf(f_tc1,"%f \t",t);
  for (i=0; i<N_TC1_X; i++)
     for (j=0; j<N_TC1_Y; j++) {
       fprintf(f_tc1,"%f \t", y_ini[i_tc1[i][j]]);
  }
  fprintf(f_tc1,"\n");

// RTC cells
  fprintf(f_tc2,"%f \t",t);
  for (i=0; i<N_TC2_X; i++)
     for (j=0; j<N_TC2_Y; j++) {
       fprintf(f_tc2,"%f \t", y_ini[i_tc2[i][j]]);
       fprintf(f_E,"%f \t", TC2_E[i][j]);
  }
  fprintf(f_tc2,"\n");
  fprintf(f_E,"\n");


 // INs
  fprintf(f_in,"%f \t",t);
  for (i=0; i<N_IN_X; i++)
     for (j=0; j<N_IN_Y; j++) {
       fprintf(f_in,"%f \t", y_ini[i_in[i][j]]);
  }
  fprintf(f_in,"\n");

 // RE cells
  fprintf(f_re,"%f \t",t);
  for (i=0; i<N_RE_X; i++)
     for (j=0; j<N_RE_Y; j++) {
        fprintf(f_re,"%f \t", y_ini[i_re[i][j]]);
  }
  fprintf(f_re,"\n");




 // Save spike times
 // HTC
 for(i=0; i<N_TC1_X; i++){
      for(j=0; j<N_TC1_Y; j++){
        sprintf(filename, "%sTC1_%d_%d", filepath, i, j);
        f_spike = fopen(filename, "w");
        fclose(f_spike);
       }
     }

 // RTC
 for(i=0; i<N_TC2_X; i++){
      for(j=0; j<N_TC2_Y; j++){
        sprintf(filename, "%sTC2_%d_%d", filepath, i, j);
        f_spike = fopen(filename, "w");
        fclose(f_spike);
       }
     }

// IN
  for(i=0; i<N_IN_X; i++){
      for(j=0; j<N_IN_Y; j++){
        sprintf(filename, "%sIN_%d_%d", filepath, i, j);
        f_spike = fopen(filename, "w");
        fclose(f_spike);
       }
    }

 // RE
   for(i=0; i<N_RE_X; i++){
       for(j=0; j<N_RE_Y; j++){
          sprintf(filename, "%sRE_%d_%d", filepath, i, j);
          f_spike = fopen(filename, "w");
          fclose(f_spike);
         }
      }

 // Configure Stimulation
  
  printf("\n Now starting simulation: t= %lf: tmax= %lf \n", t,t_end);
  
  ii = 0;
  tau = h;
  
  while( t <= t_end) {

  // Record spike times
  // HTC
	for(i=0; i<N_TC1_X; i++)
	  for(j=0; j<N_TC1_Y; j++)
	    if (y_ini[i_tc1[i][j]] > 0 && (t-TC1_SPT[i][j]) > 3) {
	 	     TC1_SPT[i][j] = t;
	         sprintf(filename, "%sTC1_%d_%d", filepath, i, j);
	         f_spike = fopen(filename, "a");
	         fprintf(f_spike,"%f\n", t);
	         fclose(f_spike);
	 	 }

 // RTC
	for(i=0; i<N_TC2_X; i++)
	  for(j=0; j<N_TC2_Y; j++)
 	    if (y_ini[i_tc2[i][j]] > 0 && (t-TC2_SPT[i][j]) > 3) {
 	    	 TC2_E[i][j] = 1 - (1 - TC2_E[i][j]*(1-0.07)) * exp(-(t-TC2_SPT[i][j])/700);
	 	     TC2_SPT[i][j] = t;
	         sprintf(filename, "%sTC2_%d_%d", filepath, i, j);
	         f_spike = fopen(filename, "a");
	         fprintf(f_spike,"%f\n", t);
	         fclose(f_spike);
	 	 }
 // IN
    for(i=0; i<N_IN_X; i++)
	  for(j=0; j<N_IN_Y; j++)
	    if (y_ini[i_in[i][j]] > 0 && (t-IN_SPT[i][j]) > 3) {
	 	    IN_SPT[i][j] = t;
	        sprintf(filename, "%sIN_%d_%d", filepath, i, j);
	        f_spike = fopen(filename, "a");
	        fprintf(f_spike,"%f\n", t);
	        fclose(f_spike);
	 	 }

 // RE
    for(i=0; i<N_RE_X; i++)
   	  for(j=0; j<N_RE_Y; j++)
   	    if (y_ini[i_re[i][j]] > 0 && (t-RE_SPT[i][j]) > 3) {
   	    	RE_E[i][j] = 1 - (1 - RE_E[i][j]*(1-0.07)) * exp(-(t-RE_SPT[i][j])/700);
   	 	    RE_SPT[i][j] = t;
   	        sprintf(filename, "%sRE_%d_%d", filepath, i, j);
   	        f_spike = fopen(filename, "a");
   	        fprintf(f_spike,"%f\n", t);
   	        fclose(f_spike);
   	 	 }

 //======================================================================

 // Constant current injection
	if (FLAG_CURRENT_INJECTION == 1)   {
       for(i=0; i<N_TC1_X; i++)
	      for(j=0; j<N_TC1_Y; j++) {
	         tc1_cell[i][j].setStim(input_tc1[i][j]);
	      }

       for(i=0; i<N_TC2_X; i++)
	      for(j=0; j<N_TC2_Y; j++) {
	         tc2_cell[i][j].setStim(input_tc2[i][j]);
	      }

	   for(i=0; i<N_IN_X; i++)
	      for(j=0; j<N_IN_Y; j++) {
	         in_cell[i][j].setStim(input_in[i][j]);
	      }
	}


 // Spindle-triggering input
   if (FLAG_SPINDLE_INPUT == 1) {

	 Stim_Spin.calc(t);

	 for(i=0; i<SPIN_X; i++)
	    for(j=0; j<SPIN_Y; j++) {
	      re_cell[i][j].setStim(Stim_Spin.I);
        }
    }


  // Square Pulse Stimulation
   if (FLAG_STIMULATION_LGN_PULSE == 1) {
	   Stim.calc(t);

	   for(i=0; i<N_TC1_X; i++)
	      for(j=0; j<N_TC1_Y; j++) {
	         tc1_cell[i][j].setStim(Stim.I);
	      }

	   for(i=0; i<N_TC2_X; i++)
	      for(j=0; j<N_TC2_Y; j++) {
	      tc2_cell[i][j].setStim(Stim.I);
	     }

	   for(i=0; i<N_IN_X; i++)
	      for(j=0; j<N_IN_Y; j++) {
	      in_cell[i][j].setStim(Stim.I);
	      }
      }


   if (FLAG_STIMULATION_TRN_PULSE == 1) {
	   Stim.calc(t);

	   for(i=0; i<N_RE_X; i++)
	     for(j=0; j<N_RE_Y; j++) {
	        if (FLAG_SPINDLE_INPUT == 1)
	           re_cell[i][j].setStim(Stim.I + Stim_Spin.I);
	        else
	    	   re_cell[i][j].setStim(Stim.I);
       }
    }



 // ==================================================================


    rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);

    t = t + tau;
    ii++;
    
    
  
   if( ((ii/(10000))*(10000) == ii)) {
      printf("\n Time: t= %lf \n", t);
    }


// Save data with a resoluton of DT = 0.2 ms

   if( ((ii/(10))*(10) == ii)){
     // HTC
	   fprintf(f_tc1,"%f \t",t);
	   for (i=0; i<N_TC1_X; i++)
	      for (j=0; j<N_TC1_Y; j++) {
	        fprintf(f_tc1,"%f \t", y_ini[i_tc1[i][j]]);
	   }
	   fprintf(f_tc1,"\n");

     // RTC
	   fprintf(f_tc2,"%f \t",t);
	   for (i=0; i<N_TC2_X; i++)
	      for (j=0; j<N_TC2_Y; j++) {
	        fprintf(f_tc2,"%f \t", y_ini[i_tc2[i][j]]);
	        fprintf(f_E,"%f \t", TC2_E[i][j]);
	   }
	   fprintf(f_tc2,"\n");
	   fprintf(f_E,"\n");

	   // IN
	    fprintf(f_in,"%f \t",t);
	    for (i=0; i<N_IN_X; i++)
	       for (j=0; j<N_IN_Y; j++) {
	         fprintf(f_in,"%f \t", y_ini[i_in[i][j]]);
	    }
	    fprintf(f_in,"\n");

		// RE
		fprintf(f_re,"%f \t",t);
		for (i=0; i<N_RE_X; i++)
		   for (j=0; j<N_RE_Y; j++) {
		      fprintf(f_re,"%f \t", y_ini[i_re[i][j]]);
		}
		fprintf(f_re,"\n");

		// STD
         fprintf(f_W_HTC_IN, "%f \t", ampa_TC1_IN[0][0][0].E);
	     fprintf(f_W_HTC_IN, "\n");

		 fprintf(f_E_TC_RE, "%f \t", ampa_TC2_RE[0][0][0].E);
		 fprintf(f_E_TC_RE, "\n");

		 fprintf(f_E_RE_TC, "%f \t", gaba_RE_TC2[0][0][0].E);
		 fprintf(f_E_RE_TC, "\n");

		// Stimulation currents
		 fprintf(f_I, "%f \t", Stim.I);
		 fprintf(f_I, "\n");

		 fprintf(f_SI, "%f \t", Stim_Spin.I);
		 fprintf(f_SI, "\n");


    }

    

}
 
 
 
 // Close files   
  fclose(f_tc1);
  fclose(f_tc2);
  fclose(f_in);
  fclose(f_re);
  fclose(f_I);
  fclose(f_SI);
  fclose(f_E);
  fclose(f_W_HTC_IN);
  fclose(f_E_TC_RE);
  fclose(f_E_RE_TC);

  
  // Print the end time of simulation
   printf("\n\nProgram ended on ");
   t_start =time(0);

   mytime = ctime(&t_start);
   printf(mytime);

}





/*******************************************************************
 // void fun()
 // Evaluate right-hand side of system of coupled ODEs
 *******************************************************************/

void fun(int unsigned NEQ, double x, double *y_ini, double *f_ini){
  
  // Local variables
  int i,j,i1,j1,I1,I2,k, id; // generic indices
  int post, pre;
 // int k_ampa, k_nmda; // local counter of number of syn per cell
  int N;
  double T_Start_Synapse;


  if (OSC == 2) {
	T_Start_Synapse = 500;
  } else  {
    T_Start_Synapse = 0; }
  
  //===================================================================
  // Begin: Cells
  
  // Thalamus
  if (TC1_EX==1) {
    for(i=0; i<N_TC1_X; i++)
       for(j=0; j<N_TC1_Y; j++)
	    tc1_cell[i][j].calc(x, y_ini+i_tc1[i][j], f_ini+i_tc1[i][j]);
  }
  
  if (TC2_EX==1) {
    for(i=0; i<N_TC2_X; i++)
      for(j=0; j<N_TC2_Y; j++)
	    tc2_cell[i][j].calc(x, y_ini+i_tc2[i][j], f_ini+i_tc2[i][j]);
  }

  if (IN_EX==1) {
   for(i=0; i<N_IN_X; i++)
     for(j=0; j<N_IN_Y; j++)
	 in_cell[i][j].calc(x, y_ini + i_in[i][j], f_ini + i_in[i][j]);
  }

  if (RE_EX==1) {
    for(i=0; i<N_RE_X; i++)
      for(j=0; j<N_RE_Y; j++) {
	  re_cell[i][j].calc(x, y_ini + i_re[i][j], f_ini + i_re[i][j]);
    }
  }

 // End: Cells
 // ==================================================================

  // External random background inputs to TC1s
  if(FLAG_RANDOM_INPUT == 1) {
  // TC1
	if (TC1_EX==1) {

      for(i=0; i<N_TC1_X; i++)
        for(j=0; j<N_TC1_Y; j++) {
        	id = i_tc1[i][j];

        	for(k=0; k<N_INPUT_TC1; k++) {
 	            ampa_ext_tc1[i][j][k].calc(g_Extern_TC1, x);
 	            f_ini[id] = f_ini[id] - (ampa_ext_tc1[i][j][k].g * y_ini[id])/(tc1_cell[i][j].S_TC*1e3);
            }
      }
	}

   // TC2
	if (TC2_EX==1) {
      for(i=0; i<N_TC2_X; i++)
        for(j=0; j<N_TC2_Y; j++) {
        	id = i_tc2[i][j];

           for(k=0; k<N_INPUT_TC2; k++){
              ampa_ext_tc2[i][j][k].calc(g_Extern_TC2, x);
              f_ini[id] = f_ini[id] - (ampa_ext_tc2[i][j][k].g * y_ini[id])/(tc2_cell[i][j].S_TC*1e3);
          }
       }
	}


  // IN
	if (IN_EX==1) {
      for(i=0; i<N_IN_X; i++)
        for(j=0; j<N_IN_Y; j++) {
          id = i_in[i][j];
          ampa_ext_in[i][j].calc(g_Extern_IN, x);
          f_ini[id] = f_ini[id] - (ampa_ext_in[i][j].g * y_ini[id])/(in_cell[i][j].S_IN*1e3);
       }
	}


  // RE
	if (RE_EX==1) {
      for(i=0; i<N_RE_X; i++)
        for(j=0; j<N_RE_Y; j++) {
          id = i_re[i][j];
          ampa_ext_re[i][j].calc(g_Extern_RE, x);
          f_ini[id] = f_ini[id] - (ampa_ext_re[i][j].g * y_ini[id])/(re_cell[i][j].S_RE*1e3);
       }
	}


  }

//==========================================
//            Synaptic Connections
//==========================================

 // HTC-HTC Gap junctional connections
  if (x >= T_Start_Synapse)  {

  if (TC1_EX==1) {
    if (FLAG_TC1_TC1_GAP == 1) {
       for(i=0; i<N_TC1_X; i++)
         for(j=0; j<N_TC1_Y; j++) {

          post = i_tc1[i][j];
          N  = conn_tc1[i][j].N_Pre_TC1;

          if (N > 0){
        	for (k=0; k<N; k++) {
     		  i1 = conn_tc1[i][j].Pre_TC1_X[k];
     		  j1 = conn_tc1[i][j].Pre_TC1_Y[k];
     		  pre = i_tc1[i1][j1];
     		  f_ini[post] = f_ini[post] - (y_ini[post] - y_ini[pre])/R_TC1_TC1/(tc1_cell[i][j].S_TC*1e3);

        	}
         }
      }
    }
  }

 // RTC-HTC Gap junctional connections

  if (TC1_EX==1 && TC2_EX==1) {
    if (FLAG_TC2_TC1_GAP == 1) {

    // For HTC CELLS
       for(i=0; i<N_TC1_X; i++)
         for(j=0; j<N_TC1_Y; j++) {

          post = i_tc1[i][j];
          N  = conn_tc1[i][j].N_Pre_TC2;

          if (N > 0){
        	for (k=0; k<N; k++) {
     		  i1 = conn_tc1[i][j].Pre_TC2_X[k];
     		  j1 = conn_tc1[i][j].Pre_TC2_Y[k];
     		  pre = i_tc2[i1][j1];
     		  f_ini[post] = f_ini[post] - (y_ini[post] - y_ini[pre])/R_TC1_TC2/(tc1_cell[i][j].S_TC*1e3);
        	}
         }
        }

       // For RTC CELLS
          for(i=0; i<N_TC2_X; i++)
            for(j=0; j<N_TC2_Y; j++) {

             post = i_tc2[i][j];
             N  = conn_tc2[i][j].N_Pre_TC1;

             if (N > 0){
           	for (k=0; k<N; k++) {
        		  i1 = conn_tc2[i][j].Pre_TC1_X[k];
        		  j1 = conn_tc2[i][j].Pre_TC1_Y[k];
        		  pre = i_tc1[i1][j1];
        		  f_ini[post] = f_ini[post] - (y_ini[post] - y_ini[pre])/R_TC1_TC2/(tc2_cell[i][j].S_TC*1e3);
           	}
            }
           }

    }
  }


  // RTC-RTC chemical synaptic connections (NOT ACTIVE!)
  if (TC2_EX==1) {
   if (FLAG_TC2_TC2 == 1) {

     for(i=0; i<N_TC2_X; i++)
       for(j=0; j<N_TC2_Y; j++){

        post = i_tc2[i][j];
        N = conn_tc2[i][j].N_Pre_TC2;

        if(N > 0) {
          for (k=0; k<N; k++) {

           i1 = conn_tc2[i][j].Pre_TC2_X[k];
           j1 = conn_tc2[i][j].Pre_TC2_Y[k];

   		   ampa_TC2_TC2[i][j][k].calc(g_AMPA_TC2_TC2, x, y_ini[post], TC2_SPT[i1][j1], STD_TC2_TC2);
   		   nmda_TC2_TC2[i][j][k].calc(g_NMDA_TC2_TC2, x, y_ini[post], TC2_SPT[i1][j1], STD_TC2_TC2);

    	   f_ini[post] = f_ini[post] - ampa_TC2_TC2[i][j][k].I/(tc2_cell[i][j].S_TC*1e3);
   	       f_ini[post] = f_ini[post] - nmda_TC2_TC2[i][j][k].I/(tc2_cell[i][j].S_TC*1e3);
        }
      }
    }
   }
  }


  // HTC-->IN
  if (TC1_EX==1 && IN_EX ==1) {
   if (FLAG_TC1_IN == 1) {

     for(i=0; i<N_IN_X; i++)
       for(j=0; j<N_IN_Y; j++){

        post = i_in[i][j];
        N = conn_in[i][j].N_Pre_TC1;

        if(N > 0) {
          for (k=0; k<N; k++) {

           i1 = conn_in[i][j].Pre_TC1_X[k];
           j1 = conn_in[i][j].Pre_TC1_Y[k];

   		   ampa_TC1_IN[i][j][k].calc(g_AMPA_TC1_IN, x, y_ini[post], TC1_SPT[i1][j1], STD_TC1_IN);
   		   nmda_TC1_IN[i][j][k].calc(g_NMDA_TC1_IN, x, y_ini[post], TC1_SPT[i1][j1], STD_TC1_IN);

   		   f_ini[post] = f_ini[post] - ampa_TC1_IN[i][j][k].I/(in_cell[i][j].S_IN*1e3);
   	       f_ini[post] = f_ini[post] - nmda_TC1_IN[i][j][k].I/(in_cell[i][j].S_IN*1e3);
         }
        }
      }
    }
  }


  // HTC-->RE
  if (TC1_EX==1 && RE_EX ==1) {
   if (FLAG_TC1_RE == 1) {

     for(i=0; i<N_RE_X; i++)
       for(j=0; j<N_RE_Y; j++){

        post = i_re[i][j];
        N = conn_re[i][j].N_Pre_TC1;

        if(N > 0) {
          for (k=0; k<N; k++) {

           i1 = conn_re[i][j].Pre_TC1_X[k];
           j1 = conn_re[i][j].Pre_TC1_Y[k];

   		   ampa_TC1_RE[i][j][k].calc(g_AMPA_TC1_RE, x, y_ini[post], TC1_SPT[i1][j1], STD_TC1_RE);
   		   nmda_TC1_RE[i][j][k].calc(g_NMDA_TC1_RE, x, y_ini[post], TC1_SPT[i1][j1], STD_TC1_RE);

   		  f_ini[post] = f_ini[post] - ampa_TC1_RE[i][j][k].I/(re_cell[i][j].S_RE*1e3);
   	      f_ini[post] = f_ini[post] - nmda_TC1_RE[i][j][k].I/(re_cell[i][j].S_RE*1e3);
         }
       }
      }
    }
   }

 // IN-RTC
  if (IN_EX==1 && TC2_EX ==1) {
   if (FLAG_IN_TC2 == 1) {

     for(i=0; i<N_TC2_X; i++)
       for(j=0; j<N_TC2_Y; j++) {

        post = i_tc2[i][j];
        N = conn_tc2[i][j].N_Pre_IN;

        if (N > 0){
   	      for (k=0; k<N; k++) {
   		   i1 = conn_tc2[i][j].Pre_IN_X[k];
   		   j1 = conn_tc2[i][j].Pre_IN_Y[k];
   		   gaba_IN_TC2[i][j][k].calc(g_GABA_IN_TC2, x, y_ini[post], IN_SPT[i1][j1], STD_IN_TC2);
   	       f_ini[post] = f_ini[post] - gaba_IN_TC2[i][j][k].I/(tc2_cell[i][j].S_TC*1e3);
        }
       }
     }
   }
  }


// RTC-RE
  if (TC2_EX==1 && RE_EX ==1) {
   if (FLAG_TC2_RE == 1) {
    for(i=0; i<N_RE_X; i++)
      for(j=0; j<N_RE_Y; j++){

          post = i_re[i][j];
          N = conn_re[i][j].N_Pre_TC2;

         if(N > 0) {
           for (k=0; k<N; k++) {

             i1 = conn_re[i][j].Pre_TC2_X[k];
             j1 = conn_re[i][j].Pre_TC2_Y[k];

             ampa_TC2_RE[i][j][k].calc(g_AMPA_TC2_RE, x, y_ini[post], TC2_SPT[i1][j1], STD_TC2_RE);
   		     nmda_TC2_RE[i][j][k].calc(g_NMDA_TC2_RE, x, y_ini[post], TC2_SPT[i1][j1], STD_TC2_RE);

   		     f_ini[post] = f_ini[post] - ampa_TC2_RE[i][j][k].I/(re_cell[i][j].S_RE*1e3);
   	         f_ini[post] = f_ini[post] - nmda_TC2_RE[i][j][k].I/(re_cell[i][j].S_RE*1e3);
           }
         }
       }
     }
   }

  // RE-HTC
  if (TC1_EX==1 && RE_EX ==1) {
     if (FLAG_RE_TC1 == 1) {
       for(i=0; i<N_TC1_X; i++)
         for(j=0; j<N_TC1_Y; j++) {

          post = i_tc1[i][j];
          N = conn_tc1[i][j].N_Pre_RE;

          if (N > 0){
     	      for (k=0; k<N; k++) {
     		   i1 = conn_tc1[i][j].Pre_RE_X[k];
     		   j1 = conn_tc1[i][j].Pre_RE_Y[k];
     		   gaba_RE_TC1[i][j][k].calc(g_GABA_RE_TC1, x, y_ini[post], RE_SPT[i1][j1], STD_RE_TC1);
               f_ini[post] = f_ini[post] - gaba_RE_TC1[i][j][k].I/(tc1_cell[i][j].S_TC*1e3);
            }
          }
         }
      }
     }

  // RE-RTC
  if (TC2_EX==1 && RE_EX ==1) {
     if (FLAG_RE_TC2 == 1) {
       for(i=0; i<N_TC2_X; i++)
         for(j=0; j<N_TC2_Y; j++) {

          post = i_tc2[i][j];
          N = conn_tc2[i][j].N_Pre_RE;

          if (N > 0){
     	      for (k=0; k<N; k++) {
     		   i1 = conn_tc2[i][j].Pre_RE_X[k];
     		   j1 = conn_tc2[i][j].Pre_RE_Y[k];

		       gaba_RE_TC2[i][j][k].calc(g_GABA_RE_TC2, x, y_ini[post], RE_SPT[i1][j1], STD_RE_TC2);
               f_ini[post] = f_ini[post] - gaba_RE_TC2[i][j][k].I/(tc2_cell[i][j].S_TC*1e3);
            }
          }
        }
      }
     }



 // RE-IN
 if (RE_EX ==1 && IN_EX ==1) {
     if (FLAG_RE_IN == 1) {

       for(i=0; i<N_IN_X; i++)
         for(j=0; j<N_IN_Y; j++) {

          post = i_in[i][j];
          N = conn_in[i][j].N_Pre_RE;

          if (N > 0){
     	      for (k=0; k<N; k++) {
     		   i1 = conn_in[i][j].Pre_RE_X[k];
     		   j1 = conn_in[i][j].Pre_RE_Y[k];
     		   gaba_RE_IN[i][j][k].calc(g_GABA_RE_IN, x, y_ini[post], RE_SPT[i1][j1], STD_RE_IN);
     	       f_ini[post] = f_ini[post] - gaba_RE_IN[i][j][k].I/(in_cell[i][j].S_IN*1e3);
           }
         }
       }
     }
   }


// RE-RE gap junctions
   if (RE_EX==1) {
     if (FLAG_RE_RE_GAP == 1) {
        for(i=0; i<N_RE_X; i++)
          for(j=0; j<N_RE_Y; j++) {

           post = i_re[i][j];
           N  = conn_re[i][j].N_Pre_RE_GAP;

          if (N > 0){
          	for (k=0; k<N; k++) {
    		    i1 = conn_re[i][j].Pre_RE_X_GAP[k];
    		    j1 = conn_re[i][j].Pre_RE_Y_GAP[k];
    		    pre = i_re[i1][j1];
    		    f_ini[post] = f_ini[post] - (y_ini[post] - y_ini[pre])/R_RE_RE/(re_cell[i][j].S_RE*1e3);
         	}
          }
        }
      }
   }


// RE-RE chemical synapses
 if (RE_EX ==1) {
     if (FLAG_RE_RE == 1) {
       for(i=0; i<N_RE_X; i++)
         for(j=0; j<N_RE_Y; j++) {

          post = i_re[i][j];
          N = conn_re[i][j].N_Pre_RE;

          if (N > 0){
     	      for (k=0; k<N; k++) {
     		   i1 = conn_re[i][j].Pre_RE_X[k];
     		   j1 = conn_re[i][j].Pre_RE_Y[k];
     		   gaba_RE_RE[i][j][k].calc(g_GABA_RE_RE, x, y_ini[post], RE_SPT[i1][j1], STD_RE_RE);
     	       f_ini[post] = f_ini[post] - gaba_RE_RE[i][j][k].I/(re_cell[i][j].S_RE*1e3);
            }
          }
        }
      }
   }


  }  // End synapse

}  // End fun





//***************************************************************************
// Solver
//***************************************************************************

void rk(unsigned n, void fun(unsigned, double, double*, double*), 
        double h, double x, double* y, double* f, double* s, double* yk)
{
  int i;
  double xk;
  double h_half = h/2.;
  
  
  fun(n, x, y, f);
  
  for(i = n-1; i >=0; --i){
    s[i] = f[i];
    yk[i] = y[i] + h_half*f[i];
  }
  
  xk = x + h_half;
  fun(n, xk, yk, f);
  
  for(i = n-1; i >=0; --i){ 
    s[i] += 2.*f[i];
    yk[i] = y[i] + h_half*f[i];
  }
  
  fun(n, xk, yk, f);
  
  for(i = n-1; i >=0; --i){
    s[i] += 2.*f[i]; yk[i] = y[i] + h*f[i];
  }
  
  xk = x + h;
  fun(n, xk, yk, f);
  
  for(i = n-1; i >=0; --i){ 
    y[i] += (h/6.)*(s[i] + f[i]);
  }
}
//***************************************************************************


double gaussrand()
{
        static double U, V;
        static int phase = 0;
        double Z;

        if(phase == 0) {
                U = (rand() + 1.) / (RAND_MAX + 2.);
                V = rand() / (RAND_MAX + 1.);
                Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
        } else
                Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

        phase = 1 - phase;

        return Z;
}




