

#define OSC 1   // 1: Delta; 2: Spindle; 3: Alpha; 4: Gamma

// Simulation time
#define	t_end  3000.0     // Delta: 3000 ms; Alpha and Gamma: 2000 ms; Spindle: 4000 ms
                          // When stimulation is applied, the simulation time is 52000 ms if 1-50 Hz is simulated


// Stimulation Input
#define FLAG_STIMULATION_LGN_PULSE  0    // 0: No stimulation; 1: Pulsatile input to LGN
#define FLAG_STIMULATION_TRN_PULSE  0    // 0: No stimulation; 1: Pulsatile input to TRN

#define STIM_START_TIME     1000         // Stimulation start time
#define STIM_DURATION       1000000      // The maximal stimulation time
#define STIM_STEP_DURATION  1000         // Stimulation step duration

#define PULSE_ON_DURATION   10           // Pulse duration
#define STIM_MAX_F          50           // Maximal stimulation frequency

#define STIM_F0             0            // Initial frequency
#define STIM_F_INCREMENT    1            // Stimulation increment step
#define STIM_AMP            0.1          // Stimulation amplitude in nA


// Spindle-triggering inputs
#define SPIN_T0   1000    // Onset time of spindle-triggering input
#define SPIN_TT    100    // Spindle-triggering input duration
#define SPIN_F     0.1    // Spindle triggering input frequency

#define SPIN_AMP  0.1    // Spindle triggering input amplitude
#define SPIN_X    10     // Subset of RE_X receiving spindle-triggering input
#define SPIN_Y    10     // Subset of RE_Y receiving spindle-triggering input

// Random seed
#define RAND_SEED 100
#define PI   3.14159

// Current injection
#define FLAG_CURRENT_INJECTION 0

// Random Background Inputs
#define FLAG_RANDOM_INPUT  1
#define RANDOM_F_TC  0.1     // Random input frequency; 0.1 per ms (100 Hz)

#define N_INPUT_TC1  1      // Number of random background input
#define N_INPUT_TC2  1


// Reversal potential
#define ENA           50
#define EK           -90
#define E_GABA_RE    -70

// Short-term depression parameters
#define TAU_STD      700
#define F_USE_SPIKE  0.07


// Connection ON/OFF
#define FLAG_TC1_TC1_GAP 1
#define FLAG_TC2_TC1_GAP 1
#define FLAG_RE_RE_GAP   1

#define FLAG_TC1_IN  1
#define FLAG_IN_TC2  1

#define FLAG_RE_IN   1
#define FLAG_RE_RE   1

#define FLAG_TC1_RE  1
#define FLAG_TC2_RE  1

#define FLAG_RE_TC1  1
#define FLAG_RE_TC2  1


#define FLAG_TC2_TC2 0    // The RTC-RTC connections are not active !


// Connection probability
#define P_TC1_TC1_GAP  0.3
#define P_TC2_TC1_GAP  0.3
#define P_RE_RE_GAP    0.3

#define P_TC2_GAP      0.2 // Percentage of RTCs that have gap junctions with HTCs
#define P_RE_GAP       0.2 // Percentage of REs that have gap junctions with other REs

#define P_TC1_IN   0.3
#define P_IN_TC2   0.3

#define P_TC1_RE   0.2
#define P_TC2_RE   0.2
#define P_RE_TC1   0.2
#define P_RE_TC2   0.2

#define P_RE_RE    0.2
#define P_RE_IN    0.05

#define P_TC2_TC2  0.2

// Turn STD ON/OFF
#define STD_TC1_IN   1
#define STD_IN_TC2   1
#define STD_RE_IN    1

#define STD_RE_RE    1

#define STD_TC1_RE   1
#define STD_TC2_RE   1

#define STD_RE_TC1   1
#define STD_RE_TC2   1

#define STD_TC2_TC2  0


// Gap junctional resistance
#define R_TC1_TC1 100   // in MOhm
#define R_TC1_TC2 300   // in MOhm
#define R_RE_RE   300   // in MOhm

#define D_TC1_TC1  2   // Radius of the gap junctions
#define D_TC2_TC1  2
#define D_TC2_TC2  2
#define D_RE_RE    2


// Maximal conductance
#define g_AMPA_TC1_IN 0.006    // (in uS)
#define g_NMDA_TC1_IN 0.003    // (in uS)

#define g_AMPA_TC1_RE 0.004
#define g_NMDA_TC1_RE 0.002

#define g_AMPA_TC2_RE 0.004
#define g_NMDA_TC2_RE 0.002

#define g_AMPA_TC2_TC2 0.002
#define g_NMDA_TC2_TC2 0.001

// GABA
#define g_GABA_RE_TC2 0.003
#define g_GABA_RE_TC1 0.003

#define g_GABA_IN_TC2 0.003

#define g_GABA_RE_IN  0.001
#define g_GABA_RE_RE  0.001



//====================================================================
//       Leakage conductance at different states
//====================================================================

// HTC
#define GL1_TC1  0.0075     // lower bound of the regular leak
#define GL2_TC1  0.0125     // upper bound of the regular leak

#define GKL_TC1_S1 0.035     // potassium leak in Delta       (Low ACh/NE)
#define GKL_TC1_S2 0.01      // potassium leak in Spindle     (Medium ACh/NE)
#define GKL_TC1_S3 0.0       // potassium leak in Alpha/Gamma (High ACh/NE)

// RTC
#define GL1_TC2  0.0075
#define GL2_TC2  0.0125

#define GKL_TC2_S1 0.035   // potassium leak in Delta
#define GKL_TC2_S2 0.01    // potassium leak in Spindle
#define GKL_TC2_S3 0.0     // potassium leak in Alpha/Gamma

// IN
#define GL1_IN  0.0075
#define GL2_IN  0.0125

#define GKL_IN_S1 0.01     // potassium leak in Delta
#define GKL_IN_S2 0.015    // potassium leak in Spindle
#define GKL_IN_S3 0.02     // potassium leak in Alpha/Gamma

// RE
#define GL1_RE  0.0075
#define GL2_RE  0.0125

#define GKL_RE_S1 0.03     // potassium leak in Delta
#define GKL_RE_S2 0.02     // potassium leak in Spindle
#define GKL_RE_S3 0.01     // potassium leak in Alpha/Gamma


//=============================================================
//            Cell numbers
//=============================================================
#define TC1_EX 1
#define TC2_EX 1
#define IN_EX  1
#define RE_EX  1

// HTC cells
#define N_TC1_X 7   // Cell number in the x direction
#define N_TC1_Y 7   // Cell number in the y direction
#define N_TC1 (N_TC1_X*N_TC1_Y)  // Total number of HTC cells

// RTC cells
#define N_TC2_X 12
#define N_TC2_Y 12
#define N_TC2 (N_TC2_X*N_TC2_Y)

// Interneurons
#define N_IN_X  8
#define N_IN_Y  8
#define N_IN (N_IN_X*N_IN_Y)

// RE cells
#define N_RE_X  10
#define N_RE_Y  10
#define N_RE (N_RE_X*N_RE_Y)


/*******************************************************************
//                Number of equations
*******************************************************************/

// Number of ODEs for each cell type
#define EQ_TC 14
#define EQ_IN 10
#define EQ_RE 9

#define N_EQ_TC1  (EQ_TC * N_TC1)
#define N_EQ_TC2  (EQ_TC * N_TC2)
#define N_EQ_IN   (EQ_IN * N_IN)
#define N_EQ_RE   (EQ_RE * N_RE)


#define N_EQ   N_EQ_TC1 + N_EQ_TC2 + N_EQ_IN + N_EQ_RE





