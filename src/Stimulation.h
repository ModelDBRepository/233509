

class STIM {

   int k;
   int ON;
   int FLAG_F;
   int FLAG_I;
   int N_Pulse;

   double F;          // Current frequency
   double Pulse_T0;
   double Pulse_T00;
   double T_Pulse;    // Total pulse duration (ON + OFF)
   double Step_T0;    // Initial time of each step


 public:
   int    FLAG_VAR_F;
   double Stim_T0;      // Starting time of stimulation
   double Stim_Dur;
   double Stim_Amp;
   double F0;            // Initial frequency
   double T_PULSE_ON;   // Pulse ON duration
   double T_Step;       // Duration of each frequency stimulation
   double STEP_F;
   double I;

   STIM();
   void init();
   void calc(double x);
};
