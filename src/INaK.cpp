#include "INaK.h" 
#include <math.h>
#include "../Constants.h"

INaK::INaK(double v) {

	G_K = 10;///////////////////////
    G_Na = 100;/////////////////////
    Vtr = -40;
    VtrK = -40;
    v2 = v - Vtr;
    v2K = v - VtrK;

    k_h = 1;
    k_n = 1;

  //  Phi = 1;
    Phi = pow(3,((Cels-36)/10));
/*
    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18);
    Beta2 = 4/(exp((40 - v2)/5) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
    Beta3 = 0.5*exp((10 - v2K)/40);
    n0 = Alpha3/(Alpha3 + Beta3);
    */
}


double INaK::E_K = -80, INaK::E_Na = 50, INaK::Cels = 36;

void INaK::init(double v) {
    v2 = v - Vtr;
    v2K = v - VtrK;

    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18);
    Beta2 = 4/(exp((40 - v2)/5) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
    Beta3 = 0.5*exp((10 - v2K)/40);
    n0 = Alpha3/(Alpha3 + Beta3);

}

void INaK::calc(double m, double h, double n, double &fm, double &fh, double &fn,
                   double v, double x){
  v2 = v - Vtr;
  v2K = v - VtrK;
  iNa = G_Na*m*m*m*h*(v - ENA);
  Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
  Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
  tau_m = 1/(Alpha1 + Beta1) / Phi;
  m_inf = Alpha1/(Alpha1 + Beta1);

  Alpha2 = 0.128*exp((17 - v2)/18);
  Beta2 = 4/(exp((40 - v2)/5) + 1);
  tau_h = k_h*1/(Alpha2 + Beta2) / Phi;
  h_inf = Alpha2/(Alpha2 + Beta2);

  fm = -(m - m_inf)/tau_m;                 
  fh = -(h - h_inf)/tau_h;                 

  iK = G_K* n*n*n*n*(v - EK);
  Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
  Beta3 = 0.5*exp((10 - v2K)/40);
  tau_n = k_n*1/(Alpha3 + Beta3) / Phi;

  n_inf = Alpha3/(Alpha3 + Beta3);
  
  fn  = -(n - n_inf)/tau_n;                 
}


