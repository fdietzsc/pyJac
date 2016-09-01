#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  // pressure dependence variable declarations
  double k0;
  double kinf;
  double Pr;

  // troe variable declarations
  double logFcent;
  double A;
  double B;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 4;
  pres_mod[0] = m + 1.5 * C[1] + 11.0 * C[5] - 0.17 * C[8];

  // reaction 5;
  pres_mod[1] = m + 1.5 * C[1] + 11.0 * C[5] - 0.17 * C[8];

  // reaction 6;
  pres_mod[2] = m + 1.5 * C[1] + 11.0 * C[5] - 0.25 * C[8];

  // reaction 7;
  pres_mod[3] = m + 1.5 * C[1] + 11.0 * C[5] - 0.62 * C[8];

  // reaction 8;
  thd = m + 0.351725400243 * C[1] - 0.35 * C[8] + 10.882917962 * C[5];
  k0 = exp(3.0395010469354865e+01 - 1.21156686 * logT);
  kinf = exp(2.2260133056545676e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(3.30000000e-01 * exp(-T / 1.00000000e-30) + 6.70000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 15;
  thd = m + 4.0 * C[5] - 0.2 * C[3] - 0.57 * C[8] + 4.13 * C[7] + 1.47 * C[1];
  k0 = exp(1.9441619737891969e+01 + 0.05525603 * logT - (-2.3574490000000001e+03 / T));
  kinf = exp(5.3697073626347098e+00 + 2.3219 * logT - (-1.7123099999999999e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

} // end get_rxn_pres_mod

