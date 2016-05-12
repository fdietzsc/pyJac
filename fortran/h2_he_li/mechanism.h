#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 11
0  H/* Species Indexes

1  H2/* Species Indexes

2  O/* Species Indexes

3  OH/* Species Indexes

4  H2O/* Species Indexes

5  CO/* Species Indexes

6  HCO/* Species Indexes

7  O2/* Species Indexes

8  HO2/* Species Indexes

9  H2O2/* Species Indexes

10  CO2/* Species Indexes

11  HE/* Species Indexes

12  N2*/

//Number of species
#define NSP 13
//Number of variables. NN = NSP + 1 (temperature)
#define NN 14
//Number of forward reactions
#define FWD_RATES 23
//Number of reversible reactions
#define REV_RATES 23
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 6

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

