#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e8;

  int low_T  = 1000;
  int high_T = 10000;

  Lt = 15;
  Ls = 15;

  fake_radius_t = 1.5;
  fake_radius_m = 1.5;
  fake_radius_b = 1.5;

  gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
  gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
  gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

  const char* bba_eq_title = "Bound binding";
  const char* bma_eq_title = "Bound motor";
  const char* ta_eq_title =  "Tail domain";
  const char* uma_eq_title = "Unbound motor";

  const char* bba_eq_fname = "bba_pe_equipartition_ratio_vs_c.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio_vs_c.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio_vs_c.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio_vs_c.txt";

  const int seeds[] = {0, 1};
  int seed_len = sizeof(seeds) / sizeof(int);

  int spring_len = 100;
  double spring_constants[spring_len];
  int min_spring = 1;
  int max_spring = 650;
  int d_spring = (max_spring - min_spring) / spring_len;
  int spring = min_spring;
  int i = 0;

  while(spring < max_spring) {
    spring_constants[i] = spring;
    spring += d_spring;
    i++;
  }

  double eq_ratios_bb[spring_len];
  double eq_ratios_bm[spring_len];
  double eq_ratios_t[spring_len];
  double eq_ratios_um[spring_len];

  for (int i = 0; i < spring_len; i++) {
    ct = spring_constants[i];
    cm = spring_constants[i];
    cb = spring_constants[i];

    onebound_data low_eq_data;
    double low_T_bb, low_T_bm, low_T_t, low_T_um;
    low_eq_data.len = 1;
    low_eq_data.bb = &low_T_bb;
    low_eq_data.bm = &low_T_bm;
    low_eq_data.t =  &low_T_t;
    low_eq_data.um = &low_T_um;

    onebound_data high_eq_data;
    double high_T_bb, high_T_bm, high_T_t, high_T_um;
    high_eq_data.len = 1;
    high_eq_data.bb = &high_T_bb;
    high_eq_data.bm = &high_T_bm;
    high_eq_data.t =  &high_T_t;
    high_eq_data.um = &high_T_um;

    onebound_forces unused_forces;
    generic_data unused_data;
    unused_data.data = &unused_forces;    

    T = low_T;
    get_onebound_equipartition_ratio(&low_eq_data, &unused_data, iterations, seeds, seed_len);
    T = high_T;
    get_onebound_equipartition_ratio(&high_eq_data, &unused_data, iterations, seeds, seed_len);

    double low_ET = (0.5*kb*low_T);
    double high_ET = (0.5*kb*high_T);
    
    eq_ratios_bb[i] = (high_T_bb*high_ET - low_T_bb*low_ET) / (0.5*kb*(high_T-low_T));
    eq_ratios_bm[i] = (high_T_bm*high_ET - low_T_bm*low_ET) / (0.5*kb*(high_T-low_T));
    eq_ratios_t[i]  = (high_T_t* high_ET - low_T_t* low_ET) / (0.5*kb*(high_T-low_T));
    eq_ratios_um[i] = (high_T_um*high_ET - low_T_um*low_ET) / (0.5*kb*(high_T-low_T));
  }

  print_data_to_file(spring_constants, eq_ratios_bb, spring_len, bba_eq_title, bba_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_bm, spring_len, bma_eq_title, bma_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_t,  spring_len, ta_eq_title,  ta_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_um, spring_len, uma_eq_title, uma_eq_fname);
  
  return 0;
}
