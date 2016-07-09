#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e7;
  
  T = 1000;
  
  Lt = 15;
  Ls = 15;
  
  fake_radius_t = 1.5;
  fake_radius_m = 1.5;
  fake_radius_b = 1.5;
  gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
  gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
  gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s
  
  const char* bba_eq_title = "Bound binding domain (bba)";
  const char* bma_eq_title = "Bound motor domain (bma)";
  const char* ta_eq_title =  "Tail domain (ta)";
  const char* uma_eq_title = "Unbound motor domain (uma)";

  const char* bba_eq_fname = "bba_pe_equipartition_ratio_vs_c.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio_vs_c.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio_vs_c.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio_vs_c.txt";

  const int seeds[] = {0, 1};
  int seed_len = sizeof(seeds) / sizeof(int);

  double spring_constants[] = {1.0, 250.0, 500.0, 750.0, 1000.0};
  double gammas[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  int num_datapoints = sizeof(spring_constants) / sizeof(double);

  double eq_ratios_bb[num_springs];
  double eq_ratios_bm[num_springs];
  double eq_ratios_t[num_springs];
  double eq_ratios_um[num_springs];

  for (int i = 0; i < num_datapoints; i++) {
    ct = spring_constants[i];
    cm = spring_constants[i];
    cb = spring_constants[i];

    onebound_data eq_data;
    eq_data.len = 1;
    eq_data.bb = &eq_ratios_bb[i];
    eq_data.bm = &eq_ratios_bm[i];
    eq_data.t =  &eq_ratios_t[i];
    eq_data.um = &eq_ratios_um[i];

    
    
    get_onebound_equipartition_ratio_and_average_force(NULL, &eq_data, 1, iterations,
	iterations+1, seeds, seed_len);
  }

  print_data_to_file(spring_constants, eq_ratios_bb, num_springs, bba_eq_title, bba_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_bm, num_springs, bma_eq_title, bma_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_t,  num_springs, ta_eq_title,  ta_eq_fname);
  print_data_to_file(spring_constants, eq_ratios_um, num_springs, uma_eq_title, uma_eq_fname);
  
  return 0;
}
