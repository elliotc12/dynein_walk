#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e8;
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

  const char* bba_eq_fname = "bba_pe_equipartition_ratio_vs_T.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio_vs_T.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio_vs_T.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio_vs_T.txt";

  const int seeds[] = {0, 1};
  int seed_len = sizeof(seeds) / sizeof(int);

  int num_Ts = 100;
  double temps[num_Ts];
  double min_temp = 150;
  double max_temp = 2000.0;
  double d_temp = (max_temp - min_temp) / num_Ts;
  double temp = min_temp;
  int i = 0;

  while(temp < max_temp) {
    temps[i] = temp;
    temp += d_temp;
    i++;
  }

  double eq_ratios_bb[num_Ts];
  double eq_ratios_bm[num_Ts];
  double eq_ratios_t [num_Ts];
  double eq_ratios_um[num_Ts];

  for (int i = 0; i < num_Ts; i++) {
    T = temps[i];

    onebound_data eq_data;
    eq_data.len = 1;
    eq_data.bb = &eq_ratios_bb[i];
    eq_data.bm = &eq_ratios_bm[i];
    eq_data.t =  &eq_ratios_t[i];
    eq_data.um = &eq_ratios_um[i];

    onebound_forces force_var_data_struct;
    generic_data unused_force_var_data;
    unused_force_var_data.len = 1;
    unused_force_var_data.data = &force_var_data_struct;

    get_onebound_equipartition_ratio(&eq_data, &unused_force_var_data, iterations, seeds, seed_len);

    eq_ratios_bb[i] = eq_data.bb[0];
    eq_ratios_bm[i] = eq_data.bm[0];
    eq_ratios_t [i] = eq_data.t[0];
    eq_ratios_um[i] = eq_data.um[0];
  }

  print_data_to_file(temps, eq_ratios_bb, num_Ts, bba_eq_title, bba_eq_fname);
  print_data_to_file(temps, eq_ratios_bm, num_Ts, bma_eq_title, bma_eq_fname);
  print_data_to_file(temps, eq_ratios_t,  num_Ts, ta_eq_title,  ta_eq_fname);
  print_data_to_file(temps, eq_ratios_um, num_Ts, uma_eq_title, uma_eq_fname);

  return 0;
}
