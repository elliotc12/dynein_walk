#include <stdio.h>
#include <string.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e5;
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

  char run_msg[512];
  const char* run_msg_base = "tempcalc (";

  prepare_data_file(bba_eq_title, bba_eq_fname);
  prepare_data_file(bma_eq_title, bma_eq_fname);
  prepare_data_file(ta_eq_title,  ta_eq_fname);
  prepare_data_file(uma_eq_title, uma_eq_fname);

  for (int i = 0; i < num_Ts; i++) {
    T = temps[i];

    onebound_data eq_data;
    double bb_double, bm_double, t_double, um_double;
    eq_data.len = 1;
    eq_data.bb = &bb_double;
    eq_data.bm = &bm_double;
    eq_data.t =  &t_double;
    eq_data.um = &um_double;

    onebound_forces force_var_data_struct;
    generic_data unused_force_var_data;
    unused_force_var_data.len = 1;
    unused_force_var_data.data = &force_var_data_struct;

    strcpy(run_msg, run_msg_base);
    char buf[50];
    sprintf(buf, "temp = %.1f, ", T);
    strcat(run_msg, buf);

    get_onebound_equipartition_ratio(
	  &eq_data, &unused_force_var_data, iterations, seeds, seed_len, run_msg);

    append_data_to_file(&T, eq_data.bb, 1, bba_eq_fname);
    append_data_to_file(&T, eq_data.bm, 1, bma_eq_fname);
    append_data_to_file(&T, eq_data.t , 1, ta_eq_fname);
    append_data_to_file(&T, eq_data.um, 1, uma_eq_fname);
  }

  return 0;
}
