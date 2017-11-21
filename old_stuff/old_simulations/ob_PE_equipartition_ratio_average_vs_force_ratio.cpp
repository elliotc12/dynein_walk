#include <cassert>
#include <stdio.h>
#include <string.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main(int argc, char** argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  
  int iterations = 1e5;
  T = 1000;

  const char* bba_eq_title = "Bound binding domain (bba)";
  const char* bma_eq_title = "Bound motor domain (bma)";
  const char* ta_eq_title =  "Tail domain (ta)";
  const char* uma_eq_title = "Unbound motor domain (uma)";

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }
  char* f_appended_name = argv[1];
  char bba_eq_fname[200];
  char bma_eq_fname[200];
  char ta_eq_fname[200];
  char uma_eq_fname[200];
  char config_eq_fname[200];

  strcpy(bba_eq_fname, "data/bba_pe_equipartition_ratio_vs_f_ratio_");
  strcpy(bma_eq_fname, "data/bma_pe_equipartition_ratio_vs_f_ratio_");
  strcpy(ta_eq_fname,  "data/ta_pe_equipartition_ratio_vs_f_ratio_");
  strcpy(uma_eq_fname, "data/uma_pe_equipartition_ratio_vs_f_ratio_");
  strcpy(config_eq_fname, "data/config_pe_equipartition_ratio_vs_f_ratio_");

  strcat(bba_eq_fname, f_appended_name);
  strcat(bma_eq_fname, f_appended_name);
  strcat(ta_eq_fname, f_appended_name);
  strcat(uma_eq_fname, f_appended_name);
  strcat(config_eq_fname, f_appended_name);

  strcat(bba_eq_fname, ".txt");
  strcat(bma_eq_fname, ".txt");
  strcat(ta_eq_fname, ".txt");
  strcat(uma_eq_fname, ".txt");
  strcat(config_eq_fname, ".txt");

  prepare_data_file(bba_eq_title, bba_eq_fname);
  prepare_data_file(bma_eq_title, bma_eq_fname);
  prepare_data_file(ta_eq_title,  ta_eq_fname);
  prepare_data_file(uma_eq_title, uma_eq_fname);

  const int seeds[] = {0, 1};
  int seed_len = sizeof(seeds) / sizeof(int);

  double spring_constants[] = {50};

  int num_gammas = 50;
  double gamma_modifiers[num_gammas];
  double min_gamma = 0.1;
  double max_gamma = 5.0;
  double d_gamma = (max_gamma - min_gamma) / num_gammas;
  double gamma = min_gamma;
  int i = 0;

  while(gamma < max_gamma) {
    gamma_modifiers[i] = gamma;
    gamma += d_gamma;
    i++;
  }

  int num_springs = sizeof(spring_constants) / sizeof(double);

  char run_msg[512];
  const char* run_msg_base = "force variance calc (";

  for (int i = 0; i < num_springs; i++) {
    ct = spring_constants[i];
    cm = spring_constants[i];
    cb = spring_constants[i];

    for (int j = 0; j < num_gammas; j++) {
      gm = fake_radius_m*6*M_PI*water_viscosity_mu*gamma_modifiers[j];
      gb = fake_radius_b*6*M_PI*water_viscosity_mu*gamma_modifiers[j];
      gt = fake_radius_t*6*M_PI*water_viscosity_mu*gamma_modifiers[j];

      onebound_data eq_data;
      double bb_double, bm_double, t_double, um_double;
      eq_data.len = 1;
      eq_data.bb = &bb_double;
      eq_data.bm = &bm_double;
      eq_data.t =  &t_double;
      eq_data.um = &um_double;

      onebound_forces force_var_data_struct;
      generic_data force_var_data;
      force_var_data.len = 1;
      force_var_data.data = &force_var_data_struct;

      strcpy(run_msg, run_msg_base);
      char buf[50];
      sprintf(buf, "c = %.1f, g_mod = %.1f ", spring_constants[i], gamma_modifiers[j]);
      strcat(run_msg, buf);

      get_onebound_equipartition_ratio(
		       &eq_data, &force_var_data, iterations, seeds, seed_len, run_msg);

      double bb_force_var_ave = (force_var_data_struct.bbx + force_var_data_struct.bby) / 2;
      double bm_force_var_ave = (force_var_data_struct.bmx + force_var_data_struct.bmy) / 2;
      double t_force_var_ave = (force_var_data_struct.tx + force_var_data_struct.ty) / 2;
      double um_force_var_ave = (force_var_data_struct.umx + force_var_data_struct.umy) / 2;

      double brownian_b_variance = gb*sqrt(2*kb*T/(gb*dt)); // should this val be sqrt'd?? var or std?
      double brownian_m_variance = gm*sqrt(2*kb*T/(gm*dt));
      double brownian_t_variance = gt*sqrt(2*kb*T/(gt*dt));

      double force_ratios_bb = brownian_b_variance / bb_force_var_ave;
      double force_ratios_bm = brownian_m_variance / bm_force_var_ave;
      double force_ratios_t  = brownian_t_variance / t_force_var_ave; 
      double force_ratios_um = brownian_m_variance / um_force_var_ave;

      append_data_to_file(eq_data.bb, &force_ratios_bb, 1, bba_eq_fname);
      append_data_to_file(eq_data.bm, &force_ratios_bm, 1, bma_eq_fname);
      append_data_to_file(eq_data.t , &force_ratios_t,  1, ta_eq_fname);
      append_data_to_file(eq_data.um, &force_ratios_um, 1, uma_eq_fname);
    }
  }

  write_config_file(config_eq_fname, NULL, NULL);
  return 0;
}
