#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e5;
  
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

  const char* bba_eq_fname = "bba_pe_equipartition_ratio_vs_f_ratio.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio_vs_f_ratio.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio_vs_f_ratio.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio_vs_f_ratio.txt";

  const int seeds[] = {0};
  int seed_len = sizeof(seeds) / sizeof(int);

  double spring_constants[] = {50, 60, 70, 80, 90, 100, 200, 300, 400};
  double gamma_modifiers[] = {0.5, 0.7, 1.0, 1.3, 1.5, 2.0, 2.5};
  int num_springs = sizeof(spring_constants) / sizeof(double);
  int num_gammas = sizeof(gamma_modifiers) / sizeof(double);

  double eq_ratios_bb[num_springs*num_gammas];
  double eq_ratios_bm[num_springs*num_gammas];
  double eq_ratios_t [num_springs*num_gammas];
  double eq_ratios_um[num_springs*num_gammas];

  double force_ratios_bb[num_springs*num_gammas];
  double force_ratios_bm[num_springs*num_gammas];
  double force_ratios_t [num_springs*num_gammas];
  double force_ratios_um[num_springs*num_gammas];

  for (int i = 0; i < num_springs; i++) {
    ct = spring_constants[i];
    cm = spring_constants[i];
    cb = spring_constants[i];

    for (int j = 0; j < num_gammas; j++) {
      gm = fake_radius_m*6*M_PI*water_viscosity_mu*gamma_modifiers[j];
      gb = fake_radius_b*6*M_PI*water_viscosity_mu*gamma_modifiers[j];
      gt = fake_radius_t*6*M_PI*water_viscosity_mu*gamma_modifiers[j];

      onebound_data eq_data;
      eq_data.len = 1;
      eq_data.bb = &eq_ratios_bb[i];
      eq_data.bm = &eq_ratios_bm[i];
      eq_data.t =  &eq_ratios_t[i];
      eq_data.um = &eq_ratios_um[i];

      onebound_forces force_var_data_struct;
      generic_data force_var_data;
      force_var_data.len = 1;
      force_var_data.data = &force_var_data_struct;

      get_onebound_equipartition_ratio(&eq_data, &force_var_data, iterations, seeds, seed_len);

      eq_ratios_bb[i*num_gammas + j] = eq_data.bb[0];
      eq_ratios_bm[i*num_gammas + j] = eq_data.bm[0];
      eq_ratios_t [i*num_gammas + j] = eq_data.t[0];
      eq_ratios_um[i*num_gammas + j] = eq_data.um[0];

      double bb_force_var_ave = (force_var_data_struct.bbx + force_var_data_struct.bby) / 2;
      double bm_force_var_ave = (force_var_data_struct.bmx + force_var_data_struct.bmy) / 2;
      double t_force_var_ave = (force_var_data_struct.tx + force_var_data_struct.ty) / 2;
      double um_force_var_ave = (force_var_data_struct.umx + force_var_data_struct.umy) / 2;

      double brownian_b_variance = sqrt(gb*sqrt(2*kb*T/(gb*dt))); // should this be sqrt?? var or std?
      double brownian_m_variance = sqrt(gm*sqrt(2*kb*T/(gm*dt)));
      double brownian_t_variance = sqrt(gt*sqrt(2*kb*T/(gt*dt)));

      force_ratios_bb[i*num_gammas + j] = brownian_b_variance / bb_force_var_ave;
      force_ratios_bm[i*num_gammas + j] = brownian_m_variance / bm_force_var_ave;
      force_ratios_t [i*num_gammas + j] = brownian_t_variance / t_force_var_ave; 
      force_ratios_um[i*num_gammas + j] = brownian_m_variance / um_force_var_ave;
    }
  }

  print_data_to_file(force_ratios_bb,eq_ratios_bb, num_springs*num_gammas, bba_eq_title, bba_eq_fname);
  print_data_to_file(force_ratios_bm,eq_ratios_bm, num_springs*num_gammas, bma_eq_title, bma_eq_fname);
  print_data_to_file(force_ratios_t, eq_ratios_t,  num_springs*num_gammas, ta_eq_title,  ta_eq_fname);
  print_data_to_file(force_ratios_um,eq_ratios_um, num_springs*num_gammas, uma_eq_title, uma_eq_fname);

  return 0;
}
