#include <cassert>
#include <limits>

#include "../default_parameters.h"
#include "../dynein_struct.h"

int main() {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  int iterations = 1e7;

  onebound_data low_T_PE_ratios;
  onebound_data high_T_PE_ratios;
  double low_T_bb_PE_ratio, low_T_bm_PE_ratio, low_T_t_PE_ratio, low_T_um_PE_ratio;
  double high_T_bb_PE_ratio, high_T_bm_PE_ratio, high_T_t_PE_ratio, high_T_um_PE_ratio;
  
  low_T_PE_ratios.bb = &low_T_bb_PE_ratio;
  low_T_PE_ratios.bm = &low_T_bm_PE_ratio;
  low_T_PE_ratios.t = &low_T_t_PE_ratio;
  low_T_PE_ratios.um = &low_T_um_PE_ratio;
  low_T_PE_ratios.len = 1;

  high_T_PE_ratios.bb = &high_T_bb_PE_ratio;
  high_T_PE_ratios.bm = &high_T_bm_PE_ratio;
  high_T_PE_ratios.t = &high_T_t_PE_ratio;
  high_T_PE_ratios.um = &high_T_um_PE_ratio;
  high_T_PE_ratios.len = 1;

  generic_data time_data;
  double t_data;
  time_data.data = &t_data;
  time_data.len = 1;

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

  // ct = 0.5;
  // cm = 0.5;
  // cb = 0.5;

  //const int seeds[] = {0, 1, 2, 3, 4};
  //int seed_len = sizeof(seeds) / sizeof(int);
  
  T = low_T;
  get_onebound_equipartition_ratio_average_per_runtime
    (&time_data, &low_T_PE_ratios, 1, iterations, iterations+1);
  
  T = high_T;
  get_onebound_equipartition_ratio_average_per_runtime
    (&time_data, &high_T_PE_ratios, 1, iterations, iterations+1);

  double d_bba_PE = high_T_bb_PE_ratio*(0.5*kb*high_T) - low_T_bb_PE_ratio*(0.5*kb*low_T);
  double d_bma_PE = high_T_bm_PE_ratio*(0.5*kb*high_T) - low_T_bm_PE_ratio*(0.5*kb*low_T);
  double d_ta_PE =  high_T_t_PE_ratio*(0.5*kb*high_T) - low_T_t_PE_ratio*(0.5*kb*low_T);
  double d_uma_PE = high_T_um_PE_ratio*(0.5*kb*high_T) - low_T_um_PE_ratio*(0.5*kb*low_T);
  double d_T = high_T - low_T;

  printf("bba PE/0.5KBT ratio: %g\n", d_bba_PE / (0.5*kb*d_T));
  printf("bma PE/0.5KBT ratio: %g\n", d_bma_PE / (0.5*kb*d_T));
  printf("ta PE/0.5KBT ratio: %g\n", d_ta_PE / (0.5*kb*d_T));
  printf("uma PE/0.5KBT ratio: %g\n", d_uma_PE / (0.5*kb*d_T));
  printf("total PE/0.5KBT ratio: %g\n",
	 (d_bba_PE + d_bma_PE + d_ta_PE + d_uma_PE) / (4*0.5*kb*d_T));
  
  return EXIT_SUCCESS;
}
