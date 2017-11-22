#ifndef PLOTTING_DEFAULTS_H
#define PLOTTING_DEFAULTS_H

const int num_corr_datapoints = 1000;
const double tau_runtime_fraction = 1e-5;
const int movie_num_frames = 1e2;

const int num_generate_pe_datapoints = 50;
const int num_generate_angle_datapoints = num_generate_pe_datapoints;
const int num_generate_force_datapoints = num_generate_pe_datapoints;

const int custom_generate_averaging_width = 1; // 0 for full-dataset averaging

#include "custom_simulation_parameters.h"

#endif
