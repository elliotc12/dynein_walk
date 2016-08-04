#pragma once

const long long iterations = 1e12;
const long long data_generation_skip_iterations = 1e8;

const int num_corr_datapoints = 1000;
const double tau_runtime_fraction = 1e-5;
const int movie_num_frames = 1e2;

const int num_generate_pe_datapoints = 250;
const int num_generate_angle_datapoints = num_generate_pe_datapoints;
const int generate_averaging_width = iterations / data_generation_skip_iterations /
  num_generate_pe_datapoints;

const int msync_after_num_writes = 10;
