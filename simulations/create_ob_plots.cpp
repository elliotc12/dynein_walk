#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"
#include "simulation_defaults.h"

void generate_correlation_fn_data(double* pe, int iters, const char* legend, char* fname_base){
  int max_tau_iter;
  if (iters*tau_runtime_fraction > num_corr_datapoints) {
    max_tau_iter = iters*tau_runtime_fraction;
  }
  else {
    max_tau_iter = num_corr_datapoints;
  }
  int d_tau_iter = max_tau_iter / num_corr_datapoints;

  char fname[200];
  strcpy(fname, fname_base);
  strcat(fname, "_correlation_fn.txt");

  char buf[256];
  sprintf(buf, "--legend='%s'", legend);
  prepare_data_file(buf, fname);

  double* taus =         (double*) malloc(num_corr_datapoints * sizeof(double));
  double* correlations = (double*) malloc(num_corr_datapoints * sizeof(double));

  double pe_ave = get_average(pe, iters);
  double pe_var = get_variance(pe, iters);

  for (int n=0; n<num_corr_datapoints; n++) {
    int tau_iter = n*d_tau_iter;
    double correlation_sum = 0;
    for (int i=0; i<(iters-tau_iter); i++) {
      correlation_sum += (pe[i] - pe_ave) * (pe[i+tau_iter] - pe_ave);
    }
    correlations[n] = correlation_sum / pe_var / (iters-tau_iter);
    taus[n] = tau_iter*dt*pe_calculation_skip_iterations;
    printf("Progress on %s: %.1f%%  \r", fname, n * 100.0 / num_corr_datapoints);
  }
  printf("Finished %s                \n", fname);

  FILE* file = fopen(fname, "a");
  append_data_to_file(taus, correlations, num_corr_datapoints, file);
  fclose(file);

  free (taus); free(correlations);
}

void generate_pe_vs_time_data(double* times, double* pe, int iters, const char* legend, char* fname_base) {
  char fname[200];
  strcpy(fname, fname_base);
  strcat(fname, ".txt");

  char buf[256];
  sprintf(buf, "--legend='%s", legend);
  prepare_data_file(buf, fname);

  double et = 0.5*kb*T;
  double* pe_local_ave = (double*) malloc(iters * sizeof(double));
  for (int i = 0; i < iters; i++) {
    if (i == 0 or i == iters-1) {
      pe_local_ave[i] = pe[i] / et;
    }
    else if (i < pe_averaging_width/2) {
      pe_local_ave[i] = get_average(pe, i*2) / et;
    }
    else if (i == pe_averaging_width/2) {
      pe_local_ave[i] = get_average(&pe[i-pe_averaging_width/2], pe_averaging_width) / et;
    }
    else if ( i > pe_averaging_width/2 and i <= (iters-pe_averaging_width/2-1)){
      pe_local_ave[i] = (pe_local_ave[i-1]*i + pe[i]/et) / (i+1);
    }
    else if (i > (iters-pe_averaging_width/2-1)) {
      pe_local_ave[i] = get_average(&pe[iters-1-(iters-1-i)*2], (iters-1-i)*2) / et;
    }
    else {
      printf("Error in PE local averaging!\n");
      exit(1);
    }
    printf("Progress for %s: %.1f%%  \r", fname, i * 100.0 / iters);
  }
  printf("Finished %s                        \n", fname);

  FILE* file = fopen(fname, "a");
  append_data_to_file(times, pe_local_ave, iters, file);
  fclose(file);

  free(pe_local_ave);
}

void generate_ave_pe_and_log_error_data(double* times, double* pe, int iters, const char* legend, char* fname_base) {
  char fname_ave[200];
  char fname_err[200];

  strcpy(fname_ave, fname_base);
  strcpy(fname_err, fname_base);

  strcat(fname_ave, "_eq_ave.txt");
  strcat(fname_err, "_log_error.txt");

  char buf[256];
  sprintf(buf, "--legend='%s'", legend);
  prepare_data_file(buf, fname_ave);
  sprintf(buf, "--legend='%s'", legend);
  prepare_data_file(buf, fname_err);

  double* pe_aves = (double*) malloc(iters * sizeof(double));
  double* log_err = (double*) malloc(iters * sizeof(double));

  double et = 0.5*kb*T;
  pe_aves[0] = pe[0]/et;
  for (int i=1; i<iters; i++) {
    pe_aves[i] = (pe_aves[i-1]*i + pe[i]/et) / (i+1);
    printf("Progress for %s: %.1f%%  \r", fname_ave, i * 100.0 / iters);
  }
  printf("Finished %s                        \n", fname_ave);

  for (int i=0; i<iters; i++) {
    log_err[i] = std::abs(pe_aves[i] - 1);
    printf("Progress for %s: %.1f%%  \r", fname_err, i * 100.0 / iters);
  }
  printf("Finished %s                       \n", fname_err);

  FILE* file_ave = fopen(fname_ave, "a");
  FILE* file_err = fopen(fname_err, "a");
  append_data_to_file(times, pe_aves, iters, file_ave);
  append_data_to_file(times, log_err, iters, file_err);
  fclose(file_ave);
  fclose(file_err);

  free(pe_aves); free(log_err);
}

void generate_angle_vs_time_data(double* times, double* angle, int iters, const char* legend, char* fname_base, double eq_angle) {
  char fname[200];
  strcpy(fname, fname_base);
  strcat(fname, ".txt");

  char buf[256];
  sprintf(buf, "--legend='%s', --hline='%g'", legend, eq_angle);
  prepare_data_file(buf, fname);

  double* angle_local_ave = (double*) malloc(iters * sizeof(double));
  for (int i = 0; i < iters; i++) {
    if (i == 0 or i == iters-1) {
      angle_local_ave[i] = angle[i];
    }
    else if (i < angle_averaging_width/2) {
      angle_local_ave[i] = get_average(angle, i*2);
    }
    else if (i == angle_averaging_width/2) {
      angle_local_ave[i] = get_average(&angle[i-angle_averaging_width/2], angle_averaging_width);
    }
    else if ( i > angle_averaging_width/2 and i <= (iters-angle_averaging_width/2-1)){
      angle_local_ave[i] = (angle_local_ave[i-1]*i + angle[i]) / (i+1);
    }
    else if (i > (iters-angle_averaging_width/2-1)) {
      angle_local_ave[i] = get_average(&angle[iters-1-(iters-1-i)*2], (iters-1-i)*2);
    }
    else {
      printf("Error in angle local averaging!\n");
      exit(1);
    }
    printf("Progress for %s: %.1f%%  \r", fname, i * 100.0 / iters);
  }
  printf("Finished %s                        \n", fname);

  FILE* file = fopen(fname, "a");
  append_data_to_file(times, angle_local_ave, iters, file);
  fclose(file);

  free(angle_local_ave);
}

int main(int argc, char** argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char data_fname[200];
  char bba_pe_fname_base[200];
  char bma_pe_fname_base[200];
  char  ta_pe_fname_base[200];
  char uma_pe_fname_base[200];
  char bba_fname_base[200];
  char bma_fname_base[200];
  char  ta_fname_base[200];
  char uma_fname_base[200];

  strcpy(data_fname, "data/onebound_data_");
  strcat(data_fname, f_appended_name);
  strcat(data_fname, ".bin");
  int data_fd = open(data_fname, O_RDONLY);

  if (errno) {
    perror("Failed opening binary data file");
    exit(errno);
  }

  strcat(bba_pe_fname_base, "data/ob_bba_pe_"); strcat(bba_pe_fname_base, f_appended_name);
  strcat(bma_pe_fname_base, "data/ob_bma_pe_"); strcat(bma_pe_fname_base, f_appended_name);
  strcat( ta_pe_fname_base, "data/ob_ta_pe_");  strcat( ta_pe_fname_base, f_appended_name);
  strcat(uma_pe_fname_base, "data/ob_uma_pe_"); strcat(uma_pe_fname_base, f_appended_name);

  strcat(bba_fname_base, "data/ob_bba_angle_"); strcat(bba_fname_base, f_appended_name);
  strcat(bma_fname_base, "data/ob_bma_angle_"); strcat(bma_fname_base, f_appended_name);
  strcat( ta_fname_base, "data/ob_ta_angle_");  strcat( ta_fname_base, f_appended_name);
  strcat(uma_fname_base, "data/ob_uma_angle_"); strcat(uma_fname_base, f_appended_name);

  struct stat data_fd_stat;
  fstat(data_fd, &data_fd_stat);
  int len = data_fd_stat.st_size / sizeof(onebound_data_generate_struct);

  onebound_data_generate_struct* data_map;
  data_map = (onebound_data_generate_struct*) mmap(NULL, len*sizeof(onebound_data_generate_struct), PROT_READ, MAP_PRIVATE, data_fd, 0);

  if (data_map == MAP_FAILED) {
    perror("Error using mmap: ");
    exit(EXIT_FAILURE);
  }

  double* time = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { time[j] = data_map[j].time; }

  /** Bound binding PE **/
  double* bba_pe = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { bba_pe[j] = data_map[j].bba_PE; }
  generate_correlation_fn_data(bba_pe, len, "Bound binding", bba_pe_fname_base);
  generate_pe_vs_time_data(time, bba_pe, len, "Bound binding", bba_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, bba_pe, len, "Bound binding", bba_pe_fname_base);
  free(bba_pe);

  /** Bound motor PE **/
  double* bma_pe = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { bma_pe[j] = data_map[j].bma_PE; }
  generate_correlation_fn_data(bma_pe, len, "Bound motor", bma_pe_fname_base);
  generate_pe_vs_time_data(time, bma_pe, len, "Bound motor", bma_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, bma_pe, len, "Bound motor", bma_pe_fname_base);
  free(bma_pe);

  /** Tail PE **/
  double* ta_pe = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { ta_pe[j] = data_map[j].ta_PE; }
  generate_correlation_fn_data(ta_pe, len, "Tail", ta_pe_fname_base);
  generate_pe_vs_time_data(time, ta_pe, len, "Tail", ta_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, ta_pe, len, "Tail", ta_pe_fname_base);
  free(ta_pe);

  /** Unbound motor PE **/
  double* uma_pe = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { uma_pe[j] = data_map[j].uma_PE; }
  generate_correlation_fn_data(uma_pe, len, "Unbound motor", uma_pe_fname_base);
  generate_pe_vs_time_data(time, uma_pe, len, "Unbound motor", uma_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, uma_pe, len, "Unbound motor", uma_pe_fname_base);
  free(uma_pe);

  /** Bound binding angle **/
  double* bba_angle = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { bba_angle[j] = data_map[j].bba; }
  generate_angle_vs_time_data(time, bba_angle, len, "Bound binding", bba_fname_base, onebound_post_powerstroke_internal_angles.bba);
  free(bba_angle);

  /** Bound motor angle **/
  double* bma_angle = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { bma_angle[j] = data_map[j].bma; }
  generate_angle_vs_time_data(time, bma_angle, len, "Bound motor", bma_fname_base, onebound_post_powerstroke_internal_angles.bma);
  free(bma_angle);

  /** Tail angle **/
  double* ta_angle = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { ta_angle[j] = data_map[j].ta; }
  generate_angle_vs_time_data(time, ta_angle, len, "Tail", ta_fname_base, onebound_post_powerstroke_internal_angles.ta);
  free(ta_angle);

  /** Unbound motor angle **/
  double* uma_angle = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { uma_angle[j] = data_map[j].uma; }
  generate_angle_vs_time_data(time, uma_angle, len, "Unbound motor", uma_fname_base, onebound_post_powerstroke_internal_angles.uma);
  free(uma_angle);
  free(time);
}
