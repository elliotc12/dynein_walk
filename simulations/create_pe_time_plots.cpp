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
  prepare_data_file(legend, fname);

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
  prepare_data_file(legend, fname);
  int pe_averaging_width = pe_averaging_width_fraction * iters;
  if (pe_averaging_width < 1) pe_averaging_width = 1;

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
  prepare_data_file(legend, fname_ave);
  prepare_data_file(legend, fname_err);

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
  char bb_fname_base[200];
  char bm_fname_base[200];
  char t_fname_base[200];
  char um_fname_base[200];
  char config_fname[200];

  char bb_bin_fname[200], bb_plot_fname[200];
  char bm_bin_fname[200], bm_plot_fname[200];
  char t_bin_fname[200],  t_plot_fname[200]; 
  char um_bin_fname[200], um_plot_fname[200];

  strcpy(bb_fname_base, "data/bba_pe_vs_time_");
  strcpy(bm_fname_base, "data/bma_pe_vs_time_");
  strcpy(t_fname_base, "data/ta_pe_vs_time_");
  strcpy(um_fname_base, "data/uma_pe_vs_time_");
  strcpy(config_fname, "data/config_pe_vs_time_");

  strcat(bb_fname_base, f_appended_name);
  strcat(bm_fname_base, f_appended_name);
  strcat(t_fname_base, f_appended_name);
  strcat(um_fname_base, f_appended_name);
  strcat(config_fname, f_appended_name);

  strcpy(bb_bin_fname, bb_fname_base); strcpy(bb_plot_fname, bb_fname_base);
  strcpy(bm_bin_fname, bm_fname_base); strcpy(bm_plot_fname, bm_fname_base);
  strcpy( t_bin_fname,  t_fname_base); strcpy( t_plot_fname,  t_fname_base);
  strcpy(um_bin_fname, um_fname_base); strcpy(um_plot_fname, um_fname_base);

  strcat(bb_bin_fname, ".bin");
  strcat(bm_bin_fname, ".bin");
  strcat(t_bin_fname, ".bin");
  strcat(um_bin_fname, ".bin");
  strcat(config_fname, ".txt");

  char* fnames[] = {bb_plot_fname, bm_plot_fname, t_plot_fname, um_plot_fname};
  const char* f_legends[] = {"Bound binding", "Bound motor", "Tail", "Unbound motor"};

  int bb_fd = open(bb_bin_fname, O_RDONLY);
  int bm_fd = open(bm_bin_fname, O_RDONLY);
  int  t_fd = open( t_bin_fname, O_RDONLY);
  int um_fd = open(um_bin_fname, O_RDONLY);
  int fds[] = {bb_fd, bm_fd, t_fd, um_fd};

  if (errno) {
    perror("Failed opening binary data files");
    exit(errno);
  }

  for (int i=0; i < 4; i++) {
    int fd = fds[i];
    struct stat my_stat;
    fstat(fd, &my_stat);
    int len = my_stat.st_size / sizeof(pe_data);

    pe_data* data;
    data = (pe_data*) mmap(NULL, len*sizeof(pe_data), PROT_READ, MAP_PRIVATE, fd, 0);
    if (errno) {
      perror("Failed to mmap a data file");
      exit(errno);
    }

    double* time_data = (double*) malloc(len * sizeof(double));
    double* pe_data = (double*) malloc(len * sizeof(double));
    for (int j=0; j < len; j++) {
      time_data[j] = data[j].time;
      pe_data[j] = data[j].pe;
    }

    generate_correlation_fn_data(pe_data, len, f_legends[i], fnames[i]);
    generate_pe_vs_time_data(time_data, pe_data, len, f_legends[i], fnames[i]);
    generate_ave_pe_and_log_error_data(time_data, pe_data, len, f_legends[i], fnames[i]);
  }
}
