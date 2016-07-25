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
  printf("%s\n", buf);
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
  char bb_pe_fname_base[200];
  char bm_pe_fname_base[200];
  char t_pe_fname_base[200];
  char um_pe_fname_base[200];
  char bb_angle_fname_base[200];
  char bm_angle_fname_base[200];
  char t_angle_fname_base[200];
  char um_angle_fname_base[200];
  char config_fname[200];

  char bb_pe_bin_fname[200], bb_pe_plot_fname[200];
  char bm_pe_bin_fname[200], bm_pe_plot_fname[200];
  char t_pe_bin_fname[200],  t_pe_plot_fname[200];
  char um_pe_bin_fname[200], um_pe_plot_fname[200];

  char bb_angle_bin_fname[200], bb_angle_plot_fname[200];
  char bm_angle_bin_fname[200], bm_angle_plot_fname[200];
  char t_angle_bin_fname[200],  t_angle_plot_fname[200];
  char um_angle_bin_fname[200], um_angle_plot_fname[200];

  strcpy(bb_pe_fname_base, "data/bba_pe_data_");
  strcpy(bm_pe_fname_base, "data/bma_pe_data_");
  strcpy(t_pe_fname_base, "data/ta_pe_data_");
  strcpy(um_pe_fname_base, "data/uma_pe_data_");
  strcpy(bb_angle_fname_base, "data/bba_angle_data_");
  strcpy(bm_angle_fname_base, "data/bma_angle_data_");
  strcpy(t_angle_fname_base, "data/ta_angle_data_");
  strcpy(um_angle_fname_base, "data/uma_angle_data_");
  strcpy(config_fname, "data/config_");

  strcat(bb_pe_fname_base, f_appended_name);
  strcat(bm_pe_fname_base, f_appended_name);
  strcat(t_pe_fname_base, f_appended_name);
  strcat(um_pe_fname_base, f_appended_name);
  strcat(bb_angle_fname_base, f_appended_name);
  strcat(bm_angle_fname_base, f_appended_name);
  strcat(t_angle_fname_base, f_appended_name);
  strcat(um_angle_fname_base, f_appended_name);
  strcat(config_fname, f_appended_name);

  strcpy(bb_pe_bin_fname, bb_pe_fname_base); strcpy(bb_pe_plot_fname, bb_pe_fname_base);
  strcpy(bm_pe_bin_fname, bm_pe_fname_base); strcpy(bm_pe_plot_fname, bm_pe_fname_base);
  strcpy( t_pe_bin_fname,  t_pe_fname_base); strcpy( t_pe_plot_fname,  t_pe_fname_base);
  strcpy(um_pe_bin_fname, um_pe_fname_base); strcpy(um_pe_plot_fname, um_pe_fname_base);
  strcpy(bb_angle_bin_fname, bb_angle_fname_base); strcpy(bb_angle_plot_fname, bb_angle_fname_base);
  strcpy(bm_angle_bin_fname, bm_angle_fname_base); strcpy(bm_angle_plot_fname, bm_angle_fname_base);
  strcpy( t_angle_bin_fname,  t_angle_fname_base); strcpy( t_angle_plot_fname,  t_angle_fname_base);
  strcpy(um_angle_bin_fname, um_angle_fname_base); strcpy(um_angle_plot_fname, um_angle_fname_base);

  strcat(bb_pe_bin_fname, ".bin");
  strcat(bm_pe_bin_fname, ".bin");
  strcat(t_pe_bin_fname, ".bin");
  strcat(um_pe_bin_fname, ".bin");
  strcat(bb_angle_bin_fname, ".bin");
  strcat(bm_angle_bin_fname, ".bin");
  strcat(t_angle_bin_fname, ".bin");
  strcat(um_angle_bin_fname, ".bin");
  strcat(config_fname, ".txt");

  char* pe_fnames[] = {bb_pe_plot_fname, bm_pe_plot_fname, t_pe_plot_fname, um_pe_plot_fname};
  char* angle_fnames[] = {bb_angle_plot_fname, bm_angle_plot_fname,
			  t_angle_plot_fname, um_angle_plot_fname};
  const char* f_legends[] = {"Bound binding", "Bound motor", "Tail", "Unbound motor"};

  int bb_pe_fd = open(bb_pe_bin_fname, O_RDONLY);
  int bm_pe_fd = open(bm_pe_bin_fname, O_RDONLY);
  int  t_pe_fd = open( t_pe_bin_fname, O_RDONLY);
  int um_pe_fd = open(um_pe_bin_fname, O_RDONLY);
  int bb_angle_fd = open(bb_angle_bin_fname, O_RDONLY);
  int bm_angle_fd = open(bm_angle_bin_fname, O_RDONLY);
  int  t_angle_fd = open( t_angle_bin_fname, O_RDONLY);
  int um_angle_fd = open(um_angle_bin_fname, O_RDONLY);
  int pe_fds[] = {bb_pe_fd, bm_pe_fd, t_pe_fd, um_pe_fd};
  int angle_fds[] = {bb_angle_fd, bm_angle_fd, t_angle_fd, um_angle_fd};

  if (errno) {
    perror("Failed opening binary data files");
    exit(errno);
  }

  for (int i=0; i < 4; i++) {
    int fd = pe_fds[i];
    struct stat my_stat;
    fstat(fd, &my_stat);
    int len = my_stat.st_size / sizeof(data_2d);

    data_2d* data;
    data = (data_2d*) mmap(NULL, len*sizeof(data_2d), PROT_READ, MAP_PRIVATE, fd, 0);
    if (errno) {
      perror("Failed to mmap a data file");
      exit(errno);
    }

    double* time_data = (double*) malloc(len * sizeof(double));
    double* pe_data = (double*) malloc(len * sizeof(double));
    for (int j=0; j < len; j++) {
      time_data[j] = data[j].time;
      pe_data[j] = data[j].d;
    }

    generate_correlation_fn_data(pe_data, len, f_legends[i], pe_fnames[i]);
    generate_pe_vs_time_data(time_data, pe_data, len, f_legends[i], pe_fnames[i]);
    generate_ave_pe_and_log_error_data(time_data, pe_data, len, f_legends[i], pe_fnames[i]);
  }

  for (int i=0; i < 4; i++) {
    int fd = angle_fds[i];
    struct stat my_stat;
    fstat(fd, &my_stat);
    int len = my_stat.st_size / sizeof(data_2d);

    data_2d* data;
    data = (data_2d*) mmap(NULL, len*sizeof(data_2d), PROT_READ, MAP_PRIVATE, fd, 0);
    if (errno) {
      perror("Failed to mmap a data file");
      exit(errno);
    }

    double* time_data = (double*) malloc(len * sizeof(double));
    double* angle_data = (double*) malloc(len * sizeof(double));
    for (int j=0; j < len; j++) {
      time_data[j] = data[j].time;
      angle_data[j] = data[j].d;
    }

    double eq_angle;
    if (i == 0) eq_angle = onebound_post_powerstroke_internal_angles.bba;
    if (i == 1) eq_angle = onebound_post_powerstroke_internal_angles.bma;
    if (i == 2) eq_angle = onebound_post_powerstroke_internal_angles.ta;
    if (i == 3) eq_angle = onebound_post_powerstroke_internal_angles.uma;
    generate_angle_vs_time_data(time_data, angle_data, len, f_legends[i], angle_fnames[i], eq_angle);
  }
}
