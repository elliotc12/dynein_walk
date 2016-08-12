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
    taus[n] = tau_iter*dt*data_generation_skip_iterations;
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

  double* pe_local_ave = (double*) malloc(num_generate_pe_datapoints * sizeof(double));
  double* sampled_times = (double*) malloc(num_generate_angle_datapoints * sizeof(double));
  int iters_per_i = iters / num_generate_pe_datapoints;
  for (int i = 0; i < num_generate_pe_datapoints; i++) {
    int iter = i*iters_per_i;
    if (iter == 0 or iter == iters-1) {
      pe_local_ave[i] = pe[iter];
    }
    else if (iter < generate_averaging_width/2) {
      pe_local_ave[i] = get_average(pe, iter*2);
    }
    else if ( iter >= generate_averaging_width/2 and iter <= (iters-generate_averaging_width/2-1)){
      pe_local_ave[i] = get_average(&pe[iter-generate_averaging_width/2],
				    generate_averaging_width);
    }
    else if (iter > (iters-generate_averaging_width/2-1)) {
      pe_local_ave[i] = get_average(&pe[iters-1-(iters-1-iter)*2], (iters-1-iter)*2);
    }
    else {
      printf("Error in PE local averaging!\n");
      exit(1);
    }
    sampled_times[i] = times[iter];
    printf("Progress for %s: %.1f%%  \r", fname, iter * 100.0 / iters);
  }
  printf("Finished %s                        \n", fname);

  FILE* file = fopen(fname, "a");
  append_data_to_file(sampled_times, pe_local_ave, num_generate_pe_datapoints, file);
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

  pe_aves[0] = pe[0];
  for (int i=1; i<iters; i++) {
    pe_aves[i] = (pe_aves[i-1]*i + pe[i]) / (i+1);
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

  double* angle_local_ave = (double*) malloc(num_generate_angle_datapoints * sizeof(double));
  double* sampled_times = (double*) malloc(num_generate_angle_datapoints * sizeof(double));
  int iters_per_i = iters / num_generate_angle_datapoints;
  for (int i = 0; i < num_generate_angle_datapoints; i++) {
    int iter = i*iters_per_i;
    if (iter == 0 or iter == iters-1) {
      angle_local_ave[i] = angle[iter];
    }
    else if (iter < generate_averaging_width/2) {
      angle_local_ave[i] = get_average(angle, iter*2);
    }
    else if (iter >= generate_averaging_width/2 and iter <= (iters-generate_averaging_width/2-1)) {
      angle_local_ave[i] = get_average(&angle[iter-generate_averaging_width/2],
				       generate_averaging_width);
    }
    else if (iter > (iters-generate_averaging_width/2-1)) {
      angle_local_ave[i] = get_average(&angle[iters-1-(iters-1-iter)*2], (iters-1-iter)*2);
    }
    else {
      printf("Error in angle local averaging!\n");
      exit(1);
    }
    sampled_times[i] = times[iter];
    printf("Progress for %s: %.1f%%  \r", fname, iter * 100.0 / iters);
  }
  printf("Finished %s                        \n", fname);

  FILE* file = fopen(fname, "a");
  append_data_to_file(sampled_times, angle_local_ave, num_generate_angle_datapoints, file);
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
  char* data_fname = new char[200];
  char* nba_pe_fname_base = new char[200];
  char* nma_pe_fname_base = new char[200];
  char*  ta_pe_fname_base = new char[200];
  char* fma_pe_fname_base = new char[200];
  char* fba_pe_fname_base = new char[200];
  char* total_pe_fname_base = new char[200];
  char* nba_fname_base = new char[200];
  char* nma_fname_base = new char[200];
  char*  ta_fname_base = new char[200];
  char* fma_fname_base = new char[200];
  char* fba_fname_base = new char[200];
  char* everything_name = new char[200];

  strcpy(data_fname, "data/bothbound_data_");
  strcat(data_fname, f_appended_name);
  strcat(data_fname, ".bin");
  int data_fd = open(data_fname, O_RDONLY);

  if (errno) {
    perror("Failed opening binary data file");
    printf("File name: %s\n", data_fname);
    exit(errno);
  }

  sprintf(everything_name, "data/everything_%s.txt", f_appended_name);

  strcpy(nba_pe_fname_base, "data/bb_nba_pe_"); strcat(nba_pe_fname_base, f_appended_name);
  strcpy(nma_pe_fname_base, "data/bb_nma_pe_"); strcat(nma_pe_fname_base, f_appended_name);
  strcpy( ta_pe_fname_base, "data/bb_ta_pe_");  strcat( ta_pe_fname_base, f_appended_name);
  strcpy(fma_pe_fname_base, "data/bb_fma_pe_"); strcat(fma_pe_fname_base, f_appended_name);
  strcpy(fba_pe_fname_base, "data/bb_fba_pe_"); strcat(fba_pe_fname_base, f_appended_name);
  strcpy(total_pe_fname_base, "data/bb_total_pe_"); strcat(total_pe_fname_base, f_appended_name);

  strcpy(nba_fname_base, "data/bb_nba_angle_"); strcat(nba_fname_base, f_appended_name);
  strcpy(nma_fname_base, "data/bb_nma_angle_"); strcat(nma_fname_base, f_appended_name);
  strcpy( ta_fname_base, "data/bb_ta_angle_");  strcat( ta_fname_base, f_appended_name);
  strcpy(fma_fname_base, "data/bb_fma_angle_"); strcat(fma_fname_base, f_appended_name);
  strcpy(fba_fname_base, "data/bb_fba_angle_"); strcat(fba_fname_base, f_appended_name);

  struct stat data_fd_stat;
  fstat(data_fd, &data_fd_stat);
  int len = data_fd_stat.st_size / sizeof(bothbound_data_generate_struct);

  bothbound_data_generate_struct* data_map;
  data_map = (bothbound_data_generate_struct*) mmap(NULL, len*sizeof(bothbound_data_generate_struct), PROT_READ, MAP_SHARED, data_fd, 0);
  if (data_map == MAP_FAILED) {
    perror("Error using mmap: ");
    exit(EXIT_FAILURE);
  }

  for (int i=0; i<len; i++) { // set the length to data size, in case simulation not complete
    if (data_map[i].nbx == 0 && data_map[i].nmx == 0 && data_map[i].tx == 0 && data_map[i].fmx == 0 && data_map[i].fbx == 0) {
      len = i;
      break;
    }
  }

  {
    FILE *f = fopen(everything_name, "w");
    fprintf(f, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            "nba_PE",
            "nma_PE",
            "ta_PE",
            "fma_PE",
	    "fba_PE",
            "nba",
            "nma",
            "ta",
            "fma",
	    "fba",
            "nbx",
            "nmx",
            "tx",
            "fmx",
            "fbx",
            "nby",
            "nmy",
            "ty",
            "fmy",
            "fby",
            "f_nbx",
            "f_nmx",
            "f_tx",
            "f_fmx",
            "f_fbx",
            "f_nby",
            "f_nmy",
            "f_ty",
            "f_fmy",
            "f_fby");
    for (int j=0;j<len;j++) {
      fprintf(f, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
              data_map[j].nba_PE,
              data_map[j].nma_PE,
              data_map[j].ta_PE,
              data_map[j].fma_PE,
	      data_map[j].fba_PE,
              data_map[j].nba,
              data_map[j].nma,
              data_map[j].ta,
              data_map[j].fma,
	      data_map[j].fba,
              data_map[j].nbx,
              data_map[j].nmx,
              data_map[j].tx,
              data_map[j].fmx,
              data_map[j].fbx,
              data_map[j].nby,
              data_map[j].nmy,
              data_map[j].ty,
              data_map[j].fmy,
              data_map[j].fby,
              data_map[j].f_nbx,
              data_map[j].f_nmx,
              data_map[j].f_tx,
              data_map[j].f_fmx,
              data_map[j].f_fbx,
              data_map[j].f_nby,
              data_map[j].f_nmy,
              data_map[j].f_ty,
              data_map[j].f_fmy,
              data_map[j].f_fby);
    }
    fclose(f);
  }

  double* time = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { time[j] = data_map[j].time; }

  double* nba_pe = (double*) malloc(len * sizeof(double));
  double* nma_pe = (double*) malloc(len * sizeof(double));
  double* ta_pe = (double*) malloc(len * sizeof(double));
  double* fma_pe = (double*) malloc(len * sizeof(double));
  double* fba_pe = (double*) malloc(len * sizeof(double));
  double* total_pe = (double*) malloc(len * sizeof(double));

  double* nba_angle = (double*) malloc(len * sizeof(double));
  double* nma_angle = (double*) malloc(len * sizeof(double));
  double* ta_angle = (double*) malloc(len * sizeof(double));
  double* fma_angle = (double*) malloc(len * sizeof(double));
  double* fba_angle = (double*) malloc(len * sizeof(double));

  double* nbx = (double*) malloc(len * sizeof(double));
  double* nmx = (double*) malloc(len * sizeof(double));
  double* tx  = (double*) malloc(len * sizeof(double));
  double* fmx = (double*) malloc(len * sizeof(double));
  double* fbx = (double*) malloc(len * sizeof(double));

  double* nby = (double*) malloc(len * sizeof(double));
  double* nmy = (double*) malloc(len * sizeof(double));
  double* ty  = (double*) malloc(len * sizeof(double));
  double* fmy = (double*) malloc(len * sizeof(double));
  double* fby = (double*) malloc(len * sizeof(double));

  double* f_nbx = (double*) malloc(len * sizeof(double));
  double* f_nmx = (double*) malloc(len * sizeof(double));
  double* f_tx  = (double*) malloc(len * sizeof(double));
  double* f_fmx = (double*) malloc(len * sizeof(double));
  double* f_fbx = (double*) malloc(len * sizeof(double));

  double* f_nby = (double*) malloc(len * sizeof(double));
  double* f_nmy = (double*) malloc(len * sizeof(double));
  double* f_ty  = (double*) malloc(len * sizeof(double));
  double* f_fmy = (double*) malloc(len * sizeof(double));
  double* f_fby = (double*) malloc(len * sizeof(double));

  for (int j=0; j < len; j++) {
    nba_pe[j] = data_map[j].nba_PE;
    nma_pe[j] = data_map[j].nma_PE;
     ta_pe[j] = data_map[j].ta_PE;
    fma_pe[j] = data_map[j].fma_PE;
    fba_pe[j] = data_map[j].fba_PE;
    nba_angle[j] = data_map[j].nba;
    nma_angle[j] = data_map[j].nma;
     ta_angle[j] = data_map[j].ta;
    fma_angle[j] = data_map[j].fma;
    fba_angle[j] = data_map[j].fba;
    nbx[j] = data_map[j].nbx;
    nmx[j] = data_map[j].nmx;
    tx [j] = data_map[j].tx;
    fmx[j] = data_map[j].fmx;
    fbx[j] = data_map[j].fbx;
    nby[j] = data_map[j].nby;
    nmy[j] = data_map[j].nmy;
    ty [j] = data_map[j].ty;
    fmy[j] = data_map[j].fmy;
    fby[j] = data_map[j].fby;
    f_nbx[j] = data_map[j].f_nbx;
    f_nmx[j] = data_map[j].f_nmx;
    f_tx[j]  = data_map[j].f_tx;
    f_fmx[j] = data_map[j].f_fmx;
    f_fbx[j] = data_map[j].f_fbx;
    f_nby[j] = data_map[j].f_nby;
    f_nmy[j] = data_map[j].f_nmy;
    f_ty[j]  = data_map[j].f_ty;
    f_fmy[j] = data_map[j].f_fmy;
    f_fby[j] = data_map[j].f_fby;
    total_pe[j] = nba_pe[j] + nma_pe[j] + ta_pe[j] + fma_pe[j] + fba_pe[j];
  }

  generate_correlation_fn_data(nba_pe, len, "Near binding", nba_pe_fname_base);
  generate_correlation_fn_data(nma_pe, len, "Near motor", nma_pe_fname_base);
  generate_correlation_fn_data(ta_pe, len,  "Tail", ta_pe_fname_base);
  generate_correlation_fn_data(fma_pe, len, "Far motor", fma_pe_fname_base);
  generate_correlation_fn_data(fba_pe, len, "Far binding", fba_pe_fname_base);

  generate_pe_vs_time_data(time, nba_pe, len, "Near binding PE", nba_pe_fname_base);
  generate_pe_vs_time_data(time, nma_pe, len, "Near motor PE", nma_pe_fname_base);
  generate_pe_vs_time_data(time, ta_pe, len, "Tail PE", ta_pe_fname_base);
  generate_pe_vs_time_data(time, fma_pe, len, "Far motor PE", fma_pe_fname_base);
  generate_pe_vs_time_data(time, fba_pe,len,"Far binding PE", fba_pe_fname_base);
  generate_pe_vs_time_data(time, total_pe, len, "Total PE", total_pe_fname_base);

  generate_ave_pe_and_log_error_data(time, nba_pe, len, "Near binding", nba_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, nma_pe, len, "Near motor", nma_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, ta_pe, len, "Tail", ta_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, fma_pe, len, "Far motor", fma_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, fba_pe, len, "Far binding", fba_pe_fname_base);
  generate_ave_pe_and_log_error_data(time, total_pe, len, "Total", total_pe_fname_base);

  generate_angle_vs_time_data(time, nba_angle, len, "Near binding angle", nba_fname_base, bothbound_pre_powerstroke_internal_angles.nba);
  generate_angle_vs_time_data(time, nma_angle, len, "Near motor angle", nma_fname_base, bothbound_pre_powerstroke_internal_angles.nma);
  generate_angle_vs_time_data(time, ta_angle, len, "Tail angle", ta_fname_base, bothbound_pre_powerstroke_internal_angles.ta);
  generate_angle_vs_time_data(time, fma_angle, len, "Far motor angle", fma_fname_base, bothbound_pre_powerstroke_internal_angles.fma);
  generate_angle_vs_time_data(time, fba_angle, len, "Far binding angle", fba_fname_base, bothbound_pre_powerstroke_internal_angles.fba);

  free(nba_pe); free(nma_pe); free(ta_pe); free(fma_pe);
  free(nba_angle); free(nma_angle); free(ta_angle); free(fma_angle);
  free(nbx); free(nmx); free(tx); free(fmx); free(fbx);
  free(nby); free(nmy); free(ty); free(fmy); free(fby);
  free(f_nbx); free(f_nmx); free(f_tx); free(f_fmx); free(f_fbx);
  free(f_nby); free(f_nmy); free(f_ty); free(f_fmy); free(f_fby);

  free(time);
}
