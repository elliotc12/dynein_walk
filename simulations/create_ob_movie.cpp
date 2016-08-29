#include <cassert>
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
#include "plotting_defaults.h"

typedef struct {
  double* bba;
  double* bma;
  double* ta;
  double* uma;
  double* bba_PE;
  double* bma_PE;
  double* ta_PE;
  double* uma_PE;
  double* bbx;     double* bby;
  double* bmx;     double* bmy;
  double* tx;      double* ty;
  double* umx;     double* umy;
  double* ubx;     double* uby;
  double* f_bbx;   double* f_bby;
  double* f_bmx;   double* f_bmy;
  double* f_tx;    double* f_ty;
  double* f_umx;   double* f_umy;
  double* f_ubx;   double* f_uby;
} movie_generate_struct;

void generate_movie(double* time, movie_generate_struct* data, int len, char* fname_base) {
  char fname[200];
  strcpy(fname, fname_base);
  strcat(fname, ".txt");

  FILE* data_file = fopen(fname, "w");
  
  int d_iter = len / movie_num_frames;
  if (d_iter < 1) d_iter = 1;

  for (int iter = 0; iter < len; iter += d_iter) {
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t"
	    "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	    "\n",
	    NEARBOUND,
	    time[iter],
	    data->bba_PE[iter], data->bma_PE[iter], data->ta_PE[iter], data->uma_PE[iter], 0.0,
	    data->bbx[iter], data->bby[iter], data->bmx[iter], data->bmy[iter],
	    data->tx[iter], data->ty[iter], data->umx[iter], data->umy[iter],
	    data->ubx[iter], data->uby[iter],
	    data->f_bbx[iter], data->f_bby[iter], data->f_bmx[iter], data->f_bmy[iter], data->f_tx[iter],
	    data->f_ty[iter], data->f_umx[iter], data->f_umy[iter], data->f_ubx[iter], data->f_uby[iter]);
    printf("Progress on %s: %.1f%%  \r", fname, iter * 100.0 / len);
  }
  printf("Finished %s                \n", fname);
  fclose(data_file);
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
  char movie_fname_base[200];

  strcpy(data_fname, "data/onebound_data_");
  strcat(data_fname, f_appended_name);
  strcat(data_fname, ".bin");
  int data_fd = open(data_fname, O_RDONLY);

  if (errno) {
    perror("Failed opening binary data file");
    printf("File name: %s\n", data_fname);
    exit(errno);
  }

  strcpy(movie_fname_base, "data/movie_"); strcat(movie_fname_base, f_appended_name);

  struct stat data_fd_stat;
  fstat(data_fd, &data_fd_stat);
  int len = data_fd_stat.st_size / sizeof(onebound_data_generate_struct);

  onebound_data_generate_struct* data_map;
  data_map = (onebound_data_generate_struct*) mmap(NULL, len*sizeof(onebound_data_generate_struct), PROT_READ, MAP_PRIVATE, data_fd, 0);

  if (data_map == MAP_FAILED) {
    perror("Error using mmap: ");
    exit(EXIT_FAILURE);
  }

  for (int i=0; i<len; i++) { // set the length to data size, in case simulation not complete
    if (data_map[i].bbx == 0 && data_map[i].bmx == 0 && data_map[i].tx == 0 && data_map[i].umx == 0 && data_map[i].ubx == 0) {
      len = i;
      break;
    }
  }

  double* time = (double*) malloc(len * sizeof(double));
  for (int j=0; j < len; j++) { time[j] = data_map[j].time; }

  double* bba_pe = (double*) malloc(len * sizeof(double));
  double* bma_pe = (double*) malloc(len * sizeof(double));
  double* ta_pe = (double*) malloc(len * sizeof(double));
  double* uma_pe = (double*) malloc(len * sizeof(double));

  double* bba_angle = (double*) malloc(len * sizeof(double));
  double* bma_angle = (double*) malloc(len * sizeof(double));
  double* ta_angle = (double*) malloc(len * sizeof(double));
  double* uma_angle = (double*) malloc(len * sizeof(double));

  double* bbx = (double*) malloc(len * sizeof(double));
  double* bmx = (double*) malloc(len * sizeof(double));
  double* tx  = (double*) malloc(len * sizeof(double));
  double* umx = (double*) malloc(len * sizeof(double));
  double* ubx = (double*) malloc(len * sizeof(double));

  double* bby = (double*) malloc(len * sizeof(double));
  double* bmy = (double*) malloc(len * sizeof(double));
  double* ty  = (double*) malloc(len * sizeof(double));
  double* umy = (double*) malloc(len * sizeof(double));
  double* uby = (double*) malloc(len * sizeof(double));

  double* f_bbx = (double*) malloc(len * sizeof(double));
  double* f_bmx = (double*) malloc(len * sizeof(double));
  double* f_tx  = (double*) malloc(len * sizeof(double));
  double* f_umx = (double*) malloc(len * sizeof(double));
  double* f_ubx = (double*) malloc(len * sizeof(double));

  double* f_bby = (double*) malloc(len * sizeof(double));
  double* f_bmy = (double*) malloc(len * sizeof(double));
  double* f_ty  = (double*) malloc(len * sizeof(double));
  double* f_umy = (double*) malloc(len * sizeof(double));
  double* f_uby = (double*) malloc(len * sizeof(double));

  for (int j=0; j < len; j++) {
    bba_pe[j] = data_map[j].bba_PE;
    bma_pe[j] = data_map[j].bma_PE;
     ta_pe[j] = data_map[j].ta_PE;
    uma_pe[j] = data_map[j].uma_PE;
    bba_angle[j] = data_map[j].bba;
    bma_angle[j] = data_map[j].bma;
     ta_angle[j] = data_map[j].ta;
    uma_angle[j] = data_map[j].uma;
    bbx[j] = data_map[j].bbx;
    bmx[j] = data_map[j].bmx;
    tx [j] = data_map[j].tx;
    umx[j] = data_map[j].umx;
    ubx[j] = data_map[j].ubx;
    bby[j] = data_map[j].bby;
    bmy[j] = data_map[j].bmy;
    ty [j] = data_map[j].ty;
    umy[j] = data_map[j].umy;
    uby[j] = data_map[j].uby;
    f_bbx[j] = data_map[j].f_bbx;
    f_bmx[j] = data_map[j].f_bmx;
    f_tx[j]  = data_map[j].f_tx;
    f_umx[j] = data_map[j].f_umx;
    f_ubx[j] = data_map[j].f_ubx;
    f_bby[j] = data_map[j].f_bby;
    f_bmy[j] = data_map[j].f_bmy;
    f_ty[j]  = data_map[j].f_ty;
    f_umy[j] = data_map[j].f_umy;
    f_uby[j] = data_map[j].f_uby;  
  }

  movie_generate_struct movie_data;
  movie_data.bba_PE = bba_pe;
  movie_data.bma_PE = bma_pe;
  movie_data.ta_PE  = ta_pe;
  movie_data.uma_PE = uma_pe;
  movie_data.bba = bba_angle;
  movie_data.bma = bma_angle;
  movie_data.ta  = ta_angle;
  movie_data.uma = uma_angle;
  movie_data.bbx = bbx;      movie_data.bby = bby;
  movie_data.bmx = bmx;      movie_data.bmy = bmy;
  movie_data.tx  = tx;	      movie_data.ty = ty;
  movie_data.umx = umx;      movie_data.umy = umy;
  movie_data.ubx = ubx;      movie_data.uby = uby;
  movie_data.f_bbx = f_bbx;
  movie_data.f_bmx = f_bmx;
  movie_data.f_tx  = f_tx;
  movie_data.f_umx = f_umx;
  movie_data.f_ubx = f_ubx;
  movie_data.f_bby = f_bby;
  movie_data.f_bmy = f_bmy;
  movie_data.f_ty  = f_ty;
  movie_data.f_umy = f_umy;
  movie_data.f_uby = f_uby;
 
  generate_movie(time, &movie_data, len, movie_fname_base);

  free(bba_pe); free(bma_pe); free(ta_pe); free(uma_pe);
  free(bba_angle); free(bma_angle); free(ta_angle); free(uma_angle);
  free(bbx); free(bmx); free(tx); free(umx); free(ubx);
  free(bby); free(bmy); free(ty); free(umy); free(uby);
  free(f_bbx); free(f_bmx); free(f_tx); free(f_umx); free(f_ubx);
  free(f_bby); free(f_bmy); free(f_ty); free(f_umy); free(f_uby);

  free(time);
}
