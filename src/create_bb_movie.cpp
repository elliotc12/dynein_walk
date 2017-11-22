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
  double* nba;
  double* nma;
  double* ta;
  double* fma;
  double* fba;
  double* nba_PE;
  double* nma_PE;
  double* ta_PE;
  double* fma_PE;
  double* fba_PE;
  double* nbx;     double* nby;
  double* nmx;     double* nmy;
  double* tx;      double* ty;
  double* fmx;     double* fmy;
  double* fbx;     double* fby;
  double* f_nbx;   double* f_nby;
  double* f_nmx;   double* f_nmy;
  double* f_tx;    double* f_ty;
  double* f_fmx;   double* f_fmy;
  double* f_fbx;   double* f_fby;
} movie_generate_struct_bothbound;

void generate_movie(double* time, movie_generate_struct_bothbound* data, int len, char* fname_base) {
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
	    BOTHBOUND,
	    time[iter],
	    data->nba_PE[iter], data->nma_PE[iter], data->ta_PE[iter], data->fma_PE[iter], data->fba_PE[iter],
	    data->nbx[iter], data->nby[iter], data->nmx[iter], data->nmy[iter],
	    data->tx[iter], data->ty[iter], data->fmx[iter], data->fmy[iter],
	    data->fbx[iter], data->fby[iter],
	    data->f_nbx[iter], data->f_nby[iter], data->f_nmx[iter], data->f_nmy[iter], data->f_tx[iter],
	    data->f_ty[iter], data->f_fmx[iter], data->f_fmy[iter], data->f_fbx[iter], data->f_fby[iter]);
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

  strcpy(data_fname, "data/bothbound_data_");
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
  int len = data_fd_stat.st_size / sizeof(bothbound_data_generate_struct);

  bothbound_data_generate_struct* data_map;
  data_map = (bothbound_data_generate_struct*) mmap(NULL, len*sizeof(bothbound_data_generate_struct), PROT_READ, MAP_PRIVATE, data_fd, 0);

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

  double* time = new double[len];
  for (int j=0; j < len; j++) { time[j] = data_map[j].time; }

  double* nba_pe = new double[len];
  double* nma_pe = new double[len];
  double* ta_pe = new double[len];
  double* fma_pe = new double[len];
  double* fba_pe = new double[len];

  double* nba_angle = new double[len];
  double* nma_angle = new double[len];
  double* ta_angle = new double[len];
  double* fma_angle = new double[len];

  double* nbx = new double[len];
  double* nmx = new double[len];
  double* tx  = new double[len];
  double* fmx = new double[len];
  double* fbx = new double[len];

  double* nby = new double[len];
  double* nmy = new double[len];
  double* ty  = new double[len];
  double* fmy = new double[len];
  double* fby = new double[len];

  double* f_nbx = new double[len];
  double* f_nmx = new double[len];
  double* f_tx  = new double[len];
  double* f_fmx = new double[len];
  double* f_fbx = new double[len];

  double* f_nby = new double[len];
  double* f_nmy = new double[len];
  double* f_ty  = new double[len];
  double* f_fmy = new double[len];
  double* f_fby = new double[len];

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
  }

  movie_generate_struct_bothbound movie_data;
  movie_data.nba_PE = nba_pe;
  movie_data.nma_PE = nma_pe;
  movie_data.ta_PE  = ta_pe;
  movie_data.fma_PE = fma_pe;
  movie_data.fba_PE = fba_pe;
  movie_data.nba = nba_angle;
  movie_data.nma = nma_angle;
  movie_data.ta  = ta_angle;
  movie_data.fma = fma_angle;
  movie_data.nbx = nbx;      movie_data.nby = nby;
  movie_data.nmx = nmx;      movie_data.nmy = nmy;
  movie_data.tx  = tx;	      movie_data.ty = ty;
  movie_data.fmx = fmx;      movie_data.fmy = fmy;
  movie_data.fbx = fbx;      movie_data.fby = fby;
  movie_data.f_nbx = f_nbx;
  movie_data.f_nmx = f_nmx;
  movie_data.f_tx  = f_tx;
  movie_data.f_fmx = f_fmx;
  movie_data.f_fbx = f_fbx;
  movie_data.f_nby = f_nby;
  movie_data.f_nmy = f_nmy;
  movie_data.f_ty  = f_ty;
  movie_data.f_fmy = f_fmy;
  movie_data.f_fby = f_fby;
 
  generate_movie(time, &movie_data, len, movie_fname_base);

  free(nba_pe); free(nma_pe); free(ta_pe); free(fma_pe);
  free(nba_angle); free(nma_angle); free(ta_angle); free(fma_angle);
  free(nbx); free(nmx); free(tx); free(fmx); free(fbx);
  free(nby); free(nmy); free(ty); free(fmy); free(fby);
  free(f_nbx); free(f_nmx); free(f_tx); free(f_fmx); free(f_fbx);
  free(f_nby); free(f_nmy); free(f_ty); free(f_fmy); free(f_fby);

  free(time);
}
