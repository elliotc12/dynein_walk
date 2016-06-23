#include "dynein_struct.h"

#include <stdlib.h>
#include <fstream>

extern double bla_init;
extern double mla_init;
extern double mra_init;
extern double bra_init;

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1) * range;
}

double square(double num) {
	return num * num;
}

double cube(double num) {
	return num * num * num;
}

double fourth(double num) {
	return num * num * num * num;
}

double fifth(double num) {
	return num * num * num * num * num;
}

double get_average(double* data, int len) {
  double sum = 0;
  for (int i = 0; i < len; i++) {
    sum += data[i];
  }
  return sum / len;
}

// void log_run(FILE* run_file, double runtime, double run_length,
// 			      double distance_traveled, int steps) {

//   float ave_step_dist = distance_traveled / steps;
//   float ave_step_time = runtime / steps;

//   printf("\n\n***********Run data**********\n");
//   printf("Run length: %g nm\n", run_length);
//   printf("Distance traveled: %g nm\n", distance_traveled);
//   printf("Steps: %d\n", steps);
//   printf("Average step length: %g nm\n", ave_step_dist);
//   printf("Average step time: %g s\n\n\n", ave_step_time);
//   fprintf(run_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
//   fprintf(run_file, "%f\t%f\t%d\t%f\t%g\n",
//   	  run_length, distance_traveled, steps, ave_step_dist, ave_step_time);
// }

void detect_nans(Dynein_bothbound* dyn_bb, const char* loc) {
  if (dyn_bb->get_nma() != dyn_bb->get_nma()) printf("get_nma returns NaN at  %s\n", loc);
  if (dyn_bb->get_fma() != dyn_bb->get_fma()) printf("get_fma returns NaN at  %s\n", loc);

  if (dyn_bb->get_nba() != dyn_bb->get_nba()) printf("get_nba returns NaN at  %s\n", loc);
  if (dyn_bb->get_fba() != dyn_bb->get_fba()) printf("get_fba returns NaN at  %s\n", loc);

  if (dyn_bb->get_nbx() != dyn_bb->get_nbx()) printf("get_nbx returns NaN at  %s\n", loc);
  if (dyn_bb->get_nmx() != dyn_bb->get_nmx()) printf("get_nmx returns NaN at  %s\n", loc);
  if (dyn_bb->get_tx()  != dyn_bb->get_tx()) printf("get_tx returns NaN at  %s\n", loc);
  if (dyn_bb->get_fmx() != dyn_bb->get_fmx()) printf("get_fmx returns NaN at  %s\n", loc);
  if (dyn_bb->get_fbx() != dyn_bb->get_fbx()) printf("get_fbx returns NaN at  %s\n", loc);

  if (dyn_bb->get_nby() != dyn_bb->get_nby()) printf("get_nby returns NaN at  %s\n", loc);
  if (dyn_bb->get_nmy() != dyn_bb->get_nmy()) printf("get_nmy returns NaN at  %s\n", loc);
  if (dyn_bb->get_ty()  != dyn_bb->get_ty()) printf("get_ty returns NaN at  %s\n", loc);
  if (dyn_bb->get_fmy() != dyn_bb->get_fmy()) printf("get_fmy returns NaN at  %s\n", loc);
  if (dyn_bb->get_fby() != dyn_bb->get_fby()) printf("get_fby returns NaN at  %s\n", loc);

  if (dyn_bb->get_d_nma() != dyn_bb->get_d_nma()) printf("get_d_nma returns NaN at  %s\n", loc);
  if (dyn_bb->get_d_fma() != dyn_bb->get_d_fma()) printf("get_d_fma returns NaN at  %s\n", loc);

  // if (dyn_bb->get_d_nba() != dyn_bb->get_d_nba()) printf("get_d_nba returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_fba() != dyn_bb->get_d_fba()) printf("get_d_fba returns NaN at  %s\n", loc);

  // if (dyn_bb->get_d_nbx() != dyn_bb->get_d_nbx()) printf("get_d_nbx returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_nmx() != dyn_bb->get_d_nmx()) printf("get_d_nmx returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_tx()  != dyn_bb->get_d_tx()) printf("get_d_tx returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_fmx() != dyn_bb->get_d_fmx()) printf("get_d_fmx returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_fbx() != dyn_bb->get_d_fbx()) printf("get_d_fbx returns NaN at  %s\n", loc);

  // if (dyn_bb->get_d_nby() != dyn_bb->get_d_nby()) printf("get_d_nby returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_nmy() != dyn_bb->get_d_nmy()) printf("get_d_nmy returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_ty()  != dyn_bb->get_d_ty()) printf("get_d_ty returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_fmy() != dyn_bb->get_d_fmy()) printf("get_d_fmy returns NaN at  %s\n", loc);
  // if (dyn_bb->get_d_fby() != dyn_bb->get_d_fby()) printf("get_d_fby returns NaN at  %s\n", loc);
}
