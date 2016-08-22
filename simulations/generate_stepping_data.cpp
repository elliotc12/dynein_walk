#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulation_defaults.h"

const int INIT_DYNARR_LEN = 10;

typedef struct {
  int num_steps;
  DynArr* step_times;
  DynArr* step_lengths;
  double dwell_time;
} stepping_data_struct;

void stepping_data_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  static State last_state = BOTHBOUND;
  static double last_nbx = ((Dynein_bothbound*) dyn)->get_nbx();
  static double last_fbx = ((Dynein_bothbound*) dyn)->get_fbx();
  static double last_bothbound_iteration = 0;

  long long max_iteration = *((long long**) job_msg)[0];
  double start_time = *((double**) job_msg)[1];
  char* run_msg = ((char**) job_msg)[2];

  stepping_data_struct* data_struct = ((stepping_data_struct**) job_msg)[3];

  if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    if (last_state == NEARBOUND) {
      double step_len = dyn_bb->get_fbx() - last_fbx;
      data_struct->step_lengths->append(step_len);
      last_fbx = dyn_bb->get_fbx();

      double step_time = (iteration - last_bothbound_iteration)*dt;
      data_struct->step_times->append(step_time);

      data_struct->num_steps++;
      printf("Switched to BB!\n");
    }
    else if (last_state == FARBOUND) {
      double step_len = dyn_bb->get_nbx() - last_nbx;
      data_struct->step_lengths->append(step_len);
      last_nbx = dyn_bb->get_nbx();

      double step_time = (iteration - last_bothbound_iteration)*dt;
      data_struct->step_times->append(step_time);

      data_struct->num_steps++;
      printf("Switched to BB!\n");
    }
    last_bothbound_iteration = iteration;
    last_state = BOTHBOUND;
  }
  else if (s == NEARBOUND) {
    if (last_state != NEARBOUND) {
      printf("Switched to NB!\n");
    }
    last_state = NEARBOUND;
  }
  else if (s == FARBOUND) {
    if (last_state != FARBOUND) {
      printf("Switched to FB!\n");
    }
    last_state = FARBOUND;
  }
  else if (s == UNBOUND) {
    if (last_state != UNBOUND) {
      data_struct->dwell_time = iteration*dt;
      printf("Switched to UB!\n");
    }
    last_state = UNBOUND;
  }

  if (iteration % (data_generation_skip_iterations*1) == 0) {
    printf("Stepping data progress (%s): %lld / %lld, %g%%                \r", run_msg,
	   iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }

  if (iteration == (max_iteration-1) - ((max_iteration-1) % data_generation_skip_iterations)) {
    printf("Finished generating stepping data (%s), process took %g seconds               \n", run_msg,
	   ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}

void make_stepping_data_file(stepping_data_struct* data, char* fname_base) {
  printf("num_steps: %d\n", data->num_steps);
  printf("dwell_time: %g\n", data->dwell_time);
  for (int i=0; i < data->step_times->get_length(); i++) {
    printf("step_time: %g\n", data->step_times->get_data()[i]);
  }
  char data_fname[200];
  sprintf(data_fname, "data/stepping_data_%s.txt", fname_base);

  FILE* data_file = fopen(data_fname, "w");

  fprintf(data_file, "#num_steps:\n");
  fprintf(data_file, "%d\n", data->num_steps);

  fprintf(data_file, "#step_lengths:\n");
  for (int i=0; i<data->step_lengths->get_length(); i++) {
    fprintf(data_file, "%g, ", data->step_lengths->get_data()[i]);
  }
  fprintf(data_file, "\n");

  fprintf(data_file, "#step_times:\n");
  for (int i=0; i<data->step_times->get_length(); i++) {
    fprintf(data_file, "%g, ", data->step_times->get_data()[i]);
  }
  fprintf(data_file, "\n");

  fclose(data_file);
}

int main(int argc, char** argv) {
  //T = 50;

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char *config_fname = new char[200];
  char *movie_config_fname = new char[200];

  sprintf(config_fname, "data/stepping_config_%s.txt", f_appended_name);
  sprintf(movie_config_fname, "data/stepping_movie_config_%s.txt", f_appended_name);

  //write_movie_config(movie_config_fname, iterations*dt);
  write_config_file(config_fname, CONFIG_INCLUDE_SKIPINFO,
		    "Initial state: onebound\nInitial conformation: equilibrium\n");

  void* job_msg[4];
  job_msg[0] = (double*) &iterations;

  double current_time = clock();
  job_msg[1] = &current_time;

  char run_msg[512];
  sprintf(run_msg, "seed = %d", (int) RAND_INIT_SEED);
  job_msg[2] = run_msg;

  stepping_data_struct data;
  data.num_steps = 0;
  data.dwell_time = 0.0;
  data.step_times = new DynArr(INIT_DYNARR_LEN);
  data.step_lengths = new DynArr(INIT_DYNARR_LEN);

  job_msg[3] = &data;

  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {eq.nma,
			    eq.fma,
			    0, 0, Ls};

  simulate(iterations*dt, RAND_INIT_SEED, BOTHBOUND, init_position,
	   stepping_data_callback, job_msg, NULL);

  make_stepping_data_file(&data, f_appended_name);
  return EXIT_SUCCESS;
}
