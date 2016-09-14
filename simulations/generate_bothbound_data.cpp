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

void write_bothbound_data_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % data_generation_skip_iterations != 0) return;
  else {
    assert(s == BOTHBOUND);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];
    int iter = iteration / data_generation_skip_iterations;

    bothbound_data_generate_struct* data_mem = ((bothbound_data_generate_struct**) job_msg)[3];

    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;

    double et = 0.5*kb*T;
    bothbound_data_generate_struct data;
    data.time = iteration*dt;
    data.nba_PE = dyn_bb->PE_nba / et;
    data.nma_PE = dyn_bb->PE_nma / et;
    data.ta_PE  = dyn_bb->PE_ta  / et;
    data.fma_PE = dyn_bb->PE_fma / et;
    data.fba_PE = dyn_bb->PE_fba / et;

    data.nba = dyn_bb->get_nba();
    data.nma = dyn_bb->get_nma();
    data.ta  = dyn_bb->get_ta();
    data.fma = dyn_bb->get_fma();
    data.fba = dyn_bb->get_fba();

    bothbound_forces f = dyn_bb->get_internal();
    data.f_nbx = f.nbx;   data.f_nby = f.nby;
    data.f_nmx = f.nmx;   data.f_nmy = f.nmy;
    data.f_tx  = f.tx;    data.f_ty  = f.ty; 
    data.f_fmx = f.fmx;   data.f_fmy = f.fmy;
    data.f_fbx = f.fbx;   data.f_fby = f.fby;

    data.nbx = dyn_bb->get_nbx();   data.nby = dyn_bb->get_nby();
    data.nmx = dyn_bb->get_nmx();   data.nmy = dyn_bb->get_nmy();
    data.tx  = dyn_bb->get_tx();    data.ty  = dyn_bb->get_ty();
    data.fmx = dyn_bb->get_fmx();   data.fmy = dyn_bb->get_fmy();
    data.fbx = dyn_bb->get_fbx();   data.fby = dyn_bb->get_fby();

    memcpy(&data_mem[iter], &data, sizeof(bothbound_data_generate_struct));

    if (iteration % (data_generation_skip_iterations*1) == 0) {
      printf("PE calculation progress (%s): %lld / %lld, %g%%                \r", run_msg,
      	     iteration, max_iteration, ((double) iteration) / max_iteration * 100);
      fflush(NULL);
    }

    if (iteration == 0) {
      msync(data_mem, sizeof(bothbound_data_generate_struct), MS_ASYNC); // asynchronously write mmap in memory to file
    }
    else if (iteration % msync_after_num_writes == 0) {
      msync(&data_mem[iter-msync_after_num_writes+1],
      msync_after_num_writes*sizeof(bothbound_data_generate_struct), MS_ASYNC);
    }

    if (iteration == (max_iteration-1) - ((max_iteration-1) % data_generation_skip_iterations)) {
      printf("Finished generating bb PE data (%s), process took %g seconds               \n", run_msg,
	     ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
}

int main(int argc, char** argv) {
  BOTHBOUND_UNBINDING_FORCE = std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  T = 100;

  int iters = iterations / data_generation_skip_iterations;

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char *config_fname = new char[200];
  char *movie_config_fname = new char[200];
  char *data_fname = new char[200];

  sprintf(data_fname, "data/bothbound_data_%s.bin", f_appended_name);
  sprintf(config_fname, "data/bb_config_%s.txt", f_appended_name);
  sprintf(movie_config_fname, "data/movie_config_%s.txt", f_appended_name);

  write_movie_config(movie_config_fname, iterations*dt);
  write_config_file(config_fname, CONFIG_INCLUDE_SKIPINFO,
		    "Initial state: bothbound\nInitial conformation: equilibrium\n");

  void* job_msg[4];
  job_msg[0] = (double*) &iterations;

  double current_time = clock();
  job_msg[1] = &current_time;

  char run_msg[512];
  sprintf(run_msg, "seed = %d", (int) RAND_INIT_SEED);
  job_msg[2] = run_msg;

  int data_fd = open(data_fname, O_RDWR | O_CREAT | O_TRUNC, S_IRWXU | S_IRGRP | S_IROTH);
  if (errno) {
    perror("Error creating data file");
    exit(errno);
  }

  ftruncate(data_fd, iters*sizeof(bothbound_data_generate_struct));
  if (errno) {
    perror("Error ftruncating data file");
    exit(errno);
  }

  void* data_mem = mmap(NULL, iters*sizeof(bothbound_data_generate_struct), PROT_WRITE, MAP_SHARED, data_fd, 0);
  if (data_mem == MAP_FAILED) {
    perror("Error using mmap: ");
    exit(EXIT_FAILURE);
  }

  job_msg[3] = data_mem;

  //bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {0.6*M_PI,
			    0.6*M_PI,
			    0, 0, Ls};

  simulate(iterations*dt, RAND_INIT_SEED, BOTHBOUND, init_position,
	   write_bothbound_data_callback, job_msg, NULL);

  munmap(data_mem, iters*sizeof(bothbound_data_generate_struct));
  close(data_fd);
  return EXIT_SUCCESS;
}
