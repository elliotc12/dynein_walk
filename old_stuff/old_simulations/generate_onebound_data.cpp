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

void write_onebound_data_callback_no_avging(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % data_generation_skip_iterations != 0) return;
  else {
    assert(s == NEARBOUND);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];
    int iter = iteration / data_generation_skip_iterations;

    onebound_data_generate_struct* data_mem = ((onebound_data_generate_struct**) job_msg)[3];

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

    double et = 0.5*kb*T;
    onebound_data_generate_struct data;
    data.time = iteration*dt;
    data.bba_PE = dyn_ob->PE_bba / et;
    data.bma_PE = dyn_ob->PE_bma / et;
    data.ta_PE  = dyn_ob->PE_ta  / et;
    data.uma_PE = dyn_ob->PE_uma / et;
    data.bba = dyn_ob->get_bba();
    data.bma = dyn_ob->get_bma() + M_PI - dyn_ob->get_bba();
    data.ta  = dyn_ob->get_uma() - dyn_ob->get_bma();
    data.uma = dyn_ob->get_uma() + M_PI - dyn_ob->get_uba();

    onebound_forces f = dyn_ob->get_internal();
    data.f_bbx = f.bbx;   data.f_bby = f.bby;
    data.f_bmx = f.bmx;   data.f_bmy = f.bmy;
    data.f_tx  = f.tx;    data.f_ty  = f.ty;
    data.f_umx = f.umx;   data.f_umy = f.umy;
    data.f_ubx = f.ubx;   data.f_uby = f.uby;

    data.bbx = dyn_ob->get_bbx();   data.bby = dyn_ob->get_bby();
    data.bmx = dyn_ob->get_bmx();   data.bmy = dyn_ob->get_bmy();
    data.tx  = dyn_ob->get_tx();    data.ty  = dyn_ob->get_ty();
    data.umx = dyn_ob->get_umx();   data.umy = dyn_ob->get_umy();
    data.ubx = dyn_ob->get_ubx();   data.uby = dyn_ob->get_uby();

    memcpy(&data_mem[iter], &data, sizeof(onebound_data_generate_struct));

    if (iteration % (data_generation_skip_iterations*1) == 0) {
      printf("PE calculation progress (%s): %lld / %lld, %g%%                \r", run_msg,
      	     iteration, max_iteration, ((double) iteration) / max_iteration * 100);
      fflush(NULL);
    }

    if (iteration == 0) {
      msync(data_mem, sizeof(onebound_data_generate_struct), MS_ASYNC); // asynchronously write mmap in memory to file
    }
    else if (iteration % msync_after_num_writes == 0) {
      msync(&data_mem[iter-msync_after_num_writes+1],
	    msync_after_num_writes*sizeof(onebound_data_generate_struct), MS_ASYNC);
    }

    if (iteration == (max_iteration-1) - ((max_iteration-1) % data_generation_skip_iterations)) {
      printf("Finished generating ob PE data (%s), process took %g seconds               \n", run_msg,
	     ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
}

void write_onebound_data_callback_with_avging(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  static onebound_data_generate_struct data_sum;
  assert(s == NEARBOUND or s == FARBOUND);
  long long max_iteration = *((long long**) job_msg)[0];
  double start_time = *((double**) job_msg)[1];
  char* run_msg = ((char**) job_msg)[2];
  int iter = iteration / data_generation_skip_iterations;

  onebound_data_generate_struct* data_mem = ((onebound_data_generate_struct**) job_msg)[3];

  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

  double et = 0.5*kb*T;
  data_sum.time += iteration*dt;
  data_sum.bba_PE += dyn_ob->PE_bba / et;
  data_sum.bma_PE += dyn_ob->PE_bma / et;
  data_sum.ta_PE  += dyn_ob->PE_ta  / et;
  data_sum.uma_PE += dyn_ob->PE_uma / et;
  data_sum.bba += dyn_ob->get_bba();
  data_sum.bma += dyn_ob->get_bma() + M_PI - dyn_ob->get_bba();
  data_sum.ta  += dyn_ob->get_uma() - dyn_ob->get_bma();
  data_sum.uma += dyn_ob->get_uma() + M_PI - dyn_ob->get_uba();

  onebound_forces f = dyn_ob->get_internal();
  data_sum.f_bbx += f.bbx;   data_sum.f_bby += f.bby;
  data_sum.f_bmx += f.bmx;   data_sum.f_bmy += f.bmy;
  data_sum.f_tx  += f.tx;    data_sum.f_ty  += f.ty;
  data_sum.f_umx += f.umx;   data_sum.f_umy += f.umy;
  data_sum.f_ubx += f.ubx;   data_sum.f_uby += f.uby;

  data_sum.bbx += dyn_ob->get_bbx();   data_sum.bby += dyn_ob->get_bby();
  data_sum.bmx += dyn_ob->get_bmx();   data_sum.bmy += dyn_ob->get_bmy();
  data_sum.tx  += dyn_ob->get_tx();    data_sum.ty  += dyn_ob->get_ty();
  data_sum.umx += dyn_ob->get_umx();   data_sum.umy += dyn_ob->get_umy();
  data_sum.ubx += dyn_ob->get_ubx();   data_sum.uby += dyn_ob->get_uby();

  if (iteration % data_generation_skip_iterations == 0) {
    data_sum.time /= data_generation_skip_iterations;
    data_sum.bba_PE /= data_generation_skip_iterations;
    data_sum.bma_PE /= data_generation_skip_iterations;
    data_sum.ta_PE  /= data_generation_skip_iterations;
    data_sum.uma_PE /= data_generation_skip_iterations;
    data_sum.bba  /= data_generation_skip_iterations;
    data_sum.bma  /= data_generation_skip_iterations;
    data_sum.ta   /= data_generation_skip_iterations;
    data_sum.uma  /= data_generation_skip_iterations;

    data_sum.f_bbx  /= data_generation_skip_iterations;   data_sum.f_bby  /= data_generation_skip_iterations;
    data_sum.f_bmx  /= data_generation_skip_iterations;   data_sum.f_bmy  /= data_generation_skip_iterations;
    data_sum.f_tx   /= data_generation_skip_iterations;   data_sum.f_ty   /= data_generation_skip_iterations;
    data_sum.f_umx  /= data_generation_skip_iterations;   data_sum.f_umy  /= data_generation_skip_iterations;
    data_sum.f_ubx  /= data_generation_skip_iterations;   data_sum.f_uby  /= data_generation_skip_iterations;

    data_sum.bbx  /= data_generation_skip_iterations;   data_sum.bby  /= data_generation_skip_iterations;
    data_sum.bmx  /= data_generation_skip_iterations;   data_sum.bmy  /= data_generation_skip_iterations;
    data_sum.tx   /= data_generation_skip_iterations;   data_sum.ty   /= data_generation_skip_iterations;
    data_sum.umx  /= data_generation_skip_iterations;   data_sum.umy  /= data_generation_skip_iterations;
    data_sum.ubx  /= data_generation_skip_iterations;   data_sum.uby  /= data_generation_skip_iterations;

    memcpy(&data_mem[iter], &data_sum, sizeof(onebound_data_generate_struct));

    data_sum.time = 0;
    data_sum.bba_PE = 0;
    data_sum.bma_PE = 0;
    data_sum.ta_PE  = 0;
    data_sum.uma_PE = 0;
    data_sum.bba = 0;
    data_sum.bma = 0;
    data_sum.ta  = 0;
    data_sum.uma = 0;

    data_sum.f_bbx = 0;   data_sum.f_bby = 0;
    data_sum.f_bmx = 0;   data_sum.f_bmy = 0;
    data_sum.f_tx  = 0;   data_sum.f_ty  = 0;
    data_sum.f_umx = 0;   data_sum.f_umy = 0;
    data_sum.f_ubx = 0;   data_sum.f_uby = 0;

    data_sum.bbx = 0;   data_sum.bby = 0;
    data_sum.bmx = 0;   data_sum.bmy = 0;
    data_sum.tx  = 0;   data_sum.ty  = 0;
    data_sum.umx = 0;   data_sum.umy = 0;
    data_sum.ubx = 0;   data_sum.uby = 0;
  }

  if (iteration % (data_generation_skip_iterations*1) == 0) {
    printf("PE calculation progress (%s): %lld / %lld, %g%%                \r", run_msg,
	   iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }

  if (iteration == 0) {
    msync(data_mem, sizeof(onebound_data_generate_struct), MS_ASYNC); // asynchronously write mmap in memory to file
  }
  else if (iteration % msync_after_num_writes == 0) {
    msync(&data_mem[iter-msync_after_num_writes+1],
	  msync_after_num_writes*sizeof(onebound_data_generate_struct), MS_ASYNC);
  }

  if (iteration == (max_iteration-1) - ((max_iteration-1) % data_generation_skip_iterations)) {
    printf("Finished generating ob PE data (%s), process took %g seconds               \n", run_msg,
	   ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}

int main(int argc, char** argv) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;
  ONEBOUND_UNBINDING_FORCE = std::numeric_limits<double>::infinity();

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

  sprintf(data_fname, "data/onebound_data_%s.bin", f_appended_name);
  sprintf(config_fname, "data/ob_config_%s.txt", f_appended_name);
  sprintf(movie_config_fname, "data/movie_config_%s.txt", f_appended_name);

  write_movie_config(movie_config_fname, iterations*dt);
  write_config_file(config_fname, CONFIG_INCLUDE_SKIPINFO,
		    "Initial state: onebound\nInitial conformation: equilibrium\n");

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

  ftruncate(data_fd, iters*sizeof(onebound_data_generate_struct));
  if (errno) {
    perror("Error ftruncating data file");
    exit(errno);
  }

  void* data_mem = mmap(NULL, iters*sizeof(onebound_data_generate_struct), PROT_WRITE, MAP_SHARED, data_fd, 0);
  if (data_mem == MAP_FAILED) {
    perror("Error using mmap: ");
    exit(EXIT_FAILURE);
  }

  job_msg[3] = data_mem;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba,
			    eq.bma + eq.bba - M_PI,
			    eq.ta + eq.bma + eq.bba - M_PI,
			    eq.ta + eq.bma + eq.bba - eq.uma,
			    0, 0};

  simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	   write_onebound_data_callback_with_avging, job_msg, NULL);

  munmap(data_mem, iters*sizeof(onebound_data_generate_struct));
  close(data_fd);
  return EXIT_SUCCESS;
}
