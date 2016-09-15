#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulation_defaults.h"

const int INIT_DYNARR_LEN = 10;
const bool logging_movie = true;
static const bool am_debugging_steps = false;

typedef struct {
  int num_steps;
  DynArr* step_times;
  DynArr* step_lengths;
  double dwell_time;
  int seed;
} stepping_data_struct;

typedef struct {
  int data;
  pthread_mutex_t* lock;
  pthread_cond_t* condvar;
} semaphore;

typedef struct {
  semaphore* sem;
  stepping_data_struct* run_data;
} stepping_thread_input_struct;

void log_stepping_data(stepping_data_struct* data_struct, void* dyn, long long iteration, long long max_iteration, State s) {
  static State last_state = BOTHBOUND;
  static double last_nbx = ((Dynein_bothbound*) dyn)->get_nbx();
  static double last_fbx = ((Dynein_bothbound*) dyn)->get_fbx();
  static double last_bothbound_iteration = 0;

  if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    if (last_state == NEARBOUND) {
      //double step_len = dyn_bb->get_fbx() - last_fbx;
      //data_struct->step_lengths->append(step_len);
      last_fbx = dyn_bb->get_fbx();

      //double step_time = (iteration - last_bothbound_iteration)*dt;
      //data_struct->step_times->append(step_time);

      data_struct->num_steps++;
      if (am_debugging_steps) printf("\nSwitched from NB to BB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    else if (last_state == FARBOUND) {
      //double step_len = dyn_bb->get_nbx() - last_nbx;
      //data_struct->step_lengths->append(step_len);
      last_nbx = dyn_bb->get_nbx();

      //double step_time = (iteration - last_bothbound_iteration)*dt;
      //data_struct->step_times->append(step_time);

      data_struct->num_steps++;
      if (am_debugging_steps) printf("\nSwitched from FB to BB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_bothbound_iteration = iteration;
    last_state = BOTHBOUND;
  }
  else if (s == NEARBOUND) {
    if (last_state != NEARBOUND) {
      if (am_debugging_steps) printf("\nSwitched from BB to NB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = NEARBOUND;
  }
  else if (s == FARBOUND) {
    if (last_state != FARBOUND) {
      if (am_debugging_steps) printf("\nSwitched from BB to FB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = FARBOUND;
  }
  else if (s == UNBOUND) {
    if (last_state != UNBOUND) {
      data_struct->dwell_time = iteration*dt;
      if (am_debugging_steps) printf("\nSwitched to UB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = UNBOUND;
  }
}

void log_stepping_movie_data(FILE* data_file, void* dyn, State s, long long iteration) {
  if (s == NEARBOUND or s == FARBOUND) {
    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
    onebound_forces dyn_ob_f = dyn_ob->get_internal();
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t"
	    "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	    "\n",
	    dyn_ob->get_state(),
	    iteration*dt,
	    dyn_ob->PE_bba, dyn_ob->PE_bma, dyn_ob->PE_ta, dyn_ob->PE_uma, 0.0,
	    dyn_ob->get_bbx(), dyn_ob->get_bby(), dyn_ob->get_bmx(), dyn_ob->get_bmy(),
	    dyn_ob->get_tx(), dyn_ob->get_ty(), dyn_ob->get_umx(), dyn_ob->get_umy(),
	    dyn_ob->get_ubx(), dyn_ob->get_uby(),
	    dyn_ob_f.bbx, dyn_ob_f.bby, dyn_ob_f.bmx, dyn_ob_f.bmy, dyn_ob_f.tx,
	    dyn_ob_f.ty, dyn_ob_f.umx, dyn_ob_f.umy, dyn_ob_f.ubx, dyn_ob_f.uby);
  }
  else if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    bothbound_forces dyn_bb_f = dyn_bb->get_internal();
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
	    iteration*dt,
	    dyn_bb->PE_nba, dyn_bb->PE_nma, dyn_bb->PE_ta, dyn_bb->PE_fma, dyn_bb->PE_fba,
	    dyn_bb->get_nbx(), dyn_bb->get_nby(), dyn_bb->get_nmx(), dyn_bb->get_nmy(),
	    dyn_bb->get_tx(), dyn_bb->get_ty(), dyn_bb->get_fmx(), dyn_bb->get_fmy(),
	    dyn_bb->get_fbx(), dyn_bb->get_fby(),
	    dyn_bb_f.nbx, dyn_bb_f.nby, dyn_bb_f.nmx, dyn_bb_f.nmy, dyn_bb_f.tx,
	    dyn_bb_f.ty, dyn_bb_f.fmx, dyn_bb_f.fmy, dyn_bb_f.fbx, dyn_bb_f.fby);
  }
  else if (s == UNBOUND) {
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t"
	    "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	    "\n",
	    UNBOUND,
	    iteration*dt,
	    0.0, 0.0, 0.0, 0.0, 0.0,
	    0.0, 0.0, 0.0, 0.0,
	    0.0, 0.0, 0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0, 0.0, 0.0, 0.0,
	    0.0, 0.0, 0.0, 0.0, 0.0);
  }
  else {
    printf("Unhandled state in stepping movie data generation!\n");
    exit(1);
  }
}

void stepping_data_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  long long max_iteration = *((long long**) job_msg)[0];
  double start_time = *((double**) job_msg)[1];
  char* run_msg = ((char**) job_msg)[2];

  stepping_data_struct* data_struct = ((stepping_data_struct**) job_msg)[3];
  FILE* data_file = ((FILE**) job_msg)[4];

  log_stepping_data(data_struct, dyn, iteration, max_iteration, s);
  if (logging_movie && iteration % data_generation_skip_iterations == 0)
    log_stepping_movie_data(data_file, dyn, s, iteration);

  if (iteration % 10 == 0) {
    printf("Stepping data progress (%s): %lld / %lld, %g%%                \r", run_msg,
	   iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }

  if (iteration == max_iteration-1) {
    printf("Finished generating stepping data (%s), process took %g seconds               \n", run_msg,
	   ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}

void stepping_data_callback_threaded(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  long long max_iteration = *((long long**) job_msg)[0];
  //double start_time = *((double**) job_msg)[1];
  //char* run_msg = ((char**) job_msg)[2];

  stepping_data_struct* data_struct = ((stepping_data_struct**) job_msg)[3];

  //printf("testing DynArrs in threads\n");
  for (int i=0; i<10; i++) {
    //data_struct->step_times->append(1.2);
  }
  // printf("data_struct: %p, step_times: %p, len: %d\n", data_struct,
  // 	 data_struct->step_times, data_struct->step_times->get_length());

  //printf("step time 7: %g, test passes\n", data_struct->step_times->get_data()[7]);

  log_stepping_data(data_struct, dyn, iteration, max_iteration, s);
}

void make_stepping_data_file(stepping_data_struct* data, char* fname_base) {
  printf("num_steps: %d\n", data->num_steps);
  printf("dwell_time: %g\n", data->dwell_time);
  for (int i=0; i < data->step_times->get_length(); i++) {
    printf("step_length: %g\n", data->step_lengths->get_data()[i]);
    printf("step_time: %g\n", data->step_times->get_data()[i]);
  }
  char data_fname[200];
  sprintf(data_fname, "data/stepping_data_%s.txt", fname_base);

  FILE* data_file = fopen(data_fname, "w");

  fprintf(data_file, "#num_steps:\n");
  fprintf(data_file, "%d\n", data->num_steps);

  fprintf(data_file, "#step_lengths:\n");
  for (int i=0; i<data->step_lengths->get_length(); i++) {
    fprintf(data_file, "%g\n", data->step_lengths->get_data()[i]);
  }
  fprintf(data_file, "\n");

  fprintf(data_file, "#step_times:\n");
  for (int i=0; i<data->step_times->get_length(); i++) {
    fprintf(data_file, "%g\n", data->step_times->get_data()[i]);
  }
  fprintf(data_file, "\n");

  fclose(data_file);
}

int old_main(int argc, char** argv) {
  T = 100;

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char *config_fname = new char[200];

  char *movie_data_fname = new char[200];
  char *movie_config_fname = new char[200];

  sprintf(config_fname, "data/stepping_config_%s.txt", f_appended_name);
  sprintf(movie_data_fname, "data/movie_%s.txt", f_appended_name);
  sprintf(movie_config_fname, "data/movie_config_%s.txt", f_appended_name);

  write_movie_config(movie_config_fname, iterations*dt);
  write_config_file(config_fname, 0, "");

  void* job_msg[5];
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
  job_msg[4] = fopen(movie_data_fname, "w");

  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {eq.nma,
			    eq.fma,
			    0, 0, Ls};

  simulate(iterations*dt, RAND_INIT_SEED, BOTHBOUND, init_position,
	   stepping_data_callback, job_msg, NULL);

  make_stepping_data_file(&data, f_appended_name);

  fclose((FILE*) job_msg[4]);
  delete data.step_times;
  delete data.step_lengths;

  return EXIT_SUCCESS;
}

void prepare_data_files(int num_args, char* f_appended_name) {
  if (num_args != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char *config_fname = new char[200];

  char *movie_data_fname = new char[200];
  char *movie_config_fname = new char[200];

  sprintf(config_fname, "data/stepping_config_%s.txt", f_appended_name);
  sprintf(movie_data_fname, "data/movie_%s.txt", f_appended_name);
  sprintf(movie_config_fname, "data/movie_config_%s.txt", f_appended_name);

  write_movie_config(movie_config_fname, iterations*dt);
  write_config_file(config_fname, 0, "");
}

void run_simulation(stepping_data_struct* run_data, int simid) {
  void* job_msg[4];
  double current_time = clock();
  char run_msg[512];
  sprintf(run_msg, "seed = %d", (int) run_data->seed);
  
  job_msg[0] = (double*) &iterations;
  job_msg[1] = &current_time;
  job_msg[2] = run_msg;
  job_msg[3] = run_data;

  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {eq.nma, eq.fma, 0, 0, Ls};

  simulate(iterations*dt, run_data->seed, BOTHBOUND, init_position,
  	   stepping_data_callback_threaded, job_msg, NULL);
}

void* stepping_worker_thread(void* arg) {
  stepping_data_struct* run_data = ((stepping_thread_input_struct*) arg)->run_data;
  semaphore* simulations_to_go = ((stepping_thread_input_struct*) arg)->sem;
  while(pthread_mutex_lock(simulations_to_go->lock) == 0) {
    if (simulations_to_go->data <= 0) {              // simulations have all run, time to exit
      pthread_mutex_unlock(simulations_to_go->lock);
      pthread_exit(NULL);
    }
    else {                                           // decrement semaphore and run a simulation
      int simulation_num = simulations_to_go->data;
      simulations_to_go->data--;
      pthread_mutex_unlock(simulations_to_go->lock);
      printf("--Running simulation %d\n", simulation_num);
      run_simulation(&run_data[simulation_num], simulation_num);
      printf("--Completed simulation %d\n", simulation_num);
      pthread_cond_signal(simulations_to_go->condvar);
    }
  }
  return NULL;
}

int main(int argc, char** argv) {
  T = 100;

  int num_threads = 4;
  semaphore* simulations_to_go = new semaphore;
  simulations_to_go->data = num_dynein_runs;
  simulations_to_go->lock = new pthread_mutex_t;
  simulations_to_go->condvar = new pthread_cond_t;

  pthread_mutex_init(simulations_to_go->lock, NULL);
  pthread_cond_init(simulations_to_go->condvar, NULL);

  stepping_thread_input_struct* stepping_data = new stepping_thread_input_struct;
  stepping_data->run_data = new stepping_data_struct[num_dynein_runs];
  stepping_data->sem = simulations_to_go;

  for (int i=0; i<num_dynein_runs; i++) {
    stepping_data->run_data[i].num_steps = 0;
    stepping_data->run_data[i].dwell_time = 0.0;
    stepping_data->run_data[i].step_times = new DynArr(INIT_DYNARR_LEN);
    stepping_data->run_data[i].step_lengths = new DynArr(INIT_DYNARR_LEN);
    stepping_data->run_data[i].seed = 0;
  }

  pthread_t* threads = new pthread_t[num_threads];

  pthread_mutex_lock(simulations_to_go->lock);

  for (int n=0; n<num_threads; n++) {
    pthread_create(&threads[n], NULL, stepping_worker_thread, (void*) stepping_data);
  }

  while (pthread_cond_wait(simulations_to_go->condvar, simulations_to_go->lock) == 0) {
    if (simulations_to_go->data <= 0) break;
  }
  pthread_mutex_unlock(simulations_to_go->lock);

  for (int n=0; n<num_threads; n++) {
    pthread_join(threads[n], NULL);
  }

  for (int m=0; m<num_dynein_runs; m++) {
    printf("num_steps for %d: %d\n", m, stepping_data->run_data[m].num_steps);
  }

  //make_stepping_data_file(&data, f_appended_name);

  // fclose((FILE*) job_msg[4]);
  // delete data.step_times;
  // delete data.step_lengths;

  return EXIT_SUCCESS;
}


// void* job_msg[5];
//   job_msg[0] = (double*) &iterations;

//   double current_time = clock();
//   job_msg[1] = &current_time;

//   char run_msg[512];
//   sprintf(run_msg, "seed = %d", (int) RAND_INIT_SEED);
//   job_msg[2] = run_msg;

//   stepping_data_struct data;
//   data.num_steps = 0;
//   data.dwell_time = 0.0;
//   data.step_times = new DynArr(INIT_DYNARR_LEN);
//   data.step_lengths = new DynArr(INIT_DYNARR_LEN);

//   job_msg[3] = &data;
//   job_msg[4] = fopen(movie_data_fname, "w");

//   bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
//   double init_position[] = {eq.nma,
// 			    eq.fma,
// 			    0, 0, Ls};

//   simulate(iterations*dt, RAND_INIT_SEED, BOTHBOUND, init_position,
// 	   stepping_data_callback, job_msg, NULL);

