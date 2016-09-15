#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulation_defaults.h"

const int INIT_DYNARR_LEN = 10;
const bool logging_movie = true;
const bool am_debugging_steps = false;
const bool display_progress = false;
const bool display_step_output = false;

typedef struct {
  int num_steps;
  DynArr* step_times;
  DynArr* step_lengths;
  double dwell_time;
} stepping_data_struct;

void log_stepping_data(stepping_data_struct* data_struct, void* dyn, long long iteration, long long max_iteration, State s) {
  static State last_state = BOTHBOUND;
  static double last_nbx = ((Dynein_bothbound*) dyn)->get_nbx();
  static double last_fbx = ((Dynein_bothbound*) dyn)->get_fbx();
  static double last_bothbound_iteration = 0;

  if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    if (last_state == NEARBOUND) {
      double step_len = dyn_bb->get_fbx() - last_fbx;
      data_struct->step_lengths->append(step_len);
      last_fbx = dyn_bb->get_fbx();

      double step_time = (iteration - last_bothbound_iteration)*dt;
      data_struct->step_times->append(step_time);

      data_struct->num_steps++;
      if (am_debugging_steps) printf("\nSwitched from NB to BB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    else if (last_state == FARBOUND) {
      double step_len = dyn_bb->get_nbx() - last_nbx;
      data_struct->step_lengths->append(step_len);
      last_nbx = dyn_bb->get_nbx();

      double step_time = (iteration - last_bothbound_iteration)*dt;
      data_struct->step_times->append(step_time);

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

  if (iteration % 10 == 0 && iteration != max_iteration && display_progress) {
    printf("Stepping data progress (%s): %lld / %lld, %g%%                \r", run_msg,
	   iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }

  if (iteration == max_iteration) {
    printf("Finished generating stepping data (%s), process took %g seconds               \n", run_msg,
	   ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}

void make_stepping_data_file(stepping_data_struct* data, char* fname_base) {
  if (display_step_output) printf("num_steps: %d\n", data->num_steps);
  if (display_step_output) printf("dwell_time: %g\n", data->dwell_time);
  for (int i=0; i < data->step_times->get_length(); i++) {
    if (display_step_output) printf("step_length: %g\n", data->step_lengths->get_data()[i]);
    if (display_step_output) printf("step_time: %g\n", data->step_times->get_data()[i]);
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

void set_input_variables(int argc, char** argv, char* run_name) {
  char c;
  *run_name = 0;

  static struct option long_options[] =
    {
      {"Ls",     required_argument,    0, 'a'},
      {"Lt",     required_argument,    0, 'b'},
      {"cb",     required_argument,    0, 'c'},
      {"cm",     required_argument,    0, 'd'},
      {"ct",     required_argument,    0, 'e'},
      {"T",      required_argument,    0, 'f'},
      {"name",   required_argument,    0, 'g'},
      {"seed",   required_argument,    0, 'h'},
      {0, 0, 0, 0}
    };

  int option_index = 0;

  while ((c = getopt_long(argc, argv, "a:b:c:d:e:f:", long_options, &option_index)) != -1) {
    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0) // option set a flag
	break;
      else {
	printf ("unknown option %s", long_options[option_index].name);
	if (optarg)
	  printf (" with arg %s", optarg);
	printf ("\n");
	break;
      }
    case 'a':
      Ls = strtod(optarg, NULL);
      break;
    case 'b':
      Lt = strtod(optarg, NULL);
      break;
    case 'c':
      cb = strtod(optarg, NULL);
      break;
    case 'd':
      cm = strtod(optarg, NULL);
      break;
    case 'e':
      ct = strtod(optarg, NULL);
      break;
    case 'f':
      T = strtod(optarg, NULL);
      break;
    case 'g':
      strcpy(run_name, optarg);
      break;
    case 'h':
      RAND_INIT_SEED = atoi(optarg);
      break;
    case '?':
      printf("Some other unknown getopt error.\n");
      exit(EXIT_FAILURE);
    default:
      printf("Default case in getopt: uh-oh!\n");
      exit(EXIT_FAILURE);
    }
  }

  if (optind != argc) {
    printf("Improper usage, all options need an option name like -ls or -T!\n");
    exit(EXIT_FAILURE);
  }
  if (run_name[0] == 0) {
    printf("name must be specified!\n");
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char** argv) {

  char* f_appended_name = new char[100];
  set_input_variables(argc, argv, f_appended_name);

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
