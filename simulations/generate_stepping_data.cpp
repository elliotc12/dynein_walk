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

const bool display_step_info = false;
const bool display_progress = false;

void log_stepping_data(FILE* data_file, void* dyn, long long iteration, long long max_iteration, State s) {
  static State  last_state = BOTHBOUND;
  static double last_bothbound_iteration = 0;

  if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    if (last_state == NEARBOUND or last_state == FARBOUND) {
      fprintf(data_file, "%.4e %.4e %.4e %.4e\n", last_bothbound_iteration*dt, iteration*dt, dyn_bb->get_nbx(), dyn_bb->get_fbx());
      if (display_step_info) printf("\nSwitched to BB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }

    last_bothbound_iteration = iteration;
    last_state = BOTHBOUND;
  }
  else if (s == NEARBOUND) {
    if (last_state != NEARBOUND) {
      if (display_step_info) printf("\nSwitched from BB to NB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = NEARBOUND;
  }
  else if (s == FARBOUND) {
    if (last_state != FARBOUND) {
      if (display_step_info) printf("\nSwitched from BB to FB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = FARBOUND;
  }
  else if (s == UNBOUND) {
    if (last_state != UNBOUND) {
      fprintf(data_file, "%.4e %.4e %.4e %.4e\n", last_bothbound_iteration*dt, iteration*dt, 0.0, 0.0);
      if (display_step_info) printf("\nSwitched to UB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
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

  FILE* stepping_data_file = ((FILE**) job_msg)[3];
  FILE* movie_data_file = ((FILE**) job_msg)[4];

  bool am_making_movie = *((bool**) job_msg)[5];

  log_stepping_data(stepping_data_file, dyn, iteration, max_iteration, s);

  if (am_making_movie && iteration % stepping_movie_framerate == 0)
    log_stepping_movie_data(movie_data_file, dyn, s, iteration);

  if (max_iteration > 0) {
    if (iteration % (int)5e5 == 0 && iteration != max_iteration && display_progress) {
      printf("Stepping data progress (%s): %lld / %lld, %g%%\n", run_msg,
  	     iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    }

    if (iteration == max_iteration and display_progress) {
      printf("Finished generating stepping data (%s), process took %g seconds\n", run_msg,
  	     ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  else if (max_iteration == 0 and iteration % (int)5e5 == 0) {
    //printf("Stepping data progress (%s): %.2g seconds\n", run_msg, iteration*dt);
    fflush(stepping_data_file);
  }
}

void set_input_variables(int argc, char** argv, char* run_name, bool* am_making_movie) {
  char c;
  *run_name = 0;

  static struct option long_options[] =
    {
      {"Ls",       required_argument,    0, 'a'},
      {"Lt",       required_argument,    0, 'b'},
      {"cb",       required_argument,    0, 'c'},
      {"cm",       required_argument,    0, 'd'},
      {"ct",       required_argument,    0, 'e'},
      {"T",        required_argument,    0, 'f'},
      {"name",     required_argument,    0, 'g'},
      {"seed",     required_argument,    0, 'h'},
      {"k_b",      required_argument,    0, 'i'},
      {"k_ub",     required_argument,    0, 'j'},
      {"k_ub_ob",  required_argument,    0, 'k'},
      {"movie",  no_argument, (int*) am_making_movie, 1},
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
    case 'i':
      low_affinity_binding_rate = strtod(optarg, NULL);
      break;
    case 'j':
      high_affinity_binding_rate = strtod(optarg, NULL);
      break;
    case 'k':
      low_affinity_unbinding_rate = strtod(optarg, NULL);
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
    printf("Improper usage, all options need an option name like --Ls or --T!\n");
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char** argv) {
  char* run_name = new char[100];
  bool am_making_movie = 0;

  set_input_variables(argc, argv, run_name, &am_making_movie);

  if (*run_name == 0) {
    sprintf(run_name, "ls-%.3e,lt-%.3e,k_b-%.3e,cb-%.3e,cm-%.3e,ct-%.3e,T-%.3e", Ls, Lt, low_affinity_binding_rate, cb, cm, ct, T);
  }

  char *stepping_data_fname = new char[200];
  char *stepping_config_fname = new char[200];

  char *movie_data_fname = new char[200];
  char *movie_config_fname = new char[200];

  sprintf(stepping_data_fname, "data/stepping_data_%s.txt", run_name);
  sprintf(stepping_config_fname, "data/stepping_config_%s.txt", run_name);
  sprintf(movie_data_fname, "data/movie_data_%s.txt", run_name);
  sprintf(movie_config_fname, "data/movie_config_%s.txt", run_name);

  write_config_file(stepping_config_fname, 0, "");

  write_movie_config(movie_config_fname, iterations*dt);

  double current_time = clock();
  int indefinite_run = 0;

  void* job_msg[6];
  job_msg[0] = &indefinite_run;
  job_msg[1] = &current_time;
  job_msg[2] = run_name;
  job_msg[3] = fopen(stepping_data_fname, "w");
  job_msg[4] = fopen(movie_data_fname, "w");
  job_msg[5] = &am_making_movie;

  printf("fname: %s\n", stepping_data_fname);
  fprintf((FILE*) job_msg[3], "#time_unbind, time_bind, nbx, fbx\n");

  if (errno) {
    perror("Error opening stepping data or movie file.\n");
    exit(errno);
  }

  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {eq.nma,
			    eq.fma,
			    0, 0, Ls};

  simulate(0, RAND_INIT_SEED, BOTHBOUND, init_position, stepping_data_callback, job_msg, NULL);

  fclose((FILE*) job_msg[3]);
  fclose((FILE*) job_msg[4]);

  return EXIT_SUCCESS;
}
