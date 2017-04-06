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

extern movie_data_struct* on_crash_old_movie_data_global_ptr;
extern movie_data_struct* on_crash_new_movie_data_global_ptr;
extern char* crash_movie_file_name_global;

bool am_making_movie = true;
bool am_debugging_onebound = false;

int num_movie_writes = 1e4;
// bytes per movie write: 213, 200mb bytes max movie size

static const long MAX_FILESIZE_PERMITTED = 1<<30;
static int NUM_STEPS = 0;

void on_crash_write_movie_buffer();

struct job_msg_t {
  long long max_iteration;
  double start_time;
  char *run_msg;
  FILE *stepping_data_file;
  FILE *movie_data_file;
};

void check_for_quitting_conditions(double time_run) {
  // P(stepping k times in t seconds) = ((t/0.06)^k*e^(t/0.06)) / k!
  if (time_run > 0.3 and NUM_STEPS == 0 and am_exiting_on_improbable_stepping) {
    printf("Zero steps in 0.3 seconds. There's a 0.67%% chance of real dynein doing this; exiting early.\n");
    exit(0);
  }
  if (time_run > 0.3 and NUM_STEPS >= 10 and am_exiting_on_improbable_stepping) {
    printf("Over 10 steps in 0.3 seconds. There's less than a 1.4%% chance of real dynein doing this; exiting early.\n");
    exit(0);
  }
  if (time_run > 0.001 and NUM_STEPS > 5 and am_exiting_on_improbable_stepping) {
    printf("Over 5 steps in 0.001 seconds. There's less than a 1e-11 chance of real dynein doing this; exiting early.\n");
    exit(0);
  }
  else if (time_run > 0.7 and am_exiting_on_improbable_stepping) {
    // printf("There's a 97%% chance of real dynein having stepped 5 times in 0.7 seconds, exiting successfully.\n");
    printf("Exiting normally after 0.7 seconds.\n");
    exit(0);
  }
}

void zero_movie_struct(movie_data_struct* data) {
  for (int i = 0; i < MOVIE_BUFFER_SIZE; i++) {
    data[i].time = -1.0;
  }
}

void log_stepping_data(FILE* data_file, void* dyn, long long iteration, long long max_iteration, State s) {
  static State  last_state = BOTHBOUND;
  static double last_bothbound_iteration = 0;

  if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    if (last_state == NEARBOUND or last_state == FARBOUND) {
      if (using_variable_timestep) {
	char temp_buffer[300];
	int num_chars_added = sprintf(temp_buffer, "%.15e %.15e %.15e %.15e\n", last_bothbound_iteration*dt, iteration*dt, dyn_bb->get_nbx(), dyn_bb->get_fbx());
	strcpy(&variable_ts_stepping_print_buffer[variable_ts_stepping_print_buffer_index], temp_buffer);
	variable_ts_stepping_print_buffer_index += num_chars_added;
      }
      else {
	fprintf(data_file, "%.15e %.15e %.15e %.15e\n", last_bothbound_iteration*dt, iteration*dt, dyn_bb->get_nbx(), dyn_bb->get_fbx());
	fflush(data_file);
      }
      NUM_STEPS++;
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
      fflush(data_file);
      if (display_step_info) printf("\nSwitched to UB at %.1f%%!\n", ((double)iteration)/max_iteration*100);
    }
    last_state = UNBOUND;
  }
}

void log_stepping_movie_data(FILE* data_file, void* dyn, State s, long long iteration) {
  static int buffer_position = 0;
  if (!am_only_writing_on_crash or (am_debugging_onebound and s != BOTHBOUND)) {
    if (--num_movie_writes > 0) {
      if (num_movie_writes == 1) printf("about to exceed movie printing line #\n");
      const char *format = "%d\t"
	"%.10g\t"
	"%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	"%g\t%g\t%g\t%g\t"
	"%g\t%g\t%g\t%g\t"
	"%g\t%g\t"
	"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	"\n";
      if (s == NEARBOUND or s == FARBOUND) {
	Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
	onebound_forces dyn_ob_f = dyn_ob->get_internal();
	fprintf(data_file, format,
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
	fprintf(data_file, format,
		BOTHBOUND,
		iteration*dt,
		dyn_bb->PE_nba, dyn_bb->PE_nma, dyn_bb->PE_ta, dyn_bb->PE_fma, dyn_bb->PE_fba,
		dyn_bb->get_nbx(), dyn_bb->get_nby(), dyn_bb->get_nmx(), dyn_bb->get_nmy(),
		dyn_bb->get_tx(), dyn_bb->get_ty(), dyn_bb->get_fmx(), dyn_bb->get_fmy(),
		dyn_bb->get_fbx(), dyn_bb->get_fby(),
		dyn_bb_f.nbx, dyn_bb_f.nby, dyn_bb_f.nmx, dyn_bb_f.nmy, dyn_bb_f.tx,
		dyn_bb_f.ty, dyn_bb_f.fmx, dyn_bb_f.fmy, dyn_bb_f.fbx, dyn_bb_f.fby);
      }
      else if (s != UNBOUND){
	printf("Unhandled state in stepping movie data generation!\n");
	exit(1);
      }
    }
  }
  else { // else write to on-crash-print buffer structs
    if (buffer_position == MOVIE_BUFFER_SIZE) {
      movie_data_struct* temp_ptr;
      temp_ptr = on_crash_old_movie_data_global_ptr;
      on_crash_old_movie_data_global_ptr = on_crash_new_movie_data_global_ptr;
      on_crash_new_movie_data_global_ptr = temp_ptr;
      zero_movie_struct(on_crash_new_movie_data_global_ptr);
      buffer_position = 0;
    }

    movie_data_struct* new_movie_buffer = on_crash_new_movie_data_global_ptr;

    if (s == NEARBOUND or s == FARBOUND) {
      Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
      onebound_forces dyn_ob_f = dyn_ob->get_internal();
      new_movie_buffer[buffer_position].state = dyn_ob->get_state();
      new_movie_buffer[buffer_position].time = iteration*dt;
      new_movie_buffer[buffer_position].PE_1 = dyn_ob->PE_bba;
      new_movie_buffer[buffer_position].PE_2 = dyn_ob->PE_bma;
      new_movie_buffer[buffer_position].PE_3 = dyn_ob->PE_uma;
      new_movie_buffer[buffer_position].PE_4 = 0.0;

      new_movie_buffer[buffer_position].x_1 = dyn_ob->get_bbx();
      new_movie_buffer[buffer_position].x_2 = dyn_ob->get_bmx();
      new_movie_buffer[buffer_position].x_3 = dyn_ob->get_tx();
      new_movie_buffer[buffer_position].x_4 = dyn_ob->get_umx();
      new_movie_buffer[buffer_position].x_5 = dyn_ob->get_ubx();
      new_movie_buffer[buffer_position].y_1 = dyn_ob->get_bby();
      new_movie_buffer[buffer_position].y_2 = dyn_ob->get_bmy();
      new_movie_buffer[buffer_position].y_3 = dyn_ob->get_ty();
      new_movie_buffer[buffer_position].y_4 = dyn_ob->get_umy();
      new_movie_buffer[buffer_position].y_5 = dyn_ob->get_uby();

      new_movie_buffer[buffer_position].fx_1 = dyn_ob_f.bbx;
      new_movie_buffer[buffer_position].fx_2 = dyn_ob_f.bmx;
      new_movie_buffer[buffer_position].fx_3 = dyn_ob_f.tx;	
      new_movie_buffer[buffer_position].fx_4 = dyn_ob_f.umx;
      new_movie_buffer[buffer_position].fx_5 = dyn_ob_f.ubx;
      new_movie_buffer[buffer_position].fy_1 = dyn_ob_f.bby;
      new_movie_buffer[buffer_position].fy_2 = dyn_ob_f.bmy;
      new_movie_buffer[buffer_position].fy_3 = dyn_ob_f.ty;	
      new_movie_buffer[buffer_position].fy_4 = dyn_ob_f.umy;
      new_movie_buffer[buffer_position].fy_5 = dyn_ob_f.uby;
    }
    else if (s == BOTHBOUND) {
      Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
      bothbound_forces dyn_bb_f = dyn_bb->get_internal();
      new_movie_buffer[buffer_position].state = BOTHBOUND;
      new_movie_buffer[buffer_position].time = iteration*dt;
      new_movie_buffer[buffer_position].PE_1 = dyn_bb->PE_nba;
      new_movie_buffer[buffer_position].PE_2 = dyn_bb->PE_nma;
      new_movie_buffer[buffer_position].PE_3 = dyn_bb->PE_fma;
      new_movie_buffer[buffer_position].PE_4 = dyn_bb->PE_fba;

      new_movie_buffer[buffer_position].x_1 = dyn_bb->get_nbx();
      new_movie_buffer[buffer_position].x_2 = dyn_bb->get_nmx();
      new_movie_buffer[buffer_position].x_3 = dyn_bb->get_tx();
      new_movie_buffer[buffer_position].x_4 = dyn_bb->get_fmx();
      new_movie_buffer[buffer_position].x_5 = dyn_bb->get_fbx();
      new_movie_buffer[buffer_position].y_1 = dyn_bb->get_nby();
      new_movie_buffer[buffer_position].y_2 = dyn_bb->get_nmy();
      new_movie_buffer[buffer_position].y_3 = dyn_bb->get_ty();
      new_movie_buffer[buffer_position].y_4 = dyn_bb->get_fmy();
      new_movie_buffer[buffer_position].y_5 = dyn_bb->get_fby();

      new_movie_buffer[buffer_position].fx_1 = dyn_bb_f.nbx;
      new_movie_buffer[buffer_position].fx_2 = dyn_bb_f.nmx;
      new_movie_buffer[buffer_position].fx_3 = dyn_bb_f.tx;	
      new_movie_buffer[buffer_position].fx_4 = dyn_bb_f.fmx;
      new_movie_buffer[buffer_position].fx_5 = dyn_bb_f.fbx;
      new_movie_buffer[buffer_position].fy_1 = dyn_bb_f.nby;
      new_movie_buffer[buffer_position].fy_2 = dyn_bb_f.nmy;
      new_movie_buffer[buffer_position].fy_3 = dyn_bb_f.ty;	
      new_movie_buffer[buffer_position].fy_4 = dyn_bb_f.fmy;
      new_movie_buffer[buffer_position].fy_5 = dyn_bb_f.fby;
    }
    buffer_position++;
  }
}

void stepping_data_callback(void* dyn, State s, void *job_msg_, data_union *job_data, long long iteration) {
  check_for_quitting_conditions(iteration*dt);
  job_msg_t job_msg = *(job_msg_t *)job_msg_;

  if (iteration % 1000 == 0) {
    if (ftell(stdout) > MAX_FILESIZE_PERMITTED) {
      printf("ERROR: stdout has gotten too big.  Exiting!\n");
      exit(1);
    }
    if (job_msg.movie_data_file && ftell(job_msg.movie_data_file) > MAX_FILESIZE_PERMITTED) {
      printf("ERROR: movie_data_file has gotten too big.  Exiting!\n");
      exit(1);
    }
  }

  log_stepping_data(job_msg.stepping_data_file, dyn, iteration, job_msg.max_iteration, s);

  if ((am_making_movie or am_debugging_onebound) && iteration % int(stepping_movie_framerate/dt) == 0)
    log_stepping_movie_data(job_msg.movie_data_file, dyn, s, iteration);

  if (iteration % (long long) (0.01 / dt) == 0) { // print time every tenth of a second.
    fprintf(job_msg.stepping_data_file, "#[time: %g]\n", iteration*dt);
    printf("#[time: %g]\n", iteration*dt);
  }

  if (job_msg.max_iteration > 0) {
    if (iteration % (int)5e5 == 0 && iteration != job_msg.max_iteration && display_progress) {
      printf("Stepping data progress (%s): %lld / %lld, %g%%\n", job_msg.run_msg,
             iteration, job_msg.max_iteration, ((double) iteration) / job_msg.max_iteration * 100);
    }

    if (iteration == job_msg.max_iteration and display_progress) {
      printf("Finished generating stepping data (%s), process took %g seconds\n", job_msg.run_msg,
             ((double) clock() - job_msg.start_time) / CLOCKS_PER_SEC);
    }
  }
  else if (job_msg.max_iteration == 0 and iteration % (int)5e5 == 0) {
    //printf("Stepping data progress (%s): %.2g seconds\n", run_msg, iteration*dt);
    fflush(job_msg.stepping_data_file);
  }
}

void set_input_variables(int argc, char** argv, char* run_name, bool* am_making_movie, double* runtime) {
  char c;
  char* label;
  *run_name = 0;
  label = 0;

  static struct option long_options[] =
    {
      {"ls",       required_argument,    0, 'a'},
      {"lt",       required_argument,    0, 'b'},
      {"cb",       required_argument,    0, 'c'},
      {"cm",       required_argument,    0, 'd'},
      {"ct",       required_argument,    0, 'e'},
      {"T",        required_argument,    0, 'f'},
      {"name",     required_argument,    0, 'g'},
      {"seed",     required_argument,    0, 'h'},
      {"k_b",      required_argument,    0, 'i'},
      {"k_ub",     required_argument,    0, 'j'},
      {"k_ub_ob",  required_argument,    0, 'k'},
      {"runtime",  required_argument,    0, 'l'},
      {"label",    required_argument,    0, 'm'},
      {"dt",       required_argument,    0, 'n'},
      // {"runtime",  required_argument,    0, 'n'},
      // {"runtime",  required_argument,    0, 'o'},
      // {"runtime",  required_argument,    0, 'p'},
      // {"runtime",  required_argument,    0, 'q'},
      // {"runtime",  required_argument,    0, 'r'},
      {"nomovie",  no_argument, (int*) am_making_movie, 0},
      {"constant-write", no_argument, (int*) &am_only_writing_on_crash, false},
      {"onebound-debugging", no_argument, (int*) &am_debugging_onebound, true},
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
      low_affinity_unbinding_rate = strtod(optarg, NULL);
      break;
    // case 'k':
    //   high_affinity_unbinding_rate = strtod(optarg, NULL);
    //   break;
    case 'l':
      *runtime = strtod(optarg, NULL);
      break;
    case 'm':
      label = optarg;
      break;
    case 'n':
      dt = strtod(optarg, NULL);
      break;
    case '?':
      printf("Some other unknown getopt error.\n");
      exit(EXIT_FAILURE);
    default:
      printf("Default case in getopt: uh-oh!\n");
      exit(EXIT_FAILURE);
    }
  }

  if (*run_name == 0) {
    if (label == 0) {
      sprintf(run_name, "k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g",
	      low_affinity_binding_rate, low_affinity_unbinding_rate, cb, cm, ct, dt);
    }
    else {
      sprintf(run_name, "%s__k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g", label,
	      low_affinity_binding_rate, low_affinity_unbinding_rate, cb, cm, ct, dt);
    }
  }

  if (optind != argc) {
    printf("Improper usage. Example: ./generate_stepping_data --label test --Ls 34.5 --T 55 --movie\n");
    exit(EXIT_FAILURE);
  }
}

void sig_handler_print_movie_buffer(int signum) {
  if (am_making_movie) {
    printf("Received sigusr1, writing to data file and turning off printing.\n");
    on_crash_write_movie_buffer();
    am_making_movie = false;
    am_only_writing_on_crash = false;
  }
  else {
    printf("Turning on printing mode.\n");
    am_making_movie = true;
    am_only_writing_on_crash = true;
  }
}

int main(int argc, char** argv) {
  setvbuf(stdout, 0, _IONBF, 0);

  char* run_name = new char[100];
  am_making_movie = true;

  crash_movie_file_name_global = new char[1000];

  double runtime = 0;

  if (am_debugging_onebound) {
    printf("turning off am_only_writing_on_crash for onebound-debugging mode\n");
    am_only_writing_on_crash = false;
  }

  set_input_variables(argc, argv, run_name, &am_making_movie, &runtime);

  if (runtime == 0 and am_making_movie and not am_only_writing_on_crash) {
    printf("error,value of am_only_writing: %d\n", (int)am_only_writing_on_crash);
    //printf("Error: run settings would cause indefinite movie data printing and fill up the disc!\n");
    exit(EXIT_FAILURE);
  }

  if (am_only_writing_on_crash) {
    struct sigaction new_action;

    new_action.sa_handler = sig_handler_print_movie_buffer;
    sigemptyset (&new_action.sa_mask);
    new_action.sa_flags = 0;

    sigaction(SIGUSR1, &new_action, NULL);
  }

  char *stepping_data_fname = new char[200];
  char *stepping_config_fname = new char[200];

  char *movie_data_fname = new char[200];
  char *movie_config_fname = new char[200];

  sprintf(stepping_data_fname, "data/stepping_data_%s.txt", run_name);
  sprintf(stepping_config_fname, "data/stepping_config_%s.txt", run_name);
  sprintf(movie_data_fname, "data/stepping_movie_data_%s.txt", run_name);
  sprintf(movie_config_fname, "data/stepping_movie_config_%s.txt", run_name);

  //technically only need this if am_only_writing_on_crash is on, but do it just in case we turn it on later
  on_crash_old_movie_data_global_ptr = new movie_data_struct[MOVIE_BUFFER_SIZE];
  on_crash_new_movie_data_global_ptr = new movie_data_struct[MOVIE_BUFFER_SIZE];
  zero_movie_struct(on_crash_old_movie_data_global_ptr);
  sprintf(crash_movie_file_name_global, "data/stepping_movie_data_%s.txt", run_name);
  write_movie_config(movie_config_fname, iterations*dt);

  write_config_file(stepping_config_fname, 0, "");

  job_msg_t job_msg;
  job_msg.max_iteration = 0;
  job_msg.start_time = clock();
  job_msg.run_msg = run_name;
  job_msg.stepping_data_file = fopen(stepping_data_fname, "w");
  job_msg.movie_data_file = 0;

  if (using_variable_timestep) {
    variable_ts_stepping_data_file = job_msg.stepping_data_file;
  }

  if (am_making_movie or am_debugging_onebound) {
    job_msg.movie_data_file = fopen(movie_data_fname, "w");
    if (!job_msg.movie_data_file) {
      printf("Error opening %s!\n", movie_data_fname);
      exit(1); 
    }
    setvbuf(job_msg.movie_data_file, NULL, _IOLBF, 0); // turn on line-buffering
    fprintf(job_msg.movie_data_file, "#State\ttime\tPE_1\tPE_2\tPE_3\tPE_4\tPE_5\t"
            "x1\ty1\tx2\ty2\tx3\ty3\tx4\ty4\tx5\ty5\t"
            "fx1\tfy1\tfx2\tfy2\tfx3\tfy3\tfx4\tfy4\tfx5\tfy5\n");
  }

  fprintf(job_msg.stepping_data_file, "# command line:");
  for (int i=0; i<argc; i++) {
    fprintf(job_msg.stepping_data_file, " %s", argv[i]);
  }

  printf("\n\n\n*********%s*********\n", run_name);
  fprintf(job_msg.stepping_data_file, "\n\n\n\n#********%s********\n", run_name);
  fprintf(job_msg.stepping_data_file, "#time_unbind, time_bind, nbx, fbx\n");

  if (errno) {
    perror("Error opening stepping data or movie file.\n");
    exit(errno);
  }

  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double init_position[] = {eq.nma,
			    eq.fma,
			    0, 0, 1.0};

  simulate(runtime, RAND_INIT_SEED, BOTHBOUND, init_position, stepping_data_callback, &job_msg, NULL);

  fclose(job_msg.stepping_data_file);
  if (job_msg.movie_data_file) fclose(job_msg.movie_data_file);

  delete[] on_crash_old_movie_data_global_ptr;
  delete[] on_crash_new_movie_data_global_ptr;

  return EXIT_SUCCESS;
}
