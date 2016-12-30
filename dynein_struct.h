#include <csignal>
#include <fenv.h>
#include <math.h>
#include <stdbool.h>

#include "MersenneTwister.h"

#ifndef DYNEIN_STRUCT_H
#define DYNEIN_STRUCT_H

typedef enum
{
  NEARBOUND,
  FARBOUND,
  BOTHBOUND,
  UNBOUND
} State;

typedef struct
{
  double bbx;   double bby;
  double bmx;   double bmy;
  double tx;    double ty;
  double umx;   double umy;
  double ubx;   double uby;
} onebound_forces;

typedef struct
{
  double nbx;   double nby;
  double nmx;   double nmy;
  double tx;    double ty;
  double fmx;   double fmy;
  double fbx;   double fby;
} bothbound_forces;

typedef struct
{
  double bba, bma, ta, uma;
} onebound_equilibrium_angles;

typedef struct
{
  double nba, nma, ta, fma, fba;
} bothbound_equilibrium_angles;

typedef struct {
  State state;
  double time;
  double PE_1, PE_2, PE_3, PE_4, PE_5;
  double x_1, x_2, x_3, x_4, x_5;
  double y_1, y_2, y_3, y_4, y_5;
  double fx_1, fx_2, fx_3, fx_4, fx_5;
  double fy_1, fy_2, fy_3, fy_4, fy_5;
} movie_data_struct;

extern double runtime, dt, kb, T, Lt, Ls, fake_radius_t,
  fake_radius_m, fake_radius_b, water_viscosity_mu, gt, gm, gb, ct, cm,
  cb, ONEBOUND_UNBINDING_FORCE, BOTHBOUND_UNBINDING_FORCE,
  MICROTUBULE_REPULSION_FORCE, MICROTUBULE_BINDING_DISTANCE,
  RAND_INIT_SEED, binding_preexponential_factor, D, tau,
  low_affinity_unbinding_rate,
  low_affinity_binding_rate,
  binding_fraction,
  DELTA_G_FORMATION_BINDING;

extern onebound_equilibrium_angles onebound_post_powerstroke_internal_angles;
extern bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles;

const bool FP_EXCEPTION_FATAL = false;
const bool am_debugging_conversions = false;
const bool am_debugging_angles = false;
const bool am_debugging_time = false;
const bool am_debugging_state_transitions = false;
const bool am_debugging_rates = true;

const bool am_only_writing_on_crash = true;

const bool crash_on_nan = true;

#ifdef __APPLE__    // OSX <fenv.h> does not have feenableexcept
void feenableexcept(int x);
#endif

const int msync_after_num_writes = 10;

class DynArr {
private:
  int len;
  int current;
  double* data;
public:
  DynArr(int init_len);
  ~DynArr();
  void append(double data);
  double* get_data();
  int get_length() {return current;}
};

/* ******************** ONEBOUND DYNEIN CLASS DEFINITION ********************** */

class Dynein_bothbound;

class Dynein_onebound {
public:
  Dynein_onebound(double bba_init, double bma_init, double fma_init, double fba_init,
                  double bbx_init, double bby_init, State s,
		  onebound_forces *internal_test,
                  onebound_forces *brownian_test,
		  onebound_equilibrium_angles* eq_angles,
                  MTRand* mtrand);

  Dynein_onebound(Dynein_bothbound* old_dynein, MTRand* rand, State s);

  MTRand *rand;

  /** Onebound functions **/
  void set_bba(double d);
  void set_bma(double d);
  void set_uma(double d);
  void set_uba(double d);

  void set_bbx(double d);
  void set_bby(double d);

  double get_bba() {return bba;}
  double get_bma() {return bma;}
  double get_uma() {return uma;}
  double get_uba() {return uba;}

  double get_bbx();
  double get_bmx();
  double get_tx();
  double get_umx();
  double get_ubx();

  double get_bby();
  double get_bmy();
  double get_ty();
  double get_umy();
  double get_uby();

  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_bba() {return d_bba;}
  double get_d_bma() {return d_bma;}
  double get_d_uma() {return d_uma;}
  double get_d_uba() {return d_uba;}

  double get_d_bbx();
  double get_d_bmx();
  double get_d_tx();
  double get_d_umx();
  double get_d_ubx();

  double get_d_bby();
  double get_d_bmy();
  double get_d_ty();
  double get_d_umy();
  double get_d_uby();

  onebound_forces get_internal() {return f;}
  onebound_forces get_brownian() {return r;}

  double get_binding_rate();
  double get_unbinding_rate();

  State get_state() {return state;}

  void update_velocities();

  double PE_bba, PE_bma, PE_ta, PE_uma;
  double get_PE() { return PE_bba + PE_bma + PE_ta + PE_uma; }

private:
  void update_brownian_forces();
  void update_internal_forces();

  onebound_equilibrium_angles eq;      //Equilibrium angles

  double bba;    //Onebound coordinates
  double bma;
  double uma;
  double uba;

  double bbx, bby;

  double d_bba;   //Onebound angular velocities
  double d_bma;
  double d_uma;
  double d_uba;

  onebound_forces r; //Brownian forces
  onebound_forces f; //Internal Forces

  onebound_forces *brownian_testcase;
  onebound_forces *internal_testcase;
  State state;
};

/* ******************* BOTHBOUND DYNEIN CLASS DEFINITION ********************** */

class Dynein_bothbound {
public:
  Dynein_bothbound(double nma_init, double fma_init, double nbx_init,
		   double nby_init, double L_init,
		   bothbound_forces* internal_test, bothbound_forces* brownian_test,
		   bothbound_equilibrium_angles* eq_angles,
		   MTRand* mtrand);

  Dynein_bothbound(Dynein_onebound* old_dynein, MTRand* rand);

  void set_nma(double d);
  void set_fma(double d);
  void set_L(double d);

  void set_dLn(double d);
  void set_dLf(double d);

  double get_nma() { return nma; }   // actual bothbound coordinates
  double get_fma() { return fma; }

  double get_nba() { return nba; }
  double get_ta()  { return ta; }
  double get_fba() { return fba; }

  double get_nbx() { return nbx; }
  double get_nmx() { return nmx; }
  double get_tx() { return tx; }
  double get_fmx() { return fmx; }
  double get_fbx() { return nbx + L; }

  double get_nby() { return nby; }
  double get_nmy() { return nmy; }
  double get_ty() { return ty; }
  double get_fmy() { return fmy; }
  double get_fby() { return nby; }

  double get_ln() { return Ln; }
  double get_lf() { return Lf; }

  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_nma();
  double get_d_fma();

  /* double get_d_nba(); // computed values */
  /* double get_d_fba(); */

  double get_d_nbx() {return 0;}
  double get_d_nmx() {return dXnm_dLn*d_Ln + dXnm_dLf*d_Lf;}
  double get_d_tx()  {return dXt_dLn*d_Ln + dXt_dLf*d_Lf;}
  double get_d_fmx() {return dXfm_dLn*d_Ln + dXfm_dLf*d_Lf;}
  double get_d_fbx() {return 0;}

  double get_d_nby() {return 0;}
  double get_d_nmy() {return dYnm_dLn*d_Ln + dYnm_dLf*d_Lf;}
  double get_d_ty()  {return dYt_dLn*d_Ln + dYt_dLf*d_Lf;}
  double get_d_fmy() {return dYfm_dLn*d_Ln + dYfm_dLf*d_Lf;}
  double get_d_fby() {return 0;}

  double get_d_Ln() { return d_Ln; };
  double get_d_Lf() { return d_Lf; };

  bothbound_forces get_internal();
  bothbound_forces get_brownian();

  double get_near_unbinding_rate();
  double get_far_unbinding_rate();

  double get_PE();
  double PE_nba, PE_nma, PE_ta, PE_fma, PE_fba;

  void update_coordinates();
  void update_velocities();

  //Variables which should be internal, but we need them public for test_bothbound
     // Various distances and angles useful in computations (see paper)
  double Ln, Lf;
  double cosAn, sinAn, cosAns, sinAns;
  double cosAf, sinAf, cosAfs, sinAfs;
  double nmx, fmx, tx;
  double nmy, fmy, ty;

  double dcosAn_dLn, dsinAn_dLn, dcosAns_dLn, dsinAns_dLn;
  double dcosAf_dLn, dsinAf_dLn, dcosAfs_dLn, dsinAfs_dLn;

  double dcosAn_dLf, dsinAn_dLf, dcosAns_dLf, dsinAns_dLf;
  double dcosAf_dLf, dsinAf_dLf, dcosAfs_dLf, dsinAfs_dLf;
  
     // Various interesting derivatives that are used in finding the
     // velocities (and are set by update_velocities).
  double dXnm_dLn, dYnm_dLn, dXnm_dLf, dYnm_dLf, dXfm_dLf, dYfm_dLf,
    dXfm_dLn, dYfm_dLn, dXt_dLn, dYt_dLn, dXt_dLf, dYt_dLf;

private:
  void update_brownian_forces();
  void update_internal_forces();

  bothbound_equilibrium_angles eq;  //Equilibrium angles
  MTRand *rand;

  double nma, fma;  //Bothbound coordinates
  double nbx, nby;
  double L;

  double d_Ln;  //Bothbound velocities
  double d_Lf;

  bothbound_forces r; //Brownian forces
  bothbound_forces f; //Internal Forces

  bothbound_forces *brownian_testcase;
  bothbound_forces *internal_testcase;

  double nba, ta, fba;
};

/* ***************************** UTILITY PROTOTYPES ****************************** */
enum WRITE_CONFIG_OMIT_FLAGS {
  CONFIG_OMIT_T = 1,
  CONFIG_OMIT_C = 2,
  CONFIG_INCLUDE_SKIPINFO = 4
};

double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
double get_average(double* data,  int len);
double get_variance(double* data, int len);
void FPE_signal_handler(int signum);
void prepare_data_file(const char* legend, char* fname);
void append_data_to_file(double* data1, double* data2, int len, FILE* fd);
void write_config_file(char* fname, int omit_flags, const char* custom_str);
void write_movie_config(char* movie_config_fname, double runtime);

/* ***************************** SIMULATION PROTOTYPES ****************************** */

typedef struct {
  long long len;
  double* bb;
  double* bm;
  double* t;
  double* um;
  double* ub;
} onebound_data;

typedef struct {
  long long len;
  double* nb;
  double* nm;
  double* t;
  double* fm;
  double* fb;
} bothbound_data;

typedef struct {
  long long len;
  void* data;
} generic_data;

typedef union {
  onebound_data ob_data;
  bothbound_data bb_data;
  generic_data g_data;
} data_union;

typedef struct {
  double* time;
  double* bb;
  double* bm;
  double* t;
  double* um;
} pe_callback_data_ptr;

typedef struct {
  double* bb;
  double* bm;
  double* t;
  double* um;
  double* f_bbx;   double* f_bby;
  double* f_bmx;   double* f_bmy;
  double* f_tx;    double* f_ty;
  double* f_umx;   double* f_umy;
  double* f_ubx;   double* f_uby;
} eq_forces_callback_data_ptr;

typedef struct {
  double time;
  double bba;
  double bma;
  double ta;
  double uma;
  double bba_PE;
  double bma_PE;
  double ta_PE;
  double uma_PE;
  double f_bbx;   double f_bby;
  double f_bmx;   double f_bmy;
  double f_tx;    double f_ty;
  double f_umx;   double f_umy;
  double f_ubx;   double f_uby;
  double bbx;   double bby;
  double bmx;   double bmy;
  double tx;    double ty;
  double umx;   double umy;
  double ubx;   double uby;
} onebound_data_generate_struct;

typedef struct {
  double time;
  double nba;
  double nma;
  double ta;
  double fma;
  double fba;
  double nba_PE;
  double nma_PE;
  double ta_PE;
  double fma_PE;
  double fba_PE;
  double f_nbx;   double f_nby;
  double f_nmx;   double f_nmy;
  double f_tx;    double f_ty;
  double f_fmx;   double f_fmy;
  double f_fbx;   double f_fby;
  double nbx;   double nby;
  double nmx;   double nmy;
  double tx;    double ty;
  double fmx;   double fmy;
  double fbx;   double fby;
} bothbound_data_generate_struct;

void store_onebound_PEs_callback(void* dyn, State s, void** job_msg, data_union* job_data, long long iteration);
void write_onebound_PEs_callback(void* dyn, State s, void** job_msg, data_union* job_data, long long iteration);
void store_onebound_PEs_and_forces_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration);

void get_onebound_PE_correlation_function(generic_data* tau_data, onebound_data* corr_data, long long num_correlations, long long iterations, long long max_tau_iter, const int* seeds, int seed_len, char* run_msg_base);
void get_onebound_equipartition_ratio_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long num_runtimes, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len, char* run_msg_base);
void get_onebound_equipartition_ratio_average_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long d_runtime_iter, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len, char* run_msg_base);
void get_onebound_equipartition_ratio(onebound_data* eq_data, generic_data* force_data, long long runtime_iter, const int* seeds, int seed_len, char* run_msg_base);

void simulate(double runtime, double rand_seed, State init_state, double* init_position,
	      void (*job)(void* dyn, State s, void *job_msg, data_union* job_data,
              long long iteration), void *job_msg, data_union* job_data);

/** create_ob/bb_plots .txt generation code **/
void generate_force_data(double* times, double* f, int len, const char* legend, char* fname_base, const char* annotation);
void generate_correlation_fn_data(double* pe, int iters, const char* legend, char* fname_base);
void generate_pe_vs_time_data(double* times, double* pe, int len, const char* legend, char* fname_base);
void generate_ave_pe_and_log_error_data(double* times, double* pe, int iters, const char* legend, char* fname_base);
void generate_angle_vs_time_data(double* times, double* angle, int len, const char* legend, char* fname_base, double eq_angle);

void on_crash_write_movie_buffer();

#endif
