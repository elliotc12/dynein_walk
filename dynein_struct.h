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

extern double runtime, dt, kb, T, Lt, Ls, fake_radius_t,
  fake_radius_m, fake_radius_b, water_viscosity_mu, gt, gm, gb, ct, cm,
  cb, ONEBOUND_UNBINDING_FORCE, BOTHBOUND_UNBINDING_FORCE,
  MICROTUBULE_REPULSION_FORCE, MICROTUBULE_BINDING_DISTANCE,
  RAND_INIT_SEED, binding_preexponential_factor,
  unbinding_preexponential_factor;

extern onebound_equilibrium_angles onebound_post_powerstroke_internal_angles;
extern bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles;

const bool FP_EXCEPTION_FATAL = false;

#ifdef __APPLE__    // OSX <fenv.h> does not have feenableexcept
void feenableexcept(int x);
#endif

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

  double get_bba();
  double get_bma();
  double get_uma();
  double get_uba();

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

  double get_d_bba();
  double get_d_bma();
  double get_d_uma();
  double get_d_uba();

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

  onebound_forces get_internal();
  onebound_forces get_brownian();

  double get_binding_rate();
  double get_unbinding_rate();

  State get_state();

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

  double nba, fba;
};

/* ***************************** UTILITY PROTOTYPES ****************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
double get_average(double* data,  int len);
double get_variance(double* data, int len);
void FPE_signal_handler(int signum);

/* ***************************** SIMULATION PROTOTYPES ****************************** */

typedef struct {
  int len;
  double* bb;
  double* bm;
  double* t;
  double* um;
  double* ub;
} onebound_data;

typedef struct {
  int len;
  double* nb;
  double* nm;
  double* t;
  double* fm;
  double* fb;
} bothbound_data;

typedef struct {
  int len;
  double* data;
} generic_data;

typedef union {
  onebound_data ob_data;
  bothbound_data bb_data;
  generic_data g_data;
} data_union;

void store_onebound_PEs_callback(void* dyn, State s, void* job_msg, data_union* job_data, int iteration);

void get_onebound_PE_correlation_function(generic_data* tau_data, onebound_data* corr_data, int num_corr_datapoints, int iterations, int max_tau_iter);
void get_onebound_equipartition_ratio_per_runtime(generic_data* runtime_data, onebound_data* eq_data, int num_eq_datapoints, int min_runtime_iter, int max_runtime_iter);
void get_onebound_equipartition_ratio_average_per_runtime(generic_data* runtime_data, onebound_data* eq_data, int num_eq_datapoints, int min_runtime_iter, int max_runtime_iter);

void print_data_to_file(double* data1, double* data2, int iterations, const char* legend, const char* fname);

void simulate(double runtime, double rand_seed, State init_state, double* init_position,
	      void (*job)(void* dyn, State s, void* job_msg, data_union* job_data, int iteration),
	      void* job_msg, data_union* job_data);

#endif
