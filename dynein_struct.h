#include <math.h>
#include "MersenneTwister.h"

extern double runtime;

const double dt = 1e-12; // s

const double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
const double T = 293.0; // K

const double lt = 10.0;   // nm, guess - not sure how DNA tail-bridge works
const double ls = 21.22; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms

// tail domain radius, not sure how to get since no info on DNA
// tail-bridge
const double fake_radius_t = 1.5;  // nm
// motor domain radius, derived from PyMol, motor radius 148.6
// angstroms
const double fake_radius_m = 1.48; // nm
// binding domain radius, derived from PyMol, binding radius 14.78
// angstroms
const double fake_radius_b = 0.14; // nm

// water viscosity is about 0.7 mPa s = 7e-4 Pa s = 7e-4 m^2 * kg/s / m^3
// mu = 7e-4 kg/s / m
//  ... thus this is ...
// mu = 7e-4 kg/s / (1e9 nm)
//    = 7e-13 kg/s / nm
const double water_viscosity_mu = 7e-13; // kg/(s*nm)

const double gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
const double gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
const double gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

const double ct = 0.01; // force*distance = energy = nm^2 * kg / s^2
const double cm = 0.1; // ???
const double cb = 0.5; // ???

const double ONEBOUND_UNBINDING_FORCE = 8e11;
const double BOTHBOUND_UNBINDING_FORCE = 8e11;

const double MICROTUBULE_REPULSION_FORCE = 30.0; // N/nm

const double MICROTUBULE_BINDING_DISTANCE = 0.2; // nm

const double RAND_INIT_SEED = 0;

typedef enum
{
  PRE_POWERSTROKE,
  POST_POWERSTROKE
} Mode;

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

const onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.0 * M_PI,
  0.6 * M_PI
};

const bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.6 * M_PI,
  1.4 * M_PI,
  0.4 * M_PI
};

/* ******************** ONEBOUND DYNEIN CLASS DEFINITION ********************** */

class Dynein_onebound {
public:
  Dynein_onebound(double bba_init, double bma_init, double fma_init, double fba_init,
                  double bbx_init, double bby_init, State s,
		  onebound_forces *internal_test,
                  onebound_forces *brownian_test,
		  onebound_equilibrium_angles* eq_angles);

  Dynein_onebound(Dynein_bothbound* old_dynein, MTRand* rand, State s);

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

  double get_PE();
  double get_KE();

  State get_state();

  void log(float t, FILE* data_file);
  void log_run(float runtime, FILE* data_file);

  void update_velocities();

private:
  void update_brownian_forces();
  void update_internal_forces();

  onebound_equilibrium_angles eq;      //Equilibrium angles
  MTRand *rand;

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

  Mode mode;
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
		   bothbound_equilibrium_angles* eq_angles);

  Dynein_bothbound(Dynein_onebound* old_dynein, MTRand* rand);

  void set_nma(double d);
  void set_fma(double d);
  void set_L(double d);

  double get_nma(); // actual bothbound coordinates
  double get_fma();

  double get_nba(); // utility fns, calculated from nma, fma
  double get_fba();

  double get_nbx();
  double get_nmx();
  double get_tx();
  double get_fmx();
  double get_fbx();

  double get_nby();
  double get_nmy();
  double get_ty();
  double get_fmy();
  double get_fby();

  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_nma();
  double get_d_fma();

  double get_d_nbx();
  double get_d_nmx();
  double get_d_tx();
  double get_d_fmx();
  double get_d_fbx();

  double get_d_nby();
  double get_d_nmy();
  double get_d_ty();
  double get_d_fmy();
  double get_d_fby();

  bothbound_forces get_internal();
  bothbound_forces get_brownian();

  double get_near_unbinding_rate();
  double get_far_unbinding_rate();

  double get_PE();
  double get_KE();

  void log(float t, FILE* data_file);
  void log_run(float runtime, FILE* data_file);

  void update_velocities();

private:
  void update_brownian_forces();
  void update_internal_forces();

  bothbound_equilibrium_angles eq;  //Equilibrium angles

  double nma, fma;  //Bothbound coordinates
  double nbx, nby;
  double L;

  double d_ln;  //Bothbound velocities
  double d_lf;

  bothbound_forces r; //Brownian forces
  bothbound_forces f; //Internal Forces

  Mode mode;
  bothbound_forces brownian_testcase;
  bothbound_forces internal_testcase;
};

/* ***************************** UTILITY PROTOTYPES ****************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
double resetLog();
