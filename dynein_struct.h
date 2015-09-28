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

const double ct = 1e-3; // force*distance = energy = nm^2 * kg / s^2
const double cm = 0.1; // ???
const double cb = 0.5; // ???

const double UNBINDING_FORCE = 8e11;

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
  double fmx;   double fmy;
  double fbx;   double fby;
} forces;

typedef struct
{
  double bba, ba, ta, fa;
} equilibrium_angles;

const equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  (108.0 / 180) * M_PI,
  (108.0 / 180) * M_PI,
  0,
  (108.0 / 180) * M_PI
};

const equilibrium_angles near_farbound_post_powerstroke_internal_angles = {
  (108.0 / 180) * M_PI,
  (108.0 / 180) * M_PI,
  (52.0 / 180) * M_PI,
  (138.0 / 180) * M_PI
};

/* ******************************** DYNEIN CLASS DEFINITION ************************************* */

class Dynein {
public:
  Dynein(double bla_init, double mla_init, double mra_init, double bra_init,
         State s, forces* internal_test, forces* brownian_test, equilibrium_angles* eq_angles);

  void set_bba(double d);
  void set_bma(double d);
  void set_fma(double d);
  void set_fba(double d);

  void set_bbx(double d);
  void set_bby(double d);

  void set_state(State s);

  double get_bba();
  double get_bma();
  double get_fma();
  double get_fba();

  double get_bbx();
  double get_bmx();
  double get_tx();
  double get_fmx();
  double get_fbx();
                  ;
  double get_bby();
  double get_bmy();
  double get_ty();
  double get_fmy();
  double get_fby();
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_bba();
  double get_d_bma();
  double get_d_fma();
  double get_d_fba();
  
  double get_d_bbx();
  double get_d_bmx();
  double get_d_tx();
  double get_d_fmx();
  double get_d_fbx();
  
  double get_d_bby();
  double get_d_bmy();
  double get_d_ty();
  double get_d_fmy();
  double get_d_fby();
   
  forces get_internal();
  forces get_brownian();

  void switch_to_bothbound();
  void unbind();

  double get_binding_probability();
  double get_unbinding_probability();

  double get_PE();
  double get_KE();

  MTRand* rand;

  State get_state();

  void log(double t, FILE* data_file);
  void resetLog();
  void update_velocities();
  
private:
  void update_brownian_forces();
  void update_internal_forces();

  void update_velocities_onebound();
  void update_velocities_bothbound();

  equilibrium_angles eq;

  double bba;
  double bma;
  double fma;
  double fba;
  
  double bbx, bby;
  
  double d_bba;   //Angular Velocities
  double d_bma;
  double d_fma;
  double d_fba;

  forces r; //Brownian forces
  forces f; //Internal Forces

  Mode mode;
  forces *brownian_testcase;
  forces *internal_testcase;
  State state;
};

/* *********************************** UTILITY PROTOTYPES ****************************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
