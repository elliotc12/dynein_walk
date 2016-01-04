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

/* ***************************** UTILITY PROTOTYPES ****************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
