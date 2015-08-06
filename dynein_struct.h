#include <math.h>

const double lt = 10.0;
const double ls = 10.0;

const double g = 1;

const double bla_init = (108.0 / 180) * M_PI;
const double mla_init = (36.0 / 180) * M_PI;
const double mra_init = (148.0 / 180) * M_PI;
const double bra_init = (72.0 / 180) * M_PI;

const double r_blx_init = 0;   const double r_bly_init = 0;
const double r_mlx_init = 0;   const double r_mly_init = 0;
const double  r_tx_init = 0;   const double  r_ty_init = 0;
const double r_mrx_init = 0;   const double r_mry_init = 0; 					
const double r_brx_init = 0;   const double r_bry_init = 0;					
                               
const double f_blx_init = 0;   const double f_bly_init = 0;
const double f_mlx_init = 0;   const double f_mly_init = 0;
const double  f_tx_init = 0;   const double  f_ty_init = 0;
const double f_mrx_init = 0;   const double f_mry_init = 0;
const double f_brx_init = 0;   const double f_bry_init = 0;

const double inctime = 0.1;
const double runtime = 1.0;

typedef enum {
	LEFTBOUND,
	RIGHTBOUND,
	BOTHBOUND
} states;

typedef struct
{
  double r_blx;   double r_bly;
  double r_mlx;   double r_mly;
  double r_tx;    double r_ty;
  double r_mrx;   double r_mry; 					
  double r_brx;   double r_bry;					

  double f_blx;   double f_bly;
  double f_mlx;   double f_mly;
  double f_tx;    double f_ty;
  double f_mrx;   double f_mry;
  double f_brx;   double f_bry;
} forces;


/* ******************************** DYNEIN CLASS DEFINITION *************************************** */

class Dynein {
public:
  Dynein(double bla_init, double mla_init, double mra_init, double bra_init, forces f);	
	
  void set_bla(double d);
  void set_mla(double d);
  void set_mra(double d);
  void set_bra(double d);

  void set_blx(double d);
  void set_bly(double d);

  void set_state(states s);
  
  void set_forces(forces f);

  double get_bla();
  double get_mla();
  double get_mra();
  double get_bra();

  double get_blx();
  double get_mlx();
  double get_tx();
  double get_mrx();
  double get_brx();
                  ;
  double get_bly();
  double get_mly();
  double get_ty();
  double get_mry();
  double get_bry();
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_bla();
  double get_d_mla();
  double get_d_mra();
  double get_d_bra();
  
  double get_d_blx();
  double get_d_mlx();
  double get_d_tx();
  double get_d_mrx();
  double get_d_brx();
  
  double get_d_bly();
  double get_d_mly();
  double get_d_ty();
  double get_d_mry();
  double get_d_bry();
   
  double get_f_blx();
  double get_f_mlx();
  double get_f_tx();
  double get_f_mrx();
  double get_f_brx();
  
  double get_f_bly();
  double get_f_mly();
  double get_f_ty();
  double get_f_mry();
  double get_f_bry();
  
  double get_r_blx();
  double get_r_mlx();
  double get_r_tx();
  double get_r_mrx();
  double get_r_brx();
  
  double get_r_bly();
  double get_r_mly();
  double get_r_ty();
  double get_r_mry();
  double get_r_bry();

  double get_PE();
  double get_KE();

  states get_state();

  void log(double t);
  void update_velocities();
  void read_init_file();
  
private:
  double bla;
  double mla;
  double mra;
  double bra;
  
  double blx, bly;
  
  double d_bla;   //Angular Velocities
  double d_mla;
  double d_mra;
  double d_bra;

  double r_blx;   double r_bly; //Brownian forces
  double r_mlx;   double r_mly;
  double r_tx;    double r_ty;
  double r_mrx;   double r_mry; 					
  double r_brx;   double r_bry;					

  double f_blx;   double f_bly;	//Internal Forces
  double f_mlx;   double f_mly;
  double f_tx;    double f_ty;
  double f_mrx;   double f_mry;
  double f_brx;   double f_bry;

  states state;
};

/* *********************************** UTILITY PROTOTYPES ****************************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
void   resetLog(Dynein* dyn);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
