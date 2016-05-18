#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

inline double sqr(double x) {
  return x*x;
}

/* ********************* BOTHBOUND DYNEIN FUNCTIONS ****************************** */

Dynein_bothbound::Dynein_bothbound(double nma_init, double fma_init, double nbx_init,
                                   double nby_init, double L_init,
                                   bothbound_forces* internal_test,
                                   bothbound_forces* brownian_test,
                                   bothbound_equilibrium_angles* eq_angles,
				   MTRand* mtrand) {
  nbx = nbx_init;
  nby = nby_init;

  nma = nma_init;
  fma = fma_init;

  L = L_init;

  internal_testcase = internal_test;
  brownian_testcase = brownian_test;

  if (eq_angles) {
    eq = *eq_angles; // use test angles
  } else {
    eq = bothbound_pre_powerstroke_internal_angles; // use experimental angles
  }

  rand = mtrand;

  update_velocities();
}

Dynein_bothbound::Dynein_bothbound(Dynein_onebound* old_dynein, MTRand* mtrand) {
  // out of old dyn
  if (old_dynein->get_state() == State::NEARBOUND) {
    nbx = old_dynein->get_bbx();
    nby = 0;

    nma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();
    fma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();

    L = old_dynein->get_ubx() - old_dynein->get_bbx();

  } else {
    nbx = old_dynein->get_ubx();
    nby = 0;

    nma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();
    fma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();

    L = old_dynein->get_bbx() - old_dynein->get_ubx();
  }

  rand = mtrand;

  update_velocities();
}

void Dynein_bothbound::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.nbx, &r.nby);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.nmx, &r.nmy);
    rand->gauss2(sqrt(2*kb*T/(gt*dt)), &r.tx, &r.ty);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.fmx, &r.fmy);
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.fbx, &r.fby);
  }
}

void Dynein_bothbound::update_internal_forces() {
  if (internal_testcase) {
    f = *internal_testcase;
  } else {
    f.nbx = 0;     f.nby = 0;     // Initialize forces to zero
    f.nmx = 0;     f.nmy = 0;
    f.tx  = 0;     f.ty  = 0;
    f.fmx = 0;     f.fmy = 0;
    f.fbx = 0;     f.fby = 0;

    double T, f1, f2, f1x, f2x, f1y, f2y;

    T = cb*(nba - eq.nba);
    f2 = T/ls;
    f2x = f2 * sin(nba);
    f2y = f2 * -cos(nba);
    f.nmx += f2x;
    f.nmy += f2y;
    f.nbx += -f2x; // Equal and opposite forces!  :)
    f.nby += -f2y; // Equal and opposite forces!  :)

    T = cm*(nma - eq.nma);
    f1 = T/ls;
    f2 = T/lt;
    f1x = f1 * sin(nba);
    f1y = f1 * -cos(nba);
    f2x = f2 * sin(nma - (M_PI - nba));
    f2y = f2 * -cos(nma - (M_PI - nba));
    f.nbx += f1x;
    f.nby += f1y;
    f.tx  += f2x;
    f.ty  += f2y;
    f.nmx += -(f1x + f2x);
    f.nmy += -(f1y + f2y);

    T = ct*(fma + fba - nma - nba - eq.ta);
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(nma + nba - M_PI);
    f1y = f1 * -cos(nma + nba - M_PI);
    f2x = f2 * -sin(fma + fba - M_PI);
    f2y = f2 * cos(fma + fba - M_PI);
    f.nmx += f1x;
    f.nmy += f1y;
    f.fmx += -f2x; // flip sign of this to fake correct force
    // (force is calculated properly, but motion equations don't take it into acct)
    f.fmy += -f2y; // flip sign of this to fake correct force
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*(fma - eq.fma);
    f1 = T / lt;
    f2 = T / ls;
    f1x = f1 * sin(fma + fba - M_PI);
    f1y = f1 * -cos(fma + fba - M_PI);
    f2x = f2 * sin(fba);
    f2y = f2 * -cos(fba);
    f.tx  += f1x;
    f.ty  += f1y;
    f.fbx += f2x;
    f.fby += f2y;
    f.fmx += (f1x + f2x); // flip sign of this to fake correct force
    f.fmy += (f1y + f2y);

    T = cb*(fba - eq.fba);
    f1 = T / ls;
    f1x = f1 * sin(fba);
    f1y = f1 * -cos(fba);
    f.fmx += -f1x; // flip sign of this to fake correct force
    f.fmy += -f1y; // flip sign of this to fake correct force
    f.fbx += -f1x;
    f.fby += -f1y;

    if (get_nmy() < 0) f.nmy += MICROTUBULE_REPULSION_FORCE * fabs(get_nmy());
    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
    if (get_fmy() < 0) f.fmy += MICROTUBULE_REPULSION_FORCE * fabs(get_fmy());
  }
}

void Dynein_bothbound::update_coordinates() {
  Ln = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(nma));
  Lf = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(fma));

  cosAn = (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln);
  sinAn = sqrt(1 - cosAn*cosAn);
  cosAns = (Ls*Ls + Ln*Ln - Lt*Lt) / (2*Ls*Ln);
  sinAns = (nma < M_PI) ? sqrt(1-cosAns*cosAns) : -sqrt(1-cosAns*cosAns);

  cosAf = -(L*L + Lf*Lf - Ln*Ln) / (2*L*Lf);
  sinAf = sqrt(1 - cosAf*cosAf);
  cosAfs = (Ls*Ls + Lf*Lf - Lt*Lt) / (2*Ls*Lf);
  sinAfs = (fma > M_PI) ? sqrt(1-cosAfs*cosAfs) : -sqrt(1-cosAfs*cosAfs);

  nmx = Ls*(cosAn*cosAns - sinAn*sinAns);
  nmy = Ls*(cosAn*sinAns + sinAn*cosAns);
  tx = Ln*cosAn;
  ty = Ln*sinAn;
  fmx = L + Ls*(cosAf*cosAfs - sinAf*sinAfs);
  fmy = Ls*(cosAf*sinAfs + sinAf*cosAfs);

  // angle of stalks from horizontal
  nba = atan2(nmy, nmx - nbx);
  fba = atan2(fmy, fmx - (nbx + L));

  assert(nma != M_PI);
  assert(fma != M_PI);
}

static const bool am_debugging_nans = false;

void Dynein_bothbound::update_velocities() {
  update_coordinates();
  update_internal_forces();
  update_brownian_forces();

  if (am_debugging_nans) printf("cosAn %g\n", cosAn);
  if (am_debugging_nans) printf("sinAn %g\n", sinAn);
  if (am_debugging_nans) printf("cosAns %g\n", cosAns);
  if (am_debugging_nans) printf("sinAns %g\n", sinAns);

  bothbound_forces rforces = r;
  bothbound_forces fforces = f;

  double dcosAn_dLn = (1/L) - (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln*Ln);
  if (am_debugging_nans) printf("dcosAn_dLn is %g\n", dcosAn_dLn);
  double dcosAn_dLf = -(Lf) / (L*Ln);
  double dsinAn_dLn = -cosAn / sqrt(1 - cosAn*cosAn) * (1/L - (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln*Ln));
  double dsinAn_dLf = cosAn / sqrt(1 - cosAn*cosAn) * Lf / (L*Ln);
  double dcosAns_dLn = 1/Ls - (Ls*Ls + Ln*Ln - Lt*Lt) / (2*Ls*Ln*Ln);
  double dcosAns_dLf = 0;
  double dsinAns_dLn = (nma < M_PI) ?
     -cosAns/sqrt(1-cosAns*cosAns) * (1/Ls - (Ls*Ls+Ln*Ln-Lt*Lt)/(2*L*Ln*Ln))
    : cosAns/sqrt(1-cosAns*cosAns) * (1/Ls - (Ls*Ls+Ln*Ln-Lt*Lt)/(2*L*Ln*Ln));
  double dsinAns_dLf = 0;

  if (am_debugging_nans) printf("dcosAn_dLn %g\n", dcosAn_dLn);
  if (am_debugging_nans) printf("dsinAn_dLn %g\n", dsinAn_dLn);
  if (am_debugging_nans) printf("dcosAns_dLn %g\n", dcosAns_dLn);
  if (am_debugging_nans) printf("dsinAns_dLn %g\n", dsinAns_dLn);

  double dcosAf_dLf = (1/L) - (L*L + Lf*Lf - Ln*Ln) / (2*L*Lf*Lf);
  double dcosAf_dLn = -(Ln) / (L*Lf);
  double dsinAf_dLf = -cosAf / sqrt(1 - cosAf*cosAf) * (1/L - (L*L + Lf*Lf - Ln*Ln) / (2*L*Lf*Lf));
  double dsinAf_dLn = cosAf / sqrt(1 - cosAf*cosAf) * Ln / (L*Lf);
  double dcosAfs_dLf = 1/Ls - (Ls*Ls + Lf*Lf - Lt*Lt) / (2*Ls*Lf*Lf);
  double dcosAfs_dLn = 0;
  double dsinAfs_dLf = (fma < M_PI) ?
     -cosAfs/sqrt(1-cosAfs*cosAfs) * (1/Ls - (Ls*Ls+Lf*Lf-Lt*Lt)/(2*L*Lf*Lf))
    : cosAfs/sqrt(1-cosAfs*cosAfs) * (1/Ls - (Ls*Ls+Lf*Lf-Lt*Lt)/(2*L*Lf*Lf));
  double dsinAfs_dLn = 0;

  double dXnm_dLn = Ls*(cosAn * dcosAn_dLn + cosAns * dcosAn_dLn
			- sinAn * dsinAns_dLn - sinAns * dsinAn_dLn);
  double dYnm_dLn = Ls*(cosAn * dsinAns_dLn + sinAns * dcosAn_dLn
			+ sinAn * dcosAns_dLn + cosAns * dsinAn_dLn);
  double dXnm_dLf = Ls*(cosAn * dcosAns_dLf + cosAns * dcosAn_dLf
			- sinAn * dsinAns_dLf - sinAns * dsinAn_dLf);
  double dYnm_dLf = Ls*(cosAn * dsinAns_dLf + sinAns * dcosAn_dLf
			+ sinAn * dcosAns_dLf + cosAns * dsinAn_dLf);

  double dXfm_dLf = Ls*(cosAf * dcosAf_dLf + cosAfs * dcosAf_dLf
			- sinAf * dsinAfs_dLf - sinAfs * dsinAf_dLf);
  double dYfm_dLf = Ls*(cosAf * dsinAfs_dLf + sinAfs * dcosAf_dLf
			+ sinAf * dcosAfs_dLf + cosAfs * dsinAf_dLf);
  double dXfm_dLn = Ls*(cosAf * dcosAfs_dLn + cosAfs * dcosAf_dLn
			- sinAf * dsinAfs_dLn - sinAfs * dsinAf_dLn);
  double dYfm_dLn = Ls*(cosAf * dsinAfs_dLn + sinAfs * dcosAf_dLn
			+ sinAf * dcosAfs_dLn + cosAfs * dsinAf_dLn);

  double dXt_dLn = cosAn;
  double dYt_dLn = sinAn;
  double dXt_dLf = 0;
  double dYt_dLf = 0;
  if (am_debugging_nans) printf("dXt_dLn is %g\n", dXt_dLn);
  if (am_debugging_nans) printf("dYt_dLn is %g\n", dYt_dLn);

  if (am_debugging_nans) printf("gm = %g, gt = %g\n", gm, gt);

  double a = -dXnm_dLn;
  double b = -dXnm_dLf;
  double c = -(get_nmx() - get_nbx()) / gm;
  double d = (get_tx() - get_nmx()) / gm;
  double e = -dXt_dLn;
  double f = -dXt_dLf;
  double g = -(get_tx() - get_nmx()) / gt;
  double h = (get_fmx() - get_tx()) / gt;
  double i = -dXfm_dLn;
  double j = -dXfm_dLf;
  double k = -(get_fmx() - get_tx()) / gm;
  double l = (get_fbx() - get_fmx()) / gm;
  double m = -dYnm_dLn;
  double n = -dYnm_dLf;
  double p = -(get_nmy() - get_nby()) / gm;
  double q = (get_ty() - get_nmy()) / gm;
  double r = -dYt_dLn;
  double s = -dYt_dLf;
  double t = -(get_ty() - get_nmy()) / gt;
  double u = (get_fmy() - get_ty()) / gt;
  double v = -dYfm_dLn;
  double w = -dYfm_dLf;
  double x = -(get_fmy() - get_ty()) / gm;
  double y = (get_fby() - get_fmy()) / gm;
  if (am_debugging_nans) {
    printf("a=%g\tb=%g\tc=%g\td=%g\te=%g\tf=%g\n",
           a,b,c,d,e,f);
    printf("g=%g\th=%g\ti=%g\tj=%g\tk=%g\tl=%g\n",
           g,h,i,j,k,l);
    printf("m=%g\tn=%g\tp=%g\tq=%g\tr=%g\ts=%g\n",
           m,n,p,q,r,s);
    printf("t=%g\tu=%g\tv=%g\tw=%g\tx=%g\ty=%g\n",
           t,u,v,w,x,y);
  }

  double x1 = -fforces.nmx / g - rforces.nmx;
  double x2 = -fforces.tx / g  - rforces.tx;
  double x3 = -fforces.fmx / g - rforces.fmx;
  double x4 = -fforces.nmy / g - rforces.nmy;
  double x5 = -fforces.ty / g  - rforces.ty;
  double x6 = -fforces.fmy / g - rforces.fmy;

  d_Ln =
    (h*l*p*t*w*x1 - g*l*p*u*w*x1 + g*l*p*s*x*x1 - f*l*p*t*x*x1 + d*l*p*u*w*x2 - c*l*q*u*w*x2 - d*l*p*s*x*x2 + 
     c*l*q*s*x*x2 - c*l*n*t*x*x2 + b*l*p*t*x*x2 - c*h*l*t*w*x4 + c*g*l*u*w*x4 - c*g*l*s*x*x4 + c*f*l*t*x*x4 - 
     d*h*l*p*w*x5 + c*h*l*q*w*x5 + c*g*l*n*x*x5 + d*f*l*p*x*x5 - b*g*l*p*x*x5 - c*f*l*q*x*x5 + d*h*l*p*s*x6 - 
     c*h*l*q*s*x6 + c*h*l*n*t*x6 - b*h*l*p*t*x6 - c*g*l*n*u*x6 - d*f*l*p*u*x6 + b*g*l*p*u*x6 + c*f*l*q*u*x6 - 
     g*k*p*s*x1*y - h*j*p*t*x1*y + f*k*p*t*x1*y + g*j*p*u*x1*y + d*k*p*s*x2*y - c*k*q*s*x2*y + c*k*n*t*x2*y - 
     b*k*p*t*x2*y - d*j*p*u*x2*y + c*j*q*u*x2*y - d*h*p*s*x3*y + c*h*q*s*x3*y - c*h*n*t*x3*y + b*h*p*t*x3*y + 
     c*g*n*u*x3*y + d*f*p*u*x3*y - b*g*p*u*x3*y - c*f*q*u*x3*y + c*g*k*s*x4*y + c*h*j*t*x4*y - c*f*k*t*x4*y - 
     c*g*j*u*x4*y - c*g*k*n*x5*y + d*h*j*p*x5*y - d*f*k*p*x5*y + b*g*k*p*x5*y - c*h*j*q*x5*y + c*f*k*q*x5*y)/
   (d*h*l*p*s*v - c*h*l*q*s*v + c*h*l*n*t*v - b*h*l*p*t*v - c*g*l*n*u*v - d*f*l*p*u*v + b*g*l*p*u*v + 
     c*f*l*q*u*v - d*h*l*p*r*w + c*h*l*q*r*w - c*h*l*m*t*w + a*h*l*p*t*w + c*g*l*m*u*w + d*e*l*p*u*w - 
     a*g*l*p*u*w - c*e*l*q*u*w + c*g*l*n*r*x + d*f*l*p*r*x - b*g*l*p*r*x - c*f*l*q*r*x - c*g*l*m*s*x - 
     d*e*l*p*s*x + a*g*l*p*s*x + c*e*l*q*s*x + c*f*l*m*t*x - c*e*l*n*t*x + b*e*l*p*t*x - a*f*l*p*t*x - 
     c*g*k*n*r*y + d*h*j*p*r*y - d*f*k*p*r*y + b*g*k*p*r*y - c*h*j*q*r*y + c*f*k*q*r*y + c*g*k*m*s*y - 
     d*h*i*p*s*y + d*e*k*p*s*y - a*g*k*p*s*y + c*h*i*q*s*y - c*e*k*q*s*y + c*h*j*m*t*y - c*f*k*m*t*y - 
     c*h*i*n*t*y + c*e*k*n*t*y + b*h*i*p*t*y - a*h*j*p*t*y - b*e*k*p*t*y + a*f*k*p*t*y - c*g*j*m*u*y + 
     c*g*i*n*u*y + d*f*i*p*u*y - b*g*i*p*u*y - d*e*j*p*u*y + a*g*j*p*u*y - c*f*i*q*u*y + c*e*j*q*u*y);

  d_Lf =
    (-(h*l*p*t*v*x1) + g*l*p*u*v*x1 - g*l*p*r*x*x1 + e*l*p*t*x*x1 - d*l*p*u*v*x2 + c*l*q*u*v*x2 + d*l*p*r*x*x2 - 
     c*l*q*r*x*x2 + c*l*m*t*x*x2 - a*l*p*t*x*x2 + c*h*l*t*v*x4 - c*g*l*u*v*x4 + c*g*l*r*x*x4 - c*e*l*t*x*x4 + 
     d*h*l*p*v*x5 - c*h*l*q*v*x5 - c*g*l*m*x*x5 - d*e*l*p*x*x5 + a*g*l*p*x*x5 + c*e*l*q*x*x5 - d*h*l*p*r*x6 + 
     c*h*l*q*r*x6 - c*h*l*m*t*x6 + a*h*l*p*t*x6 + c*g*l*m*u*x6 + d*e*l*p*u*x6 - a*g*l*p*u*x6 - c*e*l*q*u*x6 + 
     g*k*p*r*x1*y + h*i*p*t*x1*y - e*k*p*t*x1*y - g*i*p*u*x1*y - d*k*p*r*x2*y + c*k*q*r*x2*y - c*k*m*t*x2*y + 
     a*k*p*t*x2*y + d*i*p*u*x2*y - c*i*q*u*x2*y + d*h*p*r*x3*y - c*h*q*r*x3*y + c*h*m*t*x3*y - a*h*p*t*x3*y - 
     c*g*m*u*x3*y - d*e*p*u*x3*y + a*g*p*u*x3*y + c*e*q*u*x3*y - c*g*k*r*x4*y - c*h*i*t*x4*y + c*e*k*t*x4*y + 
     c*g*i*u*x4*y + c*g*k*m*x5*y - d*h*i*p*x5*y + d*e*k*p*x5*y - a*g*k*p*x5*y + c*h*i*q*x5*y - c*e*k*q*x5*y)/
   (d*h*l*p*s*v - c*h*l*q*s*v + c*h*l*n*t*v - b*h*l*p*t*v - c*g*l*n*u*v - d*f*l*p*u*v + b*g*l*p*u*v + 
     c*f*l*q*u*v - d*h*l*p*r*w + c*h*l*q*r*w - c*h*l*m*t*w + a*h*l*p*t*w + c*g*l*m*u*w + d*e*l*p*u*w - 
     a*g*l*p*u*w - c*e*l*q*u*w + c*g*l*n*r*x + d*f*l*p*r*x - b*g*l*p*r*x - c*f*l*q*r*x - c*g*l*m*s*x - 
     d*e*l*p*s*x + a*g*l*p*s*x + c*e*l*q*s*x + c*f*l*m*t*x - c*e*l*n*t*x + b*e*l*p*t*x - a*f*l*p*t*x - 
     c*g*k*n*r*y + d*h*j*p*r*y - d*f*k*p*r*y + b*g*k*p*r*y - c*h*j*q*r*y + c*f*k*q*r*y + c*g*k*m*s*y - 
     d*h*i*p*s*y + d*e*k*p*s*y - a*g*k*p*s*y + c*h*i*q*s*y - c*e*k*q*s*y + c*h*j*m*t*y - c*f*k*m*t*y - 
     c*h*i*n*t*y + c*e*k*n*t*y + b*h*i*p*t*y - a*h*j*p*t*y - b*e*k*p*t*y + a*f*k*p*t*y - c*g*j*m*u*y + 
     c*g*i*n*u*y + d*f*i*p*u*y - b*g*i*p*u*y - d*e*j*p*u*y + a*g*j*p*u*y - c*f*i*q*u*y + c*e*j*q*u*y);

  if (am_debugging_nans) printf("d_Ln is %g\n", d_Ln);
  if (am_debugging_nans) printf("d_Lf is %g\n", d_Lf);
}

double Dynein_bothbound::get_near_unbinding_rate() {
  //printf("bb.f.nby: %g, bb.r.nby: %g, both: %g", f.nby, r.nby, f.nby + r.nby);
  if (f.nby + r.nby >= BOTHBOUND_UNBINDING_FORCE) {
    return 1.0;
  } else return 0.0;
}

double Dynein_bothbound::get_far_unbinding_rate() {
  if (f.fby + r.fby >= BOTHBOUND_UNBINDING_FORCE) {
    return 1.0;
  } else return 0.0;
}

/*** Set positions, velocities and forces ***/
void Dynein_bothbound::set_nma(double d) {
    nma = d;
}

void Dynein_bothbound::set_fma(double d) {
    fma = d;
}

void Dynein_bothbound::set_L(double d) {
    L = d;
}

/*** Set velocities (for testing) ***/
void Dynein_bothbound::set_dLn(double d) {
    d_Ln = d;
}

void Dynein_bothbound::set_dLf(double d) {
    d_Lf = d;
}

/*** Angular Velocities ***/

double Dynein_bothbound::get_d_nba() {
  float d_An = -1 / sqrt(1 - (L*L+Ln*Ln-Lf*Lf)/(2*Ls*Ln)*(L*L+Ln*Ln-Lf*Lf)/(2*Ls*Ln))
    * ( (1/L - (L*L+Ln*Ln-Lf*Lf)/(2*L*Ln*Ln) )*d_Ln
	- (Lf/(L*Ln))*d_Lf );
  float d_Ans = -1 / sqrt(1 - (Ln*Ln+Ls*Ls-Lt*Lt)/(2*Ls*Ln)*(Ln*Ln+Ls*Ls-Lt*Lt)/(2*Ls*Ln))
    * ( 1/Ls - (Ln*Ln+Ls*Ls-Lt*Lt)/(2*Ls*Ln*Ls))*d_Ln;
  if (nma <= M_PI) return d_An + d_Ans;
  else return d_An - d_Ans;
}

double Dynein_bothbound::get_d_fba() {
  float d_Af = -1 / sqrt(1 - (Lf*Lf+L*L-Ln*Ln)/(2*L*Lf)*(Lf*Lf+L*L-Ln*Ln)/(2*L*Lf))
    * ( (1/L - (Lf*Lf+L*L-Ln*Ln)/(2*L*Lf*Lf) )*d_Lf
	- Ln / (L*Lf) * d_Ln);
  float d_Afs = -1 / sqrt(1 - (Ls*Ls+Lf*Lf-Lt*Lt)/(2*Ls*Lf)*(Ls*Ls+Lf*Lf-Lt*Lt)/(2*Ls*Lf))
    * ( 1/Ls - (Ls*Ls+Lf*Lf-Lt*Lt)/(2*Ls*Lf*Lf))*d_Lf;
  if (fma <= M_PI) return -(d_Af + d_Afs);
  else return -(d_Af - d_Afs);
}

double Dynein_bothbound::get_d_nma() {
  int pm = (nma > M_PI) ? -1 : 1; // sign of d_nma depends on value of nma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - Ln*Ln) / (2*lt*ls),2))
    * (Ln / (lt*ls)) * d_Ln;
}

double Dynein_bothbound::get_d_fma() {
  int pm = (fma > M_PI) ? -1 : 1; // sign of d_fma depends on value of fma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - Lf*Lf) / (2*lt*ls),2))
    * (Lf / (lt*ls)) * d_Lf;
}

/*** Get coordinates ***/
double Dynein_bothbound::get_nma() {
  return nma;
}

double Dynein_bothbound::get_fma() {
  return fma;
}

double Dynein_bothbound::get_nby() {
  return 0;
}

double Dynein_bothbound::get_nmy() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*(L*Ln)))
		+ pm*acos((pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) / (2*(Ln*ls))));
}

double Dynein_bothbound::get_ty() {
  return Ln*sqrt(1 - pow(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2), 2)
		 / (4*(pow(L, 2)*pow(Ln, 2))));
}

double Dynein_bothbound::get_fmy() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(Lf, 2) - pow(Ln, 2)) / (2*(L*Lf)))
		+ pm*acos((pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) / (2*(Lf*ls))));
}

double Dynein_bothbound::get_fby() {
  return 0;
}


/*** Get Cartesian Velocities ***/
double Dynein_bothbound::get_d_nmx() {
  return -ls * sin(nba) * get_d_nba();
}

double Dynein_bothbound::get_d_tx() {
  return get_d_nmx() + -lt * sin(get_nma()) * get_d_nma();
}

double Dynein_bothbound::get_d_fmx() {
  return -ls * sin(fba) * get_d_fba();
}

double Dynein_bothbound::get_d_nmy() {
  return ls * cos(nba) * get_d_nba();
}

double Dynein_bothbound::get_d_ty() {
  return get_d_nmy() + lt * cos(nma) * get_d_nma();
}

double Dynein_bothbound::get_d_fmy() {
  return ls * cos(fba) * get_d_fba();
}

/*** Get forces ***/
bothbound_forces Dynein_bothbound::get_internal() {
  return f;
}

bothbound_forces Dynein_bothbound::get_brownian() {
  return r;
}

/*** Get energies ***/

double Dynein_bothbound::get_PE() {
  return 0.5*cb*square(nba - eq.nba) + 0.5*cb*square(nma - eq.nma)
    + 0.5*cb*square(fma - eq.fma) + 0.5*cb*square(fba - eq.fba);
}

double Dynein_bothbound::get_KE() {
  return 0;
}

void Dynein_bothbound::log(double t, FILE* data_file) {
  fprintf(data_file, "%.2g\t%.2g\t%.2g\t%10g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	  "\t%.4f\t%.4f\t%.4f\t%.4f\t%d"
          "\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g"
          "\n",
          get_KE(), get_PE(), get_KE() + get_PE(), t,
          nbx, get_nby(), get_nmx(), get_nmy(), get_tx(), get_ty(), get_fmx(), get_fmy(), get_fbx(), get_fby(), BOTHBOUND,
          f.nbx, f.nby, f.nmx, f.nmy, f.tx, f.ty, f.fmx, f.fmy, f.fbx, f.fby);
  fprintf(data_file, "# theta_nm = %g\ttheta_fm = %g\tL = %g\n", nma, fma, L);
}
