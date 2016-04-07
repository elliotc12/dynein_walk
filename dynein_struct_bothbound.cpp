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

  Ln = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(nma));
  Lf = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(fma));

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

  Ln = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(nma));
  Lf = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(fma));

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

    double T, f1, f2, f1x, f1y, f2x, f2y;

    T = cb*(get_nba() - eq.nba);
    f2 = T/ls;
    f2x = f2 * sin(get_nba());
    f2y = f2 * -cos(get_nba());
    f.nmx += f2x;
    f.nmy += f2y;
    f.nbx += -f2x; // Equal and opposite forces!  :)
    f.nby += -f2y; // Equal and opposite forces!  :)

    T = cm*(nma - eq.nma);
    f1 = T/ls;
    f2 = T/lt;
    f1x = f1 * sin(get_nba());
    f1y = f1 * -cos(get_nba());
    f2x = f2 * sin(nma - (M_PI - get_nba()));
    f2y = f2 * -cos(nma - (M_PI - get_nba()));
    f.nbx += f1x;
    f.nby += f1y;
    f.tx  += f2x;
    f.ty  += f2y;
    f.nmx += -(f1x + f2x);
    f.nmy += -(f1y + f2y);

    T = ct*(fma + get_fba() - nma - get_nba() - eq.ta);
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(nma + get_nba() - M_PI);
    f1y = f1 * -cos(nma + get_nba() - M_PI);
    f2x = f2 * -sin(fma + get_fba() - M_PI);  // not sure if these angles are right?
    f2y = f2 * cos(fma + get_fba() - M_PI);
    f.nmx += f1x;
    f.nmy += f1y;
    f.fmx += f2x;
    f.fmy += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*(fma - eq.fma);
    f1 = T / lt;
    f2 = T / ls;
    f1x = f1 * sin(fma + get_fba() - M_PI);
    f1y = f1 * -cos(fma + get_fba() - M_PI);
    f2x = f2 * sin(get_fba());
    f2y = f2 * -cos(get_fba());
    f.tx  += f1x;
    f.ty  += f1y;
    f.fbx += f2x;
    f.fby += f2y;
    f.fmx += -(f1x + f2x);
    f.fmy += -(f1y + f2y);

    T = cb*(get_fba() - eq.fba);
    f1 = T / ls;
    f1x = f1 * sin(get_fba());
    f1y = f1 * -cos(get_fba());
    f.fmx += f1x;
    f.fmy += f1y;
    f.fbx += -f1x;
    f.fby += -f1y;

    if (get_nby() < 0) f.nby += MICROTUBULE_REPULSION_FORCE * fabs(get_nby());
    if (get_nmy() < 0) f.nmy += MICROTUBULE_REPULSION_FORCE * fabs(get_nmy());
    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
    if (get_fmy() < 0) f.fmy += MICROTUBULE_REPULSION_FORCE * fabs(get_fmy());
    if (get_fby() < 0) f.fby += MICROTUBULE_REPULSION_FORCE * fabs(get_fby());
  }
}

void Dynein_bothbound::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  double cosAn = (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln);
  double sinAn = sqrt(1 - cosAn*cosAn);
  double cosAns = (Ls*Ls + Ln*Ln - Lt*Lt) / (2*Ls*Ln);
  double sinAns = (nma < M_PI) ? sqrt(1-cosAns*cosAns) : -sqrt(1-cosAns*cosAns);

  double cosAf = (L*L + Lf*Lf - Ln*Ln) / (2*L*Lf);
  double sinAf = sqrt(1 - cosAf*cosAf);
  double cosAfs = (Ls*Ls + Lf*Lf - Lt*Lt) / (2*Ls*Lf);
  double sinAfs = (fma < M_PI) ? sqrt(1-cosAfs*cosAfs) : -sqrt(1-cosAfs*cosAfs);

  double dcosAn_dLn = (1/L) - (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln*Ln);
  double dcosAn_dLf = -(Lf) / (L*Ln);
  double dsinAn_dLn = -cosAn / sqrt(1 - cosAn*cosAn) * (1/L - (L*L + Ln*Ln - Lf*Lf) / (2*L*Ln*Ln));
  double dsinAn_dLf = cosAn / sqrt(1 - cosAn*cosAn) * Lf / (L*Ln);
  double dcosAns_dLn = 1/Ls - (Ls*Ls + Ln*Ln - Lt*Lt) / (2*Ls*Ln*Ln);
  double dcosAns_dLf = 0;
  double dsinAns_dLn = (nma < M_PI) ?
     -cosAns/sqrt(1-cosAns*cosAns) * (1/Ls - (Ls*Ls+Ln*Ln-Lt*Lt)/(2*L*Ln*Ln))
    : cosAns/sqrt(1-cosAns*cosAns) * (1/Ls - (Ls*Ls+Ln*Ln-Lt*Lt)/(2*L*Ln*Ln));
  double dsinAns_dLf = 0;

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

  // David -- emacs comand to take all below and make it 100-char lines?

  d_Ln = (dYfm_dLf*f.fmx - dXfm_dLf*f.fmy + dYfm_dLf*gm*r.fmx - dXfm_dLf*gm*r.fmy)/((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm) - ((dYfm_dLf*get_fmx() - dXfm_dLf*get_fmy() - dYfm_dLf*get_tx() + dXfm_dLf*get_ty())*
      (-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dXt_dLn*f.fmx)/gm) + (dXfm_dLn*f.tx)/gt - dXt_dLn*r.fmx + dXfm_dLn*r.tx) - 
               (-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gt) + 

        (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + 
             (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty)))/gt))/
    ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
               (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
        (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
             (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)) - 
   (((-(dYfm_dLf*get_fbx()) + dXfm_dLf*get_fby() + dYfm_dLf*get_fmx() - dXfm_dLf*get_fmy())/((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm) - 
        ((dYfm_dLf*get_fmx() - dXfm_dLf*get_fmy() - dYfm_dLf*get_tx() + dXfm_dLf*get_ty())*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*
                  (-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt) + 
             (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt))/
         ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                    (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
             (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                  (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)))*
      ((-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                  (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)*
         (-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                     ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXfm_dLn*f.nmx)/gm - (dXnm_dLn*f.fmx)/gm + dXfm_dLn*r.nmx - dXnm_dLn*r.fmx) - (-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gm - 
                  (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYnm_dLn*f.fmx)/gm) + (dXfm_dLn*f.nmy)/gm - dYnm_dLn*r.fmx + dXfm_dLn*r.nmy) - 
                       (-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gm))/gt) - 
           (-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
            (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty))) - 
        (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
              (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt))) - 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                   (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm - 
                (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm))
             /gt)*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dXt_dLn*f.fmx)/gm) + (dXfm_dLn*f.tx)/gt - dXt_dLn*r.fmx + dXfm_dLn*r.tx) - 
                  (-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty)))/gt)))/
    (-((-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt)*
         (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
              (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt))) - 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                   (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm - 
                (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm))
             /gt)) + (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + 
              (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
            (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm)) - 
         (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                 (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gm - 
              (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gm))/
          gt)*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
         (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
              (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)); // from Mathematica

  d_Lf = (-(dYfm_dLn*f.fmx) + dXfm_dLn*f.fmy - dYfm_dLn*gm*r.fmx + dXfm_dLn*gm*r.fmy)/((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm) - 
   ((-(dYfm_dLn*get_fmx()) + dXfm_dLn*get_fmy() + dYfm_dLn*get_tx() - dXfm_dLn*get_ty())*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*
             ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dXt_dLn*f.fmx)/gm) + (dXfm_dLn*f.tx)/gt - dXt_dLn*r.fmx + dXfm_dLn*r.tx) - (-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gt) + 
        (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + 
             (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty)))/gt))/
    ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
               (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
        (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
             (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)) - 
   (((dYfm_dLn*get_fbx() - dXfm_dLn*get_fby() - dYfm_dLn*get_fmx() + dXfm_dLn*get_fmy())/((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm) - 
        ((-(dYfm_dLn*get_fmx()) + dXfm_dLn*get_fmy() + dYfm_dLn*get_tx() - dXfm_dLn*get_ty())*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*
                  (-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt) + 
             (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt))/
         ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*gm*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                    (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
             (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                  (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)))*
      ((-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                  (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)*
         (-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                     ((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXfm_dLn*f.nmx)/gm - (dXnm_dLn*f.fmx)/gm + dXfm_dLn*r.nmx - dXnm_dLn*r.fmx) - (-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gm - 
                  (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYnm_dLn*f.fmx)/gm) + (dXfm_dLn*f.nmy)/gm - dYnm_dLn*r.fmx + dXfm_dLn*r.nmy) - 
                       (-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gm))/gt) - 
           (-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
            (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty))) - 
        (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
              (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt))) - 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                   (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm - 
                (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm))
             /gt)*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dXt_dLn*f.fmx)/gm) + (dXfm_dLn*f.tx)/gt - dXt_dLn*r.fmx + dXfm_dLn*r.tx) - 
                  (-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*(-((dYfm_dLn*f.fmx)/gm) + (dXfm_dLn*f.fmy)/gm - dYfm_dLn*r.fmx + dXfm_dLn*r.fmy)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-((dYt_dLn*f.fmx)/gm) + (dXfm_dLn*f.ty)/gt - dYt_dLn*r.fmx + dXfm_dLn*r.ty)))/gt)))/
    (-((-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt) + 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gt)*
         (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
              (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt))) - 
           (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                   (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm - 
                (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(-get_fmx() + get_tx()))/gm))/gm))
             /gt)) + (-((-((pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nby() - get_nmy())*(-get_nmx() + get_tx()))/pow(gm,2)) + 
              (pow(dXfm_dLn,2)*pow(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn,2)*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/pow(gm,2))*
            (-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYt_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm)) - 
         (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_nmy())*
                 (-((-(dXnm_dLn*dXfm_dLf) + dXnm_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dXnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gm - 
              (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nbx() - get_nmx())*(-((-(dYnm_dLn*dXfm_dLf) + dYnm_dLf*dXfm_dLn)*((dYfm_dLn*(get_fbx() - get_fmx()))/gm - (dXfm_dLn*(get_fby() - get_fmy()))/gm)) + (dYnm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_fbx() - get_fmx()))/gm))/gm))/
          gt)*(-((dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nby() - get_ty())*(-((-(dXt_dLn*dXfm_dLf) + dXt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
                (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dXt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmx() - get_tx()))/gt)))/gt) + 
         (dXfm_dLn*(dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*(get_nmx() - get_tx())*(-((-(dYt_dLn*dXfm_dLf) + dYt_dLf*dXfm_dLn)*((dYfm_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(-get_fmy() + get_ty()))/gm)) + 
              (dXfm_dLn*dYfm_dLf - dXfm_dLf*dYfm_dLn)*((dYt_dLn*(-get_fmx() + get_tx()))/gm - (dXfm_dLn*(get_fmy() - get_ty()))/gt)))/gt)); // from Mathematica
}

double Dynein_bothbound::get_near_unbinding_rate() {
  //printf("bb.f.nby: %g, bb.r.nby: %g, both: %g\n", f.nby, r.nby, f.nby + r.nby);
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

/*** Angular Velocities ***/

double Dynein_bothbound::get_d_nba() { // TODO: make this right
    return 1.0;
}

double Dynein_bothbound::get_d_fba() {
    return 1.0;
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

double Dynein_bothbound::get_nba() {
  double pm_n = (nma > M_PI) ? -1.0 : 1.0;
  return acos((L*L + Ln*Ln - Lf*Lf)/(2*L*Ln)) + pm_n * acos((Ls*Ls + Ln*Ln - Lt*Lt)/(2*Ls*Ln));
}

double Dynein_bothbound::get_fba() {
  double pm_f = (fma > M_PI) ? -1.0 : 1.0; // right?
  return acos((L*L + Lf*Lf - Ln*Ln)/(2*L*Lf)) + pm_f * acos((Ls*Ls + Lf*Lf - Lt*Lt)/(2*Ls*Lf));
}


double Dynein_bothbound::get_nbx() {
  return nbx;
}

double Dynein_bothbound::get_nmx() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*(L*Ln)))
         + pm*acos((pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) / (2*(Ln*ls))));
}

double Dynein_bothbound::get_tx() {
  return (pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*L);
}

double Dynein_bothbound::get_fmx() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(Lf, 2) - pow(Ln, 2)) / (2*(L*Lf)))
		+ pm*acos((pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) / (2*(Lf*ls))));;
}

double Dynein_bothbound::get_fbx() {
  return L;
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
  return -ls * sin(get_nba()) * get_d_nba();
}

double Dynein_bothbound::get_d_tx() {
  return get_d_nmx() + -lt * sin(get_nma()) * get_d_nma();
}

double Dynein_bothbound::get_d_fmx() {
  return -ls * sin(get_fba()) * get_d_fba();
}

double Dynein_bothbound::get_d_nmy() {
  return ls * cos(get_nba()) * get_d_nba();
}

double Dynein_bothbound::get_d_ty() {
  return get_d_nmy() + lt * cos(nma) * get_d_nma();
}

double Dynein_bothbound::get_d_fmy() {
  return ls * cos(get_fba()) * get_d_fba();
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
  return 0.5*cb*square(get_nba() - eq.nba) + 0.5*cb*square(nma - eq.nma)
    + 0.5*cb*square(fma - eq.fma) + 0.5*cb*square(get_fba() - eq.fba);
}

double Dynein_bothbound::get_KE() {
  return 0;
}

void Dynein_bothbound::log(int step, FILE* data_file) {
  fprintf(data_file, "%.2g\t%.2g\t%.2g\t%10d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	  "\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
          get_KE(), get_PE(), get_KE() + get_PE(), step, get_nbx(), get_nby(), get_nmx(), get_nmy(),
          get_tx(), get_ty(), get_fmx(), get_fmy(), get_fbx(), get_fby(), BOTHBOUND);
}
