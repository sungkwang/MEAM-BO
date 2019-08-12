#include "meambo.h"
#include "math_special.h"

using namespace LAMMPS_NS;

//-----------------------------------------------------------------------------
// This function finalize the intemediate density terms using the values
// received from ghost region by reverse_comm.
// The results will be copied/overwritten to the ghost region by forward_comm.
//
void
MEAMBO::meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
                      double* eatom,  int* type, int* fmap, int& errorflag)
{
  int i, elti;
  int m;
  double rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  double B, denom, rho_bkgd, rho, Z0i;

  for (i = 0; i < nlocal; i++) {
    elti = fmap[type[i]];
    if (elti >= 0) {
      rho1[i] = 0.0;
      rho2[i] = -1.0 / 3.0 * arho2b[i] * arho2b[i];
      rho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        rho1[i] = rho1[i] + arho1[i][m] * arho1[i][m];
        rho3[i] = rho3[i] - 3.0 / 5.0 * arho3b[i][m] * arho3b[i][m];
      }
      for (m = 0; m < 6; m++) {
        rho2[i] = rho2[i] + this->v2D[m] * arho2[i][m] * arho2[i][m];
      }
      for (m = 0; m < 10; m++) {
        rho3[i] = rho3[i] + this->v3D[m] * arho3[i][m] * arho3[i][m];
      }

      if (rho0[i] > 0.0) {
        if (this->ialloy == 1) {
          t_ave[i][0] = fdiv_zero(t_ave[i][0], tsq_ave[i][0]);
          t_ave[i][1] = fdiv_zero(t_ave[i][1], tsq_ave[i][1]);
          t_ave[i][2] = fdiv_zero(t_ave[i][2], tsq_ave[i][2]);
        } else if (this->ialloy == 2) {
          t_ave[i][0] = this->t1_meam[elti];
          t_ave[i][1] = this->t2_meam[elti];
          t_ave[i][2] = this->t3_meam[elti];
        } else {
          t_ave[i][0] = t_ave[i][0] / rho0[i];
          t_ave[i][1] = t_ave[i][1] / rho0[i];
          t_ave[i][2] = t_ave[i][2] / rho0[i];
        }
      }

      gamma[i] = t_ave[i][0] * rho1[i] + t_ave[i][1] * rho2[i] + t_ave[i][2] * rho3[i];

      if (rho0[i] > 0.0) {
        gamma[i] = gamma[i] / (rho0[i] * rho0[i]);
      }

      Z = this->Z_meam[elti];

      G = G_gam(gamma[i], this->ibar_meam[elti], errorflag);
      if (errorflag != 0)
        return;
      
      // get_shpfcn(this->lattce_meam[elti][elti], shp);
      get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);
      
      if (this->ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        if (this->mix_ref_t == 1) {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] + t_ave[i][2] * shp[2]) / (Z * Z);
        } else {
          gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
                (Z * Z);
        }
        Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
      }
      rho = rho0[i] * G;

      if (this->mix_ref_t == 1) {
        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] + t_ave[i][2] * shp[2]) / (Z * Z);
          Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
        }
        rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
      } else {
        if (this->bkgd_dyn == 1) {
          rho_bkgd = this->rho0_meam[elti] * Z;
        } else {
          rho_bkgd = this->rho_ref_meam[elti];
        }
      }
      rhob = rho / rho_bkgd;
      denom = 1.0 / rho_bkgd;

      G = dG_gam(gamma[i], this->ibar_meam[elti], dG);

      dgamma1[i] = (G - 2 * dG * gamma[i]) * denom;

      if (!iszero(rho0[i])) {
        dgamma2[i] = (dG / rho0[i]) * denom;
      } else {
        dgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (this->mix_ref_t == 1) {
        dgamma3[i] = rho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        dgamma3[i] = 0.0;
      }

      B = this->A_meam[elti] * this->Ec_meam[elti][elti];

      if (!iszero(rhob)) {
        if (this->emb_lin_neg == 1 && rhob <= 0) {
          frhop[i] = -B;
        } else {
          frhop[i] = B * (log(rhob) + 1.0);
        }
        if (eflag_either != 0) {
          if (eflag_global != 0) {
            if (this->emb_lin_neg == 1 && rhob <= 0) {
              *eng_vdwl = *eng_vdwl - B * rhob;
            } else {
              *eng_vdwl = *eng_vdwl + B * rhob * log(rhob);
            }
          }
          if (eflag_atom != 0) {
            if (this->emb_lin_neg == 1 && rhob <= 0) {
              eatom[i] = eatom[i] - B * rhob;
            } else {
              eatom[i] = eatom[i] + B * rhob * log(rhob);
            }
          }
        }
      } else {
        if (this->emb_lin_neg == 1) {
          frhop[i] = -B;
        } else {
          frhop[i] = B;
        }
      }
      if (isBC){
        //Z1BC needs to be copied to the ghost region to be used in nI calculaion
        for (m = 0; m < 3; m++)
          Z1BC[i] += arho1BC[i][m]*arho1BC[i][m];
        // precalculation of terms of eq. (27) in MEAM-BO paper
        Z0i = Z0BC[i] - this->z0s_meambo[BCty[DOUBLE]];
        get_delta( Z0i*Z0i, this->betaBC[DOUBLE][0], this->powerBC[DOUBLE][0], delta2_Z0[i], ddelta2_Z0[i]);
        get_delta( Z1BC[i], this->betaBC[DOUBLE][1], this->powerBC[DOUBLE][1], delta2_Z1[i], ddelta2_Z1[i]);
      }

    }
  }
}


//-----------------------------------------------------------------------------
// this function calculates the energy due to the bond order
// the half of the bond energy is assigned to the pair of atoms that forms a bond
//
void
MEAMBO::meambo_energy(int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
      double* eatom,  int* type, int* fmap, int* jnum_full, int** jlist_full, double** x)
{
  int m, jn, j, bo, itBC, BCtyp;
  int iB, i1, i2;
  int i1i2indx, jnindx;
  double Sij, Si1i2;
  double Z0i;
  double enbond, di1i2m, rIsq;
  double rI;
  double Ef, Ebond[2], fbond[2];
  double BOI, BOI_2, BOI_22, re2, eb2[3], deltaZ3[2], E_2;
  double re3, eb3[3], rIre, exp_rIre, E_3;
  double delZ0i1, delZ0i2;
  double rscrn_vdW, rscrn_diff, rscrn_diff_E_vdW;


  for (iB = 0; iB < nBC; iB++) {
    i1 = iBC[iB][0];
    i2 = iBC[iB][1];

    // Z3 density
    for (m = 0; m < 3; m++)
      Z3BC[iB] -=  0.6*arho3bBC[iB][m]*arho3bBC[iB][m]; //3/5
    for (m = 0; m < 10; m++)
      Z3BC[iB] += this->v3D[m]*arho3BC[iB][m]*arho3BC[iB][m];

    // nI: weighted local density for the bond
    for (jn = 0; jn < jnum_full[i1]; jn++) {
      j = jlist_full[i1][jn];
      if (fmap[type[j]] < 0) continue;
      jnindx = ij_index_full[i1][j];
      Sij = scrfcn_full[jnindx]*fcpair_full[jnindx]*fcpair_vdW_full[jnindx];
      if (iszero(Sij)) continue;
      nIBC[iB] += Sij*delta2_Z0[j]*delta2_Z1[j];
    }
    for (jn = 0; jn < jnum_full[i2]; jn++) {
      j = jlist_full[i2][jn];
      if (fmap[type[j]] < 0) continue;
      jnindx = ij_index_full[i2][j];
      Sij = scrfcn_full[jnindx]*fcpair_full[jnindx]*fcpair_vdW_full[jnindx];
      if (iszero(Sij)) continue;
      nIBC[iB] += Sij*delta2_Z0[j]*delta2_Z1[j];
    }

    if (nIBC[iB] < 2.0) nIBC[iB] = 2.0;
    enbond = nIBC[iB];

    i1i2indx = ij_index_full[i1][i2];
    Si1i2 = scrfcn_full[i1i2indx]*fcpair_full[i1i2indx];
    rscrn_vdW = fcpair_vdW_full[i1i2indx];
    rscrn_diff = 1.0 - rscrn_vdW;

    rIsq = 0.0;
    for (m = 0; m < 3; m++){
      di1i2m= x[i2][m] - x[i1][m];
      rIsq += di1i2m*di1i2m;
    }
    rI = sqrt(rIsq);

    //vdW contribution to energy
    itBC = this->ityBC_meambo[0]; // the first one, carbon
    if ( (this->vdW_form == 1) && (!iszero(rscrn_diff)) && (!iszero(this->evdW_96LJ_meambo[itBC]))) {
      double rc = this->svdW_96LJ_meambo[itBC];
      double rcrI = rc/rI;
      double rcrI6 = rcrI*rcrI*rcrI*rcrI*rcrI*rcrI;
      double rcrI9 = rcrI6*rcrI*rcrI*rcrI;

      double E_vdW = this->evdW_96LJ_meambo[itBC]*(2*rcrI9 - 3*rcrI6);
      rscrn_diff_E_vdW = rscrn_diff*E_vdW;
    } else {
      rscrn_diff_E_vdW = 0.0;
    }

    // triple bond contribution to energy
    bo = TRIPLE;
    BCtyp = BCty[bo];
    Z0i = Z0BC[i1] - this->z0s_meambo[BCtyp];
    delZ0i1 = get_deltaZ( Z0i*Z0i, this->betaBC[bo][0], this->powerBC[bo][0]);
    Z0i = Z0BC[i2] - this->z0s_meambo[BCtyp];
    delZ0i2 = get_deltaZ( Z0i*Z0i, this->betaBC[bo][0], this->powerBC[bo][0]);
    deltaZ3[bo] = get_deltaZ( Z3BC[iB], this->betaBC[bo][3], this->powerBC[bo][3]);
    re3 = this->re3_meambo[BCtyp];
    for (m = 0; m < 3; m++)
      eb3[m] = this->e3_meambo[BCtyp][m];
    rIre = rI/re3-1.0;
    exp_rIre = MathSpecial::fm_exp(-this->betaBC[bo][2]*rIre);
    E_3 = eb3[0]*(1.0 + rIre*eb3[1] + rIre*rIre*eb3[2])*exp_rIre;
    Ebond[bo] = E_3 + rscrn_diff_E_vdW;
    fbond[bo] = Si1i2*(delZ0i1*delZ0i2)*deltaZ3[bo];

    // double bond contribution to energy
    bo = DOUBLE;
    BCtyp = BCty[bo];
    deltaZ3[bo] = get_deltaZ( Z3BC[iB], this->betaBC[bo][3], this->powerBC[bo][3]);
    if (enbond > 2.0) {
      BOI = (enbond + 2.0)/enbond;
      BOI_2 = BOI - 2.0;
      BOI_22 = BOI_2*BOI_2;
      re2 = this->re2a_meambo[BCtyp][0]
        + this->re2a_meambo[BCtyp][1]*BOI_2
        + this->re2a_meambo[BCtyp][2]*BOI_22;
      for (m = 0; m < 3; m++)
        eb2[m] = this->e2a_0_meambo[BCtyp][m]
              + this->e2a_1_meambo[BCtyp][m]*BOI_2
              + this->e2a_2_meambo[BCtyp][m]*BOI_22;
    } else {
      re2 = this->re2a_meambo[BCtyp][0];
      for (m = 0; m < 3; m++)
        eb2[m] = this->e2a_0_meambo[BCtyp][m];
    }
    rIre = rI/re2-1.0;
    exp_rIre = MathSpecial::fm_exp(-this->betaBC[bo][2]*rIre);
    E_2 = eb2[0]*( 1.0 + rIre*eb2[1] + rIre*rIre*eb2[2] )*exp_rIre;
    Ebond[bo] = E_2 + rscrn_diff_E_vdW;
    fbond[bo] = Si1i2*(delta2_Z0[i1]*delta2_Z0[i2])*(delta2_Z1[i1]*delta2_Z1[i2])*deltaZ3[bo]
          + fbond[TRIPLE]*(1.0 - deltaZ3[TRIPLE]);

    if (eflag_either != 0) {
      Ef = Ebond[TRIPLE]*fbond[TRIPLE] + Ebond[DOUBLE]*fbond[DOUBLE];
      if (eflag_global != 0) *eng_vdwl += Ef;
      if (eflag_atom != 0) {
        //!   put the energy into the 2 atoms that forms the bond equally
        Ef = Ef*0.5;
        eatom[i1] += Ef;
        eatom[i2] += Ef;
      }
    }
  }
}


