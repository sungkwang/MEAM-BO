#include "meambo.h"
#include "math_special.h"

using namespace LAMMPS_NS;

//-----------------------------------------------------------------------------
// This function allocate/initialize the variables regarding the number of atoms.
// The arrays and matrices are reallocated when the required data size increases.
void
MEAMBO::meam_dens_setup(int atom_nmax, int nall, int n_neigh, int n_neigh_full)
{
  int i, j;

  // whenever new variables are used here, they needs to be freed in meambo_impl.cpp
  // grow local arrays if necessary
  if (atom_nmax > this->nmax) {
    memory->destroy(rho0);
    memory->destroy(rho1);
    memory->destroy(rho2);
    memory->destroy(rho3);
    memory->destroy(frhop);
    memory->destroy(gamma);
    memory->destroy(dgamma1);
    memory->destroy(dgamma2);
    memory->destroy(dgamma3);
    memory->destroy(arho2b);
    memory->destroy(arho1);
    memory->destroy(arho2);
    memory->destroy(arho3);
    memory->destroy(arho3b);
    memory->destroy(t_ave);
    memory->destroy(tsq_ave);

    this->nmax = atom_nmax;

    memory->create(rho0, this->nmax, "pair:rho0");
    memory->create(rho1, this->nmax, "pair:rho1");
    memory->create(rho2, this->nmax, "pair:rho2");
    memory->create(rho3, this->nmax, "pair:rho3");
    memory->create(frhop, this->nmax, "pair:frhop");
    memory->create(gamma, this->nmax, "pair:gamma");
    memory->create(dgamma1, this->nmax, "pair:dgamma1");
    memory->create(dgamma2, this->nmax, "pair:dgamma2");
    memory->create(dgamma3, this->nmax, "pair:dgamma3");
    memory->create(arho2b, this->nmax, "pair:arho2b");
    memory->create(arho1, this->nmax, 3, "pair:arho1");
    memory->create(arho2, this->nmax, 6, "pair:arho2");
    memory->create(arho3, this->nmax, 10, "pair:arho3");
    memory->create(arho3b, this->nmax, 3, "pair:arho3b");
    memory->create(t_ave, this->nmax, 3, "pair:t_ave");
    memory->create(tsq_ave, this->nmax, 3, "pair:tsq_ave");


    //meambo
    if (isBC) {

      memory->destroy(iBC);
      memory->destroy(Z0BC);
      memory->destroy(Z1BC);
      memory->destroy(arho1BC);
      memory->destroy(dZ0BC);
      memory->destroy(dZ1BC);
      memory->destroy(delta2_Z0);
      memory->destroy(ddelta2_Z0);
      memory->destroy(delta2_Z1);
      memory->destroy(ddelta2_Z1);
      memory->destroy(dnI);
      memory->destroy(slocal);
      memory->destroy(neighbor_check);
      memory->destroy(ij_index_full);

      memory->create(iBC, this->nmax, 2, "pair:iBC");
      memory->create(Z0BC, this->nmax, "pair:Z0BC");
      memory->create(Z1BC, this->nmax, "pair:Z1BC");
      memory->create(arho1BC, this->nmax, 3, "pair:arho1BC");
      memory->create(dZ0BC, this->nmax, 3,  "pair:dZ0BC");
      memory->create(dZ1BC, this->nmax, 3,  "pair:dZ1BC");
      memory->create(delta2_Z0, this->nmax, "pair:delta2_Z0");
      memory->create(ddelta2_Z0, this->nmax, "pair:ddelta2_Z0");
      memory->create(delta2_Z1, this->nmax, "pair:delta2_Z1");
      memory->create(ddelta2_Z1, this->nmax, "pair:ddelta2_Z1");
      memory->create(dnI, this->nmax, 3, "pair:dnI");
      memory->create(slocal, this->nmax, 3, 3, "pair:slocal");
      memory->create(neighbor_check, this->nmax, "pair:neighbor_check");
      memory->create(ij_index_full, this->nmax, this->nmax, "pair:ij_index_full");

    }
  }

  if (n_neigh > this->maxneigh) {
    this->maxneigh = n_neigh;

    // if bond order calculation is involved, then use screening terms based on
    // full neighbors not half neighbors
    if (!isBC) {
      memory->destroy(scrfcn);
      memory->destroy(dscrfcn);
      memory->destroy(fcpair);

      memory->create(scrfcn, this->maxneigh, "pair:scrfcn");
      memory->create(dscrfcn, this->maxneigh, "pair:dscrfcn");
      memory->create(fcpair, this->maxneigh, "pair:fcpair");
    }
  }


  if (isBC) {
    if (n_neigh_full > maxneigh_full) {
      maxneigh_full = n_neigh_full;

      memory->destroy(scrfcn_full);
      memory->destroy(dscrfcn_full);
      memory->destroy(fcpair_full);
      memory->destroy(fcpair_vdW_full);
      memory->destroy(dfcpair_vdW_full);
      memory->destroy(dscrfcn_vdW_full);
      memory->destroy(dZ0jBC);
      memory->destroy(dZ1jBC);

      memory->create(scrfcn_full, maxneigh_full, "pair:scrfcn_full");
      memory->create(dscrfcn_full, maxneigh_full, "pair:dscrfcn_full");
      memory->create(fcpair_full, maxneigh_full, "pair:fcpair_full");
      memory->create(fcpair_vdW_full, maxneigh_full, "pair:fcpair_vdW_full");
      memory->create(dfcpair_vdW_full, maxneigh_full, "pair:dfcpair_vdW_full");
      memory->create(dscrfcn_vdW_full, maxneigh_full, "pair:dscrfcn_vdW_full");
      memory->create(dZ0jBC, maxneigh_full,3, "pair:dZ0jBC");
      memory->create(dZ1jBC, maxneigh_full,3, "pair:dZ1jBC");
    }
  }


  // zero out local arrays
  for (i = 0; i < nall; i++) {
    rho0[i] = 0.0;
    arho2b[i] = 0.0;
    arho1[i][0] = arho1[i][1] = arho1[i][2] = 0.0;
    for (j = 0; j < 6; j++)
      arho2[i][j] = 0.0;
    for (j = 0; j < 10; j++)
      arho3[i][j] = 0.0;
    arho3b[i][0] = arho3b[i][1] = arho3b[i][2] = 0.0;
    t_ave[i][0] = t_ave[i][1] = t_ave[i][2] = 0.0;
    tsq_ave[i][0] = tsq_ave[i][1] = tsq_ave[i][2] = 0.0;
  }

  if (isBC) {
    #pragma omp simd // for vectorization with no dependency
    for (i = 0; i < nall; i++) {
      Z0BC[i] = 0.0;
      Z1BC[i] = 0.0;
      arho1BC[i][0] = arho1BC[i][1] = arho1BC[i][2] = 0.0;

      dZ0BC[i][0] = dZ0BC[i][1] = dZ0BC[i][2] =  0.0;
      dZ1BC[i][0] = dZ1BC[i][1] = dZ1BC[i][2] =  0.0;

      slocal[i][0][0] = slocal[i][0][1] = slocal[i][0][2] = 0.0;
      slocal[i][1][0] = slocal[i][1][1] = slocal[i][1][2] = 0.0;
      slocal[i][2][0] = slocal[i][2][1] = slocal[i][2][2] = 0.0;

      neighbor_check[i] = false;
    }

    #pragma omp simd // for vectorization with no dependency
    for (j = 0; j < maxneigh_full; j++){
      dZ0jBC[j][0] = dZ0jBC[j][1] = dZ0jBC[j][2] = 0.0;
      dZ1jBC[j][0] = dZ1jBC[j][1] = dZ1jBC[j][2] = 0.0;
    }
  }
}

void
MEAMBO::meam_dens_init(int i,  int* type, int* fmap, double** x,
                     int numneigh, int* firstneigh,
                     int numneigh_full, int* firstneigh_full, int fnoffset)
{
  //     Compute screening function and derivatives
  // if MEAM-BO calculation, reuse the getscreen_full value instead
  if (!isBC) getscreen(i, &scrfcn[fnoffset], &dscrfcn[fnoffset], &fcpair[fnoffset], x, numneigh,
        firstneigh, numneigh_full, firstneigh_full, type, fmap);

  //     Calculate intermediate density terms to be communicated
  calc_rho1(i, type, fmap, x, numneigh, firstneigh, fnoffset);
}


//-----------------------------------------------------------------------------
// This function calculates the screening between an atom and its half neighbors.
//
void
MEAMBO::getscreen(int i, double* scrfcn, double* dscrfcn, double* fcpair, double** x, int numneigh,
                int* firstneigh, int numneigh_full, int* firstneigh_full, int* type, int* fmap)
{
  int jn, j, kn, k;
  int elti, eltj, eltk;
  double xitmp, yitmp, zitmp, delxij, delyij, delzij, rij2, rij;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2 /*,rik*/;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2 /*,rjk*/;
  double xik, xjk, sij, fcij, sfcij, dfcij, sikj, dfikj, cikj;
  double Cmin, Cmax, delc, a, coef1, coef2;
  double dCikj;
  double rnorm, fc, dfc, drinv;

  drinv = 1.0 / this->delr_meam;
  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];

    eltj = fmap[type[j]];
    if (eltj < 0) continue;

    //     First compute screening function itself, sij
    xjtmp = x[j][0];
    yjtmp = x[j][1];
    zjtmp = x[j][2];
    delxij = xjtmp - xitmp;
    delyij = yjtmp - yitmp;
    delzij = zjtmp - zitmp;
    rij2 = delxij * delxij + delyij * delyij + delzij * delzij;
    rij = sqrt(rij2);

    const double rbound = this->ebound_meam[elti][eltj] * rij2;
    if (rij > this->rc_meam) {
      fcij = 0.0;
      dfcij = 0.0;
      sij = 0.0;
    } else {
      rnorm = (this->rc_meam - rij) * drinv;
      sij = 1.0;

      //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
      for (kn = 0; kn < numneigh_full; kn++) {
        k = firstneigh_full[kn];
        eltk = fmap[type[k]];
        if (eltk < 0) continue;
        if (k == j) continue;

        delxjk = x[k][0] - xjtmp;
        delyjk = x[k][1] - yjtmp;
        delzjk = x[k][2] - zjtmp;
        rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
        if (rjk2 > rbound) continue;

        delxik = x[k][0] - xitmp;
        delyik = x[k][1] - yitmp;
        delzik = x[k][2] - zitmp;
        rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
        if (rik2 > rbound) continue;

        xik = rik2 / rij2;
        xjk = rjk2 / rij2;
        a = 1 - (xik - xjk) * (xik - xjk);
        //     if a < 0, then ellipse equation doesn't describe this case and
        //     atom k can't possibly screen i-j
        if (a <= 0.0) continue;

        cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
        Cmax = this->Cmax_meam[elti][eltj][eltk];
        Cmin = this->Cmin_meam[elti][eltj][eltk];
        if (cikj >= Cmax) continue;
        //     note that cikj may be slightly negative (within numerical
        //     tolerance) if atoms are colinear, so don't reject that case here
        //     (other negative cikj cases were handled by the test on "a" above)
        else if (cikj <= Cmin) {
          sij = 0.0;
          break;
        } else {
          delc = Cmax - Cmin;
          cikj = (cikj - Cmin) / delc;
          sikj = fcut(cikj);
        }
        sij *= sikj;
      }

      fc = dfcut(rnorm, dfc);
      fcij = fc;
      dfcij = dfc * drinv;
    }

    //     Now compute derivatives
    dscrfcn[jn] = 0.0;
    sfcij = sij * fcij;
    if (iszero(sfcij) || iszero(sfcij - 1.0))
      goto LABEL_100;

    for (kn = 0; kn < numneigh_full; kn++) {
      k = firstneigh_full[kn];
      if (k == j) continue;
      eltk = fmap[type[k]];
      if (eltk < 0) continue;

      xktmp = x[k][0];
      yktmp = x[k][1];
      zktmp = x[k][2];
      delxjk = xktmp - xjtmp;
      delyjk = yktmp - yjtmp;
      delzjk = zktmp - zjtmp;
      rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
      if (rjk2 > rbound) continue;

      delxik = xktmp - xitmp;
      delyik = yktmp - yitmp;
      delzik = zktmp - zitmp;
      rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
      if (rik2 > rbound) continue;

      xik = rik2 / rij2;
      xjk = rjk2 / rij2;
      a = 1 - (xik - xjk) * (xik - xjk);
      //     if a < 0, then ellipse equation doesn't describe this case and
      //     atom k can't possibly screen i-j
      if (a <= 0.0) continue;

      cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
      Cmax = this->Cmax_meam[elti][eltj][eltk];
      Cmin = this->Cmin_meam[elti][eltj][eltk];
      if (cikj >= Cmax) {
        continue;
        //     Note that cikj may be slightly negative (within numerical
        //     tolerance) if atoms are colinear, so don't reject that case
        //     here
        //     (other negative cikj cases were handled by the test on "a"
        //     above)
        //     Note that we never have 0<cikj<Cmin here, else sij=0
        //     (rejected above)
      } else {
        delc = Cmax - Cmin;
        cikj = (cikj - Cmin) / delc;
        sikj = dfcut(cikj, dfikj);
        coef1 = dfikj / (delc * sikj);
        dCikj = dCfunc(rij2, rik2, rjk2);
        dscrfcn[jn] = dscrfcn[jn] + coef1 * dCikj;
      }
    }
    coef1 = sfcij;
    coef2 = sij * dfcij / rij;
    dscrfcn[jn] = dscrfcn[jn] * coef1 - coef2;

    LABEL_100:
    scrfcn[jn] = sij;
    fcpair[jn] = fcij;
  }
}

//-----------------------------------------------------------------------------
// This function calculates the intermediate density terms to be communicated
// with the surrounding ghost regions by reverse_comm (summation, from ghost to local)
// and forward_comm (copy, local to ghost).
//
void
MEAMBO::calc_rho1(int i,  int* type, int* fmap, double** x, int numneigh,
                int* firstneigh, int fnoffset)
{
  int jn, j, m, n, p, elti, eltj;
  int nv2, nv3;
  double xtmp, ytmp, ztmp, delij[3], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  // double G,Gbar,gam,shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;

  double scr;
  int jnindx;  


  elti = fmap[type[i]];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    if (isBC){
      jnindx = ij_index_full[i][j];
      scr = scrfcn_full[jnindx];
    } else {
      jnindx = fnoffset+jn;
      scr = scrfcn[jnindx];
    }

    if (!iszero(scr)) {
      delij[0] = x[j][0] - xtmp;
      delij[1] = x[j][1] - ytmp;
      delij[2] = x[j][2] - ztmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->rcutsq) {
        eltj = fmap[type[j]];
        rij = sqrt(rij2);
        ai = rij / this->re_meam[elti][elti] - 1.0;
        aj = rij / this->re_meam[eltj][eltj] - 1.0;
        ro0i = this->rho0_meam[elti];
        ro0j = this->rho0_meam[eltj];

        if (isBC)
          sij = scrfcn_full[jnindx]*fcpair_full[jnindx];
        else
          sij = scrfcn[jnindx]*fcpair[jnindx];

        rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj) * sij;
        rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj) * sij;
        rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj) * sij;
        rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj) * sij;
        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai) * sij;
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai) * sij;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai) * sij;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai) * sij;
        if (this->ialloy == 1) {
          rhoa1j = rhoa1j * this->t1_meam[eltj];
          rhoa2j = rhoa2j * this->t2_meam[eltj];
          rhoa3j = rhoa3j * this->t3_meam[eltj];
          rhoa1i = rhoa1i * this->t1_meam[elti];
          rhoa2i = rhoa2i * this->t2_meam[elti];
          rhoa3i = rhoa3i * this->t3_meam[elti];
        }
        rho0[i] = rho0[i] + rhoa0j;
        rho0[j] = rho0[j] + rhoa0i;
    
        // For ialloy = 2, use single-element value (not average)
        if (this->ialloy != 2) {
          t_ave[i][0] = t_ave[i][0] + this->t1_meam[eltj] * rhoa0j;
          t_ave[i][1] = t_ave[i][1] + this->t2_meam[eltj] * rhoa0j;
          t_ave[i][2] = t_ave[i][2] + this->t3_meam[eltj] * rhoa0j;
          t_ave[j][0] = t_ave[j][0] + this->t1_meam[elti] * rhoa0i;
          t_ave[j][1] = t_ave[j][1] + this->t2_meam[elti] * rhoa0i;
          t_ave[j][2] = t_ave[j][2] + this->t3_meam[elti] * rhoa0i;
        }
        if (this->ialloy == 1) {
          tsq_ave[i][0] = tsq_ave[i][0] + this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          tsq_ave[i][1] = tsq_ave[i][1] + this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          tsq_ave[i][2] = tsq_ave[i][2] + this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
          tsq_ave[j][0] = tsq_ave[j][0] + this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          tsq_ave[j][1] = tsq_ave[j][1] + this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          tsq_ave[j][2] = tsq_ave[j][2] + this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
        arho2b[i] = arho2b[i] + rhoa2j;
        arho2b[j] = arho2b[j] + rhoa2i;

        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
        nv2 = 0;
        nv3 = 0;
        for (m = 0; m < 3; m++) {
          arho1[i][m] = arho1[i][m] + A1j * delij[m];
          arho1[j][m] = arho1[j][m] - A1i * delij[m];
          arho3b[i][m] = arho3b[i][m] + rhoa3j * delij[m] / rij;
          arho3b[j][m] = arho3b[j][m] - rhoa3i * delij[m] / rij;
          for (n = m; n < 3; n++) {
            arho2[i][nv2] = arho2[i][nv2] + A2j * delij[m] * delij[n];
            arho2[j][nv2] = arho2[j][nv2] + A2i * delij[m] * delij[n];
            nv2 = nv2 + 1;
            for (p = n; p < 3; p++) {
              arho3[i][nv3] = arho3[i][nv3] + A3j * delij[m] * delij[n] * delij[p];
              arho3[j][nv3] = arho3[j][nv3] - A3i * delij[m] * delij[n] * delij[p];
              nv3 = nv3 + 1;
            }
          }
        }
        if (isBC){ // meambo
          double sij_vdW = scrfcn_full[jnindx]*fcpair_full[jnindx]*fcpair_vdW_full[jnindx];
          if (!iszero(sij_vdW)) {
            Z0BC[i] += sij_vdW;
            Z0BC[j] += sij_vdW;

            double A1BC = sij_vdW/rij;
            for (m = 0; m < 3; m++) {
              arho1BC[i][m] += A1BC*delij[m];
              arho1BC[j][m] -= A1BC*delij[m];
            }
          }
        }
      }
    }
  }
}


// This function allocate/initialize the variables regarding the number of bonds.
// The arrays and matrices are reallocated when the required data size increases.
void
MEAMBO::meambo_dens_setup(int bond_nmax, int n_neigh, int n_single_neigh_max)
{
  int i, j;

  bool change_any = false;
  if (bond_nmax > this->nmaxBC) {
    change_any = true;
    memory->destroy(jnumBC);
    memory->destroy(Z3BC);
    memory->destroy(nIBC);
    memory->destroy(dZ3BC);
    memory->destroy(arho3BC);
    memory->destroy(arho3bBC);

    this->nmaxBC = bond_nmax;

    // to store the neighbors of i1 and i2,
    memory->create(jnumBC, nmaxBC, "pair:jnumBC");
    memory->create(Z3BC, nmaxBC, "pair:Z3BC");
    memory->create(nIBC, nmaxBC, "pair:nIBC");
    memory->create(dZ3BC, nmaxBC, 3, "pair:dZ3BC");
    memory->create(arho3BC, nmaxBC, 10, "pair:arho3BC");
    memory->create(arho3bBC, nmaxBC, 3, "pair:arho3bBC");

  }

  if (n_single_neigh_max > this->maxneighBC_single){
    this->maxneighBC_single = n_single_neigh_max;
    change_any = true;
  }

  if (change_any){
    // to store the neighbors of i1 and i2
    memory->destroy(jlistBC);
    memory->create(jlistBC, nmaxBC, maxneighBC_single, "pair:jlistBC");
  }

  if (n_neigh > this->maxneighBC) {
    memory->destroy(scrfcnBC);
    memory->destroy(fcpairBC);
    memory->destroy(fcpair_vdWBC);
    memory->destroy(dscrfcnBC);
    memory->destroy(dZ3jBC);

    this->maxneighBC = n_neigh;

    memory->create(scrfcnBC, maxneighBC, "pair:scrfcnBC");
    memory->create(fcpairBC, maxneighBC, "pair:fcpairBC");
    memory->create(fcpair_vdWBC, maxneighBC, "pair:fcpair_vdWBC");
    memory->create(dscrfcnBC, maxneighBC, "pair:dscrfcnBC");
    memory->create(dZ3jBC, maxneighBC,3, "pair:dZ3jBC");
  }



  for (i = 0; i < this->nBC; i++) {
    Z3BC[i] = 0.0;
    nIBC[i] = 0.0;
    for (j = 0; j < 10; j++) arho3BC[i][j] = 0.0;
    arho3bBC[i][0] = arho3bBC[i][1] = arho3bBC[i][2] = 0.0;
    dZ3BC[i][0] =dZ3BC[i][1] =dZ3BC[i][2] =  0.0;
  }

  for (j = 0; j < maxneighBC; j++){
    dZ3jBC[j][0] = dZ3jBC[j][1] = dZ3jBC[j][2] = 0.0;
  }

}

//-----------------------------------------------------------------------------
// This function calculates the screening between an atom and its full neighbors.
// The results may be used instead of the results of getscreen function.
//
void
MEAMBO::meambo_get_full_screen( int* type, int* fmap, double** x, int inum, int* ilist,
      int* jnum_full, int** jlist_full)
{
  int ii,i, fnoffset;
  int jn, j, kn, k;
  int elti, eltj, eltk;
  double xitmp, yitmp, zitmp, dij[3], rij2, rij;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2;
  double xik, xjk, sij, fcij, sfcij, dfcij, sikj, dfikj, cikj;
  double Cmin, Cmax, delc, rbound, a, coef1, coef2;
  double dCikj;
  double rnorm, dfc, drinv;

  // meambo
  // two additional arrays are created, dscrfcn_vdW and fcpair_vdW
  // only if there is any possible unsaturate bond
  // therefore, the arrays are not used as arguments
  // but fnoffset is explicitly used for those conditional arrrays
  double rnorm_vdW, rcutBC, fcij_vdW, dfc_vdW, dfcij_vdW, drinv_vdW, coef3, sfcij_vdW;

  // double frcut = 1 - this->delr_meam/this->rc_meam;
  double delr_rc = this->delr_meam/this->rc_meam;
  drinv = 1.0 / this->delr_meam;

  fcij_vdW = 0.0;
  sfcij_vdW = 0.0;
  dfcij_vdW = 0.0;

  fnoffset = 0;
  for (ii = 0; ii < nall; ii++) {
    i = ilist[ii];
    elti = fmap[type[i]];
    if (elti < 0) continue;

    xitmp = x[i][0];
    yitmp = x[i][1];
    zitmp = x[i][2];

    for (jn = 0; jn < jnum_full[i]; jn++) {
      j = jlist_full[i][jn];
      eltj = fmap[type[j]];
      if (eltj < 0) continue;
      ij_index_full[i][j] = fnoffset + jn;

      //     First compute screening function itself, sij
      xjtmp = x[j][0];
      yjtmp = x[j][1];
      zjtmp = x[j][2];
      dij[0] = xjtmp - xitmp;
      dij[1] = yjtmp - yitmp;
      dij[2] = zjtmp - zitmp;
      rij2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
      rij = sqrt(rij2);

      if (rij > this->rc_meam) {
        fcij = 0.0;
        dfcij = 0.0;
        sij = 0.0;
        fcij_vdW = 0.0;
        dfcij_vdW = 0.0;
      } else {

        // 2 body screening
        rnorm = (this->rc_meam - rij) * drinv;
        fcij = dfcut(rnorm, dfc);
        dfcij = dfc * drinv;

        rcutBC = MAX(this->rcutBC_meambo[elti], this->rcutBC_meambo[eltj]);
        // drinv_vdW = 1.0/((1-frcut)*rcutBC);
        drinv_vdW = 1.0/(delr_rc*rcutBC);
        rnorm_vdW = (rcutBC - rij) * drinv_vdW;
        fcij_vdW = dfcut(rnorm_vdW, dfc_vdW);
        dfcij_vdW = dfc_vdW * drinv_vdW;
        dfcpair_vdW_full[fnoffset+jn] = -dfcij_vdW;


        // 3 body screening
        sij = 1.0;
        //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
        const double rbound = this->ebound_meam[elti][eltj] * rij2;
        for (kn = 0; kn < jnum_full[i]; kn++) {
          k = jlist_full[i][kn];
          eltk = fmap[type[k]];
          if (eltk < 0) continue;
          if (k == j) continue;

          delxjk = x[k][0] - xjtmp;
          delyjk = x[k][1] - yjtmp;
          delzjk = x[k][2] - zjtmp;
          rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
          if (rjk2 > rbound) continue;

          delxik = x[k][0] - xitmp;
          delyik = x[k][1] - yitmp;
          delzik = x[k][2] - zitmp;
          rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
          if (rik2 > rbound) continue;

          xik = rik2 / rij2;
          xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0) continue;

          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          Cmax = this->Cmax_meam[elti][eltj][eltk];
          Cmin = this->Cmin_meam[elti][eltj][eltk];
          if (cikj >= Cmax) continue;
          //     note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case here
          //     (other negative cikj cases were handled by the test on "a" above)
          else if (cikj <= Cmin) {
            sij = 0.0;
            break;
          } else {
            delc = Cmax - Cmin;
            cikj = (cikj - Cmin) / delc;
            sikj = fcut(cikj);
          }
          sij *= sikj;
        }
      }

      sfcij = sij * fcij;

      //! Now compute derivatives
      dscrfcn_full[fnoffset+jn] = 0.0;
      dscrfcn_vdW_full[fnoffset+jn] = 0.0;
      sfcij_vdW = sij * fcij * fcij_vdW;

      if ((!iszero(sfcij) && !iszero(sfcij - 1.0))
          ||(!iszero(sfcij_vdW) && !iszero(sfcij_vdW - 1.0))) {

        rbound = this->ebound_meam[elti][eltj] * rij2;
        for (kn = 0; kn < jnum_full[i]; kn++) {
          k = jlist_full[i][kn];
          if (k == j) continue;
          eltk = fmap[type[k]];
          if (eltk < 0) continue;

          xktmp = x[k][0];
          yktmp = x[k][1];
          zktmp = x[k][2];
          delxjk = xktmp - xjtmp;
          delyjk = yktmp - yjtmp;
          delzjk = zktmp - zjtmp;
          rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
          if (rjk2 > rbound) continue;

          delxik = xktmp - xitmp;
          delyik = yktmp - yitmp;
          delzik = zktmp - zitmp;
          rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
          if (rik2 > rbound) continue;

          xik = rik2 / rij2;
          xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0) continue;

          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          Cmax = this->Cmax_meam[elti][eltj][eltk];
          Cmin = this->Cmin_meam[elti][eltj][eltk];
          if (cikj >= Cmax) {
            continue;
            //     Note that cikj may be slightly negative (within numerical
            //     tolerance) if atoms are colinear, so don't reject that case
            //     here
            //     (other negative cikj cases were handled by the test on "a"
            //     above)
            //     Note that we never have 0<cikj<Cmin here, else sij=0
            //     (rejected above)
          } else {
            delc = Cmax - Cmin; //4.20b
            cikj = (cikj - Cmin) / delc; //arg for 4.20b
            sikj = dfcut(cikj, dfikj); //4.20b
            coef1 = dfikj / (delc * sikj);
            dCikj = dCfunc(rij2, rik2, rjk2); //4.17a rij in numerator will cancel out in eq 4.40
            dscrfcn_full[fnoffset+jn] += coef1*dCikj; //4.21 except sfcij
            dscrfcn_vdW_full[fnoffset+jn] += coef1*dCikj; // meambo

          }
        }
        coef1 = sfcij;
        coef2 = sij * dfcij / rij;
        dscrfcn_full[fnoffset+jn] = dscrfcn_full[fnoffset+jn]*coef1 - coef2; // 4.22a


        coef1 = sfcij_vdW;
        coef2 = sij * dfcij * fcij_vdW / rij;
        coef3 = sij * fcij * dfcij_vdW / rij;
        dscrfcn_vdW_full[fnoffset+jn] = dscrfcn_vdW_full[fnoffset+jn]*coef1 - coef2 -coef3;

      }
      scrfcn_full[fnoffset+jn] = sij;
      fcpair_full[fnoffset+jn] = fcij;
      fcpair_vdW_full[fnoffset+jn] = fcij_vdW;
    }
    fnoffset = fnoffset + jnum_full[i];
  }
}

//-----------------------------------------------------------------------------
// This function determines a bond between two atoms based on screening.
//
void
MEAMBO::meambo_get_bond(int* type, int* fmap, int inum, int* ilist,
      int* jnum_full, int** jlist_full, tagint* tag)
{
  int ii, i, elti, eltj, j, jn, fnoffset, nt;
  double Sij;

  this->nBC = 0;
  fnoffset = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    elti = fmap[type[i]];
    if (elti < 0) continue;
    for (jn = 0; jn < jnum_full[i]; jn++) {
      j = jlist_full[i][jn];
      eltj = fmap[type[j]];
      if (eltj < 0) continue;
      if (tag[i] > tag[j]) continue;

      Sij = scrfcn_full[fnoffset+jn]*fcpair_full[fnoffset+jn];
      if (iszero(Sij)) continue;

      // if C-C bond
      for (nt = 0; nt < this->ntypBC_meambo; nt++) {
        if ( (elti == this->ityBC_meambo[nt]) && (eltj == this->ityBC_meambo[nt]) ) {
          iBC[this->nBC][0] = i; // i1
          iBC[this->nBC][1] = j; // i2
          this->nBC++;
        }
      }
    }
    fnoffset = fnoffset + jnum_full[i];
  }
}

//-----------------------------------------------------------------------------
// This function merge two sets of neighbors of i1 and i2 with no duplicate.
//
void
MEAMBO::meambo_get_bond_neighborlist( int* type, int* fmap, int* jnum_full, int** jlist_full)
{
  int j, jBC, jindx, iB, i1, i2;

  //!  loop over the bonds to get distance information
  for (iB = 0; iB < nBC; iB++) {
    i1 = iBC[iB][0];
    i2 = iBC[iB][1];
    // if the number of neighbors of i2 is bigger than that of i1
    // add the neighbors of i2 first
    if (jnum_full[i1] < jnum_full[i2]){
      i1 = iBC[iB][1];
      i2 = iBC[iB][0];
    }
    jBC = 0;
    // add all neighors of i1
    for (j = 0; j < jnum_full[i1]; j++) {
      jindx = jlist_full[i1][j];
      if (fmap[type[jindx]] < 0) continue;
      jlistBC[iB][jBC] = jindx;
      neighbor_check[jindx] = true;
      jBC++;
    }

    // add if any neighor of i2 is not duplicate
    for (j = 0; j < jnum_full[i2]; j++) {
      jindx = jlist_full[i2][j];
      if (neighbor_check[jindx]) continue;
      if (fmap[type[jindx]] < 0) continue;
      jlistBC[iB][jBC] = jindx;
      jBC++;
    }
    // initialize
    for (j = 0; j < jBC; j++) {
      neighbor_check[jlistBC[iB][j]] = false;
    }
    jnumBC[iB] = jBC;
  }
}

//-----------------------------------------------------------------------------
// This function calculates the screening between the bond center to the neighbors.
// Bond center is a half way between the two atoms that form a bond.
// Also calculates intemediate terms for Z3 density at the end of the function.
//
void
MEAMBO::meambo_get_full_screenBC(double** x, int* jnumBC, int** jlistBC, int* type, int* fmap)
{
  int jn, j, kn, k, fnoffset, iB, i1, i2, m;
  int elti, eltj, eltk;
  double di1i2[3];
  double dij[3], rij2, rij;
  double dik[3];
  double rik2;
  double rjk2;
  double xik, xjk, sij, fcij, dfcij, sikj, dfikj, cikj;
  double Cmin, Cmax, delc, rbound, a, coef1, coef2;
  double dCikj;
  double rnorm, dfc, drinv, delr_rc;
  double rnorm_vdW, rcutBC, fcij_vdW, dfc_vdW, dfcij_vdW, drinv_vdW, coef3, sfcij_vdW;

  fnoffset = 0;
  for (iB = 0; iB < nBC; iB++) {
    i1 = iBC[iB][0];
    i2 = iBC[iB][1];
    for (m = 0; m < 3; m++)
      di1i2[m] = x[i2][m] - x[i1][m];
    elti = fmap[type[i1]];
    drinv = 1.0 / this->delr_meam;
    // double frcut = 1 - this->delr_meam/this->rc_meam;
    delr_rc = this->delr_meam/this->rc_meam;

    for (jn = 0; jn < jnumBC[iB]; jn++) {
      j = jlistBC[iB][jn];
      eltj = fmap[type[j]];
      if (eltj < 0) continue;
      if (j == i1 || j == i2) continue;

      //     First compute screening function itself, sij
      // distance vector from the bond center between i1 and i2 to j
      get_dijBC(i1, j, di1i2, x, dij);
      rij2 = dij[0]*dij[0] + dij[1]*dij[1] + dij[2]*dij[2];
      rij = sqrt(rij2);

      rcutBC = MAX(this->rcutBC_meambo[elti], this->rcutBC_meambo[eltj]);

      if (rij > rcutBC) {
        fcij = 0.0;
        dfcij = 0.0;
        sij = 0.0;
        fcij_vdW = 0.0;
        dfcij_vdW = 0.0;
      } else {

        // 2 body screening
        rnorm = (this->rc_meam - rij) * drinv;
        fcij = dfcut(rnorm, dfc);
        dfcij = dfc * drinv;

        // drinv_vdW = 1.0/((1-frcut)*rcutBC);
        drinv_vdW = 1.0/(delr_rc*rcutBC);
        rnorm_vdW = (rcutBC - rij) * drinv_vdW;
        fcij_vdW = dfcut(rnorm_vdW, dfc_vdW);
        dfcij_vdW = dfc_vdW * drinv_vdW;

        // 3 body screening
        sij = 1.0;
        //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
        rbound = this->ebound_meam[elti][eltj]*rij2;
        for (kn = 0; kn < jnumBC[iB]; kn++) {
          k = jlistBC[iB][kn];
          if (k == j) continue;
          if (k == i1 || k == i2) continue; // to exclude the atoms that form a bond
          eltk = fmap[type[k]];
          if (eltk < 0) continue;

          // the distance between the bond center and the atom k
          get_dijBC(i1, k, di1i2, x, dik);
          rik2 = dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2];
          if (rik2 > rbound) continue;
          rjk2 = rij2 + rik2 - 2.0*(dij[0]*dik[0] + dij[1]*dik[1] + dij[2]*dik[2]);
          if (rjk2 > rbound) continue;

          xik = rik2 / rij2;
          xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0) continue;

          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          Cmax = this->Cmax_meam[elti][eltj][eltk];
          Cmin = this->Cmin_meam[elti][eltj][eltk];
          if (cikj >= Cmax) continue;
          //     note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case here
          //     (other negative cikj cases were handled by the test on "a" above)
          else if (cikj <= Cmin) {
            sij = 0.0;
            break;
          } else {
            delc = Cmax - Cmin;
            cikj = (cikj - Cmin) / delc;
            sikj = fcut(cikj);
          }
          sij *= sikj;
        }
      }

      sfcij_vdW = sij*fcij*fcij_vdW;
    //! Now compute derivatives
      dscrfcnBC[fnoffset+jn] = 0.0;
      if (!iszero(sfcij_vdW) && !iszero(sfcij_vdW - 1.0)) {
        rbound = this->ebound_meam[elti][eltj] * rij2;
        for (kn = 0; kn < jnumBC[iB]; kn++) {
          k = jlistBC[iB][kn];
          if (k == i1 || k == i2) continue;
          if (k == j) continue;
          eltk = fmap[type[k]];
          if (eltk < 0) continue;

          get_dijBC(i1, k, di1i2, x, dik);
          rik2 = dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2];
          if (rik2 > rbound) continue;
          rjk2 = rij2 + rik2 - 2.0 * (dij[0] * dik[0] + dij[1] * dik[1] + dij[2] * dik[2]);
          if (rjk2 > rbound) continue;

          xik = rik2 / rij2;
          xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          //     if a < 0, then ellipse equation doesn't describe this case and
          //     atom k can't possibly screen i-j
          if (a <= 0.0) continue;

          cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
          Cmax = this->Cmax_meam[elti][eltj][eltk];
          Cmin = this->Cmin_meam[elti][eltj][eltk];
          if (cikj >= Cmax) {
            continue;
            //     Note that cikj may be slightly negative (within numerical
            //     tolerance) if atoms are colinear, so don't reject that case here
            //     (other negative cikj cases were handled by the test on "a" above)
            //     Note that we never have 0<cikj<Cmin here, else sij=0 (rejected above)
          } else {
            delc = Cmax - Cmin; //4.20b
            cikj = (cikj - Cmin) / delc; //arg for 4.20b
            sikj = dfcut(cikj, dfikj); //4.20b
            coef1 = dfikj / (delc * sikj);  // a part of 4.21 excluding dCijk
            dCikj = dCfunc(rij2, rik2, rjk2); //4.17a

            dscrfcnBC[fnoffset+jn] = dscrfcnBC[fnoffset+jn] + coef1 * dCikj; // meambo
          }
        }
        coef1 = sfcij_vdW;
        coef2 = sij * dfcij * fcij_vdW / rij;
        coef3 = sij * fcij * dfcij_vdW / rij;
        dscrfcnBC[fnoffset+jn] = dscrfcnBC[fnoffset+jn] * coef1 - coef2 -coef3;
      }

      scrfcnBC[fnoffset+jn] = sij;
      fcpairBC[fnoffset+jn] = fcij;
      fcpair_vdWBC[fnoffset+jn] = fcij_vdW;

      // intemediate terms for Z3 density
      if (!iszero(sfcij_vdW)){
        double A3bBC = sfcij_vdW/rij;
        double A3BC = sfcij_vdW/(rij2*rij);
        int nv3 = 0;
        for (m = 0; m < 3; m++) {
          arho3bBC[iB][m] +=  A3bBC*dij[m];
          for (int n = m; n < 3; n++) {
            for (int p = n; p < 3; p++) {
              arho3BC[iB][nv3] += A3BC*dij[m]*dij[n]*dij[p];
              nv3++;
            }
          }
        }
      }

    }
    fnoffset += jnumBC[iB];
  }
}

