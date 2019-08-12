#include "meambo.h"
#include "math_special.h"
#include <algorithm>

using namespace LAMMPS_NS;

//-----------------------------------------------------------------------------
// This function calculates forces/local stresses due to the embeding and pair energy
// pair interaction energy is also added to the final energy
//
void
MEAMBO::meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                 double* eatom,  int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                 int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom)
{
  int j, jn, k, kn, kk, m, n, p, q;
  int nv2, nv3, elti, eltj, eltk, ind;
  double xitmp, yitmp, zitmp, delij[3], rij2, rij, rij3;
  double v[6], fi[3], fj[3];
  double third, sixth;
  double pp, dUdrij, dUdsij, dUdrijm[3], force, forcem;
  double recip, phi, phip;
  double sij, scr, dscr;
  double a1, a1i, a1j, a2, a2i, a2j;
  double a3i, a3j;
  double shpi[3], shpj[3];
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  double rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  double rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  double a3, a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3], drho1drm2[3];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3], drho2drm2[3];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3], drho3drm2[3];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  double arg;
  double arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double dsij1, dsij2, force1, force2;
  double t1i, t2i, t3i, t1j, t2j, t3j;

  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  drhods1 =0;
  drhods2 =0;

  int jnindx;

  //     Compute forces atom i

  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  //     Treat each pair
  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    eltj = fmap[type[j]];
    if (isBC){
      jnindx = ij_index_full[i][j];
      scr = scrfcn_full[jnindx];
    } else {
      jnindx = fnoffset+jn;
      scr = scrfcn[jnindx];
    }

    if (!iszero(scr) && eltj >= 0) {

      delij[0] = x[j][0] - xitmp;
      delij[1] = x[j][1] - yitmp;
      delij[2] = x[j][2] - zitmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->rcutsq) {
        rij = sqrt(rij2);

        if (isBC)
          sij = scrfcn_full[jnindx]*fcpair_full[jnindx];
        else
          sij = scrfcn[jnindx]*fcpair[jnindx];

        //     Compute phi and phip
        ind = this->eltind[elti][eltj];
        pp = rij * this->rdrar;
        kk = (int)pp;
        kk = std::min(kk, this->nrar - 2);
        pp = pp - kk;
        pp = std::min(pp, 1.0);
        phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + this->phirar1[ind][kk]) * pp +
          this->phirar[ind][kk];
        phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + this->phirar4[ind][kk];
        recip = 1.0 / rij;

        if (eflag_either != 0) {
          if (eflag_global != 0)
            *eng_vdwl = *eng_vdwl + phi * sij;
          if (eflag_atom != 0) {
            eatom[i] = eatom[i] + 0.5 * phi * sij;
            eatom[j] = eatom[j] + 0.5 * phi * sij;
          }
        }

        //     Compute pair densities and derivatives
        invrei = 1.0 / this->re_meam[elti][elti];
        ai = rij * invrei - 1.0;
        ro0i = this->rho0_meam[elti];
        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai);
        drhoa0i = -this->beta0_meam[elti] * invrei * rhoa0i;
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai);
        drhoa1i = -this->beta1_meam[elti] * invrei * rhoa1i;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai);
        drhoa2i = -this->beta2_meam[elti] * invrei * rhoa2i;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai);
        drhoa3i = -this->beta3_meam[elti] * invrei * rhoa3i;

        if (elti != eltj) {
          invrej = 1.0 / this->re_meam[eltj][eltj];
          aj = rij * invrej - 1.0;
          ro0j = this->rho0_meam[eltj];
          rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj);
          drhoa0j = -this->beta0_meam[eltj] * invrej * rhoa0j;
          rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj);
          drhoa1j = -this->beta1_meam[eltj] * invrej * rhoa1j;
          rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj);
          drhoa2j = -this->beta2_meam[eltj] * invrej * rhoa2j;
          rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj);
          drhoa3j = -this->beta3_meam[eltj] * invrej * rhoa3j;
        } else {
          rhoa0j = rhoa0i;
          drhoa0j = drhoa0i;
          rhoa1j = rhoa1i;
          drhoa1j = drhoa1i;
          rhoa2j = rhoa2i;
          drhoa2j = drhoa2i;
          rhoa3j = rhoa3i;
          drhoa3j = drhoa3i;
        }

        const double t1mi = this->t1_meam[elti];
        const double t2mi = this->t2_meam[elti];
        const double t3mi = this->t3_meam[elti];
        const double t1mj = this->t1_meam[eltj];
        const double t2mj = this->t2_meam[eltj];
        const double t3mj = this->t3_meam[eltj];

        if (this->ialloy == 1) {
          rhoa1j  *= t1mj;
          rhoa2j  *= t2mj;
          rhoa3j  *= t3mj;
          rhoa1i  *= t1mi;
          rhoa2i  *= t2mi;
          rhoa3i  *= t3mi;
          drhoa1j *= t1mj;
          drhoa2j *= t2mj;
          drhoa3j *= t3mj;
          drhoa1i *= t1mi;
          drhoa2i *= t2mi;
          drhoa3i *= t3mi;
        }

        nv2 = 0;
        nv3 = 0;
        arg1i1 = 0.0;
        arg1j1 = 0.0;
        arg1i2 = 0.0;
        arg1j2 = 0.0;
        arg1i3 = 0.0;
        arg1j3 = 0.0;
        arg3i3 = 0.0;
        arg3j3 = 0.0;
        for (n = 0; n < 3; n++) {
          for (p = n; p < 3; p++) {
            for (q = p; q < 3; q++) {
              arg = delij[n] * delij[p] * delij[q] * this->v3D[nv3];
              arg1i3 = arg1i3 + arho3[i][nv3] * arg;
              arg1j3 = arg1j3 - arho3[j][nv3] * arg;
              nv3 = nv3 + 1;
            }
            arg = delij[n] * delij[p] * this->v2D[nv2];
            arg1i2 = arg1i2 + arho2[i][nv2] * arg;
            arg1j2 = arg1j2 + arho2[j][nv2] * arg;
            nv2 = nv2 + 1;
          }
          arg1i1 = arg1i1 + arho1[i][n] * delij[n];
          arg1j1 = arg1j1 - arho1[j][n] * delij[n];
          arg3i3 = arg3i3 + arho3b[i][n] * delij[n];
          arg3j3 = arg3j3 - arho3b[j][n] * delij[n];
        }

        //     rho0 terms
        drho0dr1 = drhoa0j * sij;
        drho0dr2 = drhoa0i * sij;

        //     rho1 terms
        a1 = 2 * sij / rij;
        drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1;
        drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
        for (m = 0; m < 3; m++) {
          drho1drm1[m] = a1 * rhoa1j * arho1[i][m];
          drho1drm2[m] = -a1 * rhoa1i * arho1[j][m];
        }

        //     rho2 terms
        a2 = 2 * sij / rij2;
        drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 - 2.0 / 3.0 * arho2b[i] * drhoa2j * sij;
        drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 - 2.0 / 3.0 * arho2b[j] * drhoa2i * sij;
        a2 = 4 * sij / rij2;
        for (m = 0; m < 3; m++) {
          drho2drm1[m] = 0.0;
          drho2drm2[m] = 0.0;
          for (n = 0; n < 3; n++) {
            drho2drm1[m] = drho2drm1[m] + arho2[i][this->vind2D[m][n]] * delij[n];
            drho2drm2[m] = drho2drm2[m] - arho2[j][this->vind2D[m][n]] * delij[n];
          }
          drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
          drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
        }

        //     rho3 terms
        rij3 = rij * rij2;
        a3 = 2 * sij / rij3;
        a3a = 6.0 / 5.0 * sij / rij;
        drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 - a3a * (drhoa3j - rhoa3j / rij) * arg3i3;
        drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 - a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
        a3 = 6 * sij / rij3;
        a3a = 6 * sij / (5 * rij);
        for (m = 0; m < 3; m++) {
          drho3drm1[m] = 0.0;
          drho3drm2[m] = 0.0;
          nv2 = 0;
          for (n = 0; n < 3; n++) {
            for (p = n; p < 3; p++) {
              arg = delij[n] * delij[p] * this->v2D[nv2];
              drho3drm1[m] = drho3drm1[m] + arho3[i][this->vind3D[m][n][p]] * arg;
              drho3drm2[m] = drho3drm2[m] + arho3[j][this->vind3D[m][n][p]] * arg;
              nv2 = nv2 + 1;
            }
          }
          drho3drm1[m] = (a3 * drho3drm1[m] - a3a * arho3b[i][m]) * rhoa3j;
          drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * arho3b[j][m]) * rhoa3i;
        }

        //     Compute derivatives of weighting functions t wrt rij
        t1i = t_ave[i][0];
        t2i = t_ave[i][1];
        t3i = t_ave[i][2];
        t1j = t_ave[j][0];
        t2j = t_ave[j][1];
        t3j = t_ave[j][2];

        if (this->ialloy == 1) {
          
          a1i = fdiv_zero(drhoa0j * sij, tsq_ave[i][0]);
          a1j = fdiv_zero(drhoa0i * sij, tsq_ave[j][0]);
          a2i = fdiv_zero(drhoa0j * sij, tsq_ave[i][1]);
          a2j = fdiv_zero(drhoa0i * sij, tsq_ave[j][1]);
          a3i = fdiv_zero(drhoa0j * sij, tsq_ave[i][2]);
          a3j = fdiv_zero(drhoa0i * sij, tsq_ave[j][2]);

          dt1dr1 = a1i * t1mj * (1.0 - t1i * t1mj);
          dt1dr2 = a1j * t1mi * (1.0 - t1j * t1mi);
          dt2dr1 = a2i * t2mj * (1.0 - t2i * t2mj);
          dt2dr2 = a2j * t2mi * (1.0 - t2j * t2mi);
          dt3dr1 = a3i * t3mj * (1.0 - t3i * t3mj);
          dt3dr2 = a3j * t3mi * (1.0 - t3j * t3mi); 

        } else if (this->ialloy == 2) {

          dt1dr1 = 0.0;
          dt1dr2 = 0.0;
          dt2dr1 = 0.0;
          dt2dr2 = 0.0;
          dt3dr1 = 0.0;
          dt3dr2 = 0.0;

        } else {

          ai = fdiv_zero(drhoa0j * sij, rho0[i]);
          aj = fdiv_zero(drhoa0i * sij, rho0[j]);

          dt1dr1 = ai * (t1mj - t1i);
          dt1dr2 = aj * (t1mi - t1j);
          dt2dr1 = ai * (t2mj - t2i);
          dt2dr2 = aj * (t2mi - t2j);
          dt3dr1 = ai * (t3mj - t3i);
          dt3dr2 = aj * (t3mi - t3j);
        }

        //     Compute derivatives of total density wrt rij, sij and rij(3)
        // get_shpfcn(this->lattce_meam[elti][elti], shpi);
        // get_shpfcn(this->lattce_meam[eltj][eltj], shpj);
        get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpi);
        get_shpfcn(this->lattce_meam[eltj][eltj], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpj);

        drhodr1 = dgamma1[i] * drho0dr1 +
          dgamma2[i] * (dt1dr1 * rho1[i] + t1i * drho1dr1 + dt2dr1 * rho2[i] + t2i * drho2dr1 +
                        dt3dr1 * rho3[i] + t3i * drho3dr1) -
          dgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1);
        drhodr2 = dgamma1[j] * drho0dr2 +
          dgamma2[j] * (dt1dr2 * rho1[j] + t1j * drho1dr2 + dt2dr2 * rho2[j] + t2j * drho2dr2 +
                        dt3dr2 * rho3[j] + t3j * drho3dr2) -
          dgamma3[j] * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 + shpj[2] * dt3dr2);
        for (m = 0; m < 3; m++) {
          drhodrm1[m] = 0.0;
          drhodrm2[m] = 0.0;
          drhodrm1[m] = dgamma2[i] * (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
          drhodrm2[m] = dgamma2[j] * (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
        }

        //     Compute derivatives wrt sij, but only if necessary
        if (isBC)
          dscr = dscrfcn_full[jnindx];
        else
          dscr = dscrfcn[jnindx];

        if (!iszero(dscr)) {
          drho0ds1 = rhoa0j;
          drho0ds2 = rhoa0i;
          a1 = 2.0 / rij;
          drho1ds1 = a1 * rhoa1j * arg1i1;
          drho1ds2 = a1 * rhoa1i * arg1j1;
          a2 = 2.0 / rij2;
          drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arho2b[i] * rhoa2j;
          drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arho2b[j] * rhoa2i;
          a3 = 2.0 / rij3;
          a3a = 6.0 / (5.0 * rij);
          drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
          drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

          if (this->ialloy == 1) {
            a1i = fdiv_zero(rhoa0j, tsq_ave[i][0]);
            a1j = fdiv_zero(rhoa0i, tsq_ave[j][0]);
            a2i = fdiv_zero(rhoa0j, tsq_ave[i][1]);
            a2j = fdiv_zero(rhoa0i, tsq_ave[j][1]);
            a3i = fdiv_zero(rhoa0j, tsq_ave[i][2]);
            a3j = fdiv_zero(rhoa0i, tsq_ave[j][2]);

            dt1ds1 = a1i * t1mj * (1.0 - t1i * t1mj);
            dt1ds2 = a1j * t1mi * (1.0 - t1j * t1mi);
            dt2ds1 = a2i * t2mj * (1.0 - t2i * t2mj);
            dt2ds2 = a2j * t2mi * (1.0 - t2j * t2mi);
            dt3ds1 = a3i * t3mj * (1.0 - t3i * t3mj);
            dt3ds2 = a3j * t3mi * (1.0 - t3j * t3mi);

          } else if (this->ialloy == 2) {

            dt1ds1 = 0.0;
            dt1ds2 = 0.0;
            dt2ds1 = 0.0;
            dt2ds2 = 0.0;
            dt3ds1 = 0.0;
            dt3ds2 = 0.0;

          } else {

            ai = fdiv_zero(rhoa0j, rho0[i]);
            aj = fdiv_zero(rhoa0i, rho0[j]);

            dt1ds1 = ai * (t1mj - t1i);
            dt1ds2 = aj * (t1mi - t1j);
            dt2ds1 = ai * (t2mj - t2i);
            dt2ds2 = aj * (t2mi - t2j);
            dt3ds1 = ai * (t3mj - t3i);
            dt3ds2 = aj * (t3mi - t3j);
          }

          drhods1 = dgamma1[i] * drho0ds1 +
            dgamma2[i] * (dt1ds1 * rho1[i] + t1i * drho1ds1 + dt2ds1 * rho2[i] + t2i * drho2ds1 +
                          dt3ds1 * rho3[i] + t3i * drho3ds1) -
            dgamma3[i] * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1);
          drhods2 = dgamma1[j] * drho0ds2 +
            dgamma2[j] * (dt1ds2 * rho1[j] + t1j * drho1ds2 + dt2ds2 * rho2[j] + t2j * drho2ds2 +
                          dt3ds2 * rho3[j] + t3j * drho3ds2) -
            dgamma3[j] * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 + shpj[2] * dt3ds2);
        }

        //     Compute derivatives of energy wrt rij, sij and rij[3]
        dUdrij = phip * sij + frhop[i] * drhodr1 + frhop[j] * drhodr2;
        dUdsij = 0.0;
        if (!iszero(dscr)) {
          dUdsij = phi + frhop[i] * drhods1 + frhop[j] * drhods2;
        }
        for (m = 0; m < 3; m++) {
          dUdrijm[m] = frhop[i] * drhodrm1[m] + frhop[j] * drhodrm2[m];
        }

        //     Add the part of the force due to dUdrij and dUdsij
        // 1/rij or recip for dUdsij is canceled out with rij in eq 4.17 or dCfunc
        force = dUdrij * recip + dUdsij * dscr;
        for (m = 0; m < 3; m++) {
          forcem = delij[m] * force + dUdrijm[m];
          f[i][m] = f[i][m] + forcem;
          f[j][m] = f[j][m] - forcem;
        }

        //     Tabulate per-atom virial as symmetrized stress tensor

        if (vflag_atom != 0) {
          fi[0] = delij[0] * force + dUdrijm[0];
          fi[1] = delij[1] * force + dUdrijm[1];
          fi[2] = delij[2] * force + dUdrijm[2];
          v[0] = -0.5 * (delij[0] * fi[0]);
          v[1] = -0.5 * (delij[1] * fi[1]);
          v[2] = -0.5 * (delij[2] * fi[2]);
          v[3] = -0.25 * (delij[0] * fi[1] + delij[1] * fi[0]);
          v[4] = -0.25 * (delij[0] * fi[2] + delij[2] * fi[0]);
          v[5] = -0.25 * (delij[1] * fi[2] + delij[2] * fi[1]);

          for (m = 0; m < 6; m++) {
            vatom[i][m] = vatom[i][m] + v[m];
            vatom[j][m] = vatom[j][m] + v[m];
          }
        }

        //     Now compute forces on other atoms k due to change in sij

        if (iszero(sij) || iszero(sij - 1.0)) continue; //: cont jn loop

        double dxik(0), dyik(0), dzik(0);
        double dxjk(0), dyjk(0), dzjk(0);

        for (kn = 0; kn < numneigh_full; kn++) {
          k = firstneigh_full[kn];
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            double xik, xjk, cikj, sikj, dfc, a;
            double dCikj1, dCikj2;
            double delc, rik2, rjk2;

            const double Cmax = this->Cmax_meam[elti][eltj][eltk];
            const double Cmin = this->Cmin_meam[elti][eltj][eltk];

            dsij1 = 0.0;
            dsij2 = 0.0;
            if (!iszero(sij) && !iszero(sij - 1.0)) {
              const double rbound = rij2 * this->ebound_meam[elti][eltj];
              delc = Cmax - Cmin;
              dxjk = x[k][0] - x[j][0];
              dyjk = x[k][1] - x[j][1];
              dzjk = x[k][2] - x[j][2];
              rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
              if (rjk2 <= rbound) {
                dxik = x[k][0] - x[i][0];
                dyik = x[k][1] - x[i][1];
                dzik = x[k][2] - x[i][2];
                rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
                if (rik2 <= rbound) {
                  xik = rik2 / rij2;
                  xjk = rjk2 / rij2;
                  a = 1 - (xik - xjk) * (xik - xjk);
                  if (!iszero(a)) {
                    cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                    if (cikj >= Cmin && cikj <= Cmax) {
                      cikj = (cikj - Cmin) / delc;
                      sikj = dfcut(cikj, dfc);
                      dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                      a = sij / delc * dfc / sikj;
                      dsij1 = a * dCikj1;
                      dsij2 = a * dCikj2;
                    }
                  }
                }
              }
            }

            if (!iszero(dsij1) || !iszero(dsij2)) {
              force1 = dUdsij * dsij1;
              force2 = dUdsij * dsij2;

              f[i][0] += force1 * dxik;
              f[i][1] += force1 * dyik;
              f[i][2] += force1 * dzik;
              f[j][0] += force2 * dxjk;
              f[j][1] += force2 * dyjk;
              f[j][2] += force2 * dzjk;
              f[k][0] -= force1 * dxik + force2 * dxjk;
              f[k][1] -= force1 * dyik + force2 * dyjk;
              f[k][2] -= force1 * dzik + force2 * dzjk;

              //     Tabulate per-atom virial as symmetrized stress tensor

              if (vflag_atom != 0) {
                fi[0] = force1 * dxik;
                fi[1] = force1 * dyik;
                fi[2] = force1 * dzik;
                fj[0] = force2 * dxjk;
                fj[1] = force2 * dyjk;
                fj[2] = force2 * dzjk;
                v[0] = -third * (dxik * fi[0] + dxjk * fj[0]);
                v[1] = -third * (dyik * fi[1] + dyjk * fj[1]);
                v[2] = -third * (dzik * fi[2] + dzjk * fj[2]);
                v[3] = -sixth * (dxik * fi[1] + dxjk * fj[1] + dyik * fi[0] + dyjk * fj[0]);
                v[4] = -sixth * (dxik * fi[2] + dxjk * fj[2] + dzik * fi[0] + dzjk * fj[0]);
                v[5] = -sixth * (dyik * fi[2] + dyjk * fj[2] + dzik * fi[1] + dzjk * fj[1]);

                for (m = 0; m < 6; m++) {
                  vatom[i][m] = vatom[i][m] + v[m];
                  vatom[j][m] = vatom[j][m] + v[m];
                  vatom[k][m] = vatom[k][m] + v[m];
                }
              }
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// This function calculates forces/local stresses due to the bond order
//
void
MEAMBO::meambo_force(int vflag_atom, int* type, int* fmap,
      int inum, int* ilist, int* jnum_full, int** jlist_full, double** x, double** f, double** vatom)
{
  int j, jn, k, kn, m, n, p, q;
  int bo, iB, itBC, BCtyp, i,ii, i1, i2, fnoffset, jnindx;
  int nv2, nv3, elti, eltj, eltk, ind, elti1, elti2;
  int i1i2_ind, i2i1_ind;

  double xitmp, yitmp, zitmp, dij[3], rij2, rij, rij3;
  double forcem;
  double a1;
  double dsij1, dsij2;
  double Sij, Z0_i;
  double arg1i1BC, dZ1dr1, dZ1ds1, dZ1drm1[3];
  double twooverrij;
  double di1j[3], ri1j, di2j[3], ri2j;
  double Si1j, Si2j;
  double dik[3], djk[3];
  double xik, xjk, cikj, sikj, dfc, a;
  double dCikj1, dCikj2, delc, rik2, rjk2;
  double rbound, Cmin, Cmax;
  double dscr;
  double di1i2[3], ri1i2;
  double dri1i2[3];
  double Si1i2;
  double dSi1i2_i1[3], dSi1i2_i2[3];
  double dZ0i1[3], dZ0i2[3], dZ1i1[3], dZ1i2[3];
  double ri1i2sq;
  double Ebond[2], fbond[2];
  double BOI, BOI_2, BOI_22, re2, eb2[3], deltaZ3[2], E21;
  double re3, eb3[3], rIre, exp_rIre, E_3;
  double delZ1i1, delZ1i2;
  double pE3prI, pE2prI;
  double enbond;
  double rscrn_diff, rscrn_diff_E_vdW;
  double dIjBC[3];
  double dZ0ji1[3], dZ0ji2[3], dZ1ji1[3], dZ1ji2[3];
  double dSij[3], dSji[3], drscrni1i2_vdW[3];

  // dZ0 and dZ1 calculation
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
      jnindx = ij_index_full[i][j];
      Sij = scrfcn_full[jnindx]*fcpair_full[jnindx]*fcpair_vdW_full[jnindx];
      if (iszero(Sij)) continue;

      dij[0] = x[j][0] - xitmp;
      dij[1] = x[j][1] - yitmp;
      dij[2] = x[j][2] - zitmp;
      rij2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
      rij = sqrt(rij2);
      twooverrij = 2.0/rij;

      arg1i1BC = 0.0;
      for (m = 0; m < 3; m++)
        arg1i1BC += arho1BC[i][m]*dij[m];

      dscr = dscrfcn_vdW_full[jnindx];
      if (!iszero(dscr)) {
        for (m = 0; m < 3; m++) {
          dSij[m] = dscr*dij[m];
          dSji[m] = -dSij[m];
        }
        dZ1ds1 = twooverrij*arg1i1BC;

        rbound = rij2 * this->ebound_meam[elti][eltj];
        for (kn = 0; kn < jnum_full[i]; kn++) {
          k = jlist_full[i][kn];
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            for (m = 0; m < 3; m++)
              djk[m] = x[k][m] - x[j][m];
            rjk2 = djk[0] * djk[0] + djk[1] * djk[1] + djk[2] * djk[2];
            if (rjk2 <= rbound) {
              for (m = 0; m < 3; m++)
                dik[m] = x[k][m] - x[i][m];
              rik2 = dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2];
              if (rik2 <= rbound) {
                xik = rik2 / rij2;
                xjk = rjk2 / rij2;
                a = 1 - (xik - xjk) * (xik - xjk);
                if (!iszero(a)) {
                  cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                  Cmax = this->Cmax_meam[elti][eltj][eltk];
                  Cmin = this->Cmin_meam[elti][eltj][eltk];
                  if (cikj >= Cmin && cikj <= Cmax) {
                    delc = Cmax - Cmin;
                    cikj = (cikj - Cmin) / delc;
                    sikj = dfcut(cikj, dfc);
                    dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                    a = Sij / delc * dfc / sikj;
                    dsij1 = a * dCikj1;
                    dsij2 = a * dCikj2;
                    if (!iszero(dsij1) || !iszero(dsij2)) {
                      for (m = 0; m < 3; m++) {
                        dSij[m] = dSij[m] + dsij1*dik[m];
                        dSji[m] = dSji[m] + dsij2*djk[m];
                        drscrni1i2_vdW[m] = dsij1*dik[m] + dsij2*djk[m];
                        dZ0jBC[fnoffset+kn][m] -= drscrni1i2_vdW[m];
                        dZ1jBC[fnoffset+kn][m] -= dZ1ds1*drscrni1i2_vdW[m];
                      }
                    }
                  }
                }
              }
            }
          }
        } // end of for k
        for (m = 0; m < 3; m++) {
          dZ0BC[i][m] += dSij[m];
          dZ1BC[i][m] += dZ1ds1*dSij[m];
          dZ0jBC[fnoffset+jn][m] += dSji[m];
          dZ1jBC[fnoffset+jn][m] += dZ1ds1*dSji[m];
        }
      }
      a1 = twooverrij*Sij; //2.0/rij*Sij;
      for (m = 0; m < 3; m++) {
        dZ1drm1[m] = a1*arho1BC[i][m];
      }

      a1 = a1/rij2; //2.0/rij*Sij/rij2;
      dZ1dr1 = -a1*arg1i1BC;
      for (m = 0; m < 3; m++) {
        forcem = dZ1dr1*dij[m] + dZ1drm1[m];
        dZ1BC[i][m] += forcem;
        dZ1jBC[fnoffset+jn][m] -= forcem;
      }

    }
    fnoffset =  fnoffset + jnum_full[i];
  }


  double a3, a3a, arg;
  double arg1i3BC, arg3i3BC;
  double dZ3drm[3], dZ3dr, dZ3ds;
  // dZ3BC calculation
  fnoffset = 0;
  for (iB = 0; iB < nBC; iB++) {
    i1 = iBC[iB][0];
    i2 = iBC[iB][1];
    elti = fmap[type[i1]];
    for (m = 0; m < 3; m++)
      di1i2[m] = x[i2][m] - x[i1][m];
    for (jn = 0; jn < jnumBC[iB]; jn++) {
      j = jlistBC[iB][jn];
      if (j == i1 || j == i2) continue; // to exclude the atoms that form a bond
      Sij = scrfcnBC[fnoffset+jn]*fcpairBC[fnoffset+jn]*fcpair_vdWBC[fnoffset+jn];
      if(iszero(Sij)) continue;
      eltj = fmap[type[j]];

      // distance vector from the bond center between i1 and i2 to j
      get_dijBC(i1, j, di1i2, x, dij);
      rij2 = dij[0]*dij[0] + dij[1]*dij[1] + dij[2]*dij[2];
      rij = sqrt(rij2);
      rij3 = rij*rij2;

      nv3 = 0;
      arg1i3BC = 0.0; arg3i3BC = 0.0;
      for (n = 0; n < 3; n++) {
        for (p = n; p < 3; p++) {
          for (q = p; q < 3; q++) {
            arg = dij[n]*dij[p]*dij[q]*this->v3D[nv3];
            arg1i3BC = arg1i3BC + arho3BC[iB][nv3] * arg;
            nv3 = nv3 + 1;
          }
        }
        arg3i3BC = arg3i3BC + arho3bBC[iB][n]*dij[n];
      }
      dscr = dscrfcnBC[fnoffset + jn];
      if (!iszero(dscr)) {
        for (m = 0; m < 3; m++) {
          dSij[m] = dscr*dij[m];
          dSji[m] = -dSij[m];
        }

        a3 = 2.0/rij3;
        a3a = 1.2/rij; // 6.0/5.0
        dZ3ds = a3*arg1i3BC - a3a*arg3i3BC;

        rbound = rij2 * this->ebound_meam[elti][eltj];
        for (kn = 0; kn < jnumBC[iB]; kn++) {
          k = jlistBC[iB][kn];
          if (k == i1 || k == i2) continue;
          if (k == j) continue;
          get_dijBC(i1, k, di1i2, x, dik);
          rik2 = dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2];
          if (rik2 > rbound) continue;
          for (m = 0; m < 3; m++) djk[m] = x[k][m] - x[j][m];
          rjk2 = djk[0] * djk[0] + djk[1] * djk[1] + djk[2] * djk[2];
          if (rjk2 > rbound) continue;
          xik = rik2 / rij2; xjk = rjk2 / rij2;
          a = 1 - (xik - xjk) * (xik - xjk);
          if (!iszero(a)) {
            eltk = fmap[type[k]];
            Cmax = this->Cmax_meam[elti][eltj][eltk];
            Cmin = this->Cmin_meam[elti][eltj][eltk];
            cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
            if (cikj >= Cmin && cikj <= Cmax) {
              delc = Cmax - Cmin;
              cikj = (cikj - Cmin) / delc;
              sikj = dfcut(cikj, dfc);
              dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
              a = Sij / delc * dfc / sikj;
              dsij1 = a * dCikj1;
              dsij2 = a * dCikj2;
              if (!iszero(dsij1) || !iszero(dsij2)) {
                for (m = 0; m < 3; m++) {
                  dSij[m] += dsij1*dik[m];
                  dSji[m] += dsij2*djk[m];
                  drscrni1i2_vdW[m] = dsij1*dik[m] + dsij2*djk[m];
                  dZ3jBC[fnoffset+kn][m] -= dZ3ds*drscrni1i2_vdW[m];
                }
              }
            }
          }

        } //     end of k loop
        for (m = 0; m < 3; m++) {
          dZ3BC[iB][m] += 0.5*dZ3ds*dSij[m];
          dZ3jBC[fnoffset+jn][m] += dZ3ds*dSji[m];
        }
      } // end of if

      a3 = -6.0*Sij/(rij2*rij3);
      a3a = -1.2*Sij/rij3;
      dZ3dr = (a3*arg1i3BC - a3a*arg3i3BC);

      a3 = 6.0*Sij/rij3;
      a3a = 1.2*Sij/rij;
      for (m = 0; m < 3; m++) {
        dZ3drm[m] = 0.0;
        nv2 = 0;
        for (n = 0; n < 3; n++) {
          for (p = n; p < 3; p++) {
            arg = dij[n]*dij[p]*this->v2D[nv2];
            dZ3drm[m] += arho3BC[iB][this->vind3D[m][n][p]]*arg;
            nv2 = nv2 + 1;
          }
        }
        dZ3drm[m] = a3*dZ3drm[m] - a3a*arho3bBC[iB][m];
      }
      for (m = 0; m < 3; m++) {
        forcem = dZ3dr*dij[m] + dZ3drm[m];
        dZ3BC[iB][m] += 0.5*forcem;
        dZ3jBC[fnoffset+jn][m] -= forcem;
      }
    }
    fnoffset += jnumBC[iB];
  }


  double rscrn_vdW;
  double pfI3pSi1i2, pfI3pD30i1_dD30i1, pfI3pD30i2_dD30i2;
  double pfI2pfI3, pfI3pD33I_dD3I;
  double pfI2pSi1i2;
  double pfI2pD20i1_dD20i1, pfI2pD20i2_dD20i2, pfI2pD21i1_dD21i1, pfI2pD21i2_dD21i2;
  double pfI2pD23I_dD23I, pfI2pD33I_dD33I;
  double pE2pall;
  double vdW_all;

  // final force due to bond order
  fnoffset = 0;
  for (iB = 0; iB < nBC; iB++) {
    i1 = iBC[iB][0];
    i2 = iBC[iB][1];

    elti1 = fmap[type[i1]];
    elti2 = fmap[type[i2]];

    for (m = 0; m < 3; m++)
      di1i2[m] = x[i2][m] - x[i1][m];
    ri1i2sq = di1i2[0] * di1i2[0] + di1i2[1] * di1i2[1] + di1i2[2] * di1i2[2];
    ri1i2 = sqrt(ri1i2sq);

    i1i2_ind = ij_index_full[i1][i2];
    i2i1_ind = ij_index_full[i2][i1];

    for (m = 0; m < 3; m++) {
      dZ0i1[m] = dZ0jBC[i1i2_ind][m];
      dZ1i1[m] = dZ1jBC[i1i2_ind][m];
      dZ0i2[m] = dZ0jBC[i2i1_ind][m];
      dZ1i2[m] = dZ1jBC[i2i1_ind][m];
    }

    Si1i2 = scrfcn_full[i1i2_ind]*fcpair_full[i1i2_ind];
    rscrn_vdW = fcpair_vdW_full[i1i2_ind];
    rscrn_diff = 1.0 - rscrn_vdW;

    for (m = 0; m < 3; m++) {
      dSi1i2_i1[m] = dscrfcn_full[i1i2_ind]*di1i2[m];
      dSi1i2_i2[m] = -dSi1i2_i1[m];
      dri1i2[m] = di1i2[m]/ri1i2;
    }

    // nI and dnI
    enbond = nIBC[iB];

    if (enbond > 2.0) {
      int in;
      double pnIZ0j, pnIZ1j;
      double pnIZ0i, pnIZ1i;
      double delta2_Z0Z1j, delta2_Z0Z1i;
      double Z0_j, Sji, dscrfcn_ji, pZji[3], pZi[3], drscrnji_vdW[3];
    // clear dnI every bond operation
      #pragma omp simd // for vectorization with no dependency
      for (i = 0; i < nall; i++) dnI[i][0] = dnI[i][1] = dnI[i][2] = 0.0;

    //! nI' = pnI/pSi1j*Si1j' + pnI/pSi2j*Si2j' + pnI/pZj(0)Zj(0)' + pnI/pZj(1)Zj(1)'
    //! contributions from the j neighbors of i1
      bo = DOUBLE;
      for (jn = 0; jn < jnum_full[i1]; jn++) {
        j = jlist_full[i1][jn];
        eltj = fmap[type[j]];
        if (eltj<0) continue;
        ind = ij_index_full[i1][j];
        Si1j = scrfcn_full[ind]*fcpair_full[ind]*fcpair_vdW_full[ind];
        if (iszero(Si1j)) continue;
        delta2_Z0Z1j = delta2_Z0[j]*delta2_Z1[j];
        Z0_j = 2*(Z0BC[j] - this->z0s_meambo[BCty[bo]]);
        pnIZ0j = Si1j*ddelta2_Z0[j]*Z0_j*delta2_Z1[j];
        pnIZ1j = Si1j*delta2_Z0[j]*ddelta2_Z1[j];
        for (in = 0; in < jnum_full[j]; in++) {
          i = jlist_full[j][in];
          elti = fmap[type[i]];
          if (elti<0) continue;
          ind = ij_index_full[j][i];
          Sji = scrfcn_full[ind]*fcpair_full[ind]*fcpair_vdW_full[ind];
          dscrfcn_ji = dscrfcn_vdW_full[ind]; // or dscrfcn_vdW[jindx_half]
          for (m = 0; m < 3; m++) pZji[m] = pnIZ0j*dZ0jBC[ind][m] + pnIZ1j*dZ1jBC[ind][m];

          Z0_i = 2*(Z0BC[i] - this->z0s_meambo[BCty[bo]]);
          delta2_Z0Z1i = delta2_Z0[i]*delta2_Z1[i];
          pnIZ0i = Sji*ddelta2_Z0[i]*Z0_i*delta2_Z1[i];
          pnIZ1i = Sji*delta2_Z0[i]*ddelta2_Z1[i];
          for (m = 0; m < 3; m++) pZi[m] = pnIZ0i*dZ0BC[i][m] + pnIZ1i*dZ1BC[i][m];
          if (!iszero(dscrfcn_ji)) {
            for (m = 0; m < 3; m++) dij[m] = x[j][m] - x[i][m];
            rij2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
            rbound = rij2 * this->ebound_meam[elti][eltj];
            for (m = 0; m < 3; m++) dSji[m] = -dscrfcn_ji*dij[m];
            for (kn = 0; kn < jnum_full[j]; kn++) {
              k = jlist_full[j][kn];
              eltk = fmap[type[k]];
              if (k == j || eltk < 0) continue;
              for (m = 0; m < 3; m++) djk[m] = x[k][m] - x[j][m];
              rjk2 = djk[0]*djk[0] + djk[1]*djk[1] + djk[2]*djk[2];
              if (rjk2 > rbound) continue;
              for (m = 0; m < 3; m++) dik[m] = x[k][m] - x[i][m];
              rik2 = dik[0]*dik[0] + dik[1]*dik[1] + dik[2]*dik[2];
              if (rik2 > rbound) continue;
              xik = rik2/rij2; xjk = rjk2/rij2;
              a = 1 - (xik - xjk)*(xik - xjk);
              if (iszero(a)) continue;
              cikj = (2.0*(xik + xjk) + a - 2.0)/a;
              Cmax = this->Cmax_meam[elti][eltj][eltk];
              Cmin = this->Cmin_meam[elti][eltj][eltk];
              if (cikj < Cmin || cikj > Cmax) continue;
              delc = Cmax - Cmin;
              cikj = (cikj - Cmin)/delc;
              sikj = dfcut(cikj, dfc);
              dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
              a = Sji/delc*dfc/sikj;
              dsij1 = a*dCikj1; dsij2 = a*dCikj2;
              if (!iszero(dsij1) || !iszero(dsij2)) {
                for (m = 0; m < 3; m++) dSji[m] = dSji[m] + dsij2*djk[m];
                if (i == i1){
                  for (m = 0; m < 3; m++) drscrnji_vdW[m] = dsij1*dik[m] + dsij2*djk[m];
                  if (j == i2) // i1-i2 case
                    for (m = 0; m < 3; m++) dnI[k][m] = dnI[k][m] - delta2_Z0Z1i*drscrnji_vdW[m];
                  else  // i1-j except i2 case
                    for (m = 0; m < 3; m++) dnI[k][m] = dnI[k][m] - delta2_Z0Z1j*drscrnji_vdW[m];
                }
              }
            } // end of for k
            if ( (i == i1) && (j == i2) ) { // i1-i2 case
              arg = delta2_Z0Z1i + delta2_Z0Z1j;
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + arg*dSji[m];
            } else if ( (i != i1) &&  (j == i2) ) { // j-i2 (excluding i1) case
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + delta2_Z0Z1i*dSji[m];
            } else if ( (i == i1) && (j != i2) ){ // i1-j(excluding i2) case
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + delta2_Z0Z1j*dSji[m];
            }
          }

          if ( (i == i1) && (j == i2) ) { // i1-i2 case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZi[m] + pZji[m];
          } else if ( (i != i1) &&  (j == i2) ) { // j-i2 (excluding i1) case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZi[m] + pZji[m];
          } else if ( (i == i1) && (j != i2) ){ // i1-j(excluding i2) case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZji[m];
          } else if ( (i != i1) && (j != i2) ) { // j neighbors of bond iB
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZji[m];
          }
        } // ii
      } // j

      //! contributions from the j neighbors of i2
      for (jn = 0; jn < jnum_full[i2]; jn++) {
        j = jlist_full[i2][jn];
        eltj = fmap[type[j]];
        if (eltj<0) continue;
        ind = ij_index_full[i2][j];
        Si2j = scrfcn_full[ind]*fcpair_full[ind]*fcpair_vdW_full[ind];
        if (iszero(Si2j)) continue;
        delta2_Z0Z1j = delta2_Z0[j]*delta2_Z1[j];
        Z0_j = 2*(Z0BC[j] - this->z0s_meambo[BCty[bo]]);
        pnIZ0j = Si2j*ddelta2_Z0[j]*Z0_j*delta2_Z1[j];
        pnIZ1j = Si2j*delta2_Z0[j]*ddelta2_Z1[j];
        for (in = 0; in < jnum_full[j]; in++) {
          i = jlist_full[j][in];
          elti = fmap[type[i]];
          if (elti<0) continue;
          ind = ij_index_full[j][i];
          Sji = scrfcn_full[ind]*fcpair_full[ind]*fcpair_vdW_full[ind];
          dscrfcn_ji = dscrfcn_vdW_full[ind];
          for (m = 0; m < 3; m++) pZji[m] = pnIZ0j*dZ0jBC[ind][m] + pnIZ1j*dZ1jBC[ind][m];
          Z0_i = 2*(Z0BC[i] - this->z0s_meambo[BCty[bo]]);
          delta2_Z0Z1i = delta2_Z0[i]*delta2_Z1[i];
          pnIZ0i = Sji*ddelta2_Z0[i]*Z0_i*delta2_Z1[i];
          pnIZ1i = Sji*delta2_Z0[i]*ddelta2_Z1[i];
          for (m = 0; m < 3; m++) pZi[m] = pnIZ0i*dZ0BC[i][m] + pnIZ1i*dZ1BC[i][m];
          if (!iszero(dscrfcn_ji)) {
            for (m = 0; m < 3; m++) dij[m] = x[j][m] - x[i][m];
            rij2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
            rbound = rij2 * this->ebound_meam[elti][eltj];
            for (m = 0; m < 3; m++) dSji[m] = -dscrfcn_ji*dij[m];
            for (kn = 0; kn < jnum_full[j]; kn++) {
              k = jlist_full[j][kn];
              eltk = fmap[type[k]];
              if (k == j || eltk < 0) continue;
              for (m = 0; m < 3; m++) djk[m] = x[k][m] - x[j][m];
              rjk2 = djk[0]*djk[0] + djk[1]*djk[1] + djk[2]*djk[2];
              if (rjk2 > rbound) continue;
              for (m = 0; m < 3; m++) dik[m] = x[k][m] - x[i][m];
              rik2 = dik[0]*dik[0] + dik[1]*dik[1] + dik[2]*dik[2];
              if (rik2 > rbound) continue;
              xik = rik2/rij2; xjk = rjk2/rij2;
              a = 1 - (xik - xjk)*(xik - xjk);
              if (iszero(a)) continue;
              cikj = (2.0*(xik + xjk) + a - 2.0)/a;
              Cmax = this->Cmax_meam[elti][eltj][eltk];
              Cmin = this->Cmin_meam[elti][eltj][eltk];
              if (cikj < Cmin || cikj > Cmax) continue;
              delc = Cmax - Cmin;
              cikj = (cikj - Cmin)/delc;
              sikj = dfcut(cikj, dfc);
              dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
              a = Sji/delc*dfc/sikj;
              dsij1 = a*dCikj1; dsij2 = a*dCikj2;
              if (!iszero(dsij1) || !iszero(dsij2)) {
                for (m = 0; m < 3; m++) dSji[m] = dSji[m] + dsij2*djk[m];
                if (i == i2){
                  for (m = 0; m < 3; m++) drscrnji_vdW[m] = dsij1*dik[m] + dsij2*djk[m];
                  if (j == i1) // i2-i1 case
                    for (m = 0; m < 3; m++) dnI[k][m] = dnI[k][m] - delta2_Z0Z1i*drscrnji_vdW[m];
                  else  // i2-j except i1 case
                    for (m = 0; m < 3; m++) dnI[k][m] = dnI[k][m] - delta2_Z0Z1j*drscrnji_vdW[m];
                }
              }
            } // end of for k
            if ( (i == i2) && (j == i1) ) { // i2-i1 case
              double temp = delta2_Z0Z1i + delta2_Z0Z1j;
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + temp*dSji[m];
            } else if ( (i != i2) &&  (j == i1) ) { // j-i1 (excluding i2) case
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + delta2_Z0Z1i*dSji[m];
            } else if ( (i == i2) && (j != i1) ){ // i2-j(excluding i1) case
              for (m = 0; m < 3; m++) dnI[j][m] = dnI[j][m] + delta2_Z0Z1j*dSji[m];
            }
          }
          if ( (i == i2) && (j == i1) ) { // i1-i2 case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZi[m] + pZji[m];
          } else if ( (i != i2) &&  (j == i1) ) { // j-i2 (excluding i1) case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZi[m] + pZji[m];
          } else if ( (i == i2) && (j != i1) ){ // i1-j(excluding i2) case
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZji[m];
          } else if ( (i != i2) && (j != i1) ) { // j neighbors of bond iB
            for (m = 0; m < 3; m++) dnI[i][m] = dnI[i][m] + pZji[m];
          }
        } // ii
      } // j
    }

    // vdW contribution to energy
    itBC = this->ityBC_meambo[0]; // the first one, carbon


    if ((this->vdW_form == 1) && (!iszero(rscrn_diff)) && (!iszero(this->evdW_96LJ_meambo[itBC]))) {
      double rc = this->svdW_96LJ_meambo[itBC];
      double rcrI = rc/ri1i2;
      double rcrI6 = rcrI*rcrI*rcrI*rcrI*rcrI*rcrI;
      double rcrI9 = rcrI6 * rcrI*rcrI*rcrI;
      double E_vdW = this->evdW_96LJ_meambo[itBC]*(2*rcrI9 - 3*rcrI6);
      double dE_vdW = -18*(this->evdW_96LJ_meambo[itBC]/ri1i2)*(rcrI9 - rcrI6);
      double dfci1i2_vdW = dfcpair_vdW_full[i1i2_ind]; //for vdW interaction

      rscrn_diff_E_vdW = rscrn_diff*E_vdW;

      double pEpSi1i2 = -E_vdW;
      double pEpE_vdWpE_vdWprI = rscrn_diff*dE_vdW;
      vdW_all = pEpSi1i2*dfci1i2_vdW + pEpE_vdWpE_vdWprI;

    } else {
      rscrn_diff_E_vdW = 0.0;
      vdW_all = 0.0;
    }

    double Z0i1, Z0i2, delZ0i1, delZ0i2, ddelZ0i1, ddelZ0i2, ddeltaZ3;
    /* triple bond */
    bo = TRIPLE;
    BCtyp = BCty[bo];
    Z0i1 = Z0BC[i1] - this->z0s_meambo[BCtyp];
    Z0i2 = Z0BC[i2] - this->z0s_meambo[BCtyp];
    get_delta( Z0i1*Z0i1, this->betaBC[bo][0], this->powerBC[bo][0], delZ0i1, ddelZ0i1);
    get_delta( Z0i2*Z0i2, this->betaBC[bo][0], this->powerBC[bo][0], delZ0i2, ddelZ0i2);
    get_delta( Z3BC[iB], this->betaBC[bo][3], this->powerBC[bo][3], deltaZ3[bo], ddeltaZ3);
    re3 = this->re3_meambo[BCtyp];
    for (m = 0; m < 3; m++) eb3[m] = this->e3_meambo[BCtyp][m];
    rIre = ri1i2/re3-1.0;
    exp_rIre = MathSpecial::fm_exp(-this->betaBC[bo][2]*rIre);
    E_3 = eb3[0]*(1.0 + rIre*eb3[1] + rIre*rIre*eb3[2])*exp_rIre;
    Ebond[bo] = E_3 + rscrn_diff_E_vdW;
    fbond[bo] = Si1i2*(delZ0i1*delZ0i2)*deltaZ3[bo];

// fI3' = pfI3/pSi1i2*Si1i2' + pfI3/pD3(0)i1*D3(0)i1' + pfI3/pD3(0)i2*D3(0)i2' + pfI3/pD3(3)I*D3(3)I'
    pE3prI = -this->betaBC[bo][2]/re3*E_3 + eb3[0]/re3*(eb3[1] + 2.0*eb3[2]*rIre) * exp_rIre;
    pfI3pSi1i2 = delZ0i1*delZ0i2*deltaZ3[bo];
    double pfI3pD30i1 = Si1i2*delZ0i2*deltaZ3[bo];
    double pfI3pD30i2 = Si1i2*delZ0i1*deltaZ3[bo];
    double dD30i1 = ddelZ0i1*2.0*Z0i1;
    double dD30i2 = ddelZ0i2*2.0*Z0i2;
    pfI3pD30i1_dD30i1 = pfI3pD30i1*dD30i1;
    pfI3pD30i2_dD30i2 = pfI3pD30i2*dD30i2;
    pfI2pfI3 = 1.0 - deltaZ3[bo];
    double pfI3pD33I = Si1i2*delZ0i1*delZ0i2;
    double dD33I = ddeltaZ3;
    pfI3pD33I_dD3I = pfI3pD33I*dD33I;
    /* end of triple bond */

    /* double bond */
    bo = DOUBLE;
    BCtyp = BCty[bo];

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
      BOI_2 = 0.0;
      re2 = this->re2a_meambo[BCtyp][0];
      for (m = 0; m < 3; m++)
        eb2[m] = this->e2a_0_meambo[BCtyp][m];
    }

    delZ0i1 = delta2_Z0[i1];
    delZ0i2 = delta2_Z0[i2];
    get_delta( Z3BC[iB], this->betaBC[bo][3], this->powerBC[bo][3], deltaZ3[bo], ddeltaZ3);

    delZ1i1 = delta2_Z1[i1];
    delZ1i2 = delta2_Z1[i2];

    rIre = ri1i2/re2-1.0;
    exp_rIre = MathSpecial::fm_exp(-this->betaBC[bo][2]*rIre);
    E21 = ( 1.0 + rIre*eb2[1] + rIre*rIre*eb2[2] )*exp_rIre;
    Ebond[bo] = eb2[0]*E21 + rscrn_diff_E_vdW;
    fbond[bo] = Si1i2*(delZ0i1*delZ0i2)*deltaZ3[bo]*(delZ1i1*delZ1i2)
                + fbond[TRIPLE]*(1.0 - deltaZ3[TRIPLE]);
// fI2' = pfI2/pSi1i2*Si1i2' + pfI2/pD2(0)i1*D2(0)i1' + pfI2/pD2(0)i2*D2(0)i2'
//      + pfI2/pD2(1)i1*D2(1)i1' + pfI2/pD2(1)i2*D2(1)i2' + pfI2/pD2(3)I*D2(3)I'
//      + pfI2/pfI3*fI3' + pfI2/pD3(3)I*D3(3)I'

    pfI2pSi1i2 = delZ0i1*delZ0i2*delZ1i1*delZ1i2*deltaZ3[bo];
    double pfI2pD20i1 = Si1i2*delZ0i2*delZ1i1*delZ1i2*deltaZ3[bo];
    double pfI2pD20i2 = Si1i2*delZ0i1*delZ1i1*delZ1i2*deltaZ3[bo];
    double pfI2pD21i1 = Si1i2*delZ0i1*delZ0i2*delZ1i2*deltaZ3[bo];
    double pfI2pD21i2 = Si1i2*delZ0i1*delZ0i2*delZ1i1*deltaZ3[bo];
    Z0i1 = Z0BC[i1] - this->z0s_meambo[BCtyp];
    Z0i2 = Z0BC[i2] - this->z0s_meambo[BCtyp];
    double dD20i1 = ddelta2_Z0[i1]*2.0*Z0i1;
    double dD20i2 = ddelta2_Z0[i2]*2.0*Z0i2;
    double dD21i1 = ddelta2_Z1[i1];
    double dD21i2 = ddelta2_Z1[i2];

    pfI2pD20i1_dD20i1 = pfI2pD20i1*dD20i1;
    pfI2pD20i2_dD20i2 = pfI2pD20i2*dD20i2;
    pfI2pD21i1_dD21i1 = pfI2pD21i1*dD21i1;
    pfI2pD21i2_dD21i2 = pfI2pD21i2*dD21i2;

    double pfI2pD23I = Si1i2*delZ0i1*delZ0i2*delZ1i1*delZ1i2;
    double pfI2pD33I = -fbond[TRIPLE];
    double dD23I = ddeltaZ3;
    pfI2pD23I_dD23I = pfI2pD23I*dD23I;
    pfI2pD33I_dD33I = pfI2pD33I*dD33I;

    double beta2 = this->betas_meambo[BCtyp][2];
    pE2prI = eb2[0]*(eb2[1]/re2 + 2.0*eb2[2]*rIre/re2)*exp_rIre + eb2[0]*E21*(-beta2/re2);
    if (enbond > 2.0) {
  //! E21(iB) = (1 + rIre*eb2(1,iB) + rIr22*eb2(2,iB))*exp_rIre
      double pE2pe0 = E21;
      double pE2pe1 = eb2[0]*rIre*exp_rIre;
      double pE2pe2 = eb2[0]*rIre*rIre*exp_rIre;
      double rIr2_2 = ri1i2/(re2*re2);
      double pE2pr2 = eb2[0]*E21*beta2*rIr2_2 - eb2[0]*(eb2[1]*rIr2_2 + 2.0*eb2[2]*rIre*rIr2_2)*exp_rIre;

  //! ek' = pekpBOI*BOI'
      double pe0pBOI = this->e2a_1_meambo[BCtyp][0] + 2.0*this->e2a_2_meambo[BCtyp][0]*BOI_2;
      double pe1pBOI = this->e2a_1_meambo[BCtyp][1] + 2.0*this->e2a_2_meambo[BCtyp][1]*BOI_2;
      double pe2pBOI = this->e2a_1_meambo[BCtyp][2] + 2.0*this->e2a_2_meambo[BCtyp][2]*BOI_2;

  //! re2' = pr2/pBOI*BOI'
      double pr2pBOI = this->re2a_meambo[BCtyp][1] + 2.0*this->re2a_meambo[BCtyp][2]*BOI_2;
  //! BOI' = pBOI/pnI*nI'
      double pBOIpnI = -2.0/(enbond*enbond);
      pE2pall = pBOIpnI*(pE2pe0*pe0pBOI + pE2pe1*pe1pBOI + pE2pe2*pe2pBOI + pE2pr2*pr2pBOI);
    } else {
      pE2pall = 0;
    }

    /* end of double bond */

    // i1-j iteration
    for (jn = 0; jn < jnumBC[iB]; jn++) {
      j = jlistBC[iB][jn];
      if (j == i1 || j == i2) continue; // to exclude the atoms that form a bond
      eltj = fmap[type[j]];

      for (m = 0; m < 3; m++) {
        di1j[m] = x[j][m] - x[i1][m];
        di2j[m] = x[j][m] - x[i2][m];
      }

      ri1j = di1j[0]*di1j[0] + di1j[1]*di1j[1] + di1j[2]*di1j[2];
      ri2j = di2j[0]*di2j[0] + di2j[1]*di2j[1] + di2j[2]*di2j[2];

      get_dijBC(i1, j, di1i2, x, dIjBC);

      // derivative Si1i2 w.r.t. j
      double dSi1i2j[3] = {0};
      if (!iszero(Si1i2) && !iszero(Si1i2 - 1.0)) {
        rbound = ri1i2sq * this->ebound_meam[elti1][elti2];
        for (m = 0; m < 3; m++) djk[m] = x[j][m] - x[i2][m];
        rjk2 = djk[0]*djk[0] + djk[1]*djk[1] + djk[2]*djk[2];
        if (rjk2 <= rbound) {
          for (m = 0; m < 3; m++) dik[m] = x[j][m] - x[i1][m];
          rik2 = dik[0]*dik[0] + dik[1]*dik[1] + dik[2]*dik[2];
          if (rik2 <= rbound) {
            xik = rik2/ri1i2sq; xjk = rjk2/ri1i2sq;
            a = 1 - (xik - xjk)*(xik - xjk);
            if (!iszero(a)) {
              cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
              Cmax = this->Cmax_meam[elti1][elti2][eltj];
              Cmin = this->Cmin_meam[elti1][elti2][eltj];
              if (cikj >= Cmin && cikj <= Cmax) {
                delc = Cmax - Cmin;
                cikj = (cikj - Cmin)/delc;
                sikj = dfcut(cikj, dfc);
                dCfunc2(ri1i2sq, rik2, rjk2, dCikj1, dCikj2);
                a = Si1i2/delc*dfc/sikj;
                dsij1 = a*dCikj1; dsij2 = a*dCikj2;
                if (!iszero(dsij1) || !iszero(dsij2)) {
                  for (m = 0; m < 3; m++) {
                    dSi1i2_i1[m] = dSi1i2_i1[m] + dsij1*dik[m];
                    dSi1i2_i2[m] = dSi1i2_i2[m] + dsij2*djk[m];
                    dSi1i2j[m] = -(dsij1*dik[m] + dsij2*djk[m]);
                  }
                }
              }
            }
          }
        }
      }

      if (ri1j < this->rcutsq){
        int ind = ij_index_full[i1][j];
        for (m = 0; m < 3; m++)  dZ0ji1[m] = dZ0jBC[ind][m];
        for (m = 0; m < 3; m++)  dZ1ji1[m] = dZ1jBC[ind][m];
      } else {
        for (m = 0; m < 3; m++)  dZ0ji1[m] = 0.0;
        for (m = 0; m < 3; m++)  dZ1ji1[m] = 0.0;
      }
      if (ri2j < this->rcutsq){
        int ind = ij_index_full[i2][j];
        for (m = 0; m < 3; m++)  dZ0ji2[m] = dZ0jBC[ind][m];
        for (m = 0; m < 3; m++)  dZ1ji2[m] = dZ1jBC[ind][m];
      } else {
        for (m = 0; m < 3; m++)  dZ0ji2[m] = 0.0;
        for (m = 0; m < 3; m++)  dZ1ji2[m] = 0.0;
      }

      double dfbond, df3j[3], dEbond, pfI2pD21i1i2, pfI2pD20i1i2, pfI3pD30i1i2;
      bo = TRIPLE;
      for (m = 0; m < 3; m++) {
        dEbond = 0;
        pfI3pD30i1i2 = pfI3pD30i1_dD30i1*dZ0ji1[m] + pfI3pD30i2_dD30i2*dZ0ji2[m];
        dfbond = pfI3pSi1i2*dSi1i2j[m] + pfI3pD30i1i2 + pfI3pD33I_dD3I*dZ3jBC[fnoffset+jn][m];
        df3j[m] = dfbond;
        forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
        f[j][m] += forcem;
        for (n = 0; n < 3; n++) slocal[j][m][n] += forcem*dIjBC[n];
      }
      bo = DOUBLE;
      for (m = 0; m < 3; m++) {
        dEbond = 0;
        if (enbond > 2) dEbond = pE2pall*dnI[j][m];
        pfI2pD20i1i2 = pfI2pD20i1_dD20i1*dZ0ji1[m] + pfI2pD20i2_dD20i2*dZ0ji2[m];
        pfI2pD21i1i2 = pfI2pD21i1_dD21i1*dZ1ji1[m] + pfI2pD21i2_dD21i2*dZ1ji2[m];
        dfbond = pfI2pSi1i2*dSi1i2j[m] + pfI2pD20i1i2 + pfI2pD21i1i2 + pfI2pfI3*df3j[m]
                + pfI2pD23I_dD23I*dZ3jBC[fnoffset+jn][m] + pfI2pD33I_dD33I*dZ3jBC[fnoffset+jn][m];
        forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
        f[j][m] += forcem;
        for (n = 0; n < 3; n++) slocal[j][m][n] += forcem*dIjBC[n];
      }

    }

    // assign force due to the bond between atom i1-i2
    double dfbond, df3i1[3], df3i2[3], dEbond, pfI2pD21i1i2, pfI2pD20i1i2, pfI3pD30i1i2;
    bo = TRIPLE;
    arg = pE3prI + vdW_all;
    for (m = 0; m < 3; m++) {
      dEbond = arg*dri1i2[m];
      pfI3pD30i1i2 = pfI3pD30i1_dD30i1*dZ0BC[i1][m] + pfI3pD30i2_dD30i2*dZ0i2[m];
      dfbond = pfI3pSi1i2*dSi1i2_i1[m] + pfI3pD30i1i2 + pfI3pD33I_dD3I*dZ3BC[iB][m];
      df3i1[m] = dfbond;
      forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
      f[i1][m] += forcem;
      for (n = 0; n < 3; n++) slocal[i1][m][n] -= 0.5*forcem*di1i2[n];
    }
    for (m = 0; m < 3; m++) {
      dEbond = -arg*dri1i2[m];
      pfI3pD30i1i2 = pfI3pD30i2_dD30i2*dZ0BC[i2][m] + pfI3pD30i1_dD30i1*dZ0i1[m];
      dfbond = pfI3pSi1i2*dSi1i2_i2[m] + pfI3pD30i1i2 + pfI3pD33I_dD3I*dZ3BC[iB][m];
      df3i2[m] = dfbond;
      forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
      f[i2][m] += forcem;
      for (n = 0; n < 3; n++) slocal[i2][m][n] += 0.5*forcem*di1i2[n]; // di2i1
    }
    bo = DOUBLE;
    arg = pE2prI + vdW_all;
    for (m = 0; m < 3; m++) {
      dEbond = arg*dri1i2[m];
      if (enbond > 2) dEbond = dEbond + pE2pall*dnI[i1][m];
      pfI2pD20i1i2 = pfI2pD20i1_dD20i1*dZ0BC[i1][m] + pfI2pD20i2_dD20i2*dZ0i2[m];
      pfI2pD21i1i2 = pfI2pD21i1_dD21i1*dZ1BC[i1][m] + pfI2pD21i2_dD21i2*dZ1i2[m];
      dfbond = pfI2pSi1i2*dSi1i2_i1[m] + pfI2pD20i1i2 + pfI2pD21i1i2 + pfI2pfI3*df3i1[m]
            + pfI2pD23I_dD23I*dZ3BC[iB][m] + pfI2pD33I_dD33I*dZ3BC[iB][m];
      forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
      f[i1][m] += forcem;
      for (n = 0; n < 3; n++) slocal[i1][m][n] -= 0.5*forcem*di1i2[n];
    }
    for (m = 0; m < 3; m++) {
      dEbond = -arg*dri1i2[m];
      if (enbond > 2) dEbond = dEbond + pE2pall*dnI[i2][m];
      pfI2pD20i1i2 = pfI2pD20i2_dD20i2*dZ0BC[i2][m] + pfI2pD20i1_dD20i1*dZ0i1[m];
      pfI2pD21i1i2 = pfI2pD21i2_dD21i2*dZ1BC[i2][m] + pfI2pD21i1_dD21i1*dZ1i1[m];
      dfbond = pfI2pSi1i2*dSi1i2_i2[m] + pfI2pD20i1i2 + pfI2pD21i1i2 + pfI2pfI3*df3i2[m]
            + pfI2pD23I_dD23I*dZ3BC[iB][m] + pfI2pD33I_dD33I*dZ3BC[iB][m];
      forcem = dfbond*Ebond[bo] + fbond[bo]*dEbond;
      f[i2][m] += forcem;
      for (n = 0; n < 3; n++) slocal[i2][m][n] += 0.5*forcem*di1i2[n]; // di2i1
    }
    fnoffset =  fnoffset + jnumBC[iB];
  }

  if (vflag_atom != 0) {     // Copy slocal(3,3,.) to vatom(6,.)
    for (ii = 0; ii < nall; ii++) {
      vatom[ii][0] = vatom[ii][0] + slocal[ii][0][0];
      vatom[ii][1] = vatom[ii][1] + slocal[ii][1][1];
      vatom[ii][2] = vatom[ii][2] + slocal[ii][2][2];
      vatom[ii][3] = vatom[ii][3] + 0.5 * (slocal[ii][0][1] + slocal[ii][1][0]);
      vatom[ii][4] = vatom[ii][4] + 0.5 * (slocal[ii][0][2] + slocal[ii][2][0]);
      vatom[ii][5] = vatom[ii][5] + 0.5 * (slocal[ii][1][2] + slocal[ii][2][1]);
    }
  }
}

