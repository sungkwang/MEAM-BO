#include "meambo.h"
#include <cmath>
using namespace LAMMPS_NS;

//-----------------------------------------------------------------------------
// This function sets paramters values read from the MEAM paramter files
// Also initialize the paramters with default values
//
void
MEAMBO::meam_setup_global(int nelt, lattice_t* lat, double* z, int* ielement, double* atwt, double* alpha,
                        double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                        double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                        int* ibar)
{

  int i;
  double tmplat[maxelt];

  this->neltypes = nelt;

  for (i = 0; i < nelt; i++) {
    this->lattce_meam[i][i] = lat[i];

    this->Z_meam[i] = z[i];
    this->ielt_meam[i] = ielement[i];
    this->alpha_meam[i][i] = alpha[i];
    this->beta0_meam[i] = b0[i];
    this->beta1_meam[i] = b1[i];
    this->beta2_meam[i] = b2[i];
    this->beta3_meam[i] = b3[i];
    tmplat[i] = alat[i];
    this->Ec_meam[i][i] = esub[i];
    this->A_meam[i] = asub[i];
    this->t0_meam[i] = t0[i];
    this->t1_meam[i] = t1[i];
    this->t2_meam[i] = t2[i];
    this->t3_meam[i] = t3[i];
    this->rho0_meam[i] = rozero[i];
    this->ibar_meam[i] = ibar[i];

    if (this->lattce_meam[i][i] == FCC)
      this->re_meam[i][i] = tmplat[i] / sqrt(2.0);
    else if (this->lattce_meam[i][i] == BCC)
      this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 2.0;
    else if (this->lattce_meam[i][i] == HCP)
      this->re_meam[i][i] = tmplat[i];
    else if (this->lattce_meam[i][i] == DIM)
      this->re_meam[i][i] = tmplat[i];
    else if ((this->lattce_meam[i][i] == DIA) || (this->lattce_meam[i][i] == DIA3))
      this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 4.0;
    else if (this->lattce_meam[i][i] == LIN) 
      this->re_meam[i][i] = tmplat[i];  
    else if (this->lattce_meam[i][i] == ZIG) 
      this->re_meam[i][i] = tmplat[i];  
    else if (this->lattce_meam[i][i] == TRI) 
      this->re_meam[i][i] = tmplat[i];        
    else {
      //           error
    }
  }

  // Set some defaults
  this->rc_meam = 4.0;
  this->delr_meam = 0.1;
  setall2d(this->attrac_meam, 0.0);
  setall2d(this->repuls_meam, 0.0);
  setall3d(this->Cmax_meam, 2.8);
  setall3d(this->Cmin_meam, 2.0);
  setall2d(this->ebound_meam, pow(2.8, 2) / (4.0 * (2.8 - 1.0)));
  setall2d(this->delta_meam, 0.0);
  setall2d(this->nn2_meam, 0);
  setall2d(this->zbl_meam, 1);
  this->gsmooth_factor = 99.0;
  this->augt1 = 1;
  this->ialloy = 0;
  this->mix_ref_t = 0;
  this->emb_lin_neg = 0;
  this->bkgd_dyn = 0;
  this->erose_form = 0;
  // for trimer, zigzag, line refernece structure, sungkwang
  setall2d(this->stheta_meam, 1.0); // stheta = sin(theta/2*pi/180) where theta is 180, so 1.0
  setall2d(this->ctheta_meam, 0.0); // stheta = cos(theta/2*pi/180) where theta is 180, so 0

  // meam-bo paramters
  this->ntypBC_meambo = 0; // number of bond types
  this->vdW_form = 0; // van der Waal interaction form
  // bond data defaults
  for (i = 0; i < nelt; i++) {
    this->ityBC_meambo[i] = 0; // the first one carbon
    this->rcutBC_meambo[i] = this->rc_meam;
    this->evdW_96LJ_meambo[i] = 0.0;
    this->svdW_96LJ_meambo[i] = 0.0;
  }
}

//-----------------------------------------------------------------------------
// This function sets paramters values read from the MEAM bond order paramter file
// Also initialize the paramters with default values
//
void
MEAMBO::meambo_setup_bond_params(int nelt, int lenbc, int lenarr,
           int* i1elt, int* i2elt, int* nb, int* z0s, int* z1s,
           double* beta0, double* beta1, double* beta2, double* beta3,
           double* p0, double* p1, double* p2, double* p3,
           double* e2a_0_0, double* e2a_0_1, double* e2a_0_2,
           double* e2a_1_0, double* e2a_1_1, double* e2a_1_2,
           double* e2a_2_0, double* e2a_2_1, double* e2a_2_2,
           double* re2a_0, double* re2a_1, double* re2a_2,
           double* e3_0, double* e3_1, double* e3_2,
           double* re3)

{

  int i,j;

  // defaults
  for (i = 0; i < maxbondtype; i++) {
    this->i1elt_meambo[i] = 0;
    this->nb_meambo[i] = 0;
    this->z0s_meambo[i] = 0;
    this->z1s_meambo[i] = 0;
    for (j = 0; j < 4; j++) {
      this->betas_meambo[i][j] = 0.0;
      this->ps_meambo[i][j] = 0.0;
    }
    for (j = 0; j < 3; j++) {
      this->e2a_0_meambo[i][j] = 0.0;
      this->e2a_1_meambo[i][j] = 0.0;
      this->e2a_2_meambo[i][j] = 0.0;
      this->re2a_meambo[i][j] = 0.0;
      this->e3_meambo[i][j] = 0.0;
    }
    this->re3_meambo[i] = 0.0;
  }

  this->nbondrec = lenbc;

  // for (i = nelt; i < (nelt+lenbc); i++) {
  for (i = 0; i < lenbc; i++) {
    this->i1elt_meambo[i] = i1elt[i];
    this->i2elt_meambo[i] = i2elt[i];    
    
    this->nb_meambo[i] = nb[i];
    this->z0s_meambo[i] = z0s[i];
    this->z1s_meambo[i] = z0s[i];    

    this->betas_meambo[i][0] = beta0[i];
    this->betas_meambo[i][1] = beta1[i];
    this->betas_meambo[i][2] = beta2[i];
    this->betas_meambo[i][3] = beta3[i];

    this->ps_meambo[i][0] = p0[i];
    this->ps_meambo[i][1] = p1[i];
    this->ps_meambo[i][2] = p2[i];
    this->ps_meambo[i][3] = p3[i];

    this->e2a_0_meambo[i][0] = e2a_0_0[i];
    this->e2a_0_meambo[i][1] = e2a_0_1[i];
    this->e2a_0_meambo[i][2] = e2a_0_2[i];

    this->e2a_1_meambo[i][0] = e2a_1_0[i];
    this->e2a_1_meambo[i][1] = e2a_1_1[i];
    this->e2a_1_meambo[i][2] = e2a_1_2[i];

    this->e2a_2_meambo[i][0] = e2a_2_0[i];
    this->e2a_2_meambo[i][1] = e2a_2_1[i];
    this->e2a_2_meambo[i][2] = e2a_2_2[i];

    this->re2a_meambo[i][0] = re2a_0[i];
    this->re2a_meambo[i][1] = re2a_1[i];
    this->re2a_meambo[i][2] = re2a_2[i];

    this->e3_meambo[i][0] = e3_0[i];
    this->e3_meambo[i][1] = e3_1[i];
    this->e3_meambo[i][2] = e3_2[i];

    this->re3_meambo[i] = re3[i];
  }

  int BCtyp;
  for (int bo = TRIPLE; bo >= DOUBLE; --bo) {
    // BCty[bo] = ntype-1 + bo+1;
    // BCty[bo] = nelt-1 + bo+1;
    BCty[bo] = bo;
    BCtyp = BCty[bo];
    for (int m = 0; m < 4; m++)
      betaBC[bo][m] = this->betas_meambo[BCtyp][m];

    powerBC[bo][0] = this->ps_meambo[BCtyp][0];
    powerBC[bo][1] = this->ps_meambo[BCtyp][1];
    powerBC[bo][3] = this->ps_meambo[BCtyp][3];
  }

}
