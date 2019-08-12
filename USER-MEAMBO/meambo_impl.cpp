
#include "meambo.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MEAMBO::MEAMBO(Memory* mem)
  : memory(mem)
{
  phir = phirar = phirar1 = phirar2 = phirar3 = phirar4 = phirar5 = phirar6 = NULL;

  nmax = 0;
  rho0 = rho1 = rho2 = rho3 = frhop = NULL;
  gamma = dgamma1 = dgamma2 = dgamma3 = arho2b = NULL;
  arho1 = arho2 = arho3 = arho3b = t_ave = tsq_ave = NULL;

  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = NULL;

  neltypes = 0;
  for (int i = 0; i < maxelt; i++) {
      Z_meam[i] = A_meam[i] = rho0_meam[i] = beta0_meam[i] =
      beta1_meam[i]= beta2_meam[i] = beta3_meam[i] =
      t0_meam[i] = t1_meam[i] = t2_meam[i] = t3_meam[i] =
      rho_ref_meam[i] = ibar_meam[i] = ielt_meam[i] = 0.0;
    for (int j = 0; j < maxelt; j++) {
      lattce_meam[i][j] = FCC;
      Ec_meam[i][j] = re_meam[i][j] = alpha_meam[i][j] = delta_meam[i][j] = 
      ebound_meam[i][j] = attrac_meam[i][j] = repuls_meam[i][j] = 0.0;
      nn2_meam[i][j] = zbl_meam[i][j] = eltind[i][j] = 0;
    }
  }


  //meambo
  maxneigh_full = 0;
  nmaxBC = 0;
  maxneighBC = 0;
  maxneighBC_single  =0;
  nBC = 0;
  iBC = NULL;
  isBC = false;
  Z0BC = Z1BC = NULL;
  delta2_Z0 = ddelta2_Z0 = delta2_Z1 = ddelta2_Z1 = NULL;
  Z3BC = NULL;
  arho1BC = arho3BC = arho3bBC = NULL;
  nIBC = NULL;

  scrfcn_full = dscrfcn_full = fcpair_full = NULL;
  fcpair_vdW_full = dfcpair_vdW_full = dscrfcn_vdW_full = NULL;

  scrfcnBC = fcpairBC = fcpair_vdWBC = dscrfcnBC = NULL;

  slocal = NULL;

  dZ0BC = dZ1BC = NULL;
  dZ0jBC = dZ1jBC = NULL;
  dZ3BC = dZ3jBC = NULL;
  dnI = NULL;

  jnumBC = NULL;
  jlistBC = NULL;
  ij_index_full = NULL;
  neighbor_check = NULL;

}

MEAMBO::~MEAMBO()
{
  memory->destroy(this->phirar6);
  memory->destroy(this->phirar5);
  memory->destroy(this->phirar4);
  memory->destroy(this->phirar3);
  memory->destroy(this->phirar2);
  memory->destroy(this->phirar1);
  memory->destroy(this->phirar);
  memory->destroy(this->phir);

  memory->destroy(this->rho0);
  memory->destroy(this->rho1);
  memory->destroy(this->rho2);
  memory->destroy(this->rho3);
  memory->destroy(this->frhop);
  memory->destroy(this->gamma);
  memory->destroy(this->dgamma1);
  memory->destroy(this->dgamma2);
  memory->destroy(this->dgamma3);
  memory->destroy(this->arho2b);

  memory->destroy(this->arho1);
  memory->destroy(this->arho2);
  memory->destroy(this->arho3);
  memory->destroy(this->arho3b);
  memory->destroy(this->t_ave);
  memory->destroy(this->tsq_ave);

  memory->destroy(this->scrfcn);
  memory->destroy(this->dscrfcn);
  memory->destroy(this->fcpair);


  //meambo
  if (isBC) {
    memory->destroy(this->iBC);
    memory->destroy(this->jnumBC);
    memory->destroy(this->jlistBC);

    memory->destroy(this->Z0BC);
    memory->destroy(this->Z1BC);

    memory->destroy(this->delta2_Z0);
    memory->destroy(this->ddelta2_Z0);
    memory->destroy(this->delta2_Z1);
    memory->destroy(this->ddelta2_Z1);

    memory->destroy(this->Z3BC);
    memory->destroy(this->arho1BC);
    memory->destroy(this->arho3BC);
    memory->destroy(this->arho3bBC);

    memory->destroy(this->nIBC);
    memory->destroy(this->slocal);

    memory->destroy(this->scrfcn_full);
    memory->destroy(this->dscrfcn_full);
    memory->destroy(this->fcpair_full);
    memory->destroy(this->fcpair_vdW_full);
    memory->destroy(this->dfcpair_vdW_full);
    memory->destroy(this->dscrfcn_vdW_full);

    memory->destroy(this->scrfcnBC);
    memory->destroy(this->fcpairBC);
    memory->destroy(this->fcpair_vdWBC);
    memory->destroy(this->dscrfcnBC);

    memory->destroy(this->dZ0BC);
    memory->destroy(this->dZ1BC);
    memory->destroy(this->dZ0jBC);
    memory->destroy(this->dZ1jBC);

    memory->destroy(this->dZ3BC);
    memory->destroy(this->dZ3jBC);

    memory->destroy(this->dnI);
    memory->destroy(this->ij_index_full);
    memory->destroy(this->neighbor_check);

  }
}
