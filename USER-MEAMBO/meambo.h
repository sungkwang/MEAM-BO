#ifndef LMP_MEAMBO_H
#define LMP_MEAMBO_H

#include "memory.h"
#include <cmath>
#include <cstdlib>
#include "math_special.h"
#include "math_const.h"

#define maxelt 5
#define maxbondtype maxelt*5

namespace LAMMPS_NS {

typedef enum { FCC, BCC, HCP, DIM, DIA, DIA3, B1, C11, L12, B2, CH4, LIN, ZIG, TRI  } lattice_t;
typedef enum {DOUBLE, TRIPLE} BondOrder; // index for bond order  


class MEAMBO
{
public:
  MEAMBO(Memory* mem);
  ~MEAMBO();

private:
  Memory* memory;

  // rcut = force cutoff
  // rcutsq = force cutoff squared

  double rcut, rcutsq;

  // Ec_meam = cohesive energy
  // re_meam = nearest-neighbor distance
  // B_meam = bulk modulus
  // Z_meam = number of first neighbors for reference structure
  // ielt_meam = atomic number of element
  // A_meam = adjustable parameter
  // alpha_meam = sqrt(9*Omega*B/Ec)
  // rho0_meam = density scaling parameter
  // delta_meam = heat of formation for alloys
  // beta[0-3]_meam = electron density constants
  // t[0-3]_meam = coefficients on densities in Gamma computation
  // rho_ref_meam = background density for reference structure
  // ibar_meam(i) = selection parameter for Gamma function for elt i,
  // lattce_meam(i,j) = lattce configuration for elt i or alloy (i,j)
  // neltypes = maximum number of element type defined
  // eltind = index number of pair (similar to Voigt notation; ij = ji)
  // phir = pair potential function array
  // phirar[1-6] = spline coeffs
  // attrac_meam = attraction parameter in Rose energy
  // repuls_meam = repulsion parameter in Rose energy
  // nn2_meam = 1 if second nearest neighbors are to be computed, else 0
  // zbl_meam = 1 if zbl potential for small r to be use, else 0
  // emb_lin_neg = 1 if linear embedding function for rhob to be used, else 0
  // bkgd_dyn = 1 if reference densities follows Dynamo, else 0
  // Cmin_meam, Cmax_meam = min and max values in screening cutoff
  // rc_meam = cutoff distance for meam
  // delr_meam = cutoff region for meam
  // ebound_meam = factor giving maximum boundary of sceen fcn ellipse
  // augt1 = flag for whether t1 coefficient should be augmented
  // ialloy = flag for newer alloy formulation (as in dynamo code)
  // mix_ref_t = flag to recover "old" way of computing t in reference config
  // erose_form = selection parameter for form of E_rose function
  // gsmooth_factor = factor determining length of G smoothing region
  // vind[23]D = Voight notation index maps for 2 and 3D
  // v2D,v3D = array of factors to apply for Voight notation

  // nr,dr = pair function discretization parameters
  // nrar,rdrar = spline coeff array parameters

  // theta = angle between three atoms in line, zigzag, and trimer reference structures
  // stheta_meam = sin(theta/2) in radian used in line, zigzag, and trimer reference structures
  // ctheta_meam = cos(theta/2) in radian used in line, zigzag, and trimer reference structures

  double Ec_meam[maxelt][maxelt], re_meam[maxelt][maxelt];
  double Z_meam[maxelt];
  double A_meam[maxelt], alpha_meam[maxelt][maxelt], rho0_meam[maxelt];
  double delta_meam[maxelt][maxelt];
  double beta0_meam[maxelt], beta1_meam[maxelt];
  double beta2_meam[maxelt], beta3_meam[maxelt];
  double t0_meam[maxelt], t1_meam[maxelt];
  double t2_meam[maxelt], t3_meam[maxelt];
  double rho_ref_meam[maxelt];
  int ibar_meam[maxelt], ielt_meam[maxelt];
  lattice_t lattce_meam[maxelt][maxelt];
  int nn2_meam[maxelt][maxelt];
  int zbl_meam[maxelt][maxelt];
  int eltind[maxelt][maxelt];
  int neltypes;

  double** phir;

  double **phirar, **phirar1, **phirar2, **phirar3, **phirar4, **phirar5, **phirar6;

  double attrac_meam[maxelt][maxelt], repuls_meam[maxelt][maxelt];

  double Cmin_meam[maxelt][maxelt][maxelt];
  double Cmax_meam[maxelt][maxelt][maxelt];
  double rc_meam, delr_meam, ebound_meam[maxelt][maxelt];
  int augt1, ialloy, mix_ref_t, erose_form;
  int emb_lin_neg, bkgd_dyn;
  double gsmooth_factor;
  int vind2D[3][3], vind3D[3][3][3];
  int v2D[6], v3D[10];

  int nr, nrar;
  double dr, rdrar;
  
  //angle for trimer, zigzag, line reference structures
  double stheta_meam[maxelt][maxelt];
  double ctheta_meam[maxelt][maxelt];  

  /* parameters for unsaturated bond and non-bonded interaction */
 
  // nbondrec = number of records in MEAM-BO parameter file
  // ityBC_meambo = type of bond for each pair of atoms that forms a bond
  // rcutBC_meambo = cutoff radius for bond formation such as van der Waals radious
  // i1elt_meambo = type of 1st element that forms a bond, i.e. i1-i2
  // i2elt_meambo = type of the other element that forms a bond 
  // nb_meambo = bond type, double or triple
  // z0s_meambo = parameter to satisfy geometic condition of Z0 depending on bond/element type
  // z1s_meambo = parameter to satisfy geometic condition of Z1 depending on bond/element type
  // betas_meambo = parameters for gaussian-like function to satisfy geometric condition
  // ps_meambo = parameters for gaussian-like function to satisfy geometric condition
  // re3_meambo = experimental distance for triple bond
  // e3_meambo = polynomial parameters to fit experimental energy for a triple bond state
  // re2a_meambo = polynomial parameters to fit experimental distances for multiple double bond states
  // e2a_0_meambo = polynomial parameters to fit experimental energy for multiple double bond states 
  // e2a_1_meambo = polynomial parameters to fit experimental energy for multiple double bond states 
  // e2a_2_meambo = polynomial parameters to fit experimental energy for multiple double bond states 
  // vdW_form = selection parameter for form of vdW function
  // evdW_96LJ_meambo = epsilon (energy) parameter in 9-6 Lennard Jones equation for vdW interaction
  // svdW_96LJ_meambo = sigma (distance) paramters in 9-6 Lennard Jones equation for vdW interaction
  
  int nbondrec;  
  int ityBC_meambo[maxelt+1];
  double rcutBC_meambo[maxelt];
  
  int i1elt_meambo[maxbondtype];
  int i2elt_meambo[maxbondtype];  
  int nb_meambo[maxbondtype];
  
  int z0s_meambo[maxbondtype];
  int z1s_meambo[maxbondtype];
  double betas_meambo[maxbondtype][4];
  double ps_meambo[maxbondtype][4];
  double e3_meambo[maxbondtype][3];
  double re3_meambo[maxbondtype];
  double e2a_0_meambo[maxbondtype][3];
  double e2a_1_meambo[maxbondtype][3];
  double e2a_2_meambo[maxbondtype][3];
  double re2a_meambo[maxbondtype][3];

  int vdW_form;
  double evdW_96LJ_meambo[maxelt];
  double svdW_96LJ_meambo[maxelt];

  // intermediate variable for bond order calculation
  int BCty[maxbondtype]; 
  double betaBC[2][4], powerBC[2][4]; 

public:
  int nmax; // max number of atoms  for allocation
  int maxneigh; // max number of neighbors for allocation

  // intermediate density terms
  double *rho, *rho0, *rho1, *rho2, *rho3, *frhop;
  double *gamma, *dgamma1, *dgamma2, *dgamma3, *arho2b;
  double **arho1, **arho2, **arho3, **arho3b, **t_ave, **tsq_ave;

  // intermediate screening terms
  double *scrfcn, *dscrfcn, *fcpair;

  // meambo variables
  int maxneigh_full; // maximum number of neighbors for allocation
  int maxneighBC; // maximum number of neighors of all bonds for allocation
  int maxneighBC_single; // maximum number of neighors of a single bond for allocation
  int nmaxBC; // maximum number of bonds for allocation
  int nall; // number of local + number of ghosts

  bool isBC; // determine to execute bond order calculation

  int ntypBC_meambo; // number of elements that have unsaturated bonds
  int nBC; // number of bonds  
  int **iBC; // atom index for each bond, i1 = iBC[0], i2 = iBC[1]
  int **jlistBC, *jnumBC; // neighbor list of each bond excluding the atoms that form a bond
  bool * neighbor_check; // matrix to determine duplicates of neighbors of bond atoms
  
  int ** ij_index_full; // index matrix to convert [i][j] to [in*offset+jn]

  // intermediate screening terms for full i-j pairs
  double *scrfcn_full, *dscrfcn_full, *fcpair_full, *fcpair_vdW_full, *dfcpair_vdW_full, *dscrfcn_vdW_full;
  // intermediate screening terms for full I(bond cernter)-j pairs
  double *scrfcnBC, *fcpairBC, *fcpair_vdWBC, *dscrfcnBC ;

  // intermediate bond density terms
  double **arho1BC, **arho3BC, **arho3bBC;  
  double *Z0BC, *Z1BC, *Z3BC, *nIBC; 
  double **dZ0BC, **dZ1BC, **dZ3BC, **dnI;
  double **dZ0jBC, **dZ1jBC, **dZ3jBC;
  double *delta2_Z0, *ddelta2_Z0, *delta2_Z1, *ddelta2_Z1;
  
  double ***slocal; // local stress

protected:
  // meam_funcs.cpp

  //-----------------------------------------------------------------------------
  // Cutoff function
  //
  static double fcut(const double xi) {
    double a;
    if (xi >= 1.0)
      return 1.0;
    else if (xi <= 0.0)
      return 0.0;
    else {
      a = 1.0 - xi;
      a *= a; a *= a;
      a = 1.0 - a;
      return a * a;
    }
  }

  //-----------------------------------------------------------------------------
  // Cutoff function and derivative
  //
  static double dfcut(const double xi, double& dfc) {
    double a, a3, a4, a1m4;
    if (xi >= 1.0) {
      dfc = 0.0;
      return 1.0;
    } else if (xi <= 0.0) {
      dfc = 0.0;
      return 0.0;
    } else {
      a = 1.0 - xi;
      a3 = a * a * a;
      a4 = a * a3;
      a1m4 = 1.0-a4;

      dfc = 8 * a1m4 * a3;
      return a1m4*a1m4;
    }
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rij
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static double dCfunc(const double rij2, const double rik2, const double rjk2) {
    double rij4, a, asq, b,denom;

    rij4 = rij2 * rij2;
    a = rik2 - rjk2;
    b = rik2 + rjk2;
    asq = a*a;
    denom = rij4 - asq;
    denom = denom * denom;
    return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rik and rjk
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static void dCfunc2(const double rij2, const double rik2, const double rjk2,
               double& dCikj1, double& dCikj2) {
    double rij4, rik4, rjk4, a, denom;

    rij4 = rij2 * rij2;
    rik4 = rik2 * rik2;
    rjk4 = rjk2 * rjk2;
    a = rik2 - rjk2;
    denom = rij4 - a * a;
    denom = denom * denom;
    dCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
    dCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;
  }

  double G_gam(const double gamma, const int ibar, int &errorflag) const;
  double dG_gam(const double gamma, const int ibar, double &dG) const;
  static double zbl(const double r, const int z1, const int z2);
  static double erose(const double r, const double re, const double alpha, const double Ec, const double repuls, const double attrac, const int form);

  // static void get_shpfcn(const lattice_t latt, double (&s)[3]);
  static void get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3]);
  
  static int get_Zij(const lattice_t latt);
  // static int get_Zij2(const lattice_t latt, const double cmin, const double cmax, double &a, double &S);
  static int get_Zij2(const lattice_t latt, const double cmin, const double cmax, 
                      const double sthe, double &a, double &S);
  static int get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S);
 
  void meam_checkindex(int, int, int, int*, int*);
  void getscreen(int i, double* scrfcn, double* dscrfcn, double* fcpair, double** x, int numneigh,
      int* firstneigh, int numneigh_full, int* firstneigh_full,  int* type, int* fmap);
  void calc_rho1(int i,  int* type, int* fmap, double** x, int numneigh, int* firstneigh, int fnoffset);

  void alloyparams();
  void compute_pair_meam();
  double phi_meam(double, int, int);
  void compute_reference_density();
  void get_tavref(double*, double*, double*, double*, double*, double*, double, double, double, double,
                  double, double, double, int, int, lattice_t);
  void get_sijk(double, int, int, int, double*);
  void get_densref(double, int, int, double*, double*, double*, double*, double*, double*, double*, double*);
  void interpolate_meam(int);

  /* meambo functions */
  //-----------------------------------------------------------------------------
  // Delta(Z) function
  static double get_deltaZ(const double z, const double beta, const double power){
    // delt = exp(-beta*(abs(z)**power))
    return MathSpecial::fm_exp(-beta*pow(fabs(z), power));
  }

  //-----------------------------------------------------------------------------
  // derivative of Delta(Z) function
  static void get_delta(const double z, const double beta, const double power, double& delt, double& ddelt) {
    //      delt = exp(-beta*(abs(z)**power))
    //      ddelt = -beta*power*abs(z)**(power-1.0)*delt
    double zabs = fabs(z);
    delt = MathSpecial::fm_exp(-beta*pow(zabs, power));
    ddelt = -beta*power*pow(zabs, (power-1.0))*delt;
  }

  //-----------------------------------------------------------------------------
  // this gives vector between the bond center, formed by i1 and i2, to the its neighbor, j.
  static void get_dijBC(const int i1, const int j, const double *di1i2, double **x, double* dijBC) {
  double di1j;
  for (int m = 0; m < 3; m++){
    di1j = x[j][m] - x[i1][m];
    dijBC[m] = di1j - 0.5*di1i2[m];
  }
}

public:
  void meam_setup_global(int nelt, lattice_t* lat, double* z, int* ielement, double* atwt, double* alpha,
      double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
      double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
      int* ibar);
  void meam_setup_param(int which, double value, int nindex, int* index /*index(3)*/, int* errorflag);
  void meam_setup_done(double* cutmax);
  void meam_dens_setup(int atom_nmax, int nall, int n_neigh, int n_neigh_full);
  void meam_dens_init(int i, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
      int numneigh_full, int* firstneigh_full, int fnoffset);
  void meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
      double* eatom,  int* type, int* fmap, int& errorflag);
  void meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
      double* eatom,  int* type, int* fmap, double** x, int numneigh, int* firstneigh,
      int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom);

  // meambo functions
  void meambo_get_full_screenBC(double** x, int* jnumBC, int** jlistBC, int* type, int* fmap);
  void meambo_get_full_screen( int* type, int* fmap, double** x, int inum_half,
      int* ilist_half, int* numneigh_full, int** firstneigh_full);
  void meambo_get_bond(int* type, int* fmap, int inum_half,
      int* ilist_half, int* numneigh_full, int** firstneigh_full, tagint* tag);
  void meambo_get_bond_neighborlist( int* type, int* map, int* jnum, int** jlist);
  void meambo_dens_setup(int atom_nmax, int n_neigh, int n_single_neigh_max); 
  void meambo_energy(int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
      double* eatom,  int* type, int* fmap, int *numneigh_full, int** firstneigh_full, double** x);
  void meambo_force(int vflag_atom, int* type, int* fmap, int inum_full, int* ilist_full,
      int *numneigh_full, int** firstneigh_full, double** x, double** f, double** vatom);
  void meambo_setup_bond_params(int nelt, int lenbc, int lenarr,int* i1elt, int* i2elt, int* nb, 
      int* z0s, int* z1s, double* beta0, double* beta1, double* beta2, double* beta3,
      double* p0, double* p1, double* p2, double* p3,
      double* e2a_0_0, double* e2a_0_1, double* e2a_0_2,
      double* e2a_1_0, double* e2a_1_1, double* e2a_1_2,
      double* e2a_2_0, double* e2a_2_1, double* e2a_2_2,
      double* re2a_0, double* re2a_1, double* re2a_2,
      double* e3_0, double* e3_1, double* e3_2, double* re3);
};

// Helper functions

static inline bool iszero(const double f) {
  return fabs(f) < 1e-20;
}

template <typename TYPE, size_t maxi, size_t maxj>
static inline void setall2d(TYPE (&arr)[maxi][maxj], const TYPE v) {
  for (size_t i = 0; i < maxi; i++)
    for (size_t j = 0; j < maxj; j++)
      arr[i][j] = v;
}

template <typename TYPE, size_t maxi, size_t maxj, size_t maxk>
static inline void setall3d(TYPE (&arr)[maxi][maxj][maxk], const TYPE v) {
  for (size_t i = 0; i < maxi; i++)
    for (size_t j = 0; j < maxj; j++)
      for (size_t k = 0; k < maxk; k++)
        arr[i][j][k] = v;
}

static inline double fdiv_zero(const double n, const double d) {
  if (iszero(d))
    return 0.0;
  return n / d;
}

}
#endif
