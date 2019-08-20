/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Michael Baskes (Mississippi State University)
   Contributing author: Ricolindo L Carino (Mississippi State University)
   Contributing author: Sungkwang Mun (Mississippi State University)
   based on the work of Greg Wagner (SNL) and Sebastian Hütter (OvGU)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "meambo.h"
#include "pair_meambo.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

static const int nkeywords = 27;
static const char *keywords[] = {
  "Ec","alpha","rho0","delta","lattce",
  "attrac","repuls","nn2","Cmin","Cmax","rc","delr",
  "augt1","gsmooth_factor","re","ialloy",
  "mixture_ref_t","erose_form","zbl",
  "emb_lin_neg","bkgd_dyn", "theta",
  "ntypBC","rcutBC","vdW_form","evdW_96LJ","svdW_96LJ"}; 



/* ---------------------------------------------------------------------- */

PairMEAMBO::PairMEAMBO(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  allocated = 0;

  nelements = 0;
  elements = NULL;
  mass = NULL;
  meambo_inst = new MEAMBO(memory);

  // set comm size needed by this Pair
  // the number of variables per atom to be communicated

  comm_forward = 47;
  comm_reverse = 34;

}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMEAMBO::~PairMEAMBO()
{
  delete meambo_inst;

  for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  delete [] mass;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairMEAMBO::compute(int eflag, int vflag)
{
  int i, ii, num_all_half_neighbors, num_all_full_neighbors, num_max_neighbors,  num_all_neighbors, errorflag;

  int inum_half, *ilist_half, *numneigh_half, **firstneigh_half;
  int inum_full, *ilist_full, *numneigh_full, **firstneigh_full;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // neighbor list info

  inum_half = listhalf->inum;
  ilist_half = listhalf->ilist;
  numneigh_half = listhalf->numneigh;
  firstneigh_half = listhalf->firstneigh;

  inum_full = listfull->inum;
  ilist_full = listfull->ilist;
  numneigh_full = listfull->numneigh;
  firstneigh_full = listfull->firstneigh;

  // check size of scrfcn based on half neighbor list

  int nlocal = atom->nlocal;
  meambo_inst->nall = nlocal + atom->nghost;

  // strip neighbor lists of any special bond flags before using with MEAM
  // necessary before doing neigh_f2c and neigh_c2f conversions each step

  if (neighbor->ago == 0) {
    neigh_strip(inum_half,ilist_half,numneigh_half,firstneigh_half);
    neigh_strip(meambo_inst->nall,ilist_full,numneigh_full,firstneigh_full);
  }

  num_all_half_neighbors = 0;
  for (ii = 0; ii < inum_half; ii++)
    num_all_half_neighbors += numneigh_half[ilist_half[ii]];

  num_all_full_neighbors = 0;
  if (meambo_inst->isBC) {
    for (ii = 0; ii < meambo_inst->nall; ii++) {
      num_all_full_neighbors += numneigh_full[ilist_full[ii]];
    }
  }

  meambo_inst->meam_dens_setup(atom->nmax, meambo_inst->nall, num_all_half_neighbors, num_all_full_neighbors);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  // int ntype = atom->ntypes;
  tagint *tag = atom->tag;

  // 3 stages of MEAM calculation
  // loop over my atoms followed by communication

  errorflag = 0;

  if (meambo_inst->isBC) {
    meambo_inst->meambo_get_full_screen(type,map,x,inum_full, ilist_full,
            numneigh_full,firstneigh_full);
  }

  int offset = 0;
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meambo_inst->meam_dens_init(i, type, map, x, numneigh_half[i], firstneigh_half[i],
                    numneigh_full[i], firstneigh_full[i], offset);
    offset += numneigh_half[i];

  }

  if (meambo_inst->isBC) {
    meambo_inst->meambo_get_bond(type,map,inum_full, ilist_full,
          numneigh_full,firstneigh_full, tag);
  }

  if (meambo_inst->isBC && meambo_inst->nBC>0) {
    // need to allocate bond order terms that needs the number of bonds, meambo_inst->nBC
    num_all_neighbors = 0;
    num_max_neighbors = 0;
    for (int iB = 0; iB < meambo_inst->nBC; iB++) {
      int i1 = meambo_inst->iBC[iB][0];
      num_all_neighbors += numneigh_full[i1];
      if (num_max_neighbors < numneigh_full[i1])
        num_max_neighbors = numneigh_full[i1];
    }
    num_all_neighbors = num_all_neighbors*2;
    num_max_neighbors = num_max_neighbors*2;

    // if number of the bonds is larger than the max number of atoms,
    // then use number of bonds to relocate. Otherwise use nmax.
    double bond_nmax;
    if (atom->nmax<meambo_inst->nBC)
      bond_nmax = meambo_inst->nBC;
    else
      bond_nmax = atom->nmax;

    meambo_inst->meambo_dens_setup(bond_nmax, num_all_neighbors, num_max_neighbors);
    meambo_inst->meambo_get_bond_neighborlist(type, map, numneigh_full, firstneigh_full);

     //     Compute screening function and derivatives for unsaturated bonds
    meambo_inst->meambo_get_full_screenBC(x, meambo_inst->jnumBC, meambo_inst->jlistBC, type, map);
  }

  //Reverse communicate (accumulate) partial ghost atom values to local atoms
  comm->reverse_comm_pair(this);

  meambo_inst->meam_dens_final(nlocal,eflag_either,eflag_global,eflag_atom,
              &eng_vdwl,eatom,type,map,errorflag);

  if (errorflag) {
    char str[128];
    sprintf(str,"MEAM library error %d",errorflag);
    error->one(FLERR,str);
  }

  //Forward communicate (copy) full local atom values to ghost atoms
  comm->forward_comm_pair(this);

  // vptr is first value in vatom if it will be used by meam_force()
  // else vatom may not exist, so pass dummy ptr
  double **vptr;
  if (vflag_atom) vptr = vatom;
  else vptr = NULL;

  offset = 0;
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meambo_inst->meam_force(i,eflag_either,eflag_global,eflag_atom,
                vflag_atom,&eng_vdwl,eatom,type,map,x,
                numneigh_half[i],firstneigh_half[i],
                numneigh_full[i],firstneigh_full[i],
                offset,f,vptr);
    offset += numneigh_half[i];
  }

  if (meambo_inst->isBC && meambo_inst->nBC>0) {
    meambo_inst->meambo_energy(eflag_either, eflag_global, eflag_atom,
                &eng_vdwl, eatom, type, map, numneigh_full, firstneigh_full, x);
    meambo_inst->meambo_force(vflag_atom, type, map, inum_full, ilist_full,
                numneigh_full, firstneigh_full, x,f,vptr);
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMEAMBO::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMEAMBO::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMEAMBO::coeff(int narg, char **arg)
{
  int i,j,m,n;

  if (!allocated) allocate();

  // Add argument for bond file info (can be NULL for MEAM-only)
  if (narg < 7) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read MEAM element names between 2 filenames
  // nelements = # of MEAM elements
  // elements = list of unique element names

  if (nelements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
    delete [] mass;
  }
  nelements = narg - 5 - atom->ntypes;

  // printf("narg:%d, ntype:%d, nelt:%d\n", narg, atom->ntypes, nelements);

  if (nelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");
  elements = new char*[nelements];
  mass = new double[nelements];

  for (i = 0; i < nelements; i++) {
    n = strlen(arg[i+3]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],arg[i+3]);
  }

  // read MEAM library (1st), parameter (2nd) files, and BO paramters.
  // pass all parameters to MEAM-BO package
  // tell MEAM-BO package that setup is done

  read_files(arg[2],arg[2+nelements+1],arg[narg-1]);
  meambo_inst->meam_setup_done(&cutmax);

  // read args that map atom types to MEAM elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (i = 4 + nelements; i < (narg-1); i++) {
    m = i - (4+nelements) + 1;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass for i,i in atom class

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,mass[map[i]]);
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMEAMBO::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style MEAMBO requires newton pair on");

  // need full and half neighbor list

  int irequest_full = neighbor->request(this,instance_me);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
  neighbor->requests[irequest_full]->ghost = 1;
  int irequest_half = neighbor->request(this,instance_me); // no ghost for half neighbor
  neighbor->requests[irequest_half]->id = 2;

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */

void PairMEAMBO::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = ptr;
  else if (id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMEAMBO::init_one(int i, int j)
{
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairMEAMBO::read_files(char *globalfile, char *userfile, char *BCfile)
{
  // open global meamf file on proc 0

  FILE *fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(globalfile);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open 1st MEAM paramter file %s",globalfile);
      error->one(FLERR,str);
    }
  }

  // allocate parameter arrays

  int params_per_line = 19;

  lattice_t *lat = new lattice_t[nelements];
  int *ielement = new int[nelements];
  int *ibar = new int[nelements];
  double *z = new double[nelements];
  double *atwt = new double[nelements];
  double *alpha = new double[nelements];
  double *b0 = new double[nelements];
  double *b1 = new double[nelements];
  double *b2 = new double[nelements];
  double *b3 = new double[nelements];
  double *alat = new double[nelements];
  double *esub = new double[nelements];
  double *asub = new double[nelements];
  double *t0 = new double[nelements];
  double *t1 = new double[nelements];
  double *t2 = new double[nelements];
  double *t3 = new double[nelements];
  double *rozero = new double[nelements];

  bool *found = new bool[nelements];
  for (int i = 0; i < nelements; i++) found[i] = false;

  // read each set of params from global MEAM file
  // one set of params can span multiple lines
  // store params if element name is in element list
  // if element name appears multiple times, only store 1st entry

  int i,n,nwords;
  char **words = new char*[80];
  char line[MAXLINE],copy[MAXLINE],*ptr;
  int eof = 0;

  int nset = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in 1st MEAM paramter file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    nwords = 0;
    words[nwords++] = strtok(line,"' \t\n\r\f");
    while ((words[nwords++] = strtok(NULL,"' \t\n\r\f"))) continue;

    // skip if element name isn't in element list

    for (i = 0; i < nelements; i++)
      if (strcmp(words[0],elements[i]) == 0) break;
    if (i >= nelements) continue;

    // skip if element already appeared

    if (found[i] == true) continue;
    found[i] = true;

    // map lat string to an integer
    // reference structures availble for unary interaction
    if (strcmp(words[1],"fcc") == 0) lat[i] = FCC;
    else if (strcmp(words[1],"bcc") == 0) lat[i] = BCC;
    else if (strcmp(words[1],"hcp") == 0) lat[i] = HCP;
    else if (strcmp(words[1],"dim") == 0) lat[i] = DIM;
    else if (strcmp(words[1],"dia") == 0) lat[i] = DIA;
    else if (strcmp(words[1],"dia3") == 0) lat[i] = DIA3;
    else if (strcmp(words[1],"lin") == 0) lat[i] = LIN; 
    else if (strcmp(words[1],"zig") == 0) lat[i] = ZIG; 
    else if (strcmp(words[1],"tri") == 0) lat[i] = TRI; 
    else {
      char str[512];
      snprintf(str,512, "Unrecognized lattice type '%s' in 1st MEAM paramter file. " 
                        "Available types are fcc, bcc, hcp, dim, dia, dia3, lin, zig, and tri.", 
                        words[1]);
      error->all(FLERR, str);
    }

    // store parameters

    z[i] = atof(words[2]);
    ielement[i] = atoi(words[3]);
    atwt[i] = atof(words[4]);
    alpha[i] = atof(words[5]);
    b0[i] = atof(words[6]);
    b1[i] = atof(words[7]);
    b2[i] = atof(words[8]);
    b3[i] = atof(words[9]);
    alat[i] = atof(words[10]);
    esub[i] = atof(words[11]);
    asub[i] = atof(words[12]);
    t0[i] = atof(words[13]);
    t1[i] = atof(words[14]);
    t2[i] = atof(words[15]);
    t3[i] = atof(words[16]);
    rozero[i] = atof(words[17]);
    ibar[i] = atoi(words[18]);

    nset++;
  }

  // error if didn't find all elements in file

  if (nset != nelements)
    error->all(FLERR,"Did not find all elements in 1st MEAM paramter file");

  // pass element parameters to MEAM package

  meambo_inst->meam_setup_global(nelements,lat,z,ielement,atwt,alpha,b0,b1,b2,b3,
                       alat,esub,asub,t0,t1,t2,t3,rozero,ibar);

  // set element masses

  for (i = 0; i < nelements; i++) mass[i] = atwt[i];

  // clean-up memory


  delete [] lat;
  delete [] ielement;
  delete [] ibar;
  delete [] z;
  delete [] atwt;
  delete [] alpha;
  delete [] b0;
  delete [] b1;
  delete [] b2;
  delete [] b3;
  delete [] alat;
  delete [] esub;
  delete [] asub;
  delete [] t0;
  delete [] t1;
  delete [] t2;
  delete [] t3;
  delete [] rozero;
  delete [] found;

  // done if user param file is NULL

  if (strcmp(userfile,"NULL") == 0) return;

  // open user param file on proc 0

  if (comm->me == 0) {
    fp = force->open_potential(userfile);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open 2nd MEAM paramter file %s",userfile);
      error->one(FLERR,str);
    }
  }

  // read settings
  // pass them one at a time to MEAM package
  // match strings to list of corresponding ints

  int which;
  double value;
  int nindex,index[3];
  int maxparams = 6;
  char **params = new char*[maxparams];
  int nparams;

  eof = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nparams = atom->count_words(line);
    if (nparams == 0) continue;

    // words = ptrs to all words in line

    nparams = 0;
    params[nparams++] = strtok(line,"=(), '\t\n\r\f");
    while (nparams < maxparams &&
           (params[nparams++] = strtok(NULL,"=(), '\t\n\r\f")))
      continue;
    nparams--;

    for (which = 0; which < nkeywords; which++)
      if (strcmp(params[0],keywords[which]) == 0) break;
    if (which == nkeywords) {
      char str[128];
      snprintf(str,128,"Keyword '%s' in 2nd MEAM paramter file not recognized",
               params[0]);
      error->all(FLERR,str);
    }
    for (i = 0; i < 3; i++) index[i] = 0;
    nindex = nparams - 2;
    for (i = 0; i < nindex; i++) index[i] = atoi(params[i+1]) - 1;

    // map lattce_meam value to an integer
    value = 0;
    if (which == 4) {
      if (strcmp(params[nparams-1],"fcc") == 0) value = FCC;
      else if (strcmp(params[nparams-1],"bcc") == 0) value = BCC;
      else if (strcmp(params[nparams-1],"hcp") == 0) value = HCP;
      else if (strcmp(params[nparams-1],"dim") == 0) value = DIM;
      else if (strcmp(params[nparams-1],"dia") == 0) value = DIA;
      else if (strcmp(params[nparams-1],"b1")  == 0) value = B1;
      else if (strcmp(params[nparams-1],"c11") == 0) value = C11;
      else if (strcmp(params[nparams-1],"l12") == 0) value = L12;
      else if (strcmp(params[nparams-1],"b2")  == 0) value = B2;
      else if (strcmp(params[nparams-1],"ch4")  == 0) value = CH4; 
      else if (strcmp(params[nparams-1],"lin")  == 0) value = LIN; 
      else if (strcmp(params[nparams-1],"zig")  == 0) value = ZIG; 
      else if (strcmp(params[nparams-1],"tri")  == 0) value = TRI; 
      else {
        char str[512];
        snprintf(str,512,"Unrecognized lattice type '%s' in 2nd MEAM parameter file. " 
                         "Available types are fcc, bcc, hcp, dim, dia, b1, c11, l12, b2, ch4, lin, zig, and tri.",
                         params[nparams-1]);
        error->all(FLERR,str);
      }
    } else {
      value = atof(params[nparams-1]);
    }
    // pass single setting to MEAM package

    int errorflag = 0;
    meambo_inst->meam_setup_param(which,value,nindex,index,&errorflag);
    if (errorflag) {
      char str[128];
      snprintf(str,128,"MEAM library error %d",errorflag);
      error->all(FLERR,str);
    }
  }

  delete [] params;


  // =====================================================================
  // begin section to read BCfile
  // =====================================================================

  if (strcmp(BCfile,"NULL") == 0){
    if (comm->me == 0) {
      error->warning(FLERR,"No MEAMBO BCfile is provided. Will run regular MEAM instead");
    }
  }else {
    if (comm->me == 0) {
      fp = force->open_potential(BCfile);
      if (fp == NULL) {
        fclose(fp);
        error->one(FLERR, "Cannot open MEAMBO BCfile. To run original MEAM, set NULL");
      }
    }

    if (meambo_inst->ntypBC_meambo > 0){
      meambo_inst->isBC = true;
    } else {
      meambo_inst->isBC = false;
      error->warning(FLERR,"ntypBC (number of bond types) parameter is not found or set to 0. Run MEAM instead");
    }

    // assume max number of bond types are the square of the number of elements 
    int len_arrays = nelements*nelements;
    int *i1elt = new int[len_arrays];
    int *i2elt = new int[len_arrays];    
    int *nb = new int[len_arrays];
    int *z0s = new int[len_arrays];
    int *z1s = new int[len_arrays];    

    double *betas0 = new double[len_arrays];
    double *betas1 = new double[len_arrays];
    double *betas2 = new double[len_arrays];
    double *betas3 = new double[len_arrays];

    double *p0 = new double[len_arrays];
    double *p1 = new double[len_arrays];
    double *p2 = new double[len_arrays];
    double *p3 = new double[len_arrays];

    double *e3_0 = new double[len_arrays];
    double *e3_1 = new double[len_arrays];
    double *e3_2 = new double[len_arrays];

    double *e2a_0_0 = new double[len_arrays];
    double *e2a_0_1 = new double[len_arrays];
    double *e2a_0_2 = new double[len_arrays];

    double *e2a_1_0 = new double[len_arrays];
    double *e2a_1_1 = new double[len_arrays];
    double *e2a_1_2 = new double[len_arrays];

    double *e2a_2_0 = new double[len_arrays];
    double *e2a_2_1 = new double[len_arrays];
    double *e2a_2_2 = new double[len_arrays];

    double *re2a_0 = new double[len_arrays];
    double *re2a_1 = new double[len_arrays];
    double *re2a_2 = new double[len_arrays];

    double *re3 = new double[len_arrays];

    // initialize
    for (int i = 0; i < len_arrays; i++) {
      i1elt[i] = 0;
      i2elt[i] = 0;
      nb[i] = 0;
      z0s[i] = 0;
      z1s[i] = 0;
      betas0[i] = betas1[i] = betas2[i] = betas3[i] = 0.0;
      p0[i] = p1[i] = p2[i] = p3[i] = 0.0;
      e3_0[i] = e3_1[i] = e3_2[i] = 0.0;
      e2a_0_0[i] = e2a_0_1[i] = e2a_0_2[i] = 0.0;
      e2a_1_0[i] = e2a_1_1[i] = e2a_1_2[i] = 0.0;
      e2a_2_0[i] = e2a_2_1[i] = e2a_2_2[i] = 0.0;
      re2a_0[i] = re2a_1[i] = re2a_2[i] = 0.0;
      re3[i] = 0.0;
    }

    nset = 0;
    int bond = 0;

    eof = 0;
    params_per_line = 13;  // the initial line for a single entry

    // BCfile structure is:
    // # DATE: YYYY-MM-DD meambc file for lammps package meambo
    // # elt_i1 elt_i2 bond Z(0) Z(1) beta0 beta1 beta2 beta3 p1 p2 p3 p3
    // # double bond: e2a_0(1:3), e2a_1(1:3), e2a_2(1:3), re2a(1:3)
    // # triple bond: e3(1:3), re3

    // cycle through all lines
    while (1) {

      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        }
        else {
          n = strlen(line) + 1;
        }
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);

      // strip comment, skip line if blank

      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
      if (nwords == 0) continue;      
      
      // tokenize elt & bond
      strcpy(copy,line);
      nwords = 0;
      words[nwords++] = strtok(copy,"' \t\n\r\f");
      while ((words[nwords++] = strtok(NULL,"' \t\n\r\f"))) continue;

      bond = atoi(words[2]);
      if (bond == 2) params_per_line = 25; // 13 + 4*3
      else if (bond == 3) params_per_line = 17; // 13 + 4
      else error->all(FLERR,"Incorrect bond order in MEAMBO BCfile, which must be 2 or 3.");

      while (nwords < params_per_line) {
        n = strlen(line);
        if (comm->me == 0) {
          ptr = fgets(&line[n],MAXLINE-n,fp);
          if (ptr == NULL) {
            eof = 1;
            fclose(fp);
          }
          else {
            n = strlen(line) + 1;
          }
        }
        MPI_Bcast(&eof,1,MPI_INT,0,world);
        if (eof) break;
        MPI_Bcast(&n,1,MPI_INT,0,world);
        MPI_Bcast(line,n,MPI_CHAR,0,world);
        if ((ptr = strchr(line,'#'))) *ptr = '\0';

        nwords = atom->count_words(line);
      }

      if (nwords != params_per_line)
        error->all(FLERR,"Incorrect format in MEAMBO BCfile. "
                         "MEAMBO BCfile should have the following form.\n"
                         "  elt_i1 elt_i2 bond Z(0) Z(1) beta0 beta1 beta2 beta3 p1 p2 p3 p3\n"
                         "** Additional parameters for double bond, i.e. (bond = 2)\n"
                         "  e2a_0(1:3)\n"
                         "  e2a_1(1:3)\n"
                         "  e2a_2(1:3)\n"
                         "  re2a(1:3)\n"
                         "** Additional parameters for triple bond, i.e. (bond = 3)\n"
                         "  e3(1:3), re3\n");

      nwords = 0;
      words[nwords++] = strtok(line,"' \t\n\r\f");
      while ((words[nwords++] = strtok(NULL,"' \t\n\r\f"))) continue;

      // search element i for bond order calculation from the element list
      int ind;
      for (ind = 0; ind < nelements; ind++) { 
        if (strcmp(words[0],elements[ind]) == 0) break;
      }
      if (ind >= nelements) {
        char str[128];
        sprintf(str,"Cannot find element i1 '%s' for bond in the element list.", words[0]);
        error->one(FLERR,str);
      }
      i1elt[nset]  = ind;
      
      // search element j for bond order calculation from the element list
      for (ind = 0; ind < nelements; ind++) { 
        if (strcmp(words[1],elements[ind]) == 0) break;
      }
      if (ind >= nelements) {
        char str[128];
        sprintf(str,"Cannot find element i2 '%s' for bond in the element list.", words[1]);
        error->one(FLERR,str);
      }      
      i2elt[nset]  = ind;
      nb[nset]     = atoi(words[2]);
      z0s[nset]     = atoi(words[3]); 
      z1s[nset]     = atoi(words[4]); 
      betas0[nset] = atof(words[5]);
      betas1[nset] = atof(words[6]);
      betas2[nset] = atof(words[7]);
      betas3[nset] = atof(words[8]);
      p0[nset]    = atof(words[9]);
      p1[nset]    = atof(words[10]);
      p2[nset]    = atof(words[11]);
      p3[nset]    = atof(words[12]);

      if (bond == 2) { // # double bond: e2a_0, e2a_1, e2a_2, re2a
        e2a_0_0[nset] = atof(words[13]);
        e2a_0_1[nset] = atof(words[14]);
        e2a_0_2[nset] = atof(words[15]);
        e2a_1_0[nset] = atof(words[16]);
        e2a_1_1[nset] = atof(words[17]);
        e2a_1_2[nset] = atof(words[18]);
        e2a_2_0[nset] = atof(words[19]);
        e2a_2_1[nset] = atof(words[20]);
        e2a_2_2[nset] = atof(words[21]);
        re2a_0[nset]    = atof(words[22]);
        re2a_1[nset]    = atof(words[23]);
        re2a_2[nset]    = atof(words[24]);
      }
      else if (bond == 3) { // # triple bond: e3, re3
        e3_0[nset] = atof(words[13]);
        e3_1[nset] = atof(words[14]);
        e3_2[nset] = atof(words[15]);
        re3[nset]   = atof(words[16]);
      }
      nset++;
    }

    // store parameters
    meambo_inst->meambo_setup_bond_params(nelements,nset,len_arrays,
        i1elt,i2elt,nb,z0s,z1s,betas0,betas1,betas2,betas3,
        p0,p1,p2,p3, e2a_0_0, e2a_0_1, e2a_0_2,
        e2a_1_0, e2a_1_1, e2a_1_2, e2a_2_0, e2a_2_1, e2a_2_2,
        re2a_0, re2a_1, re2a_2, e3_0, e3_1, e3_2, re3);

    // clean-up memory

    delete [] i1elt;
    delete [] i2elt;    
    delete [] nb;
    delete [] z0s;
    delete [] z1s;    
    delete [] betas0;
    delete [] betas1;
    delete [] betas2;
    delete [] betas3;
    delete [] p0;
    delete [] p1;
    delete [] p2;
    delete [] p3;
    delete [] e3_0;
    delete [] e3_1;
    delete [] e3_2;
    delete [] e2a_0_0;
    delete [] e2a_0_1;
    delete [] e2a_0_2;
    delete [] e2a_1_0;
    delete [] e2a_1_1;
    delete [] e2a_1_2;
    delete [] e2a_2_0;
    delete [] e2a_2_1;
    delete [] e2a_2_2;
    delete [] re2a_0;
    delete [] re2a_1;
    delete [] re2a_2;
    delete [] re3;

    // =====================================================================
    // end   section to read BCfile
    // =====================================================================
  }

  delete [] words;

}

/* ---------------------------------------------------------------------- */

int PairMEAMBO::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = meambo_inst->rho0[j];
    buf[m++] = meambo_inst->rho1[j];
    buf[m++] = meambo_inst->rho2[j];
    buf[m++] = meambo_inst->rho3[j];
    buf[m++] = meambo_inst->frhop[j];
    buf[m++] = meambo_inst->gamma[j];
    buf[m++] = meambo_inst->dgamma1[j];
    buf[m++] = meambo_inst->dgamma2[j];
    buf[m++] = meambo_inst->dgamma3[j];
    buf[m++] = meambo_inst->arho2b[j];
    buf[m++] = meambo_inst->arho1[j][0];
    buf[m++] = meambo_inst->arho1[j][1];
    buf[m++] = meambo_inst->arho1[j][2];
    buf[m++] = meambo_inst->arho2[j][0];
    buf[m++] = meambo_inst->arho2[j][1];
    buf[m++] = meambo_inst->arho2[j][2];
    buf[m++] = meambo_inst->arho2[j][3];
    buf[m++] = meambo_inst->arho2[j][4];
    buf[m++] = meambo_inst->arho2[j][5];
    for (k = 0; k < 10; k++) buf[m++] = meambo_inst->arho3[j][k];
    buf[m++] = meambo_inst->arho3b[j][0];
    buf[m++] = meambo_inst->arho3b[j][1];
    buf[m++] = meambo_inst->arho3b[j][2];
    buf[m++] = meambo_inst->t_ave[j][0];
    buf[m++] = meambo_inst->t_ave[j][1];
    buf[m++] = meambo_inst->t_ave[j][2];
    buf[m++] = meambo_inst->tsq_ave[j][0];
    buf[m++] = meambo_inst->tsq_ave[j][1];
    buf[m++] = meambo_inst->tsq_ave[j][2];

    if (meambo_inst->isBC) {  //meambo
      buf[m++] = meambo_inst->Z0BC[j];
      buf[m++] = meambo_inst->Z1BC[j];

      buf[m++] = meambo_inst->delta2_Z0[j];
      buf[m++] = meambo_inst->ddelta2_Z0[j];
      buf[m++] = meambo_inst->delta2_Z1[j];
      buf[m++] = meambo_inst->ddelta2_Z1[j];

      buf[m++] = meambo_inst->arho1BC[j][0];
      buf[m++] = meambo_inst->arho1BC[j][1];
      buf[m++] = meambo_inst->arho1BC[j][2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairMEAMBO::unpack_forward_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    meambo_inst->rho0[i] = buf[m++];
    meambo_inst->rho1[i] = buf[m++];
    meambo_inst->rho2[i] = buf[m++];
    meambo_inst->rho3[i] = buf[m++];
    meambo_inst->frhop[i] = buf[m++];
    meambo_inst->gamma[i] = buf[m++];
    meambo_inst->dgamma1[i] = buf[m++];
    meambo_inst->dgamma2[i] = buf[m++];
    meambo_inst->dgamma3[i] = buf[m++];
    meambo_inst->arho2b[i] = buf[m++];
    meambo_inst->arho1[i][0] = buf[m++];
    meambo_inst->arho1[i][1] = buf[m++];
    meambo_inst->arho1[i][2] = buf[m++];
    meambo_inst->arho2[i][0] = buf[m++];
    meambo_inst->arho2[i][1] = buf[m++];
    meambo_inst->arho2[i][2] = buf[m++];
    meambo_inst->arho2[i][3] = buf[m++];
    meambo_inst->arho2[i][4] = buf[m++];
    meambo_inst->arho2[i][5] = buf[m++];
    for (k = 0; k < 10; k++) meambo_inst->arho3[i][k] = buf[m++];
    meambo_inst->arho3b[i][0] = buf[m++];
    meambo_inst->arho3b[i][1] = buf[m++];
    meambo_inst->arho3b[i][2] = buf[m++];
    meambo_inst->t_ave[i][0] = buf[m++];
    meambo_inst->t_ave[i][1] = buf[m++];
    meambo_inst->t_ave[i][2] = buf[m++];
    meambo_inst->tsq_ave[i][0] = buf[m++];
    meambo_inst->tsq_ave[i][1] = buf[m++];
    meambo_inst->tsq_ave[i][2] = buf[m++];

    if (meambo_inst->isBC) {  //meambo
      meambo_inst->Z0BC[i] = buf[m++];
      meambo_inst->Z1BC[i] = buf[m++];

      meambo_inst->delta2_Z0[i] = buf[m++];
      meambo_inst->ddelta2_Z0[i] = buf[m++];
      meambo_inst->delta2_Z1[i] = buf[m++];
      meambo_inst->ddelta2_Z1[i] = buf[m++];

      meambo_inst->arho1BC[i][0] = buf[m++];
      meambo_inst->arho1BC[i][1] = buf[m++];
      meambo_inst->arho1BC[i][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairMEAMBO::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = meambo_inst->rho0[i];
    buf[m++] = meambo_inst->arho2b[i];
    buf[m++] = meambo_inst->arho1[i][0];
    buf[m++] = meambo_inst->arho1[i][1];
    buf[m++] = meambo_inst->arho1[i][2];
    buf[m++] = meambo_inst->arho2[i][0];
    buf[m++] = meambo_inst->arho2[i][1];
    buf[m++] = meambo_inst->arho2[i][2];
    buf[m++] = meambo_inst->arho2[i][3];
    buf[m++] = meambo_inst->arho2[i][4];
    buf[m++] = meambo_inst->arho2[i][5];
    for (k = 0; k < 10; k++) buf[m++] = meambo_inst->arho3[i][k];
    buf[m++] = meambo_inst->arho3b[i][0];
    buf[m++] = meambo_inst->arho3b[i][1];
    buf[m++] = meambo_inst->arho3b[i][2];
    buf[m++] = meambo_inst->t_ave[i][0];
    buf[m++] = meambo_inst->t_ave[i][1];
    buf[m++] = meambo_inst->t_ave[i][2];
    buf[m++] = meambo_inst->tsq_ave[i][0];
    buf[m++] = meambo_inst->tsq_ave[i][1];
    buf[m++] = meambo_inst->tsq_ave[i][2];

    if (meambo_inst->isBC) { //meambo
      buf[m++] = meambo_inst->Z0BC[i];
      buf[m++] = meambo_inst->arho1BC[i][0];
      buf[m++] = meambo_inst->arho1BC[i][1];
      buf[m++] = meambo_inst->arho1BC[i][2];
    }

  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMEAMBO::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    meambo_inst->rho0[j] += buf[m++];
    meambo_inst->arho2b[j] += buf[m++];
    meambo_inst->arho1[j][0] += buf[m++];
    meambo_inst->arho1[j][1] += buf[m++];
    meambo_inst->arho1[j][2] += buf[m++];
    meambo_inst->arho2[j][0] += buf[m++];
    meambo_inst->arho2[j][1] += buf[m++];
    meambo_inst->arho2[j][2] += buf[m++];
    meambo_inst->arho2[j][3] += buf[m++];
    meambo_inst->arho2[j][4] += buf[m++];
    meambo_inst->arho2[j][5] += buf[m++];
    for (k = 0; k < 10; k++) meambo_inst->arho3[j][k] += buf[m++];
    meambo_inst->arho3b[j][0] += buf[m++];
    meambo_inst->arho3b[j][1] += buf[m++];
    meambo_inst->arho3b[j][2] += buf[m++];
    meambo_inst->t_ave[j][0] += buf[m++];
    meambo_inst->t_ave[j][1] += buf[m++];
    meambo_inst->t_ave[j][2] += buf[m++];
    meambo_inst->tsq_ave[j][0] += buf[m++];
    meambo_inst->tsq_ave[j][1] += buf[m++];
    meambo_inst->tsq_ave[j][2] += buf[m++];

    if (meambo_inst->isBC) {      //meambo
      meambo_inst->Z0BC[j] += buf[m++];
      meambo_inst->arho1BC[j][0] += buf[m++];
      meambo_inst->arho1BC[j][1] += buf[m++];
      meambo_inst->arho1BC[j][2] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairMEAMBO::memory_usage()
{
  double bytes = 11 * meambo_inst->nmax * sizeof(double);
  // bytes += (3 + 6 + 10 + 3 + 3 + 3) * meambo_inst->nmax * sizeof(double);
  // bytes += 3 * meambo_inst->maxneigh * sizeof(double);
  bytes += (3 + 6 + 10 + 3 + 3 + 3 + 1 + 6 + 100) * meambo_inst->nmax * sizeof(double);
  bytes += 10 * meambo_inst->maxneigh * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   strip special bond flags from neighbor list entries
   are not used with MEAM
   need to do here so Fortran lib doesn't see them
   done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
------------------------------------------------------------------------- */

void PairMEAMBO::neigh_strip(int inum, int *ilist,
                           int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j] &= NEIGHMASK;
  }
}
