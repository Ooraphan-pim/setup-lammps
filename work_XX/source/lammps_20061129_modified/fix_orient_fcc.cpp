/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Koenraad Janssens and David Olmsted (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "fix_orient_fcc.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "neighbor.h"
#include "comm.h"
#include "output.h"
#include "error.h"

#define BIG 1000000000

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixOrientFCC::FixOrientFCC(int narg, char **arg) : Fix(narg, arg)
{
  MPI_Comm_rank(world,&me);

  if (narg != 11) error->all("Illegal fix orient/fcc command");

  neigh_full_every = 1;

  nstats = atoi(arg[3]);
  direction_of_motion = atoi(arg[4]);
  a = atof(arg[5]);
  Vxi = atof(arg[6]);
  uxif_low = atof(arg[7]);
  uxif_high = atof(arg[8]);

  if (direction_of_motion == 0) {
    int n = strlen(arg[9]) + 1;
    chifilename = new char[n];
    strcpy(chifilename,arg[9]);
    n = strlen(arg[10]) + 1;
    xifilename = new char[n];
    strcpy(xifilename,arg[10]);
  } else if (direction_of_motion == 1) {
    int n = strlen(arg[9]) + 1;
    xifilename = new char[n];
    strcpy(xifilename,arg[9]);
    n = strlen(arg[10]) + 1;
    chifilename = new char[n];
    strcpy(chifilename,arg[10]);
  } else error->all("Illegal fix orient/fcc command");

  // initializations

  PI = 4.0*atan(1.0);
  half_fcc_nn = 6;
  use_xismooth = false;
  double xicutoff = 1.57;
  xicutoffsq = xicutoff * xicutoff;
  cutsq = 0.5 * a*a*xicutoffsq;
  nmax = 0;

  // read xi and chi reference orientations from files

  if (me == 0) {
    char line[512];
    char *result;
    int count;

    FILE *infile = fopen(xifilename,"r");
    if (infile == NULL) error->one("Fix orient/fcc file open failed");
    for (int i = 0; i < 6; i++) {
      result = fgets(line,512,infile);
      if (!result) error->one("Fix orient/fcc file read failed");
      count = sscanf(line,"%lg %lg %lg",&Rxi[i][0],&Rxi[i][1],&Rxi[i][2]);
      if (count != 3) error->one("Fix orient/fcc file read failed");
    }
    fclose(infile);

    infile = fopen(chifilename,"r");
    if (infile == NULL) error->one("Fix orient/fcc file open failed");
    for (int i = 0; i < 6; i++) {
      result = fgets(line,512,infile);
      if (!result) error->one("Fix orient/fcc file read failed");
      count = sscanf(line,"%lg %lg %lg",&Rchi[i][0],&Rchi[i][1],&Rchi[i][2]);
      if (count != 3) error->one("Fix orient/fcc file read failed");
    }
    fclose(infile);
  }

  MPI_Bcast(&Rxi[0][0],18,MPI_DOUBLE,0,world);
  MPI_Bcast(&Rchi[0][0],18,MPI_DOUBLE,0,world);

  // make copy of the reference vectors

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 3; j++) {
      half_xi_chi_vec[0][i][j] = Rxi[i][j];
      half_xi_chi_vec[1][i][j] = Rchi[i][j];
    }

  // compute xiid,xi0,xi1 for all 12 neighbors
  // xi is the favored crystal
  // want order parameter when actual is Rchi

  double xi_sq,dxi[3],rchi[3];

  xiid = 0.0;
  for (int i = 0; i < 6; i++) {
    rchi[0] = Rchi[i][0];
    rchi[1] = Rchi[i][1];
    rchi[2] = Rchi[i][2];
    find_best_ref(rchi,0,xi_sq,dxi);
    xiid += sqrt(xi_sq);
    for (int j = 0; j < 3; j++) rchi[j] = -rchi[j];
    find_best_ref(rchi,0,xi_sq,dxi);
    xiid += sqrt(xi_sq);
  }

  xiid /= 12.0;
  xi0 = uxif_low * xiid;
  xi1 = uxif_high * xiid;
}

/* ---------------------------------------------------------------------- */

FixOrientFCC::~FixOrientFCC()
{
  delete [] xifilename;
  delete [] chifilename;
  if (nmax) delete [] nbr;
}

/* ---------------------------------------------------------------------- */

int FixOrientFCC::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= THERMO;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::init()
{
  if (use_xismooth) comm->maxforward_fix = MAX(comm->maxforward_fix,62);
  else comm->maxforward_fix = MAX(comm->maxforward_fix,50);

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if (thermo_print || thermo_energy) thermo_flag = 1;
  else thermo_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::setup()
{
  eflag_on = 1;
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  eflag_on = 0;
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::post_force(int vflag)
{
  int i,j,k,m,n,nn,nsort,id_self;
  int *neighs;
  double edelta,added_energy,omega;
  double dx,dy,dz,rsq,xismooth,xi_sq,duxi,duxi_other;
  double dxi[3];
  double *dxiptr;
  bool found_myself;

  int eflag = false;
  if (thermo_flag) {
    if (eflag_on) eflag = true;
    else if (output->next_thermo == update->ntimestep) eflag = true;
  }

  // set local ptrs

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  // insure nbr data structure is adequate size

  if (nall > nmax) {
    if (nmax) delete [] nbr;
    nbr = new Nbr[nall];
    nmax = nall;
  }

  // loop over owned atoms and build Nbr data structure of neighbors
  // use full neighbor list

  if (eflag) added_energy = 0.0;
  int count = 0;
  int mincount = BIG;
  int maxcount = 0;

  for (i = 0; i < nlocal; i++) {
    neighs = neighbor->firstneigh_full[i];
    n = neighbor->numneigh_full[i];

    if (n < mincount) mincount = n;
    if (n > maxcount) {
      if (maxcount) delete [] sort;
      sort = new Sort[n];
      maxcount = n;
    }

    // loop over all neighbors of atom i
    // for those within cutsq, build sort data structure
    // store local id, rsq, delta vector, xismooth (if included)

    nsort = 0;
    for (k = 0; k < n; k++) {
      j = neighs[k];
      count++;

      dx = x[i][0] - x[j][0];
      dy = x[i][1] - x[j][1];
      dz = x[i][2] - x[j][2];
      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < cutsq) {
	sort[nsort].id = j;
	sort[nsort].rsq = rsq;
	sort[nsort].delta[0] = dx;
	sort[nsort].delta[1] = dy;
	sort[nsort].delta[2] = dz;
	if (use_xismooth) {
	  xismooth = (xicutoffsq - 2.0*rsq/(a*a)) / (xicutoffsq - 1.0);
	  sort[nsort].xismooth = 1.0 - fabs(1.0-xismooth);
	}
	nsort++;
      }
    }

    // sort neighbors by rsq distance
    // no need to sort if nsort <= 12
    
    if (nsort > 12) qsort(sort,nsort,sizeof(Sort),compare);

    // copy up to 12 nearest neighbors into nbr data structure
    // operate on delta vector via find_best_ref() to compute dxi
    
    n = MIN(12,nsort);
    nbr[i].n = n;
    if (n == 0) continue;

    double xi_total = 0.0;
    for (j = 0; j < n; j++) {
      find_best_ref(sort[j].delta,0,xi_sq,dxi);
      xi_total += sqrt(xi_sq);
      nbr[i].id[j] = sort[j].id;
      nbr[i].dxi[j][0] = dxi[0]/n;
      nbr[i].dxi[j][1] = dxi[1]/n;
      nbr[i].dxi[j][2] = dxi[2]/n;
      if (use_xismooth) nbr[i].xismooth[j] = sort[j].xismooth;
    }
    xi_total /= n;

    // compute potential derivative to xi

    if (xi_total < xi0) {
      nbr[i].duxi = 0.0;
      if (eflag) edelta = 0.0;
    } else if (xi_total > xi1) {
      nbr[i].duxi = 0.0;
      if (eflag) edelta = Vxi;
    } else {
      omega = (0.5*PI)*(xi_total-xi0) / (xi1-xi0);
      nbr[i].duxi = PI*Vxi*sin(2.0*omega) / (2.0*(xi1-xi0));
      if (eflag) edelta = Vxi*(1 - cos(2.0*omega)) / 2.0;
    }
    if (eflag) added_energy += edelta;
  }
  
  if (maxcount) delete [] sort;
  
  // communicate to acquire nbr data for ghost atoms

  comm->comm_fix(this);

  // compute grain boundary force on each owned atom
  // skip atoms not in group

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    n = nbr[i].n;
    duxi = nbr[i].duxi;

    for (j = 0; j < n; j++) {
      dxiptr = &nbr[i].dxi[j][0];
      if (use_xismooth) {
	xismooth = nbr[i].xismooth[j];
        f[i][0] += duxi * dxiptr[0] * xismooth;
        f[i][1] += duxi * dxiptr[1] * xismooth;
        f[i][2] += duxi * dxiptr[2] * xismooth;
      } else {
	f[i][0] += duxi * dxiptr[0];
        f[i][1] += duxi * dxiptr[1];
        f[i][2] += duxi * dxiptr[2];
      }

      // m = local index of neighbor
      // id_self = ID for atom I in atom M's neighbor list
      // if M is local atom, id_self will be local ID of atom I
      // if M is ghost atom, id_self will be global ID of atom I

      m = nbr[i].id[j]; 
      if (m < nlocal) id_self = i;
      else id_self = tag[i];
      found_myself = false;
      nn = nbr[m].n;

      for (k = 0; k < nn; k++) {
	if (id_self == nbr[m].id[k]) {
	  if (found_myself) error->one("Fix orient/fcc found self twice");
	  found_myself = true;
	  duxi_other = nbr[m].duxi;
	  dxiptr = &nbr[m].dxi[k][0];
	  if (use_xismooth) {
	    xismooth = nbr[m].xismooth[k];
	    f[i][0] -= duxi_other * dxiptr[0] * xismooth;
	    f[i][1] -= duxi_other * dxiptr[1] * xismooth;
	    f[i][2] -= duxi_other * dxiptr[2] * xismooth;
	  } else {
	    f[i][0] -= duxi_other * dxiptr[0];
	    f[i][1] -= duxi_other * dxiptr[1];
	    f[i][2] -= duxi_other * dxiptr[2];
	  }
	}
      }
    }
  }

  // sum energy across procs

  if (eflag)
    MPI_Allreduce(&added_energy,&total_added_e,1,MPI_DOUBLE,MPI_SUM,world);

  // print statistics every nstats timesteps

  if (nstats && update->ntimestep % nstats == 0) {
    int total;
    MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,world);
    double ave = total/atom->natoms;

    int min,max;
    MPI_Allreduce(&mincount,&min,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&maxcount,&max,1,MPI_INT,MPI_MAX,world);

    if (me == 0) { 
      if (screen) {
	fprintf(screen,"orient step %d: %g atoms have %d neighbors\n",
		update->ntimestep,atom->natoms,total);
	fprintf(screen,"  neighs: min = %d, max = %d, ave = %g\n",
		min,max,ave);
      }
      if (logfile) {
	fprintf(logfile,"orient step %d: %g atoms have %d neighbors\n",
		update->ntimestep,atom->natoms,total);
	fprintf(logfile,"  neighs: min = %d, max = %d, ave = %g\n",
		min,max,ave);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixOrientFCC::thermo_fields(int n, int *flags, char **keywords)
{
  if (n == 0) return 1;
  flags[0] = 3;
  strcpy(keywords[0],"Orient");
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixOrientFCC::thermo_compute(double *values)
{
  values[0] = total_added_e;
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixOrientFCC::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,k,id,num;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int m = 0;

  for (i = 0; i < n; i++) {
    k = list[i];
    num = nbr[k].n;
    buf[m++] = static_cast<double> (num);
    buf[m++] = nbr[k].duxi;

    for (j = 0; j < num; j++) {
      if (use_xismooth) buf[m++] = nbr[m].xismooth[j];
      buf[m++] = nbr[k].dxi[j][0];
      buf[m++] = nbr[k].dxi[j][1];
      buf[m++] = nbr[k].dxi[j][2];

      // id stored in buf needs to be global ID
      // if k is a local atom, it stores local IDs, so convert to global
      // if k is a ghost atom (already comm'd), its IDs are already global

      id = nbr[k].id[j];
      if (k < nlocal) id = tag[id];
      buf[m++] = static_cast<double> (id);
    }

    m += (12-num) * 3;
    if (use_xismooth) m += 12-num;
  }

  if (use_xismooth) return 62;
  return 50;
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::unpack_comm(int n, int first, double *buf)
{
  int i,j,num;
  int last = first + n;
  int m = 0;

  for (i = first; i < last; i++) {
    nbr[i].n = num = static_cast<int> (buf[m++]);
    nbr[i].duxi = buf[m++];

    for (j = 0; j < num; j++) {
      if (use_xismooth) nbr[i].xismooth[j] = buf[m++];
      nbr[i].dxi[j][0] = buf[m++];
      nbr[i].dxi[j][1] = buf[m++];
      nbr[i].dxi[j][2] = buf[m++];
      nbr[i].id[j] = static_cast<int> (buf[m++]);
    }

    m += (12-num) * 3;
    if (use_xismooth) m += 12-num;
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientFCC::find_best_ref(double *displs, int which_crystal, 
				 double &xi_sq, double *dxi)
{
  int i;
  double dot,tmp;

  double  best_dot  = -1.0;         // best is biggest (smallest angle)
  int     best_i    = -1;
  int     best_sign = 0;

  for (i = 0; i < half_fcc_nn; i++) {
    dot = displs[0] * half_xi_chi_vec[which_crystal][i][0] +
      displs[1] * half_xi_chi_vec[which_crystal][i][1] +
      displs[2] * half_xi_chi_vec[which_crystal][i][2];
    if (fabs(dot) > best_dot) {
      best_dot = fabs(dot);
      best_i = i;
      if (dot < 0.0) best_sign = -1;
      else best_sign = 1;
    }
  }

  xi_sq = 0.0;
  for (i = 0; i < 3; i++) {
    tmp = displs[i] - best_sign * half_xi_chi_vec[which_crystal][best_i][i];
    xi_sq += tmp*tmp;
  }

  if (xi_sq > 0.0) {
    double xi = sqrt(xi_sq);
    for (i = 0; i < 3; i++)
      dxi[i] = (best_sign * half_xi_chi_vec[which_crystal][best_i][i] - 
		displs[i]) / xi;
  } else dxi[0] = dxi[1] = dxi[2] = 0.0;
}

/* ----------------------------------------------------------------------
   compare two neighbors I and J in sort data structure
   called via qsort in post_force() method
   is a static method so can't access sort data structure directly
   return -1 if I < J, 0 if I = J, 1 if I > J
   do comparison based on rsq distance
------------------------------------------------------------------------- */

int FixOrientFCC::compare(const void *pi, const void *pj)
{
  FixOrientFCC::Sort *ineigh = (FixOrientFCC::Sort *) pi;
  FixOrientFCC::Sort *jneigh = (FixOrientFCC::Sort *) pj;

  if (ineigh->rsq < jneigh->rsq) return -1;
  else if (ineigh->rsq > jneigh->rsq) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

int FixOrientFCC::memory_usage()
{
  int bytes = atom->nmax * sizeof(Nbr);
  return bytes;
}
