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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_buck_coul_long.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "memory.h"
#include "neighbor.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairBuckCoulLong::~PairBuckCoulLong()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut_lj);
    memory->destroy_2d_double_array(cut_ljsq);
    memory->destroy_2d_double_array(a);
    memory->destroy_2d_double_array(rho);
    memory->destroy_2d_double_array(c);
    memory->destroy_2d_double_array(rhoinv);
    memory->destroy_2d_double_array(buck1);
    memory->destroy_2d_double_array(buck2);
    memory->destroy_2d_double_array(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulLong::compute(int eflag, int vflag)
{
  int i,j,k,numneigh,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r2inv,r6inv,forcecoul,forcebuck,fforce,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double factor,phicoul,phibuck,r,rexp;
  int *neighs;
  double **f;

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq) {
	  r = sqrt(rsq);
	  grij = g_ewald * r;
	  expm2 = exp(-grij*grij);
	  t = 1.0 / (1.0 + EWALD_P*grij);
	  erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	  prefactor = qqrd2e * qtmp*q[j]/r;
	  forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	  if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
	} else forcecoul = 0.0;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
          r = sqrt(rsq);
	  rexp = exp(-r*rhoinv[itype][jtype]);
	  forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
	} else forcebuck = 0.0;

	fforce = (forcecoul + factor_lj*forcebuck) * r2inv;

	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}

	if (eflag) {
	  if (newton_pair || j < nlocal) factor = 1.0;
	  else factor = 0.5;
	  if (rsq < cut_coulsq) {
	    phicoul = prefactor*erfc;
	    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
	    eng_coul += factor*phicoul;
	  }
	  if (rsq < cut_ljsq[itype][jtype]) {
	    phibuck = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
	      offset[itype][jtype];
	    eng_vdwl += factor*factor_lj*phibuck;
	  }
	}

	if (vflag == 1) {
	  if (newton_pair || j < nlocal) {
	    virial[0] += delx*delx*fforce;
	    virial[1] += dely*dely*fforce;
	    virial[2] += delz*delz*fforce;
	    virial[3] += delx*dely*fforce;
	    virial[4] += delx*delz*fforce;
	    virial[5] += dely*delz*fforce;
	  } else {
	    virial[0] += 0.5*delx*delx*fforce;
	    virial[1] += 0.5*dely*dely*fforce;
	    virial[2] += 0.5*delz*delz*fforce;
	    virial[3] += 0.5*delx*dely*fforce;
	    virial[4] += 0.5*delx*delz*fforce;
	    virial[5] += 0.5*dely*delz*fforce;
	  }
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBuckCoulLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  cut_lj = memory->create_2d_double_array(n+1,n+1,"pair:cut_lj");
  cut_ljsq = memory->create_2d_double_array(n+1,n+1,"pair:cut_ljsq");
  a = memory->create_2d_double_array(n+1,n+1,"pair:a");
  rho = memory->create_2d_double_array(n+1,n+1,"pair:rho");
  c = memory->create_2d_double_array(n+1,n+1,"pair:c");
  rhoinv = memory->create_2d_double_array(n+1,n+1,"pair:rhoinv");
  buck1 = memory->create_2d_double_array(n+1,n+1,"pair:buck1");
  buck2 = memory->create_2d_double_array(n+1,n+1,"pair:buck2");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBuckCoulLong::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all("Illegal pair_style command");

  cut_lj_global = atof(arg[0]);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = atof(arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuckCoulLong::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a_one = atof(arg[2]);
  double rho_one = atof(arg[3]);
  double c_one = atof(arg[4]);

  double cut_lj_one = cut_lj_global;
  if (narg == 6) cut_lj_one = atof(arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuckCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  rhoinv[i][j] = 1.0/rho[i][j];
  buck1[i][j] = a[i][j]/rho[i][j];
  buck2[i][j] = 6.0*c[i][j];
     
  if (offset_flag) {
    double rexp = exp(-cut_lj[i][j]/rho[i][j]);
    offset[i][j] = a[i][j]*rexp - c[i][j]/pow(cut_lj[i][j],6.0);
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBuckCoulLong::init_style()
{
  // require an atom style with charge defined

  if (atom->charge_allow == 0)
    error->all("Must use charged atom style with this pair style");

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

 if (force->kspace == NULL) 
    error->all("Pair style is incompatible with KSpace style");
  else if (strcmp(force->kspace_style,"ewald") == 0)
    g_ewald = force->kspace->g_ewald;
  else if (strcmp(force->kspace_style,"pppm") == 0)
    g_ewald = force->kspace->g_ewald;
  else error->all("Pair style is incompatible with KSpace style");
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&a[i][j],sizeof(double),1,fp);
	fwrite(&rho[i][j],sizeof(double),1,fp);
	fwrite(&c[i][j],sizeof(double),1,fp);
	fwrite(&cut_lj[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&a[i][j],sizeof(double),1,fp);
	  fread(&rho[i][j],sizeof(double),1,fp);
	  fread(&c[i][j],sizeof(double),1,fp);
	  fread(&cut_lj[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulLong::single(int i, int j, int itype, int jtype,
			     double rsq, double factor_coul, double factor_lj,
			     int eflag, One &one)
{
  double r2inv,r6inv,r,rexp,grij,expm2,t,erfc,prefactor;
  double forcecoul,forcebuck,phicoul,phibuck;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    t = 1.0 / (1.0 + EWALD_P*grij);
    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    r = sqrt(rsq);
    rexp = exp(-r*rhoinv[itype][jtype]);
    forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
  } else forcebuck = 0.0;
  one.fforce = (forcecoul + factor_lj*forcebuck) * r2inv;
  
  if (eflag) {
    if (rsq < cut_coulsq) {
      phicoul = prefactor*erfc;
      if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
      one.eng_coul = phicoul;
    } else one.eng_coul = 0.0;
    if (rsq < cut_ljsq[itype][jtype]) {
      phibuck = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
	offset[itype][jtype];
      one.eng_vdwl = factor_lj*phibuck;
    } else one.eng_vdwl = 0.0;
  }
}
