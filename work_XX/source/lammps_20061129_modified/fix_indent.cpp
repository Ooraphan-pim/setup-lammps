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
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_indent.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

#define NONE     0
#define SPHERE   1
#define CYLINDER 2

/* ---------------------------------------------------------------------- */

FixIndent::FixIndent(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 4) error->all("Illegal fix indent command");
  k = atof(arg[3]);

  // set input line defaults

  istyle = NONE;
  vx = vy = vz = 0.0;
  scaleflag = 1;
  radflag = 0;
  r0_start = 0.0;

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of fix indent with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling to indenter force constant, geometry, and velocity

  k /= xscale;
  k3 = k/3.0;
  vx *= xscale;
  vy *= yscale;
  vz *= zscale;

  if (istyle == SPHERE) {
    x0 *= xscale;
    y0 *= yscale;
    z0 *= zscale;
    r0_stop *= xscale;
    r0_start *= xscale;
  } else if (istyle == CYLINDER) {
    if (cdim == 0) {
      c1 *= yscale;
      c2 *= zscale;
      r0_stop *= xscale;
      r0_start *= xscale;
    } else if (cdim == 1) {
      c1 *= xscale;
      c2 *= zscale;
      r0_stop *= yscale;
      r0_start *= yscale;
    } else if (cdim == 2) {
      c1 *= xscale;
      c2 *= yscale;
      r0_stop *= zscale;
      r0_start *= zscale;
    }
  } else error->all("Illegal fix indent command");

  // time 0 for the indenter movement

  ntimestep_initial = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixIndent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndent::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if (thermo_print || thermo_energy) thermo_flag = 1;
  else thermo_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixIndent::setup()
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

void FixIndent::min_setup()
{
  eflag_on = 1;
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force(int vflag)
{
  bool eflag = false;
  if (thermo_flag) {
    if (eflag_on) eflag = true;
    else if (output->next_thermo == update->ntimestep) eflag = true;
  }
  if (eflag) eng = 0.0;

  // set current r0
  // for minimization, always set to r0_stop

  double r0;
  if (!radflag || update->whichflag) r0 = r0_stop;
  else {
    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;
    r0 = r0_start + delta * (r0_stop-r0_start);
  }

  // spherical indenter

  if (istyle == SPHERE) {

    // x1,y1,z1 = current position of indenter from original x0,y0,z0

    double delta = (update->ntimestep - ntimestep_initial) * update->dt;
    double x1 = x0 + delta*vx;
    double y1 = y0 + delta*vy;
    double z1 = z0 + delta*vz;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	delx = x[i][0] - x1;
	dely = x[i][1] - y1;
	delz = x[i][2] - z1;
	r = sqrt(delx*delx + dely*dely + delz*delz);
	dr = r - r0;
	if (dr >= 0.0) continue;
	fmag = k*dr*dr;
	f[i][0] += delx*fmag/r;
	f[i][1] += dely*fmag/r;
	f[i][2] += delz*fmag/r;
	if (eflag) eng -= k3 * dr*dr*dr;
      }

  // cylindrical indenter

  } else {

    // c1new,c2new = current coords of indenter axis from original c1,c2
	      
    double delta = (update->ntimestep - ntimestep_initial) * update->dt;
    double c1new,c2new;
    if (cdim == 0) {
      c1new = c1 + delta*vy;
      c2new = c2 + delta*vz;
    } else if (cdim == 1) {
      c1new = c1 + delta*vx;
      c2new = c2 + delta*vz;
    } else {
      c1new = c1 + delta*vx;
      c2new = c2 + delta*vy;
    }
    
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    double delx,dely,delz,r,dr,fmag;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (cdim == 0) {
	  delx = 0;
	  dely = x[i][1] - c1new;
	  delz = x[i][2] - c2new;
	} else if (cdim == 1) {
	  delx = x[i][0] - c1new;
	  dely = 0;
	  delz = x[i][2] - c2new;
	} else {
	  delx = x[i][0] - c1new;
	  dely = x[i][1] - c2new;
	  delz = 0;
	}
	r = sqrt(delx*delx + dely*dely + delz*delz);
	dr = r - r0;
	if (dr >= 0.0) continue;
	fmag = k*dr*dr;
	f[i][0] += delx*fmag/r;
	f[i][1] += dely*fmag/r;
	f[i][2] += delz*fmag/r;
	if (eflag) eng -= k3 * dr*dr*dr;
      }
  }

  if (eflag) MPI_Allreduce(&eng,&etotal,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixIndent::thermo_fields(int n, int *flags, char **keywords)
{
  if (n == 0) return 1;
  flags[0] = 3;
  strcpy(keywords[0],"Indent");
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixIndent::thermo_compute(double *values)
{
  values[0] = etotal;
  return 1;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of fix indent input line 
------------------------------------------------------------------------- */

void FixIndent::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal fix indent command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      x0 = atof(arg[iarg+1]);
      y0 = atof(arg[iarg+2]);
      z0 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = SPHERE;
      iarg += 5;
    } else if (strcmp(arg[iarg],"cylinder") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"x") == 0) cdim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) cdim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) cdim = 2;
      else error->all("Illegal fix indent command");
      c1 = atof(arg[iarg+2]);
      c2 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = CYLINDER;
      iarg += 5;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix indent command");
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rstart") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      radflag = 1;
      r0_start = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal fix indent command");
      iarg += 2;
    } else error->all("Illegal fix indent command");
  }
}
