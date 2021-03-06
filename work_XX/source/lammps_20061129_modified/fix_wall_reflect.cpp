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

#include "string.h"
#include "fix_wall_reflect.h"
#include "atom.h"
#include "domain.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FixWallReflect::FixWallReflect(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 4) error->all("Illegal fix wall/reflect command");

  xloflag = xhiflag = yloflag = yhiflag = zloflag = zhiflag = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"xlo") == 0) xloflag = 1;
    else if (strcmp(arg[iarg],"xhi") == 0) xhiflag = 1;
    else if (strcmp(arg[iarg],"ylo") == 0) yloflag = 1;
    else if (strcmp(arg[iarg],"yhi") == 0) yhiflag = 1;
    else if (strcmp(arg[iarg],"zlo") == 0) zloflag = 1;
    else if (strcmp(arg[iarg],"zhi") == 0) zhiflag = 1;
    else error->all("Illegal fix wall/reflect command");
  }
}

/* ---------------------------------------------------------------------- */

int FixWallReflect::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::initial_integrate()
{
  double xlo = domain->boxxlo;
  double xhi = domain->boxxhi;
  double ylo = domain->boxylo;
  double yhi = domain->boxyhi;
  double zlo = domain->boxzlo;
  double zhi = domain->boxzhi;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (xloflag && x[i][0] < xlo) {
	x[i][0] = xlo + (xlo - x[i][0]);
	v[i][0] = -v[i][0];
      }
      if (xhiflag && x[i][0] > xhi) {
	x[i][0] = xhi - (x[i][0] - xhi);
	v[i][0] = -v[i][0];
      }
      if (yloflag && x[i][1] < ylo) {
	x[i][1] = ylo + (ylo - x[i][1]);
	v[i][1] = -v[i][1];
      }
      if (yhiflag && x[i][1] > yhi) {
	x[i][1] = yhi - (x[i][1] - yhi);
	v[i][1] = -v[i][1];
      }
      if (zloflag && x[i][2] < zlo) {
	x[i][2] = zlo + (zlo - x[i][2]);
	v[i][2] = -v[i][2];
      }
      if (zhiflag && x[i][2] > zhi) {
	x[i][2] = zhi - (x[i][2] - zhi);
	v[i][2] = -v[i][2];
      }
    }
  }
}
