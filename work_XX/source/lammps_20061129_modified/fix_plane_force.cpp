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
#include "stdlib.h"
#include "fix_plane_force.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FixPlaneForce::FixPlaneForce(int narg, char **arg) : Fix(narg, arg)
{
  if (narg != 6) error->all("Illegal fix planeforce command");
  xdir = atof(arg[3]);
  ydir = atof(arg[4]);
  zdir = atof(arg[5]);
}

/* ---------------------------------------------------------------------- */

int FixPlaneForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlaneForce::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(1,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPlaneForce::min_setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixPlaneForce::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dot;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dot = f[i][0]*xdir + f[i][1]*ydir + f[i][2]*zdir;
      f[i][0] -= dot * xdir;
      f[i][1] -= dot * ydir;
      f[i][2] -= dot * zdir;
    }
}

/* ---------------------------------------------------------------------- */

void FixPlaneForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPlaneForce::min_post_force(int vflag)
{
  post_force(vflag);
}
