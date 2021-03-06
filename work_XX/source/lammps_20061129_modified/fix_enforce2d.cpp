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
#include "fix_enforce2d.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FixEnforce2D::FixEnforce2D(int narg, char **arg) : Fix(narg, arg)
{
  if (narg != 3) error->all("Illegal fix enforce2d command");
}

/* ---------------------------------------------------------------------- */

int FixEnforce2D::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::init()
{
  granular = atom->check_style("granular");
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::setup()
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

void FixEnforce2D::min_setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force(int vflag)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][2] = 0.0;
      f[i][2] = 0.0;
    }

  // for granular systems, zero xy rotational componenets

  if (granular) {
    double **phiv = atom->phiv;
    double **phia = atom->phia;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	phiv[i][0] = 0.0;
	phiv[i][1] = 0.0;
	phia[i][0] = 0.0;
	phia[i][1] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::min_post_force(int vflag)
{
  post_force(vflag);
}
