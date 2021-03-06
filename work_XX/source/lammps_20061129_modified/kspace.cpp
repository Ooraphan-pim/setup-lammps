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

#include "stdlib.h"
#include "string.h"
#include "kspace.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

KSpace::KSpace(int narg, char **arg)
{
  energy = 0.0;
  order = 5;
  gridflag = 0;
  gewaldflag = 0;
  slabflag = 0;
  slab_volfactor = 1;
}

/* ----------------------------------------------------------------------
   modify parameters of the KSpace style 
------------------------------------------------------------------------- */

void KSpace::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mesh") == 0) {
      if (iarg+4 > narg) error->all("Illegal kspace_modify command");
      nx_pppm = atoi(arg[iarg+1]);
      ny_pppm = atoi(arg[iarg+2]);
      nz_pppm = atoi(arg[iarg+3]);
      if (nx_pppm == 0 && ny_pppm == 0 && nz_pppm == 0) gridflag = 0;
      else gridflag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"order") == 0) {
      if (iarg+2 > narg) error->all("Illegal kspace_modify command");
      order = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald") == 0) {
      if (iarg+2 > narg) error->all("Illegal kspace_modify command");
      g_ewald = atof(arg[iarg+1]);
      if (g_ewald == 0.0) gewaldflag = 0;
      else gewaldflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"slab") == 0) {
      if (iarg+2 > narg) error->all("Illegal kspace_modify command");
      slab_volfactor = atof(arg[iarg+1]);
      iarg += 2;
      if (slab_volfactor < 2.0) error->all("Bad slab parameter");
      slabflag = 1;
    } else error->all("Illegal kspace_modify command");
  }
}
