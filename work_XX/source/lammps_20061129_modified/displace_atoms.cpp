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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "displace_atoms.h"
#include "system.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "group.h"
#include "error.h"

#define MOVE 1
#define RAMP 2

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

void DisplaceAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all("Illegal displace_atoms command");

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for displace_atoms ...\n");
  sys->init();

  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");

  // group and style

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all("Could not find displace_atoms group ID");
  int groupbit = group->bitmask[igroup];

  int style;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else error->all("Illegal displace_atoms command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of displace_atoms with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms directly by specified 3-vector distance

  if (style == MOVE) {

    double delx,dely,delz;
    delx = xscale*atof(arg[2]);
    dely = yscale*atof(arg[3]);
    delz = zscale*atof(arg[4]);

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	x[i][0] += delx;
	x[i][1] += dely;
	x[i][2] += delz;
      }
    }

    // move atoms in ramped fashion

  } else if (style == RAMP) {

    int d_dim;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all("Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*atof(arg[3]);
      d_hi = xscale*atof(arg[4]);
    } else if (d_dim == 1) {
      d_lo = yscale*atof(arg[3]);
      d_hi = yscale*atof(arg[4]);
    } else if (d_dim == 2) {
      d_lo = zscale*atof(arg[3]);
      d_hi = zscale*atof(arg[4]);
    }

    int coord_dim;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all("Illegal velocity ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*atof(arg[6]);
      coord_hi = xscale*atof(arg[7]);
    } else if (coord_dim == 1) {
      coord_lo = yscale*atof(arg[6]);
      coord_hi = yscale*atof(arg[7]);
    } else if (coord_dim == 2) {
      coord_lo = zscale*atof(arg[6]);
      coord_hi = zscale*atof(arg[7]);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double fraction,dramp;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
	fraction = MAX(fraction,0.0);
	fraction = MIN(fraction,1.0);
	dramp = d_lo + fraction*(d_hi - d_lo);
	x[i][d_dim] += dramp;
      }
    }
  }

  // move atoms to new processors
  // enforce PBC before in case atoms are outside box

  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();

  // check if any atoms were lost

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (natoms != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via displacement: original %.15g current %.15g",
	    atom->natoms,natoms);
    error->all(str);
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_atoms input line 
------------------------------------------------------------------------- */

void DisplaceAtoms::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal displace_atoms command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal displace_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal displace_atoms command");
      iarg += 2;
    } else error->all("Illegal displace_atoms command");
  }
}
