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
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "velocity.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "temperature.h"
#include "temp_full.h"
#include "random_park.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#define CREATE   1
#define SET      2
#define SCALE    3
#define RAMP     4
#define ZERO     5

#define ALL      1
#define LOCAL    2
#define GEOM     3

#define WARMUP 100
#define SMALL  0.001

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

void Velocity::command(int narg, char **arg)
{
  // require atom masses to all be set

  if (domain->box_exist == 0) 
    error->all("Velocity command before simulation box is defined");
  if (atom->natoms == 0)
    error->all("Velocity command with no atoms existing");
  atom->check_mass();

  if (narg < 2) error->all("Illegal velocity command");

  // identify group

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all("Could not find velocity group ID");
  groupbit = group->bitmask[igroup];

  // identify style

  if (strcmp(arg[1],"create") == 0) style = CREATE;
  else if (strcmp(arg[1],"set") == 0) style = SET;
  else if (strcmp(arg[1],"scale") == 0) style = SCALE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"zero") == 0) style = ZERO;
  else error->all("Illegal velocity command");

  // set defaults

  tempwhich = -1;
  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 1;
  rotation_flag = 0;
  loop_flag = ALL;
  scale_flag = 1;

  // read options from end of input line
  // change defaults as options specify

  if (style == CREATE) options(narg-4,&arg[4]);
  else if (style == SET) options(narg-5,&arg[5]);
  else if (style == SCALE) options(narg-3,&arg[3]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == ZERO) options(narg-3,&arg[3]);

  // set scaling for SET and RAMP styles

  if (style == SET || style == RAMP) {
    if (scale_flag && domain->lattice == NULL)
      error->all("Use of velocity with undefined lattice");

    if (scale_flag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;
  }

  // initialize velocities based on style

  if (style == CREATE) create(narg-2,&arg[2]);
  else if (style == SET) set(narg-2,&arg[2]);
  else if (style == SCALE) scale(narg-2,&arg[2]);
  else if (style == RAMP) ramp(narg-2,&arg[2]);
  else if (style == ZERO) zero(narg-2,&arg[2]);
}

/* ---------------------------------------------------------------------- */

void Velocity::create(int narg, char **arg)
{
  int i;

  double t_desired = atof(arg[0]);
  int seed = atoi(arg[1]);

  // if tempwhich = -1, create a new temperature full style with the vel group
  // else use pre-defined temperature

  Temperature *temperature;
  if (tempwhich == -1) {
    char **arg = new char*[3];
    arg[0] = "temp";
    arg[1] = group->names[igroup];
    arg[2] = "full";
    temperature = new TempFull(3,arg);
    delete [] arg;
  } else temperature = force->templist[tempwhich];

  // initialize temperature computation
  // warn if groups don't match

  if (igroup != temperature->igroup && comm->me == 0)
    error->warning("Mismatch between velocity and temperature groups");
  temperature->init();

  // store a copy of current velocities

  double **v = atom->v;
  int nlocal = atom->nlocal;
  double **vhold = memory->create_2d_double_array(nlocal,3,"velocity:vnew");

  for (i = 0; i < nlocal; i++) {
    vhold[i][0] = v[i][0];
    vhold[i][1] = v[i][1];
    vhold[i][2] = v[i][2];
  }

  // create new velocities, in uniform or gaussian distribution
  // loop option determines looping style, ALL is default
  //   ALL = loop over all natoms, only set those I own via atom->map
  //    cannot do this if atom IDs do not span 1-Natoms (some were deleted)
  //    will produce same V, independent of P, if atoms were read-in
  //    will NOT produce same V, independent of P, if used create_atoms
  //   LOCAL = only loop over my atoms, adjust RNG to be proc-specific
  //    will never produce same V, independent of P
  //   GEOM = only loop over my atoms
  //    choose RNG for each atom based on its xyz coord (geometry)
  //    will always produce same V, independent of P
  // adjust by factor for atom mass
  // for 2d, set Vz to 0.0

  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int dimension = force->dimension;
  int mass_require = atom->mass_require;

  int m;
  double vx,vy,vz,factor;
  RanPark *random;
  
  if (loop_flag == ALL) {

    // create an atom map if one doesn't exist already

    int mapflag = 0;
    if (atom->map_style == 0) {
      mapflag = 1;
      atom->map_style = 1;
      atom->map_init();
      atom->map_set();
    }

    random = new RanPark(seed);

    if (atom->tag_enable == 0)
      error->all("Cannot use velocity create loop all unless atoms have IDs");
    int natoms = static_cast<int> (atom->natoms);

    // check that atom IDs span range from 1 to natoms

    int *tag = atom->tag;
    int idmin = natoms;
    int idmax = 0;

    for (i = 0; i < nlocal; i++) {
      idmin = MIN(idmin,tag[i]);
      idmax = MAX(idmax,tag[i]);
    }
    int idminall,idmaxall;
    MPI_Allreduce(&idmin,&idminall,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&idmax,&idmaxall,1,MPI_INT,MPI_MAX,world);

    if (idminall != 1 || idmaxall != natoms) {
      char *str = "Cannot use velocity create loop all with non-contiguous atom IDs";
      error->all(str);
    }
    
    // loop over all atoms in system
    // generate RNGs for all atoms, only assign to ones I own
    // use either per-type mass or per-atom rmass

    for (i = 1; i <= natoms; i++) {
      if (dist_flag == 0) {
	vx = random->uniform();
	vy = random->uniform();
	vz = random->uniform();
      } else {
	vx = random->gaussian();
	vy = random->gaussian();
	vz = random->gaussian();
      }
      m = atom->map(i);
      if (m >= 0 && m < nlocal) {
	if (mask[m] & groupbit) {
	  if (mass_require) factor = 1.0/sqrt(mass[type[m]]);
	  else factor = 1.0/sqrt(rmass[m]);
	  v[m][0] = vx * factor;
	  v[m][1] = vy * factor;
	  if (dimension == 3) v[m][2] = vz * factor;
	  else v[m][2] = 0.0;
	}
      }
    }

    // delete temporary atom map

    if (mapflag) {
      atom->map_delete();
      atom->map_style = 0;
    }

  } else if (loop_flag == LOCAL) {
    random = new RanPark(seed + comm->me);
    for (i = 0; i < WARMUP; i++) random->uniform();

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (dist_flag == 0) {
	  vx = random->uniform();
	  vy = random->uniform();
	  vz = random->uniform();
	} else {
	  vx = random->gaussian();
	  vy = random->gaussian();
	  vz = random->gaussian();
	}
	if (mass_require) factor = 1.0/sqrt(mass[type[i]]);
	else factor = 1.0/sqrt(rmass[i]);
	v[i][0] = vx * factor;
	v[i][1] = vy * factor;
	if (dimension == 3) v[i][2] = vz * factor;
	else v[i][2] = 0.0;
      }
    }

  } else if (loop_flag == GEOM) {
    random = new RanPark(seed);
    double **x = atom->x;
    
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	triple(x[i][0],x[i][1],x[i][2],&vx,&vy,&vz,seed,random);
	if (mass_require) factor = 1.0/sqrt(mass[type[i]]);
	else factor = 1.0/sqrt(rmass[i]);
	v[i][0] = vx * factor;
	v[i][1] = vy * factor;
	if (dimension == 3) v[i][2] = vz * factor;
	else v[i][2] = 0.0;
      }
    }
  }

  // apply momentum and rotation zeroing

  if (momentum_flag) zero_momentum();
  if (rotation_flag) zero_rotation();

  // scale temp to desired value

  double t = temperature->compute();
  rescale(t,t_desired);

  // if sum_flag set, add back in previous velocities 

  if (sum_flag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	v[i][0] += vhold[i][0];
	v[i][1] += vhold[i][1];
	v[i][2] += vhold[i][2];
      }
    }
  }

  // free local memory
  // if temperature was created, delete it

  memory->destroy_2d_double_array(vhold);
  delete random;
  if (tempwhich == -1) delete temperature;
}

/* ---------------------------------------------------------------------- */

void Velocity::set(int narg, char **arg)
{
  int xflag,yflag,zflag;
  double vx,vy,vz;

  if (strcmp(arg[0],"NULL") == 0) xflag = 0;
  else {
    xflag = 1;
    vx = xscale * atof(arg[0]);
  }
  if (strcmp(arg[1],"NULL") == 0) yflag = 0;
  else {
    yflag = 1;
    vy = yscale * atof(arg[1]);
  }
  if (strcmp(arg[2],"NULL") == 0) zflag = 0;
  else {
    zflag = 1;
    vz = zscale * atof(arg[2]);
  }

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int dimension = force->dimension;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (sum_flag == 0) {
	if (xflag) v[i][0] = vx;
	if (yflag) v[i][1] = vy;
	if (zflag && dimension == 3) v[i][2] = vz;
      } else {
	if (xflag) v[i][0] += vx;
	if (yflag) v[i][1] += vy;
	if (zflag && dimension == 3) v[i][2] += vz;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   rescale velocities of a group after computing its temperature 
------------------------------------------------------------------------- */

void Velocity::scale(int narg, char **arg)
{
  double t_desired = atof(arg[0]);

  // if tempwhich = -1, create a new temperature full style with the vel group
  // else use pre-defined temperature

  Temperature *temperature;
  if (tempwhich == -1) {
    char **arg = new char*[3];
    arg[0] = "temp";
    arg[1] = group->names[igroup];
    arg[2] = "full";
    temperature = new TempFull(3,arg);
    delete [] arg;
  } else temperature = force->templist[tempwhich];

  // initialize temperature computation
  // warn if groups don't match

  if (igroup != temperature->igroup && comm->me == 0)
    error->warning("Mismatch between velocity and temperature groups");
  temperature->init();

  // scale temp to desired value

  double t = temperature->compute();
  rescale(t,t_desired);

  // if temperature was created, delete it

  if (tempwhich == -1) delete temperature;
}

/* ----------------------------------------------------------------------
   apply a ramped set of velocities 
------------------------------------------------------------------------- */

void Velocity::ramp(int narg, char **arg)
{
  int v_dim;
  if (strcmp(arg[0],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[0],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[0],"vz") == 0) v_dim = 2;
  else error->all("Illegal velocity command");

  if (v_dim == 2 && force->dimension == 2) 
    error->all("Velocity ramp in z for a 2d problem");

  double v_lo,v_hi;
  if (v_dim == 0) {
    v_lo = xscale*atof(arg[1]);
    v_hi = xscale*atof(arg[2]);
  } else if (v_dim == 1) {
    v_lo = yscale*atof(arg[1]);
    v_hi = yscale*atof(arg[2]);
  } else if (v_dim == 0) {
    v_lo = zscale*atof(arg[1]);
    v_hi = zscale*atof(arg[2]);
  }

  int coord_dim;
  if (strcmp(arg[3],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[3],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[3],"z") == 0) coord_dim = 2;
  else error->all("Illegal velocity command");

  double coord_lo,coord_hi;
  if (coord_dim == 0) {
    coord_lo = xscale*atof(arg[4]);
    coord_hi = xscale*atof(arg[5]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*atof(arg[4]);
    coord_hi = yscale*atof(arg[5]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*atof(arg[4]);
    coord_hi = zscale*atof(arg[5]);
  }

  // vramp = ramped velocity component for v_dim
  // add or set based on sum_flag
  
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fraction,vramp;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      if (sum_flag) v[i][v_dim] += vramp;
      else v[i][v_dim] = vramp;
    }
}

/* ----------------------------------------------------------------------
   zero linear or angular momentum of a group
------------------------------------------------------------------------- */

void Velocity::zero(int narg, char **arg)
{
  if (strcmp(arg[0],"linear") == 0) zero_momentum();
  else if (strcmp(arg[0],"angular") == 0) zero_rotation();
  else error->all("Illegal velocity command");
}

/* ----------------------------------------------------------------------
   rescale velocities of group atoms to t_new from t_old 
------------------------------------------------------------------------- */

void Velocity::rescale(double t_old, double t_new)
{
  if (t_old == 0.0) error->all("Attempting to rescale a 0.0 temperature");

  double factor = sqrt(t_new/t_old);

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] *= factor;
      v[i][1] *= factor;
      v[i][2] *= factor;
    }
}

/* ----------------------------------------------------------------------
   zero the linear momentum of a group of atoms by adjusting v by -Vcm
------------------------------------------------------------------------- */

void Velocity::zero_momentum()
{
  // cannot have 0 atoms in group

  if (group->count(igroup) == 0.0)
    error->all("Cannot zero momentum of 0 atoms");

  // compute velocity of center-of-mass of group

  double masstotal = group->mass(igroup);
  double vcm[3];
  group->vcm(igroup,masstotal,vcm);

  // adjust velocities by vcm to zero linear momentum

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vcm[0];
      v[i][1] -= vcm[1];
      v[i][2] -= vcm[2];
    }
}

/* ----------------------------------------------------------------------
   zero the angular momentum of a group of atoms by adjusting v by -(w x r)
------------------------------------------------------------------------- */

void Velocity::zero_rotation()
{
  int i,j;

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0.0)
    error->all("Cannot zero momentum of 0 atoms");

  // compute omega (angular velocity) of group around center-of-mass

  double xcm[3],angmom[3],inertia[3][3],omega[3];
  double masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  group->angmom(igroup,xcm,angmom);
  group->inertia(igroup,xcm,inertia);
  group->omega(igroup,angmom,inertia,omega);

  // adjust velocities to zero omega
  // vnew_i = v_i - w x r_i
  // must use unwrapped coords to compute r_i correctly

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;
  
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = (x[i][0] + xbox*xprd) - xcm[0];
      dy = (x[i][1] + ybox*yprd) - xcm[1];
      dz = (x[i][2] + zbox*zprd) - xcm[2];
      v[i][0] -= omega[1]*dz - omega[2]*dy;
      v[i][1] -= omega[2]*dx - omega[0]*dy;
      v[i][2] -= omega[0]*dy - omega[1]*dx;
    }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of velocity input line 
------------------------------------------------------------------------- */

void Velocity::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal velocity command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dist") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"uniform") == 0) dist_flag = 0;
      else if (strcmp(arg[iarg+1],"gaussian") == 0) dist_flag = 1;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) sum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) sum_flag = 1;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mom") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) momentum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) momentum_flag = 1;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rot") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) rotation_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) rotation_flag = 1;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      for (tempwhich = 0; tempwhich < force->ntemp; tempwhich++)
	if (strcmp(arg[iarg+1],force->templist[tempwhich]->id) == 0) break;
      if (tempwhich == force->ntemp) 
	error->all("Could not find velocity temperature ID");
      iarg += 2;
    } else if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"all") == 0) loop_flag = ALL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop_flag = LOCAL;
      else if (strcmp(arg[iarg+1],"geom") == 0) loop_flag = GEOM;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal velocity command");
      if (strcmp(arg[iarg+1],"box") == 0) scale_flag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scale_flag = 1;
      else error->all("Illegal velocity command");
      iarg += 2;
    } else error->all("Illegal velocity command");
  }
}

/* ---------------------------------------------------------------------- */

#define IA1 1366
#define IC1 150889
#define IM1 714025
#define IA2 8121
#define IC2 28411
#define IM2 134456
#define IA3 7141
#define IC3 54773
#define IM3 259200

void Velocity::triple(double x, double y, double z, 
		      double *vx, double *vy, double *vz,
		      int seed, RanPark *random)
{
  // seed 1,2,3 = combination of atom coord in each dim and user-input seed
  // map geometric extent into range of each of 3 RNGs
  // warm-up each RNG by calling it twice

  double fraction;
  int seed1,seed2,seed3;

  fraction = (x - domain->boxxlo) / domain->xprd;
  seed1 = static_cast<int> (fraction * IM1);
  seed1 = (seed1+seed) % IM1;
  seed1 = (seed1*IA1+IC1) % IM1;
  seed1 = (seed1*IA1+IC1) % IM1;

  fraction = (y - domain->boxylo) / domain->yprd;
  seed2 = static_cast<int> (fraction * IM2);
  seed2 = (seed2+seed) % IM2;
  seed2 = (seed2*IA2+IC2) % IM2;
  seed2 = (seed2*IA2+IC2) % IM2;

  fraction = (z - domain->boxzlo) / domain->zprd;
  seed3 = static_cast<int> (fraction * IM3);
  seed3 = (seed3+seed) % IM3;
  seed3 = (seed3*IA3+IC3) % IM3;
  seed3 = (seed3*IA3+IC3) % IM3;

  // fraction = 0-1 with giving each dim an equal weighting
  // use fraction to reset Park/Miller RNG seed

  fraction = 1.0*seed1/(3*IM1) + 1.0*seed2/(3*IM2) + 1.0*seed3/(3*IM3);
  random->reset(fraction);

  // use RNG to set velocities after warming up twice

  random->uniform();
  random->uniform();

  if (dist_flag == 0) {
    *vx = random->uniform();
    *vy = random->uniform();
    *vz = random->uniform();
  } else {
    *vx = random->gaussian();
    *vy = random->gaussian();
    *vz = random->gaussian();
  }
}
