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

#ifndef FIX_H
#define FIX_H

#include "lammps.h"

class Fix : public LAMMPS {
 public:
  char *id,*style;
  int igroup,groupbit;
  int restart_global;            // 1 if Fix saves global state, 0 if not
  int restart_peratom;           // 1 if Fix saves peratom state, 0 if not
  int force_reneighbor;          // 1 if Fix forces reneighboring, 0 if not
  int next_reneighbor;           // next timestep to force a reneighboring
  int thermo_print;              // 1 if Fix prints info during thermo, 0 no
  int thermo_energy;             // 1 if Fix adds to thermo energy, 0 if not
  int nevery;                    // how often to call an end_of_step fix
  int neigh_half_once;           // 0/1 if needs half neigh list occasionally
  int neigh_half_every;          // 0/1 if needs half neigh list every step
  int neigh_full_once;           // 0/1 if needs full neigh list occasionally
  int neigh_full_every;          // 0/1 if needs full neigh list every step
  double virial[6];              // fix contribution to pressure virial

  int INITIAL_INTEGRATE,PRE_EXCHANGE,PRE_NEIGHBOR;    // mask settings
  int POST_FORCE,FINAL_INTEGRATE,END_OF_STEP,THERMO;
  int INITIAL_INTEGRATE_RESPA,POST_FORCE_RESPA,FINAL_INTEGRATE_RESPA;
  int MIN_POST_FORCE;

  Fix(int, char **);
  virtual ~Fix();
  void modify_params(int, char **);

  virtual int setmask() = 0;

  virtual void init() {}
  virtual void setup() {}
  virtual void min_setup() {}
  virtual void initial_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual void end_of_step() {}
  virtual void write_restart(FILE *) {}
  virtual void restart(char *) {}

  virtual int memory_usage() {return 0;}
  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int) {}
  virtual int pack_exchange(int, double *) {return 0;}
  virtual int unpack_exchange(int, double *) {return 0;}
  virtual int pack_restart(int, double *) {return 0;}
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) {return 0;}
  virtual int maxsize_restart() {return 0;}

  virtual void initial_integrate_respa(int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int) {}

  virtual void min_post_force(int) {}

  virtual int pack_comm(int, int *, double *, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual int thermo_fields(int, int *, char **) {return 0;}
  virtual int thermo_compute(double *) {return 0;}
  virtual void dump() {}

  virtual int dof(int) {return 0;}
  virtual void dilate(int, double, double, double, double) {}

  virtual int modify_param(int, char **) {return 0;}
};

#endif
