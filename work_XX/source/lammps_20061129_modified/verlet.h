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

#ifndef VERLET_H
#define VERLET_H

#include "integrate.h"

class Verlet : public Integrate {
 public:
  Verlet(int, char **);
  ~Verlet() {}
  void init();
  void setup();
  void iterate(int);

 private:
  int virial_every;       // what vflag should be on every timestep (0,1,2)
  int virial_thermo;      // what vflag should be on thermo steps (1,2)
  int pairflag,torqueflag,granflag;

  int maxpair;            // copies of Update quantities
  double **f_pair;

  void force_clear(int);
};

#endif
