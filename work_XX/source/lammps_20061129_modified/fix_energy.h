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

#ifndef FIX_ENERGY_H
#define FIX_ENERGY_H

#include "fix.h"

class FixEnergy : public Fix {
  friend class DumpCustom;
  friend class MinCG;

 public:
  FixEnergy(int, char **);
  ~FixEnergy();
  int setmask();
  void init();
  void dump();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int memory_usage();

 private:
  int nmax,eamstyle;
  double *energy;
};

#endif
