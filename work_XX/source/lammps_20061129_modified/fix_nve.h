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

#ifndef FIX_NVE_H
#define FIX_NVE_H

#include "fix.h"

class FixNVE : public Fix {
 public:
  FixNVE(int, char **);
  ~FixNVE() {}
  int setmask();
  void init();
  void initial_integrate();
  void final_integrate();
  void initial_integrate_respa(int, int);
  void final_integrate_respa(int);

 private:
  double dtv,dtf;
  double *step_respa;
  int mass_require;
};

#endif
