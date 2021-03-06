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

#ifndef FIX_MOMENTUM_H
#define FIX_MOMENTUM_H

#include "fix.h"

class FixMomentum : public Fix {
 public:
  FixMomentum(int, char **);
  ~FixMomentum() {}
  int setmask();
  void init();
  void end_of_step();

 private:
  int linear,angular;
  int xflag,yflag,zflag;
  double masstotal;
};

#endif
