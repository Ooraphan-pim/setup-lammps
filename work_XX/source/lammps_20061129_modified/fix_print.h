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

#ifndef FIX_PRINT_H
#define FIX_PRINT_H

#include "fix.h"

class FixPrint : public Fix {
 public:
  FixPrint(int, char **);
  ~FixPrint();
  int setmask();
  void end_of_step();

 private:
  int me;
  char *line,*copy,*work;
};

#endif
