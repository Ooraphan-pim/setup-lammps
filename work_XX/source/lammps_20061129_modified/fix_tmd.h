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

#ifndef FIX_TMD_H
#define FIX_TMD_H

#include "stdio.h"
#include "fix.h"

class FixTMD : public Fix {
 public:
  FixTMD(int, char **);
  ~FixTMD();
  int setmask();
  void init();
  void initial_integrate();
  void initial_integrate_respa(int,int);

  int memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:
  int me;
  int nfileevery,previous_stat;
  FILE *fp;
  double rho_start,rho_stop,rho_old,masstotal;
  double dtv,dtf,dtfm;
  double work_lambda,work_analytical;
  double **xf,**xold;

  void readfile(char *);
  void open(char *);
};

#endif
