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

#ifndef PAIR_SOFT_H
#define PAIR_SOFT_H

#include "pair.h"

class PairSoft : public Pair {
  friend class Pair;

 public:
  PairSoft();
  ~PairSoft();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void single(int, int, int, int, double, double, double, int, One &);

 private:
  double PI;
  double cut_global;
  double **prestart,**prestop;
  double **prefactor;
  double **cut;

  void allocate();
};

#endif
