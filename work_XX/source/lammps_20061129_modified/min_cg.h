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

#ifndef MIN_CG_H
#define MIN_CG_H

#include "min.h"

class MinCG : public Min {
 public:
  MinCG() {}
  virtual ~MinCG() {}
  void init();
  void run();

  virtual void iterate(int);

 protected:
  int virial_thermo;      // what vflag should be on thermo steps (1,2)
  int pairflag,torqueflag,granflag;
  int neigh_every,neigh_delay,neigh_dist_check;   // copies of reneigh criteria

  int maxpair;            // copies of Update quantities
  double **f_pair;

  int ifix_minimize;      // fix that stores gradient vecs
  double ecurrent;        // current potential energy
  double mindist,maxdist; // min/max dist for any coord delta in line search

  int ndof;               // # of degrees-of-freedom on this proc
  double *g,*h;           // local portion of gradient, searchdir vectors

  typedef int (MinCG::*FnPtr)(int, double *, double *, double,
			      double, double, double &, int &);
  FnPtr linemin;          // ptr to linemin functions

  int linemin_scan(int, double *, double *, double,
		   double, double, double &, int &);
  int linemin_secant(int, double *, double *, double,
		     double, double, double &, int &);

  void setup();
  void setup_vectors();
  void eng_force(int *, double **, double **, double *);
  void force_clear(int);
  void write_energy_by_type(int);
};

#endif
