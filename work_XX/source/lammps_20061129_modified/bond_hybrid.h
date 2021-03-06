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

#ifndef BOND_HYBRID_H
#define BOND_HYBRID_H

#include "stdio.h"
#include "bond.h"

class BondHybrid : public Bond {
  friend class Force;

 public:
  BondHybrid();
  ~BondHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void single(int, double, int, int, double, int, double &, double &);
  int memory_usage();

 private:
  int nstyles;                  // # of different bond styles
  Bond **styles;                // class list for each Bond style
  char **keywords;              // keyword for each Bond style
  int *map;                     // which style each bond type points to

  int *nbondlist;               // # of bonds in sub-style bondlists
  int *maxbond;                 // max # of bonds sub-style lists can store
  int ***bondlist;              // bondlist for each sub-style
  
  void allocate();
};

#endif
