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

#ifndef TIMER_H
#define TIMER_H

#include "lammps.h"

#define TIME_TOTAL       0
#define TIME_LOOP        1
#define TIME_PAIR        2
#define TIME_BOND        3
#define TIME_KSPACE      4
#define TIME_NEIGHBOR    5
#define TIME_COMM        6
#define TIME_OUTPUT      7

#define TIME_N           8

class Timer : public LAMMPS {
 public:
  double *array;

  Timer();
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

 private:
  double previous_time;
};

#endif
