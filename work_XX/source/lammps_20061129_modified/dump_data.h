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

#ifndef DUMP_DATA_H
#define DUMP_DATA_H

#include "dump.h"

class DumpData : public Dump {
 public:
  DumpData(int, char**);
  ~DumpData() {}
  void init();

 private:
  int scale_flag;            // 1 if atom coords are scaled, 0 if no
  int image_flag;            // 1 if append box count to atom coords, 0 if no

  int modify_param(int, char **);
  void write_header(int);
  void write_vel_header(int);
  int count();
  int pack();
  void write_data(int, double *);
  void write_vel(int, double*);

  typedef void (DumpData::*FnPtrHeader)(int);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int);
  void header_item(int);
  void header_dlo(int);

  typedef int (DumpData::*FnPtrPack)();
  FnPtrPack pack_choice;               // ptr to pack functions
  int pack_scale_image();
  int pack_scale_noimage();
  int pack_noscale_image();
  int pack_noscale_noimage();

  typedef void (DumpData::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_image(int, double *);
  void write_noimage(int, double *);
};

#endif
