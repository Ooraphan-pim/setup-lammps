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

#ifndef DUMP_H
#define DUMP_H

#include "stdio.h"
#include "lammps.h"

class Dump : public LAMMPS {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  int igroup,groupbit;       // group that Dump is performed on
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int vel_flag;              // 0 if no velocity section, 1 if velocity section

  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format;              // format string for the file write
  double *buf;               // memory for atom quantities
  int maxbuf;                // size of buf
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one atom

  Dump(int, char **);
  virtual ~Dump();
  virtual void init() {}
  void write();
  void modify_params(int, char **);
  virtual int memory_usage();

 protected:
  virtual void openfile();
  virtual int modify_param(int, char **) {return 0;}

  virtual void write_header(int) = 0;
  virtual int count() = 0;
  virtual int pack() = 0;
  virtual void write_data(int, double *) = 0;
  virtual void write_vel_header(int) = 0;
  virtual void write_vel(int, double *) = 0;
};

#endif
