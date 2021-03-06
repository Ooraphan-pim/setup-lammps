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

#include "mpi.h"
#include "fft3d_wrap.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FFT3d::FFT3d(MPI_Comm comm, int nfast, int nmid, int nslow,
	     int in_ilo, int in_ihi, int in_jlo, int in_jhi,
	     int in_klo, int in_khi,
	     int out_ilo, int out_ihi, int out_jlo, int out_jhi,
	     int out_klo, int out_khi,
	     int scaled, int permute, int *nbuf)
{
  plan = fft_3d_create_plan(comm,nfast,nmid,nslow,
			    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
			    scaled,permute,nbuf);
  if (plan == NULL) error->one("Could not create 3d FFT plan");
}

/* ---------------------------------------------------------------------- */

FFT3d::~FFT3d()
{
  fft_3d_destroy_plan(plan);
}

/* ---------------------------------------------------------------------- */

void FFT3d::compute(double *in, double *out, int flag)
{
  fft_3d((FFT_DATA *) in,(FFT_DATA *) out,flag,plan);
}

/* ---------------------------------------------------------------------- */

void FFT3d::timing1d(double *in, int nsize, int flag)
{
  fft_1d_only((FFT_DATA *) in,nsize,flag,plan);
}
