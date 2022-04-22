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

/* ----------------------------------------------------------------------
   Contributing author: Mathias Puetz (SNL) and friends
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "dihedral_multi_harmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define TOLERANCE 0.05
#define SMALL     0.001

/* ----------------------------------------------------------------------
   free all arrays 
------------------------------------------------------------------------- */

DihedralMultiHarmonic::~DihedralMultiHarmonic()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(a1);
    memory->sfree(a2);
    memory->sfree(a3);
    memory->sfree(a4);
    memory->sfree(a5);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralMultiHarmonic::compute(int eflag, int vflag)
{
  int n,i1,i2,i3,i4,type,factor;
  double rfactor;
  double vb1x,vb1y,vb1z,vb2x,vb2y;
  double vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1;
  double sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22;
  double a33,a12,a13,a23,sx1,sx2,sx12,sy1,sy2,sy12;
  double sz1,sz2,sz12,s2,sin2;

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {

    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    if (newton_bond) factor = 4;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
      if (i3 < nlocal) factor++;
      if (i4 < nlocal) factor++;
      }
    rfactor = 0.25 * factor;

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
    domain->minimum_image(&vb1x,&vb1y,&vb1z);

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    domain->minimum_image(&vb2x,&vb2y,&vb2z);

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;
    domain->minimum_image(&vb2xm,&vb2ym,&vb2zm);

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    domain->minimum_image(&vb3x,&vb3y,&vb3z);

    // c0 calculation
        
    sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
    sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
    sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
        
    rb1 = sqrt(sb1);
    rb3 = sqrt(sb3);
        
    c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

    // 1st and 2nd angle
        
    b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    b1mag = sqrt(b1mag2);
    b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
    b2mag = sqrt(b2mag2);
    b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    b3mag = sqrt(b3mag2);

    ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
    r12c1 = 1.0 / (b1mag*b2mag);
    c1mag = ctmp * r12c1;

    ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
    r12c2 = 1.0 / (b2mag*b3mag);
    c2mag = ctmp * r12c2;

    // cos and sin of 2 angles and final c

    sin2 = MAX(1.0 - c1mag*c1mag,0.0);
    sc1 = sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0/sc1;

    sin2 = MAX(1.0 - c2mag*c2mag,0.0);
    sc2 = sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0/sc2;

    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0 + c1mag*c2mag) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      if (screen) {
	fprintf(screen,"Dihedral problem: %d %d %d %d %d %d\n",
		comm->me,update->ntimestep,
		atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	fprintf(screen,"  1st atom: %d %g %g %g\n",
		comm->me,x[i1][0],x[i1][1],x[i1][2]);
	fprintf(screen,"  2nd atom: %d %g %g %g\n",
		comm->me,x[i2][0],x[i2][1],x[i2][2]);
	fprintf(screen,"  3rd atom: %d %g %g %g\n",
		comm->me,x[i3][0],x[i3][1],x[i3][2]);
	fprintf(screen,"  4th atom: %d %g %g %g\n",
		comm->me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy
    // p = sum (i=1,5) a_i * c**(i-1)
    // pd = dp/dc

    p = a1[type] + c*(a2[type] + c*(a3[type] + c*(a4[type] + c*a5[type])));
    pd = a2[type] + c*(2.0*a3[type] + c*(3.0*a4[type] + c*4.0*a5[type]));

    if (eflag) energy += rfactor * p;

    a = pd;
    c = c * a;
    s12 = s12 * a;
    a11 = (-c*sb1*s1);
    a22 = sb2*(2.0*c0*s12 - c*(s1+s2));
    a33 = (-c*sb3*s2);
    a12 = r12c1*(c1mag*c*s1 + c2mag*s12);
    a13 = rb1*rb3*s12;
    a23 = r12c2*(-c2mag*c*s2 - c1mag*s12);

    sx1  = a11*vb1x + a12*vb2x + a13*vb3x;
    sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
    sx12 = a13*vb1x + a23*vb2x + a33*vb3x;
    sy1  = a11*vb1y + a12*vb2y + a13*vb3y;
    sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
    sy12 = a13*vb1y + a23*vb2y + a33*vb3y;
    sz1  = a11*vb1z + a12*vb2z + a13*vb3z;
    sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
    sz12 = a13*vb1z + a23*vb2z + a33*vb3z;

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= sx1;
      f[i1][1] -= sy1;
      f[i1][2] -= sz1;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += sx2 + sx1;
      f[i2][1] += sy2 + sy1;
      f[i2][2] += sz2 + sz1;
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += sx12 - sx2;
      f[i3][1] += sy12 - sy2;
      f[i3][2] += sz12 - sz2;
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] -= sx12;
      f[i4][1] -= sy12;
      f[i4][2] -= sz12;
    }

    // virial contribution

    if (vflag) {
      virial[0] -= rfactor * (vb1x*sx1 + vb2x*sx2 + vb3x*sx12);
      virial[1] -= rfactor * (vb1y*sy1 + vb2y*sy2 + vb3y*sy12);
      virial[2] -= rfactor * (vb1z*sz1 + vb2z*sz2 + vb3z*sz12);
      virial[3] -= rfactor * (vb1x*sy1 + vb2x*sy2 + vb3x*sy12);
      virial[4] -= rfactor * (vb1x*sz1 + vb2x*sz2 + vb3x*sz12);
      virial[5] -= rfactor * (vb1y*sz1 + vb2y*sz2 + vb3y*sz12);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DihedralMultiHarmonic::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  a1 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:a1");
  a2 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:a2");
  a3 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:a3");
  a4 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:a4");
  a5 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:a5");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralMultiHarmonic::coeff(int which, int narg, char **arg)
{
  if (which != 0) error->all("Invalid coeffs for this dihedral style");
  if (narg != 6) error->all("Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  double a1_one = atof(arg[1]);
  double a2_one = atof(arg[2]);
  double a3_one = atof(arg[3]);
  double a4_one = atof(arg[4]);
  double a5_one = atof(arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    a1[i] = a1_one;
    a2[i] = a2_one;
    a3[i] = a3_one;
    a4[i] = a4_one;
    a5[i] = a5_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralMultiHarmonic::write_restart(FILE *fp)
{
  fwrite(&a1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a3[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a4[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a5[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralMultiHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&a1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&a2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&a3[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&a4[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&a5[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&a1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a3[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a4[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a5[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}