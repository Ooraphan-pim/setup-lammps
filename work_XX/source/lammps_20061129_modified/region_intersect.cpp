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

#include "stdlib.h"
#include "string.h"
#include "region_intersect.h"
#include "domain.h"
#include "error.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

RegIntersect::RegIntersect(int narg, char **arg) : Region(narg, arg)
{
  if (narg < 5) error->all("Illegal region command");
  int n = atoi(arg[2]);
  if (n < 2) error->all("Illegal region command");
  options(narg-(n+3),&arg[n+3]);

  // build list of regions to intersect

  list = new int[n];
  nregion = 0;

  int iregion;
  for (int iarg = 0; iarg < n; iarg++) {
    for (iregion = 0; iregion < domain->nregion; iregion++)
      if (strcmp(arg[iarg+3],domain->regions[iregion]->id) == 0) break;
    if (iregion == domain->nregion)
      error->all("Region union region ID does not exist");
    list[nregion++] = iregion;
  }

  // extent of intersection of regions

  Region **regions = domain->regions;

  extent_xlo = regions[list[0]]->extent_xlo;
  extent_ylo = regions[list[0]]->extent_ylo;
  extent_zlo = regions[list[0]]->extent_zlo;
  extent_xhi = regions[list[0]]->extent_xhi;
  extent_yhi = regions[list[0]]->extent_yhi;
  extent_zhi = regions[list[0]]->extent_zhi;

  for (int ilist = 1; ilist < nregion; ilist++) {
    extent_xlo = MAX(extent_xlo,regions[list[ilist]]->extent_xlo);
    extent_ylo = MAX(extent_ylo,regions[list[ilist]]->extent_ylo);
    extent_zlo = MAX(extent_zlo,regions[list[ilist]]->extent_zlo);
    extent_xhi = MIN(extent_xhi,regions[list[ilist]]->extent_xhi);
    extent_yhi = MIN(extent_yhi,regions[list[ilist]]->extent_yhi);
    extent_zhi = MIN(extent_zhi,regions[list[ilist]]->extent_zhi);
  }
}

/* ---------------------------------------------------------------------- */

RegIntersect::~RegIntersect()
{
  delete [] list;
}

/* ---------------------------------------------------------------------- */

int RegIntersect::match(double x, double y, double z)
{
  int ilist;
  Region **regions = domain->regions;
  for (ilist = 0; ilist < nregion; ilist++)
    if (!regions[list[ilist]]->match(x,y,z)) break;

  int inside;                          // inside if matched all regions
  if (ilist == nregion) inside = 1;
  else inside = 0;

  return !(inside ^ interior);         // 1 if same, 0 if different
}
