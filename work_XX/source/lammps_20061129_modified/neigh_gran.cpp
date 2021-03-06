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

#include "neighbor.h"
#include "atom.h"
#include "modify.h"
#include "fix_shear_history.h"
#include "error.h"

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_no_newton()
{
  int i,j,m,n,nn;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;
  int *npartner;
  int **partner;
  double ***shearpartner;

  if (history >= 0) {
    npartner = ((FixShearHistory *) modify->fix[history])->npartner;
    partner = ((FixShearHistory *) modify->fix[history])->partner;
    shearpartner = 
      ((FixShearHistory *) modify->fix[history])->shearpartner;
  }

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) {
	add_pages(npage);
	if (history >= 0) add_pages_history(npage);
      }
    }

    n = 0;
    neighptr = &pages[npage][npnt];
    if (history >= 0) {
      nn = 0;
      touchptr = &pages_touch[npage][npnt];
      shearptr = &pages_shear[npage][3*npnt];
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (exclude && exclusion(i,j,type,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) {
	neighptr[n] = j;

	if (history >= 0) {
	  if (rsq < radsum*radsum) {
	    for (m = 0; m < npartner[i]; m++)
	      if (partner[i][m] == tag[j]) break;
	    if (m < npartner[i]) {
	      touchptr[n] = 1;
	      shearptr[nn++] = shearpartner[i][m][0];
	      shearptr[nn++] = shearpartner[i][m][1];
	      shearptr[nn++] = shearpartner[i][m][2];
	    } else {
	      touchptr[n] = 0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	    }
	  } else {
	    touchptr[n] = 0;
	    shearptr[nn++] = 0.0;
	    shearptr[nn++] = 0.0;
	    shearptr[nn++] = 0.0;
	  }
	}

	n++;
      }
    }	       

    if (history >= 0) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
    }
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   no shear history is allowed for this option
   pair added to list if atoms i and j are both owned and i < j
   if j is ghost only me or other proc adds pair
   decision based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_newton()
{
  int i,j,n,itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }

    n = 0;
    neighptr = &pages[npage][npnt];

    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (j >= nlocal) {
	jtag = tag[j];
	if (itag > jtag) {
	  if ((itag+jtag) % 2 == 0) continue;
	} else if (itag < jtag) {
	  if ((itag+jtag) % 2 == 1) continue;
	} else {
	  if (x[j][2] < ztmp) continue;
	  else if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	  else if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp)
	    continue;
	}
      }

      if (exclude && exclusion(i,j,type,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);
      
      if (rsq <= cutsq) neighptr[n++] = j;
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   each owned atom i checks own bin and surrounding bins in non-Newton stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_bin_no_newton()
{
  int i,j,k,m,n,nn,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;
  int *npartner;
  int **partner;
  double ***shearpartner;

  if (history >= 0) {
    npartner = ((FixShearHistory *) modify->fix[history])->npartner;
    partner = ((FixShearHistory *) modify->fix[history])->partner;
    shearpartner = 
      ((FixShearHistory *) modify->fix[history])->shearpartner;
  }

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) {
	add_pages(npage);
	if (history >= 0) add_pages_history(npage);
      }
    }

    n = 0;
    neighptr = &pages[npage][npnt];
    if (history >= 0) {
      nn = 0;
      touchptr = &pages_touch[npage][npnt];
      shearptr = &pages_shear[npage][3*npnt];
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    ibin = coord2bin(x[i]);

    // loop over all atoms in surrounding bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    for (k = 0; k < nstencil; k++) {
      j = binhead[ibin+stencil[k]];
      while (j >= 0) {
	if (j <= i) {
	  j = bins[j];
	  continue;
	}

	if (exclude && exclusion(i,j,type,mask,molecule)) {
	  j = bins[j];
	  continue;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radi + radius[j];
	cutsq = (radsum+skin) * (radsum+skin);

	if (rsq <= cutsq) {
	  neighptr[n] = j;

	  if (history >= 0) {
	    if (rsq < radsum*radsum) {
	      for (m = 0; m < npartner[i]; m++)
		if (partner[i][m] == tag[j]) break;
	      if (m < npartner[i]) {
		touchptr[n] = 1;
		shearptr[nn++] = shearpartner[i][m][0];
		shearptr[nn++] = shearpartner[i][m][1];
		shearptr[nn++] = shearpartner[i][m][2];
	      } else {
		touchptr[n] = 0;
		shearptr[nn++] = 0.0;
		shearptr[nn++] = 0.0;
		shearptr[nn++] = 0.0;
	      }
	    } else {
	      touchptr[n] = 0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	    }
	  }

	  n++;
	}

	j = bins[j];
      }
    }

    if (history >= 0) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
    }
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with full Newton's 3rd law
   no shear history is allowed for this option
   every pair stored exactly once by some processor
   each owned atom i checks its own bin and other bins in Newton stencil
------------------------------------------------------------------------- */

void Neighbor::granular_bin_newton()
{
  int i,j,k,n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }

    n = 0;
    neighptr = &pages[npage][npnt];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i

    j = bins[i];
    while (j >= 0) {
      if (j >= nlocal) {
	if ((x[j][2] < ztmp) || (x[j][2] == ztmp && x[j][1] < ytmp) ||
	    (x[j][2] == ztmp && x[j][1]  == ytmp && x[j][0] < xtmp)) {
	  j = bins[j];
	  continue;
	}
      }

      if (exclude && exclusion(i,j,type,mask,molecule)) {
	j = bins[j];
	continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) neighptr[n++] = j;

      j = bins[j];
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      j = binhead[ibin+stencil[k]];
      while (j >= 0) {
	if (exclude && exclusion(i,j,type,mask,molecule)) {
	  j = bins[j];
	  continue;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radi + radius[j];
	cutsq = (radsum+skin) * (radsum+skin);

	if (rsq <= cutsq) neighptr[n++] = j;

	j = bins[j];
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}
