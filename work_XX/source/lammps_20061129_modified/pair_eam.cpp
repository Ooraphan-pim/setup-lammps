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
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eam.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "neighbor.h"
#include "memory.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairEAM::PairEAM()
{
  nmax = 0;
  rho = NULL;
  fp = NULL;

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAM::~PairEAM()
{
  memory->sfree(rho);
  memory->sfree(fp);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    delete [] map;
    delete [] type2frho;
    memory->destroy_2d_int_array(type2rhor);
    memory->destroy_2d_int_array(type2z2r);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete [] funcfl[i].file;
      memory->sfree(funcfl[i].frho);
      memory->sfree(funcfl[i].rhor);
      memory->sfree(funcfl[i].zr);
    }
    memory->sfree(funcfl);
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy_2d_double_array(setfl->frho);
    memory->destroy_2d_double_array(setfl->rhor);
    memory->destroy_3d_double_array(setfl->z2r);
    delete setfl;
  }

  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy_2d_double_array(fs->frho);
    memory->destroy_3d_double_array(fs->rhor);
    memory->destroy_3d_double_array(fs->z2r);
    delete fs;
  }

  memory->destroy_2d_double_array(frho);
  memory->destroy_2d_double_array(rhor);
  memory->destroy_2d_double_array(z2r);

  memory->destroy_3d_double_array(frho_spline);
  memory->destroy_3d_double_array(rhor_spline);
  memory->destroy_3d_double_array(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairEAM::compute(int eflag, int vflag)
{
  int i,j,k,m,numneigh,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,p,fforce,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  double *coeff;
  int *neighs;
  double **f;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(rho);
    memory->sfree(fp);
    nmax = atom->nmax;
    rho = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho");
    fp = (double *) memory->smalloc(nmax*sizeof(double),"pair:fp");
  }

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	if (newton_pair || j < nlocal) {
	  coeff = rhor_spline[type2rhor[itype][jtype]][m];
	  rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	}
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  for (i = 0; i < nlocal; i++) {
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) eng_vdwl += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  }

  // communicate derivative of embedding function

  comm->comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	r = sqrt(rsq);
	p = r*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);

	// rhoip = derivative of (density at atom j due to atom i)
	// rhojp = derivative of (density at atom i due to atom j)
	// phi = pair potential energy
	// phip = phi'
	// z2 = phi * r
	// z2p = (phi * r)' = (phi' r) + phi
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

	coeff = rhor_spline[type2rhor[itype][jtype]][m];
	rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = z2r_spline[type2z2r[itype][jtype]][m];
	z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	
	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	fforce = psip*recip;

	f[i][0] -= delx*fforce;
	f[i][1] -= dely*fforce;
	f[i][2] -= delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] += delx*fforce;
	  f[j][1] += dely*fforce;
	  f[j][2] += delz*fforce;
	}

	if (eflag) {
	  if (newton_pair || j < nlocal) eng_vdwl += phi;
	  else eng_vdwl += 0.5*phi;
	}

	if (vflag == 1) {
	  if (newton_pair || j < nlocal) {
	    virial[0] -= delx*delx*fforce;
	    virial[1] -= dely*dely*fforce;
	    virial[2] -= delz*delz*fforce;
	    virial[3] -= delx*dely*fforce;
	    virial[4] -= delx*delz*fforce;
	    virial[5] -= dely*delz*fforce;
	  } else {
	    virial[0] -= 0.5*delx*delx*fforce;
	    virial[1] -= 0.5*dely*dely*fforce;
	    virial[2] -= 0.5*delz*delz*fforce;
	    virial[3] -= 0.5*delx*dely*fforce;
	    virial[4] -= 0.5*delx*delz*fforce;
	    virial[5] -= 0.5*dely*delz*fforce;
	  }
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[n] = -1;

  type2frho = new int[n+1];
  type2rhor = memory->create_2d_int_array(n+1,n+1,"pair:type2rhor");
  type2z2r = memory->create_2d_int_array(n+1,n+1,"pair:type2z2r");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEAM::settings(int narg, char **arg)
{
  if (narg > 0) error->all("Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairEAM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all("Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2],funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *) 
      memory->srealloc(funcfl,nfuncfl*sizeof(Funcfl),"pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file,arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i = j

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (i == j) {
	setflag[i][i] = 1;
	map[i] = ifuncfl;
	atom->set_mass(i,funcfl[ifuncfl].mass);
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAM::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax,funcfl[m].cut);
  } else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;

  return cutmax;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAM::init_style()
{
  // set communication sizes in comm class

  comm->maxforward_pair = MAX(comm->maxforward_pair,1);
  comm->maxreverse_pair = MAX(comm->maxreverse_pair,1);

  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  cutforcesq = cutmax*cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairEAM::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl-1];

  int me = comm->me;
  FILE *fp;
  char line[MAXLINE];

  if (me == 0) {
    fp = fopen(filename,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(str);
    }
  }

  int tmp;
  if (me == 0) {
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg",&tmp,&file->mass);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %d %lg %lg",
	   &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->mass,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  file->frho = (double *) memory->smalloc((file->nrho+1)*sizeof(double),
					  "pair:frho");
  file->rhor = (double *) memory->smalloc((file->nr+1)*sizeof(double),
					  "pair:rhor");
  file->zr = (double *) memory->smalloc((file->nr+1)*sizeof(double),
					"pair:zr");

  if (me == 0) grab(fp,file->nrho,&file->frho[1]);
  MPI_Bcast(&file->frho[1],file->nrho,MPI_DOUBLE,0,world);

  if (me == 0) grab(fp,file->nr,&file->zr[1]);
  MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) grab(fp,file->nr,&file->rhor[1]);
  MPI_Bcast(&file->rhor[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairEAM::file2array()
{
  int i,j,k,m,n;
  int ntypes = atom->ntypes;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax,rhomax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1) * file->dr);
    rhomax = MAX(rhomax,(file->nrho-1) * file->drho);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array
  
  nfrho = nfuncfl + 1;
  memory->destroy_2d_double_array(frho);
  frho = (double **) memory->create_2d_double_array(nfrho,nrho+1,"pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r,p,cof1,cof2,cof3,cof4;
  
  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = 0.166666667*p*(p*p-1.0);
      frho[n][m] = cof1*file->frho[k-1] + cof2*file->frho[k] + 
	cof3*file->frho[k+1] + cof4*file->frho[k+2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy_2d_double_array(rhor);
  rhor = (double **) memory->create_2d_double_array(nrhor,nr+1,"pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = 0.166666667*p*(p*p-1.0);
      rhor[n][m] = cof1*file->rhor[k-1] + cof2*file->rhor[k] +
	cof3*file->rhor[k+1] + cof4*file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;
  memory->destroy_2d_double_array(z2r);
  z2r = (double **) memory->create_2d_double_array(nz2r,nr+1,"pair:z2r");

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  double zri,zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
	r = (m-1)*dr;

	p = r/ifile->dr + 1.0;
	k = static_cast<int> (p);
	k = MIN(k,ifile->nr-2);
	k = MAX(k,2);
	p -= k;
	p = MIN(p,2.0);
	cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
	cof2 = 0.5*(p*p-1.0)*(p-2.0);
	cof3 = -0.5*p*(p+1.0)*(p-2.0);
	cof4 = 0.166666667*p*(p*p-1.0);
	zri = cof1*ifile->zr[k-1] + cof2*ifile->zr[k] +
	  cof3*ifile->zr[k+1] + cof4*ifile->zr[k+2];

	p = r/jfile->dr + 1.0;
	k = static_cast<int> (p);
	k = MIN(k,jfile->nr-2);
	k = MAX(k,2);
	p -= k;
	p = MIN(p,2.0);
	cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
	cof2 = 0.5*(p*p-1.0)*(p-2.0);
	cof3 = -0.5*p*(p+1.0)*(p-2.0);
	cof4 = 0.166666667*p*(p*p-1.0);
	zrj = cof1*jfile->zr[k-1] + cof2*jfile->zr[k] +
	  cof3*jfile->zr[k+1] + cof4*jfile->zr[k+2];

	z2r[n][m] = 27.2*0.529 * zri*zrj;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2z2r not used

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    irow = map[i];
    if (irow == -1) continue;
    for (j = 1; j <= ntypes; j++) {
      icol = map[j];
      if (icol == -1) continue;
      if (irow < icol) {
	irow = map[j];
	icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairEAM::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy_3d_double_array(frho_spline);
  memory->destroy_3d_double_array(rhor_spline);
  memory->destroy_3d_double_array(z2r_spline);

  frho_spline = memory->create_3d_double_array(nfrho,nrho+1,7,"pair:frho");
  rhor_spline = memory->create_3d_double_array(nrhor,nr+1,7,"pair:rhor");
  z2r_spline = memory->create_3d_double_array(nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];
  
  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) + 
		    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;
  
  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) - 
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] - 
      2.0*(spline[m+1][6]-spline[m][6]);
  }
  
  spline[n][4] = 0.0;
  spline[n][3] = 0.0;
  
  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairEAM::grab(FILE *fp, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fp);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while (ptr = strtok(NULL," \t\n\r\f")) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

void PairEAM::single(int i, int j, int itype, int jtype,
		      double rsq, double factor_coul, double factor_lj,
		      int eflag, One &one)
{
  int m;
  double r,p,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  double *coeff;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);
  
  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp[i]*rhojp + fp[j]*rhoip + phip;
  one.fforce = -psip*recip;

  if (eflag) {
    one.eng_vdwl = phi;
    one.eng_coul = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void PairEAM::single_embed(int i, int itype, double &fpi,
			   int eflag, double &phi)
{
  double p = rho[i]*rdrho + 1.0;
  int m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  
  double *coeff = frho_spline[type2frho[itype]][m];
  fpi = (coeff[0]*p + coeff[1])*p + coeff[2];
  if (eflag) phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
}

/* ---------------------------------------------------------------------- */

int PairEAM::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAM::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairEAM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

int PairEAM::memory_usage()
{
  int bytes = 2 * nmax * sizeof(double);
  return bytes;
}
