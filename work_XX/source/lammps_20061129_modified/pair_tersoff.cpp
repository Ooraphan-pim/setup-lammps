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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tersoff.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "memory.h"
#include "neighbor.h"
#include "memory.h"
#include "error.h"

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairTersoff::PairTersoff()
{
  neigh_half_every = 0;
  neigh_full_every = 1;
  single_enable = 0;
  one_coeff = 1;

  PI2 = 2.0*atan(1.0);
  PI4 = atan(1.0);

  nelements = 0;
  elements = NULL;
  nparams = 0;
  maxparam = 0;
  params = NULL;
  elem2param = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoff::~PairTersoff()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  delete [] params;
  memory->destroy_3d_int_array(elem2param);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairTersoff::compute(int eflag, int vflag)
{
  int i,j,k,m,n,itag,jtag,itype,jtype,ktype,iparam_ij,iparam_ijk,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,rsq1,rsq2,eng,fforce;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *neighs;
  double **f;

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // loop over full neighbor list of my atoms

  for (i = 0; i < nlocal; i++) {
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // two-body interactions, skip half of them

    neighs = neighbor->firstneigh_full[i];
    numneigh = neighbor->numneigh_full[i];

    for (m = 0; m < numneigh; m++) {
      j = neighs[m];
      jtag = tag[j];

      if (itag > jtag) {
	if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
	if ((itag+jtag) % 2 == 1) continue;
      } else {
	if (x[j][2] < x[i][2]) continue;
	else if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	else if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp)
	    continue;
      }

      jtype = map[type[j]];

      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem2param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fforce,eflag,eng);

      if (eflag) eng_vdwl += eng;
      f[i][0] += fforce*delx;
      f[i][1] += fforce*dely;
      f[i][2] += fforce*delz;
      f[j][0] -= fforce*delx;
      f[j][1] -= fforce*dely;
      f[j][2] -= fforce*delz;
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff

    for (m = 0; m < numneigh; m++) {
      j = neighs[m];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;
      for (n = 0; n < numneigh; n++) {
	if (n == m) continue;
	k = neighs[n];
	ktype = map[type[k]];
	iparam_ijk = elem2param[itype][jtype][ktype];

	delr2[0] = x[k][0] - xtmp;
	delr2[1] = x[k][1] - ytmp;
	delr2[2] = x[k][2] - ztmp;
	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	if (rsq2 > params[iparam_ijk].cutsq) continue;

	zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fforce,prefactor,eflag,eng);

      if (eflag) eng_vdwl += eng;
      f[i][0] += fforce*delr1[0];
      f[i][1] += fforce*delr1[1];
      f[i][2] += fforce*delr1[2];
      f[j][0] -= fforce*delr1[0];
      f[j][1] -= fforce*delr1[1];
      f[j][2] -= fforce*delr1[2];

      // attractive term via loop over k

      for (n = 0; n < numneigh; n++) {
	if (n == m) continue;
	k = neighs[n];
	ktype = map[type[k]];
	iparam_ijk = elem2param[itype][jtype][ktype];

	delr2[0] = x[k][0] - xtmp;
	delr2[1] = x[k][1] - ytmp;
	delr2[2] = x[k][2] - ztmp;
	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	if (rsq2 > params[iparam_ijk].cutsq) continue;

	attractive(&params[iparam_ijk],prefactor,
		   rsq1,rsq2,delr1,delr2,fi,fj,fk);

	f[i][0] += fi[0];
	f[i][1] += fi[1];
	f[i][2] += fi[2];
	f[j][0] += fj[0];
	f[j][1] += fj[1];
	f[j][2] += fj[2];
	f[k][0] += fk[0];
	f[k][1] += fk[1];
	f[k][2] += fk[2];
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairTersoff::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairTersoff::settings(int narg, char **arg)
{
  if (narg != 0) error->all("Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTersoff::coeff(int narg, char **arg)
{
  int i,j,m,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all("Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters
  
  read_file(arg[2]);
  setup();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoff::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::init_style()
{
  if (atom->tag_enable == 0)
    error->all("Pair style Tersoff requires atom IDs");
  if (force->newton_pair == 0)
    error->all("Pair style Tersoff requires newton pair on");
}

/* ---------------------------------------------------------------------- */

void PairTersoff::read_file(char *file)
{
  int params_per_line = 15;
  char **words = new char*[params_per_line+1];

  if (params) delete [] params;
  params = NULL;
  nparams = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Tersoff potential file %s",file);
      error->one(str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int i,n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
	eof = 1;
	fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if (ptr = strchr(line,'#')) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
	  eof = 1;
	  fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if (ptr = strchr(line,'#')) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all("Incorrect format in Tersoff potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
					  "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].lam3 = atof(words[3]);
    params[nparams].c = atof(words[4]);
    params[nparams].d = atof(words[5]);
    params[nparams].h = atof(words[6]);
    params[nparams].powern = atof(words[7]);
    params[nparams].beta = atof(words[8]);
    params[nparams].lam2 = atof(words[9]);
    params[nparams].bigb = atof(words[10]);
    params[nparams].bigr = atof(words[11]);
    params[nparams].bigd = atof(words[12]);
    params[nparams].lam1 = atof(words[13]);
    params[nparams].biga = atof(words[14]);

    if (params[nparams].lam3 < 0.0 || params[nparams].c <= 0.0 || 
	params[nparams].d <= 0.0 || params[nparams].powern <= 0.0 || 
	params[nparams].beta < 0.0 || params[nparams].lam2 < 0.0 || 
	params[nparams].bigb < 0.0 || params[nparams].bigr < 0.0 ||
	params[nparams].bigd < 0.0 ||
	params[nparams].bigd >= params[nparams].bigr ||
	params[nparams].lam3 < 0.0 || params[nparams].biga < 0.0)
      error->all("Illegal Tersoff parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::setup()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  if (elem2param) memory->destroy_3d_int_array(elem2param);
  elem2param = memory->create_3d_int_array(nelements,nelements,nelements,
					   "pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
	n = -1;
	for (m = 0; m < nparams; m++) {
	  if (i == params[m].ielement && j == params[m].jelement && 
	      k == params[m].kelement) {
	    if (n >= 0) error->all("Potential file has duplicate entry");
	    n = m;
	  }
	}
	if (n < 0) error->all("Potential file is missing an entry");
	elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}  

/* ---------------------------------------------------------------------- */

void PairTersoff::repulsive(Param *param, double rsq, double &fforce,
			    int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;

  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::zeta(Param *param, double rsqij, double rsqik,
			 double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] + 
	      delrij[2]*delrik[2]) / (rij*rik);

  arg = pow(param->lam3 * (rij-rik),3);
  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);
  
  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::force_zeta(Param *param, double rsq, double zeta_ij,
			     double &fforce, double &prefactor,
			     int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairTersoff::attractive(Param *param, double prefactor,
			     double rsqij, double rsqik,
			     double *delrij, double *delrik,
			     double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);
  
  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  
  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  
  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(PI4/ters_D) * cos(PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}   

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) * 
			  (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) * 
			   pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);
			  
  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_gijk(double costheta, Param *param)
{
  double ters_c = param->c;
  double ters_d = param->d;

  return 1.0 + pow(ters_c/ters_d,2.0) -
    pow(ters_c,2.0) / (pow(ters_d,2.0) + pow(param->h - costheta,2.0));
};

/* ---------------------------------------------------------------------- */

double PairTersoff::ters_gijk_d(double costheta, Param *param)
{
  double numerator = -2.0 * pow(param->c,2) * (param->h - costheta);
  double denominator = pow(pow(param->d,2.0) + 
			   pow(param->h - costheta,2.0),2.0);
  return numerator/denominator;
}

/* ---------------------------------------------------------------------- */

void PairTersoff::ters_zetaterm_d(double prefactor,
				  double *rij_hat, double rij,
				  double *rik_hat, double rik,
				  double *dri, double *drj, double *drk,
				  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  tmp = pow(param->lam3 * (rij-rik),3.0);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  ex_delr_d = (3.0*pow(param->lam3,3.0)) * pow(rij-rik,2.0)*ex_delr;
  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairTersoff::costheta_d(double *rij_hat, double rij,
			     double *rik_hat, double rik,
			     double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
};

/* ----------------------------------------------------------------------
   vector routines
------------------------------------------------------------------------- */

/*
double PairTersoff::vec3_dot(double *x, double *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void PairTersoff::vec3_add(double *x, double *y, double *z)
{
  z[0] = x[0]+y[0];
  z[1] = x[1]+y[1];
  z[2] = x[2]+y[2];
}

void PairTersoff::vec3_scale(double k, double *x, double *y)
{
  y[0] = k*x[0];
  y[1] = k*x[1];
  y[2] = k*x[2];
}

void PairTersoff::vec3_scaleadd(double k, double *x, double *y, double *z)
{
  z[0] = k*x[0]+y[0];
  z[1] = k*x[1]+y[1];
  z[2] = k*x[2]+y[2];
}

*/
