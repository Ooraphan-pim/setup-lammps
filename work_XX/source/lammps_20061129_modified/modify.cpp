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

#include "stdio.h"
#include "string.h"
#include "modify.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#define FixInclude
#include "style.h"
#undef FixInclude

#define DELTA 1

// mask settings - same as mask settings in fix.cpp

#define INITIAL_INTEGRATE  1
#define PRE_EXCHANGE       2
#define PRE_NEIGHBOR       4
#define POST_FORCE         8
#define FINAL_INTEGRATE   16
#define END_OF_STEP       32
#define THERMO            64
#define INITIAL_INTEGRATE_RESPA 128
#define POST_FORCE_RESPA        256
#define FINAL_INTEGRATE_RESPA   512
#define MIN_POST_FORCE         1024

enum {NEITHER,PRINT,ENERGY,BOTH};   // same as thermo.cpp

/* ---------------------------------------------------------------------- */

Modify::Modify()
{
  nfix = maxfix = 0;
  n_initial_integrate = 0;
  n_pre_exchange = n_pre_neighbor = 0;
  n_post_force = n_final_integrate = n_end_of_step = n_thermo = 0;
  n_initial_integrate_respa = n_post_force_respa = n_final_integrate_respa = 0;
  n_min_post_force = 0;

  fix = NULL;
  fmask = NULL;
  list_initial_integrate = NULL;
  list_pre_exchange = list_pre_neighbor = NULL;
  list_post_force = list_final_integrate = list_end_of_step = NULL;
  list_thermo = NULL;
  list_initial_integrate_respa = list_post_force_respa = NULL;
  list_final_integrate_respa = NULL;
  list_min_post_force = NULL;

  end_of_step_every = NULL;

  nfix_restart_global = 0;
  id_restart_global = style_restart_global = state_restart_global = NULL;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = NULL;
  index_restart_peratom = NULL;
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes
  // don't invoke atom->update_callback() since atom has already been deleted

  for (int i = 0; i < nfix; i++) delete fix[i];
  memory->sfree(fix);
  memory->sfree(fmask);

  delete [] list_initial_integrate;
  delete [] list_pre_exchange;
  delete [] list_pre_neighbor;
  delete [] list_post_force;
  delete [] list_final_integrate;
  delete [] list_end_of_step;
  delete [] list_thermo;
  delete [] list_initial_integrate_respa;
  delete [] list_post_force_respa;
  delete [] list_final_integrate_respa;
  delete [] list_min_post_force;

  delete [] end_of_step_every;

  restart_deallocate();
}

/* ----------------------------------------------------------------------
   initialize all fixes and lists of fixes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i;

  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate();

  // init each fix

  comm->maxforward_fix = comm->maxreverse_fix = 0;
  for (i = 0; i < nfix; i++) fix[i]->init();

  // create lists of fixes to call at each stage of run

  list_init(INITIAL_INTEGRATE,n_initial_integrate,list_initial_integrate);
  list_init(PRE_EXCHANGE,n_pre_exchange,list_pre_exchange);
  list_init(PRE_NEIGHBOR,n_pre_neighbor,list_pre_neighbor);
  list_init(POST_FORCE,n_post_force,list_post_force);
  list_init(FINAL_INTEGRATE,n_final_integrate,list_final_integrate);
  list_init_end_of_step(END_OF_STEP,n_end_of_step,list_end_of_step);
  list_init_thermo(THERMO,n_thermo,list_thermo);

  list_init(INITIAL_INTEGRATE_RESPA,
	    n_initial_integrate_respa,list_initial_integrate_respa);
  list_init(POST_FORCE_RESPA,
	    n_post_force_respa,list_post_force_respa);
  list_init(FINAL_INTEGRATE_RESPA,
	    n_final_integrate_respa,list_final_integrate_respa);

  list_init(MIN_POST_FORCE,n_min_post_force,list_min_post_force);
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes
------------------------------------------------------------------------- */

void Modify::setup()
{
  if (update->whichflag == 0)
    for (int i = 0; i < nfix; i++) fix[i]->setup();
  else
    for (int i = 0; i < nfix; i++) fix[i]->min_setup();
}

/* ----------------------------------------------------------------------
   1st half of integrate call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate()
{
  for (int i = 0; i < n_initial_integrate; i++)
    fix[list_initial_integrate[i]]->initial_integrate();
}

/* ----------------------------------------------------------------------
   pre_exchange call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_exchange()
{
  for (int i = 0; i < n_pre_exchange; i++)
    fix[list_pre_exchange[i]]->pre_exchange();
}

/* ----------------------------------------------------------------------
   pre_neighbor call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_neighbor()
{
  for (int i = 0; i < n_pre_neighbor; i++)
    fix[list_pre_neighbor[i]]->pre_neighbor();
}

/* ----------------------------------------------------------------------
   force adjustment call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_force(int vflag)
{
  for (int i = 0; i < n_post_force; i++)
    fix[list_post_force[i]]->post_force(vflag);
}

/* ----------------------------------------------------------------------
   2nd half of integrate call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++)
    fix[list_final_integrate[i]]->final_integrate();
}

/* ----------------------------------------------------------------------
   end-of-timestep call only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0)
      fix[list_end_of_step[i]]->end_of_step();
}

/* ----------------------------------------------------------------------
   thermo_fields call only for relevant fixes
   called with n = 0, just query how many fields each fix returns
   called with n > 0, get the field names
------------------------------------------------------------------------- */

int Modify::thermo_fields(int n, int *flags, char **keywords)
{
  int i,j,nfirst,flag_print,flag_energy;

  int m = 0;
  if (n == 0)
    for (i = 0; i < n_thermo; i++)
      m += fix[list_thermo[i]]->thermo_fields(0,NULL,NULL);
  else
    for (i = 0; i < n_thermo; i++) {
      nfirst = m;
      m += fix[list_thermo[i]]->thermo_fields(n,&flags[m],&keywords[m]);
      for (j = nfirst; j < m; j++) {
	if (fix[list_thermo[i]]->thermo_print && 
	    (flags[j] == PRINT || flags[j] == BOTH)) flag_print = PRINT;
	else flag_print = NEITHER;
	if (fix[list_thermo[i]]->thermo_energy && 
	    (flags[j] == ENERGY || flags[j] == BOTH)) flag_energy = ENERGY;
	else flag_energy = NEITHER;
	flags[j] = flag_print + flag_energy;
      }
    }
  return m;
}

/* ----------------------------------------------------------------------
   thermo_compute call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::thermo_compute(double *values)
{
  int m = 0;
  for (int i = 0; i < n_thermo; i++)
    m += fix[list_thermo[i]]->thermo_compute(&values[m]);
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate_respa(int ilevel, int flag)
{
  for (int i = 0; i < n_initial_integrate_respa; i++)
    fix[list_initial_integrate_respa[i]]->initial_integrate_respa(ilevel,flag);
}

/* ----------------------------------------------------------------------
   rRESPA force adjustment call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_post_force_respa; i++)
    fix[list_post_force_respa[i]]->post_force_respa(vflag,ilevel,iloop);
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate_respa(int ilevel)
{
  for (int i = 0; i < n_final_integrate_respa; i++)
    fix[list_final_integrate_respa[i]]->final_integrate_respa(ilevel);
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_force(int vflag)
{
  for (int i = 0; i < n_min_post_force; i++)
    fix[list_min_post_force[i]]->min_post_force(vflag);
}

/* ----------------------------------------------------------------------
   add a new fix or replace one with same ID
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Fix command before simulation box is defined");
  if (narg < 3) error->all("Illegal fix command");

  // find group ID

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find fix group ID");

  // if fix ID exists:
  //   set newflag = 0 so create new fix in same location in fix list
  //   error if new style does not match old style
  //     since can't replace it (all when-to-invoke ptrs would be invalid)
  //   warn if new group != old group
  //   delete old fix
  //   set ptr to NULL in case new fix scans list of fixes
  // if fix ID does not exist:
  //   set newflag = 1 so create new fix
  //   extend fix and fmask lists as necessary

  int ifix,newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;

  if (ifix < nfix) {
    newflag = 0;
    if (strcmp(arg[2],fix[ifix]->style) != 0)
      error->all("Replacing a fix, but new style != old style");
    if (fix[ifix]->igroup != igroup && comm->me == 0)
      error->warning("Replacing a fix, but new group != old group");
    delete fix[ifix];
    atom->update_callback(ifix);
    fix[ifix] = NULL;
  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix,maxfix*sizeof(Fix *),"modify:fix");
      fmask = (int *) 
	memory->srealloc(fmask,maxfix*sizeof(int),"modify:fmask");
    }
  }

  // create the Fix

  if (0) return;         // dummy line to enable else-if macro expansion

#define FixClass
#define FixStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) fix[ifix] = new Class(narg,arg);
#include "style.h"
#undef FixClass

  else error->all("Invalid fix style");

  // if fix is new, set it's mask values and increment nfix

  if (newflag) {
    fmask[ifix] = fix[ifix]->setmask();
    nfix++;
  }

  // check if Fix is in restart_global list
  // if yes, pass state info to the Fix so it can reset itself

  for (int i = 0; i < nfix_restart_global; i++)
    if (strcmp(id_restart_global[i],fix[ifix]->id) == 0 &&
	strcmp(style_restart_global[i],fix[ifix]->style) == 0) {
      fix[ifix]->restart(state_restart_global[i]);
      if (comm->me == 0) {
	char *str = "Resetting global state of Fix %s Style %s "
	  "from restart file info\n";
	if (screen) fprintf(screen,str,fix[ifix]->id,fix[ifix]->style);
	if (logfile) fprintf(logfile,str,fix[ifix]->id,fix[ifix]->style);
      }
    }

  // check if Fix is in restart_peratom list
  // if yes, loop over atoms so they can extract info from atom->extra array

  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strcmp(id_restart_peratom[i],fix[ifix]->id) == 0 &&
	strcmp(style_restart_peratom[i],fix[ifix]->style) == 0) {
      for (int j = 0; j < atom->nlocal; j++)
	fix[ifix]->unpack_restart(j,index_restart_peratom[i]);
      if (comm->me == 0) {
	char *str = "Resetting per-atom state of Fix %s Style %s "
	  "from restart file info\n";
	if (screen) fprintf(screen,str,fix[ifix]->id,fix[ifix]->style);
	if (logfile) fprintf(logfile,str,fix[ifix]->id,fix[ifix]->style);
      }
    }
}

/* ----------------------------------------------------------------------
   modify a Fix's parameters
------------------------------------------------------------------------- */

void Modify::modify_fix(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal fix_modify command");

  // lookup Fix ID

  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;
  if (ifix == nfix) error->all("Could not find fix_modify ID");
  
  fix[ifix]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(char *id)
{
  int ifix;

  // find which fix it is and delete it

  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) error->all("Could not find unfix ID");

  delete fix[ifix];
  if (atom) atom->update_callback(ifix);

  // move other Fixes and fmask down in list one slot

  for (int i = ifix+1; i < nfix; i++) fix[i-1] = fix[i];
  for (int i = ifix+1; i < nfix; i++) fmask[i-1] = fmask[i];
  nfix--;
}

/* ----------------------------------------------------------------------
   write to restart file for all Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
------------------------------------------------------------------------- */

void Modify::write_restart(FILE *fp)
{
  int me = comm->me;

  int count = 0;
  for (int i = 0; i < nfix; i++) 
    if (fix[i]->restart_global) count++;

  if (me == 0) fwrite(&count,sizeof(int),1,fp);

  int n;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_global) {
      if (me == 0) {
	n = strlen(fix[i]->id) + 1;
	fwrite(&n,sizeof(int),1,fp);
	fwrite(fix[i]->id,sizeof(char),n,fp);
	n = strlen(fix[i]->style) + 1;
	fwrite(&n,sizeof(int),1,fp);
	fwrite(fix[i]->style,sizeof(char),n,fp);
      }
      fix[i]->write_restart(fp);
    }

  count = 0;
  for (int i = 0; i < nfix; i++) 
    if (fix[i]->restart_peratom) count++;

  if (me == 0) fwrite(&count,sizeof(int),1,fp);

  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_peratom) {
      if (me == 0) {
	n = strlen(fix[i]->id) + 1;
	fwrite(&n,sizeof(int),1,fp);
	fwrite(fix[i]->id,sizeof(char),n,fp);
	n = strlen(fix[i]->style) + 1;
	fwrite(&n,sizeof(int),1,fp);
	fwrite(fix[i]->style,sizeof(char),n,fp);
	n = fix[i]->maxsize_restart();
	fwrite(&n,sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   read in restart file data on all previously defined Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
   return maxsize of extra info that will be stored with any atom
------------------------------------------------------------------------- */

int Modify::read_restart(FILE *fp)
{
  // nfix_restart_global = # of restart entries with global state info

  int me = comm->me;
  if (me == 0) fread(&nfix_restart_global,sizeof(int),1,fp);
  MPI_Bcast(&nfix_restart_global,1,MPI_INT,0,world);

  // allocate space for each entry

  if (nfix_restart_global) {
    id_restart_global = new char*[nfix_restart_global];
    style_restart_global = new char*[nfix_restart_global];
    state_restart_global = new char*[nfix_restart_global];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, chunk of state data

  int n;
  for (int i = 0; i < nfix_restart_global; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id_restart_global[i] = new char[n];
    if (me == 0) fread(id_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(id_restart_global[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    style_restart_global[i] = new char[n];
    if (me == 0) fread(style_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(style_restart_global[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    state_restart_global[i] = new char[n];
    if (me == 0) fread(state_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(state_restart_global[i],n,MPI_CHAR,0,world);
  }

  // nfix_restart_peratom = # of restart entries with peratom info

  int maxsize = 0;

  if (me == 0) fread(&nfix_restart_peratom,sizeof(int),1,fp);
  MPI_Bcast(&nfix_restart_peratom,1,MPI_INT,0,world);

  // allocate space for each entry

  if (nfix_restart_peratom) {
    id_restart_peratom = new char*[nfix_restart_peratom];
    style_restart_peratom = new char*[nfix_restart_peratom];
    index_restart_peratom = new int[nfix_restart_peratom];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, maxsize of one atom's data
  // set index = which set of extra data this fix represents

  for (int i = 0; i < nfix_restart_peratom; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id_restart_peratom[i] = new char[n];
    if (me == 0) fread(id_restart_peratom[i],sizeof(char),n,fp);
    MPI_Bcast(id_restart_peratom[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    style_restart_peratom[i] = new char[n];
    if (me == 0) fread(style_restart_peratom[i],sizeof(char),n,fp);
    MPI_Bcast(style_restart_peratom[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    maxsize += n;

    index_restart_peratom[i] = i;
  }

  return maxsize;
}

/* ----------------------------------------------------------------------
   delete all lists of restart file Fix info
------------------------------------------------------------------------- */

void Modify::restart_deallocate()
{
  if (nfix_restart_global) {
    for (int i = 0; i < nfix_restart_global; i++) {
      delete [] id_restart_global[i];
      delete [] style_restart_global[i];
      delete [] state_restart_global[i];
    }
    delete [] id_restart_global;
    delete [] style_restart_global;
    delete [] state_restart_global;
  }

  if (nfix_restart_peratom) {
    for (int i = 0; i < nfix_restart_peratom; i++) {
      delete [] id_restart_peratom[i];
      delete [] style_restart_peratom[i];
    }
    delete [] id_restart_peratom;
    delete [] style_restart_peratom;
    delete [] index_restart_peratom;
  }

  nfix_restart_global = nfix_restart_peratom = 0;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for end_of_step fixes
   also create end_of_step_every[]
------------------------------------------------------------------------- */

void Modify::list_init_end_of_step(int mask, int &n, int *&list)
{
  delete [] list;
  delete [] end_of_step_every;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    }
}

/* ----------------------------------------------------------------------
   create list of fix indices for thermo fixes
   also must have thermo print/energy flag set via fix_modify
------------------------------------------------------------------------- */

void Modify::list_init_thermo(int mask, int &n, int *&list)
{
  delete [] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && (fix[i]->thermo_print || fix[i]->thermo_energy))
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && (fix[i]->thermo_print || fix[i]->thermo_energy))
      list[n++] = i;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory from all fixes
------------------------------------------------------------------------- */

int Modify::memory_usage()
{
  int bytes = 0;
  for (int i = 0; i < nfix; i++) bytes += fix[i]->memory_usage();
  return bytes;
}
