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
#include "string.h"
#include "stdlib.h"
#include "sys/types.h"
#include "dirent.h"
#include "read_restart.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "special.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

#define AtomInclude
#include "style.h"
#undef AtomInclude

#define LB_FACTOR 1.1

/* ---------------------------------------------------------------------- */

void ReadRestart::command(int narg, char **arg)
{
  if (domain->box_exist) 
    error->all("Cannot read_restart after simulation box is defined");

  if (narg != 1) error->all("Illegal read_restart command");

  MPI_Comm_rank(world,&me);

  // if filename contains "*", search dir for latest restart file

  char *file = new char[strlen(arg[0]) + 16];
  if (strchr(arg[0],'*')) file_search(arg[0],file);
  else strcpy(file,arg[0]);

  // check if filename contains "%"

  int multiproc;
  if (strchr(file,'%')) multiproc = 1;
  else multiproc = 0;

  // open single restart file or base file for multiproc case

  if (me == 0) {
    if (screen) fprintf(screen,"Reading restart file ...\n");
    char *hfile;
    if (multiproc) {
      char *hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"rb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(str);
    }
    if (multiproc) delete [] hfile;
  }

  // read header info and create atom style and simulation box

  header();
  domain->box_exist = 1;

  // problem setup using info from header

  int n;
  if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);

  atom->allocate_type_arrays();
  atom->grow(n);

  domain->set_initial_box();
  domain->set_global_box();
  comm->set_procs();
  domain->set_local_box();

  // read groups, ntype-length arrays, force field, fix info from file
  // nextra = max # of extra quantities stored with each atom

  group->read_restart(fp);
  if (atom->mass_require) mass();
  if (atom->dipole_require) dipole();
  force_fields();

  int nextra = modify->read_restart(fp);
  atom->nextra_store = nextra;
  atom->extra = memory->create_2d_double_array(n,nextra,"atom:extra");

  // if single file:
  //   proc 0 reads atoms from file, one chunk per proc (nprocs_file)
  // else if one file per proc:
  //   proc 0 reads chunks from series of files (nprocs_file)
  // proc 0 bcasts each chunk to other procs
  // each proc unpacks the atoms, saving ones in it's sub-domain

  int maxbuf = 0;
  double *buf = NULL;

  double subxlo = domain->subxlo;
  double subxhi = domain->subxhi;
  double subylo = domain->subylo;
  double subyhi = domain->subyhi;
  double subzlo = domain->subzlo;
  double subzhi = domain->subzhi;

  int m;
  double xtmp,ytmp,ztmp;
  char *perproc = new char[strlen(file) + 16];
  char *ptr = strchr(file,'%');

  for (int iproc = 0; iproc < nprocs_file; iproc++) {
    if (me == 0) {
      if (multiproc) {
	fclose(fp);
	*ptr = '\0';
	sprintf(perproc,"%s%d%s",file,iproc,ptr+1);
	*ptr = '%';
	fp = fopen(perproc,"rb");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open restart file %s",perproc);
	  error->one(str);
	}
      }
      fread(&n,sizeof(int),1,fp);
    }

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n > maxbuf) {
      maxbuf = n;
      delete [] buf;
      buf = new double[maxbuf];
    }

    if (me == 0) fread(buf,sizeof(double),n,fp);
    MPI_Bcast(buf,n,MPI_DOUBLE,0,world);

    m = 0;
    while (m < n) {
      xtmp = buf[m+1];
      ytmp = buf[m+2];
      ztmp = buf[m+3];
      if (xtmp >= subxlo && xtmp < subxhi &&
	  ytmp >= subylo && ytmp < subyhi &&
	  ztmp >= subzlo && ztmp < subzhi)
	m += atom->unpack_restart(&buf[m]);
      else m += static_cast<int> (buf[m]);
    }
  }

  // close restart file and clean-up memory
  
  if (me == 0) fclose(fp);
  delete [] buf;
  delete [] file;
  delete [] perproc;

  // check that all atoms were assigned to procs

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %.15g atoms\n",natoms);
    if (logfile) fprintf(logfile,"  %.15g atoms\n",natoms);
  }

  if (natoms != atom->natoms) error->all("Did not assign all atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  %d bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  %d bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  %d angles\n",atom->nangles);
      if (logfile) fprintf(logfile,"  %d angles\n",atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  %d dihedrals\n",atom->ndihedrals);
      if (logfile) fprintf(logfile,"  %d dihedrals\n",atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  %d impropers\n",atom->nimpropers);
      if (logfile) fprintf(logfile,"  %d impropers\n",atom->nimpropers);
    }
  }

  // check if tags are being used
  // create global mapping and bond topology now that system is defined

  int flag = 0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->tag[i] > 0) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,world);
  if (flag_all == 0) atom->tag_enable = 0;

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
  if (atom->molecular) {
    Special special;
    special.build();
  }
}

/* ----------------------------------------------------------------------
   search for all files in dir matching infile which contains "*"
   replace "*" with latest timestep value to create outfile name
   if infile also contains "%", need to use "base" when search directory
------------------------------------------------------------------------- */

void ReadRestart::file_search(char *infile, char *outfile)
{
  // if filename contains "%" replace "%" with "base"

  int multiproc;
  char *pattern = new char[strlen(infile) + 16];
  char *ptr;

  if (ptr = strchr(infile,'%')) {
    multiproc = 1;
    *ptr = '\0';
    sprintf(pattern,"%s%s%s",infile,"base",ptr+1);
    *ptr = '%';
  } else {
    multiproc = 0;
    strcpy(pattern,infile);
  }

  // scan all files in directory, searching for files that match pattern
  // maxnum = largest int that matches "*"

  int n = strlen(pattern) + 16;
  char *begin = new char[n];
  char *middle = new char[n];
  char *end = new char[n];

  ptr = strchr(pattern,'*');
  *ptr = '\0';
  strcpy(begin,pattern);
  strcpy(end,ptr+1);
  int nbegin = strlen(begin);
  int maxnum = -1;

  if (me == 0) {
    struct dirent *ep;
    DIR *dp = opendir("./");
    if (dp == NULL) 
      error->one("Cannot open dir to search for restart file");
    while (ep = readdir(dp)) {
      if (strstr(ep->d_name,begin) != ep->d_name) continue;
      if ((ptr = strstr(&ep->d_name[nbegin],end)) == NULL) continue;
      if (strlen(end) == 0) ptr = ep->d_name + strlen(ep->d_name);
      *ptr = '\0';
      if (strlen(&ep->d_name[nbegin]) < n) {
	strcpy(middle,&ep->d_name[nbegin]);
	if (atoi(middle) > maxnum) maxnum = atoi(middle);
      }
    }
    closedir(dp);
    if (maxnum < 0) error->one("Found no restart file matching pattern");
  }

  delete [] pattern;
  delete [] begin;
  delete [] middle;
  delete [] end;

  // create outfile with maxint substituted for "*"
  // use original infile, not pattern, since need to retain "%" in filename

  ptr = strchr(infile,'*');
  *ptr = '\0';
  sprintf(outfile,"%s%d%s",infile,maxnum,ptr+1);
  *ptr = '*';
}

/* ----------------------------------------------------------------------
   read header of restart file 
------------------------------------------------------------------------- */

void ReadRestart::header()
{
  int px,py,pz;
  int xperiodic,yperiodic,zperiodic;
  int boundary[3][2];

  // read flags and values until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    // check restart file version, compare to LAMMPS version

    if (flag == 0) {
      char *version = read_char();
      if (strcmp(version,universe->version) != 0) {
	error->warning("Restart file version does not match LAMMPS version");
	if (screen) fprintf(screen,"  restart file = %s, LAMMPS = %s\n",
			    version,universe->version);
      }
      delete [] version;

      // reset unit_style only if different
      // so that timestep,neighbor-skin are not changed

    } else if (flag == 1) {
      char *style = read_char();
      if (strcmp(style,update->unit_style) != 0 && me == 0)
	error->warning("Resetting unit_style to restart file value");
      if (strcmp(style,update->unit_style) != 0) update->set_units(style);
      delete [] style;

    } else if (flag == 2) {
      update->ntimestep = read_int();

      // set dimension from restart file, warn if different

    } else if (flag == 3) {
      int dimension = read_int();
      if (dimension != force->dimension && me == 0)
	error->warning("Resetting dimension to restart file value");
      force->dimension = dimension;
      if (force->dimension == 2 && domain->zperiodic == 0)
	error->all("Cannot run 2d simulation with nonperiodic Z dimension");

      // read nprocs_file from restart file, warn if different

    } else if (flag == 4) {
      nprocs_file = read_int();
      if (nprocs_file != comm->nprocs && me == 0)
	error->warning("Restart file used different # of processors");

      // don't set procgrid, warn if different

    } else if (flag == 5) {
      px = read_int();
    } else if (flag == 6) {
      py = read_int();
    } else if (flag == 7) {
      pz = read_int();
      if (comm->user_procgrid[0] != 0 && 
	  (px != comm->user_procgrid[0] || py != comm->user_procgrid[1] || 
	   pz != comm->user_procgrid[2]) && me == 0)
	error->warning("Restart file used different 3d processor grid");

      // don't set newton_pair, warn if different
      // set newton_bond from restart file, warn if different

    } else if (flag == 8) {
      int newton_pair = read_int();
      if (newton_pair != force->newton_pair && me == 0)
	error->warning("Restart file used different newton pair setting");
    } else if (flag == 9) {
      int newton_bond = read_int();
      if (newton_bond != force->newton_bond && me == 0)
	error->warning("Resetting newton bond to restart file value");
      force->newton_bond = newton_bond;
      if (force->newton_pair || force->newton_bond) force->newton = 1;
      else force->newton = 0;

      // set boundary settings from restart file, warn if different

    } else if (flag == 10) {
      xperiodic = read_int();
    } else if (flag == 11) {
      yperiodic = read_int();
    } else if (flag == 12) {
      zperiodic = read_int();
    } else if (flag == 13) {
      boundary[0][0] = read_int();
    } else if (flag == 14) {
      boundary[0][1] = read_int();
    } else if (flag == 15) {
      boundary[1][0] = read_int();
    } else if (flag == 16) {
      boundary[1][1] = read_int();
    } else if (flag == 17) {
      boundary[2][0] = read_int();
    } else if (flag == 18) {
      boundary[2][1] = read_int();

      int flag = 0;
      if ((xperiodic != domain->xperiodic || yperiodic != domain->yperiodic ||
	   zperiodic != domain->zperiodic)) flag = 1;
      if (boundary[0][0] != domain->boundary[0][0] || 
	  boundary[0][1] != domain->boundary[0][1] ||
	  boundary[1][0] != domain->boundary[1][0] || 
	  boundary[1][1] != domain->boundary[1][1] ||
	  boundary[2][0] != domain->boundary[2][0] || 
	  boundary[2][1] != domain->boundary[2][1]) flag = 1;
      if (flag && me == 0) 
	error->warning("Resetting boundary settings to restart file values");

      domain->xperiodic = xperiodic;
      domain->yperiodic = yperiodic;
      domain->zperiodic = zperiodic;
      domain->boundary[0][0] = boundary[0][0];
      domain->boundary[0][1] = boundary[0][1];
      domain->boundary[1][0] = boundary[1][0];
      domain->boundary[1][1] = boundary[1][1];
      domain->boundary[2][0] = boundary[2][0];
      domain->boundary[2][1] = boundary[2][1];
  
      domain->nonperiodic = 0;
      if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
	domain->nonperiodic = 1;
	if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
	    boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
	    boundary[2][0] >= 2 || boundary[2][1] >= 2)
	  domain->nonperiodic = 2;
      }

      // create atom class
      // if style = hybrid, read additional sub-class arguments

    } else if (flag == 19) {
      char *style = read_char();

      int nwords,n;
      char **words;

      if (strcmp(style,"hybrid") != 0) {
	nwords = 1;
	words = new char*[nwords];
	words[0] = style;
      } else {
	if (me == 0) fread(&nwords,sizeof(int),1,fp);
	MPI_Bcast(&nwords,1,MPI_INT,0,world);
	nwords++;
	words = new char*[nwords];
	words[0] = style;
	for (int i = 1; i < nwords; i++) {
	  if (me == 0) fread(&n,sizeof(int),1,fp);
	  MPI_Bcast(&n,1,MPI_INT,0,world);
	  words[i] = new char[n];
	  if (me == 0) fread(words[i],sizeof(char),n,fp);
	  MPI_Bcast(words[i],n,MPI_CHAR,0,world);
	}
      }

      Atom *old = atom;
	
      if (0) return;         // dummy line to enable else-if macro expansion

#define AtomClass
#define AtomStyle(key,Class) \
      else if (strcmp(style,#key) == 0) atom = new Class(nwords,words);
#include "style.h"
#undef AtomClass

      else error->all("Unknown atom style in restart file");

      if (old) {
	atom->settings(old);
	delete old;
      }

      for (int i = 0; i < nwords; i++) delete [] words[i];
      delete [] words;

    } else if (flag == 20) {
      atom->natoms = read_double();
    } else if (flag == 21) {
      atom->ntypes = read_int();
    } else if (flag == 22) {
      atom->nbonds = read_int();
    } else if (flag == 23) {
      atom->nbondtypes = read_int();
    } else if (flag == 24) {
      atom->bond_per_atom = read_int();
    } else if (flag == 25) {
      atom->nangles = read_int();
    } else if (flag == 26) {
      atom->nangletypes = read_int();
    } else if (flag == 27) {
      atom->angle_per_atom = read_int();
    } else if (flag == 28) {
      atom->ndihedrals = read_int();
    } else if (flag == 29) {
      atom->ndihedraltypes = read_int();
    } else if (flag == 30) {
      atom->dihedral_per_atom = read_int();
    } else if (flag == 31) {
      atom->nimpropers = read_int();
    } else if (flag == 32) {
      atom->nimpropertypes = read_int();
    } else if (flag == 33) {
      atom->improper_per_atom = read_int();

    } else if (flag == 34) {
      domain->boxxlo = read_double();
    } else if (flag == 35) {
      domain->boxxhi = read_double();
    } else if (flag == 36) {
      domain->boxylo = read_double();
    } else if (flag == 37) {
      domain->boxyhi = read_double();
    } else if (flag == 38) {
      domain->boxzlo = read_double();
    } else if (flag == 39) {
      domain->boxzhi = read_double();

    } else if (flag == 40) {
      force->special_lj[1] = read_double();
    } else if (flag == 41) {
      force->special_lj[2] = read_double();
    } else if (flag == 42) {
      force->special_lj[3] = read_double();
    } else if (flag == 43) {
      force->special_coul[1] = read_double();
    } else if (flag == 44) {
      force->special_coul[2] = read_double();
    } else if (flag == 45) {
      force->special_coul[3] = read_double();

    } else error->all("Invalid flag in header of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadRestart::mass()
{
  double *mass = new double[atom->ntypes+1];
  if (me == 0) fread(&mass[1],sizeof(double),atom->ntypes,fp);
  MPI_Bcast(&mass[1],atom->ntypes,MPI_DOUBLE,0,world);
  atom->set_mass(mass);
  delete [] mass;
}

/* ---------------------------------------------------------------------- */

void ReadRestart::dipole()
{
  double *dipole = new double[atom->ntypes+1];
  if (me == 0) fread(&dipole[1],sizeof(double),atom->ntypes,fp);
  MPI_Bcast(&dipole[1],atom->ntypes,MPI_DOUBLE,0,world);
  atom->set_dipole(dipole);
  delete [] dipole;
}

/* ---------------------------------------------------------------------- */

void ReadRestart::force_fields()
{
  int n;
  char *style;

  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n) {
    style = (char *) memory->smalloc(n*sizeof(char),"read_restart:style");
    if (me == 0) fread(style,sizeof(char),n,fp);
    MPI_Bcast(style,n,MPI_CHAR,0,world);

    if (force->pair == NULL || strcmp(style,force->pair_style)) {
      if (force->pair) {
	if (me == 0)
	  error->warning("Resetting pair_style to restart file value");
	delete force->pair;
      }
      force->create_pair(style);
    }
      
    memory->sfree(style);
    force->pair->read_restart(fp);
  } else force->create_pair("none");

  if (atom->bonds_allow) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n) {
      style = (char *) memory->smalloc(n*sizeof(char),"read_restart:style");
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      if (force->bond == NULL || strcmp(style,force->bond_style)) {
	if (force->bond) {
	  if (me == 0)
	    error->warning("Resetting bond_style to restart file value");
	  delete force->bond;
	}
	force->create_bond(style);
      }
      
      memory->sfree(style);
      force->bond->read_restart(fp);
    } else force->create_bond("none");
  }

  if (atom->angles_allow) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n) {
      style = (char *) memory->smalloc(n*sizeof(char),"read_restart:style");
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      if (force->angle == NULL || strcmp(style,force->angle_style)) {
	if (force->angle) {
	  if (me == 0)
	    error->warning("Resetting angle_style to restart file value");
	  delete force->angle;
	}
	force->create_angle(style);
      }

      memory->sfree(style);
      force->angle->read_restart(fp);
    } else force->create_angle("none");
  }

  if (atom->dihedrals_allow) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n) {
      style = (char *) memory->smalloc(n*sizeof(char),"read_restart:style");
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      if (force->dihedral == NULL || strcmp(style,force->dihedral_style)) {
	if (force->dihedral) {
	  if (me == 0)
	    error->warning("Resetting dihedral_style to restart file value");
	  delete force->dihedral;
	}
	force->create_dihedral(style);
      }

      memory->sfree(style);
      force->dihedral->read_restart(fp);
    } else force->create_dihedral("none");
  }

  if (atom->impropers_allow) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n) {
      style = (char *) memory->smalloc(n*sizeof(char),"read_restart:style");
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      if (force->improper == NULL || strcmp(style,force->improper_style)) {
	if (force->improper) {
	  if (me == 0)
	    error->warning("Resetting improper_style to restart file value");
	  delete force->improper;
	}
	force->create_improper(style);
      }

      memory->sfree(style);
      force->improper->read_restart(fp);
    } else force->create_improper("none");
  }
}

/* ----------------------------------------------------------------------
   read an int from restart file 
------------------------------------------------------------------------- */

int ReadRestart::read_int()
{
  int value;
  if (me == 0) fread(&value,sizeof(int),1,fp);
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file 
------------------------------------------------------------------------- */

double ReadRestart::read_double()
{
  double value;
  if (me == 0) fread(&value,sizeof(double),1,fp);
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char str from restart file 
------------------------------------------------------------------------- */

char *ReadRestart::read_char()
{
  int n;
  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  char *value = new char[n];
  if (me == 0) fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}
