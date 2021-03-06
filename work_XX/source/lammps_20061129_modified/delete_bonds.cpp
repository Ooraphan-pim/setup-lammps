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
#include "stdlib.h"
#include "string.h"
#include "delete_bonds.h"
#include "system.h"
#include "atom.h"
#include "domain.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "group.h"
#include "special.h"
#include "error.h"

#define MULTI    1
#define ATOM     2
#define BOND     3
#define ANGLE    4
#define DIHEDRAL 5
#define IMPROPER 6
#define STATS    7

/* ---------------------------------------------------------------------- */

void DeleteBonds::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Delete_bonds command before simulation box is defined");
  if (atom->natoms == 0)
    error->all("Delete_bonds command with no atoms existing");
  if (atom->molecular == 0)
    error->all("Cannot use delete_bonds with non-molecular system");
  if (narg < 2) error->all("Illegal delete_bonds command");

  // init entire system since comm->borders is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for delete_bonds ...\n");
  sys->init();

  if (comm->me == 0 && screen) fprintf(screen,"Deleting bonds ...\n");

  // identify group

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all("Cannot find delete_bonds group ID");
  int groupbit = group->bitmask[igroup];
  
  // set style and which = type value

  int style;
  if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else if (strcmp(arg[1],"atom") == 0) style = ATOM;
  else if (strcmp(arg[1],"bond") == 0) style = BOND;
  else if (strcmp(arg[1],"angle") == 0) style = ANGLE;
  else if (strcmp(arg[1],"dihedral") == 0) style = DIHEDRAL;
  else if (strcmp(arg[1],"improper") == 0) style = IMPROPER;
  else if (strcmp(arg[1],"stats") == 0) style = STATS;
  else error->all("Illegal delete_bonds command");

  int iarg = 2;
  int which;
  if (style != MULTI && style != STATS) {
    if (narg < 3) error->all("Illegal delete_bonds command");
    which = atoi(arg[2]);
    iarg++;
  }

  // grab optional keywords

  int undo_flag = 0;
  int remove_flag = 0;
  int special_flag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"undo") == 0) undo_flag = 1;
    else if (strcmp(arg[iarg],"remove") == 0) remove_flag = 1;
    else if (strcmp(arg[iarg],"special") == 0) special_flag = 1;
    else error->all("Illegal delete_bonds command");
    iarg++;
  }

  // border swap to insure type and mask is current for off-proc atoms
  // enforce PBC before in case atoms are outside box

  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();

  // set topology interactions either off or on
  // criteria for an interaction to potentially be changed (set flag = 1)
  //   all atoms in interaction must be in group
  //   for style = MULTI, no other criteria
  //   for style = ATOM, at least one atom is specified type
  //   for style = BOND/ANGLE/DIHEDRAL/IMPROPER, interaction is specified type
  //   for style = STATS only compute stats, flag is always 0
  // if flag = 1
  //   set interaction type negative if undo_flag = 0
  //   set interaction type positive if undo_flag = 1

  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int i,m,n,flag;
  int atom1,atom2,atom3,atom4;

  if (atom->bonds_allow) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_bond[i]; m++) {
	atom1 = atom->map(atom->bond_atom[i][m]);
	if (atom1 == -1) error->one("Bond atom missing in delete_bonds");
	if (mask[i] & groupbit && mask[atom1] & groupbit) {
	  flag = 0;
	  if (style == MULTI) flag = 1;
	  if (style == ATOM && 
	      (type[i] == which || type[atom1] == which)) flag = 1;
	  if (style == BOND && (bond_type[i][m] == which)) flag = 1;
	  if (flag) {
	    if (undo_flag == 0 && bond_type[i][m] > 0)
	      bond_type[i][m] = -bond_type[i][m]; 
	    if (undo_flag == 1 && bond_type[i][m] < 0)
	      bond_type[i][m] = -bond_type[i][m]; 
	  }
	}
      }
    }
  }

  if (atom->angles_allow) {
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_angle[i]; m++) {
	atom1 = atom->map(atom->angle_atom1[i][m]);
	atom2 = atom->map(atom->angle_atom2[i][m]);
	atom3 = atom->map(atom->angle_atom3[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1)
	  error->one("Angle atom missing in delete_bonds");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit) {
	  flag = 0;
	  if (style == MULTI) flag = 1;
	  if (style == ATOM && 
	      (type[atom1] == which || type[atom2] == which ||
	       type[atom3] == which)) flag = 1;
	  if (style == ANGLE && (angle_type[i][m] == which)) flag = 1;
	  if (flag) {
	    if (undo_flag == 0 && angle_type[i][m] > 0)
	      angle_type[i][m] = -angle_type[i][m]; 
	    if (undo_flag == 1 && angle_type[i][m] < 0)
	      angle_type[i][m] = -angle_type[i][m]; 
	  }
	}
      }
    }
  }

  if (atom->dihedrals_allow) {
    int *num_dihedral = atom->num_dihedral;
    int **dihedral_type = atom->dihedral_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_dihedral[i]; m++) {
	atom1 = atom->map(atom->dihedral_atom1[i][m]);
	atom2 = atom->map(atom->dihedral_atom2[i][m]);
	atom3 = atom->map(atom->dihedral_atom3[i][m]);
	atom4 = atom->map(atom->dihedral_atom4[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
	  error->one("Dihedral atom missing in delete_bonds");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit && mask[atom4] & groupbit) {
	  flag = 0;
	  if (style == MULTI) flag = 1;
	  if (style == ATOM && 
	      (type[atom1] == which || type[atom2] == which ||
	       type[atom3] == which || type[atom4] == which)) flag = 1;
	  if (style == DIHEDRAL && (dihedral_type[i][m] == which)) flag = 1;
	  if (flag) {
	    if (undo_flag == 0 && dihedral_type[i][m] > 0)
	      dihedral_type[i][m] = -dihedral_type[i][m]; 
	    if (undo_flag == 1 && dihedral_type[i][m] < 0)
	      dihedral_type[i][m] = -dihedral_type[i][m]; 
	  }
	}
      }
    }
  }

  if (atom->impropers_allow) {
    int *num_improper = atom->num_improper;
    int **improper_type = atom->improper_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_improper[i]; m++) {
	atom1 = atom->map(atom->improper_atom1[i][m]);
	atom2 = atom->map(atom->improper_atom2[i][m]);
	atom3 = atom->map(atom->improper_atom3[i][m]);
	atom4 = atom->map(atom->improper_atom4[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
	  error->one("Improper atom missing in delete_bonds");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit && mask[atom4] & groupbit) {
	  flag = 0;
	  if (style == MULTI) flag = 1;
	  if (style == ATOM && 
	      (type[atom1] == which || type[atom2] == which ||
	       type[atom3] == which || type[atom4] == which)) flag = 1;
	  if (style == IMPROPER && (improper_type[i][m] == which)) flag = 1;
	  if (flag) {
	    if (undo_flag == 0 && improper_type[i][m] > 0)
	      improper_type[i][m] = -improper_type[i][m]; 
	    if (undo_flag == 1 && improper_type[i][m] < 0)
	      improper_type[i][m] = -improper_type[i][m]; 
	  }
	}
      }
    }
  }

  // remove interactions if requested
  // only if all atoms in bond, angle, etc are in the delete_bonds group

  if (remove_flag) {

    if (atom->bonds_allow) {
      for (i = 0; i < nlocal; i++) {
	m = 0;
	while (m < atom->num_bond[i]) {
	  if (atom->bond_type[i][m] <= 0) {
	    atom1 = atom->map(atom->bond_atom[i][m]);
	    if (mask[i] & groupbit && mask[atom1] & groupbit) {
	      n = atom->num_bond[i];
	      atom->bond_type[i][m] = atom->bond_type[i][n-1];
	      atom->bond_atom[i][m] = atom->bond_atom[i][n-1];
	      atom->num_bond[i]--;
	    } else m++;
	  } else m++;
	}
      }
    }

    if (atom->angles_allow) {
      for (i = 0; i < nlocal; i++) {
	m = 0;
	while (m < atom->num_angle[i]) {
	  if (atom->angle_type[i][m] <= 0) {
	    atom1 = atom->map(atom->angle_atom1[i][m]);
	    atom2 = atom->map(atom->angle_atom2[i][m]);
	    atom3 = atom->map(atom->angle_atom3[i][m]);
	    if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
		mask[atom3] & groupbit) {
	      n = atom->num_angle[i];
	      atom->angle_type[i][m] = atom->angle_type[i][n-1];
	      atom->angle_atom1[i][m] = atom->angle_atom1[i][n-1];
	      atom->angle_atom2[i][m] = atom->angle_atom2[i][n-1];
	      atom->angle_atom3[i][m] = atom->angle_atom3[i][n-1];
	      atom->num_angle[i]--;
	    } else m++;
	  } else m++;
	}
      }
    }

    if (atom->dihedrals_allow) {
      for (i = 0; i < nlocal; i++) {
	m = 0;
	while (m < atom->num_dihedral[i]) {
	  if (atom->dihedral_type[i][m] <= 0) {
	    atom1 = atom->map(atom->dihedral_atom1[i][m]);
	    atom2 = atom->map(atom->dihedral_atom2[i][m]);
	    atom3 = atom->map(atom->dihedral_atom3[i][m]);
	    atom4 = atom->map(atom->dihedral_atom4[i][m]);
	    if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
		mask[atom3] & groupbit && mask[atom4] & groupbit) {
	      n = atom->num_dihedral[i];
	      atom->dihedral_type[i][m] = atom->dihedral_type[i][n-1];
	      atom->dihedral_atom1[i][m] = atom->dihedral_atom1[i][n-1];
	      atom->dihedral_atom2[i][m] = atom->dihedral_atom2[i][n-1];
	      atom->dihedral_atom3[i][m] = atom->dihedral_atom3[i][n-1];
	      atom->dihedral_atom4[i][m] = atom->dihedral_atom4[i][n-1];
	      atom->num_dihedral[i]--;
	    } else m++;
	  } else m++;
	}
      }
    }

    if (atom->impropers_allow) {
      for (i = 0; i < nlocal; i++) {
	m = 0;
	while (m < atom->num_improper[i]) {
	  if (atom->improper_type[i][m] <= 0) {
	    atom1 = atom->map(atom->improper_atom1[i][m]);
	    atom2 = atom->map(atom->improper_atom2[i][m]);
	    atom3 = atom->map(atom->improper_atom3[i][m]);
	    atom4 = atom->map(atom->improper_atom4[i][m]);
	    if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
		mask[atom3] & groupbit && mask[atom4] & groupbit) {
	      n = atom->num_improper[i];
	      atom->improper_type[i][m] = atom->improper_type[i][n-1];
	      atom->improper_atom1[i][m] = atom->improper_atom1[i][n-1];
	      atom->improper_atom2[i][m] = atom->improper_atom2[i][n-1];
	      atom->improper_atom3[i][m] = atom->improper_atom3[i][n-1];
	      atom->improper_atom4[i][m] = atom->improper_atom4[i][n-1];
	      atom->num_improper[i]--;
	    } else m++;
	  } else m++;
	}
      }
    }

  }

  // if interactions were removed, recompute global counts

  if (remove_flag) {

    if (atom->bonds_allow) {
      int nbonds = 0;
      for (i = 0; i < nlocal; i++) nbonds += atom->num_bond[i];
      MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_INT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->nbonds /= 2;
    }

    if (atom->angles_allow) {
      int nangles = 0;
      for (i = 0; i < nlocal; i++) nangles += atom->num_angle[i];
      MPI_Allreduce(&nangles,&atom->nangles,1,MPI_INT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->nangles /= 3;
    }

    if (atom->dihedrals_allow) {
      int ndihedrals = 0;
      for (i = 0; i < nlocal; i++) ndihedrals += atom->num_dihedral[i];
      MPI_Allreduce(&ndihedrals,&atom->ndihedrals,1,MPI_INT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->ndihedrals /= 4;
    }

    if (atom->impropers_allow) {
      int nimpropers = 0;
      for (i = 0; i < nlocal; i++) nimpropers += atom->num_improper[i];
      MPI_Allreduce(&nimpropers,&atom->nimpropers,1,MPI_INT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->nimpropers /= 4;
    }

  }

  // compute and print stats

  int tmp;
  int bond_on,bond_off;
  int angle_on,angle_off;
  int dihedral_on,dihedral_off;
  int improper_on,improper_off;

  if (atom->bonds_allow) {
    bond_on = bond_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++)
	if (atom->bond_type[i][m] > 0) bond_on++;
	else bond_off++;
    MPI_Allreduce(&bond_on,&tmp,1,MPI_INT,MPI_SUM,world);
    bond_on = tmp;
    MPI_Allreduce(&bond_off,&tmp,1,MPI_INT,MPI_SUM,world);
    bond_off = tmp;
    if (force->newton_bond == 0) {
      bond_on /= 2;
      bond_off /= 2;
    }
  }

  if (atom->angles_allow) {
    angle_on = angle_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++)
	if (atom->angle_type[i][m] > 0) angle_on++;
	else angle_off++;
    MPI_Allreduce(&angle_on,&tmp,1,MPI_INT,MPI_SUM,world);
    angle_on = tmp;
    MPI_Allreduce(&angle_off,&tmp,1,MPI_INT,MPI_SUM,world);
    angle_off = tmp;
    if (force->newton_bond == 0) {
      angle_on /= 3;
      angle_off /= 3;
    }
  }

  if (atom->dihedrals_allow) {
    dihedral_on = dihedral_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++)
	if (atom->dihedral_type[i][m] > 0) dihedral_on++;
	else dihedral_off++;
    MPI_Allreduce(&dihedral_on,&tmp,1,MPI_INT,MPI_SUM,world);
    dihedral_on = tmp;
    MPI_Allreduce(&dihedral_off,&tmp,1,MPI_INT,MPI_SUM,world);
    dihedral_off = tmp;
    if (force->newton_bond == 0) {
      dihedral_on /= 4;
      dihedral_off /= 4;
    }
  }

  if (atom->impropers_allow) {
    improper_on = improper_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++)
	if (atom->improper_type[i][m] > 0) improper_on++;
	else improper_off++;
    MPI_Allreduce(&improper_on,&tmp,1,MPI_INT,MPI_SUM,world);
    improper_on = tmp;
    MPI_Allreduce(&improper_off,&tmp,1,MPI_INT,MPI_SUM,world);
    improper_off = tmp;
    if (force->newton_bond == 0) {
      improper_on /= 4;
      improper_off /= 4;
    }
  }

  if (comm->me == 0) {
    if (screen) {
      if (atom->bonds_allow)
	fprintf(screen,"  %d total bonds, %d turned on, %d turned off\n",
		atom->nbonds,bond_on,bond_off);
      if (atom->angles_allow)
	fprintf(screen,"  %d total angles, %d turned on, %d turned off\n",
		atom->nangles,angle_on,angle_off);
      if (atom->dihedrals_allow)
	fprintf(screen,"  %d total dihedrals, %d turned on, %d turned off\n",
		atom->ndihedrals,dihedral_on,dihedral_off);
      if (atom->impropers_allow)
	fprintf(screen,"  %d total impropers, %d turned on, %d turned off\n",
		atom->nimpropers,improper_on,improper_off);
    }
    if (logfile) {
      if (atom->bonds_allow)
	fprintf(logfile,"  %d total bonds, %d turned on, %d turned off\n",
		atom->nbonds,bond_on,bond_off);
      if (atom->angles_allow)
	fprintf(logfile,"  %d total angles, %d turned on, %d turned off\n",
		atom->nangles,angle_on,angle_off);
      if (atom->dihedrals_allow)
	fprintf(logfile,"  %d total dihedrals, %d turned on, %d turned off\n",
		atom->ndihedrals,dihedral_on,dihedral_off);
      if (atom->impropers_allow)
	fprintf(logfile,"  %d total impropers, %d turned on, %d turned off\n",
		atom->nimpropers,improper_on,improper_off);
    }
  }

  // re-compute special list if requested

  if (special_flag) {
    Special special;
    special.build();
  }
}
