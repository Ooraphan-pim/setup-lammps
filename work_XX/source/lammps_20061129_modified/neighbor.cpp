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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "limits.h"
#include "neighbor.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "update.h"
#include "respa.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#define PGDELTA 1
#define LB_FACTOR 1.5
#define SMALL 1.0e-6
#define EXDELTA 1

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor()
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  style = 1;
  every = 1;
  delay = 10;
  dist_check = 1;
  pgsize = 10000;
  oneatom = 2000;

  maxlocal = 0;
    
  cutneighsq = NULL;
  fixchecklist = NULL;

  // last neighbor info

  maxhold = 0;
  xhold = NULL;

  // pair exclusion list info

  nex_type = maxex_type = 0;
  ex1_type = ex2_type = NULL;
  ex_type = NULL;

  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = NULL;

  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = NULL;

  // bin info

  maxhead = 0;
  binhead = NULL;
  maxbin = 0;
  bins = NULL;

  maxstencil = 0;
  stencil = NULL;
  maxstencil_full = 0;
  stencil_full = NULL;

  // half neighbor list info

  half = half_command = 0;
  maxpage = 0;
  numneigh = NULL;
  firstneigh = NULL;
  pages = NULL;

  // full neighbor list info

  full = 0;
  maxpage_full = 0;
  numneigh_full = NULL;
  firstneigh_full = NULL;
  pages_full = NULL;

  // shear history neighbor list info

  history = -1;
  firsttouch = NULL;
  firstshear = NULL;
  pages_touch = NULL;
  pages_shear = NULL;

  // multiple respa neighbor list info

  respa = 0;
  maxpage_inner = 0;
  maxpage_middle = 0;
  numneigh_inner = NULL;
  firstneigh_inner = NULL;
  pages_inner = NULL;
  numneigh_middle = NULL;
  firstneigh_middle = NULL;
  pages_middle = NULL;

  // bond list info

  maxbond = 0;
  bondlist = NULL;
  maxangle = 0;
  anglelist = NULL;
  maxdihedral = 0;
  dihedrallist = NULL;
  maximproper = 0;
  improperlist = NULL;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{
  memory->destroy_2d_double_array(cutneighsq);
  delete [] fixchecklist;
  memory->destroy_2d_double_array(xhold);

  memory->sfree(ex1_type);
  memory->sfree(ex2_type);
  memory->destroy_2d_int_array(ex_type);

  memory->sfree(ex1_group);
  memory->sfree(ex2_group);
  delete [] ex1_bit;
  delete [] ex2_bit;

  memory->sfree(ex_mol_group);
  delete [] ex_mol_bit;

  memory->sfree(binhead);
  memory->sfree(bins);
  memory->sfree(stencil);
  memory->sfree(stencil_full);

  memory->destroy_2d_int_array(bondlist);
  memory->destroy_2d_int_array(anglelist);
  memory->destroy_2d_int_array(dihedrallist);
  memory->destroy_2d_int_array(improperlist);

  memory->sfree(numneigh);
  memory->sfree(firstneigh);
  for (int i = 0; i < maxpage; i++) memory->sfree(pages[i]);
  memory->sfree(pages);

  memory->sfree(numneigh_full);
  memory->sfree(firstneigh_full);
  for (int i = 0; i < maxpage_full; i++) memory->sfree(pages_full[i]);
  memory->sfree(pages_full);

  if (history >= 0) {
    memory->sfree(firsttouch);
    memory->sfree(firstshear);
    for (int i = 0; i < maxpage; i++) memory->sfree(pages_touch[i]);
    for (int i = 0; i < maxpage; i++) memory->sfree(pages_shear[i]);
    memory->sfree(pages_touch);
    memory->sfree(pages_shear);
  }

  if (respa) {
    memory->sfree(numneigh_inner);
    memory->sfree(firstneigh_inner);
    for (int i = 0; i < maxpage_inner; i++) memory->sfree(pages_inner[i]);
    memory->sfree(pages_inner);
    memory->sfree(numneigh_middle);
    memory->sfree(firstneigh_middle);
    for (int i = 0; i < maxpage_middle; i++) memory->sfree(pages_middle[i]);
    memory->sfree(pages_middle);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,m,n;

  ncalls = ndanger = 0;

  // error check

  if (delay > 0 && (delay % every) != 0)
    error->all("Neighbor delay must be 0 or multiple of every setting");

  // ------------------------------------------------------------------
  // settings

  // set cutneigh and trigger distance for reneighboring

  if (force->pair) cutneigh = force->pair->cutforce + skin;
  else cutneigh = skin;
  triggersq = 0.25*skin*skin;
  if (cutneighsq == NULL)
    cutneighsq = memory->create_2d_double_array(atom->ntypes+1,atom->ntypes+1,
						"neigh:cutneighsq");

  // set neighbor cutoffs with skin included
  // if no pair defined, cutneigh is just skin

  n = atom->ntypes;

  if (force->pair) {
    double cutoff;
    double **cutsq = force->pair->cutsq;
    for (i = 1; i <= n; i++)
      for (j = i; j <= n; j++) {
	cutoff = sqrt(cutsq[i][j]);
	cutneighsq[i][j] = (cutoff+skin) * (cutoff+skin);
	cutneighsq[j][i] = cutneighsq[i][j];
      }
  } else {
    for (i = 1; i <= n; i++)
      for (j = i; j <= n; j++) {
	cutneighsq[i][j] = skin*skin;
	cutneighsq[j][i] = cutneighsq[i][j];
      }
  }

  // check other classes that can induce reneighboring in decide()

  restart_check = 0;
  if (output->restart_every) restart_check = 1;

  delete [] fixchecklist;
  fixchecklist = NULL;
  fixchecklist = new int[modify->nfix];

  fix_check = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;

  must_check = 0;
  if (restart_check || fix_check) must_check = 1;

  // set special_flag for 1-2, 1-3, 1-4 neighbors
  // flag[0] is not used, flag[1] = 1-2, flag[2] = 1-3, flag[3] = 1-4
  // flag = 0 if both LJ/Coulomb special values are 0.0
  // flag = 1 if both LJ/Coulomb special values are 1.0
  // flag = 2 otherwise or if KSpace solver is enabled
  // (pairwise portion of KSpace solver uses all 1-2,1-3,1-4 neighbors)

  if (force->special_lj[1] == 0.0 && force->special_coul[1] == 0.0) 
    special_flag[1] = 0;
  else if (force->special_lj[1] == 1.0 && force->special_coul[1] == 1.0) 
    special_flag[1] = 1;
  else special_flag[1] = 2;

  if (force->special_lj[2] == 0.0 && force->special_coul[2] == 0.0) 
    special_flag[2] = 0;
  else if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0) 
    special_flag[2] = 1;
  else special_flag[2] = 2;

  if (force->special_lj[3] == 0.0 && force->special_coul[3] == 0.0) 
    special_flag[3] = 0;
  else if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) 
    special_flag[3] = 1;
  else special_flag[3] = 2;

  if (force->kspace) special_flag[1] = special_flag[2] = special_flag[3] = 2;

  // ------------------------------------------------------------------
  // memory management

  // free xhold and bins if not needed for this run

  if (dist_check == 0) {
    memory->destroy_2d_double_array(xhold);
    maxhold = 0;
  }

  if (style == 0) {
    memory->sfree(bins);
    memory->sfree(binhead);
    memory->sfree(stencil);
    memory->sfree(stencil_full);
    maxbin = maxhead = maxstencil = maxstencil_full = 0;
  }

  // 1st time allocation of xhold and bins

  if (dist_check) {
    if (maxhold == 0) {
      maxhold = atom->nmax;
      xhold = memory->create_2d_double_array(maxhold,3,"neigh:xhold");
    }
  }

  if (style == 1) {
    if (maxbin == 0) {
      maxbin = atom->nmax;
      bins = (int *) memory->smalloc(maxbin*sizeof(int),"bins");
    }
  }
    
  // exclusion lists for type, group, molecule settings from neigh_modify

  n = atom->ntypes;

  if (nex_type == 0 && nex_group == 0 && nex_mol == 0) exclude = 0;
  else exclude = 1;

  if (nex_type) {
    memory->destroy_2d_int_array(ex_type);
    ex_type = (int **) memory->create_2d_int_array(n+1,n+1,"neigh:ex_type");

    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
	ex_type[i][j] = 0;

    for (i = 0; i < nex_type; i++) {
      if (ex1_type[i] <= 0 || ex1_type[i] > n || 
	  ex2_type[i] <= 0 || ex2_type[i] > n)
	error->all("Invalid atom type in neighbor exclusion list");
      ex_type[ex1_type[i]][ex2_type[i]] = 1;
      ex_type[ex2_type[i]][ex1_type[i]] = 1;
    }
  }

  if (nex_group) {
    delete [] ex1_bit;
    delete [] ex2_bit;
    ex1_bit = new int[nex_group];
    ex2_bit = new int[nex_group];

    for (i = 0; i < nex_group; i++) {
      ex1_bit[i] = group->bitmask[ex1_group[i]];
      ex2_bit[i] = group->bitmask[ex2_group[i]];
    }
  }

  if (nex_mol) {
    delete [] ex_mol_bit;
    ex_mol_bit = new int[nex_mol];

    for (i = 0; i < nex_mol; i++)
      ex_mol_bit[i] = group->bitmask[ex_mol_group[i]];
  }

  // ------------------------------------------------------------------
  // half and full pairwise neighbor lists

  // determine whether to build half and full lists

  maxlocal = atom->nmax;
  int half_previous = half;
  int full_previous = full;

  half_once = full_once = 0;
  half_every = full_every = 0;
  if (force->pair) half_every = force->pair->neigh_half_every;
  if (force->pair) full_every = force->pair->neigh_full_every;

  for (i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->neigh_half_every) half_every = 1;
    if (modify->fix[i]->neigh_full_every) full_every = 1;
    if (modify->fix[i]->neigh_half_once) half_once = 1;
    if (modify->fix[i]->neigh_full_once) full_once = 1;
  }

  half = full = 0;
  if (half_every || half_once || half_command) half = 1;
  if (full_every || full_once) full = 1;
  half = 1;

  // setup/delete memory for half and full lists

  if (half == 0 && half_previous) {
    memory->sfree(numneigh);
    memory->sfree(firstneigh);
    for (i = 0; i < maxpage; i++) memory->sfree(pages[i]);
    memory->sfree(pages);
    pages = NULL;
    maxpage = 0;
  } else if (half && half_previous == 0) {
    numneigh =
      (int *) memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh");
    firstneigh =
      (int **) memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh");
    add_pages(0);
  } else if (half && half_previous) {
    memory->sfree(numneigh);
    memory->sfree(firstneigh);
    numneigh =
      (int *) memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh");
    firstneigh =
      (int **) memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh");
  }

  if (full == 0 && full_previous) {
    memory->sfree(numneigh_full);
    memory->sfree(firstneigh_full);
    for (i = 0; i < maxpage_full; i++) memory->sfree(pages_full[i]);
    memory->sfree(pages_full);
    pages_full = NULL;
    maxpage_full = 0;
  } else if (full && full_previous == 0) {
    numneigh_full =
      (int *) memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_full");
    firstneigh_full =
      (int **) memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_full");
    add_pages_full(0);
  } else if (full && full_previous) {
    memory->sfree(numneigh_full);
    memory->sfree(firstneigh_full);
    numneigh_full =
      (int *) memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_full");
    firstneigh_full =
      (int **) memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_full");
  }

  // setup/delete memory for shear history neighbor lists
  // history = index of granular shear history fix if it exists
  
  int history_previous = history;
  history = -1;
  if (force->pair_match("gran/history") || force->pair_match("gran/hertzian"))
    for (i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"SHEAR_HISTORY") == 0) history = i;

  if (history == -1 && history_previous >= 0) {
    memory->sfree(firsttouch);
    memory->sfree(firstshear);
    for (i = 0; i < maxpage; i++) memory->sfree(pages_touch[i]);
    for (i = 0; i < maxpage; i++) memory->sfree(pages_shear[i]);
    memory->sfree(pages_touch);
    memory->sfree(pages_shear);
    pages_touch = NULL;
    pages_shear = NULL;
  } else if (history >= 0 && history_previous == -1) {
    firsttouch = (int **) memory->smalloc(maxlocal*sizeof(int *),"firsttouch");
    firstshear = (double **)
      memory->smalloc(maxlocal*sizeof(double *),"firstshear");
    add_pages_history(0);
  } else if (history >= 0 && history_previous >= 0) {
    memory->sfree(firsttouch);
    memory->sfree(firstshear);
    firsttouch = (int **) memory->smalloc(maxlocal*sizeof(int *),"firsttouch");
    firstshear = (double **)
      memory->smalloc(maxlocal*sizeof(double *),"firstshear");
  }

  // setup/delete memory for rRESPA neighbor lists
  // respa = 1 if rRESPA requires extra neighbor lists
  // set neighbor cutoffs for multiple lists

  int respa_previous = respa;
  respa = 0;
  if (update->whichflag == 0 && strcmp(update->integrate_style,"respa") == 0) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  if (respa && half_every == 0)
    error->all("Cannot use rRESPA with full neighbor lists");

  if (respa == 0 && respa_previous) {
    memory->sfree(numneigh_inner);
    memory->sfree(firstneigh_inner);
    for (i = 0; i < maxpage; i++) memory->sfree(pages_inner[i]);
    memory->sfree(pages_inner);
    pages_inner = NULL;
    maxpage_inner = 0;
    if (respa_previous == 2) {
      memory->sfree(numneigh_middle);
      memory->sfree(firstneigh_middle);
      for (i = 0; i < maxpage; i++) memory->sfree(pages_middle[i]);
      memory->sfree(pages_middle);
      pages_middle = NULL;
      maxpage_middle = 0;
    }
  } else if (respa && respa_previous == 0) {
    numneigh_inner = (int *) 
      memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_inner");
    firstneigh_inner = (int **) 
      memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_inner");
    add_pages_inner(0);
    if (respa == 2) {
      numneigh_middle = (int *) 
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_middle");
      firstneigh_middle = (int **) 
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_middle");
      add_pages_middle(0);
    }
  } else if (respa && respa_previous) {
    memory->sfree(numneigh_inner);
    memory->sfree(firstneigh_inner);
    numneigh_inner = (int *) 
      memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_inner");
    firstneigh_inner = (int **) 
      memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_inner");
    if (respa == 2) {
      memory->sfree(numneigh_middle);
      memory->sfree(firstneigh_middle);
      numneigh_middle = (int *) 
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_middle");
      firstneigh_middle = (int **) 
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_middle");
    }
  }

  if (respa) {
    double *cut_respa = ((Respa *) update->integrate)->cutoff;
    cut_inner_sq = (cut_respa[1] + skin) * (cut_respa[1] + skin);
    cut_middle_sq = (cut_respa[3] + skin) * (cut_respa[3] + skin);
    cut_middle_inside_sq = (cut_respa[0] - skin) * (cut_respa[0] - skin);
  }

  // set ptrs to correct half and full build functions
  // cannot combine granular and rRESPA

  if (half) {
    if (strcmp(atom->style,"granular") == 0) {
      if (style == 0) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::granular_nsq_no_newton;
	else half_build = &Neighbor::granular_nsq_newton;
      } else if (style == 1) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::granular_bin_no_newton;
	else half_build = &Neighbor::granular_bin_newton;
      }
    } else if (respa) {
      if (style == 0) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::respa_nsq_no_newton;
	else half_build = &Neighbor::respa_nsq_newton;
      } else if (style == 1) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::respa_bin_no_newton;
	else half_build = &Neighbor::respa_bin_newton;
      }
    } else {
      if (style == 0) {
	if (force->newton_pair == 0) {
	  if (full_every) half_build = &Neighbor::half_full_no_newton;
	  else half_build = &Neighbor::half_nsq_no_newton;
	} else {
	  if (full_every) half_build = &Neighbor::half_full_newton;
	  else half_build = &Neighbor::half_nsq_newton;
	}
      } else if (style == 1) {
	if (force->newton_pair == 0) {
	  if (full_every) half_build = &Neighbor::half_full_no_newton;
	  else half_build = &Neighbor::half_bin_no_newton;
	} else {
	  if (full_every) half_build = &Neighbor::half_full_newton;
	  else half_build = &Neighbor::half_bin_newton;
	}
      }
    }
  } else half_build = NULL;

  if (full) {
    if (style == 0) full_build = &Neighbor::full_nsq;
    else full_build = &Neighbor::full_bin;
  } else full_build = NULL;

  // ------------------------------------------------------------------
  // bond neighbor lists

  // 1st time allocation of bond lists

  if (atom->molecular && atom->nbonds && maxbond == 0) {
    if (nprocs == 1) maxbond = atom->nbonds;
    else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / nprocs);
    bondlist = memory->create_2d_int_array(maxbond,3,"neigh:bondlist");
  }

  if (atom->molecular && atom->nangles && maxangle == 0) {
    if (nprocs == 1) maxangle = atom->nangles;
    else maxangle = static_cast<int> (LB_FACTOR * atom->nangles / nprocs);
    anglelist =  memory->create_2d_int_array(maxangle,4,"neigh:anglelist");
  }

  if (atom->molecular && atom->ndihedrals && maxdihedral == 0) {
    if (nprocs == 1) maxdihedral = atom->ndihedrals;
    else maxdihedral = static_cast<int> 
	   (LB_FACTOR * atom->ndihedrals / nprocs);
    dihedrallist = 
      memory->create_2d_int_array(maxdihedral,5,"neigh:dihedrallist");
  }

  if (atom->molecular && atom->nimpropers && maximproper == 0) {
    if (nprocs == 1) maximproper = atom->nimpropers;
    else maximproper = static_cast<int>
	   (LB_FACTOR * atom->nimpropers / nprocs);
    improperlist = 
      memory->create_2d_int_array(maximproper,5,"neigh:improperlist");
  }

  // set flags that determine which bond neighboring routines to use
  // SHAKE sets bonds and angles negative
  // bond_quartic sets bonds to 0
  // delete_bonds sets all interactions negative


  int bond_off = 0;
  int angle_off = 0;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"shake") == 0)
      bond_off = angle_off = 1;
  if (force->bond && force->bond_match("quartic")) bond_off = 1;

  if (atom->bonds_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (bond_off) break;
      for (m = 0; m < atom->num_bond[i]; m++)
	if (atom->bond_type[i][m] <= 0) bond_off = 1;
    }
  }

  if (atom->angles_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (angle_off) break;
      for (m = 0; m < atom->num_angle[i]; m++)
	if (atom->angle_type[i][m] <= 0) angle_off = 1;
    }
  }

  int dihedral_off = 0;
  if (atom->dihedrals_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (dihedral_off) break;
      for (m = 0; m < atom->num_dihedral[i]; m++)
	if (atom->dihedral_type[i][m] <= 0) dihedral_off = 1;
    }
  }

  int improper_off = 0;
  if (atom->impropers_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (improper_off) break;
      for (m = 0; m < atom->num_improper[i]; m++)
	if (atom->improper_type[i][m] <= 0) improper_off = 1;
    }
  }

  // set ptrs to correct intra-molecular build functions

  if (bond_off) bond_build = &Neighbor::bond_partial;
  else bond_build = &Neighbor::bond_all;

  if (angle_off) angle_build = &Neighbor::angle_partial;
  else angle_build = &Neighbor::angle_all;

  if (dihedral_off) dihedral_build = &Neighbor::dihedral_partial;
  else dihedral_build = &Neighbor::dihedral_all;

  if (improper_off) improper_build = &Neighbor::improper_partial;
  else improper_build = &Neighbor::improper_all;

  // set intra-molecular neighbor list counts to 0
  // in case all are turned off but potential is still defined

  nbondlist = nanglelist = ndihedrallist = nimproperlist = 0;
}

/* ---------------------------------------------------------------------- */

int Neighbor::decide()
{
  if (must_check) {
    int n = update->ntimestep;
    if (restart_check && n == output->next_restart) return 1;
    for (int i = 0; i < fix_check; i++)
      if (n == modify->fix[fixchecklist[i]]->next_reneighbor) return 1;
  }

  ago++;
  if (ago >= delay && ago % every == 0) {
    if (dist_check == 0) return 1;
    else return check_distance();
  } else return 0;
}

/* ---------------------------------------------------------------------- */

int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;

  int nlocal = atom->nlocal;
  double **x = atom->x;
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - xhold[i][0];
    dely = x[i][1] - xhold[i][1];
    delz = x[i][2] - xhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > triggersq) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

/* ----------------------------------------------------------------------
   build all needed neighbor lists every few timesteps
   half, full, bond lists are created as needed
------------------------------------------------------------------------- */

void Neighbor::build()
{
  ago = 0;
  ncalls++;

  // store current nlocal used on this build (used by fix shear/history)

  nlocal_neighbor = atom->nlocal;

  // store current atom positions if needed

  if (dist_check) {
    double **x = atom->x;
    int nlocal = atom->nlocal;
    if (nlocal > maxhold) {
      maxhold = atom->nmax;
      memory->destroy_2d_double_array(xhold);
      xhold = memory->create_2d_double_array(maxhold,3,"neigh:xhold");
    }
    for (int i = 0; i < nlocal; i++) {
      xhold[i][0] = x[i][0];
      xhold[i][1] = x[i][1];
      xhold[i][2] = x[i][2];
    }
  }

  // extend atom arrays if necessary
  // check half/full instead of half_every/full_every so memory will be
  //   allocated correctly whenever build_half() and build_full() are called

  if (atom->nlocal > maxlocal) {
    maxlocal = atom->nmax;

    if (half) {
      memory->sfree(numneigh);
      memory->sfree(firstneigh);
      numneigh = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh");
      firstneigh = (int **)
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh");
    }

    if (full) {
      memory->sfree(numneigh_full);
      memory->sfree(firstneigh_full);
      numneigh_full = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_full");
      firstneigh_full = (int **)
      memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_full");
    }

    if (history >= 0) {
      memory->sfree(firsttouch);
      memory->sfree(firstshear);
      firsttouch = (int **) 
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firsttouch");
      firstshear = (double **)
	memory->smalloc(maxlocal*sizeof(double *),"neigh:firstshear");
    }

    if (respa) {
      memory->sfree(numneigh_inner);
      memory->sfree(firstneigh_inner);
      numneigh_inner = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_inner");
      firstneigh_inner = (int **)
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_inner");
      if (respa == 2) {
	memory->sfree(numneigh_middle);
	memory->sfree(firstneigh_middle);
	numneigh_middle = (int *)
	  memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_middle");
	firstneigh_middle = (int **)
	  memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_middle");
      }
    }
  }

  // extend bin list if necessary

  if (style && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->sfree(bins);
    bins = (int *) memory->smalloc(maxbin*sizeof(int),"bins");
  }

  // list construction for pairs and bonds
  // full comes first in case half is built from full

  if (full_every) (this->*full_build)();
  if (half_every) (this->*half_build)();

  if (atom->molecular) {
    if (atom->nbonds) (this->*bond_build)();
    if (atom->nangles) (this->*angle_build)();
    if (atom->ndihedrals) (this->*dihedral_build)();
    if (atom->nimpropers) (this->*improper_build)();
  }
}

/* ----------------------------------------------------------------------
   one-time call to build a half neighbor list made by other classes
------------------------------------------------------------------------- */

void Neighbor::build_half()
{
  (this->*half_build)();
}

/* ----------------------------------------------------------------------
   one-time call to build a full neighbor list made by other classes
------------------------------------------------------------------------- */

void Neighbor::build_full()
{
  (this->*full_build)();
}

/* ----------------------------------------------------------------------
   setup neighbor binning parameters
   bin numbering is global: 0 = 0.0 to binsize, 1 = binsize to 2*binsize
			    nbin-1 = prd-binsize to binsize
			    nbin = prd to prd+binsize
                            -1 = -binsize to 0.0
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   prd must be filled exactly by integer # of bins
     so procs on both sides of PBC see same bin boundary
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
   stencil() = bin offsets in 1d sense for stencil of surrounding bins
   stencil_full() = bin offsets in 1d sense for stencil for full neighbor list
------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  double cutneighinv = 1.0/cutneigh;

  // test for too many global bins in any dimension due to huge domain

  if (2.0*domain->xprd*cutneighinv > INT_MAX ||
      2.0*domain->yprd*cutneighinv > INT_MAX ||
      2.0*domain->zprd*cutneighinv > INT_MAX)
    error->all("Domain too large for neighbor bins");

  // divide box into bins
  // optimal size is roughly 1/2 the cutoff

  nbinx = static_cast<int> (2.0*domain->xprd*cutneighinv);
  nbiny = static_cast<int> (2.0*domain->yprd*cutneighinv);
  if (force->dimension == 3)
    nbinz = static_cast<int> (2.0*domain->zprd*cutneighinv);
  else nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  binsizex = domain->xprd/nbinx;
  binsizey = domain->yprd/nbiny;
  binsizez = domain->zprd/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  //         add in SMALL for round-off safety

  double coord;
  int mbinxhi,mbinyhi,mbinzhi;

  coord = domain->subxlo - cutneigh - SMALL*domain->xprd;
  mbinxlo = static_cast<int> ((coord-domain->boxxlo)*bininvx);
  if (coord < 0.0) mbinxlo = mbinxlo - 1;
  coord = domain->subxhi + cutneigh + SMALL*domain->xprd;
  mbinxhi = static_cast<int> ((coord-domain->boxxlo)*bininvx);

  coord = domain->subylo - cutneigh - SMALL*domain->yprd;
  mbinylo = static_cast<int> ((coord-domain->boxylo)*bininvy);
  if (coord < 0.0) mbinylo = mbinylo - 1;
  coord = domain->subyhi + cutneigh + SMALL*domain->yprd;
  mbinyhi = static_cast<int> ((coord-domain->boxylo)*bininvy);

  coord = domain->subzlo - cutneigh - SMALL*domain->zprd;
  mbinzlo = static_cast<int> ((coord-domain->boxzlo)*bininvz);
  if (coord < 0.0) mbinzlo = mbinzlo - 1;
  coord = domain->subzhi + cutneigh + SMALL*domain->zprd;
  mbinzhi = static_cast<int> ((coord-domain->boxzlo)*bininvz);

  // extend bins by 1 to insure stencil extent is included

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  mbinz = mbinzhi - mbinzlo + 1;

  // test for too many total local bins due to huge domain

  if (1.0*mbinx*mbiny*mbinz > INT_MAX)
    error->all("Domain too large for neighbor bins");

  // memory for bin ptrs

  mbins = mbinx*mbiny*mbinz;
  if (mbins > maxhead) {
    maxhead = mbins;
    memory->sfree(binhead);
    binhead = (int *) memory->smalloc(maxhead*sizeof(int),"neigh:binhead");
  }

  // create stencil of bins whose closest corner to central bin
  //   is within neighbor cutoff
  // next(xyz) = how far the stencil could possibly extend
  // for partial Newton (newton = 0)
  //   stencil is all surrounding bins including self
  // for full Newton (newton = 1)
  //   stencil is bins to the "upper right" of central bin
  //   stencil does NOT include self
  // for full neighbor list (full = 1)
  //   stencil is all surrounding bins including self, regardless of Newton
  //   stored in stencil_full
  // 3d creates xyz stencil, 2d is only xy

  int nextx = static_cast<int> (cutneigh*bininvx);
  if (nextx*binsizex < cutneigh) nextx++;
  int nexty = static_cast<int> (cutneigh*bininvy);
  if (nexty*binsizey < cutneigh) nexty++;
  int nextz = static_cast<int> (cutneigh*bininvz);
  if (nextz*binsizez < cutneigh) nextz++;

  int nmax = (2*nextz+1) * (2*nexty+1) * (2*nextx+1);
  if (nmax > maxstencil) {
    maxstencil = nmax;
    memory->sfree(stencil);
    stencil = (int *) memory->smalloc(maxstencil*sizeof(int),"neigh:stencil");
  }

  int i,j,k;
  nstencil = 0;
  double cutsq = cutneigh*cutneigh;

  if (force->dimension == 3) {
    if (force->newton_pair == 0) {
      for (k = -nextz; k <= nextz; k++)
	for (j = -nexty; j <= nexty; j++)
	  for (i = -nextx; i <= nextx; i++)
	    if (bin_distance(i,j,k) < cutsq)
	      stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
    } else {
      for (k = 0; k <= nextz; k++)
	for (j = -nexty; j <= nexty; j++)
	  for (i = -nextx; i <= nextx; i++)
	    if (k > 0 || j > 0 || (j == 0 && i > 0))
	      if (bin_distance(i,j,k) < cutsq)
		stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
    }
  } else {
    if (force->newton_pair == 0) {
      for (j = -nexty; j <= nexty; j++)
	for (i = -nextx; i <= nextx; i++)
	  if (bin_distance(i,j,0) < cutsq)
	    stencil[nstencil++] = j*mbinx + i;
    } else {
	for (j = 0; j <= nexty; j++)
	  for (i = -nextx; i <= nextx; i++)
	    if (j > 0 || (j == 0 && i > 0))
	      if (bin_distance(i,j,0) < cutsq)
		stencil[nstencil++] = j*mbinx + i;
    }
  }

  if (full) {
    if (nmax > maxstencil_full) {
      maxstencil_full = nmax;
      memory->sfree(stencil_full);
      stencil_full = (int *) memory->smalloc(maxstencil_full*sizeof(int),
					     "neigh:stencil_full");
    }
    nstencil_full = 0;
    if (force->dimension == 3) {
      for (k = -nextz; k <= nextz; k++)
	for (j = -nexty; j <= nexty; j++)
	  for (i = -nextx; i <= nextx; i++)
	    if (bin_distance(i,j,k) < cutsq)
	      stencil_full[nstencil_full++] = k*mbiny*mbinx + j*mbinx + i;
    } else {
      for (j = -nexty; j <= nexty; j++)
	for (i = -nextx; i <= nextx; i++)
	  if (bin_distance(i,j,0) < cutsq)
	    stencil_full[nstencil_full++] = j*mbinx + i;
    }
  }
}
      
/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double Neighbor::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0)
    delx = (i-1)*binsizex;
  else if (i == 0)
    delx = 0.0;
  else
    delx = (i+1)*binsizex;

  if (j > 0)
    dely = (j-1)*binsizey;
  else if (j == 0)
    dely = 0.0;
  else
    dely = (j+1)*binsizey;

  if (k > 0)
    delz = (k-1)*binsizez;
  else if (k == 0)
    delz = 0.0;
  else
    delz = (k+1)*binsizez;
 
  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   modify parameters of the pair-wise neighbor build
------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      every = atoi(arg[iarg+1]);
      if (every <= 0) error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      delay = atoi(arg[iarg+1]);
      if (delay < 0) error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) dist_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) dist_check = 0;
      else error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      pgsize = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      oneatom = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");

      if (strcmp(arg[iarg+1],"type") == 0) {
	if (iarg+4 > narg) error->all("Illegal neigh_modify command");
	if (nex_type == maxex_type) {
	  maxex_type += EXDELTA;
	  ex1_type = (int *) memory->srealloc(ex1_type,maxex_type*sizeof(int),
					      "neigh:ex1_type");
	  ex2_type = (int *) memory->srealloc(ex2_type,maxex_type*sizeof(int),
					      "neigh:ex2_type");
	}
	ex1_type[nex_type] = atoi(arg[iarg+2]);
	ex2_type[nex_type] = atoi(arg[iarg+3]);
	nex_type++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"group") == 0) {
	if (iarg+4 > narg) error->all("Illegal neigh_modify command");
	if (nex_group == maxex_group) {
	  maxex_group += EXDELTA;
	  ex1_group = 
	    (int *) memory->srealloc(ex1_group,maxex_group*sizeof(int),
				     "neigh:ex1_group");
	  ex2_group = 
	    (int *) memory->srealloc(ex2_group,maxex_group*sizeof(int),
				     "neigh:ex2_group");
	}
	ex1_group[nex_group] = group->find(arg[iarg+2]);
	ex2_group[nex_group] = group->find(arg[iarg+3]);
	if (ex1_group[nex_group] == -1 || ex2_group[nex_group] == -1)
	  error->all("Invalid group ID in neigh_modify command");
	nex_group++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"molecule") == 0) {
	if (iarg+3 > narg) error->all("Illegal neigh_modify command");
	if (atom->molecular == 0) {
	  char *str = "Must use molecular atom style with neigh_modify exclude molecule";
	  error->all(str);
	}
	if (nex_mol == maxex_mol) {
	  maxex_mol += EXDELTA;
	  ex_mol_group = 
	    (int *) memory->srealloc(ex_mol_group,maxex_mol*sizeof(int),
				     "neigh:ex_mol_group");
	}
	ex_mol_group[nex_mol] = group->find(arg[iarg+2]);
	if (ex_mol_group[nex_mol] == -1)
	  error->all("Invalid group ID in neigh_modify command");
	nex_mol++;
	iarg += 3;

      } else if (strcmp(arg[iarg+1],"none") == 0) {
	nex_type = nex_group = nex_mol = 0;
	iarg += 2;
      } else error->all("Illegal neigh_modify command");

    } else error->all("Illegal neigh_modify command");
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

int Neighbor::memory_usage()
{
  int bytes = 0;

  bytes += maxhold*3 * sizeof(double);

  if (style == 0) {
    bytes += maxbin * sizeof(int);
    bytes += maxhead * sizeof(int);
    bytes += maxstencil * sizeof(int);
    bytes += maxstencil_full * sizeof(int);
  }

  if (half) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage*pgsize * sizeof(int);
  }

  if (full) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage_full*pgsize * sizeof(int);
  }

  if (history >= 0) {
    bytes += maxlocal * sizeof(int *);
    bytes += maxlocal * sizeof(double *);
    bytes += maxpage*pgsize * sizeof(int);
    bytes += maxpage*pgsize*3 * sizeof(double);
  }

  if (respa) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage_inner*pgsize * sizeof(int);
    if (respa == 2) {
      bytes += maxlocal * sizeof(int);
      bytes += maxlocal * sizeof(int *);
      bytes += maxpage_middle*pgsize * sizeof(int);
    }
  }

  bytes += maxbond*3 * sizeof(int);
  bytes += maxangle*4 * sizeof(int);
  bytes += maxdihedral*5 * sizeof(int);
  bytes += maximproper*5 * sizeof(int);

  return bytes;
}

/* ----------------------------------------------------------------------
   add pages to half or full neighbor list, starting at npage
------------------------------------------------------------------------- */

void Neighbor::add_pages(int npage)
{
  maxpage += PGDELTA;
  pages = (int **) 
    memory->srealloc(pages,maxpage*sizeof(int *),"neigh:pages");
  for (int i = npage; i < maxpage; i++)
    pages[i] = (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages[i]");
}

void Neighbor::add_pages_full(int npage)
{
  maxpage_full += PGDELTA;
  pages_full = (int **) 
    memory->srealloc(pages_full,maxpage_full*sizeof(int *),"neigh:pages_full");
  for (int i = npage; i < maxpage_full; i++)
    pages_full[i] =
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_full[i]");
}

/* ----------------------------------------------------------------------
   add pages to granular neighbor list, starting at npage
------------------------------------------------------------------------- */

void Neighbor::add_pages_history(int npage)
{
  pages_touch = (int **)
    memory->srealloc(pages_touch,maxpage*sizeof(int *),"neigh:pages_touch");
  pages_shear = (double **)
    memory->srealloc(pages_shear,maxpage*sizeof(double *),
		     "neigh:pages_shear");
  for (int i = npage; i < maxpage; i++) {
    pages_touch[i] = (int *)
      memory->smalloc(pgsize*sizeof(int),"neigh:pages_touch[i]");
    pages_shear[i] = (double *)
      memory->smalloc(3*pgsize*sizeof(double),"neigh:pages_shear[i]");
  }
}

/* ----------------------------------------------------------------------
   add pages to rRESPA inner neighbor list, starting at npage_inner
------------------------------------------------------------------------- */

void Neighbor::add_pages_inner(int npage_inner)
{
  maxpage_inner += PGDELTA;
  pages_inner = (int **) 
    memory->srealloc(pages_inner,maxpage_inner*sizeof(int *),
		     "neigh:pages_inner");
  for (int i = npage_inner; i < maxpage_inner; i++)
    pages_inner[i] = 
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_inner[i]");
}

/* ----------------------------------------------------------------------
   add pages to rRESPA middle neighbor list, starting at npage_middle
------------------------------------------------------------------------- */

void Neighbor::add_pages_middle(int npage_middle)
{
  maxpage_middle += PGDELTA;
  pages_middle = (int **) 
    memory->srealloc(pages_middle,maxpage_middle*sizeof(int *),
		     "neigh:pages_middle");
  for (int i = npage_middle; i < maxpage_middle; i++)
    pages_middle[i] = 
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_middle[i]");
}

/* ----------------------------------------------------------------------
   determine if atom j is in special list of atom i
   if it is not, return 0
   if it is and special flag is 0 (both coeffs are 0.0), return -1
   if it is and special flag is 1 (both coeffs are 1.0), return 0
   if it is and special flag is 2 (otherwise), return 1,2,3
     for which neighbor it is (and which coeff it maps to)
------------------------------------------------------------------------- */

int Neighbor::find_special(int i, int j)
{
  int *list = atom->special[i];
  int n1 = atom->nspecial[i][0];
  int n2 = atom->nspecial[i][1];
  int n3 = atom->nspecial[i][2];
  int tag = atom->tag[j];

  for (int i = 0; i < n3; i++) {
    if (list[i] == tag) {
      if (i < n1) {
	if (special_flag[1] == 0) return -1;
	else if (special_flag[1] == 1) return 0;
	else return 1;
      } else if (i < n2) {
	if (special_flag[2] == 0) return -1;
	else if (special_flag[2] == 1) return 0;
	else return 2;
      } else {
	if (special_flag[3] == 0) return -1;
	else if (special_flag[3] == 1) return 0;
	else return 3;
      }
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void Neighbor::bin_atoms()
{
  int i,ibin,nlocal,nall;
  double **x;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  x = atom->x;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin ghost atoms 1st, so will be at end of linked list
  // then bin owned atoms

  for (i = nlocal; i < nall; i++) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  for (i = 0; i < nlocal; i++) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }
}

/* ----------------------------------------------------------------------
   convert atom coords into local bin #
   only ghost atoms will have coord >= boxhi or coord < boxlo
   take special care to insure ghosts are put in correct bins
   this is necessary so that both procs on either side of PBC
     treat a pair of atoms straddling the PBC in a consistent way
------------------------------------------------------------------------- */

int Neighbor::coord2bin(double *x)
{
  int ix,iy,iz;

  if (x[0] >= domain->boxxhi)
    ix = static_cast<int> ((x[0]-domain->boxxhi)*bininvx) + nbinx - mbinxlo;
  else if (x[0] >= domain->boxxlo)
    ix = static_cast<int> ((x[0]-domain->boxxlo)*bininvx) - mbinxlo;
  else
    ix = static_cast<int> ((x[0]-domain->boxxlo)*bininvx) - mbinxlo - 1;
  
  if (x[1] >= domain->boxyhi)
    iy = static_cast<int> ((x[1]-domain->boxyhi)*bininvy) + nbiny - mbinylo;
  else if (x[1] >= domain->boxylo)
    iy = static_cast<int> ((x[1]-domain->boxylo)*bininvy) - mbinylo;
  else
    iy = static_cast<int> ((x[1]-domain->boxylo)*bininvy) - mbinylo - 1;
  
  if (x[2] >= domain->boxzhi)
    iz = static_cast<int> ((x[2]-domain->boxzhi)*bininvz) + nbinz - mbinzlo;
  else if (x[2] >= domain->boxzlo)
    iz = static_cast<int> ((x[2]-domain->boxzlo)*bininvz) - mbinzlo;
  else
    iz = static_cast<int> ((x[2]-domain->boxzlo)*bininvz) - mbinzlo - 1;

  return (iz*mbiny*mbinx + iy*mbinx + ix + 1);
}

/* ----------------------------------------------------------------------
   test if atom pair i,j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
------------------------------------------------------------------------- */

int Neighbor::exclusion(int i, int j, int *type, int *mask, int *molecule)
{
  int m;

  if (nex_type && ex_type[type[i]][type[j]]) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask[i] & ex1_bit[m] && mask[j] & ex2_bit[m]) return 1;
      if (mask[i] & ex2_bit[m] && mask[j] & ex1_bit[m]) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)
      if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
	  molecule[i] == molecule[j]) return 1;
  }

  return 0;
}
