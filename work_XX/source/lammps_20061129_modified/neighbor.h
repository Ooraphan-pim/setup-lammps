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

#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "lammps.h"

class Neighbor : public LAMMPS {
 public:
  int style;                       // 0 = nsq, 1 = binned
  int every;                       // build every this many steps
  int delay;                       // delay build for this many steps
  int dist_check;                  // 0 = always build, 1 = only if 1/2 dist
  int ago;                         // how many steps ago neighboring occurred
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom

  double skin;                     // skin distance
  double cutneigh;                 // neighbor cutoff

  int ncalls;                      // # of times build has been called
  int ndanger;                     // # of dangerous builds
  int nlocal_neighbor;             // nlocal at last build

  int half;                        // 0/1 if half pair list ever built
  int full;                        // 0/1 if full pair list ever built
  int half_every;                  // 0/1 if half list built every step
  int full_every;                  // 0/1 if full list built every step
  int half_once;                   // 0/1 if half pair list built occasionally
  int full_once;                   // 0/1 if full pair list built occasionally
  int half_command;                // 0/1 if command requires half list

  int *numneigh;                   // # of half neighbors for each atom
  int **firstneigh;                // ptr to 1st half neighbor of each atom

  int *numneigh_full;              // # of full neighbors for each atom
  int **firstneigh_full;           // ptr to 1st full neighbor of each atom

  int nbondlist;                   // list of bonds to compute
  int **bondlist;

  int nanglelist;                  // list of angles to compute
  int **anglelist;

  int ndihedrallist;               // list of dihedrals to compute
  int **dihedrallist;

  int nimproperlist;               // list of impropers to compute
  int **improperlist;
                                   // granular neighbor list
  int **firsttouch;                // ptr to first touch flag for each atom
  double **firstshear;             // ptr to first shear values for each atom

                                   // rRESPA neighbor lists
  int *numneigh_inner;             // # of inner pair neighbors for each atom
  int **firstneigh_inner;          // ptr to inner 1st neigh of each atom
  int *numneigh_middle;            // # of middle pair neighbors for each atom
  int **firstneigh_middle;         // ptr to middle 1st neigh of each atom

  Neighbor();
  ~Neighbor();
  void init();
  int decide();               // decide whether to build or not
  int check_distance();       // check max distance moved since last build
  void setup_bins();          // setup bins based on box and cutoff
  void build();               // create all neighbor lists (half,full,bond)
  void build_half();          // one-time creation of half neighbor list
  void build_full();          // one-time creation of full neighbor list
  void modify_params(int, char**);  // modify parameters of neighbor build
  int memory_usage();         // tally memory usage
  
 private:
  int me,nprocs;

  int maxlocal;                    // size of numneigh, firstneigh arrays
  int maxbond,maxangle,maxdihedral,maximproper;   // size of bond lists

  int must_check;                  // 1 if must check other classes to reneigh
  int restart_check;               // 1 if restart enabled, 0 if no
  int fix_check;                   // # of fixes that induce reneigh
  int *fixchecklist;               // which fixes to check

  double **cutneighsq;             // neighbor cutoff squared for type pairs
  double triggersq;                // sq distance to trigger build on

  double **xhold;                  // atom coords at last neighbor build
  int maxhold;                     // size of xhold array

  int nbinx,nbiny,nbinz;           // # of global bins
  int *bins;                       // ptr to next atom in each bin
  int maxbin;                      // size of bins array

  int *binhead;                    // ptr to 1st atom in each bin
  int maxhead;                     // size of binhead array

  int nstencil;                    // # of bins in stencil
  int *stencil;                    // stencil list of bin offsets
  int maxstencil;                  // size of stencil

  int mbins;                       // # of local bins and offset
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int **pages;                     // half neighbor list pages
  int maxpage;                     // # of half pages currently allocated

  int **pages_full;                // full neighbor list pages
  int maxpage_full;                // # of full pages currently allocated

  int nstencil_full;               // # of bins in full neighbor stencil
  int *stencil_full;               // full neighbor stencil list of bin offsets
  int maxstencil_full;             // size of full neighbor stencil

                                   // granular neighbor list
  int history;                     // -1 if history not needed
                                   // else stores fix shear history index
  int **pages_touch;               // pages of touch flags
  double **pages_shear;            // pages of shear values

                                   // rRESPA neighbor lists
  int respa;                       // 0 = single neighbor list
                                   // 1 = additional inner list
                                   // 2 = additional inner and middle list
  double inner[2],middle[2];       // rRESPA cutoffs for extra lists
  double cut_inner_sq;		   // outer cutoff for inner neighbor list
  double cut_middle_sq;            // outer cutoff for middle neighbor list
  double cut_middle_inside_sq;     // inner cutoff for middle neighbor list

  int **pages_inner;               // inner neighbor list pages
  int maxpage_inner;               // # of inner pages currently allocated
  int **pages_middle;              // middle neighbor list pages
  int maxpage_middle;              // # of middle pages currently allocated

  int special_flag[4];             // flags for 1-2, 1-3, 1-4 neighbors

  int exclude;                     // 0 if no type/group exclusions, 1 if yes

  int nex_type;                    // # of entries in type exclusion list
  int maxex_type;                  // max # in type list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int maxex_group;                 // max # in group list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  int nex_mol;                     // # of entries in molecule exclusion list
  int maxex_mol;                   // max # in molecule list
  int *ex_mol_group;               // molecule group #'s to exclude
  int *ex_mol_bit;                 // molecule group bits to exclude

  typedef void (Neighbor::*FnPtr)();
  FnPtr half_build;                   // ptr to half pair list functions
  FnPtr full_build;                   // ptr to full pair list functions

  void half_nsq_no_newton();          // via nsq w/ pair newton off
  void half_nsq_newton();             // via nsq w/ pair newton on
  void half_bin_no_newton();          // via binning w/ pair newton off
  void half_bin_newton();             // via binning w/ pair newton on
  void half_full_no_newton();         // via full list w/ pair newton off
  void half_full_newton();            // via full list w/ pair newton on

  void full_nsq();                    // 2 fns for full neighbor lists
  void full_bin();

  void granular_nsq_no_newton();      // 4 fns for granular systems
  void granular_nsq_newton();
  void granular_bin_no_newton();
  void granular_bin_newton();

  void respa_nsq_no_newton();         // 4 fns for multiple respa lists
  void respa_nsq_newton();
  void respa_bin_no_newton();
  void respa_bin_newton();

  FnPtr bond_build;                   // ptr to bond list functions
  void bond_all();                    // bond list with all bonds
  void bond_partial();                // exclude certain bonds

  FnPtr angle_build;                  // ptr to angle list functions
  void angle_all();                   // angle list with all angles
  void angle_partial();               // exclude certain angles

  FnPtr dihedral_build;               // ptr to dihedral list functions
  void dihedral_all();                // dihedral list with all dihedrals
  void dihedral_partial();            // exclude certain dihedrals

  FnPtr improper_build;               // ptr to improper list functions
  void improper_all();                // improper list with all impropers
  void improper_partial();            // exclude certain impropers

  void add_pages(int);                // add pages to half neigh list
  void add_pages_full(int);           // add pages to full neigh list
  void add_pages_history(int);        // add pages to granular neigh list
  void add_pages_inner(int);          // add pages to respa inner list
  void add_pages_middle(int);         // add pages to respa middle list

  void bin_atoms();                     // bin all atoms
  double bin_distance(int, int, int);   // distance between binx
  int coord2bin(double *);              // mapping atom coord to a bin
  int find_special(int, int);           // look for special neighbor
  int exclusion(int, int, int *, int *, int *);  // test for pair exclusion
};

#endif
