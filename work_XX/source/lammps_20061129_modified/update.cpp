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

#include "string.h"
#include "update.h"
#include "neighbor.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#define IntegrateInclude
#define MinimizeInclude
#include "style.h"
#undef IntegrateInclude
#undef MinimizeInclude

/* ---------------------------------------------------------------------- */

Update::Update()
{
  int n;
  char *str;

  whichflag = -1;
  ntimestep = 0;
  first_update = 0;

  maxpair = 0;
  f_pair = NULL;

  unit_style = NULL;
  set_units("lj");

  str = "verlet";
  n = strlen(str) + 1;
  integrate_style = new char[n];
  strcpy(integrate_style,str);
  integrate = new Verlet(0,NULL);

  str = "cg";
  n = strlen(str) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,str);
  minimize = new MinCG();
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  memory->destroy_2d_double_array(f_pair);

  delete [] unit_style;

  delete [] integrate_style;
  delete integrate;

  delete [] minimize_style;
  delete minimize;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // init the appropriate integrate or minimize class
  // if neither (e.g. from write_restart) then just return

  if (whichflag == -1) return;
  else if (whichflag == 0) integrate->init();
  else if (whichflag == 1) minimize->init();

  // only set first_update if a run or minimize is being performed

  first_update = 1;
}

/* ---------------------------------------------------------------------- */

void Update::set_units(char *style)
{
  if (strcmp(style,"lj") == 0) {
    force->boltz = 1.0;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    dt = 0.005;
    neighbor->skin = 0.3;

  } else if (strcmp(style,"real") == 0) {
    force->boltz = 0.001987191;
    force->mvv2e = 48.88821 * 48.88821;
    force->ftm2v = 1.0 / 48.88821 / 48.88821;
    force->nktv2p = 68589.796;
    force->qqr2e = 332.0636;
    force->qe2f = 23.0602;
    dt = 1.0;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"metal") == 0) {
    force->boltz = 8.617e-5;
    force->mvv2e = 1.0365e-4;
    force->ftm2v = 1 / 1.0365e-4;
    force->nktv2p = 1.602e6;
    force->qqr2e = 14.40659;
    force->qe2f = 1.0;
    dt = 0.001;
    neighbor->skin = 2.0;
  } else error->all("Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  if (0) return;      // dummy line to enable else-if macro expansion

#define IntegrateClass
#define IntegrateStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) integrate = new Class(narg-1,&arg[1]);
#include "style.h"
#undef IntegrateClass

  else error->all("Illegal run_style command");

  int n = strlen(arg[0]) + 1;
  integrate_style = new char[n];
  strcpy(integrate_style,arg[0]);
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal min_style command");

  delete [] minimize_style;
  delete minimize;

  if (0) return;      // dummy line to enable else-if macro expansion

#define MinimizeClass
#define MinimizeStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) minimize = new Class();
#include "style.h"
#undef MinimizeClass

  else error->all("Illegal min_style command");

  int n = strlen(arg[0]) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,arg[0]);
}

/* ----------------------------------------------------------------------
   memory usage of update and integrate/minimize
------------------------------------------------------------------------- */

int Update::memory_usage()
{
  int bytes = maxpair*3 * sizeof(double);
  if (whichflag == 0) bytes += integrate->memory_usage();
  else if (whichflag == 1) bytes += minimize->memory_usage();
  return bytes;
}
