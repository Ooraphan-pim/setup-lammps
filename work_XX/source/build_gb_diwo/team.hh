

#ifndef DLO_TEAM_HH_LOADED
#define DLO_TEAM_HH_LOADED

#include <mpi.h>


class Everybody;


class Team {
public:

  Team( Everybody &  everybody,  int procs_per_team);
  Team( const Team &); // not defined

//data

  MPI_Comm  team_comm;
  int       team_name;           // color for MPI_Comm_split
  int       team_me;             // rank
  int       team_num_procs;

};





#endif   // do not add anything after this #endif
 
