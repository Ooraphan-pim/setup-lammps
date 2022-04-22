



#ifndef DLO_EVERYBODY_HH_DEFINED
#define DLO_EVERYBODY_HH_DEFINED

#include <vector>
#include <mpi.h>
#include <cstdlib>

// class Team;  // did not work
#include "team.hh"


class Everybody {

  friend  Team::Team( Everybody &, int );

// Nested class defining per-team data in everybody
private:
  struct Team_data {
    friend class Everybody;

    int  num_procs;
    int  root_proc;
  };

public:

  Everybody(int *  argc,  char ***  argv);
  Everybody(const Everybody &);               // Not defined
  ~Everybody();

// data
  MPI_Comm  all_comm;
  int       all_me;
  int       all_num_procs;
  int       num_teams;
  std::vector< Team_data >  team_table;

private:
  static int  count_to_prevent_dup;

};






#endif     // Do not put anything after this #endif
