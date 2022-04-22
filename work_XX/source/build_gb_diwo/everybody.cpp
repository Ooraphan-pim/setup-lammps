
#include "everybody.hh"
#include <mpi.h>
#include <iostream>
#include <cstdlib>
using std::abort;

int Everybody::count_to_prevent_dup = 0;


Everybody::Everybody(int *  argc,  char ***  argv)
  : all_comm( MPI_COMM_WORLD ),
    all_me(0),
    all_num_procs(0),
    num_teams(0),
    team_table() {

  if( Everybody::count_to_prevent_dup == 0 ) {
    Everybody::count_to_prevent_dup = 1;
  } else {
    std::cerr << "Attempt to create a second Everybody." << std::endl;
    abort();
  }

  MPI_Init(argc, argv);

  MPI_Comm_size( all_comm, &all_num_procs );
  MPI_Comm_rank( all_comm, &all_me );
}

Everybody::~Everybody()
{
  MPI_Finalize();
}
