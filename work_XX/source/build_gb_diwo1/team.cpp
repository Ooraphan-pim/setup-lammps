

#include "team.hh"
#include "everybody.hh"
#include <mpi.h>
#include <iostream>

using std::cerr;
using std::endl;


Team::Team( Everybody &  everybody,
	    int  procs_per_team      )
{
  int  mpi_initialized;
  int  result = MPI_Initialized( &(mpi_initialized) );
  if( ! mpi_initialized ) {
    cerr << "Team::Team  MPI_Initialized returned false" << endl;
    abort();
  }
  if( procs_per_team == 1 ) {
    if( everybody.all_me == 0 ) {
      team_name = MPI_UNDEFINED;
    } else {
      team_name = everybody.all_me;
    }
    everybody.num_teams = everybody.all_num_procs - 1;
    for( int  ii = 0; ii < everybody.num_teams; ++ii ) {
      Everybody::Team_data  team_data;
      team_data.num_procs = 1;
      team_data.root_proc = 1 + ii;
      everybody.team_table.push_back( team_data );
    }
  } else if( procs_per_team == 2 ) {
    everybody.num_teams = ( everybody.all_num_procs / 2 ) - 1;
    if( 2 * everybody.num_teams + 2 != everybody.all_num_procs ) {
      cerr << "procs_per_team is 2, num_procs is not even" << endl;
      abort();
    }
    team_name = everybody.all_me / 2;
    if( team_name == 0 ) {
      team_name = MPI_UNDEFINED;
    }
    for( int  ii = 0; ii < everybody.num_teams; ++ii ) {
      Everybody::Team_data  team_data;
      team_data.num_procs = 2;
      team_data.root_proc = 2 * ii;
      everybody.team_table.push_back( team_data );
    }
    
  } else {
    cerr << "procs_per_team > 2 not supported (yet)." << endl;
    abort();
  }

  // debugging
  //  cerr << "rank  " << everybody.all_me << " before MPI_Comm_split.  team_name  " << team_name << endl;
  // end debugging

  MPI_Comm_split( everybody.all_comm,
		  team_name,
		  0,
		  &(team_comm)            );

  // debugging
  //  cerr << "rank  " << everybody.all_me << " after MPI_Comm_split.  team_comm  " << team_comm << endl;
  // end debugging

  if( team_name == MPI_UNDEFINED ) {
    team_me = -1;
    team_num_procs = -1;
  } else {
    MPI_Comm_size( team_comm, &team_num_procs );
    MPI_Comm_rank( team_comm, &team_me );


    // debugging
    if( team_num_procs != procs_per_team ) {
      cerr << "Team::Team  internal error team_num_procs" << endl;
      abort();
    }
    // end debugging
  }
  
  // debugging
  //  cerr << "rank  " << everybody.all_me << " team_num_procs  " << team_num_procs << "  team_me  "  << team_me << endl;
  // end debugging

}
