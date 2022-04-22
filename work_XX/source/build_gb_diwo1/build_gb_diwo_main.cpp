#include "dispatcher_worker_params.hh"
#include "everybody.hh"
#include "team.hh"
#include "worker.hh"
#include "dispatcher.hh"
#include "drone.hh"
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// extern "C" void _PSC_ftn_init(int argc, char **argv);   // This was for tribe

const int  procs_per_team = 1;


int main(int argc, char **  argv)
{
  //  _PSC_ftn_init(argc, argv);

  Everybody  everybody( &argc, &argv );

  // debugging
  //   if( everybody.all_me == 0 ) {
  //     cerr << "everybody.all_num_procs  " << everybody.all_num_procs << endl;
  //   }
  // end debugging

  Team  team(everybody, procs_per_team);  // second arg is procs_per_team

  if( everybody.num_teams == 0 ) {
    cerr << "not enough processors" << endl;
    abort();
  }

  if( everybody.all_me == 0 ) {
    dispatcher(everybody );
  } else if( procs_per_team == 2 && everybody.all_me == 1 ) {
    drone( everybody );
  } else {
    worker(everybody, team);
  }

  return(0);
}
