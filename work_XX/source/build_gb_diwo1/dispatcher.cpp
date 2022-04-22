#include "dispatcher_worker_params.hh"
#include "dispatcher.hh"
#include "everybody.hh"
#include "team.hh"
#include <mpi.h>
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// Changed from two underscores to none for fenrir

extern "C"  void cpp_read_line(char *, int *);
extern "C"  void cpp_write_line(char *);


void dispatcher( const Everybody &  everybody )
{

  bool  done_flag = false;
  int   our_num_teams = everybody.num_teams;
  char  send_buf[DLO_CHAR_BUF_SIZE];

  while(1) {
    // determine the next assignment  (unless we are done)
    int  count;
    if( ! done_flag ) {
      cpp_read_line(send_buf, &(count) );
      if( count > DLO_CHAR_BUF_SIZE ) {
	cerr << "dispatcher  input line too long  count" << count << endl;
	abort();
      }
    }
    if( count == -1 ) {
      done_flag = true;
    } else {
      for( size_t  ii = count; ii < DLO_CHAR_BUF_SIZE; ++ii ) {
	send_buf[ii] = ' ';
      }
    }

    // debugging
    //    cerr << "dispatcher  after read assignment.  done_flag  " << done_flag << endl;
    // end debugging

    // wait for a worker to be available
    char recv_buf[DLO_CHAR_BUF_SIZE];
    MPI_Status  recv_status;

    while(1) {
      MPI_Recv( recv_buf,
		DLO_CHAR_BUF_SIZE,
		MPI_CHAR,
		MPI_ANY_SOURCE,
		MPI_ANY_TAG,
		everybody.all_comm,
		&(recv_status)        );

      // debugging
//             cerr << "dispatcher  after MPI_Recv.  " << endl;
//             cerr << "     dispatcher  recv_status.MPI_TAG  " << recv_status.MPI_TAG << endl;
//             cerr << "     dispatcher  recv_status.MPI_SOURCE  " << recv_status.MPI_SOURCE << endl;
      // end debugging

      if( recv_status.MPI_TAG == DLO_OUTPUT_TAG ) {
	cpp_write_line( recv_buf );
      } else if ( recv_status.MPI_TAG == DLO_READY_TAG ) {
	break;
      } else if ( recv_status.MPI_TAG == DLO_ERROR_ABORT_TAG ) {
	cerr << "received error abort message from worker" << endl;
	abort();
      } else {
	cerr << "dispatcher received invalid tag from worker" << endl;
	abort();
      }
    }       // wait until a worker (or team ) is ready

    int  worker = recv_status.MPI_SOURCE;
    if( done_flag ) {

      // debugging
//             cerr << "dispatcher  before sending DLO_EXIT_TAG" << endl;
//             cerr << "     dispatcher  worker  " << worker << endl;
      // end debugging

      MPI_Send( send_buf,
	        0,
	        MPI_CHAR,
	        worker,
	        DLO_EXIT_TAG,
	        everybody.all_comm );
      --our_num_teams;

      // debugging
//             cerr << "dispatcher  after sending DLO_EXIT_TAG and decrementing our_num_teams" << endl;
//             cerr << "     dispatcher  our_num_teams  " << our_num_teams << endl;
      // end debugging

      if( our_num_teams == 0 ) {
	break;  // we are done
      }
    } else {

      // debugging
//             cerr << "dispatcher  before sending assignment " << endl;
//             cerr << "     dispatcher  send_buf  " << send_buf << endl;
//             cerr << "     dispatcher  worker  " << worker << endl;
      // end debugging

      MPI_Send( send_buf,
	        count,
	        MPI_CHAR,
	        worker,
	        DLO_ASSIGNMENT_TAG,
	        everybody.all_comm );

      // debugging
//             cerr << "dispatcher  after sending assignment" << endl;
      // end debugging

    }
  }    // main loop

  // If two processors per team need to tell the drone to exit.
  if( everybody.team_table[0].num_procs == 2 ) {
    MPI_Send( send_buf,
	      0,
	      MPI_CHAR,
	      1,                 // drone
	      DLO_EXIT_TAG,
	      everybody.all_comm );
  }

  return;
}
