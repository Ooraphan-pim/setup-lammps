#include "dispatcher_worker_params.hh"
#include "worker.hh"
#include "everybody.hh"
#include "team.hh"
#include "build_gb_diwo.hh"
#include <mpi.h>
#include "library.h"  // Lammps
#include <iostream>
#include <string>
//#include <unistd.h>    // debugging
using std::cerr;
using std::endl;
using std::string;

void  worker( const Everybody &  everybody,
	     const Team &  team            )
{


  char         recv_buf[DLO_CHAR_BUF_SIZE];   // our assignment
  char         send_buf[DLO_CHAR_BUF_SIZE];
  int          result;
  //  MPI_Request  recv_request;
  MPI_Status   recv_status;



  // Open lammps

//   int team_leader_global_rank = 0;    // Not enough, need PID or somesuch
//   if( team.team_me == 0 ) {
//     team_leader_global_rank = everybody.all_me;
//   }
//   MPI_Bcast( &(team_leader_global_rank),
// 	     1,
// 	     MPI_INT,
// 	     0,
// 	     team.team_comm     );
//   ostringstream  ost;
//   ost << team_leader_global_rank;
//   string lammps_log = 


  string lammps_log = "none";          // No easy way to create unique name yet.
                                       // The actual log file is created later.
  int num_strings = 3;
  char * strings[3];

  char * lammps_dummy = "bad_juju_here";
  strings[0] = lammps_dummy;

  char * lammps_log_label = "-log";
  strings[1] = lammps_log_label;

  char  lammps_log_char[lammps_log.size()+1];
  strcpy(lammps_log_char, lammps_log.c_str());
  strings[2] = lammps_log_char;

//  debugging
//  if( me == 0 ) {    team.team_me is needed 
//    (*run_out_file) << "lammps_log_char  " << lammps_log_char << endl;
//     (*run_out_file) << "num_strings  " << num_strings << endl;
//     (*run_out_file) << "strings[0]  >>" << strings[0] << "<<" << endl;
//     (*run_out_file) << "strings[1]  >>" << strings[1] << "<<" << endl;
//     (*run_out_file) << "strings[2]  >>"  << strings[2] << "<<" << endl;
//  } 
//end debugging
 
  lammps_open(num_strings, strings, team.team_comm);

  //debugging
  //  (*run_out_file) << "after lammps open" << endl;
  //  exit(3);
  //end debugging

  int  task_tag;
  int  num_char;

  for( size_t  ii = 0; ii < DLO_CHAR_BUF_SIZE; ++ii ) {
    send_buf[ii] = ' ';
  }

  while(1) {
    if( team.team_me == 0 ) {

      // debugging
      //      cerr << "worker  " << everybody.all_me << " before MPI_Send  everybody.all_comm  " << everybody.all_comm << endl;
      // end debugging

      result = MPI_Send( send_buf,
			 DLO_CHAR_BUF_SIZE,                // really, nothing
			 MPI_CHAR,
			 0,                // dispatcher
			 DLO_READY_TAG,
			 everybody.all_comm );

      // debugging
      //      cerr << "worker  " << everybody.all_me << "  after MPI_Send  result  " << result << endl;
      // end debugging

      if( result != 0 ) {
	abort();
      }

      // debugging
      //      cerr << "worker  " << everybody.all_me << "  before MPI_Recv  everybody.all_comm  " << everybody.all_comm << endl;
      // end debugging

      result = MPI_Recv(recv_buf,
			DLO_CHAR_BUF_SIZE, 
			MPI_CHAR,
			0,                       // dispatcher
			MPI_ANY_TAG, 
			everybody.all_comm,
			&(recv_status)          );

      // debugging
      //      cerr << "worker  " << everybody.all_me << "  after MPI_Recv  result  " << result << "  recv_status.MPI_TAG  " << recv_status.MPI_TAG << "  recv_buf  " << recv_buf << endl;
      // end debugging

      if( result != 0 ) {
	abort();
      }
      result = MPI_Get_count( &(recv_status), MPI_CHAR, &(num_char) );
      if( result != 0 ) {
	abort();
      }
      if( num_char > DLO_CHAR_BUF_SIZE ) {
	cerr << "ERROR  worker   num_char > DLO_CHAR_BUF_SIZE" << endl;
	abort();
      }


      task_tag = recv_status.MPI_TAG;
    }        // if( team leader )

    // debugging
    //    cerr << "worker  " << everybody.all_me << "  before MPI_Bcast task_tag  team.team_comm  " << team.team_comm << endl;
    // end debugging

    MPI_Bcast( &(task_tag),
	       1,
	       MPI_INT,
	       0,
	       team.team_comm     );

    if( task_tag == DLO_EXIT_TAG ) {
      break;
    } else if( task_tag == DLO_ASSIGNMENT_TAG ) {
      // what we expect
    } else {
      cerr << "worker  bad tag received  " << task_tag << endl;
      abort();
    }

// We have our assignment


    MPI_Bcast( &(num_char),
	       1,
	       MPI_INT,
	       0,
	       team.team_comm     );



    // debugging
//     cerr << "worker  " << everybody.all_me << "  before MPI_Bcast recv_buf team.team_comm  " << team.team_comm << endl;
    // end debugging


    MPI_Bcast( recv_buf,
	       DLO_CHAR_BUF_SIZE,
	       MPI_CHAR,
	       0,
	       team.team_comm     );

    // debugging
//     cerr << "worker  " << everybody.all_me << "  after MPI_Bcast" << endl;
    // end debugging


    // Do our work

    string  name_in;
    for( size_t  ii = 0;  ii < num_char; ++ii ) {
      name_in.push_back(recv_buf[ii]);
    }
    build_gb_diwo( name_in, team.team_comm );

    // debugging
//     cerr << "worker  " << everybody.all_me << "  before MPI_Barrier  team.team_comm  " << team.team_comm << endl;
    // end debugging

    MPI_Barrier( team.team_comm );

    // debugging
//     cerr << "worker  " << everybody.all_me << "  after MPI_Barrier" << endl;
    // end debugging


    // Report back
//     if( team.team_me == 0 ) {
//       int result;
//       result = MPI_Send( recv_buf,
// 			 DLO_CHAR_BUF_SIZE, 
// 			 MPI_CHAR,
// 			 DLO_MANAGER_RANK,
// 			 DLO_OUTPUT_TAG,
// 			 everybody.all_comm  );
//     }


  }

  // debugging
  //  cerr << "worker  " << everybody.all_me << "  before lammps_close()" << endl;
  // end debugging

  lammps_close();

  // debugging
  //  cerr << "worker  " << everybody.all_me << "  after lammps_close()" << endl;
  // end debugging

  return;
}
