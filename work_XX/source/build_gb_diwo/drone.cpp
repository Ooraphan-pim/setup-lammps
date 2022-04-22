#include "dispatcher_worker_params.hh"
#include "drone.hh"
#include "everybody.hh"
#include <mpi.h>
#include <iostream>
#include <cstdlib>
using std::cerr;
using std::endl;
using std::string;
using std::abort;

void  drone( const Everybody &  everybody )
{
  char  recv_buf[DLO_CHAR_BUF_SIZE];
  MPI_Status  recv_status;
  int  result;

  result = MPI_Recv(recv_buf,
		    DLO_CHAR_BUF_SIZE, 
		    MPI_CHAR,
		    0,                       // dispatcher
		    MPI_ANY_TAG, 
		    everybody.all_comm,
		    &(recv_status)          );

  if( recv_status.MPI_TAG != DLO_EXIT_TAG ) {
    cerr << "drone  recv_status.MPI_TAG not DLO_EXIT_TAG  " << recv_status.MPI_TAG << endl;
    abort();
  }

  return;
}
