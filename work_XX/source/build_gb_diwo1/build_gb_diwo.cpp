#include "build_gb_diwo.hh"
#include "build_gb_diwo_stuff.hh"
#include <mpi.h>
#include <string>
#include <iostream>
#include <fstream>
using std::string;
using std::cin;
using std::cerr;
using std::cout;
using std::endl;

// Changed from two underscores to none for fenrir

extern "C"  int  diwo_start_build(const char *, const int *,
				  char *, int *);

//std::ofstream *  run_out_file;

int  build_gb_diwo( const string &  name_in,
		    MPI_Comm  our_comm        )
{
  int me;
  MPI_Comm_rank(our_comm, &me);

  // Now we call diwo_start_build which will create the first data file, etc.
  bool delete_by_e_or_p;
  int int_delete_by_e_or_p;  // 1 for true, 0 for false
  int name_len;
  char  name_array[256];
  if( me == 0 ) {

    // debugging
    cerr << "build_gb_diwo  name_in:  " << name_in << endl;
    // end debugging

    int  name_in_len = name_in.size();
    if( name_in_len > 254 ) {
      cerr << "ERROR  name_in_len  " << name_in_len << endl;
      abort();
    }
    char  name_in_array[256];
    for( size_t  ii = 0; ii < name_in_len; ++ii ) {
      name_in_array[ii] = name_in[ii];
    }
    for( size_t  ii = name_in_len; ii < 256; ++ii ) {
      name_in_array[ii] = ' ';
    }
    int_delete_by_e_or_p = diwo_start_build(name_in_array, &(name_in_len),
					     name_array, &(name_len)         );
    if( name_len > 256 ) {
      cerr << "ERROR  name_len > 256  name_len  " << name_len << endl;
      abort();
    }
  }

  // debugging
  //  if( me == 0 ) {
  //    cerr << "build_gb_diwo  before:  int_delete_by_e_or_p  " << int_delete_by_e_or_p << endl;
  //  }
  // end debugging

  MPI_Bcast(&int_delete_by_e_or_p, 1, MPI_INT, 0, our_comm);

  // debugging
  //  if( me == 0 ) {
  //    cerr << "build_gb_diwo  before:  name_len  " << name_len << endl;
  //  }
  // end debugging

  MPI_Bcast(&name_len, 1, MPI_INT, 0, our_comm);

  for( size_t  ii = name_len; ii < 256; ++ii ) {
    name_array[ii] = ' ';
  }

  // debugging
  //  if( me == 0 ) {
  //    cerr << "build_gb_diwo  before:  name_array  " << name_array << endl;
  //  }
  // end debugging

  MPI_Bcast(name_array, 256, MPI_CHAR, 0, our_comm);

  // debugging
  //  if( me == 0 ) {
  //    cerr << "build_gb_diwo  after" << endl;
  //  }
  // end debugging

  string  name;
  for( size_t  ii = 0; ii < name_len; ++ii ) {
    name.push_back(name_array[ii]);
  }

  // debugging
  //  if( me == 0 ) {
  //    cerr << "build_gb_diwo  name:  " << name << endl;
  //  }
  // end debugging

  //  if( me == 0 ) {
  //    string run_out_name = name + ".cpp.run.out";
  //    run_out_file = new std::ofstream;
  //    (*run_out_file).open(run_out_name.c_str());
  //  }

  if( int_delete_by_e_or_p == 1 ) {
    delete_by_e_or_p = true;
  } else if( int_delete_by_e_or_p == 0 ) {
    delete_by_e_or_p = false;
  } else {
    cerr << " internal error, int_delete_by_e_or_p, stopping" << endl;
    abort();
  }

  build_gb_diwo_lammps_setup(name, our_comm);

  if( delete_by_e_or_p ) {
    cerr << " delete_by_e_or_p is true.  coming real soon now" << endl;
    abort();
  } else {
    diwo_lammps_loop_not_e_or_p( our_comm );
  }

  // Now in worker.cpp  lammps_close();
  //  (*run_out_file).close();
  //  delete run_out_file;

  return(0);
}
