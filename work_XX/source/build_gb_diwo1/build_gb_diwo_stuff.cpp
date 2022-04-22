

#include "build_gb_diwo_stuff.hh"
// #include "build_gb_diwo.hh"
#include <mpi.h>
#include "library.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
using std::string;
//using std::cin;
using std::cerr;
using std::cout;
using std::endl;
using std::getline;
using std::vector;

const int  MAXLINE = 1024;

// Changed from two underscores to none for fenrir

extern "C"  int  diwo_continue_build_not_e_or_p();

std::ifstream  control_file;

vector<string>  loop_lammps_commands;

// Do everything else for the not_e_or_p_case.
//    (1) read in the read_data and minimize lines from control_file
//    (2) call lammps with them
//    (3) call diwo_continue_build_not_e_or_p
//         repeat (2) and (3) until diwo_continue_build_not_e_or_p says "done"
void  diwo_lammps_loop_not_e_or_p( MPI_Comm  our_comm )
{
  int me;
  MPI_Comm_rank(our_comm, &me);

  loop_lammps_commands.clear();
  while(1) {
    char  char_str_in[MAXLINE];
    int str_size;       // -2 means error; -1 eof, 0 means found "*", so return
    read_control_line_diwo(str_size, char_str_in, our_comm);
    if( str_size < 0 ) { 
      abort();
    } else   if( str_size == 0 ) {
      break;
    }
    loop_lammps_commands.push_back(char_str_in);
  }

//   if( loop_lammps_commands.size() != 2 ) {
//     if( me == 0 ) {
//       cerr << "diwo_lammps_loop_not_e_or_p  loop_lammps_commands.size() not 2:  "  
//         << loop_lammps_commands.size() << endl;
//     }
//     exit(2);
//   }

  control_file.close();

  // Now the lammps loop
  while(1) {
    for( size_t ii = 0; ii < loop_lammps_commands.size();  ++ii ) {
      char temp_c_str[MAXLINE];
      strcpy( temp_c_str, loop_lammps_commands[ii].c_str() );

      //       if( me == 0 ) {
      //         cerr << "temp_c_str  " << temp_c_str << endl;
      //       }

      string str_out = lammps_command( temp_c_str );

      //       if( me == 0 ) {
      //         cerr << "str_out  " << str_out << endl;
      //       }

    }
    int  finished;
    if( me == 0 ) {
      finished = diwo_continue_build_not_e_or_p();
    }
    MPI_Bcast(&finished, 1, MPI_INT, 0, our_comm);
    if( finished == 1 ) {
      break;
    }
  }

  return;
}         // diwo_lammps_loop_not_e_or_p



// First read the original input and do setup
void  build_gb_diwo_lammps_setup(const string &  name_in,
                                 MPI_Comm        our_comm )
{
  int me;
  MPI_Comm_rank(our_comm, &me);


//   string lammps_log = "temp." + name_in + ".min.log";
//   int num_strings = 3;
//   char * strings[3];
//   char * lammps_dummy = "bad_juju_here";
//   strings[0] = lammps_dummy;
//   char * lammps_log_label = "-log";
//   strings[1] = lammps_log_label;
//   char  lammps_log_char[lammps_log.size()+1];
//   strcpy(lammps_log_char, lammps_log.c_str());
//   strings[2] = lammps_log_char;
// //  debugging
// //  if( me == 0 ) {
// //    cerr << "lammps_log_char  " << lammps_log_char << endl;
// //     cerr << "num_strings  " << num_strings << endl;
// //     cerr << "strings[0]  >>" << strings[0] << "<<" << endl;
// //     cerr << "strings[1]  >>" << strings[1] << "<<" << endl;
// //     cerr << "strings[2]  >>"  << strings[2] << "<<" << endl;
// //  } 
// //end debugging
//   lammps_open(num_strings, strings, our_comm);
//   //debugging
//   //  cerr << "after lammps open" << endl;
//   //  exit(3);
//   //end debugging


  if( me == 0 ) {
    string control_name = "temp." + name_in + ".min.in";

    // debugging
    //    cerr << "control_name  " << control_name << endl;
    // end debugging

    control_file.open(control_name.c_str());

    // debugging
    //    cerr << "control_file.good() " <<  control_file.good()  << endl;
    // end debugging
  }

  while(1) {
    char  char_str_in[MAXLINE];
    int str_size;       // -2 means error; -1 eof, 0 means found "*", so return

    // debugging
    //    cerr << "before read_control_line_diwo" << endl;
    // end debugging

    read_control_line_diwo(str_size, char_str_in, our_comm);
    if( str_size < 0 ) {
      abort();
    } else if( str_size == 0 ) {
      break;
    }

// For input file in use with * at begining, this point is never reached

    // debugging
    //    cerr << "after read_control_line_diwo  " << char_str_in << endl;
    // end debugging

    string str_out = lammps_command(char_str_in);

    // debugging
    //    if( me == 0 ) {
    //      cerr << "char_str_in  " << char_str_in << endl;
    //      cerr << "str_out  "  << str_out << endl;
    //    }
    // end debugging

  }       // while(1) which is reading the initialization control file

  return;
}          // build_gb_diwo_lammps_setup


void  read_control_line_diwo( int &   str_size,
                              char *  char_str_in,
                              MPI_Comm  our_comm   )
{
  int me;
  MPI_Comm_rank(our_comm, &me);

  if( me == 0 ) {
    if( ! control_file.good() ) {
      cerr << "not control_file.good() at entrance to read_control_line_diwo" << endl;
      str_size = -2;
    } else {
      while(1) {         // To skip over blank lines, which lammps does not like
        string  str_in;
        getline(control_file, str_in);
        if( ! control_file.good() ) {
          cerr << "read_control_line_diwo: not control_file.good() after read. Unexpected eof?" 
               << endl;
          str_size = -1;           // end of file
          break;    // process as end of file
        } else {
          str_size = 2 + str_in.size();         // +1 for newline, +1 for null
          if( str_size == 2 ) {
            continue;                 // do not process blank lines
          }
        }
        if( str_size <= MAXLINE ) {
          str_in.push_back('\n');
          strcpy(char_str_in, str_in.c_str());
        }
        break;    // process this line
      }
      if( str_size > 0 ) {
        if( char_str_in[0] == '*' ) {     // This says return 
          str_size = 0;
        }
      }
    }
  }

  MPI_Bcast(&str_size, 1, MPI_INT, 0, our_comm);
  if( str_size < 0 ) {
    exit(2);
  } else if( str_size == 0 ) {
    return;
  } else if( str_size > MAXLINE ) {
    if( me == 0 ) {
      cerr << "str_in to long  str_size " << str_size << endl;
      cerr << "char_str_in  " << char_str_in << endl;
    }
    exit(2);
  }
    
  MPI_Bcast(char_str_in, str_size, MPI_CHAR, 0, our_comm);

  return;
}             // read_control_line_diwo
