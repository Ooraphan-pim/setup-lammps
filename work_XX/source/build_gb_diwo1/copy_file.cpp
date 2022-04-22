#include <fstream>
#include <iostream>
#include "copy_file.hh"
using std::cerr;
using std::endl;

// WARNING, will overwrite output file

// void copy_file_error(const char * const  part_one,
// 		     const char * const  part_two = " ")
// {
//   cerr << "copy_file: " << part_one << " " << part_two << endl;
//   exit(1);
// }

// Changed from two to none underscore for fenrir

void copy_file( const char * const  infile_name,
		const char * const  outfile_name )
{
  std::ifstream infile(infile_name);
  if( ! infile ) {
    //    copy_file_error("open failed for input file", infile_name);
    cerr << "copy_file:  open failed for input file  "  << infile_name << endl;
    exit(1);
  }
  std::ofstream outfile(outfile_name);
  if( ! outfile ) {
    //    copy_file_error("open failed for output file", outfile_name);
    cerr << "copy_file:  open failed for output file  " << outfile_name << endl;
    exit(1);
  }

  char char_in;
  while( infile.get( char_in ) ) {
    outfile.put( char_in  );
  }

  if( ! infile.eof() ) {
    //    copy_file_error("read error?");
    cerr << "copy_file:  read error?" << endl;
    exit(1);
  }
  if( ! outfile ) {
    //    copy_file_error("write error?");
    cerr << "copy_file:  write_error?" << endl;
    exit(1);
  }

}    // copy_file


