#include <fstream>
#include <iostream>
#include "copy_file.hh"
using std::cerr;
using std::endl;

// WARNING, will overwrite output file


// Changed from two to none underscore for fenrir

#ifdef BGB_NO_UNDERSCORES
void copy_file( const char * const  infile_name,
		const char * const  outfile_name )
#else
void copy_file__( const char * const  infile_name,
		  const char * const  outfile_name )
#endif

{
  std::ifstream infile(infile_name);
  if( ! infile ) {
    cerr << "copy_file:  open failed for input file  "  << infile_name << endl;
    exit(1);
  }
  std::ofstream outfile(outfile_name);
  if( ! outfile ) {
    cerr << "copy_file:  open failed for output file  " << outfile_name << endl;
    exit(1);
  }

  char char_in;
  while( infile.get( char_in ) ) {
    outfile.put( char_in  );
  }

  if( ! infile.eof() ) {
    cerr << "copy_file:  read error?" << endl;
    exit(1);
  }
  if( ! outfile ) {
    cerr << "copy_file:  write_error?" << endl;
    exit(1);
  }

}    // copy_file
