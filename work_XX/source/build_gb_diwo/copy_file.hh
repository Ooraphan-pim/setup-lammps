#ifndef DLO_COPY_FILE_HH_LOADED
#define DLO_COPY_FILE_HH_LOADED

#include <cstdlib>
// Changed from two underscores to none for fenrir


extern "C" {


#ifdef BGB_NO_UNDERSCORES
  void copy_file( const char * const  infile_name,
		  const char * const  outfile_name );
#else
  void copy_file__( const char * const  infile_name,
		    const char * const  outfile_name );
#endif

}


#endif  // do not add anything after this #endif
