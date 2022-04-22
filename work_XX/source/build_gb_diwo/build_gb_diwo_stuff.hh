#ifndef DLO_BUILD_GB_DIWO_STUFF_HH_LOADED
#define DLO_BUILD_GB_DIWO_STUFF_HH_LOADED

#include <mpi.h>
//class string;
#include <string>
#include <cstdlib>

void  build_gb_diwo_lammps_setup(const std::string &  name_in, MPI_Comm  our_comm);
void  read_control_line_diwo(int &   str_size,
			     char *  char_str_in,
			     MPI_Comm  our_comm   );
void  diwo_lammps_loop_not_e_or_p(MPI_Comm  our_comm);



#endif    // Do not put anything after this #endif
