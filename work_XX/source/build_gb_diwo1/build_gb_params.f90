!! build_gb_params.f90

module build_gb_param_mod
  implicit none


!! values for (deletion) method
  integer, parameter :: delete_by_averaging = 11001
  integer, parameter :: fixed_crystal       = 11002
  integer, parameter :: single_crystal      = 11003
  integer, parameter :: no_deletion         = 11004

!! values for island_cut_type
  integer, parameter :: island_circle       = 12001


end module build_gb_param_mod
