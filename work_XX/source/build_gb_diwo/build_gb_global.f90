!! build_gb_global.f90

module build_gb_global_mod
  use material_mod
  implicit none

  logical :: double_edge = .false.

  logical :: jiggle_atoms = .false.
  real*8  :: jiggle_amount = -1.0d0          !! maximum random jiggle
  real*8  :: jiggle_region = -1.0d0          !! distance from nominal boundary to jiggle atoms
  integer :: jiggle_seed = -1   !! gb_energy_cc_diwo increments the original seed
                                !! to give a different seed in each .init file

!! For island
  logical :: island_flag
  integer :: island_cut_type  !! values are in build_gb_param_mod
  real*8  :: island_radius    !! for island_cut_type .eq. island_circle   (alat/2)
  real*8  :: island_min_in_plane          !! (alat/2)
  real*8  :: island_min_thickness         !! (alat/2)

  logical :: fixed_x_len
  real*8  :: left_x_len
  real*8  :: right_x_len

  integer :: left_pbv_ii(3, 3)    !! first row is left_xx_pbv, etc.
  integer :: right_pbv_ii(3, 3)

  real*8  :: left_trans(3, 3)           !! rotation matrix
  real*8  :: right_trans(3, 3)

  real*8  :: left_len(3)          !! box size.  left_len(2) = right_len(2),
  real*8  :: right_len(3)         !!    and the same for z.

  character(256) :: name

  logical :: min_dist_flag        !! If min_dist_flag is true, then
  real*8  :: min_dist             !! atoms closer than min_dist are removed

  real*8  :: full_min_dist        !! distance to use for nbr list

  real*8  :: low_min_dist         !! If we are doing incremental deletion for
                                  !! enegy calculation where to start.
                                  !! (In that case min_dist is where to end)
  real*8  :: e_per_a              !! energy per atom (bulk crystal)

  integer :: crystal_to_delete
  integer :: method               !! deletion_method

  real*8  :: delta  !! called nbr_delta in input files.

  type(material) :: mat

  real*8  :: xi_chi_alat

  real*8  :: perlb(3), perub(3), perlen(3), box_center(3)   !! simulation box

  logical :: vacuum_gap_flag = .false.   !! for build_gb_for_mob 
  real*8  :: half_gap = 0.0d0       !!    accessed by build_gb_atom_code_dsc

  logical :: nbrs_by_dist
  logical :: full_nbr_list
  logical :: two_boundaries
  logical :: sigma_flag

  real*8  :: concentration

  integer :: sigma_number
  real*8  :: global_current_min
  logical :: lammps_flag
  logical :: output_flag

  logical :: lammps_loop_finished
  logical :: one_at_a_time = .false.

  logical :: fix_yy, fix_zz
  real*8  :: fix_yy_width, fix_zz_width

end module build_gb_global_mod
