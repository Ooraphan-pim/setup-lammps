!! build_gb_control.f90

module build_gb_control_mod
  use build_gb_param_mod
  implicit none

  type build_gb_control_type

!! For building island grain
     logical        :: island_flag
     integer        :: island_cut_type      !! values are in build_gb_param_mod
     real*8         :: island_radius        !! for island_cut_type .eq. island_circle
     real*8         :: island_min_in_plane  !! units are alat/2 (radius, in_plane, thickness)
     real*8         :: island_min_thickness

     logical        :: fixed_x_len       !! if true, use these rather than PBV for length 
     real*8         :: left_x_len        !! normal to the boundary
     real*8         :: right_x_len

     real*8         :: alat
     real*8         :: xi_chi_alat  !! different from alat for build_gb_for_mob
     character(3)   :: crystal
!! concentration is optional.  If it is absent or zero, pure material.
!! Otherwise the concentration of species 2.   
     real*8         :: concentration          !! default is zero
     integer        :: left_pbv_ii(3,3)       !! first row is left_xx_pbv, etc.
     integer        :: right_pbv_ii(3,3)
     character(256) :: name
     logical        :: nbrs_by_dist           !! default is .false.
     logical        :: two_boundaries         !! default is .false.
     logical        :: min_dist_flag          !! If min_dist_flag is true, then
     real*8         :: min_dist               !! atoms closer than min_dist are removed
     real*8         :: full_min_dist          !! Distance for nbr list
     real*8         :: low_min_dist           !! If we are doing incremental deletion for
                                              !! enegy calculation where to start.
                                              !! (In that case min_dist is where to end)
     real*8         :: e_per_a                !! energy per bulk atom
   !! If min_dist_flag is false, set min_dist to -1
     integer        :: crystal_to_delete
     integer        :: method                 !! deletion_method
     logical        :: sigma_flag             !! calculate sigma number if possible
     logical        :: lammps_flag
     logical        :: output_flag            !! if false, delete .data, .vectors files
     logical        :: one_at_a_time          !! delete one atom at a time
     real*8         :: nbr_delta              !! How different do nbr distances
                                              !! need to be to be different passes.
                                              !! Default is 0.1 angstrom.
!! stuff for _dsc  ( single boundary version )     
     logical        :: dsc_flag               !! true if we are doing _dsc
     real*8         :: edge_size
     real*8         :: dsc_vec(3)
     real*8         :: dsc_xx_adj
     real*8         :: xx_boundary_adj
     character(256) :: main_name
     integer        :: pbv_num
     integer        :: dsc_num
     integer        :: dsc_int_vec(3)
     integer        :: place_index
     integer        :: method_num
     character(256) :: lammps_control_file
     logical        :: delete_by_energy
     logical        :: show_pressure
     logical        :: delete_by_pressure
     character(256) :: deletion_control_file
     character(256) :: lammps_command         !! e.g.  ./lmp_serial
            !! or mpiexec -n 24 /home/dolmste/lammps/lammps-30Jan06/src/lmp_liberty
     character(256) :: scr_file
     character(256) :: sub_sh_file 

     logical        :: double_edge
     logical        :: jiggle_atoms
     real*8         :: jiggle_amount
     real*8         :: jiggle_region
     integer        :: jiggle_seed

!! stuff for bow_out
     logical :: fix_yy, fix_zz
     real*8  :: fix_yy_width, fix_zz_width

  end type build_gb_control_type


contains

  subroutine write_gb_init( init, out_name )
    use int_to_string2_mod
    implicit none
    type(build_gb_control_type), intent(in) :: init
    character(256), intent(out) :: out_name

   character(256) :: temp_name, out_file
   integer :: dim_ii

   temp_name = trim(init%main_name)
   temp_name = trim(temp_name)//"."//trim(int_to_string2(init%pbv_num))
   temp_name = trim(temp_name)//"."//trim(int_to_string2(init%dsc_num))
   do dim_ii = 1, 3
      temp_name = trim(temp_name)//"."//trim(int_to_string2(init%dsc_int_vec(dim_ii)))
   end do
   temp_name = trim(temp_name)//"."//trim(int_to_string2(init%place_index))
   temp_name = trim(temp_name)//"."//trim(int_to_string2(init%method_num))
   out_name = temp_name
   temp_name = trim(temp_name)//".init"
   out_file = temp_name

   open(99, file=out_file, status="new")

   write(99, "('$material_card')")
   write(99, "(' alat = ', d25.15)") init%alat
   if( init%crystal .eq. "FCC" ) then
      write(99, "(' crystal = ""FCC""')")
   else if( init%crystal .eq. "BCC" ) then
      write(99, "(' crystal = ""BCC""')")
   else
      write(0,*) "write_gb_init  crystal not FCC, BCC?  ", init%crystal
      call abort()
   end if
   write(99, "(' /')")

   write(99, "('$pbv_card')")
   write(99, "('  left_xx_pbv = ', 3i15)") init%left_pbv_ii(1,1:3)
   write(99, "('  left_yy_pbv = ', 3i15)") init%left_pbv_ii(2,1:3)
   write(99, "('  left_zz_pbv = ', 3i15)") init%left_pbv_ii(3,1:3)
   write(99, "('  right_xx_pbv = ', 3i15)") init%right_pbv_ii(1,1:3)
   write(99, "('  right_yy_pbv = ', 3i15)") init%right_pbv_ii(2,1:3)
   write(99, "('  right_zz_pbv = ', 3i15)") init%right_pbv_ii(3,1:3)
   write(99, "(' /')")

   write(99, "('$control_card')")
   write(99,*) " main_name = """, trim(init%main_name), """"
   write(99, "('  pbv_num = ', i15)") init%pbv_num
   write(99, "('  edge_size = ', d25.15)") init%edge_size
   write(99, "('  min_dist = ', d25.15)") init%min_dist
   write(99, "('  low_min_dist = ', d25.15)") init%low_min_dist
   write(99, "('  e_per_a = ', d25.15)") init%e_per_a
   write(99, "('  dsc_num = ', i15)") init%dsc_num
   write(99, "('  dsc_int_vec = ', 3i15)") init%dsc_int_vec
   write(99, "('  place_index = ', i15)") init%place_index
   write(99, "('  method_num = ', i15)") init%method_num
   write(99, "('  dsc_vec = ', 3d25.15)") init%dsc_vec
   write(99, "('  dsc_xx_adj = ', d25.15)") init%dsc_xx_adj
   write(99, "('  xx_boundary_adj = ', d25.15)") init%xx_boundary_adj
   if( init%method .eq. delete_by_averaging ) then
      write(99, "('  deletion_method = ""delete_by_averaging""')")
   else if( init%method .eq. fixed_crystal ) then
      write(99, "('  deletion_method = ""fixed_crystal""')")
   else if( init%method .eq. single_crystal ) then
      write(99, "('  deletion_method = ""single_crystal""')")
   else if( init%method .eq. no_deletion ) then
      write(0,*) "write_gb_init  method is no_deletion"
      write(0,*) "  currently not supported?  how did we get here?"
      call abort()
   else
      write(0,*) "write_gb_init  invalid init%method  ", init%method
      call abort()
   end if
   if( init%method .eq. fixed_crystal ) then
      write(99, "('  crystal_to_delete = ',i15)") init%crystal_to_delete
   end if
   write(99, "('  nbr_delta = ', d25.15)") init%nbr_delta
   if( init%one_at_a_time ) then
      write(99, "('  one_at_a_time = .true.')")
   else
      write(99, "('  one_at_a_time = .false.')")
   end if
   write(99, "('  lammps_flag = .true.')")
   write(99, "('  output_flag = .false.')")
   write(99,*) ' lammps_control_file = "'//trim(init%lammps_control_file)//'"'
   if( init%delete_by_energy ) then
      write(99, "('  delete_by_energy = .true.')")
   else
      write(99, "('  delete_by_energy = .false.')")
   end if
   if( init%delete_by_pressure ) then
      write(99, "('  delete_by_pressure = .true.')")
   else
      write(99, "('  delete_by_pressure = .false.')")
   end if
   if( init%delete_by_energy .or. init%delete_by_pressure ) then
      if( init%show_pressure ) then
         write(99, "('  show_pressure = .true.')")
      end if
      write(99,*) ' deletion_control_file = "'//trim(init%deletion_control_file)//'"'
   end if
   write(99,*) ' lammps_command = "'//trim(init%lammps_command)//'"'
   write(99,*) ' scr_file = "'//trim(init%scr_file)//'"'
   write(99,*) ' sub_sh_file = "'//trim(init%sub_sh_file)//'"'
   if( init%fixed_x_len ) then
      write(99, "('  fixed_x_len = .true.')")
      write(99, "('  left_x_len = ', d25.15)") init%left_x_len
      write(99, "('  right_x_len = ', d25.15)") init%right_x_len
   else
      write(99, "('  fixed_x_len = .false.')")
   end if

   if( init%double_edge ) then
      write(99, "('  double_edge = .true.')")
   else
      write(99, "('  double_edge = .false.')")
   end if
   if( init%jiggle_atoms ) then
      write(99, "('  jiggle_atoms = .true.')")
      write(99, "('  jiggle_amount = ', d25.15)") init%jiggle_amount
      write(99, "('  jiggle_region = ', d25.15)") init%jiggle_region
      write(99, "('  jiggle_seed = ', i15)") init%jiggle_seed
   else
      write(99, "('  jiggle_atoms = .false.')")
   end if

   write(99, "(' /')")

   close(99)

   return
 end subroutine write_gb_init
   





  subroutine initialize_build_gb_control(init)
    implicit none
    type(build_gb_control_type), intent(out) :: init

    init%name = ""
    init%alat = -1.0d0
    init%xi_chi_alat = -1.0d0
    init%crystal = ""
    init%left_pbv_ii = 0
    init%right_pbv_ii = 0
    init%method = -1
    init%crystal_to_delete = -1
    init%nbrs_by_dist = .false.
    init%two_boundaries = .false.
    init%concentration = 0.0d0
    init%sigma_flag = .true.
    init%lammps_flag = .false.
    init%output_flag = .true.
    init%one_at_a_time = .false.
    init%low_min_dist = -1.0d0
    init%full_min_dist = -1.0d0
    init%e_per_a = 1.0d0
    init%edge_size = -1.0d0   !! build_gb_dsc must check this itself!
    init%nbr_delta = 0.1d0

    init%fixed_x_len = .false.
    init%left_x_len = -1.0d0
    init%right_x_len = -1.0d0

!! _dsc stuff
    init%edge_size = -1.0d0   !! build_gb_dsc must check this itself!
!! the other _dsc stuff is "internal", so just trust it
    init%delete_by_energy = .false.  !! default
    init%show_pressure = .false.
    init%delete_by_pressure = .false.
    init%double_edge = .false.
    init%jiggle_atoms = .false.
    init%jiggle_amount = -1.0d0
    init%jiggle_region = -1.0d0
    init%jiggle_seed = -1

    return
  end subroutine initialize_build_gb_control


  subroutine check_build_gb_control(init)
    implicit none
    type(build_gb_control_type), intent(in) :: init

    if( init%alat .lt. 0 ) then
       write(0,*) 'alat missing or negative'
       call exit(1)
    end if
    if( init%name .eq. "" ) then
       if( (.not. init%dsc_flag) .or. init%main_name .eq. " " ) then
          write(0,*) "name missing (or blank)"
          call exit(1)
       end if
    end if
    if( init%method .eq. -1 ) then
       write(0,*) "method missing or -1"
       call exit(1)
    end if
    if(       init%min_dist .ne. -1.0d0             &
        .and. ( .not. init%min_dist_flag ) ) then
       write(0,*) "not min_dist_flag, but min_dist not -1"
       call exit(1)
    end if
    if( init%lammps_flag ) then
       if( init%low_min_dist .eq. -1.0d0 ) then
          write(0,*) "lammps_flag, but low_min_dist is -1"
          call exit(1)
       end if
       if( init%low_min_dist .gt. init%min_dist ) then
          write(0,*) "low_min_dist .gt. min_dist"
          call exit(1)
       end if
       if( init%e_per_a .eq. 1.0d0 ) then
          write(0,*) "lammps_flag, but e_per_a missing (or +1)"
          call exit()
       end if
    end if
    if( init%lammps_flag .and. init%output_flag ) then
       write(0,*) "lammps_flag and output_flag both on"
       call exit(1)
    end if
    if( init%one_at_a_time .and. ( .not. init%lammps_flag ) ) then
       write(0,*) "one_at_a_time, but not lammps_flag"
       call exit(1)
    end if

    return
  end subroutine check_build_gb_control


end module build_gb_control_mod
