!! build_gb_sub.f90

!! (1) Rotate crystals
!! (2) Cut crystals
!! (2a)  Move by point in DSC primitive cell
!! (3) Worry about boundaries
!! (4) For alloy, randomize type
!! (5) Write data file

module build_gb_sub_dsc_mod
  implicit none


  real*8 :: zero_vec(3) = (/0, 0, 0/)


!!  interface
!!     integer*4 function system(command_line)
!!       character(*), intent(in) :: command_line
!!     end function system
!!  end interface


contains


!! Returns from the middle of the code when doing single_crystal.
!!
subroutine build_gb_dsc(init,                                  &
                        num_atoms_out, area_out,               &
                        energy_out,                            &
                        dist_out,                              &
                        xi_chi_flag)
  use build_gb_param_mod
  use build_gb_control_mod
  use material_mod
  use build_gb_global_mod
  use build_gb_input_mod
  use build_gb_stuff_dsc_mod
  use build_gb_atom_data_mod
  use build_gb_atom_code_dsc_mod
  use build_gb_boundary_mod
  use call_lammps_dsc_mod
  use int_to_string2_mod
  use delete_file_mod
  use system_mod
  implicit none


!! for two boundaries, area_out is the total 2 * L_y * L_z
  type(build_gb_control_type), intent(in):: init
  integer, intent(out) :: num_atoms_out
  real*8, intent(out)  :: area_out
  real*8, intent(out)  :: energy_out
  real*8, intent(out)  :: dist_out
  logical, intent(in)  :: xi_chi_flag      !! if true, create xi/chi files

  logical :: delete_by_e_or_p

  real*8 :: dsc_vec(3)
  real*8 :: dsc_xx_adj
  real*8 :: xx_boundary_adj

  character(256) :: temp_name

  integer :: number_to_delete, save_num_atoms
  logical :: save_lammps_loop_finished         !! not needed?

!!  integer :: r_step
!!  real*8  :: r_e_1, r_e_2

  real*8  :: best_energy, energy, b_e_2
  integer :: b_num_atoms

!!  real*8  :: r_p_xx, r_p_yy, r_p_zz
!!  real*8  :: b_p_xx, b_p_yy, b_p_zz

  real*8  :: b_dist
  integer :: counter, result

  integer :: num_free
  real*8  :: energy_free

  character(512) :: line

  integer :: left_max_ii, right_max_ii

!! If alloy flag is true, the data file has the actual species.  (OBSOLETE INFO)
!! If alloy flag is false, the species in the data file is the crystal the
!!   atom was built in.  (OBSOLETE INFO)

  logical :: alloy_flag

!!!   logical :: just_one_crystal = .false.

  if( init%delete_by_energy .or. init%delete_by_pressure ) then
     delete_by_e_or_p = .true.
  else
     delete_by_e_or_p = .false.
  end if

  dsc_vec = init%dsc_vec
  dsc_xx_adj = init%dsc_xx_adj
  xx_boundary_adj = init%xx_boundary_adj

  call check_build_gb_control(init)     !! Does not return on failure

  call edit_input(init)

!! Three basic sections here
!!    normal = planar
!!    single crystal
!!    island grain
!! single crystal is marked by "method", but island grain has its own flag

  if( method .eq. single_crystal ) then
!!!      just_one_crystal = .true.
     if( xi_chi_flag ) then
        call xi_chi_files()
     end if
     perlb(1) = 0.0d0
     call compute_max_index( right_len, right_max_ii )
     num_atoms_built = 0
     call build_atoms_dsc( right_max_ii, right_trans, right_len,         &
                           -1.00d0, 1, zero_vec, 0.0d0 )
     actual_num_atoms = num_atoms_built
     num_atoms_deleted = 0
     alloy_flag = .false.
     call write_atoms_dsc(alloy_flag, .false.)
     num_atoms_out = actual_num_atoms
     return

  else   !! actually building a boundary
     if( xi_chi_flag ) then
        call xi_chi_files()
     end if

     if( island_flag ) then
        right_len = left_len
     else
        if( init%edge_size .eq. -1.0 ) then
           write(0,*) "edge_size not initialized"
           call exit(1)
        end if
     end if

     call compute_max_index( left_len, left_max_ii )
     call compute_max_index( right_len, right_max_ii )
     num_atoms_built = 0
     call build_atoms_dsc( left_max_ii, left_trans, left_len,       &
                           init%edge_size, -1,                      &   
                           dsc_vec, dsc_xx_adj + xx_boundary_adj )
     call build_atoms_dsc( right_max_ii, right_trans, right_len,    &
                           init%edge_size, 1,                       &
                           zero_vec, xx_boundary_adj )
     actual_num_atoms = num_atoms_built
     num_atoms_deleted = 0

!! debugging
!!$     if( island_flag ) then
!!$        call print_atoms()
!!$        return
!!$     end if
!! end debugging

     lammps_loop_finished = .false.

!! worry about boundary
     global_current_min = 0.0d0  !! Hereafter controlled by enforce_min_dist
     if( min_dist_flag ) then
        call enforce_min_dist( ( .false. ) )
        if( delete_by_e_or_p .or. init%delete_by_pressure ) then
           call mark_orig_pos()
        end if
     end if

  end if      !! else of if( method .eq. single_crystal )

  if( concentration .ne. 0.0 ) then
     alloy_flag = .true.
  else
     alloy_flag = .false.
  end if

  if( two_boundaries ) then
     write(0,*) "build_gb_sub_dsc  dsc version should not have two boundaries?"
     call exit(1)
     area_out = 2.0d0 * ( perub(2) - perlb(2) ) * ( perub(3) - perlb(3) )
!!       length_out = perub(1) - perlb(1)
  else
     area_out = ( perub(2) - perlb(2) ) * ( perub(3) - perlb(3) )
!!       length_out = -1.0d0
  end if

  if( delete_by_e_or_p .and. ( .not. method .eq. single_crystal ) ) then
     if( .not. min_dist_flag ) then
        write(0,*) "build_gb_sub_dsc  internal error  min_dist_flag(2)"
        call abort()
     end if

     save_lammps_loop_finished = lammps_loop_finished
     save_num_atoms = actual_num_atoms

     trial_delete_loop:  &
     do
        if( lammps_loop_finished ) then
           exit trial_delete_loop
        end if
        call enforce_min_dist(.true.)
     end do trial_delete_loop
     number_to_delete = num_atoms_deleted - orig_atoms_deleted
     lammps_loop_finished = save_lammps_loop_finished
     actual_num_atoms = save_num_atoms
     call write_atoms_dsc(alloy_flag, .true.)
  end if

  counter = 0
  lammps_loop:  &
  do
!! output
!!$     if( counter .gt. 0 .and. lammps_loop_finished ) then
!!$        exit lammps_loop     !! we are done
!!$     end if
     counter = counter + 1

     if( .not. delete_by_e_or_p ) then
        call write_atoms_dsc(alloy_flag, .false.)

!! debugging
!!        write(0,*) "actual_num_atoms  ", actual_num_atoms
!!        stop
!! end debugging

     end if

     if( .not. lammps_flag ) then
        if( output_flag ) then
           call exit(1)     !! We want to leave the files, so cannot write over them.
        else
           temp_name = trim(name)//".data"
!!           call delete_file(trim(temp_name), len_trim(temp_name))  !! not on fenrir
           call delete_file(temp_name, len_trim(temp_name))
        end if

        exit lammps_loop
     end if

     if( delete_by_e_or_p .and. ( counter .gt. 1 ) ) then
        call setup_deletion_control_file(init, actual_num_atoms)
        actual_num_atoms = actual_num_atoms - 1
     end if

     call call_lammps_dsc(name, init%lammps_command, num_free, energy_free)
     energy = energy_free - e_per_a * num_free
     if( ( counter .eq. 1 ) .or. ( energy .lt. best_energy ) ) then
        best_energy = energy
        b_e_2 = energy_free
        b_num_atoms = num_free
        b_dist = global_current_min
        line = "mv -f "                            &
               //trim(name)//".temp.restart "      &
               //trim(name)//".min.restart"

!! debugging
!!        write(0,*) trim(line)
!! end debugging
        result = system(trim(line)//char(0))
        if( result .ne. 0 ) then
           write(0,*) "Warning:  system->mv failed"
!! Depending on the control file, the file many not exist  !! call abort()
        end if

!! debugging
        line = "cp -f "                            &
               //trim(name)//".data "              &
               //trim(name)//".data.keep"
        result = system(trim(line)//char(0))
        if( result .ne. 0 ) then
           write(0,*) "Warning:  system->mv .data (debugging) failed"
           call abort()
        end if
!! end debugging

     end if

!! debugging
        line = "cp -f "                            &
               //trim(name)//".data "              &
               //trim(name)//".data."              &
               //trim(int_to_string2(counter))//".keep"
        result = system(trim(line)//char(0))
        if( result .ne. 0 ) then
           write(0,*) "Warning:  system->mv .data (debugging) failed"
           call abort()
        end if
!! end debugging

     write(0,"('energy lammps  ', 3i4, f7.3, i2, f8.4, i8, f9.3)")     &
                init%dsc_int_vec,                                      &
                init%xx_boundary_adj,                                  &
                init%method_num,                                       &
                global_current_min,                                    &
                num_free,                                              &
                1000.0 * energy / area_out

     temp_name = trim(name)//".data"
     call delete_file(temp_name, len_trim(temp_name))
     if( .not. min_dist_flag ) then
        write(0,*) "build_gb_sub_dsc  internal error  min_dist_flag"
        call abort()
     end if
     if( delete_by_e_or_p ) then
        temp_name = trim(name)//".min.in"
        call delete_file(temp_name, len_trim(temp_name))
        if( counter > number_to_delete ) then
           exit lammps_loop
        end if
     else
        call enforce_min_dist(.true.)
        if( lammps_loop_finished ) then
           exit lammps_loop                !! we are done
        end if
     end if
  end do lammps_loop

!! debugging
  if( lammps_flag ) then
     write(0,*) 'num energy calcs  ', counter
  end if
!! end debugging

  if( lammps_flag ) then
     energy_out = b_e_2
     num_atoms_out = b_num_atoms
     dist_out = b_dist
  else
     energy_out = 0.0d0
     num_atoms_out = 0
     dist_out = 0.0d0
  end if
        
  return
end subroutine build_gb_dsc

end module build_gb_sub_dsc_mod
