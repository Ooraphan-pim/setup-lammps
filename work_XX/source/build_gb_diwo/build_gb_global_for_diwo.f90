!! Probably the routines to be called by C++ had better be outside the module






module build_gb_global_for_diwo_mod
!!  use system_mod
  use build_gb_control_mod
  use build_gb_run_out_mod
  implicit none


!!$  interface
!!$     integer*4 function system(command_line)
!!$       character(*), intent(in) :: command_line
!!$     end function system
!!$  end interface


!!$ Apparently iargc does not work if fortran subroutine called from C++
!!$
!!$  interface
!!$     integer*4 function iargc()
!!$     end function iargc
!!$  end interface
!!$
!!$  interface
!!$     subroutine getarg(arg_no, arg)
!!$       integer*4, intent(in) :: arg_no
!!$       character(*), intent(out) :: arg
!!$     end subroutine getarg
!!$  end interface




  type(build_gb_control_type) :: global_init
  character(256)              :: global_temp_name




  logical :: global_dsc_flag
  logical :: global_delete_by_e_or_p
  real*8  :: global_dsc_vec(3)
  real*8  :: global_dsc_xx_adj
  real*8  :: global_xx_boundary_adj
  integer :: global_left_max_ii
  integer :: global_right_max_ii
  logical :: global_alloy_flag
  real*8  :: global_area_out
  integer :: global_number_to_delete
  real*8  :: global_best_energy
  real*8  :: global_b_e_2
  real*8  :: global_b_dist
  integer :: global_b_num_atoms
  real*8  :: global_b_energy_with
  real*8  :: global_b_e_with_2
  integer :: global_b_num_with
  logical :: global_first_time
  integer :: global_counter

contains
  subroutine diwo_read_build_gb_input_dsc(name_in)
    use build_gb_input_mod
    implicit none
    character(256), intent(in) :: name_in

!!$    integer,parameter :: in_unit = 5

    integer,parameter :: in_unit = 23 

!!$    integer :: num_args
!!$    character(256) :: arg_in

!! Apparently iargc() does not work when fortan is called from c++.
!!$    num_args = iargc()
!!$    if( num_args .ne. 1 ) then
!!$       write(0,*) "one arg, the .init file   iargc   ", num_args
!!$       call exit(1)
!!$    end if
!!$    call getarg(1, arg_in)
!!$    open(in_unit, file=trim(arg_in), status="old", action="read")

    open(in_unit, file=trim(name_in), status="old", action="read")

    global_dsc_flag = .true.
    call read_build_gb_input(in_unit, global_init, global_dsc_flag)

!!$    close(in_unit)

    return
  end subroutine diwo_read_build_gb_input_dsc


  subroutine diwo_temp_name()
    use int_to_string2_mod
    implicit none

    integer :: dim_ii

    global_temp_name = trim(global_init%main_name)
    global_temp_name = trim(global_temp_name)//"."          &
         //trim(int_to_string2(global_init%pbv_num))
    global_temp_name = trim(global_temp_name)//"."          &
         //trim(int_to_string2(global_init%dsc_num))
    do dim_ii = 1, 3
       global_temp_name = trim(global_temp_name)//"."       &
            //trim(int_to_string2(global_init%dsc_int_vec(dim_ii)))
    end do
    global_temp_name = trim(global_temp_name)//"."          &
         //trim(int_to_string2(global_init%place_index))
    global_temp_name = trim(global_temp_name)//"."          &
         //trim(int_to_string2(global_init%method_num))

    global_init%name = "temp."//trim(global_temp_name)

    return
  end subroutine diwo_temp_name

!!$  subroutine diwo_sed_control_file
!!$    implicit none
!!$
!!$    integer :: result
!!$    character(256) :: line
!!$
!!$    line = "cat "//trim(global_init%lammps_control_file)         &
!!$         //" | sed ""s/Q/"//"temp."//trim(global_temp_name)    &
!!$         //"/"" > temp."//trim(global_temp_name)//".min.in"
!!$    result = system(line)
!!$    if( result .ne. 0 ) then
!!$       write(run_out_unit,*) 'diwo_sed_control_file: system->sed failed'
!!$       call abort()
!!$    end if
!!$
!!$    return
!!$  end subroutine diwo_sed_control_file


  subroutine diwo_start_build_sub(delete_by_e_or_p)
    use build_gb_global_mod
    use build_gb_input_mod
    use build_gb_atom_data_mod
    use build_gb_atom_code_dsc_mod
    use build_gb_boundary_mod
    use build_gb_stuff_dsc_mod
    use random_numbers_mod
!!    use call_lammps_dsc_mod
!!    use delete_file_mod
    implicit none

    logical, intent(out) :: delete_by_e_or_p

    real*8, parameter :: zero_vec(3) = (/0, 0, 0/)    !! constant

!! temporaries
    logical :: save_lammps_loop_finished
    integer :: save_num_atoms                 

!!     real*8  :: energy
!!     real*8  :: energy_free
!!     integer :: num_free
!!     character(256) :: local_temp_name

    global_first_time = .true.
    global_counter = 0

    if( global_init%delete_by_energy .or. global_init%delete_by_pressure ) then
       global_delete_by_e_or_p = .true.
    else
       global_delete_by_e_or_p = .false.
    end if
    delete_by_e_or_p = global_delete_by_e_or_p

    global_dsc_vec = global_init%dsc_vec
    global_dsc_xx_adj = global_init%dsc_xx_adj
    global_xx_boundary_adj = global_init%xx_boundary_adj

    call check_build_gb_control(global_init)    !! Does not return on failure
    call edit_input(global_init)

    if( jiggle_atoms ) then
       call init_random_from_seed( jiggle_seed )
    end if

    if( method .eq. single_crystal ) then
       write(run_out_unit,*) "diwo_start_build_sub  did not expect single_crystal"
       call abort()
    end if

    if( global_init%edge_size .eq. -1.0 ) then
       write(run_out_unit,*) "diwo_start_build_sub  edge_size not initialized"
       call exit(1)
    end if

    call compute_max_index( left_len, global_left_max_ii )
    call compute_max_index( right_len, global_right_max_ii )

    num_atoms_built = 0
    call build_atoms_dsc( global_left_max_ii, left_trans, left_len,       &
                          global_init%edge_size, -1,                      &   
                          global_dsc_vec,                                 &
                          global_dsc_xx_adj + global_xx_boundary_adj )
    call build_atoms_dsc( global_right_max_ii, right_trans, right_len,    &
                          global_init%edge_size, 1,                       &
                          zero_vec,                                       &
                          global_xx_boundary_adj )
    actual_num_atoms = num_atoms_built
    num_atoms_deleted = 0

    lammps_loop_finished = .false.        !! in build_gb_global_mod

!! worry about boundary the first time
    global_current_min = 0.0d0   !! After this it is handled by enforce_min_dist
    if( min_dist_flag ) then
       call enforce_min_dist( .false. )
       if( global_delete_by_e_or_p ) then
          call mark_orig_pos()
       end if
    end if

    if( concentration .ne. 0.0 ) then
       global_alloy_flag = .true.
    else
       global_alloy_flag = .false.
    end if

    if( two_boundaries ) then
       write(run_out_unit,*) "diwo_start_build_sub  not expecting two boundaries"
       call abort()
    end if

    global_area_out = ( perub(2) - perlb(2) ) * ( perub(3) - perlb(3) )

    if( global_delete_by_e_or_p ) then
       if( .not. min_dist_flag ) then
          write(run_out_unit,*) "diwo_start_build_sub  internal error  min_dist_flag(2)"
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
       global_number_to_delete = num_atoms_deleted - orig_atoms_deleted
       lammps_loop_finished = save_lammps_loop_finished
       actual_num_atoms = save_num_atoms
       call write_atoms_dsc(global_alloy_flag, .true.)
    end if

    if( .not. global_delete_by_e_or_p ) then
       call write_atoms_dsc(global_alloy_flag, .false.)
    end if
     
    if( .not. lammps_flag ) then
       write(run_out_unit,*) "diwo_start_build_sub  lammps_flag is false ?"
       call abort()
    end if


    return
  end subroutine diwo_start_build_sub


!! Because of problems compiling with Intel compilers, the change is made:
!!       call delete_file(trim(_name_), len_trim(_name_))
!!  -->  call delete_file(_name_, len_trim(_name_))
 
  subroutine diwo_continue_build_not_e_or_p_sub(result)
    use build_gb_global_mod
    use build_gb_boundary_mod
    use build_gb_atom_code_dsc_mod
    use get_lammps_result_mod
    use delete_file_mod
    implicit none

    integer, intent(out) :: result   !! 1 means finished, otherwise zero

    integer :: num_free
    real*8  :: energy_free, energy
    integer :: num_with                  !! values including "inner" edge atoms
    real*8  :: total_energy_with, energy_with
    character(256) :: local_temp_name
    character(256) :: copy_in_name, copy_out_name
    integer        :: open_status

    call get_lammps_result( global_init%name, .false., double_edge,                 &
                            num_free, energy_free,                                  &
                            num_with, total_energy_with )

    energy = energy_free - e_per_a * num_free
    energy_with = total_energy_with - e_per_a * num_with

    if( double_edge ) then
       write(run_out_unit,"('energy lammps  ', 3i4, f7.3, i2, f8.4, i8, f9.3, i8, f9.3)")     &
               global_init%dsc_int_vec,                               &
               global_init%xx_boundary_adj,                           &
               global_init%method_num,                                &
               global_current_min,                                    &
               num_free,                                              &
               1000.0d0 * energy / global_area_out,                   &
               num_with,                                              &
               1000.0d0 * energy_with / global_area_out
    else
        write(run_out_unit,"('energy lammps  ', 3i4, f7.3, i2, f8.4, i8, f9.3)")     &
               global_init%dsc_int_vec,                               &
               global_init%xx_boundary_adj,                           &
               global_init%method_num,                                &
               global_current_min,                                    &
               num_free,                                              &
               1000.0d0 * energy / global_area_out
     end if

    global_counter = global_counter + 1

    if( global_first_time .or. (energy .lt. global_best_energy) ) then
       global_first_time = .false.
       global_best_energy = energy
       global_b_e_2 = energy_free
       global_b_num_atoms = num_free
       global_b_dist = global_current_min
       global_b_energy_with = energy_with             !! This is *not* what is being minimized
       global_b_e_with_2 = total_energy_with
       global_b_num_with = num_with

!! If we wrote a restart file, it is because we want to keep it.
!! If not, it will not exist.
!!$       line = "cp -f "                            &
!!$            //trim(global_init%name)//".temp.restart "              &
!!$            //trim(global_init%name)//".min.restart"
!!$       result = system(trim(line))

       copy_in_name = trim(global_init%name)//".temp.restart"
       open(99, file=copy_in_name, iostat = open_status, status="old", action="read")
       close(99)
       if( open_status .eq. 0 ) then
          copy_out_name = trim(global_init%name)//".min.restart"
          call copy_file( trim(copy_in_name)//achar(0),                     &
                          trim(copy_out_name)//achar(0) )
       end if


!!$       line = "cp -f "                            &
!!$            //trim(global_init%name)//".data "              &
!!$            //trim(global_init%name)//".data.keep"
!!$       result = system(trim(line))
!!$       if( result .ne. 0 ) then
!!$          write(run_out_unit,*) "system->cp .data (debugging) failed"
!!$          call abort()
!!$       end if

!! for figuring out excess volume
       copy_in_name = trim(global_init%name)//".data"
       copy_out_name = trim(global_init%name)//".data.keep"
       call copy_file( trim(copy_in_name)//achar(0),                     &
                       trim(copy_out_name)//achar(0) )
    end if

    local_temp_name = trim(global_init%name)//".data"
    call delete_file(local_temp_name, len_trim(local_temp_name))
    local_temp_name = trim(global_init%name)//".temp.restart"
    call delete_file(local_temp_name, len_trim(local_temp_name))

!!     local_temp_name = trim(global_temp_name)//".half_xi.vectors"
!!     call delete_file(trim(local_temp_name), len_trim(local_temp_name))
!!     local_temp_name = trim(global_temp_name)//".half_chi.vectors"
!!     call delete_file(trim(local_temp_name), len_trim(local_temp_name))

    if( lammps_loop_finished ) then
       result = 1           !! we are done

       local_temp_name = trim(global_temp_name)//".result"
       open(99, file=local_temp_name, status="new")
       if( double_edge ) then
          write(99, "(f18.14, i10, f9.3, i10, f9.3)")                   &
                     global_b_dist,                                     &
                     global_b_num_atoms,                                &
                     1000.0d0 * global_best_energy / global_area_out,   &
                     global_b_num_with,                                 &
                     1000.0d0 * global_b_energy_with / global_area_out
       else
          write(99, "(f18.14, i10, f9.3)")                              &
                     global_b_dist,                                     &
                     global_b_num_atoms,                                &
                     1000.0d0 * global_best_energy / global_area_out
       end if

       close(99)

       write(run_out_unit,*) "num energy calcs  ", global_counter

       local_temp_name = "temp."//trim(global_temp_name)//".min.in"
       call delete_file(local_temp_name, len_trim(local_temp_name))
       local_temp_name = "temp."//trim(global_temp_name)//".min.log"
       call delete_file(local_temp_name, len_trim(local_temp_name))
    else
       call enforce_min_dist(.true.)
       call write_atoms_dsc(global_alloy_flag, .false.) 
       result = 0
    end if


    return;
  end subroutine diwo_continue_build_not_e_or_p_sub


end module build_gb_global_for_diwo_mod
