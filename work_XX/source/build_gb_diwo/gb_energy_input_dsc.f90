

module gb_energy_input_dsc_mod
  implicit none

  logical :: island_flag
  character(256) :: island_cut_type
  real*8  :: island_radius
  real*8  :: island_min_in_plane
  real*8  :: island_min_thickness
  logical :: build_xi_chi_files     !! just for island.  if true, just xi_chi_files

  logical :: fixed_x_len  !! if true, min_len(and left and right) are actual, repeat dist ignored
  logical :: require_two_repeats
  real*8 :: edge_size
  real*8 :: desired_len         !! includes edge
  real*8 :: min_len             !! includes edge
  real*8 :: min_len_left, min_len_right         !! includes edge
  real*8 :: min_width
  real*8 :: min_dist
  real*8 :: low_min_dist
  real*8 :: e_per_a
  real*8 :: alat
  real*8 :: zero_t_alat       !! Only used by gb_build_for_mob
  real*8 :: max_angle
  real*8 :: nbr_delta         !! What to consider different in nbr distances
  logical :: adj_name_flag
  character(256) :: adj_name
  logical :: one_at_a_time
  integer :: method
  logical :: do_method(3)
  logical :: lammps_flag
  logical :: output_flag
  logical :: mob_shell_flag    !! Just for build_gb_for_mob, not in "init" record
  integer :: dsc_num
  logical :: dsc_half_flag
  logical :: do_now_flag    !! if true, do minimizations.  if false write out
                            !! init information for later
  character(256) :: main_name
  character(256) :: project_task
  character(256) :: lammps_control_file
  logical :: delete_by_energy
  logical :: show_pressure
  logical :: delete_by_pressure
  character(256) :: deletion_control_file
  character(256) :: lammps_command
  character(256) :: scr_file
  character(256) :: sub_sh_file
  character(256) :: head_file
  character(256) :: line_file
  character(256) :: tail_file
  logical        :: double_edge
  logical        :: jiggle_atoms
  real*8         :: jiggle_amount
  real*8         :: jiggle_region
  integer        :: jiggle_seed

  character(3) :: crystal

contains

  subroutine get_gb_energy_input_dsc(single_crystal_flag, control_file_argnum)
    implicit none
    logical, intent(in) :: single_crystal_flag
    integer, intent(in) :: control_file_argnum

    integer :: iargc

    integer :: num_args
    character(256) :: arg_in
    integer, parameter :: in_unit = 99

    namelist /control_card/ desired_len, min_len, min_width,            & 
                            min_len_left, min_len_right,                &
                            min_dist, low_min_dist,                     &
                            e_per_a, max_angle,                         &
                            alat,                                       &
                            zero_t_alat,                                &
                            one_at_a_time, adj_name,                    &
                            lammps_flag, edge_size,                     &
                            output_flag,                                &
                            mob_shell_flag,                             &
                            dsc_num,                                    &
                            dsc_half_flag,                              &
                            nbr_delta,                                  &
                            main_name,                                  &
                            do_now_flag,                                &
                            method,                                     &
                            project_task,                               &
                            lammps_control_file,                        &
                            delete_by_energy,                           &
                            show_pressure,                              &
                            delete_by_pressure,                         &
                            deletion_control_file,                      &
                            lammps_command,                             &
                            scr_file,                                   &
                            sub_sh_file,                                &
                            head_file,                                  &
                            line_file,                                  &
                            tail_file,                                  &
                            crystal,                                    &
                            fixed_x_len,                                &
                            require_two_repeats,                        &
                            double_edge,                                &
                            jiggle_atoms,                               &
                            jiggle_amount,                              &
                            jiggle_region,                              &
                            jiggle_seed,                                &
!! For building island grain
                            island_flag,                                &
                            island_cut_type,                            &
                            island_radius,                              &
                            island_min_in_plane,                        &
                            island_min_thickness,                       &
                            build_xi_chi_files        !! If true, just xi chi files


    build_xi_chi_files = .false.
    island_flag = .false.        !! default
    island_cut_type = ""
    island_radius = -1.0d0               !! units are alat/2
    island_min_in_plane = -1.0d0         !! units are alat/2
    island_min_thickness = -1.0d0        !! units are alat/2

    fixed_x_len = .false.          !! default
    require_two_repeats = .true.   !! defaule
    crystal = "FCC"                !! as default
    main_name = ""
    project_task = ""
    lammps_control_file = "badjuju here"
    delete_by_energy = .false.
    show_pressure = .false.
    delete_by_pressure = .false.
    deletion_control_file = ""
    lammps_command = "./lmp_serial"
    scr_file = "badjuju here"
    sub_sh_file = "badjuju here"
    head_file = "badjuju here"
    line_file = "badjuju here"
    tail_file = "badjuju here"
    double_edge = .false.
    jiggle_atoms = .false.
    jiggle_amount = -1.0d0
    jiggle_region = -1.0d0
    jiggle_seed = -1

    do_now_flag = .true.
    method = -1              !! do all three methods

    edge_size = -1.0d0
    desired_len = -1.0d0
    min_len = 0.0d0        !! default
    min_len_left = -1.0d0
    min_len_right = -1.0d0
    min_width = -1.0d0    
    min_dist = -1.0d0
    low_min_dist = -1.0d0
    e_per_a = 1.0d0
    alat = -1.0d0
    zero_t_alat = -1.0d0
    max_angle = 361.0d0         !! If no max_angle given, do all angles.
    adj_name = ""
    one_at_a_time = .false.
    lammps_flag = .true.
    output_flag = .false.
    mob_shell_flag = .false.
    dsc_num = -1
    dsc_half_flag = .false.
    nbr_delta = 0.1d0

    num_args = iargc()
    if( num_args .ne. control_file_argnum ) then
       write(0,*) "One argument, the name of the control file, except ___start.prl."
       write(0,*) "Then two arguments, the max jobs per qsub, then the control file."
       call exit(1)
    end if

    call getarg(control_file_argnum, arg_in)
    open(in_unit, file=trim(arg_in), status="old", action="read")
    read(in_unit, control_card)

    if( build_xi_chi_files ) then
       if( .not. island_flag ) then
          write(0,*) "build_xi_chi_files, but not island_flag, stopping"
          call exit(1)
       end if
    end if

    if( island_flag ) then
       if( island_cut_type .eq. "" ) then
          write(0,*) "island_cut_type missing"
          call exit(1)
       end if
       if( island_cut_type .eq. "circle" ) then
          if( island_radius .le. 0 ) then
             write(0,*) "island_radius missing (or non-positive)"
             call exit(1)
          end if
       end if
       if( island_min_in_plane .le. 0 ) then
          write(0,*) "island_min_in_plane missing (or non-positive)"
          call exit(1)
       end if
       if( island_min_thickness .le. 0 ) then
          write(0,*) "island_min_thickness missing (or non-positive)"
          call exit(1)
       end if
    end if

!! debugging
!!$    write(0,*) "lammps_control_file"
!!$    write(0,*) trim(lammps_control_file)
!! end debugging

    if( trim(crystal) .ne. "FCC"                &
        .and. trim(crystal) .ne. "BCC" ) then
       write(0,*) "unsupported crystal  ", trim(crystal)
       call exit(1)
    end if
    if( method .eq. -1 ) then
       do_method = .true.
    else
       do_method = .false.
       if( method .lt. 1 .or. method .gt. 3 ) then
          write(0,*) "invalid method: ", method
          call exit(1)
       end if
       do_method(method) = .true.
    end if

    if( .not. do_now_flag ) then
       if( main_name .eq. "" ) then
          write(0,*) "do_now_flag is .false., but main_name missing or ''"
          call exit(1)
       end if
    end if

    if( delete_by_energy .and. delete_by_pressure ) then
       write(0,*) "delete_by_energy and delete_by_pressure both true"
       call exit(1)
    end if
    if( delete_by_energy .or. delete_by_pressure ) then
       if( deletion_control_file .eq. "" ) then
          write(0,*)                 &
           "delete_by_energy or delete_by_pressure is .true., but deletion_control_file is missing or ''"
          call exit(1)
       end if
    end if
    if( delete_by_pressure ) then
       show_pressure = .true.
    end if


    if( adj_name .eq. "" ) then
       adj_name_flag = .false.
    else
       adj_name_flag = .true.
    end if

    if( desired_len .ne. -1.0d0 ) then
       write(0,*) "WARNING: desired_len is no longer supported"
    end if

    if( .not. island_flag ) then
       if( min_len_left .eq. -1.0d0 ) then
          min_len_left = min_len
       end if
       if( min_len_right .eq. -1.0d0 ) then
          min_len_right = min_len
       end if
       if( min_len_left .eq. 0.0d0   .or.  min_len_right .eq. 0.0d0 ) then
          write(0,*) "problem with min_len, or min_len_left, or min_len_right"
          call exit(1)
       end if
       if( min_width .eq. -1.0d0 ) then
          write(0,*) "min_width missing or -1"
          call exit(1)
       end if
    end if

    if( .not. single_crystal_flag ) then
       if( dsc_num .lt. 1 ) then
          write(0,*) "dsc_num missing or .lt. 1"
          call exit(1)
       end if
       if( min_dist .eq. -1.0d0 ) then
          write(0,*) "min_dist missing or -1"
          call exit(1)
       end if
       if( low_min_dist .eq. -1.0d0 ) then
          write(0,*) "low_min_dist missing or -1"
          call exit(1)
       end if
       if( .not. island_flag ) then
          if( edge_size .eq. -1.0d0 ) then
             write(0,*) "edge_size missing or -1"
             call exit(1)
          end if
       end if
       if( e_per_a .eq. 1.0d0 ) then
          write(0,*) "e_per_a missing or -1"
          call exit(1)
       end if
    end if

    if( jiggle_atoms ) then
       if( jiggle_amount .le. 0.0d0 ) then
          write(0,*) "jiggle_atoms, but jiggle_amount missing or not positive"
          call exit(1)
       end if
       if( jiggle_region .le. 0.0d0 ) then
          write(0,*) "jiggle_atoms, but jiggle_region missing or not positive"
          call exit(1)
       end if
       if( jiggle_seed .le. 0 ) then
          write(0,*) "jiggle_atoms, but jiggle_seed missing or not positive"
          call exit(1)
       end if
    else
       if( jiggle_amount .ne. -1.0d0 ) then
          write(0,*) "jiggle_amount given, but not jiggle_atoms"
          call exit(1)
       end if
       if( jiggle_region .ne. -1.0d0 ) then
          write(0,*) "jiggle_region given, but not jiggle_atoms"
          call exit(1)
       end if
       if( jiggle_seed .ne. -1 ) then
          write(0,*) "jiggle_seed given, but not jiggle_atoms"
          call exit(1)
       end if
    end if

    close(in_unit)

    return
  end subroutine get_gb_energy_input_dsc

  subroutine compute_dsc_vec( dsc_mat,           &
                              dsc_int_vec,       &
                              dsc_vec )
    implicit none
    real*8, intent(in)  :: dsc_mat(3,3)
    integer, intent(in) :: dsc_int_vec(3)
    real*8, intent(out) :: dsc_vec(3)

    integer :: ii, jj
    real*8  :: dsc_apply(3)

    dsc_apply = dsc_int_vec / dble( dsc_num )

    dsc_vec = 0.0d0
    do ii = 1, 3
       do jj = 1, 3
          dsc_vec(ii) = dsc_vec(ii) + dsc_apply(jj) * dsc_mat(jj, ii)
       end do
    end do

    return
  end subroutine compute_dsc_vec


end module gb_energy_input_dsc_mod
