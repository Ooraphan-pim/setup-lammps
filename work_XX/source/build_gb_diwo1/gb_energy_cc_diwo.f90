!! for the intel compiler on fenrir, trying for_flush instead of flush (no good)
 
program gb_energy_cc_diwo
  use gb_energy_input_dsc_mod
  use build_gb_run_out_mod
  use build_gb_param_mod
  use build_gb_control_mod
  use build_gb_sub_dsc_mod
  use call_lammps_dsc_mod
  use delete_file_mod
  use system_mod
  implicit none

  type(build_gb_control_type) init

  integer, parameter :: in_unit = 21
  integer, parameter :: task_unit = 41

  real*8  :: island_gb_energy_result(3)
  real*8  :: island_gb_dist(3)

  real*8  :: dsc_mat(3,3)

  integer :: dsc_int_vec(3)
  real*8  :: dsc_vec(3)
  integer :: dsc_ii, dsc_jj, dsc_kk

  integer :: num_init
  logical :: program_first_time

  integer :: dim_ii, read_status
!!  character(256) :: name
!!  character(256) :: name_in
  character(256) :: adjusted_name
  character(256) :: temp_name
  integer :: len_sq
  real*8  :: len(3)
  real*8  :: len_right
  real*8  :: min_isle(3)
  integer :: temp_factor
  integer :: min_factor
  integer :: num_atoms
  real*8  :: area
  character(256) :: line
  integer :: result
  character(256) :: ignore_string_one, ignore_string_two
  character(256) :: ignore_string_three, ignore_string_four
  integer :: dsc_denom
  integer :: dsc_num_in(3)
  integer :: pbv_num
  real*8  :: csl_x_step
  integer :: csl_by_pp    !! pp is left
  integer :: csl_by_qq
  integer :: csl_by_dsc
  integer :: csl_by_lcm
  integer :: dsc_xx_in_steps(3)
  integer :: dsc_xx
  real*8  :: dsc_xx_adj
  real*8  :: xx_boundary_adj

  real*8  :: r_e_2
  real*8  :: r_dist
  real*8  :: gb_energy_out

  real*8  :: gb_best_energy(3)
  real*8  :: gb_best_dist(3)
  real*8  :: gb_best_place(3)
  integer :: gb_best_dsc(3, 3)
  real*8  :: best_energy
  real*8  :: best_dist
  real*8  :: best_place
  integer :: best_dsc(3)
  integer :: best_method

  integer :: method_ii
  integer :: place_index

  run_out_unit = 0

  call get_gb_energy_input_dsc(.false., 1)      !! single_crystal_flag
  if( alat .eq. -1.0d0 ) then
     write(0,*) "alat missing or -1"
     call exit(1)
  end if

  if( adj_name_flag ) then
     adjusted_name = trim(adj_name)
  else
     write(0,*) "adj_name  missing"
     call exit(1)
  end if
  
  if( do_now_flag ) then
     temp_name = trim(main_name) // "." // trim(adjusted_name) // ".energy"

     line = "cat temp.min.in | sed ""s/Q/temp."     &
            // trim(main_name)                      &
            // "."                                  &
            // trim(adjusted_name)                  &
            // "/"" > temp."                        &
            // trim(main_name)                      &
            // "."                                  &
            // trim(adjusted_name)                  &
            // ".min.in"
     result = system(trim(line)//char(0))
     if( result .ne. 0 ) then
        write(0,*) 'gb_energy: do_now_flag: system->sed failed'
        call abort()
     end if
  else
     temp_name = trim(main_name)//"."//trim(adjusted_name)//".task_list"
     open(task_unit, file=temp_name, status="new")
     temp_name = trim(main_name)//"."//trim(adjusted_name)//".control"
  end if

  open(6, file=trim(temp_name), status="new")

  call initialize_build_gb_control(init)

  init%fixed_x_len = fixed_x_len

  init%island_flag = island_flag
  if( island_flag ) then
     if( island_cut_type .eq. "circle" ) then
        init%island_cut_type = island_circle
     else
        write(0,*) "island_cut_type unknown  ", trim(island_cut_type)
        call exit(1)
     end if
     if( init%island_cut_type .eq. island_circle ) then
        init%island_radius = island_radius
     end if
     init%island_min_in_plane = island_min_in_plane
     init%island_min_thickness = island_min_thickness
  end if

  init%alat = alat
  init%crystal = crystal
  init%concentration = 0.0d0
  temp_name = "temp." // trim(main_name) // "." // trim(adjusted_name)
  init%name = trim(temp_name)
  init%min_dist = min_dist
  init%low_min_dist = low_min_dist
  init%e_per_a = e_per_a
  init%min_dist_flag = .true.
  init%crystal_to_delete = 1
  init%method = fixed_crystal
  init%two_boundaries = .false.
  init%output_flag = output_flag
  if( output_flag ) then
     init%lammps_flag = .false.
  else
     init%lammps_flag = lammps_flag
  end if
  init%one_at_a_time = one_at_a_time
  init%edge_size = edge_size
  init%nbr_delta = nbr_delta
  init%lammps_control_file = trim(lammps_control_file)
  init%delete_by_energy = delete_by_energy
  init%show_pressure = show_pressure
  init%delete_by_pressure = delete_by_pressure
  init%deletion_control_file = trim(deletion_control_file)
  init%lammps_command = trim(lammps_command)
  init%scr_file = trim(scr_file)
  init%sub_sh_file = trim(sub_sh_file)

  init%dsc_flag = .true.

  program_first_time = .true.

  read_loop:   &
  do
     read(5,*,iostat=read_status) ignore_string_one, pbv_num
     if( read_status .gt. 0 ) then
        write(0,*) 'read error'
        call abort()
     else if( read_status .lt. 0 ) then
        exit read_loop
     end if
     write(6,*) 'PBV:  ', pbv_num

!!     call flush(6)
!!     call for_flush(6)

!!      write(0,*) 'PBV:  ', pbv_num

     do dim_ii = 1, 3
        read(5,*) init%left_pbv_ii(dim_ii, 1:3), init%right_pbv_ii(dim_ii, 1:3)
        write(6,"(3i4,'  ', 3i4)") init%left_pbv_ii(dim_ii, 1:3),                &
                                   init%right_pbv_ii(dim_ii, 1:3)
     end do

     if( island_flag ) then

!! The lengths for the cell for island grain are given by the left pbv (the "matrix grain")
        min_isle(1) = island_min_in_plane
        min_isle(2) = island_min_in_plane
        min_isle(3) = island_min_thickness
        do dim_ii = 1, 3
           len_sq = dot_product( init%left_pbv_ii(dim_ii, 1:3), init%left_pbv_ii(dim_ii, 1:3))
           len(dim_ii) = sqrt( dble(len_sq) )
           min_factor = ceiling( min_isle(dim_ii) / len(dim_ii) )
           init%left_pbv_ii(dim_ii, 1:3) =  min_factor * init%left_pbv_ii(dim_ii, 1:3)
        end do

!! debugging 
        do dim_ii = 1, 3
           write(0,"(3i4,'  ', 3i4)") init%left_pbv_ii(dim_ii, 1:3),         &
                                      init%right_pbv_ii(dim_ii, 1:3)
        end do
!!$        stop
!! end debugging

        num_init = 0

!! To start, call do_one_placement directly

        dsc_vec = 0.0d0
        dsc_xx_adj = 0.0d0
                !!        xx_boundary_adjust = 0.0d0
        dsc_int_vec = 0
        place_index = 0

        call do_one_placement( island_gb_energy_result, island_gb_dist )

!! WORKING HERE
        write(0,*) "after do_one_placement"
        stop

     else     !! if( island_flag )

        read(5,*)   !! blank line
        read(5,*) ignore_string_one, csl_x_step
        read(5,*) ignore_string_one, csl_by_pp,           &
                  ignore_string_two, csl_by_qq,           &
                  ignore_string_three, csl_by_dsc,        &
                  ignore_string_four, csl_by_lcm
        read(5,*) ignore_string_one,  dsc_xx_in_steps

!! debugging
!!     write(6,*) "csl_x_step  ", csl_x_step
!!     write(6,*) "csl_by_pp   ", csl_by_pp
!!     write(6,*) "csl_by_qq   ", csl_by_qq
!!     write(6,*) "csl_by_dsc  ", csl_by_dsc
!! end debugging

        csl_x_step = abs(csl_x_step)
        csl_by_pp = abs(csl_by_pp)
        csl_by_qq = abs(csl_by_qq)
        csl_by_dsc = abs(csl_by_dsc)

        read(5,*)   !! "dsc"
        do dim_ii = 1, 3
           read(5,*)   !! dsc as rational line
        end do
        read(5,*)  !! blank line

        read(5,*)  ignore_string_one, dsc_denom
        do dim_ii = 1, 3
           read(5,*)  dsc_num_in
           dsc_mat(dim_ii,1:3) = dsc_num_in / dble( dsc_denom )
        end do


!! The lengths in the y and z directions are the same right and left
        do dim_ii = 1, 3
           len_sq = dot_product( init%left_pbv_ii(dim_ii, 1:3), init%left_pbv_ii(dim_ii, 1:3))
           len(dim_ii) = sqrt( dble(len_sq) )
        end do
        min_factor = ceiling( min_len_left / len(1) )
        if( require_two_repeats ) then
           temp_factor = max( 2, min_factor )
        else
           temp_factor = min_factor
        end if
        init%left_pbv_ii(1, 1:3) = temp_factor * init%left_pbv_ii(1, 1:3)
        do dim_ii = 2, 3
           min_factor = ceiling( min_width / len(dim_ii) )
           if( require_two_repeats ) then
              min_factor = max( min_factor, 2 )
           end if
           init%left_pbv_ii(dim_ii, 1:3) = min_factor * init%left_pbv_ii(dim_ii, 1:3)
           init%right_pbv_ii(dim_ii, 1:3) = min_factor * init%right_pbv_ii(dim_ii, 1:3)
        end do
     
        len_sq = dot_product( init%right_pbv_ii(1, 1:3), init%right_pbv_ii(1, 1:3))
        len_right = sqrt( dble(len_sq) )
        min_factor = ceiling( min_len_right / len_right )
        if( require_two_repeats ) then
           temp_factor = max( 2, min_factor )
        else
           temp_factor = min_factor
        end if
        init%right_pbv_ii(1, 1:3) = temp_factor * init%right_pbv_ii(1, 1:3)

        do dim_ii = 1, 3
           write(0,"(3i4,'  ', 3i4)") init%left_pbv_ii(dim_ii, 1:3),                &
                                      init%right_pbv_ii(dim_ii, 1:3)
        end do

        if( fixed_x_len ) then
           init%left_x_len = min_len_left
           init%right_x_len = min_len_right
        end if

        num_init = 0

        call loop_over_dsc_vectors()

     end if   !! else of if( island_flag )

     if( do_now_flag ) then
        do method_ii = 1, 3
           if( do_method(method_ii) ) then
              write(6,"('PBV ', i6, '  Method ', i2, 3i6, f8.3, f8.4, f9.3)")    &
                   pbv_num, method_ii, gb_best_dsc(1:3, method_ii),               &
                   gb_best_place(method_ii),                                      &
                   gb_best_dist(method_ii),                                       &
                   gb_best_energy(method_ii)
              if( program_first_time ) then
                 program_first_time = .false.
                 best_energy = gb_best_energy(method_ii)
                 best_dist = gb_best_dist(method_ii)
                 best_place = gb_best_place(method_ii)
                 best_dsc = gb_best_dsc(1:3, method_ii)
                 best_method = method_ii
              else
                 if( gb_best_energy(method_ii) < best_energy ) then
                    best_energy = gb_best_energy(method_ii)
                    best_dist = gb_best_dist(method_ii)
                    best_place = gb_best_place(method_ii)
                    best_dsc = gb_best_dsc(1:3, method_ii)
                    best_method = method_ii
                 end if
              end if
           end if   !! if( do_method(method_ii)
        end do      !! method_ii
        write(6,*)
        write(6,"('PBV ', i6, ' best ', i2, 3i6, f8.3, f18.14, f9.3)")    &
             pbv_num, best_method, best_dsc,                          &
             best_place, best_dist, best_energy
     else
        write(6,*) "num_init  ", num_init
     end if     !! else of if( do_now_flag )

     read(5,*)   !! blank line

  end do read_loop

  if( do_now_flag ) then
     temp_name = "temp." // trim(main_name) // "." // trim(adjusted_name) // ".min.in"
     call delete_file(trim(temp_name), len_trim(temp_name))
  end if

  stop

contains

  subroutine do_one_placement( gb_energy_result, gb_dist )
    implicit none

    real*8, intent(out) :: gb_energy_result(3)
    real*8, intent(out) :: gb_dist(3)       !! minimum nbr distance used (i.e. best)

    integer :: mm_ii
    character(256) :: init_file_name, name_from_init, min_in_name

!! Try the three methods
    do mm_ii = 1, 3
       if( do_method(mm_ii) ) then
          if( mm_ii .eq. 1 ) then
             init%method = fixed_crystal
             init%crystal_to_delete = 1
          else if( mm_ii .eq. 2 ) then
             init%method = fixed_crystal
             init%crystal_to_delete = 2
          else if( mm_ii .eq. 3 ) then
             init%method = delete_by_averaging
             init%crystal_to_delete = -1
          else
             write(0,*) "gb_energy_cc_diwo:  internal error  mm_ii"
             call abort()
          end if

          init%dsc_vec = dsc_vec
          init%dsc_xx_adj = dsc_xx_adj
          if( island_flag ) then
             init%xx_boundary_adj = 0.0d0
          else
             init%xx_boundary_adj = xx_boundary_adj
          end if
          init%dsc_int_vec = dsc_int_vec
          init%method_num = mm_ii
          init%place_index = place_index

          if( do_now_flag ) then 
!! If output_flag is true, build_gb_dsc will not return.
             call build_gb_dsc(init,                                 & 
                               num_atoms, area,                      &
                               r_e_2, r_dist,                        &
                               .false.)             !! no xi/chi files, thank you
             gb_energy_out = r_e_2 - e_per_a * num_atoms
             gb_energy_out = 1000.0 * gb_energy_out / (area )
             gb_energy_result(mm_ii) = gb_energy_out
             gb_dist(mm_ii) = r_dist
             write(6,"('GB energy(',i1,')  ', 4i6, f8.3, f8.4, f9.3)")      &
                  mm_ii, pbv_num, dsc_int_vec, xx_boundary_adj, r_dist,       &
                  gb_energy_out

!!             call flush(6)
!!             call for_flush(6)

             write(0,"('GB energy(',i1,')  ', 4i6, f8.3, f8.4, f9.3)")      &
                  mm_ii, pbv_num, dsc_int_vec, xx_boundary_adj, r_dist,       &
                  gb_energy_out
          else
             init%main_name = main_name
             init%pbv_num = pbv_num
             init%dsc_num = dsc_num
             call write_gb_init( init, name_from_init )
             init_file_name = trim(name_from_init)//".init"
             write(task_unit, "(a)") trim(init_file_name)
             min_in_name = "temp."//trim(name_from_init)//".min.in"

!! debugging
!!$             write(0,*) "gb_energy_cc_diwo:  lammps_control_file"
!!$             write(0,*) trim(lammps_control_file)
!! end debugging

             line = "cat "//trim(lammps_control_file)                              &
                    //" | sed ""s/Q/"//"temp."//trim(name_from_init)               &
                    //"/"" > "//trim(min_in_name)

!! debugging
!!             write(0,*) "gb_energy_cc_diwo:  line (for system(line))"
!!             write(0,"(a)") line
!! end debugging


             result = system(trim(line)//char(0))
             if( result .ne. 0 ) then
                write(0,*) 'gb_energy_cc_diwo: do_one_placement: system->sed failed'
                call abort()
             end if
             num_init = num_init + 1
          end if
       end if
    end do

    return
  end subroutine do_one_placement


  subroutine loop_over_boundary_placements(best_gb_energy,         &
                                           best_gb_dist,           &
                                           best_placement)
    implicit none

    real*8, intent(out)  :: best_gb_energy(3)
    real*8, intent(out)  :: best_gb_dist(3)
    real*8, intent(out)  :: best_placement(3)

    logical, allocatable :: plane_exists(:)
    integer :: csl_by_dsc_step
    integer :: our_lcm
    integer :: lcm_test_ii, step_ii
    integer :: test_lcm
    integer :: total_steps
    integer :: qq_steps, pp_steps
    integer :: dsc_temp, dsc_steps
    integer :: qq_plane, pp_plane
    integer :: low_plane, high_plane
    integer :: ii
    real*8  :: average_plane
    real*8  :: gb_energy_one_placement(3)
    real*8  :: gb_dist_one_placement(3)
    logical :: first_time

    csl_by_dsc_step = csl_by_dsc * dsc_num
    if( dsc_half_flag ) then
       csl_by_dsc_step = 2 * csl_by_dsc_step
    end if
    our_lcm = -1
    lcm_loop:  &
    do lcm_test_ii = 1, 48
       test_lcm = lcm_test_ii * csl_by_lcm
       if( modulo(test_lcm, csl_by_dsc_step) .eq. 0 ) then
          our_lcm = test_lcm
          exit lcm_loop
       end if
    end do lcm_loop
    if( our_lcm .eq. -1 ) then
       write(0,*) "loop_over_boundary_placements  did not find an lcm, extend limit?"
       call abort()
    end if

    total_steps = our_lcm
    allocate( plane_exists(0:total_steps) )
    plane_exists = .false.

    qq_steps = total_steps / csl_by_qq
    pp_steps = total_steps / csl_by_pp
    dsc_temp = total_steps / csl_by_dsc_step
    dsc_steps = dsc_xx * dsc_temp
    if( dsc_half_flag ) then
       dsc_steps = dsc_steps + dsc_temp / 2
    end if

    do step_ii = 0, csl_by_qq
       qq_plane = step_ii * qq_steps
       plane_exists( qq_plane ) = .true.
    end do
    do step_ii = 1, csl_by_pp
       pp_plane = step_ii * pp_steps + dsc_steps
       pp_plane = modulo( pp_plane, total_steps )
       plane_exists( pp_plane ) = .true.
    end do

    !! Note that plane_exists(0) and plane_exists(total_steps) are true,
    !! because of the qq_plane loop.

    place_index = 0  !! for writing out init for "do later"
    first_time = .true.
    low_plane = 0
    outer_loop:  &
    do
       high_plane = low_plane + 1
       inner_loop:  &
       do
          if( high_plane > total_steps ) then
             exit outer_loop   !! only happens when low_plane .eq. total_steps
          end if
          if( plane_exists( high_plane ) ) then  !! loop until it exists, or is too big
             exit inner_loop
          else
             high_plane = high_plane + 1  
          end if
       end do inner_loop
       !! Here we have high_plane > low_plane, and both exist.
       average_plane = ( dble(low_plane) + dble(high_plane) ) / 2.0d0
       xx_boundary_adj = average_plane * csl_x_step / dble(total_steps)

       place_index = place_index + 1
       call do_one_placement( gb_energy_one_placement, gb_dist_one_placement )
       if( do_now_flag ) then
          if( first_time ) then
             first_time = .false.
             best_gb_energy = gb_energy_one_placement
             best_gb_dist = gb_dist_one_placement
             best_placement = xx_boundary_adj
          else
             do ii = 1, 3
                if( do_method(ii) ) then
                   if( gb_energy_one_placement(ii) < best_gb_energy(ii) ) then
                      best_gb_energy(ii) = gb_energy_one_placement(ii)
                      best_gb_dist(ii)   = gb_dist_one_placement(ii)
                      best_placement(ii) = xx_boundary_adj
                   end if
                end if
             end do
          end if
       end if

       low_plane = high_plane
    end do outer_loop

    if( do_now_flag ) then
       do ii = 1, 3
          if( do_method(ii) ) then
             write(6,"('PBV ', i6, '  Method ', i2, '  DSC ', 3i6, f8.3, f8.4, f9.3)")    &
                  pbv_num, ii, dsc_int_vec,    &
                  best_placement(ii), best_gb_dist(ii), best_gb_energy(ii)
          end if
       end do
    end if

    deallocate( plane_exists )

    return
  end subroutine loop_over_boundary_placements


  subroutine loop_over_dsc_vectors()
    implicit none

    logical :: first_time
    real*8 :: gb_energy_one_dsc(3)
    real*8 :: gb_dist_one_dsc(3)
    real*8 :: gb_place_one_dsc(3)


    first_time = .true.
    do dsc_ii = 1, dsc_num
       dsc_int_vec(1) = dsc_ii - 1
       do dsc_jj = 1, dsc_num
          dsc_int_vec(2) = dsc_jj - 1
          do dsc_kk = 1, dsc_num
             dsc_int_vec(3) = dsc_kk - 1
             call compute_dsc_vec( dsc_mat,       &
                                   dsc_int_vec,   &
!!   now global                                   dsc_num,       &
                                   dsc_vec )
             if( dsc_half_flag ) then
                dsc_xx_adj = csl_x_step / ( 2 * csl_by_dsc * dsc_num )
             else
                dsc_xx_adj = 0.0d0
             end if
             dsc_xx = dot_product( dsc_xx_in_steps, dsc_int_vec )
             call loop_over_boundary_placements(gb_energy_one_dsc,   &
                                                gb_dist_one_dsc,     &
                                                gb_place_one_dsc)
             if( do_now_flag ) then
                if( first_time ) then
                   first_time = .false.
                   gb_best_energy = gb_energy_one_dsc
                   gb_best_dist = gb_dist_one_dsc
                   gb_best_place = gb_place_one_dsc
                   do method_ii = 1, 3
                      gb_best_dsc(1:3, method_ii) = dsc_int_vec
                   end do
                else
                   do method_ii = 1, 3
                      if( do_method(method_ii) ) then
                         if( gb_energy_one_dsc(method_ii) < gb_best_energy(method_ii) ) then
                            gb_best_energy(method_ii) = gb_energy_one_dsc(method_ii)
                            gb_best_dist(method_ii) = gb_dist_one_dsc(method_ii)
                            gb_best_place(method_ii) = gb_place_one_dsc(method_ii)
                            gb_best_dsc(1:3, method_ii) = dsc_int_vec
                         end if
                      end if
                   end do
                end if    !! else of if( first_time )
             end if       !! if( do_now_flag )
          end do
       end do
    end do

    return
  end subroutine loop_over_dsc_vectors


end program gb_energy_cc_diwo
