module build_gb_boundary_mod
  implicit none

  integer, parameter :: max_passes = 100
  integer, parameter :: max_one_at_a_time = 500

contains


!!
!! Warning, returns from within code if method = no_deletion
!!
!! For each boundary, global_current_min must be set to zero before
!! calling enforce_min_dist for the first time.  (If you use
!! global_current_min.)
!!
  subroutine enforce_min_dist(one_pass_flag)
    use build_gb_param_mod
    use build_gb_global_mod
    use gb_neighbor_mod
    use build_gb_run_out_mod
    use smallish_real_8_mod
    implicit none
    logical, intent(in) :: one_pass_flag

    integer, save :: num_passes

    real*8  :: shortest_nbr_dist, current_min
    logical :: deletion_finished, no_nbr_flag
    integer :: pass_high

    if( .not. one_pass_flag ) then
       num_passes = 0
       pass_high = max_passes
       deletion_finished = .false.
    else
       pass_high = num_passes + 1
       deletion_finished = .true.    !! we just do one pass at a time
    end if

!! debugging
    write(run_out_unit,*) 'enforce: ', one_pass_flag, num_passes, pass_high
!! end debugging

    deletion_loop: &
    do num_passes = num_passes + 1, pass_high
!! debugging
       write(run_out_unit,*) 'nbr deletion pass: ', num_passes
!!       write(82,*) 'nbr deletion pass: ', num_passes
!!       write(83,*) 'nbr deletion pass: ', num_passes
!! end debugging
       call build_nbr_list(full_min_dist)
       if( num_passes .eq. 1 ) then
          if( two_boundaries .and. sigma_flag ) then
             call get_sigma_number()   !! without two_boundaries, no go
          else
             sigma_number = -1
          end if
          if( method .eq. no_deletion ) then

!! debugging
!!             write(run_out_unit,*) "return from within deletion_loop; method is no_deletion"
!! end debugging

             return        !! nothing for us to do
          end if
       end if

!! debugging
!!       write(run_out_unit,*) "before find_shortest_nbr call"
!! end debugging

       call find_shortest_nbr( no_nbr_flag, shortest_nbr_dist )
       if( no_nbr_flag ) then
          deletion_finished = .true.
          lammps_loop_finished = .true.
          exit deletion_loop       !! no more to do
       end if
       if( (.not. one_pass_flag) .and.(shortest_nbr_dist .gt. low_min_dist) ) then

!! Not valid, could be second time through loop, etc !! global_current_min = 0.0

          deletion_finished = .true.

!! debugging
!!          write(run_out_unit,*) "exiting_deletion loop; shortest_nbr_dist .gt. low_min_dist"
!! end debugging

          if( shortest_nbr_dist .gt. min_dist ) then
             lammps_loop_finished = .true.

!! debugging
!!          write(run_out_unit,*) "  and  shortest_nbr_dist .gt. min_dist"
!! end debugging

          end if

          exit deletion_loop
       end if
       if( one_at_a_time ) then
          current_min = shortest_nbr_dist + smallish_real_8
       else
          current_min = shortest_nbr_dist + delta
       end if
       current_min = min( current_min, min_dist )  !! Strictly speaking, no effect
!!        global_current_min = shortest_nbr_dist
       global_current_min = current_min

!! debugging
       write(run_out_unit,*) "global_current_min  ", global_current_min
!! end debugging

       call delete_nbrs( current_min, (one_at_a_time .and. one_pass_flag) )
    end do deletion_loop


    if( .not. deletion_finished ) then
       write(run_out_unit,*) 'enforce_min_dist: overran max_passes'
       call abort()
    end if
    num_passes = num_passes - 1          !! Thd do loop increments past pass_high
!! debugging
    write(run_out_unit,*) 'delete_nbrs: num_passes: ', num_passes
!! end debugging
!! debugging
    write(run_out_unit,*) 'enforce: shortest_nbr_dist', shortest_nbr_dist
    write(run_out_unit,*) 'global_current_min  ', global_current_min
!! end debugging

    return
  end subroutine enforce_min_dist


end module build_gb_boundary_mod
