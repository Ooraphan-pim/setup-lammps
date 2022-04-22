module gb_neighbor_mod
!!  use fcc_mod
  use build_gb_atom_data_mod
  use build_gb_run_out_mod
  implicit none

!! Global debugging parameter to allow analyze_dump_file to control which config to debug.

  logical :: global_neighbor_debug_flag = .false.
  logical :: global_show_nbrs_flag = .false.

!! Global parameters (intended to be private to neighbor_mod)

!!  integer, parameter :: max_nbrs = 20   !! now in atoms.f90
  integer, parameter :: max_bins = 1600000
  integer, parameter :: max_depth = 100

  real*8, parameter :: fraction = 7.0d0/11.0d0

!!   logical, parameter :: full_nbr_list = .false. !!nbr pairs appear in both atoms lists

!! Types (intended to be private to neighbor_mod)

!!$  type nbr_entry_type
!!$     real*8 :: dist_sq
!!$!!     real*8 :: xi_chi_error(2)
!!$  end type nbr_entry_type

!! Global variables (intended to be private to neighbor_mod)

  real*8  :: max_nbr_dist
  real*8  :: max_nbr_dist_sq
!!    real*8  :: max_error            !! the square of a distance

  integer :: num_slices(3)
  real*8  :: slice_len(3)
  integer :: num_bins


  integer :: atom_bin(3, max_atoms, max_crystals)

  real*8  :: atom_rel_pos(3, max_atoms, max_crystals)    !! position, but put into box

  integer :: bin_head(max_bins)

!!   integer :: next_atom(max_atoms)               !! in atoms.f90
!!   integer :: num_nbrs(max_atoms)

!!$  type(nbr_entry_type) :: nbr_list(max_nbrs, max_atoms) 

!! debugging
  integer :: high_bin_count
!! end debugging

!! debugging/analysis
!!  real*8 :: max_actual_dist_sq
!!

contains


!! debugging subroutine
!***************************************************************************
!                   show_nbr_pairs                   
!***************************************************************************

  subroutine show_nbr_pairs()
    use build_gb_global_mod
    implicit none

    integer atom_ii, nbr_ii, nbr_atom_ii

    do atom_ii = 1, num_atoms_built
       do nbr_ii = 1, atom(atom_ii)%num_nbrs
          nbr_atom_ii = atom(atom_ii)%nbr_list(nbr_ii)%nbr_index
          write(82,"(i10,i4,5x,i10,i4)")                               &
                   atom_ii, atom(atom_ii)%crystal,                     &
                   nbr_atom_ii, atom(nbr_atom_ii)%crystal
          write(83,"(3f8.2,5x,3f8.2)") atom(atom_ii)%pos,              &
                                       atom(nbr_atom_ii)%pos
          if( .not. full_nbr_list ) then   !! print reverse as well
             write(82,"(i10,i4,5x,i10,i4)")                            &
                      nbr_atom_ii, atom(nbr_atom_ii)%crystal,          &
                      atom_ii, atom(atom_ii)%crystal
             write(83,"(3f8.2,5x,3f8.2)") atom(nbr_atom_ii)%pos,       &
                                          atom(atom_ii)%pos
          end if
       end do
    end do

    return
  end subroutine show_nbr_pairs
!! end debugging subroutine



!! debugging subroutine
!***************************************************************************
!                   show_num_nbrs
!***************************************************************************

  subroutine show_num_nbrs()
    use build_gb_global_mod
    implicit none

    integer :: max_num_nbrs

    if( full_nbr_list ) then
       max_num_nbrs = maxval( atom(1:num_atoms_built)%num_nbrs )
       call sub_show_num_nbrs()
    end if

    return

  contains

    subroutine sub_show_num_nbrs
      implicit none

      integer :: atoms_with_num_nbrs(0:max_num_nbrs)
      integer :: atom_ii, temp_num_nbrs, total_atoms
      integer :: ii

!! debugging debugging
!!$      integer :: debug_index, debug_jj
!!$      real*8  :: debug_pos_diff(3)
!!$      real*8  :: debug_sq_dist
!! end debugging debugging


      atoms_with_num_nbrs = 0
      do atom_ii = 1, num_atoms_built
         temp_num_nbrs = atom(atom_ii)%num_nbrs

!! debugging debugging
!!$         if( temp_num_nbrs .eq. 0 ) then
!!$            write(run_out_unit,*) 'no neighbors'
!!$            write(run_out_unit,*) 'atom_ii: ', atom_ii
!!$            write(run_out_unit,*) 'atom_rel_pos(1:3, atom_ii): ', atom_rel_pos(1:3, atom_ii)
!!$            write(run_out_unit,*) 'atom_pos(1:3, atom_ii): ', atom_pos(1:3, atom_ii)
!!$            write(run_out_unit,*) 'atom, distance'
!!$            do debug_index = 1, num_atoms_built
!!$               if( debug_index .eq. atom_ii ) cycle
!!$               debug_pos_diff =                                           &
!!$                  atom_rel_pos(1:3,atom_ii) - atom_rel_pos(1:3,debug_index)
!!$               do debug_jj = 1, 3
!!$                  debug_pos_diff(debug_jj) = debug_pos_diff(debug_jj)   &
!!$                    - perlen(debug_jj)                                  &
!!$                         * nint( debug_pos_diff(debug_jj) / perlen(debug_jj) )
!!$               end do
!!$               debug_sq_dist = dot_product(debug_pos_diff, debug_pos_diff)
!!$               if( sqrt(debug_sq_dist) .lt. max_nbr_dist ) then
!!$                  write(run_out_unit,"(i5, 4f14.4)") debug_index, sqrt(debug_sq_dist), atom_rel_pos(1:3, debug_index)
!!$               end if
!!$            end do
!!$            stop
!!$         end if
!! end debugging debugging

         atoms_with_num_nbrs(temp_num_nbrs) =                            &
              atoms_with_num_nbrs(temp_num_nbrs) + 1
      end do

      write(81,*)
      write(81,*) 'atoms by number of neighbors'
      write(81,*) 'max nbr distance: ', max_nbr_dist 
      total_atoms = 0
      do ii = 0, max_num_nbrs
         write(81,*) ii, atoms_with_num_nbrs(ii)
         total_atoms = total_atoms + atoms_with_num_nbrs(ii)
      end do
      write(81,*)
      if( total_atoms .ne. num_atoms_built ) then
         write(run_out_unit,*) 'problem: num_atoms_built, total_atoms: ', num_atoms_built, total_atoms
      else
         write(81,*) 'total_atoms: ', total_atoms
      end if

      return
    end subroutine sub_show_num_nbrs

  end subroutine show_num_nbrs
!! end debugging subroutine


!***************************************************************************
!                   find_shortest_nbr
!***************************************************************************

  subroutine find_shortest_nbr( no_nbr_flag, shortest_dist )
    implicit none

    logical, intent(out) :: no_nbr_flag      !! set to true if no nbr pairs at all
    real*8, intent(out)  :: shortest_dist

    integer :: atom_ii, nbr_ii
    real*8  :: shortest_dist_sq

    no_nbr_flag = .true.

    atom_loop: &
    do atom_ii = 1, num_atoms_built
       if( atom(atom_ii)%deleted_atom ) then
          cycle atom_loop
       end if
       do nbr_ii = 1, atom(atom_ii)%num_nbrs
          if( no_nbr_flag ) then
             no_nbr_flag = .false.
             shortest_dist_sq = atom(atom_ii)%nbr_list(nbr_ii)%dist_sq
          else
             shortest_dist_sq =                                                   &
                  min( shortest_dist_sq, atom(atom_ii)%nbr_list(nbr_ii)%dist_sq )
          end if
       end do
    end do atom_loop

    if( .not. no_nbr_flag ) then
       if( shortest_dist_sq .lt. 1.0d-8 ) then
          shortest_dist = 0.0
       else
          shortest_dist = sqrt( shortest_dist_sq ) 
       end if

!! debugging/analysis
       write(run_out_unit,*) 'find_shortest_nbr: ', shortest_dist
!! end debugging/analysis
    else
!! debugging
       write(run_out_unit,*) 'find_shortest_nbr; no_nbr_flag'
!! end debugging

    end if

    return

  end subroutine find_shortest_nbr



!***************************************************************************
!                   delete_nbrs
!***************************************************************************
!
! This is for gb.  Neighbors are atoms that are too close and need deletion.
! Two methods are available, either atoms are delete from a specific crystal, 
! or we delete both atoms and replace them by one at the average position.
!

  subroutine delete_nbrs( min_dist_allowed, one_atom_flag )
    use build_gb_global_mod
    implicit none

    real*8, intent(in) :: min_dist_allowed
    logical, intent(in) :: one_atom_flag

    integer, save :: prior_atom_ii = -1     !! used if one atom at a time

    integer :: atom_ii, nbr_ii
    real*8  :: min_allowed_sq
    integer :: counter, atom_ii_start

    min_allowed_sq = min_dist_allowed * min_dist_allowed

    if( one_atom_flag ) then
       if( prior_atom_ii .eq. -1 ) then
          atom_ii_start = 1
       else
          atom_ii_start = prior_atom_ii + fraction * num_atoms_built
          if( atom_ii_start .gt. num_atoms_built ) then
             atom_ii_start = atom_ii_start - num_atoms_built
          end if
       end if
    else
       atom_ii_start = 1
    end if
    atom_ii = atom_ii_start - 1

    atom_loop: &
    do counter = 1, num_atoms_built
       atom_ii = atom_ii + 1
       if( atom_ii .gt. num_atoms_built ) then
          atom_ii = atom_ii - num_atoms_built
       end if
       if( atom(atom_ii)%deleted_atom ) then
          cycle atom_loop
       end if
       do nbr_ii = 1, atom(atom_ii)%num_nbrs
          if( atom(atom_ii)%nbr_list(nbr_ii)%dist_sq < min_allowed_sq ) then
             call delete_one_atom()
             if( one_atom_flag ) then
                prior_atom_ii = atom_ii
                exit atom_loop            !! we deleted one
             else
                cycle atom_loop            !! no longer the same atom, nbrs not defined
             end if
          end if
       end do
    end do atom_loop

    !! debugging
    write(run_out_unit,*) "actual_num_atoms  ", actual_num_atoms

    return

  contains

!! Warning, returns early if the nbr atom was already deleted during this pass.
!! (If we are only deleting one per pass, that is not a problem.)
    subroutine delete_one_atom()
      use build_gb_param_mod
      implicit none

      integer :: dim_ii, nbr_atom, crystal_one, crystal_two
      integer :: atom_to_delete
      real*8  :: pos_one(3), pos_two(3), rel_pos_two(3), ave_pos(3)

      nbr_atom = atom(atom_ii)%nbr_list(nbr_ii)%nbr_index

      if( atom(nbr_atom)%deleted_atom ) then
!! This is a situation where an atom with lower atom_ii had this as a nbr and deleted
!! it already.  We cannot do anything.
         return
      end if

      if( method .eq. delete_by_averaging ) then
         pos_one = atom(atom_ii)%pos
         pos_two = atom(nbr_atom)%pos
         rel_pos_two = pos_two - pos_one
         do dim_ii = 1, 3
            rel_pos_two(dim_ii) = rel_pos_two(dim_ii)                             &
                 - perlen(dim_ii) * nint( rel_pos_two(dim_ii) / perlen(dim_ii) )
         end do
         ave_pos = pos_one + 0.5d0 * rel_pos_two
!! Put the average atom over atom_one, delete atom_two
!! atom_one is now a new atom, its nbrs are undefined.
         atom(nbr_atom)%deleted_atom = .true.
         atom(atom_ii)%num_nbrs = 0
         atom(atom_ii)%pos = ave_pos
         atom(atom_ii)%crystal = 3
      else if( method .eq. fixed_crystal ) then
         crystal_one = atom(atom_ii)%crystal
         crystal_two = atom(nbr_atom)%crystal
         if( crystal_one .eq. crystal_two ) then
            write(run_out_unit,*) 'delete_one_atom: both are same crystal'
            call abort()
         end if
         if( crystal_one .eq. crystal_to_delete ) then
            atom_to_delete = atom_ii
         else if( crystal_two .eq. crystal_to_delete ) then
            atom_to_delete = nbr_atom
         else
            write(run_out_unit,*) 'atom_to_delete: internal error: neither is crystal_to_delete'
            write(run_out_unit,*) '  atom_ii: ', atom_ii
            write(run_out_unit,*) '  nbr_atom: ', nbr_atom
            write(run_out_unit,*) '  crystal_one: ', crystal_one
            write(run_out_unit,*) '  crystal_two: ', crystal_two
            write(run_out_unit,*) '  crystal_to_delete: ', crystal_to_delete
            call abort()
         end if
         atom(atom_to_delete)%deleted_atom = .true.
!! debugging or information
!!         write(run_out_unit,"('deleted: ', 3f14.4)") atom(atom_to_delete)%pos
!! end debugging
      else
         write(run_out_unit,*) 'delete_one_atom: invalid method'
         call abort()
      end if

      num_atoms_deleted = num_atoms_deleted + 1
      actual_num_atoms = actual_num_atoms - 1

      return
    end subroutine delete_one_atom

  end subroutine delete_nbrs

    


!***************************************************************************
!                   build_nbr_list
!***************************************************************************
!
! We want to find each nbr pair once, then put the information in both
! atoms nbr lists.
!
! To find each nbr pair once, we go through each atom.  We look for neighbors
! that are further down the chain in the current bin.  And we look for
! neighbors in other bins with equal or higher bin_int numbers.
!



  subroutine build_nbr_list(max_nbr_dist_in)
    use build_gb_global_mod
    implicit none

    real*8, intent(in) :: max_nbr_dist_in
!!    real*8, intent(in) :: nn_dist_in

    integer :: atom_ii
    integer :: bin_ii, bin_jj, bin_kk, int_bin
    integer :: search_ii, search_jj, search_kk, int_search
!!     integer :: link_one, link_two
    integer :: atom_one, atom_two
    integer :: bin(3), search(3)
    integer :: adj_perlen(3)

!!    real*8  :: atom_one_pos(3), atom_two_pos(3)  !! from atom_rel_pos
!!    type(nbr_entry_type) :: nbr_entry
!!    logical :: are_they_nbrs

!! debugging
    integer :: depth
!! end debugging

!!    max_error = (max_nbr_dist_in - nn_dist_in)**2

!! bin_atoms calls setup_bins which sets max_nbr_dist to max_nbr_dist_in
    call bin_atoms(max_nbr_dist_in)

!! debugging
!!      call count_binning()
!!      stop
!! end if

    do atom_ii = 1, num_atoms_built
       atom(atom_ii)%num_nbrs = 0
    end do

    do bin_ii = 0, num_slices(1) - 1
       do bin_jj = 0, num_slices(2) - 1
          bin_kk_loop: &
          do bin_kk = 0, num_slices(3) - 1
             bin(1) = bin_ii
             bin(2) = bin_jj
             bin(3) = bin_kk
             call bin_to_int_bin(bin, int_bin)
             if( bin_head(int_bin) .eq. 0 ) then
                cycle bin_kk_loop
             end if
             atom_one = bin_head(int_bin)
!! debugging
             depth = 0
!! end debugging
             chain_one_loop: &
             do
!! debugging
                depth = depth + 1
                if( depth .gt. high_bin_count ) then
                   write(run_out_unit,*) 'internal error, chain_one_loop depth'
                   write(run_out_unit,*) "depth  ", depth
                   write(run_out_unit,*) "high_bin_count  ", high_bin_count
                   call abort()
                end if
!! end debugging
                call find_atom_nbrs()
                if( atom(atom_one)%next_atom .eq. 0 ) then
                   exit chain_one_loop
                else
                   atom_one = atom(atom_one)%next_atom
                end if
             end do chain_one_loop
          end do bin_kk_loop
       end do                 !! bin_jj 
    end do                    !! bin_ii

!! debugging
!!    if( global_neighbor_debug_flag .and. .true. ) then

    if( nbrs_by_dist ) then
       call show_num_nbrs()
    end if

    if( global_show_nbrs_flag ) then
       call show_nbr_pairs()
    end if

!!       stop
!!    end if
!! end debugging

    return

  contains

    subroutine find_atom_nbrs()
      implicit none

!! debugging
      integer :: depth_two
!! end debugging

!!      atom_one_pos = atom_rel_pos(1:3, atom_one)

!! first look for neighbors in this bin
      adj_perlen = 0

      if( atom(atom_one)%next_atom .ne. 0 ) then

         atom_two = atom(atom_one)%next_atom

!! debugging
         depth_two = 0
!! end debugging

         this_bin_chain_two_loop: &
         do
!! debugging
            depth_two = depth_two + 1
            if( depth_two .gt. high_bin_count ) then
               write(run_out_unit,*) 'internal error, this_bin_chain_two_loop depth'
               call abort()
            end if
!! end debugging

!!                 atom_two_pos = atom_rel_pos(1:3, atom_two)
            call check_are_they_nbrs(atom_one,                &
                                     atom_two,                &
                                     adj_perlen)

            if( atom(atom_two)%next_atom .eq. 0 ) then
               exit this_bin_chain_two_loop
            else
               atom_two = atom(atom_two)%next_atom
            end if
         end do this_bin_chain_two_loop

      end if   !! if atom(atom_one)%next_atom .ne. 0 

!! Now check the other bins
      search_ii_loop: &
      do search_ii = -1, 1
         search(1) = bin(1) + search_ii
         if( search(1) .ge. num_slices(1) ) then
            if( two_boundaries ) then
               search(1) = 0
               adj_perlen(1) = 1
            else
               cycle search_ii_loop
            end if
         else if( search(1) .eq. -1 ) then
            if( two_boundaries ) then
               search(1) = num_slices(1) - 1
               adj_perlen(1) = -1
            else
               cycle search_ii_loop
            end if
         else
            adj_perlen(1) = 0
         end if
         do search_jj = -1, 1
            search(2) = bin(2) + search_jj
            if( search(2) .ge. num_slices(2) ) then
               search(2) = 0
               adj_perlen(2) = 1
            else if( search(2) .eq. -1 ) then
               search(2) = num_slices(2) - 1
               adj_perlen(2) = -1
            else
               adj_perlen(2) = 0
            end if
            search_kk_loop: &
            do search_kk = -1, 1
               search(3) = bin(3) + search_kk
               if( search(3) .ge. num_slices(3) ) then
                  search(3) = 0
                  adj_perlen(3) = 1
               else if( search(3) .eq. -1 ) then
                  search(3) = num_slices(3) - 1
                  adj_perlen(3) = -1
               else
                  adj_perlen(3) = 0
               end if

               call bin_to_int_bin(search, int_search)

               if( int_search .le. int_bin ) then
                  cycle search_kk_loop    !! we already did our own bin
                                          !! and we will do lower bins
                                          !! with int_search as bin_int
               end if
               if( bin_head(int_search) .eq. 0 ) then
                  cycle search_kk_loop    !! empty bin
               end if

!! debugging
               depth_two = 0
!! end debugging

               atom_two = bin_head(int_search)
               chain_two_loop: &
               do
!! debugging
                  depth_two = depth_two + 1
                  if( depth_two .gt. high_bin_count ) then
                     write(run_out_unit,*) 'internal error: regular chain_two_loop depth'
                     call abort()
                  end if
!! end debugging

                  call check_are_they_nbrs(atom_one, atom_two, adj_perlen)

                  if( atom(atom_two)%next_atom .eq. 0 ) then
                     exit chain_two_loop
                  else
                     atom_two = atom(atom_two)%next_atom
                  end if
               end do chain_two_loop
            end do search_kk_loop
         end do                     !! search_jj
      end do search_ii_loop         !! search_ii

      return
    end subroutine find_atom_nbrs

  end subroutine build_nbr_list


!***************************************************************************
!                   check_are_they_nbrs
!***************************************************************************

  subroutine check_are_they_nbrs(atom_one, atom_two, adj_perlen)
    use build_gb_atom_data_mod
    use build_gb_global_mod
    implicit none

    integer, intent(in) :: atom_one
    integer, intent(in) :: atom_two
    integer, intent(in) :: adj_perlen(3)

    real*8  :: rel_pos(3)
    integer :: ii
    type(nbr_entry_type) :: nbr_entry
    real*8  :: dist_sq

!!    real*8  :: xi_error, chi_error

!!$    if( atom_one .eq. 0 ) then
!!$       write(run_out_unit,*) 'check_are_they_nbrs: atom_one is zero'
!!$       call abort()
!!$    end if
!!$    if( atom_two .eq. 0 ) then
!!$       write(run_out_unit,*) 'check_are_they_nbrs: atom_two is zero'
!!$       call abort()
!!$    end if

    if( atom(atom_one)%deleted_atom ) then
       write(run_out_unit,*) 'check_are_they_nbrs: atom_one is a deleted atom'
       call abort()
    end if
    if( atom(atom_two)%deleted_atom ) then
       write(run_out_unit,*) 'check_are_they_nbrs: atom_two is a deleted atom'
       call abort()
    end if

    rel_pos = atom(atom_two)%rel_pos - atom(atom_one)%rel_pos
    dist_sq = 0.0
    do ii = 1, 3
       rel_pos(ii) = rel_pos(ii) + adj_perlen(ii) * perlen(ii)
       dist_sq = dist_sq + rel_pos(ii)**2
    end do

    if( dist_sq .lt. max_nbr_dist_sq ) then
       nbr_entry%dist_sq = dist_sq
       atom(atom_one)%num_nbrs = atom(atom_one)%num_nbrs + 1
       if( atom(atom_one)%num_nbrs .gt. max_nbrs ) then
          write(run_out_unit,*) 'overran max_nbrs, atom_one'
          call abort()
       end if
       nbr_entry%nbr_index = atom_two
       atom(atom_one)%nbr_list(atom(atom_one)%num_nbrs) = nbr_entry
       if( full_nbr_list ) then
          atom(atom_two)%num_nbrs = atom(atom_two)%num_nbrs + 1
          if( atom(atom_two)%num_nbrs .gt. max_nbrs ) then
             write(run_out_unit,*) 'overran max_nbrs, atom_two'
             call abort()
          end if
          nbr_entry%nbr_index = atom_one
          atom(atom_two)%nbr_list(atom(atom_two)%num_nbrs) = nbr_entry
       end if

!! debugging/analysis
!!       if( dist_sq .gt. max_actual_dist_sq ) then
!!          max_actual_dist_sq = dist_sq
!!       end if
!! end debugging/analysis

    end if

    return
  end subroutine check_are_they_nbrs
    


!! debugging subroutine
!***************************************************************************
!                   count_binning
!***************************************************************************

  subroutine count_binning()
    implicit none

    integer :: how_many(high_bin_count)
    integer :: number_empty
    integer :: total_atoms
    integer :: total_bins
    integer :: ii, depth, curr_slot

    total_bins = num_slices(1) * num_slices(2) * num_slices(3)

    how_many = 0
    number_empty = 0

    do ii = 1, total_bins
       if( bin_head(ii) .eq. 0 ) then
          number_empty = number_empty + 1
       else
          depth = 0
          curr_slot = bin_head(ii)
          chain_loop: &
          do
             depth = depth + 1
             if( depth .gt. high_bin_count ) then
                write(run_out_unit,*) 'count_binning: internal error'
                call abort()
             end if
             if( atom(curr_slot)%next_atom .eq. 0 ) then
                exit chain_loop
             else
                curr_slot = atom(curr_slot)%next_atom
             end if
          end do chain_loop
          how_many(depth) = how_many(depth) + 1
       end if
    end do

    total_atoms = 0

    write(run_out_unit,*) 'number of bins: ',       total_bins
    write(run_out_unit,*) 'number of empty bins: ', number_empty
    write(run_out_unit,*) 'bin count, number of bins: '
    do ii = 1, high_bin_count
       write(run_out_unit,*) '   ', ii, how_many(ii)
       total_atoms = total_atoms + ii * how_many(ii)
    end do

    if( total_atoms .ne. num_atoms_built ) then
       write(run_out_unit,*) 'Problem: num_atoms_built, total_atoms: ', num_atoms_built, total_atoms
    else
       write(run_out_unit,*)
       write(run_out_unit,*) 'total atoms: ', total_atoms
    end if

    return
  end subroutine count_binning
!! debugging subroutine



!***************************************************************************
!                   bin_atoms
!***************************************************************************

  subroutine bin_atoms(max_nbr_dist_in)
    use build_gb_global_mod
    implicit none

    real*8, intent(in) :: max_nbr_dist_in

    integer :: atom_ii, jj, int_bin, curr_slot, depth
!!      integer :: cry_ii
    integer :: bin(3)
    real*8  :: pos(3), rel_pos(3)

    call setup_bins(max_nbr_dist_in)

    do atom_ii = 1, num_atoms_built
       if( .not. atom(atom_ii)%deleted_atom ) then
          pos = atom(atom_ii)%pos
          do jj = 1, 3
             rel_pos(jj) = pos(jj)                                         &
                - perlen(jj) * nint( (pos(jj)-box_center(jj))/perlen(jj) )
             bin(jj) = floor( ( rel_pos(jj) - perlb(jj)) / slice_len(jj) ) 
             if( bin(jj) .lt. -1 ) then
                write(run_out_unit,*) 'bin < -1'
                write(run_out_unit,*) 'jj: ', jj
                write(run_out_unit,*) 'pos: ', pos
                write(run_out_unit,*) 'perlen: ', perlen
                write(run_out_unit,*) 'box_center: ', box_center
                write(run_out_unit,*) 'slice_len: ', slice_len
                write(run_out_unit,*) 'num_slices: ', num_slices
                write(run_out_unit,*) 'rel_pos(jj): ', rel_pos(jj)
                write(run_out_unit,*) 'bin(jj): ', bin(jj)
                call abort()
             else if( bin(jj) .eq. -1 ) then
                bin(jj) = 0
             else if( bin(jj) .eq. num_slices(jj) ) then
                bin(jj) = num_slices(jj) - 1
             else if( bin(jj) .gt. num_slices(jj) ) then
                write(run_out_unit,*) 'bin > num_slices'
                write(run_out_unit,*) 'jj: ', jj
                write(run_out_unit,*) 'bin: ', bin
                write(run_out_unit,*) 'num_slices: ', num_slices
                write(run_out_unit,*) 'rel_pos: ', rel_pos
                write(run_out_unit,*) 'atom_ii: ', atom_ii
                call abort()
             end if
          end do
          atom(atom_ii)%rel_pos = rel_pos

          atom(atom_ii)%bin = bin

          call bin_to_int_bin(bin, int_bin)

!! debugging
!!       write(run_out_unit,"(i5, 3i4, i10)") atom_ii, bin, int_bin
!! end debugging


!! put atom_ii into the link list for bin "int_bin"
          if( bin_head(int_bin) .eq. 0 ) then
             bin_head(int_bin) = atom_ii
             if( 1 .gt. high_bin_count ) then
                high_bin_count = 1
             end if
          else
             curr_slot = bin_head(int_bin)
             depth = 1
             chain_loop: &
             do
                depth = depth + 1
                if( depth .gt. high_bin_count ) then
                   high_bin_count = depth
                end if
!! max_depth is intended to prevent infinite loops while debugging
                if( depth .gt. max_depth ) then
                   write(run_out_unit,*) 'overran max_depth'
                   call abort()
                end if
                if( atom(curr_slot)%next_atom .eq. 0 ) then
                   atom(curr_slot)%next_atom = atom_ii
                   exit chain_loop
                else
                   curr_slot = atom(curr_slot)%next_atom
                end if
             end do chain_loop
          end if
       end if      !! if atom_ii has not been deleted
    end do         !! atom_ii = 1, num_atoms_built


    return
  end subroutine bin_atoms

!***************************************************************************
!                   bin_to_int_bin
!***************************************************************************

  subroutine bin_to_int_bin(bin, int_bin)
    implicit none

    integer, intent(in)  :: bin(3)
    integer, intent(out) :: int_bin

    int_bin = 1                                                &
                 + bin(1)*num_slices(2)*num_slices(3)          &
                 + bin(2)*num_slices(3)                        &
                 + bin(3)
!! debugging
    if( int_bin .le. 0 ) then
       write(run_out_unit,*) 'bin_to_int_bin: internal error int_bin < 1'
       call abort()
    else if( int_bin .gt. num_bins ) then
       write(run_out_unit,*) 'bin_to_int_bin: internal error int_bin > num_bins'
       call abort()
    end if
!! end debugging

    return
  end subroutine bin_to_int_bin

!***************************************************************************
!                   setup_bins
!***************************************************************************

  subroutine setup_bins(max_nbr_dist_in)
    use build_gb_global_mod
    implicit none

    real*8, intent(in) :: max_nbr_dist_in

    integer :: ii, total_bins, atom_ii

    max_nbr_dist = max_nbr_dist_in
    max_nbr_dist_sq = max_nbr_dist**2

    total_bins = 1
    do ii = 1, 3
       num_slices(ii) = floor( perlen(ii) / max_nbr_dist )
       if( num_slices(ii) .le. 0 ) then
          write(run_out_unit,*) 'perlen < max_nbr_dist, bailing out.'
          call abort()
       end if
       slice_len(ii)  = perlen(ii) / num_slices(ii)
       if( slice_len(ii) .lt. max_nbr_dist_in ) then
          write(run_out_unit,*) 'internal error: slice_len'
          call abort()
       end if
       total_bins     = total_bins * num_slices(ii)
    end do

    if( total_bins .gt. max_bins ) then
       write(run_out_unit,*) 'overran max_bins'
       call abort()
    end if

    bin_head(1:total_bins) = 0
    do atom_ii = 1, num_atoms_built
       atom(atom_ii)%next_atom = 0
    end do
    num_bins = total_bins


    return
  end subroutine setup_bins


  subroutine get_sigma_number()
    use build_gb_global_mod
    use x_dir_to_crystal_mod
    implicit none

    real*8, parameter :: tiny_fudge = 1.0d-14
    real*8, parameter :: edge_fudge = 1.0d-12
    real*8, parameter :: small_fudge = 1.0d-10

!! debugging
!!$    logical :: atom_ii_matched
!!$    real*8  :: atom_ii_best_dist
!!$    integer :: atom_ii_best_match
!! end debugging

    integer :: atom_ii
!!    integer :: bin_ii, bin_jj, bin_kk
!!    integer :: int_bin
    integer :: search_ii, search_jj, search_kk, int_search
    integer :: atom_two
    integer :: bin(3), search(3)
    integer :: adj_perlen(3)
    real*8  :: test_len, adj_test_len
!!    real*8  :: other_adj_test_len

    integer :: loop_x_dir, other_x_dir
    integer :: loop_type, other_type
    integer :: num_loop_type_atoms
    integer :: num_other_type_atoms
    integer :: num_matches

    real*8  :: pos_xx, adj_x
    real*8  :: pos(3), rel_pos(3)
    integer :: jj
    real*8  :: our_perlen(3), our_center(3)

    real*8  :: largest_accepted = -1.0d0
    real*8  :: smallest_rejected = 9999.0d0

    if( left_len(1) .lt. right_len(1) ) then
       call x_dir_to_crystal_ii(-1, loop_type)
       test_len = left_len(1)
       loop_x_dir = -1
       other_x_dir = 1
       adj_x = test_len
    else
       call x_dir_to_crystal_ii(1, loop_type)
       test_len = right_len(1)
       loop_x_dir = 1
       other_x_dir = -1
       adj_x = -1 * test_len
    end if
    adj_test_len = test_len - edge_fudge
!!    other_adj_test_len = test_len + edge_fudge
    if( loop_type .eq. 1 ) then
       other_type = 2
    else if( loop_type .eq. 2 ) then
       other_type = 1
    else
       write(run_out_unit,*) 'get_sigma_number: case_trap: loop_type'
       call abort()
    end if

!! debugging
!!$    write(run_out_unit,*) 'loop_type: ', loop_type
!!$    write(run_out_unit,*) 'other_type: ', other_type
!!$    write(run_out_unit,*) 'test_len: ', test_len
!!$    write(run_out_unit,*) 'loop_x_dir: ', loop_x_dir
!!$    write(run_out_unit,*) 'other_x_dir: ', other_x_dir
!!$    write(run_out_unit,*) 'adj_x: ', adj_x
!!$    write(run_out_unit,*) 'adj_test_len: ', adj_test_len
!! end debugging


    our_perlen(1) = test_len
    our_perlen(2) = perlen(2)
    our_perlen(3) = perlen(3)
    our_center(1) = 0.5d0 * adj_x
    our_center(2) = box_center(2)
    our_center(3) = box_center(3)

    num_matches = 0
    num_loop_type_atoms = 0
    num_other_type_atoms = 0

    atom_loop:  &
    do atom_ii = 1, num_atoms_built
       if( atom(atom_ii)%crystal .eq. other_type ) then
          pos_xx = other_x_dir * atom(atom_ii)%pos(1)
          if( pos_xx .le. adj_test_len ) then
!!            if(       pos_xx .gt. edge_fudge                  &
!!                .and. pos_xx .le. other_adj_test_len ) then
             num_other_type_atoms = num_other_type_atoms + 1
          end if
          cycle atom_loop
       end if
!! Now we know it is a loop type atom.  But for the two boundary setup,
!! we may want to exclude it.
       pos_xx = loop_x_dir * atom(atom_ii)%pos(1)
       if( .false. ) then                         !! this is old (change
                                                  !! atoms_two treatment if
                                                  !! you want this to work.)
          if( pos_xx .gt. adj_test_len ) then
             cycle atom_loop
          end if
       else                                     !! this is new
          if( pos_xx .lt. edge_fudge ) then
             cycle atom_loop
          end if
       end if

!! Now we want to try to match it.
       num_loop_type_atoms = num_loop_type_atoms + 1

!! debugging
!!$       atom_ii_matched = .false.
!!$       atom_ii_best_dist = 9999.0d0
!!$       atom_ii_best_match = -1
!! end debugging

       pos = atom(atom_ii)%pos
       pos(1) = pos(1) + adj_x


!! new version follows

       rel_pos(1) = pos(1)
       do jj = 2, 3
          rel_pos(jj) = pos(jj)                                         &
                - our_perlen(jj)                                        &
                  * nint( (pos(jj)-our_center(jj))/our_perlen(jj) )
       end do
       do jj = 1, 3

!! end of new version, orig follows commented out

!!$ was      do jj = 1, 3
!!$ was         rel_pos(jj) = pos(jj)                                         &
!!$ was               - our_perlen(jj)                                        &
!!$ was                 * nint( (pos(jj)-our_center(jj))/our_perlen(jj) )

          bin(jj) = floor( ( rel_pos(jj) - perlb(jj)) / slice_len(jj) ) 
          if( bin(jj) .lt. -1 ) then
             write(run_out_unit,*) 'get_sigma_number: bin < -1'
             write(run_out_unit,*) 'jj: ', jj
             write(run_out_unit,*) 'pos: ', pos
             write(run_out_unit,*) 'our_perlen: ', our_perlen
             write(run_out_unit,*) 'our_center: ', our_center
             write(run_out_unit,*) 'slice_len: ', slice_len
             write(run_out_unit,*) 'num_slices: ', num_slices
             write(run_out_unit,*) 'rel_pos(jj): ', rel_pos(jj)
             write(run_out_unit,*) 'bin(jj): ', bin(jj)
             call abort()
          else if( bin(jj) .eq. -1 ) then
             bin(jj) = 0
          else if( bin(jj) .eq. num_slices(jj) ) then
             bin(jj) = num_slices(jj) - 1
          else if( bin(jj) .gt. num_slices(jj) ) then
             write(run_out_unit,*) 'get_sigma_number: bin > num_slices'
             write(run_out_unit,*) 'jj: ', jj
             write(run_out_unit,*) 'bin: ', bin
             write(run_out_unit,*) 'num_slices: ', num_slices
             write(run_out_unit,*) 'rel_pos: ', rel_pos
             write(run_out_unit,*) 'atom_ii: ', atom_ii
             call abort()
          end if
       end do

       do search_ii = -1, 1
          search(1) = bin(1) + search_ii
          if( search(1) .ge. num_slices(1) ) then
             search(1) = 0
             adj_perlen(1) = 1
          else if( search(1) .eq. -1 ) then
             search(1) = num_slices(1) - 1
             adj_perlen(1) = -1
          else
             adj_perlen(1) = 0
          end if
          do search_jj = -1, 1
             search(2) = bin(2) + search_jj
             if( search(2) .ge. num_slices(2) ) then
                search(2) = 0
                adj_perlen(2) = 1
             else if( search(2) .eq. -1 ) then
                search(2) = num_slices(2) - 1
                adj_perlen(2) = -1
             else
                adj_perlen(2) = 0
             end if
             search_kk_loop: &
             do search_kk = -1, 1
               search(3) = bin(3) + search_kk
               if( search(3) .ge. num_slices(3) ) then
                  search(3) = 0
                  adj_perlen(3) = 1
               else if( search(3) .eq. -1 ) then
                  search(3) = num_slices(3) - 1
                  adj_perlen(3) = -1
               else
                  adj_perlen(3) = 0
               end if

               call bin_to_int_bin(search, int_search)

               if( bin_head(int_search) .eq. 0 ) then
                  cycle search_kk_loop    !! empty bin
               end if

               atom_two = bin_head(int_search)
               chain_two_loop: &
               do
                  if( atom(atom_two)%crystal .eq. other_type ) then
                     pos_xx = other_x_dir * atom(atom_two)%pos(1)
!!                       if(       pos_xx .gt. edge_fudge                  &
!!                           .and. pos_xx .le. other_adj_test_len ) then
                     if( pos_xx .le. adj_test_len ) then
                        call check_do_they_match()
                     end if
                  end if

                  if( atom(atom_two)%next_atom .eq. 0 ) then
                     exit chain_two_loop
                  else
                     atom_two = atom(atom_two)%next_atom
                  end if
               end do chain_two_loop
            end do search_kk_loop
         end do                     !! search_jj
      end do                        !! search_ii

!! debugging
!!$      if( .not. atom_ii_matched ) then
!!$         write(run_out_unit,*) 'atom_ii_best_dist: ', atom_ii_best_dist
!!$         if( atom_ii_best_dist .eq. 9999.0d0 ) then
!!$            write(run_out_unit,"(6f12.4)") atom(atom_ii)%pos, rel_pos
!!$         else
!!$            write(run_out_unit,"(6f12.4)") atom(atom_ii)%pos, atom(atom_ii_best_match)%pos
!!$         end if
!!$      end if
!! end debugging

   end do atom_loop

   write(run_out_unit,*) 'num_loop_type_atoms:  ', num_loop_type_atoms
   write(run_out_unit,*) 'num_other_type_atoms: ', num_other_type_atoms
   write(run_out_unit,*) 'num_matches:          ', num_matches

   if( num_other_type_atoms .ne. num_loop_type_atoms ) then
      write(run_out_unit,*) 'get_sigma_number: internal error: num atoms loop/other'
      call abort()
   end if

   if( modulo(num_loop_type_atoms, num_matches) .ne. 0 ) then
      write(run_out_unit,*) 'Warning: get_sigma_number: bad juju: sigma not an integer'
      sigma_number = -1
   else
      sigma_number = num_loop_type_atoms/num_matches
      write(run_out_unit,*) 'Sigma: ', sigma_number
   end if

   write(run_out_unit,*) 'smallest_rejected: ', smallest_rejected
   write(run_out_unit,*) 'largest_accepted:  ', largest_accepted

   return

 contains
   subroutine check_do_they_match()
     implicit none

     integer :: dim_ii
     logical :: match
     real*8  :: rel_pos_two(3)
     real*8  :: diff

!! debugging
!!$     real*8  :: temp_dist
!! end debugging

     rel_pos_two = atom(atom_two)%rel_pos
     
     if( atom(atom_ii)%deleted_atom ) then
        write(run_out_unit,*) 'check_do_they_match: atom_ii is a deleted atom'
        call abort()
     end if
     if( atom(atom_two)%deleted_atom ) then
        write(run_out_unit,*) 'check_do_they_match: atom_two is a deleted atom'
        call abort()
     end if

     match = .true.
     dim_loop:  &
     do dim_ii = 1, 3
        diff = rel_pos_two(dim_ii) - rel_pos(dim_ii)
        diff = diff + adj_perlen(dim_ii) * perlen(dim_ii)
        if( abs(diff) .gt. small_fudge ) then
           match = .false.
           if( abs(diff) .lt. smallest_rejected ) then
              smallest_rejected = abs(diff)
           end if
           exit dim_loop
        else
           if( abs(diff) .gt. largest_accepted ) then
              largest_accepted = abs(diff)
           end if
        end if
     end do dim_loop

     if( match ) then

!! debugging
!!$        write(run_out_unit,"(3f12.4,'  ',3f12.4)") rel_pos, rel_pos_two
!! end debugging

!! debugging
!!$        if( atom_ii_matched ) then
!!$           write(run_out_unit,*) 'check_do_they_match: already matched'
!!$           call abort()
!!$        end if
!!$        atom_ii_matched = .true.
!! end debugging

        num_matches = num_matches + 1

!! debugging
!!$     else
!!$        temp_dist = 0.0d0
!!$        do dim_ii = 1, 3
!!$           diff = rel_pos_two(dim_ii) - rel_pos(dim_ii)
!!$           diff = diff + adj_perlen(dim_ii) * perlen(dim_ii)
!!$           temp_dist = temp_dist + diff**2
!!$        end do
!!$        temp_dist = sqrt(temp_dist)
!!$        if( temp_dist .lt. atom_ii_best_dist ) then
!!$           atom_ii_best_dist = temp_dist
!!$           atom_ii_best_match = atom_two
!!$        end if
!! end debugging

     end if

     return
   end subroutine check_do_they_match

 end subroutine get_sigma_number

end module gb_neighbor_mod
