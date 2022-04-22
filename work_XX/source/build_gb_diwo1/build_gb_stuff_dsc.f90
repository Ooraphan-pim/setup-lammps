!! build_gb_stuff_dsc.f90

module build_gb_stuff_dsc_mod
  use system_mod
  implicit none

  real*8, parameter, private :: tiny_fudge = 1.0d-14

!!   interface
!!      integer*4 function system(command_line)
!!        character(*), intent(in) :: command_line
!!      end function system
!!   end interface
 
 

contains

!!**********************************************************
!!               setup_deletion_control_file
!!**********************************************************

  subroutine setup_deletion_control_file(init, num_atoms_expected)
    use build_gb_control_mod
    use delete_file_mod
    use int_to_string2_mod
    implicit none

    type(build_gb_control_type), intent(in) :: init
    integer, intent(in) :: num_atoms_expected
    
    character(256) :: line
    character(256) :: temp_name
    integer :: read_status
    integer :: tag_in, type_in
    real*8  :: energy_in
    real*8  :: ss_diag_in(3)
    real*8  :: pressure
    real*8  :: pos_in(3)
    logical :: first_time
    real*8  :: worst_energy = 9999.0d0
    real*8  :: worst_pressure = -9999.0d0
    real*8  :: worst_pos(3) = 9999.0d0
    integer :: worst_tag = -1
    real*8  :: actual_worst_energy = 9999.0d0
    real*8  :: actual_worst_pressure = -9999.0d0
    integer :: num_atoms_in
    integer :: skip_ii
    character(256) :: string_one, string_two
    integer :: atoms_read
    integer :: result, our_type

    temp_name = trim(init%name)//".dump"
    open(99, file=temp_name, status='old', action='read')
    do skip_ii = 1, 3
       read(99,*)
    end do
    read(99,*) num_atoms_in
    if( num_atoms_in .ne. num_atoms_expected ) then
       write(0,*) "delete_by_energy  num_atoms mismatch"
       write(0,*) "num_atoms_expected  ", num_atoms_expected
       write(0,*) "num_atoms_in        ", num_atoms_in
       call abort()
    end if
    do skip_ii = 1, 4
       read(99,*)
    end do
    read(99,*) string_one, string_two
    if( string_two .ne. "ATOMS" ) then
       write(0,*) "delete_by_energy  problem finding ATOMS header"
       write(0,*) "string_one  ", string_one
       write(0,*) "string_two  ", string_two
       call abort()
    end if

    first_time = .true.
    atoms_read = 0
    read_loop:  &
    do
       if( init%show_pressure ) then
          read(99, *, iostat=read_status) tag_in, type_in, pos_in,      &
                                          energy_in, ss_diag_in
       else
          read(99, *, iostat=read_status) tag_in, type_in, pos_in, energy_in
       end if
       if( read_status .gt. 0 ) then
          write(0,*) "setup_deletion_control_file  read error on dump file"
          call abort()
       else if( read_status .lt. 0 ) then
          exit read_loop
       end if
       atoms_read = atoms_read + 1
       if( init%show_pressure ) then
          pressure = -sum(ss_diag_in) / 3.0d0
       end if
       our_type = type_in
       if( our_type .gt. 5 ) then
          our_type = our_type - 5
       end if
       if( our_type .gt. 3 ) then
          cycle read_loop           !! edge atom, not grain boundary
       end if
       if( first_time ) then
          first_time = .false.
          worst_tag = tag_in
          worst_energy = energy_in
          actual_worst_energy = energy_in
          worst_pos = pos_in
          if( init%show_pressure ) then
             worst_pressure = pressure
             actual_worst_pressure = pressure
          end if
       else
          if( init%delete_by_energy ) then
             if( energy_in .gt. worst_energy ) then
                worst_tag = tag_in
                worst_energy = energy_in
                worst_pos = pos_in
                if( init%show_pressure ) then
                   worst_pressure = pressure
                end if
             end if
             if( init%show_pressure ) then
                if( pressure .gt. actual_worst_pressure ) then
                   actual_worst_pressure = pressure
                end if
             end if
          else if( init%delete_by_pressure ) then
             if( pressure .gt. worst_pressure ) then
                worst_tag = tag_in
                worst_pressure = pressure
                worst_energy = energy_in
                worst_pos = pos_in
             end if
             if( energy_in .gt. actual_worst_energy ) then
                actual_worst_energy = energy_in
             end if
          else
             write(0,*) "setup_deleteion_control_file  internal error"
             call abort()
          end if
       end if
    end do read_loop
    close(99)

    if( atoms_read .ne. num_atoms_in ) then
       write(0,*) "setup_deletion_control_file  num_atoms mismatch(2)"
       write(0,*) "num_atoms_in  ", num_atoms_in
       write(0,*) "atoms_read    ", atoms_read
       call abort()
    end if

!! debugging
    if( init%delete_by_energy ) then
       if( init%show_pressure ) then
          write(0,"('worst: tag, pos, energy, pressure, actual_worst_pressure  ',"    &
                    //"i5, 3f14.4, f9.4, 2g12.4)")  worst_tag,                          &
                    worst_pos, worst_energy, worst_pressure, actual_worst_pressure
       else
          write(0,"('worst: tag, pos, energy  ',"        &
                    //"i5, 3f14.4, f9.4)")  worst_tag,    &
                    worst_pos, worst_energy
       end if
    else
       write(0,"('worst: tag, pos, pressure, energy, actual_worst_energy  ',"    &
                 //"i5, 3f14.4, g12.4, 2f9.4)")  worst_tag,                           &
                 worst_pos, worst_pressure, worst_energy, actual_worst_energy
    end if
!! end debugging

    line = "cat "//trim(init%deletion_control_file)                 &
           //" | sed ""s/Q/"//trim(init%name)                       &
           //"/;s/Z/"//trim(int_to_string2(worst_tag))              &
           //"/"" > "//trim(init%name)//".min.in"
    result = system(trim(line)//char(0))
    if( result .ne. 0 ) then
       write(0,*) "setup_deletion_control_file:  system->sed  failed"
       call abort()
    end if

    temp_name = trim(init%name)//".dump"
!!    call delete_file(trim(temp_name), len_trim(temp_name))    !! not working on fenrir
    call delete_file(temp_name, len_trim(temp_name))

!! debugging
!!    stop
!! end debugging

    return
  end subroutine setup_deletion_control_file



!!**********************************************************
!!               xi_chi_files
!!**********************************************************

  subroutine xi_chi_files()
    use build_gb_global_mod
    implicit none

    character(256) xi_or_chi_file

    xi_or_chi_file = trim(name)//".half_xi.vectors"
    call do_xi_or_chi_file(xi_or_chi_file, left_trans)
    xi_or_chi_file = trim(name)//".half_chi.vectors"
    call do_xi_or_chi_file(xi_or_chi_file, right_trans)

    return
  end subroutine xi_chi_files


!!**********************************************************
!!               do_xi_or_chi_file
!!**********************************************************

  subroutine do_xi_or_chi_file( filename, trans_matrix )
    use build_gb_global_mod
    implicit none

    character(256), intent(in) :: filename
    real*8, intent(in)         :: trans_matrix(3,3)

    integer :: ii
    real*8  :: trans_vec(3)

    open(99, file=filename, status='new')
    do ii = 1, mat%crystal%num_half_nn
       trans_vec = matmul( trans_matrix, mat%crystal%half_nn(1:3, ii) )
       trans_vec = trans_vec * xi_chi_alat
       write(99,"(3f21.15)") trans_vec
    end do
    close(99)

    return
  end subroutine do_xi_or_chi_file


!!**********************************************************
!!               compute_max_index
!!**********************************************************

  subroutine compute_max_index( length_vector, max_index )
    use material_mod
    use build_gb_global_mod
    implicit none

    real*8, intent(in)   :: length_vector(3)
    integer, intent(out) :: max_index

    real*8  :: diagonal = 0.0d0
    integer :: ii

    do ii = 1, 3
       diagonal = diagonal + length_vector(ii)*length_vector(ii)
    end do

    diagonal = sqrt( diagonal )
    max_index = ceiling( diagonal / mat%alat ) + 1

    return
  end subroutine compute_max_index

  subroutine build_atoms_dsc( max_index_in,       &
                              trans_matrix,       &
                              length_vector,      &
                              edge_size,          &
                              x_dir,              &
                              dsc_vec_in,         &
                              xx_adj_in )
    use material_mod
    use build_gb_global_mod
    use build_gb_atom_data_mod
    use x_dir_to_crystal_mod
    implicit none

    integer, intent(in)  :: max_index_in
    real*8, intent(in)   :: trans_matrix(3, 3)
    real*8, intent(in)   :: length_vector(3)
    real*8, intent(in)   :: edge_size
    integer, intent(in)  :: x_dir                 !!  1 or -1
    real*8, intent(in)   :: dsc_vec_in(3)
    real*8, intent(in)   :: xx_adj_in

    logical :: keep_flag
    integer :: max_index
    integer :: ii_xx, ii_yy, ii_zz, ii_bb
    integer :: crystal_ii
    real*8  :: cell_base(3)
    real*8  :: fcc_pos(3)              !! position before rotation
    real*8  :: trans_pos(3)            !! position after rotation

    real*8  :: x_low_cutoff, x_high_cutoff
    real*8  :: dsc_vec(3)
    real*8  :: dsc_vec_len
    real*8  :: xx_adj

    dsc_vec_len = sqrt( dot_product(dsc_vec_in, dsc_vec_in) )
    max_index = max_index_in + ceiling( dsc_vec_len / 2.0d0 ) + 1
    max_index = max_index + ceiling( abs(xx_adj_in) / 2.0d0 ) + 1
    dsc_vec = dsc_vec_in * mat%alat / 2.0d0
    xx_adj = xx_adj_in * mat%alat / 2.0d0

!! debugging
!!    write(0,*) "WARNING  using twice dsc_vec for debugging"
!!    dsc_vec = dsc_vec_in * mat%alat
!! end debugging

!! debugging
!!    write(6,*) "dsc_vec  ", dsc_vec
!! end debugging

    x_low_cutoff = -tiny_fudge * length_vector(1)
    if( two_boundaries ) then
       x_high_cutoff = (1 + tiny_fudge) * length_vector(1)
    else
       x_high_cutoff = (1 - tiny_fudge) * length_vector(1)
    end if

    do ii_xx = -max_index, max_index
       cell_base(1) = ii_xx * mat%alat
       do ii_yy = -max_index, max_index
          cell_base(2) = ii_yy * mat%alat
          do ii_zz = -max_index, max_index
             cell_base(3) = ii_zz * mat%alat
             do ii_bb = 1, mat%crystal%num_basis
                fcc_pos = cell_base + mat%alat * mat%crystal%basis(1:3, ii_bb)
                fcc_pos = fcc_pos + dsc_vec
                call do_atom()
             end do
          end do
       end do
    end do

    return

  contains
    subroutine do_atom()
      implicit none

      trans_pos = matmul(trans_matrix, fcc_pos)
      if(       ( trans_pos(3) .gt. -tiny_fudge * length_vector(3) )         &
          .and. ( trans_pos(3) .lt. (1 - tiny_fudge) * length_vector(3) )    &
          .and. ( trans_pos(2) .gt. -tiny_fudge * length_vector(2) )         &
          .and. ( trans_pos(2) .lt. (1 - tiny_fudge) * length_vector(2) ) )     then
         !! not eliminated yet
      else
         return
      end if

      call x_dir_to_crystal_ii( x_dir, crystal_ii )

      if( island_flag ) then
         if( ( trans_pos(1) .gt. -tiny_fudge * length_vector(1) )        &
              .and. ( trans_pos(1) .lt. (1 - tiny_fudge) * length_vector(1) ) ) then
            !! inside simulation cell
         else
            return
         end if
         call cut_out_island( trans_pos, crystal_ii, keep_flag )
         if( .not. keep_flag ) then
            return
         end if
      else
         trans_pos(1) = trans_pos(1) - xx_adj

         if(       ( x_dir * trans_pos(1) .gt. x_low_cutoff ) &
             .and. ( x_dir * trans_pos(1) .lt. x_high_cutoff ) ) then
            !! not eliminated
         else
            return
         end if

         if( x_dir * trans_pos(1) .ge. x_high_cutoff - edge_size ) then
            crystal_ii = crystal_ii + 3    !! It is an edge atom, to be "fixed"
         end if
      end if

      num_atoms_built = num_atoms_built + 1
      if( num_atoms_built .gt. max_atoms ) then
         write(0,*) 'Overran max_atoms.'
         call exit(1)
      end if
      atom(num_atoms_built)%crystal = crystal_ii
      atom(num_atoms_built)%pos = trans_pos
      atom(num_atoms_built)%deleted_atom = .false.

      return
    end subroutine do_atom


  end subroutine build_atoms_dsc

  
  subroutine cut_out_island( position, crystal, keep_flag )
    use build_gb_param_mod
    use material_mod
    use build_gb_global_mod
    implicit none

    real*8, intent(in)   :: position(3)
    integer, intent(in)  :: crystal
    logical, intent(out) :: keep_flag

    real*8  :: radius_sq, radius_sq_limit
    integer :: dim_ii

    if( island_cut_type .eq. island_circle ) then
       radius_sq_limit = ((mat%alat/2.0d0) * island_radius)**2
       radius_sq = 0.0d0
       do dim_ii = 1, 2
          radius_sq = radius_sq + (position(dim_ii) - box_center(dim_ii))**2
       end do

       if( crystal .eq. 1 ) then
          if( radius_sq .gt. radius_sq_limit ) then
             keep_flag = .true.
          else
             keep_flag = .false.
          end if
       else if( crystal .eq. 2 ) then
          if( radius_sq .lt. radius_sq_limit ) then
             keep_flag = .true.
          else
             keep_flag = .false.
          end if
       else
          write(0,*) "cut_out_island  internal error  crystal  ", crystal
          call exit(1)
       end if
    else
       write(0,*) "cut_out_island  internal error  island_cut_type  ", island_cut_type
    end if

    return
  end subroutine cut_out_island


end module build_gb_stuff_dsc_mod
