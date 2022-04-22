!! randomize_type_dsc.f90

module randomize_type_dsc_mod
  use build_gb_global_mod
  use build_gb_atom_data_mod
  use lammps_data_file_mod
  use scatter_randomly_g95_mod
  implicit none

  integer, parameter :: type_diff = 5

contains

  subroutine randomize_type_dsc(atom_table)
    implicit none

    type( lammps_df_atom ), intent(inout)  :: atom_table( actual_num_atoms )

    integer :: total_sites
    integer :: solute_sites
    integer :: matrix_type
    integer :: solute_type
    integer :: type_two_sites
    logical :: flipped_types

    integer, allocatable :: solute_list(:)

    integer :: tag_ii, solute_ii
    integer :: this_solute, this_species
    integer :: num_solute, num_matrix

    total_sites = actual_num_atoms
    type_two_sites = nint( concentration * total_sites )
    if( dble(type_two_sites) .gt. 0.5d0 * total_sites ) then
       solute_sites = total_sites - type_two_sites
       flipped_types = .true.
       matrix_type = type_diff
       solute_type = 0
    else
       solute_sites = type_two_sites
       flipped_types = .false.
       matrix_type = 0
       solute_type = type_diff
    end if

    if( dble(solute_sites) .gt. 0.5d0 * total_sites ) then
       write(0,*) 'randomize_type: internal error, solute_sites'
       call abort()
    end if

    if( solute_sites .le. 0 ) then 
       write(0,*) 'concentration too small, no solutes'
       write(0,*) 'solute_sites: ', solute_sites
       call abort()
    end if

!! debugging
    write(0,*) 'randomize_type'
    write(0,*) 'concentration: ', concentration
    write(0,*) 'total_sites: ', total_sites
    write(0,*) 'solute_sites: ', solute_sites
    write(0,*) 'flipped_types: ', flipped_types
    write(0,*) 'matrix_type: ', matrix_type
    write(0,*) 'solute_type: ', solute_type
!! end debugging

    allocate( solute_list(1:solute_sites) )
    call scatter_randomly_g95(solute_sites, total_sites, 0, solute_list)
    write(0,*) 'Warning: seed not yet implemented'
    do tag_ii = 1, total_sites
       atom_table(tag_ii)%species =                           &
            atom_table(tag_ii)%species + matrix_type 
    end do
    do solute_ii = 1, solute_sites
       this_solute = solute_list(solute_ii)
!! debugging
!!       write(0,*) this_solute
!! end debugging
       atom_table(this_solute)%species =                       &
            atom_table(this_solute)%species + solute_type - matrix_type
    end do
    deallocate( solute_list )

!! sanity check
    num_matrix = 0
    num_solute = 0
    do tag_ii = 1, total_sites
       this_species = atom_table(tag_ii)%species
       if( flipped_types ) then
          if( this_species .ge. type_diff ) then
             num_matrix = num_matrix + 1
          else
             num_solute = num_solute + 1
          end if
       else
          if( this_species .ge. type_diff ) then
             num_solute = num_solute + 1
          else
             num_matrix = num_matrix + 1
          end if
       end if
    end do
    if( num_solute .ne. solute_sites ) then
       write(0,*) 'randomize_type: internal error, num_solute'
       call abort()
    end if
    if( num_solute + num_matrix .ne. total_sites ) then
       write(0,*) 'randomize_type: internal error, num_matrix'
       call abort()
    end if
!! end sanity check

    return
  end subroutine randomize_type_dsc


end module randomize_type_dsc_mod
