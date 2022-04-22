!! build_gb_atom_code.f90

module build_gb_atom_code_dsc_mod
  use build_gb_atom_data_mod
  use build_gb_run_out_mod
  implicit none


contains

!! !!**********************************************************
!! !!               print_atoms
!! !!**********************************************************
!! 
!!   subroutine print_atoms()
!!     implicit none
!! 
!!     integer atom_ii, tag_ii
!! 
!! !! debugging
!!     write(0,*) 'actual_num_atoms: ', actual_num_atoms
!! !! end debugging
!! 
!!     tag_ii = -1
!!     do atom_ii = 1, num_atoms_built
!!        if( .not. atom(atom_ii)%deleted_atom ) then
!!           tag_ii = tag_ii + 1
!!           write(0,"(i10,i3,3f12.6)")                             & 
!!                   tag_ii,                                        &
!!                   atom(atom_ii)%crystal - 1,                     &
!!                   atom(atom_ii)%pos
!!        end if
!!     end do
!!     if( tag_ii .ne. actual_num_atoms - 1 ) then
!!        write(0,*) 'internal error: tag_ii'
!!        write(0,*) 'tag_ii+1: ', tag_ii + 1
!!        write(0,*) 'actual_num_atoms: ', actual_num_atoms
!!        call abort()
!!     end if
!! 
!!     return
!!   end subroutine print_atoms




!!**********************************************************
!!               write_atoms_dsc
!!**********************************************************

  subroutine write_atoms_dsc(alloy_flag, use_orig_data_flag)
    use build_gb_global_mod
    use randomize_type_dsc_mod
    use lammps_data_file_mod
    implicit none

    logical, intent(in) :: alloy_flag
    logical, intent(in) :: use_orig_data_flag  !! .true. for delete_by_energy

    type( lammps_df_header ) :: data_header
    type( lammps_df_atom )   :: atom_table( actual_num_atoms )

    integer atom_ii, tag_ii
    character(256) data_file
    logical :: our_deleted_flag

    data_file = trim(name)//".data"

    data_header%header = "build_gb"
    data_header%num_atoms = actual_num_atoms
    if( alloy_flag ) then
       if( double_edge ) then
          data_header%num_species = 14
       else
          data_header%num_species = 10
       end if
    else
       if( double_edge ) then
          data_header%num_species = 7
       else
          data_header%num_species = 5
       end if
    end if
    data_header%perlb = perlb
    data_header%perub = perub

!! debugging
    write(run_out_unit,*) 'actual_num_atoms: ', actual_num_atoms
!! end debugging

    tag_ii = 0
    do atom_ii = 1, num_atoms_built
       if( use_orig_data_flag ) then
          our_deleted_flag = atom(atom_ii)%deleted_first_pass
       else
          our_deleted_flag = atom(atom_ii)%deleted_atom
       end if
       if( .not. our_deleted_flag ) then
          tag_ii = tag_ii + 1
          atom_table(tag_ii)%tag = tag_ii - 1   !! C numbering from zero
!! write_lammps data file adds back the "1" to both tag and species! :-) 
          atom_table(tag_ii)%species = atom(atom_ii)%crystal - 1
          if( use_orig_data_flag ) then
             atom_table(tag_ii)%pos = atom(atom_ii)%orig_pos
          else
             atom_table(tag_ii)%pos = atom(atom_ii)%pos
          end if
!! debugging
!!$          write(0,"(i10,i3,3f12.6)")                             & 
!!$                  tag_ii - 1,                                    &
!!$                  atom(atom_ii)%crystal - 1,                     &
!!$                  atom(atom_ii)%pos
!! end debugging

       end if
    end do
    if( tag_ii .ne. actual_num_atoms ) then
       write(run_out_unit,*) 'internal error: tag_ii'
       call abort()
    end if

    if( alloy_flag ) then
       call randomize_type_dsc( atom_table )
    end if

    call write_lammps_data_file( data_header, atom_table, data_file )

    return
  end subroutine write_atoms_dsc


  subroutine mark_orig_pos()
    implicit none

    integer atom_ii

    do atom_ii = 1, num_atoms_built
       atom(atom_ii)%deleted_first_pass = atom(atom_ii)%deleted_atom
       atom(atom_ii)%orig_pos = atom(atom_ii)%pos
    end do

    orig_atoms_deleted = num_atoms_deleted

    return
  end subroutine mark_orig_pos


end module build_gb_atom_code_dsc_mod
