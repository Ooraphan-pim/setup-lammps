module lammps_data_file_mod
  implicit none

  type :: lammps_df_atom
     integer :: tag             !! tag, species, mol numbered from zero
     integer :: species
     real*8  :: pos(3)
     integer :: mol             !! optional (molecule)
     real*8  :: charge          !! optional
     real*8  :: vel(3)          !! optional
  end type lammps_df_atom

  type :: lammps_df_header
     character(80) :: header
     integer :: num_atoms
     integer :: num_species
     real*8  :: perlb(3)
     real*8  :: perub(3)
     integer :: num_bonds
     integer :: num_angles
     integer :: num_bond_types
     integer :: num_angle_types
  end type lammps_df_header

  type :: lammps_df_bond
     integer :: bond_tag          !! all counted from 1
     integer :: bond_type
     integer :: atom_one
     integer :: atom_two
  end type lammps_df_bond

  type :: lammps_df_angle
     integer :: angle_tag
     integer :: angle_type
     integer :: atom_one
     integer :: atom_two
     integer :: atom_three
  end type lammps_df_angle


!!$  type lammps_config
!!$     type( lammps_df_header ) data_header
!!$     type( lammps_df_atom )   atom( data_header%num_atoms)
!!$  end type lammps_config

contains

!!$  subroutine write_lammps_data_file( config, filename )
!!$    implicit none
!!$
!!$    type( lammps_config), intent(in)        :: config
!!$    character(256), intent(in)              :: filename
!!$
!!$    open(99, file=filename, status='new')
!!$    !! write the file here
!!$    close(99)
!!$
!!$    return
!!$  end subroutine write_lammps_data_file


!!    filename = "-" means stdin
  subroutine write_lammps_data_file( header,                  &
                                     atom_table,              &
                                     filename,                &
                                     mol_flag,                &
                                     charge_flag,             &
                                     bond_flag,               &
                                     bond_table,              &
                                     angle_flag,              &
                                     angle_table,             &
                                     mass_flag,               &
                                     mass,                    &
                                     vel_flag )
    implicit none

    type( lammps_df_header), intent(in)     :: header
    type( lammps_df_atom), intent(in)       :: atom_table(header%num_atoms)
    character(256), intent(in)              :: filename
    logical, optional, intent(in)           :: mol_flag
    logical, optional, intent(in)           :: charge_flag
    logical, optional, intent(in)           :: bond_flag
    type( lammps_df_bond ), optional, intent(in) :: bond_table(header%num_bonds)
    logical, optional, intent(in)           :: angle_flag
    type( lammps_df_angle ), optional, intent(in) :: angle_table(header%num_angles)
    logical, optional, intent(in)           :: mass_flag
    real*8, optional, intent(in)            :: mass(header%num_species)
    logical, optional, intent(in)           :: vel_flag

    integer :: header_length
    integer :: atom_ii, type_ii, bond_ii, angle_ii
    logical :: our_mol_flag
    logical :: our_charge_flag
    logical :: our_vel_flag
    logical :: our_bond_flag
    logical :: our_angle_flag
    logical :: our_mass_flag

    integer :: out_unit

!! debugging
!!    write(0,*) "write_lammps_data_file  filename:  ", trim(filename)
!! end debugging

    our_mol_flag = .false.
    if( present(mol_flag) ) then
       if( mol_flag ) then
          our_mol_flag = .true.
       end if
    end if

    our_charge_flag = .false.
    if( present(charge_flag) ) then
       if( charge_flag ) then
          our_charge_flag = .true.
       end if
    end if
    
    our_bond_flag = .false.
    if( present(bond_flag) ) then
       if( bond_flag ) then
          our_bond_flag = .true.
          if( .not. present( bond_table ) ) then
             write(0,*) "write_lammps_data_file  bond_flag, but not bond_table"
             call exit(1)
          end if
       end if
    end if
    if( .not. our_bond_flag ) then
       if( present(bond_table) ) then
          if( header%num_bonds .ne. 0 ) then
             write(0,*) "write_lammps_data_file  bond_table, but not bond_flag=.true."
             call exit(1)
          end if
       end if
    end if

    our_angle_flag = .false.
    if( present(angle_flag) ) then
       if( angle_flag ) then
          our_angle_flag = .true.
          if( .not. present(angle_table) ) then
             write(0,*) "write_lammps_data_file  angle_flag, but not angle_table"
             call exit(1)
          end if
       end if
    end if
    if( .not. our_angle_flag ) then
       if( present(angle_table) ) then
          if( header%num_angles .ne. 0 ) then
             write(0,*) "write_lammps_data_file  angle_table, but not angle_flag=.true."
             call exit(1)
          end if
       end if
    end if

    our_mass_flag = .false.
    if( present(mass_flag) ) then
       if( mass_flag ) then
          our_mass_flag = .true.
          if( .not. present(mass) ) then
             write(0,*) "write_lammps_data_file  mass_flag, but not mass"
             call exit(1)
          end if
       end if
    end if
    if( .not. our_mass_flag ) then
       if( present(mass) ) then
          write(0,*) "write_lammps_data_file  mass, but not mass_flag=.true."
       end if
    end if

    our_vel_flag = .false.
    if( present( vel_flag ) ) then
       if( vel_flag ) then
          our_vel_flag = .true.
       end if
    end if

    if( filename .eq. "-" ) then
       out_unit = 6
    else
       out_unit = 99
       open(out_unit, file=trim(filename), status='new')
    end if

    header_length = len_trim(header%header)
!!    write(out_unit,"(a<header_length>)") trim(header%header)   g95 refuses
!!    write(out_unit,"(a80)") trim(header%header)
    write(out_unit,*) trim(header%header)
    write(out_unit,*)
    write(out_unit,"(i10, ' atoms')") header%num_atoms
    write(out_unit,"(i3, ' atom types')") header%num_species
    write(out_unit,*)

    if( our_bond_flag ) then
       write(out_unit,"(i10, ' bonds')") header%num_bonds
       write(out_unit,"(i3, ' bond types')") header%num_bond_types
       write(out_unit,*)
    end if
    if( our_angle_flag ) then
       write(out_unit,"(i10, ' angles')") header%num_angles
       write(out_unit,"(i3, ' angle types')") header%num_angle_types
       write(out_unit,*)
    end if


    write(out_unit,"(2f18.8,' xlo xhi')") header%perlb(1), header%perub(1)
    write(out_unit,"(2f18.8,' ylo yhi')") header%perlb(2), header%perub(2)
    write(out_unit,"(2f18.8,' zlo zhi')") header%perlb(3), header%perub(3)
    write(out_unit,*)

    if( our_mass_flag ) then
       write(out_unit,"('Masses')")
       write(out_unit,*)
       do type_ii = 1, header%num_species
          write(out_unit,"(i3,f18.8)") type_ii, mass(type_ii)
       end do
       write(out_unit,*)
    end if

    write(out_unit,"('Atoms')")
    write(out_unit,*)

    do atom_ii = 1, header%num_atoms
       if( our_mol_flag ) then
          if( our_charge_flag ) then
             write(out_unit,"(i10,i10,i3,4f18.8)")                  &
                   atom_table(atom_ii)%tag + 1,               &
                   atom_table(atom_ii)%mol + 1,               &
                   atom_table(atom_ii)%species + 1,           &
                   atom_table(atom_ii)%charge,                &
                   atom_table(atom_ii)%pos
          else
             write(out_unit,"(i10,i10,i3,3f18.8)")                  &
                   atom_table(atom_ii)%tag + 1,               &
                   atom_table(atom_ii)%mol + 1,               &
                   atom_table(atom_ii)%species + 1,           &
                   atom_table(atom_ii)%pos
          end if
       else
          if( our_charge_flag ) then
             write(out_unit,"(i10,i3,4f18.8)")                      &
                   atom_table(atom_ii)%tag + 1,               &
                   atom_table(atom_ii)%species + 1,           &
                   atom_table(atom_ii)%charge,                &
                   atom_table(atom_ii)%pos
          else
             write(out_unit,"(i10,i3,3f18.8)")                      &
                   atom_table(atom_ii)%tag + 1,               &
                   atom_table(atom_ii)%species + 1,           &
                   atom_table(atom_ii)%pos
          end if
       end if
    end do

    if( our_vel_flag ) then
       write(out_unit,*)
       write(out_unit,"('Velocities')")
       write(out_unit,*)
       do atom_ii = 1, header%num_atoms
          write(out_unit,"(i10,3f18.8)") atom_table(atom_ii)%tag + 1, atom_table(atom_ii)%vel
       end do
    end if

    if( our_bond_flag ) then
       write(out_unit,*)
       write(out_unit,"('Bonds')")
       write(out_unit,*)
       do bond_ii = 1, header%num_bonds
          write(out_unit,"(i10,i3,3i10)")                          &
                bond_table(bond_ii)%bond_tag,                &
                bond_table(bond_ii)%bond_type,               &
                bond_table(bond_ii)%atom_one,                &
                bond_table(bond_ii)%atom_two
       end do
    end if

    if( our_angle_flag ) then
       write(out_unit,*)
       write(out_unit,"('Angles')")
       write(out_unit,*)
       do angle_ii = 1, header%num_angles
          write(out_unit,"(i10,i3,3i10)")                          &
                angle_table(angle_ii)%angle_tag,             &
                angle_table(angle_ii)%angle_type,            &
                angle_table(angle_ii)%atom_one,              &
                angle_table(angle_ii)%atom_two,              &
                angle_table(angle_ii)%atom_three
       end do
    end if

    if( out_unit .eq. 99 ) then
       close(99)
    end if

    return
  end subroutine write_lammps_data_file


end module lammps_data_file_mod
