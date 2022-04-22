!! build_gb_atom_data.f90

module build_gb_atom_data_mod
  implicit none

  integer, parameter :: max_atoms = 600000
  integer, parameter :: max_crystals = 2
  integer, parameter :: num_crystals = 2
  integer, parameter :: max_nbrs = 20

  type nbr_entry_type
     integer :: nbr_index
     real*8  :: dist_sq
  end type nbr_entry_type

  type :: atom_type
     integer :: crystal
     integer :: atom_type     !! species, only used for alloys
     real*8  :: pos(3)
     real*8  :: orig_pos(3)
!! for neighbors
     real*8  :: rel_pos(3)
     integer :: bin(3)
     integer :: next_atom        !! next atom in same bin
     integer :: num_nbrs
     type( nbr_entry_type) :: nbr_list(max_nbrs)
!! for cleaning up boundary
     logical :: deleted_atom   !! if true, it was deleted
     logical :: deleted_first_pass
  end type atom_type

  type( atom_type ) :: atom(max_atoms)


  integer num_atoms_built
  integer num_atoms_deleted
  integer actual_num_atoms

  integer orig_atoms_deleted


!!  integer num_atoms


end module build_gb_atom_data_mod
