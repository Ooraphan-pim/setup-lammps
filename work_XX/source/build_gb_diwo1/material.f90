!! material.f90

module material_mod
  use crystal_mod
  implicit none

  type material
     real*8 :: alat         ! In angstroms
     type(crystal_type) :: crystal      ! just "FCC" for now
  end type material

end module material_mod
