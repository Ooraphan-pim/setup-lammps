module x_dir_to_crystal_mod
contains

!!**********************************************************
!!               x_dir_to_crystal_ii
!!**********************************************************

  subroutine x_dir_to_crystal_ii( x_dir, crystal_ii )
    implicit none

    integer, intent(in)  :: x_dir          !  -1 or 1
    integer, intent(out) :: crystal_ii     !   1 or 2

    if( x_dir .eq. -1 ) then
       crystal_ii = 1
    else if( x_dir .eq. 1 ) then
       crystal_ii = 2
    else
       write(0,*) 'x_dir_to_crystal_ii: internal error'
       call abort()
    end if

    return
  end subroutine x_dir_to_crystal_ii

end module x_dir_to_crystal_mod
