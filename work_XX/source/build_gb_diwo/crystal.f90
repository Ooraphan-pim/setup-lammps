!! crystal.f90

module crystal_mod
  implicit none

  integer, private, parameter :: max_num_basis = 4
  integer, private, parameter :: max_half_nn   = 6

  type crystal_type
     character(3) :: name
     integer      :: num_basis
     real*8       :: basis(3, max_num_basis)
     integer      :: num_nn
!! This assumes the crystal has inversion symmetry.
     integer      :: num_half_nn
     real*8       :: half_nn(3, max_half_nn)
     integer      :: num_half_recip
     real*8       :: half_recip(3, max_half_nn)
  end type crystal_type
     
  real*8, private, parameter ::                        &
     fcc_basis(3,4) = reshape(                         &
            (/  0.0d0,  0.0d0,  0.0d0,                 &
                0.0d0,  0.5d0,  0.5d0,                 &
                0.5d0,  0.0d0,  0.5d0,                 &
                0.5d0,  0.5d0,  0.0d0   /),            &
            (/ 3, 4 /) )

         !! The other 6 are the opposites of these
  real*8, private, parameter ::                           &
     fcc_half_nn(3,6) = reshape(                          &
            (/  0.0d0,  0.5d0,  0.5d0,                    &
                0.0d0,  0.5d0, -0.5d0,                    &
                0.5d0,  0.0d0,  0.5d0,                    &
                0.5d0,  0.0d0, -0.5d0,                    &
                0.5d0,  0.5d0,  0.0d0,                    &
                0.5d0, -0.5d0,  0.0d0   /),               &
            (/ 3, 6 /) )
  real*8, private, parameter ::                           &
       fcc_half_recip(3,4) = reshape(                    &
            (/  1.0d0,  1.0d0,  1.0d0,                   &
               -1.0d0,  1.0d0,  1.0d0,                   &
                1.0d0, -1.0d0,  1.0d0,                   &
                1.0d0,  1.0d0, -1.0d0   /),              &
            (/ 3, 4 /) )

  real*8, private, parameter ::                        &
     bcc_basis(3,2) = reshape(                         &
            (/  0.0d0,  0.0d0,  0.0d0,                 &
                0.5d0,  0.5d0,  0.5d0   /),            &
            (/ 3, 2 /) )

         !! The other 4 are the opposites of these
  real*8, private, parameter ::                           &
     bcc_half_nn(3,4) = reshape(                          &
            (/  0.5d0,  0.5d0,  0.5d0,                    &
                0.5d0,  0.5d0, -0.5d0,                    &
                0.5d0, -0.5d0,  0.5d0,                    &
                0.5d0, -0.5d0, -0.5d0   /),               &
            (/ 3, 4 /) )

  real*8, private, parameter ::                            &
       bcc_half_recip(3, 6) = 2.0d0 * fcc_half_nn


contains

  subroutine check_crystal_name(crystal, result)
    implicit none

    type(crystal_type) :: crystal
    !!    character(3), intent(in) :: crystal
    logical, intent(out)     :: result

!! debugging
!!$    integer ii
!! end debugging

    if( crystal%name .eq. 'FCC' ) then
       crystal%num_basis = 4;
       crystal%basis(:,1:4) = fcc_basis;
       crystal%num_nn = 12;
       crystal%num_half_nn = 6;
       crystal%half_nn(:,1:6) = fcc_half_nn;
       crystal%num_half_recip = 4;
       crystal%half_recip(:,1:4) = fcc_half_recip
!! debugging
!!$       do ii = 1, crystal%num_basis
!!$          write(0,*) crystal%basis(1:3, ii)
!!$       end do
!! end debugging
    else if( crystal%name .eq. 'BCC' ) then
       crystal%num_basis = 2;
       crystal%basis(:,1:2) = bcc_basis;
       crystal%num_nn = 8;
       crystal%num_half_nn = 4;
       crystal%half_nn(:,1:4) = bcc_half_nn;
       crystal%num_half_recip = 6;
       crystal%half_recip(:,1:6) = bcc_half_recip
!! debugging
!!$       do ii = 1, crystal%num_basis
!!$          write(0,*) crystal%basis(1:3, ii)
!!$       end do
!! end debugging
    else
       write(0,*) 'check_crystal_name: invalid crystal, should be "FCC"'
       call abort()
    end if

    result = .true.

    return
  end subroutine check_crystal_name


  subroutine crystal_in_lattice( vector, crystal, result )
    implicit none

    integer, intent(in)      :: vector(3)
    character(3), intent(in) :: crystal
    logical, intent(out)     :: result

    integer :: temp, temp_sum, ii

    if( crystal .eq. "FCC" ) then
       temp = vector(1) + vector(2) + vector(3)
       if( modulo(temp, 2) .eq. 0 ) then
          result = .true.
       else
          result = .false.
       end if
    else if( crystal .eq. "BCC" ) then
       temp_sum = 0
       do ii = 1, 3
          temp = modulo( vector(ii), 2 )   !! either zero or positive one
          temp_sum = temp_sum + temp
       end do
       if( ( temp_sum .eq. 0 )                     &
            .or. ( temp_sum .eq. 3 ) ) then
          result = .true.
       else
          result = .false.
       end if
    else
       write(0,*) 'crystal_in_lattice: invalid crystal, should be "FCC"'
       call abort()
    end if

    return
  end subroutine crystal_in_lattice

end module crystal_mod
