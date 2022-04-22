module random_direction_mod
  implicit none


contains

  subroutine random_direction(direction)
    use random_numbers_mod
    implicit none

    real*8, intent(out) :: direction(3)

    real*8  :: rand_zero_one
    real*8  :: zz, temp, remainder, theta
    integer :: dim_ii

    real*8  :: pi, two_pi

    pi = acos(-1.0d0)
    two_pi = 2.0d0 * pi

    rand_zero_one = random_real()

    zz = 2.0d0 * rand_zero_one - 1.0d0

    temp = 1.0d0 - zz**2
    if( abs( temp ) .lt. 1.0d-12 ) then
       remainder = 0.0d0
    else
       remainder = sqrt( temp )
    end if

    rand_zero_one = random_real()
    theta = two_pi * rand_zero_one

    direction(1) = remainder * cos( theta )
    direction(2) = remainder * sin( theta )
    direction(3) = zz

    temp = 0.0d0
    do dim_ii = 1, 3
       temp = temp + direction(dim_ii)**2
    end do

    if( abs( 1.0d0 - temp ) .gt. 1.0d-8 ) then
       write(0,*) "not unit  ", direction, temp
       call exit(1)
    end if

    return
  end subroutine random_direction

end module random_direction_mod
