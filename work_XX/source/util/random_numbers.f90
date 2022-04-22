module random_numbers_mod
  implicit none

  private

  public :: random_real, init_random_seed, init_random_from_seed

  logical :: initialized = .false.

  integer, parameter :: big = 2147483647    !! (2**31-1)


contains

  function  random_integer()
    implicit none

    integer :: random_integer

    real*8  :: random_num

    if( .not. initialized ) then
       write(0,*) "call init_random_seed or init_random_from_seed before calling random_integer"
       call exit(1)
    end if

    random_num = random_real()
    random_integer = random_num * big

  end function random_integer


!!     returns a pseudo-ranom number 0 <= random_real <= 1
!!  either init_random_seed()  or  init_random_from_seed  should be called first
!!
!!     random_real will call init_random_seed to do initialization to support existing code

  function  random_real()     
    implicit none

    real*8  :: random_real

    if( .not. initialized ) then
       call init_random_seed()
    end if

    call random_number( random_real )

  end function random_real


  SUBROUTINE init_random_from_seed( seed_in )
    implicit none

    integer, intent(in) :: seed_in

    INTEGER :: ii, n
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    real*8  :: ignore;

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    do ii = 1, n
       seed(ii) = seed_in + 37 * ii
    end do

    CALL RANDOM_SEED(PUT = seed)

!!    write(0,*) "seed  ", seed

    DEALLOCATE(seed)

    initialized = .true.

!!     prime by calling twice.  (Presumably not needed, but no harm either.)
    ignore = random_real()
    ignore = random_real()
    

    return
  END SUBROUTINE init_random_from_seed


  SUBROUTINE init_random_seed()
    implicit none

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

!!    write(0,*) "seed  ", seed

    DEALLOCATE(seed)

    initialized = .true.

    return
  END SUBROUTINE init_random_seed

end module random_numbers_mod
