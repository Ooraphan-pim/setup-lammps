module jiggle_position_mod
  implicit none

contains

!! changes position by amount times a random direction

  subroutine jiggle_position( position, amount )
    use random_direction_mod
    implicit none

    real*8, intent(inout) :: position(3)
    real*8, intent(in)    :: amount

    real*8  :: direction(3)

    call random_direction( direction )

    position = position + amount * direction

    return
  end subroutine jiggle_position

end module jiggle_position_mod
