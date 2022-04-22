module scatter_randomly_g95_mod
  implicit none

contains

!! ********************************************************************
!!                         scatter_randomly_g95
!! ********************************************************************
!!      Picks _m_ items randomly out of _n_, without replacement.
!!
!!      If you wish to scatter _m_ solutes, say, in a matrix with
!!      _n_ sites, call as
!!          call scatter_randomly(_m_, _n_, _seed_, _array_name_(1:_m_))
!!      What is returned is a list of _m_ integers, each between 1 and _n_.
!!      If _m_ is greater than _n_/2, flip your problem when calling.
!!
!!      If _seed_ is not zero, it is used as a seed for the random number
!!      generator.
!!
!! A purely O(change) algorithm using only change calls to random_sub
!! is available from 
!! http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/59883
!! As implemented here:
!!    Allocate a list change_table(1:change)
!!    Make a list table(1:total) with table(j) = j.
!!    curr_len = total
!!    do ii = 1, change
!!      use random_sub to pick a number kk, 1 <= kk <= curr_len
!!      change_table(ii) = kk
!!      table(kk) = table(curr_len)
!!      curr_len = curr_len - 1
!!      !! (Notice that above is ok for kk = curr_len.)
!!    end do
!!

subroutine scatter_randomly_g95(change_count, matrix_count, seed, change_array)
  implicit none
  integer, intent(in)  :: change_count, matrix_count, seed
  integer, intent(out) :: change_array(change_count)

  integer :: matrix_list(matrix_count)
  integer :: ii, current_list_length, random_int

!! to be fixed
  if( seed .ne. 0 ) then
     write(0,*) 'scatter_randomly_g95: seed not yet implemented'
     call abort()
  end if


  if( change_count .gt. matrix_count / 2 ) then
     write(0,*) 'scatter_randomly: more solute than matrix, flip it.'
     call abort()
  end if

  do ii = 1, matrix_count
     matrix_list(ii) = ii
  end do
  current_list_length = matrix_count

  do ii = 1, change_count
     if( ii .eq. 1 ) then
        random_int = random_sub_g95(current_list_length, seed)
     else
        random_int = random_sub_g95(current_list_length, 0)
     end if
     change_array(ii) = matrix_list(random_int)
     matrix_list(random_int) = matrix_list(current_list_length)
     current_list_length = current_list_length - 1
  end do

  return

  contains
!! ================================================================
!!     If seed is non-zero, initializes seed as well as giving a random_number.
!!     After initialization, call with seed = 0

    function random_sub_g95(total, seed)
      implicit none
      
      integer, intent(in)   :: total
      integer, intent(in)   :: seed
      integer               :: random_sub_g95

      integer, parameter :: big = 2147483647    !! (2**31-1)
!!       logical, save :: first_call = .true.

!!      integer  pass

      integer :: temp_2
      real*4  :: random 

!!       integer, external :: irandm

!! to be fixed
      if( seed .ne. 0 ) then
         write(0,*) 'random_sub_g95: seed not yet implemented'
         call abort()
      end if
!! end of to be fixed

!!       pass = seed

      call random_number(random)

      temp_2 = ceiling(random * total)

      if( temp_2 .le. 0 ) then
         if( random .eq. 0.0 ) then
            temp_2 = 1
         else
            write(0,*) 'random_sub_g95: temp_2 <= 0 problem'
            call abort()
         end if
      else if( temp_2 .gt. total ) then
         write(0,*) 'random_sub_g95: temp_2 > total problem'
         call abort()
      end if

      random_sub_g95 = temp_2
      return
    end function random_sub_g95
  end subroutine scatter_randomly_g95

end module scatter_randomly_g95_mod
