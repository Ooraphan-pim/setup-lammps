!! ********************************************************************
!!                             module delete_file
!! ********************************************************************


module delete_file_mod
implicit none

contains


!! ********************************************************************
!!                  delete_file
!! ********************************************************************
!!
!!  call as "call delete_file(__name__, len_trim(__name__))
!!

subroutine delete_file(delete_name, name_len)
  implicit none
  
  character*(*), intent(in) :: delete_name
  integer, intent(in)       :: name_len

!!  character*256             :: temp_name
!! debugging
!!  write(0,*) "delete_file  entry  delete_name:  ", trim(delete_name)
!! end debugging
!!  if( name_len .gt. 256 ) then
!!     write(0,*) "delete_file  name to long"
!!     call exit(1)
!!  end if
!!  temp_name = delete_name(1:name_len)
!! debugging
!!  write(0,*) "delete_file  before open  temp_name:  ", trim(temp_name)
!! end debugging

  open(99, file=trim(delete_name))
  close(99, status='delete')

!! debugging
!!  write(0,*) "delete_file  after close"
!! end debugging

  return
end subroutine delete_file


end module delete_file_mod

