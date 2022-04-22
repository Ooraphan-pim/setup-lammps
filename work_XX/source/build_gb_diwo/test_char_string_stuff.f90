

subroutine  cpp_read_line( line, count )
  implicit none
  character(256), intent(out) :: line
  integer, intent(out) :: count

  integer :: read_status

  read(5,"(a256)",iostat=read_status) line
  if( read_status .gt. 0 ) then
     write(0,*) "cpp_read_line  read_error"
     count = -1
  else if( read_status .lt. 0 ) then
     count = -1
  else
!! debugging
!!    write(0,*) "cpp_read_line  ", trim(line)
!! end debugging
     count = len_trim(line)
     if( count .eq. 0 ) then
        count = -1
     end if
  end if

  return
end subroutine cpp_read_line


subroutine cpp_write_line(line)
  implicit none
  character(256), intent(in) :: line

  write(6,*) trim(line)

  return
end subroutine cpp_write_line

