module get_lammps_result_mod
  implicit none
 
!!$  interface
!!$     integer*4 function system(command_line)
!!$       character(*), intent(in) :: command_line
!!$     end function system
!!$  end interface


!! Because of problems compiling with Intel compilers, the change was made:
!!       call delete_file(trim(_name_), len_trim(_name_))
!!  -->  call delete_file(_name_, len_trim(_name_))


 
contains

  subroutine get_lammps_result( name, delete_log_flag, double_edge_flag,                    &
                                num_out, energy_out,                                        &
                                num_with_out, energy_with_out              )
    use delete_file_mod
    implicit none
    character(256), intent(in) :: name
    logical, intent(in) :: delete_log_flag
    logical, intent(in) :: double_edge_flag
    integer, intent(out) :: num_out
    real*8, intent(out) :: energy_out
    integer, intent(out) :: num_with_out          !! values including the inner edge atoms
    real*8, intent(out) :: energy_with_out

    character(256) :: temp_name
    character(512) :: line
    integer :: result

    integer :: type_ii
    integer :: type_in, num_in
    real*8  :: energy_in
    character(256) :: label_in
    integer :: num_free, num_lines_added, read_status
    integer :: num_with
    real*8  :: energy_free
    real*8  :: energy_with

    temp_name = trim(name)//".min.log"
    open(99, file=temp_name, status="old", action="read")
    num_free = 0
    energy_free = 0.0d0
    num_with = 0
    energy_with = 0.0d0
    num_lines_added = 0
!! There are five lines labelled 'Final_energy_by_type' we want the first three
!! (If double_edge_flag there are 7 and we need the first three and the first five.)
    read_loop: do
       read(99,*,iostat=read_status) label_in, type_in, num_in, energy_in 
       if( read_status .gt. 0 ) then
          cycle read_loop         !! not our type of line
       else if( read_status .lt. 0 ) then
          write(0,*) "unexpected eof on ___.min.log"
          call abort()
       end if
       if( label_in .eq. "Final_energy_by_type" ) then
          num_lines_added = num_lines_added + 1
          if( type_in .ne. num_lines_added ) then
             write(0,*) "get_lammps_result_diwo: unexpected type_in"
             call abort()
          end if
          if( num_lines_added .le. 3 ) then 
             num_free = num_free + num_in
             energy_free = energy_free + energy_in
          end if
          if( double_edge_flag .and. num_lines_added .le. 5 ) then
             num_with = num_with + num_in
             energy_with = energy_with + energy_in
          end if
       end if
       if( double_edge_flag ) then
          if( num_lines_added .eq. 5 ) then
             exit read_loop
          end if
       else
          if( num_lines_added .eq. 3 ) then
             exit read_loop
          end if
       end if
    end do read_loop

    close(99)
    num_out = num_free
    energy_out = energy_free
    num_with_out = num_with
    energy_with_out = energy_with

    if( delete_log_flag ) then
       temp_name = trim(name)//".min.log"
       call delete_file(temp_name, len_trim(temp_name))
    end if

    return
  end subroutine get_lammps_result

end module get_lammps_result_mod
