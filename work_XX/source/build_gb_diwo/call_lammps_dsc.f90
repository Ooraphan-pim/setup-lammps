module call_lammps_dsc_mod
  use system_mod
  implicit none
 
!!  interface
!!     integer*4 function system(command_line)
!!       character(*), intent(in) :: command_line
!!     end function system
!!  end interface
 
contains

  subroutine call_lammps_dsc(name, lammps_command, num_out, energy_out)
    use delete_file_mod
    implicit none
    character(256), intent(in) :: name
    character(256), intent(in) :: lammps_command   !! e.g.  ./lmp_serial
    integer, intent(out) :: num_out
    real*8, intent(out) :: energy_out

    character(256) :: temp_name
    character(512) :: line
    integer :: result

    integer :: type_ii
    integer :: type_in, num_in
    real*8  :: energy_in
    character(256) :: label_in
    integer :: num_free
    real*8  :: energy_free

!!    line = "./lmp_serial -log "//trim(name)//".min.log -in "//trim(name)//".min.in"
    line = trim(lammps_command)                           &
           //" -log "//trim(name)//".min.log -in "//trim(name)//".min.in"
!! debugging
    write(0,*) trim(line)
!! end debugging
    result = system(trim(line)//char(0))
    if( result .ne. 0 ) then
       write(0,*) 'system->lmp_serial failed'
       call abort()
    end if

!!    line ="grep -A2 Step "//trim(name)//".min.log | tail -1 > "//trim(name)//".min.result"

    line ="grep 'Final_energy_by_type' "//trim(name)//".min.log > "//trim(name)//".min.result"
    result = system(trim(line)//char(0))
    if( result .ne. 0 ) then
       write(0,*) 'system->grep failed'
       call abort()
    end if
    temp_name = trim(name)//".min.result"
    open(99, file=trim(temp_name), status="old", action="read")
!!    read(99,*) r_step, r_e_1, r_e_2, r_p_xx, r_p_yy, r_p_zz
    num_free = 0
    energy_free = 0.0d0
    do type_ii = 1, 3  !! There are 5, we want the first 3
       read(99,*) label_in, type_in, num_in, energy_in 
       if( type_in .ne. type_ii ) then
          write(0,*) "call_lammps_dsc: unexpected type_in"
          call abort()
       end if
       num_free = num_free + num_in
       energy_free = energy_free + energy_in
    end do
    close(99)
    num_out = num_free
    energy_out = energy_free

    temp_name = trim(name)//".min.log"
    call delete_file(temp_name, len_trim(temp_name))
    temp_name = trim(name)//".min.result"
    call delete_file(temp_name, len_trim(temp_name))

    return
  end subroutine call_lammps_dsc

end module call_lammps_dsc_mod
