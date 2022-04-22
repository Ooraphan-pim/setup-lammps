module system_mod
  implicit none

   interface
      integer*4 function system(command_line)
        character(*), intent(in) :: command_line
      end function system
   end interface

end module system_mod

