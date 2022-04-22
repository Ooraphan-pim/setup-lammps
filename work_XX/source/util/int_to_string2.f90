!! ********************************************************************
!!                             module int_to_string2_mod
!! ********************************************************************


module int_to_string2_mod
  implicit none


contains




!! ************************************************************
!!           int_to_string2
!! ************************************************************
  function int_to_string2( int_in )
    implicit none

    character(256) :: int_to_string2
    integer, intent(in) :: int_in

    character(256) :: temp_string
    
    write(temp_string, *) int_in
    read(temp_string, *)  int_to_string2
    return
  end function int_to_string2
   





end module int_to_string2_mod

