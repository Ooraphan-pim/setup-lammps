!! Probably the routines to be called by C++ had better be outside the module

function diwo_start_build(name_in, name_in_len, out_name, out_name_len)
  use build_gb_global_for_diwo_mod
  use build_gb_run_out_mod
  implicit none
  character(256), intent(in) :: name_in
  integer, intent(in) :: name_in_len
  character(256), intent(out) :: out_name
  integer, intent(out) :: out_name_len
  integer :: diwo_start_build

  logical :: delete_by_e_or_p
  character(256) :: our_name
  character(256) :: run_out_name
  integer :: name_len

  if( name_in_len > 256 ) then
     write(0,*) "diwo_start_build  name_in_len > 256  name_in_len  ", name_in_len
     call abort()
     name_len = 256        !! liberty is putting a null in the string?
  else
     name_len = name_in_len
  end if
  our_name = name_in(1:name_len)

  run_out_name = trim(our_name)//".f90.run.out"
  open(run_out_unit, file=trim(run_out_name), status = "new")
  call diwo_read_build_gb_input_dsc(our_name)
  call diwo_temp_name()
!!    call diwo_sed_control_file()    !! task moved to gb_energy_cc_diwo
  call diwo_start_build_sub(delete_by_e_or_p)

  if( delete_by_e_or_p ) then
     diwo_start_build = 1
  else
     diwo_start_build = 0
  end if

  out_name = " "
  out_name = trim(global_temp_name)
  out_name_len = len_trim(global_temp_name)

  return
end function diwo_start_build


!! Result of 1 means finished, 0 otherwise
function diwo_continue_build_not_e_or_p()
  use build_gb_global_for_diwo_mod
  implicit none

  integer :: diwo_continue_build_not_e_or_p

  integer :: result
  
  call diwo_continue_build_not_e_or_p_sub(result)

  diwo_continue_build_not_e_or_p = result

  return
end function diwo_continue_build_not_e_or_p
