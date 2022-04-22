!! Gather up the results from a set of runs



program gb_energy_gather
  use gb_energy_input_dsc_mod
  use int_to_string2_mod
  use delete_file_mod
  use system_mod
  implicit none

  character(256) :: control_name
  integer :: read_status
  integer :: pbv_num
  character(256) :: ignore_string
  integer :: num_init_in
  integer :: num_result_read
  character(256) :: name_with_pbv
  character(256) :: name_with_dsc_num
  character(256) :: name_with_dsc_vec
  character(256) :: name_with_place
  character(256) :: name_with_method
  character(256) :: final_restart_name
  character(256) :: method_restart_name
  character(256) :: dsc_restart_name
  character(256) :: final_data_name
  character(256) :: method_data_name
  character(256) :: dsc_data_name
  character(256) :: line
  integer :: system_result, open_status

  logical :: pbv_first_time
  integer :: method_ii, skip_ii
  logical :: dist_valid

  real*8  :: gb_best_energy(3)
  real*8  :: gb_best_dist(3)
  integer :: gb_best_atoms(3)
  integer :: gb_best_place(3)
  integer :: gb_best_dsc(3, 3)
  real*8  :: best_energy
  real*8  :: best_dist
  integer :: best_place
  integer :: best_atoms
  integer :: best_dsc(3)
  integer :: best_method

  integer :: dsc_int_vec(3)

  call get_gb_energy_input_dsc(.false., 1)      !! single_crystal_flag

  if( delete_by_energy .or. delete_by_pressure ) then
     dist_valid = .false.
  else
     dist_valid = .true.
  end if


  control_name = trim(main_name)//"."//trim(adj_name)//".control"
  open(22, file=control_name, status="old", action="read")

  pbv_loop:  &
  do
     read(22, *, iostat = read_status) ignore_string, pbv_num 
     if( read_status .gt. 0 ) then
        write(0,*) "gb_energy_start  read error"
        call abort()
     else if( read_status .lt. 0 ) then
        exit pbv_loop
     end if
     do skip_ii = 1, 3
        read(22,*)  !!  The pbv lines
     end do
     read(22,*) ignore_string, num_init_in
     name_with_pbv = trim(main_name)//"."//trim(int_to_string2(pbv_num))
     name_with_dsc_num = trim(name_with_pbv)//"."//trim(int_to_string2(dsc_num))  
     num_result_read = 0
     pbv_first_time = .true.
     call loop_over_dsc_vectors()
     if( num_result_read .ne. num_init_in ) then
        write(0,*) "gb_energy_gather  did not gather the expected number of runs"
        write(0,*) "  pbv_num        ", pbv_num
        write(0,*) "  num_init_in    ", num_init_in
        write(0,*) "  num_result_read  ", num_result_read
        call abort()
     end if
     do method_ii = 1, 3
        if( do_method(method_ii) ) then
           if( dist_valid ) then
              write(6,"('PBV ', i6, '  Method ', i2, 3i6, i8, f8.4, i9, f9.3)")    &
                   pbv_num, method_ii, gb_best_dsc(1:3, method_ii),               &
                   gb_best_place(method_ii),                                      &
                   gb_best_dist(method_ii),                                       &
                   gb_best_atoms(method_ii),                                      &
                   gb_best_energy(method_ii)
           else
              write(6,"('PBV ', i6, '  Method ', i2, 3i6, i8, i9, f9.3)")        &
                   pbv_num, method_ii, gb_best_dsc(1:3, method_ii),               &
                   gb_best_place(method_ii),                                      &
                   gb_best_atoms(method_ii),                                      &
                   gb_best_energy(method_ii)
           end if
           method_restart_name = trim(name_with_dsc_num)               &
                                 //trim(int_to_string2(method_ii))     &
                                 //".method.restart"
           method_data_name = trim(name_with_dsc_num)               &
                              //trim(int_to_string2(method_ii))     &
                              //".method.data"
           if( pbv_first_time .or. (gb_best_energy(method_ii) < best_energy ) ) then
              pbv_first_time = .false.
              best_energy = gb_best_energy(method_ii)
              best_dist = gb_best_dist(method_ii)
              best_place = gb_best_place(method_ii)
              best_atoms = gb_best_atoms(method_ii)
              best_dsc = gb_best_dsc(1:3, method_ii)
              best_method = method_ii
              final_restart_name = trim(name_with_dsc_num)               &
                                   //".min.restart"
              final_data_name = trim(name_with_dsc_num)               &
                                //".min.data"
              line = "cp -f "                              &
                     //trim(method_restart_name)//" "      &
                     //trim(final_restart_name)
              system_result = system(trim(line)//char(0))
              line = "cp -f "                              &
                     //trim(method_data_name)//" "      &
                     //trim(final_data_name)
              system_result = system(trim(line)//char(0))
           end if
           open(99, file=method_restart_name, status="old", action="read", iostat=open_status)
           close(99)
           if( open_status .eq. 0 ) then
              call delete_file(trim(method_restart_name), len_trim(method_restart_name))
           end if
           open(99, file=method_data_name, status="old", action="read", iostat=open_status)
           close(99)
           if( open_status .eq. 0 ) then
              call delete_file(trim(method_data_name), len_trim(method_data_name))
           end if
        end if   !! if( do_method(method_ii)
     end do      !! method_ii
     write(6,*)
     if( dist_valid ) then
        write(6,"('PBV', i6, ' DSC', i2, ' best', i2, '   ', 3i3, i8, f18.14, i9, f9.3)")    &
          pbv_num, dsc_num, best_method, best_dsc,                          &
          best_place, best_dist, best_atoms, best_energy
     else
        write(6,"('PBV', i6, ' DSC', i2, ' best', i2, '   ', 3i3, i8, i9, f9.3)")    &
          pbv_num, dsc_num, best_method, best_dsc,                          &
          best_place, best_atoms, best_energy
     end if

     

  end do pbv_loop

  stop

contains

  subroutine do_one_placement(place_label_ii,              &
                              success, gb_energy_result,   &
                              gb_dist, gb_atoms )
    implicit none
    integer, intent(in)  :: place_label_ii
    logical, intent(out) :: success
    real*8, intent(out)  :: gb_energy_result(3)
    real*8, intent(out)  :: gb_dist(3)
    integer, intent(out) :: gb_atoms(3)

    character(256) :: init_name
    character(256) :: result_name
    character(256) :: scr_name
    character(256) :: sub_sh_name
    character(256) :: temp_restart_name
    character(256) :: temp_name
    integer :: open_status
    integer :: method_ii
    integer :: jobid_in

    logical :: success_all
    logical :: success_some

    success_some = .false.
    success_all = .true.

    method_loop:  &
    do method_ii = 1, 3
       if( do_method(method_ii) ) then
          name_with_method = trim(name_with_place)//"."//trim(int_to_string2(method_ii))
          init_name = trim(name_with_method)//".init"
          result_name = trim(name_with_method)//".result"
          scr_name = "temp."//trim(name_with_method)//".scr"
          sub_sh_name = "temp."//trim(name_with_method)//".sub.sh"
          temp_restart_name = "temp."//trim(name_with_method)//".temp.restart"
          !! debugging
!!$          write(0,*) trim(init_name)
!!$          stop
          !! end debugging

          open(99, file=init_name, status="old", action="read", iostat=open_status)
          close(99)
          if( open_status .ne. 0 ) then
             success_all = .false.
             cycle method_loop
          end if
          open(99, file=result_name, status="old", action="read", iostat=open_status)
          if( open_status .ne. 0 ) then
             success_all = .false.
             cycle method_loop
          end if
          success_some = .true.
          read(99, *) gb_dist(method_ii), gb_atoms(method_ii), gb_energy_result(method_ii) 
          close(99)
          num_result_read = num_result_read + 1
          if( dist_valid ) then
             write(6,"('GB energy(',i1,')  ', 4i6, i8, f8.4, i9, f9.3)")      &
                  method_ii, pbv_num, dsc_int_vec, place_label_ii,            &
                  gb_dist(method_ii),                                         &
                  gb_atoms(method_ii),                                        &
                  gb_energy_result(method_ii)
          else
             write(6,"('GB energy(',i1,')  ', 4i6, i8, i9, f9.3)")            &
                  method_ii, pbv_num, dsc_int_vec, place_label_ii,            &
                  gb_atoms(method_ii),                                        &
                  gb_energy_result(method_ii)
          end if

          open(99, file=scr_name, status="old", action="read", iostat=open_status)
          close(99)
          if( open_status .eq. 0 ) then
             call delete_file(trim(scr_name), len_trim(scr_name))
          end if
          open(99, file=sub_sh_name, status="old", action="read", iostat=open_status)
          close(99)
          if( open_status .eq. 0 ) then
             call delete_file(trim(sub_sh_name), len_trim(sub_sh_name))
          end if
!!          open(99, file=temp_restart_name, status="old", action="read", iostat=open_status)
!!          close(99)
!!          if( open_status .eq. 0 ) then
!!             call delete_file(trim(temp_restart_name), len_trim(temp_restart_name))
!!          end if

          temp_name = "temp."//trim(name_with_method)//".jobid"
          open(99, file=temp_name, status="old", iostat=open_status)
          if( open_status .eq. 0 ) then
             read(99, *) jobid_in
             close(99, status='delete')
             temp_name = "test_"//trim(name_with_method)//".o"         &
                                //trim(int_to_string2(jobid_in))
             open(99, file=temp_name, status="old",iostat=open_status)
             if( open_status .eq. 0 ) then
                close(99, status='delete')
             else
                close(99)
             end if
             temp_name = "test_"//trim(name_with_method)//".e"         &
                                //trim(int_to_string2(jobid_in))
             open(99, file=temp_name, status="old",iostat=open_status)
             if( open_status .eq. 0 ) then
                close(99, status='delete')
             else
                close(99)
             end if
          else
             close(99)
          end if
       end if    !! if( do_method(method_ii) )
    end do method_loop

    if( success_all ) then
       if( .not. success_some ) then
          write(0,*) "gb_energy_start  did not do any methods"
          call abort()
       else
          success = .true.
       end if
    else
       if( success_some ) then
          write(0,*) "gb_energy_start  found some methods, but not all"
          call abort()
       else
          success = .false. 
       end if
    end if

    return
  end subroutine do_one_placement
       

  subroutine loop_over_placements(best_gb_energy,         &
                                  best_gb_dist,           &
                                  best_gb_atoms,          &
                                  best_place_index)
    implicit none

    real*8,  intent(out)  :: best_gb_energy(3)
    real*8,  intent(out)  :: best_gb_dist(3)
    integer, intent(out)  :: best_gb_atoms(3)
    integer, intent(out)  :: best_place_index(3)

    real*8  :: gb_energy_one_place(3)
    real*8  :: gb_dist_one_place(3)
    integer :: gb_atoms_one_place(3)

    character(256) :: min_restart_name
    character(256) :: min_data_name
    integer place_index, mm_ii
    logical result
    logical first_time

    first_time = .true.
    place_index = 0
    place_loop:  &
    do
       place_index = place_index + 1
       name_with_place = trim(name_with_dsc_vec)//"."//trim(int_to_string2(place_index))
       call do_one_placement(place_index,             &
                             result,                  &
                             gb_energy_one_place,     &
                             gb_dist_one_place,       &
                             gb_atoms_one_place )
       if( .not. result ) then
          exit place_loop
       end if

       do mm_ii = 1, 3
          if( do_method(mm_ii) ) then
             min_restart_name = "temp."//trim(name_with_place)//"."       &
                                //trim(int_to_string2(mm_ii))     &
                                //".min.restart"
             min_data_name = "temp."//trim(name_with_place)//"."       &
                                //trim(int_to_string2(mm_ii))     &
                                //".data.keep"
             if( first_time .or. (gb_energy_one_place(mm_ii) < best_gb_energy(mm_ii)) ) then
                best_gb_energy(mm_ii) = gb_energy_one_place(mm_ii)
                best_gb_dist(mm_ii) = gb_dist_one_place(mm_ii)
                best_gb_atoms(mm_ii) = gb_atoms_one_place(mm_ii)
                best_place_index(mm_ii) = place_index
                dsc_restart_name = trim(name_with_dsc_vec)               &
                                      //trim(int_to_string2(mm_ii))     &
                                      //".dsc.restart"
                dsc_data_name = trim(name_with_dsc_vec)               &
                                      //trim(int_to_string2(mm_ii))     &
                                      //".dsc.data"
                line = "cp -f "                            &
                     //trim(min_restart_name)//" "              &
                     //trim(dsc_restart_name)
                system_result = system(trim(line)//char(0))
                line = "cp -f "                            &
                     //trim(min_data_name)//" "              &
                     //trim(dsc_data_name)
                system_result = system(trim(line)//char(0))
             end if
!!             open(99, file=min_restart_name, status="old", action="read", iostat=open_status)
!!             close(99)
!!             if( open_status .eq. 0 ) then
!!                call delete_file(trim(min_restart_name), len_trim(min_restart_name))
!!             end if
             open(99, file=min_data_name, status="old", action="read", iostat=open_status)
             close(99)
             if( open_status .eq. 0 ) then
                call delete_file(trim(min_data_name), len_trim(min_data_name))
             end if
          end if
       end do

       if( first_time ) then
          first_time = .false.
       end if
       
    end do place_loop

    do mm_ii = 1, 3
       if( do_method(mm_ii) ) then
          if( dist_valid ) then
             write(6,"('PBV ', i6, '  Method ', i2, '  DSC ', 3i6, i8, f8.4, i9, f9.3)") &
                  pbv_num, mm_ii, dsc_int_vec,                           &
                  best_place_index(mm_ii), best_gb_dist(mm_ii),          &
                  best_gb_atoms(mm_ii), best_gb_energy(mm_ii)
          else
             write(6,"('PBV ', i6, '  Method ', i2, '  DSC ', 3i6, i8, i9, f9.3)") &
                  pbv_num, mm_ii, dsc_int_vec,                           &
                  best_place_index(mm_ii),                               &
                  best_gb_atoms(mm_ii), best_gb_energy(mm_ii)
          end if
       end if
    end do


    return
  end subroutine loop_over_placements


  subroutine loop_over_dsc_vectors()
    implicit none

    logical :: first_time
    real*8  :: gb_energy_one_dsc(3)
    real*8  :: gb_dist_one_dsc(3)
    integer :: gb_atoms_one_dsc(3)
    integer :: gb_place_one_dsc(3)

    integer :: dsc_ii, dsc_jj, dsc_kk

    first_time = .true.
    do dsc_ii = 0, dsc_num - 1
       dsc_int_vec(1) = dsc_ii
       do dsc_jj = 0, dsc_num - 1
          dsc_int_vec(2) = dsc_jj
          do dsc_kk = 0, dsc_num - 1
             dsc_int_vec(3) = dsc_kk
             name_with_dsc_vec = trim(name_with_dsc_num)                &
                                 //"."//trim(int_to_string2(dsc_ii))    &
                                 //"."//trim(int_to_string2(dsc_jj))    &
                                 //"."//trim(int_to_string2(dsc_kk))

             !! debugging
!!$             write(0,*) trim(name_with_dsc_vec)
!!$             stop
             !! end debugging

             call loop_over_placements(gb_energy_one_dsc,     &
                                       gb_dist_one_dsc,       &
                                       gb_atoms_one_dsc,      &
                                       gb_place_one_dsc)

             

             do method_ii = 1, 3
                if( do_method(method_ii) ) then
                   dsc_restart_name = trim(name_with_dsc_vec)               &
                                      //trim(int_to_string2(method_ii))     &
                                      //".dsc.restart"
                   dsc_data_name = trim(name_with_dsc_vec)               &
                                      //trim(int_to_string2(method_ii))     &
                                      //".dsc.data"
                   if( first_time .or. (gb_energy_one_dsc(method_ii) < gb_best_energy(method_ii)) ) then
                      gb_best_energy(method_ii) = gb_energy_one_dsc(method_ii)
                      gb_best_dist(method_ii) = gb_dist_one_dsc(method_ii)
                      gb_best_atoms(method_ii) = gb_atoms_one_dsc(method_ii)
                      gb_best_place(method_ii) = gb_place_one_dsc(method_ii)
                      gb_best_dsc(1:3, method_ii) = dsc_int_vec
                      method_restart_name = trim(name_with_dsc_num)               &
                                            //trim(int_to_string2(method_ii))     &
                                            //".method.restart"
                      method_data_name = trim(name_with_dsc_num)               &
                                         //trim(int_to_string2(method_ii))     &
                                         //".method.data"
                      line = "cp -f "                            &
                           //trim(dsc_restart_name)//" "              &
                           //trim(method_restart_name)
                      system_result = system(trim(line)//char(0))
                      line = "cp -f "                            &
                           //trim(dsc_data_name)//" "              &
                           //trim(method_data_name)
                      system_result = system(trim(line)//char(0))
                   end if
                   open(99, file=dsc_restart_name, status="old", action="read", iostat=open_status)
                   close(99)
                   if( open_status .eq. 0 ) then
                      call delete_file(trim(dsc_restart_name), len_trim(dsc_restart_name))
                   end if
                   open(99, file=dsc_data_name, status="old", action="read", iostat=open_status)
                   if( open_status .eq. 0 ) then
                      call delete_file(trim(dsc_data_name), len_trim(dsc_data_name))
                   end if
                end if
             end do

             if( first_time ) then
                first_time = .false.
             end if    !! else of if( first_time )

          end do       !! do dsc_kk
       end do          !! do dsc_jj
    end do             !! do dsc_ii

    return
  end subroutine loop_over_dsc_vectors

end program gb_energy_gather
