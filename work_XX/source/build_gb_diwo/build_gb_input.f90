!! build_gb_input.f90


module build_gb_input_mod
  use build_gb_run_out_mod
  implicit none

contains

  subroutine read_build_gb_input(in_unit, init, dsc_flag)
    use build_gb_param_mod
    use build_gb_control_mod
    implicit none

    integer, intent(in) :: in_unit
    type(build_gb_control_type), intent(out) :: init
    logical, intent(in) :: dsc_flag

    character(256) :: name
    real*8  :: min_dist             !! atoms closer than min_dist are removed
    integer :: crystal_to_delete

    real*8 :: alat
    character(3) :: crystal
    character(80) :: deletion_method
    integer :: left_xx_pbv(3), left_yy_pbv(3), left_zz_pbv(3)
    integer :: right_xx_pbv(3), right_yy_pbv(3), right_zz_pbv(3)

    logical :: fixed_x_len
    real*8  :: left_x_len
    real*8  :: right_x_len

    real*8  :: edge_size  !! width of fixed regions, for dsc version only

    real*8  :: nbr_delta  !! How different do nbr distances need to be to be
                          !! different deletion passes.

    logical :: nbrs_by_dist
    logical :: two_boundaries
    logical :: sigma_flag

    real*8  :: concentration

!! concentration is optional.  If it is absent or zero, pure material.
!! Otherwise the concentration of species 2.   

!! These entries are for _dsc, "do_later"
    character(256) :: main_name
    integer :: pbv_num
    real*8  :: low_min_dist
    real*8  :: e_per_a
    integer :: dsc_num
    integer :: dsc_int_vec(3)
    integer :: place_index
    integer :: method_num
    real*8  :: dsc_vec(3)
    real*8  :: dsc_xx_adj
    real*8  :: xx_boundary_adj
    logical :: one_at_a_time
    logical :: lammps_flag
    logical :: output_flag
    character(256) :: lammps_control_file
    logical :: delete_by_energy
    logical :: show_pressure
    logical :: delete_by_pressure
    character(256) :: deletion_control_file
    character(256) :: lammps_command
    character(256) :: scr_file
    character(256) :: sub_sh_file
    logical        :: double_edge
    logical        :: jiggle_atoms
    real*8         :: jiggle_amount
    real*8         :: jiggle_region
    integer        :: jiggle_seed

    namelist /material_card/  alat,  crystal
    namelist /pbv_card/  left_xx_pbv, left_yy_pbv, left_zz_pbv,     &
                         right_xx_pbv, right_yy_pbv, right_zz_pbv
    namelist /control_card/ name,                                   &
                            min_dist, deletion_method,              &
                            crystal_to_delete,                      &
                            nbrs_by_dist,                           &
                            two_boundaries,                         &
                            concentration,                          &
                            sigma_flag,                             &
                            edge_size,                              &
                            nbr_delta,                              &
                            fixed_x_len,                            &
                            left_x_len,                             &
                            right_x_len,                            &
!! These entries are for _dsc, "do_later"
                            main_name,                              &
                            pbv_num,                                &
                            low_min_dist,                           &
                            e_per_a,                                &
                            dsc_num,                                &
                            dsc_int_vec,                            &
                            place_index,                            &
                            method_num,                             &
                            dsc_vec,                                &
                            dsc_xx_adj,                             &
                            xx_boundary_adj,                        &
                            one_at_a_time,                          &
                            lammps_flag,                            &
                            output_flag,                            &
                            lammps_control_file,                    &
                            delete_by_energy,                       &
                            show_pressure,                          &
                            delete_by_pressure,                     &
                            deletion_control_file,                  &
                            lammps_command,                         &
                            scr_file,                               &
                            sub_sh_file,                            &
                            double_edge,                            &
                            jiggle_atoms,                           &
                            jiggle_amount,                          &
                            jiggle_region,                          &
                            jiggle_seed


    call initialize_build_gb_control(init)
    name = init%name
    alat = init%alat
    crystal = init%crystal
    crystal_to_delete = init%crystal_to_delete
    nbrs_by_dist = init%nbrs_by_dist
    two_boundaries = init%two_boundaries
    concentration = init%concentration
    sigma_flag = init%sigma_flag
    edge_size = init%edge_size
    nbr_delta = init%nbr_delta
    lammps_flag = init%lammps_flag         !! .false.
    output_flag = init%output_flag         !! .true.

    fixed_x_len = init%fixed_x_len
    left_x_len = init%left_x_len
    right_x_len = init%right_x_len

    double_edge = init%double_edge
    jiggle_atoms = init%jiggle_atoms
    jiggle_amount = init%jiggle_amount
    jiggle_region = init%jiggle_region
    jiggle_seed = init%jiggle_seed

    min_dist = -1.0
    left_xx_pbv = 0
    left_yy_pbv = 0
    left_zz_pbv = 0
    right_xx_pbv = 0
    right_yy_pbv = 0
    right_zz_pbv = 0
    deletion_method = ""

!! These entries are for _dsc, "do_later"
!! They are written by gb_energy_cc_dsc, so do not contain user typos
    main_name = ""
    pbv_num = -1
    low_min_dist = -1.0d0
    e_per_a = -1.0d0
    dsc_num = -1
    place_index = -1
    dsc_int_vec = -1
    method_num = -1
    dsc_vec(3) = 0.0d0
    dsc_xx_adj = 0.0d0
    xx_boundary_adj = 0.0d0
    one_at_a_time = .false.
    lammps_control_file = "temp.min.in"
    delete_by_energy = .false.
    show_pressure = .false.
    delete_by_pressure = .false.
    deletion_control_file = ""
    lammps_command = "./lmp_serial"
    scr_file = "temp.energy_dsc.scr"
    sub_sh_file = "temp.energy_dsc.sub.sh"


    write(run_out_unit,*) 'reading material_card'
    read(in_unit, material_card)
    if( alat .lt. 0 ) then
       write(run_out_unit,*) 'alat missing or negative'
       call exit(1)
    end if
    init%alat = alat
    init%crystal = crystal

    write(run_out_unit,*) 'reading pbv_card'
    read(in_unit, pbv_card)
    
    init%fixed_x_len = fixed_x_len
    if( fixed_x_len ) then
       init%left_x_len = left_x_len
       init%right_x_len = right_x_len
       if( left_x_len .le. 0.0d0 ) then
          write(0,*) "fixed_x_len, but left_x_len missing or not positive"
          call exit(1)
       end if
       if( right_x_len .le. 0.0d0 ) then
          write(0,*) "fixed_x_len, but right_x_len missing or not positive"
          call exit(1)
       end if
    end if

    init%left_pbv_ii(1, 1:3) = left_xx_pbv
    init%left_pbv_ii(2, 1:3) = left_yy_pbv
    init%left_pbv_ii(3, 1:3) = left_zz_pbv
    init%right_pbv_ii(1, 1:3) = right_xx_pbv
    init%right_pbv_ii(2, 1:3) = right_yy_pbv
    init%right_pbv_ii(3, 1:3) = right_zz_pbv

    write(run_out_unit,*) 'reading control_card'
    read(in_unit, control_card)
    if( .not. dsc_flag ) then
       if( name .eq. "" ) then
          write(run_out_unit,*) "name missing (or blank)"
          call exit(1)
       end if
    end if
    init%name = name

    if( jiggle_atoms ) then
       if( jiggle_amount .le. 0.0d0 ) then
          write(run_out_unit,*) "jiggle_atoms, but jiggle_amount missing or not positive"
          call exit(1)
       end if
       if( jiggle_region .le. 0.0d0 ) then
          write(run_out_unit,*) "jiggle_atoms, but jiggle_region missing or not positive"
          call exit(1)
       end if
       if( jiggle_seed .le. 0 ) then
          write(run_out_unit,*) "jiggle_atoms, jiggle_seed missing or not positive"
          call exit(1)
       end if
    else
       if( jiggle_amount .ne. -1.0d0 ) then
          write(run_out_unit,*) "jiggle_amount given, but not jiggle_atoms"
          call exit(1)
       end if
       if( jiggle_region .ne. -1.0d0 ) then
          write(run_out_unit,*) "jiggle_region given, but not jiggle_atoms"
          call exit(1)
       end if
       if( jiggle_seed .ne. -1 ) then
          write(run_out_unit,*) "jiggle_seed give, but not jiggle_atoms"
          call exit(1)
       end if
    end if

    if( min_dist .eq. -1.0 ) then
       init%min_dist_flag = .false.
       init%min_dist = min_dist
    else
       init%min_dist_flag = .true.
       init%min_dist = min_dist
       if( min_dist .lt. 0 ) then
          write(run_out_unit,*) 'min_dist negative.  Omit entirely for no min_dist.'
          call exit(1)
       end if
    end if

    if( low_min_dist .eq. -1.0d0 ) then
       init%low_min_dist = init%min_dist
    else
       init%low_min_dist = low_min_dist
    end if

    init%crystal_to_delete = crystal_to_delete

    if( deletion_method .eq. "delete_by_averaging" ) then
       init%method = delete_by_averaging
       if( crystal_to_delete .ne. -1 ) then
          write(run_out_unit,*) 'delete_by_averaging, but crystal_to_delete given'
          call exit(1)
       end if
    else if( deletion_method .eq. "fixed_crystal" ) then
       init%method = fixed_crystal
       if( crystal_to_delete .ne. 1 .and. crystal_to_delete .ne. 2 ) then
          write(run_out_unit,*) 'crystal_to_delete must be 1 or 2'
          call exit(1)
       end if
    else if( deletion_method .eq. "single_crystal" ) then
       init%method = single_crystal
    else if( deletion_method .eq. "" ) then
       write(run_out_unit,*) 'deletion_method missing or blank'
       call exit(1)
    else
       write(run_out_unit,*) 'invalid deletion_method: ', deletion_method 
       write(run_out_unit,*) '  sb delete_by_averaging or fixed_crystal'
       call exit(1)
    end if

    init%nbr_delta = nbr_delta

    init%nbrs_by_dist = nbrs_by_dist
    init%two_boundaries = two_boundaries
    init%concentration = concentration
    init%sigma_flag = sigma_flag

    init%edge_size = edge_size

    init%dsc_flag = dsc_flag

    init%main_name = main_name
    init%pbv_num = pbv_num
    init%e_per_a = e_per_a
    init%dsc_num = dsc_num
    init%dsc_int_vec = dsc_int_vec
    init%place_index = place_index
    init%method_num = method_num
    init%dsc_vec = dsc_vec
    init%dsc_xx_adj = dsc_xx_adj
    init%xx_boundary_adj = xx_boundary_adj
    init%one_at_a_time = one_at_a_time
    init%lammps_flag = lammps_flag
    init%output_flag = output_flag
    init%lammps_control_file = lammps_control_file

    if( delete_by_energy .and. delete_by_pressure ) then
       write(run_out_unit,*) "both delete_by_energy and delete_by_pressure"
       call exit(1)
    end if
    if( delete_by_pressure ) then
       show_pressure = .true.
    end if
    if( delete_by_energy .or. delete_by_pressure ) then
       if( deletion_control_file .eq. "" ) then
          write(run_out_unit,*) "read_build_gb_input  delete_by_energy ",          &
                     "or delete_by_pressure,",                            &
                     "but deletion control file missing or ''"
          call exit(1)
       end if
    end if

    init%delete_by_energy = delete_by_energy
    init%show_pressure = show_pressure
    init%delete_by_pressure = delete_by_pressure
    init%deletion_control_file = deletion_control_file
    init%lammps_command = lammps_command
    init%scr_file = scr_file
    init%sub_sh_file = sub_sh_file

    init%double_edge = double_edge
    init%jiggle_atoms = jiggle_atoms
    init%jiggle_amount = jiggle_amount
    init%jiggle_region = jiggle_region
    init%jiggle_seed = jiggle_seed

    return
  end subroutine read_build_gb_input



  subroutine edit_input(init)
    use build_gb_control_mod
    use build_gb_global_mod
    use crystal_mod
    implicit none

    type(build_gb_control_type), intent(in) :: init

    logical :: result

    fixed_x_len = init%fixed_x_len
    if( fixed_x_len ) then
       left_x_len = init%left_x_len
       right_x_len = init%right_x_len
    end if

!! For island
    island_flag = init%island_flag
    if( island_flag ) then
       island_cut_type = init%island_cut_type
       if( island_cut_type .eq. island_circle ) then
          island_radius = init%island_radius
       end if
       island_min_in_plane = init%island_min_in_plane
       island_min_thickness = init%island_min_thickness
    end if


    if( init%lammps_flag .and. (.not. init%dsc_flag) ) then
       write(run_out_unit,*) 'build_gb_input  edit_input  lammps_flag, but not dsc_flag. ??'
       call exit(1)
    end if

    mat%alat = init%alat
    mat%crystal%name = init%crystal
    call check_crystal_name(mat%crystal, result)
    if( .not. result ) then
       write(run_out_unit,*) 'check_crystal_name returned false'
       call exit(1)
    end if

    if( init%xi_chi_alat .eq. -1.0d0 ) then
       xi_chi_alat = init%alat
    else
       xi_chi_alat = init%xi_chi_alat
    end if

!! For bow-out
    fix_yy = init%fix_yy
    fix_yy_width = init%fix_yy_width
    fix_zz = init%fix_zz
    fix_zz_width = init%fix_zz_width


    concentration = init%concentration

    left_pbv_ii = init%left_pbv_ii
    right_pbv_ii = init%right_pbv_ii

    call validate_pbv( left_pbv_ii, right_pbv_ii, init%crystal )
    call build_trans_matrix( left_pbv_ii, left_trans )
    call build_trans_matrix( right_pbv_ii, right_trans )
    

!! debugging
    write(run_out_unit,*) 'left_trans: '
    write(run_out_unit,*) left_trans
    write(run_out_unit,*) 'right_trans: '
    write(run_out_unit,*) right_trans
!! end debugging 

    if( fixed_x_len ) then
       left_len(1) = 0.5d0 * mat%alat * left_x_len
       right_len(1) = 0.5d0 * mat%alat * right_x_len
    end if
    name = init%name
    nbrs_by_dist = init%nbrs_by_dist
    if( nbrs_by_dist ) then
       full_nbr_list = .true.
    else
       full_nbr_list = .false.
    end if
    two_boundaries = init%two_boundaries
    sigma_flag = init%sigma_flag
    min_dist_flag = init%min_dist_flag
    min_dist = init%min_dist
    if( init%full_min_dist .eq. -1.0d0 ) then
       full_min_dist = min_dist
    else
       full_min_dist = init%full_min_dist
    end if
    low_min_dist = init%low_min_dist
    e_per_a = init%e_per_a
    lammps_flag = init%lammps_flag
    output_flag = init%output_flag
    one_at_a_time = init%one_at_a_time
    method = init%method
    crystal_to_delete = init%crystal_to_delete
    delta = init%nbr_delta

    double_edge = init%double_edge
    jiggle_atoms = init%jiggle_atoms
    jiggle_region = init%jiggle_region
    jiggle_amount = init%jiggle_amount
    jiggle_seed = init%jiggle_seed

    return
  end subroutine edit_input


  subroutine validate_pbv( left, right , crystal )
    use build_gb_global_mod
    implicit none

    integer, intent(in)      :: left(3, 3)
    integer, intent(in)      :: right(3, 3)
    character(3), intent(in) :: crystal

    integer :: ii, left_len_sq, right_len_sq
    integer :: left_vec(3), right_vec(3)

!! check non-zero
!! check that they are lattice vectors
!! check that they are orthogonal
!! check that y and z are the same length between left and right

    call validate_one_pbv( left, crystal )
    call validate_one_pbv( right, crystal )

    do ii = 1, 3
       left_vec = left(ii, 1:3)
       left_len_sq = dot_product( left_vec, left_vec )
       left_len(ii) = 0.5d0 * mat%alat * sqrt( dble( left_len_sq ) )
       if( .not. island_flag ) then
          right_vec = right(ii, 1:3)
          right_len_sq = dot_product( right_vec, right_vec )
          right_len(ii) = 0.5d0 * mat%alat * sqrt( dble( right_len_sq ) )
          if( ii .ne. 1 ) then
             if( left_len_sq .ne. right_len_sq ) then
                write(run_out_unit,*) 'y or z lengths not compatible: ii: ', ii
                call exit(1)
             end if
          end if
       end if
    end do

    if( island_flag ) then
       perlb = 0.0d0
       perub = left_len
    else
       perlb(2:3) = 0.0d0
       perub(2:3) = left_len(2:3)     !! left and right are identical here
       perlb(1) = -left_len(1)
       perub(1) = right_len(1)
    end if
    perlen = perub - perlb
    box_center = ( perlb + perub ) / 2.0d0

    return
  end subroutine validate_pbv



  subroutine validate_one_pbv( pbv_ii, crystal )
    use crystal_mod
    implicit none

    integer, intent(in) :: pbv_ii(3, 3)
    character(3), intent(in) :: crystal

    integer :: ii, jj
    integer :: vect_ii(3), vect_jj(3)
    logical :: result
    integer :: len_sq, prod

    do ii = 1, 3
       vect_ii = pbv_ii(ii, 1:3)
       len_sq = dot_product( vect_ii, vect_ii )
       if( len_sq .eq. 0 ) then
          write(run_out_unit,*) 'some pbv vector is missing or zero'
          call exit(1)
       end if
       call crystal_in_lattice( vect_ii, crystal, result )
       if( .not. result ) then
          write(run_out_unit,*) 'some pbv not a lattice vector: pbv: ', vect_ii
          call exit(1)
       end if
       do jj = ii + 1, 3
          vect_jj = pbv_ii(jj, 1:3)
          prod = dot_product( vect_jj, vect_ii )
          if( prod .ne. 0 ) then
             write(run_out_unit,*) 'two pbv not orthogonal'
             write(run_out_unit,*) '  ii, pbv: ', ii, vect_ii
             write(run_out_unit,*) '  jj, pbv: ', jj, vect_jj
             call exit(1)
          end if
       end do
    end do

    return
  end subroutine validate_one_pbv


  subroutine build_trans_matrix( matrix_ii, matrix_real )
    implicit none

    integer, intent(in) :: matrix_ii(3, 3)
    real*8, intent(out) :: matrix_real(3, 3)

    integer :: ii
    integer :: vec_ii(3)
    real*8  :: vec_real(3)
    real*8  :: length

    !! convert to rotation matrix by normalizing each vector

    do ii = 1, 3
       vec_ii = matrix_ii(ii, 1:3)
       length = dot_product( vec_ii, vec_ii )
       length = sqrt( length )
       vec_real = vec_ii / length
       matrix_real(ii, 1:3) = vec_real
    end do

    return
  end subroutine build_trans_matrix

end module build_gb_input_mod
