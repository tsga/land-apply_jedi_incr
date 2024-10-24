 program apply_incr_noahmp_snow

 use netcdf

 use NoahMPdisag_module, only : noahmp_type, UpdateAllLayers

 implicit none

 include 'mpif.h'

 type(noahmp_type)      :: noahmp_state

 integer :: res, len_land_vec 
 character(len=8) :: date_str 
 character(len=2) :: hour_str
 logical :: frac_grid

 ! index to map between tile and vector space 
 integer, allocatable :: tile2vector(:,:) 
 double precision, allocatable :: increment(:) 
 double precision, allocatable :: swe_back(:) 
 double precision, allocatable :: snow_depth_back(:) 

 integer :: ierr, irank, nprocs, myrank, lunit, ncid, n
 integer :: ntiles, ens_size, ens_mem, tile_num
 character(len=3) :: ens_str
 logical :: file_exists

 ! restart variables that apply to full grid cell 
 ! (cf those that are land only)
 type grid_type
     double precision, allocatable :: land_frac          (:)
     double precision, allocatable :: swe                (:)
     double precision, allocatable :: snow_depth         (:)
     character(len=10)  :: name_snow_depth
     character(len=10)  :: name_swe
 endtype 
 type(grid_type) :: grid_state

 character(len=512) :: orog_path, rst_path_full, inc_path_full
 character(len=256) :: rst_path, inc_path
 character*20       :: otype ! orography filename stub. For atm only, oro_C${RES}, for atm/ocean oro_C${RES}.mx100

 character(len=512) :: restart_file
 character(len=1)   :: tilech
 character(len=512) :: ioerrmsg

 namelist /noahmp_snow/ date_str, hour_str, res, frac_grid, rst_path, inc_path, orog_path, otype, ntiles, ens_size
!
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    print*
    print*,"starting apply_incr_noahmp_snow program on rank ", myrank, ' of ', nprocs, ' procs'

    ! SET NAMELIST DEFAULTS
    frac_grid = .false.
    rst_path = './'
    inc_path = './'
    ntiles = 6
    ens_size = 1

    ! READ NAMELIST 

    inquire (file='apply_incr_nml', exist=file_exists) 

    if (.not. file_exists) then
        write (6, *) 'ERROR: apply_incr_nml does not exist'
        call mpi_abort(mpi_comm_world, 10)
    end if

    open (action='read', file='apply_incr_nml', iostat=ierr, newunit=lunit, iomsg=ioerrmsg)
    read (nml=noahmp_snow, iostat=ierr, unit=lunit)
    close (lunit)
    if (ierr /= 0) then
        print*, "Error code from namelist read", ierr
        write(6,*) trim(ioerrmsg)         
        call mpi_abort(mpi_comm_world, 10)
    end if
    ! uncommenting this may help catch a namelist error
    ! if (myrank==0) then
    !    write (6, noahmp_snow)
    ! !    print*, 'ens_size ', ens_size, ' ntiles ', ntiles
    ! end if

    ! SET VARIABLE NAMES FOR SNOW OVER LAND AND GRID
    if (frac_grid) then 
        noahmp_state%name_snow_depth =  'snodl     '
        noahmp_state%name_swe =         'weasdl    '
        grid_state%name_snow_depth =    'snwdph    '
        grid_state%name_swe =           'sheleg    '
    else
        noahmp_state%name_snow_depth =  'snwdph    '
        noahmp_state%name_swe =         'sheleg    '
        grid_state%name_snow_depth =    'snwdph    '
        grid_state%name_swe =           'sheleg    '
    endif

    do irank=myrank, ntiles*ens_size - 1, nprocs
        ens_mem = irank/ntiles + 1            !ensemble member
        tile_num = MOD(irank, ntiles) + 1      !tile number

        write(ens_str, '(I3.3)') ens_mem

!TBCL: keep the default for ens_size=1
        if(ens_size > 1) then 
            rst_path_full = trim(rst_path)//"/mem"//ens_str//"/"
            inc_path_full = trim(inc_path)//"/mem"//ens_str//"/"
        else
            rst_path_full = trim(rst_path)      !//"/mem000/"
            inc_path_full = trim(inc_path)      !//"/mem000/"
        endif

        print*
        print*, "Proc ", myrank, " ensemble member ", ens_mem, " tile ", tile_num

        ! GET MAPPING INDEX (see subroutine comments re: source of land/sea mask)

        call get_fv3_mapping(myrank, ens_mem, tile_num, rst_path_full, date_str, hour_str, res, len_land_vec, frac_grid, tile2vector)
    
        ! SET-UP THE NOAH-MP STATE  AND INCREMENT
        
        ! The allocations are inside the loop because different ensemble members could have different len_land_vec
        allocate(noahmp_state%swe                (len_land_vec)) ! values over land only
        allocate(noahmp_state%snow_depth         (len_land_vec)) ! values over land only 
        allocate(noahmp_state%active_snow_layers (len_land_vec)) 
        allocate(noahmp_state%swe_previous       (len_land_vec))
        allocate(noahmp_state%snow_soil_interface(len_land_vec,7))
        allocate(noahmp_state%temperature_snow   (len_land_vec,3))
        allocate(noahmp_state%snow_ice_layer     (len_land_vec,3))
        allocate(noahmp_state%snow_liq_layer     (len_land_vec,3))
        allocate(noahmp_state%temperature_soil   (len_land_vec))
        allocate(increment   (len_land_vec)) ! increment to snow depth over land

        if (frac_grid) then
            allocate(grid_state%land_frac          (len_land_vec)) 
            allocate(grid_state%swe                (len_land_vec)) ! values over full grid
            allocate(grid_state%snow_depth         (len_land_vec)) ! values over full grid
            allocate(swe_back                      (len_land_vec)) ! save background 
            allocate(snow_depth_back               (len_land_vec)) !
        endif

        ! READ RESTART FILE 
        write(tilech, '(i1.1)') (tile_num)
        restart_file = trim(rst_path_full)//"/"//date_str//"."//hour_str//"0000.sfc_data.tile"//tilech//".nc"

        call   read_fv3_restart(trim(restart_file), res, ncid, &          !tile_num, rst_path_full, date_str, hour_str, 
                    len_land_vec, tile2vector, frac_grid, noahmp_state, grid_state)

        ! READ SNOW DEPTH INCREMENT

        call   read_fv3_increment(tile_num, inc_path_full, date_str, hour_str, res, &
                    len_land_vec, tile2vector, noahmp_state%name_snow_depth, increment)
    
        if (frac_grid) then ! save background
            swe_back = noahmp_state%swe
            snow_depth_back = noahmp_state%snow_depth
        endif 

        ! ADJUST THE SNOW STATES OVER LAND

!TBCL: return and check error code from this call (for now assume it is well handled inside function)
        call UpdateAllLayers(len_land_vec, increment, noahmp_state)

        ! IF FRAC GRID, ADJUST SNOW STATES OVER GRID CELL

        if (frac_grid) then

            ! get the land frac 
            call  read_fv3_orog(tile_num, res, orog_path, otype, len_land_vec, tile2vector, & 
                    grid_state)

            do n=1,len_land_vec 
                    grid_state%swe(n) = grid_state%swe(n) + & 
                                    grid_state%land_frac(n)* ( noahmp_state%swe(n) - swe_back(n)) 
                    grid_state%snow_depth(n) = grid_state%snow_depth(n) + & 
                                    grid_state%land_frac(n)* ( noahmp_state%snow_depth(n) - snow_depth_back(n)) 
            enddo

        endif

        ! WRITE OUT ADJUSTED RESTART

        call   write_fv3_restart(trim(restart_file), noahmp_state, grid_state, res, ncid, len_land_vec, & 
                    frac_grid, tile2vector) 

        ! CLOSE RESTART FILE 
        ! print*
        ! print*,"apply_incr_noahmp_snow, closing restart on proc ", myrank, " ensemble member ", ens_mem, " tile ", tile_num
        ierr = nf90_close(ncid)
        call netcdf_err( ierr, "closing restart file "//trim(restart_file) )
        
        ! Deallocate. These are required incase a single process loops through multiple tiles with different mapping     
        if (allocated(tile2vector)) deallocate(tile2vector)   
        deallocate(noahmp_state%swe) ! values over land only
        deallocate(noahmp_state%snow_depth) ! values over land only 
        deallocate(noahmp_state%active_snow_layers) 
        deallocate(noahmp_state%swe_previous)
        deallocate(noahmp_state%snow_soil_interface)
        deallocate(noahmp_state%temperature_snow)
        deallocate(noahmp_state%snow_ice_layer)
        deallocate(noahmp_state%snow_liq_layer)
        deallocate(noahmp_state%temperature_soil)
        deallocate(increment) ! increment to snow depth over land

        if (frac_grid) then
            deallocate(grid_state%land_frac) 
            deallocate(grid_state%swe) ! values over full grid
            deallocate(grid_state%snow_depth) ! values over full grid
            deallocate(swe_back) ! save background 
            deallocate(snow_depth_back) !
        endif

        ! print*
        ! print*, "finisheed loop proc ", myrank, " ensemble member ", ens_mem, " tile ", tile_num

    enddo

    print*, "Finisheed on proc ", myrank
    call mpi_finalize(ierr)

 contains 

!--------------------------------------------------------------
! if at netcdf call returns an error, print out a message
! and stop processing.
!--------------------------------------------------------------
 subroutine netcdf_err( err, string )

        implicit none

        include 'mpif.h'

        integer, intent(in) :: err
        character(len=*), intent(in) :: string
        character(len=80) :: errmsg

        if( err == nf90_noerr )return
        errmsg = nf90_strerror(err)
        print*,''
        print*,'fatal error: ', trim(string), ': ', trim(errmsg)
        print*,'stop.'
        call mpi_abort(mpi_comm_world, 999)

        return
 end subroutine netcdf_err


!--------------------------------------------------------------
! Get land sea mask from fv3 restart, and use to create 
! index for mapping from tiles (FV3 UFS restart) to vector
!  of land locations (offline Noah-MP restart)
! NOTE: slmsk in the restarts counts grid cells as land if 
!       they have a non-zero land fraction. Excludes grid 
!       cells that are surrounded by sea (islands). The slmsk 
!       in the oro_grid files (used by JEDI for screening out 
!       obs is different, and counts grid cells as land if they 
!       are more than 50% land (same exclusion of islands). If 
!       we want to change these definitations, may need to use 
!       land_frac field from the oro_grid files.
!--------------------------------------------------------------

 subroutine get_fv3_mapping(myrank, ens_mem, tile_num, rst_path, date_str, hour_str, res, & 
                len_land_vec, frac_grid, tile2vector)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: myrank, ens_mem, tile_num, res
 character(len=*), intent(in) :: rst_path
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 logical, intent(in) :: frac_grid
 integer, allocatable, intent(out) :: tile2vector(:,:)
 integer :: len_land_vec

 character(len=512) :: restart_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr,  ncid
 integer :: id_dim, id_var, fres
 integer :: slmsk(res,res) ! saved as double in the file, but i think this is OK
 integer :: vtype(res,res) ! saved as double in the file, but i think this is OK
 integer, parameter :: vtype_landice=15
 double precision :: fice(res,res)
 double precision, parameter :: fice_fhold = 0.00001
 integer :: i, j, nn

    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    restart_file = trim(rst_path)//"/"//date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ! write (6, *) 'calculate mapping from land mask in ', trim(restart_file)

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ! READ MASK 
    ierr=nf90_inq_varid(ncid, "slmsk", id_var)
    call netcdf_err(ierr, 'reading slmsk id' )
    ierr=nf90_get_var(ncid, id_var, slmsk)
    call netcdf_err(ierr, 'reading slmsk' )
 
    ! REMOVE GLACIER GRID POINTS
    ierr=nf90_inq_varid(ncid, "vtype", id_var)
    call netcdf_err(ierr, 'reading vtype id' )
    ierr=nf90_get_var(ncid, id_var, vtype)
    call netcdf_err(ierr, 'reading vtype' )

    ! remove land grid cells if glacier land type
    do i = 1, res 
        do j = 1, res  
            if ( vtype(i,j) ==  vtype_landice)  slmsk(i,j)=0 ! vtype is integer, but stored as double
        enddo 
    enddo
 
    if (frac_grid) then 

        write (6, *) 'fractional grid: ammending mask to exclude sea ice from', trim(restart_file)

        ierr=nf90_inq_varid(ncid, "fice", id_var)
        call netcdf_err(ierr, 'reading fice id' )
        ierr=nf90_get_var(ncid, id_var, fice)
        call netcdf_err(ierr, 'reading fice' )

        ! remove land grid cells if ice is present
        do i = 1, res 
            do j = 1, res  
                if (fice(i,j) > fice_fhold ) slmsk(i,j)=0
            enddo 
        enddo


    endif
 
    ! get number of land points
    len_land_vec = 0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk(i,j) == 1)  len_land_vec = len_land_vec+ 1  
        enddo 
    enddo
    
    write(6,*) 'Number of land points on proc ', myrank, 'ensmem ', ens_mem, ' tilenum ', tile_num, ' :',  len_land_vec

    allocate(tile2vector(len_land_vec,2)) 

    nn=0
    do i = 1, res 
        do j = 1, res 
             if ( slmsk(i,j) == 1)   then 
                nn=nn+1
                tile2vector(nn,1) = i 
                tile2vector(nn,2) = j 
             endif
        enddo 
    enddo

end subroutine get_fv3_mapping


!--------------------------------------------------------------
! open fv3 restart, and read in required variables
! file is opened as read/write and remains open
!--------------------------------------------------------------
 subroutine read_fv3_restart(restart_file, res, ncid, &                      !tile_num, rst_path, date_str, hour_str, 
                len_land_vec,tile2vector, frac_grid, noahmp_state, grid_state)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: res, len_land_vec   !tile_num, 
 character(len=*), intent(in) :: restart_file  !rst_path
!  character(len=8), intent(in) :: date_str 
!  character(len=2), intent(in) :: hour_str 
 integer, intent(in) :: tile2vector(len_land_vec,2)
 logical, intent(in) :: frac_grid

 integer, intent(out) :: ncid
 type(noahmp_type), intent(inout)  :: noahmp_state
 type(grid_type), intent(inout)  :: grid_state

!  character(len=512) :: restart_file
!  character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr, id_dim, fres
 integer :: nn

    ! OPEN FILE

    inquire(file=trim(restart_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'restart_file does not exist, ', &
                    trim(restart_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ! write (6, *) 'opening ', trim(restart_file)

    ierr=nf90_open(trim(restart_file),nf90_write,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(restart_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
    call netcdf_err(ierr, 'reading xaxis_1 in '//trim(restart_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading xaxis_1 in '//trim(restart_file) )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong in '//trim(restart_file)
       call mpi_abort(mpi_comm_world, ierr)
    endif

   ! read swe over land (file name: sheleg, vert dim 1) 
    ! this call has file name added for inspection. 
    !In case of failure, filename has all info about proc rank, ens member, and tile number
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                        noahmp_state%name_swe, noahmp_state%swe)

    ! read snow_depth over land (file name: snwdph, vert dim 1)
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                        noahmp_state%name_snow_depth, noahmp_state%snow_depth)

   if (frac_grid) then 
       ! read swe over grid cell  (file name: sheleg, vert dim 1) 
        call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                            grid_state%name_swe, grid_state%swe)

        ! read snow_depth  over grid cell (file name: snwdph, vert dim 1)
        call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                            grid_state%name_snow_depth, grid_state%snow_depth)
    endif

    ! read active_snow_layers (file name: snowxy, vert dim: 1) 
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                        'snowxy    ', noahmp_state%active_snow_layers)

    ! read swe_previous (file name: sneqvoxy, vert dim: 1) 
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 0, & 
                        'sneqvoxy  ', noahmp_state%swe_previous)

    ! read snow_soil_interface (file name: zsnsoxy, vert dim: 7) 
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, 7,  tile2vector, & 
                        'zsnsoxy   ', noahmp_state%snow_soil_interface)

    ! read temperature_snow (file name: tsnoxy, vert dim: 3) 
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, 3, tile2vector, & 
                        'tsnoxy    ', noahmp_state%temperature_snow)

    ! read snow_ice_layer (file name:  snicexy, vert dim: 3) 
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, 3, tile2vector, & 
                        'snicexy    ', noahmp_state%snow_ice_layer)

    ! read snow_liq_layer (file name: snliqxy, vert dim: 3) 
    call read_nc_var3D(ncid, trim(restart_file), len_land_vec, res, 3, tile2vector, & 
                        'snliqxy   ', noahmp_state%snow_liq_layer)

    ! read temperature_soil (file name: stc, use layer 1 only, vert dim: 1) 
    call read_nc_var2D(ncid, trim(restart_file), len_land_vec, res, tile2vector, 4, & 
                        'stc       ', noahmp_state%temperature_soil)

end subroutine read_fv3_restart


!--------------------------------------------------------------
! open fv3 orography file, and read in land fraction
!--------------------------------------------------------------
 subroutine read_fv3_orog(tile_num, res, orog_path, otype, len_land_vec, tile2vector, & 
                grid_state)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: tile_num, res, len_land_vec
 character(len=*), intent(in)    :: orog_path
 character(len=20), intent(in)   :: otype
 integer, intent(in) :: tile2vector(len_land_vec,2)
 type(grid_type), intent(inout) :: grid_state

 character(len=250) :: filename
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ncid, id_dim, id_var, ierr, fres

    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    filename =trim(orog_path)//"/"//trim(otype)//".tile"//rankch//".nc"

    inquire(file=trim(filename), exist=file_exists)

    if (.not. file_exists) then
            print *, 'filename does not exist, ', &
                    trim(filename) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ! write (6, *) 'opening ', trim(filename)

    ierr=nf90_open(trim(filename),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(filename) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'lon', id_dim)
    call netcdf_err(ierr, 'reading lon id from '//trim(filename) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading lon from '//trim(filename) )

    if ( fres /= res) then
       print*,'fatal error: dimensions wrong in file '//trim(filename)
       call mpi_abort(mpi_comm_world, ierr)
    endif

   ! read swe over grid cell  (file name: sheleg, vert dim 1) 
    call read_nc_var2D(ncid, trim(filename), len_land_vec, res, tile2vector, 0, & 
                        'land_frac  ', grid_state%land_frac)

    ! close file 
    ! write (6, *) 'closing ', trim(filename)

    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(filename) )

end subroutine read_fv3_orog

!--------------------------------------------------------------
!  read in snow depth increment from jedi increment file
!  file format is same as restart file
!--------------------------------------------------------------
 subroutine read_fv3_increment(tile_num, inc_path, date_str, hour_str, res, & 
                len_land_vec,tile2vector, control_var, increment)

 implicit none 

 include 'mpif.h'

 integer, intent(in) :: tile_num, res, len_land_vec
 character(len=*), intent(in) :: inc_path
 character(len=8), intent(in) :: date_str 
 character(len=2), intent(in) :: hour_str 
 integer, intent(in) :: tile2vector(len_land_vec,2)
 character(len=10), intent(in)  :: control_var
 double precision, intent(out) :: increment(len_land_vec)     ! snow depth increment

 character(len=512) :: incr_file
 character(len=1) :: rankch
 logical :: file_exists
 integer :: ierr 
 integer :: id_dim, id_var, fres, ncid
 integer :: nn


 
    ! OPEN FILE
    write(rankch, '(i1.1)') (tile_num)
    incr_file = trim(inc_path)//"/"//"snowinc."//date_str//"."//hour_str//"0000.sfc_data.tile"//rankch//".nc"

    inquire(file=trim(incr_file), exist=file_exists)

    if (.not. file_exists) then
            print *, 'incr_file does not exist, ', &
                    trim(incr_file) , ' exiting'
            call mpi_abort(mpi_comm_world, 10) 
    endif

    ! write (6, *) 'opening ', trim(incr_file)

    ierr=nf90_open(trim(incr_file),nf90_nowrite,ncid)
    call netcdf_err(ierr, 'opening file: '//trim(incr_file) )

    ! CHECK DIMENSIONS
    ierr=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
    call netcdf_err(ierr, 'reading xaxis_1 '//trim(incr_file) )
    ierr=nf90_inquire_dimension(ncid,id_dim,len=fres)
    call netcdf_err(ierr, 'reading xaxis_1 '//trim(incr_file) )

    if ( fres /= res) then
       print*,'fatal error: dimension fres ',fres, ' in '//trim(incr_file), ' not equal to res ',res
       call mpi_abort(mpi_comm_world, ierr)
    endif

    ! read snow_depth (file name: snwdph, vert dim 1)
    call read_nc_var2D(ncid, trim(incr_file), len_land_vec, res, tile2vector, 0, & 
                        control_var, increment)

    ! close file 
    ! write (6, *) 'closing ', trim(incr_file)

    ierr=nf90_close(ncid)
    call netcdf_err(ierr, 'closing file: '//trim(incr_file) )

end subroutine  read_fv3_increment

!--------------------------------------------------------
! Subroutine to read in a 2D variable from netcdf file, 
! and save to noahmp vector
!--------------------------------------------------------

subroutine read_nc_var2D(ncid, file_name, len_land_vec, res, tile2vector, in3D_vdim,  & 
                         var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res 
    ! file name added for inspection. In case of failure, filename has all info about proc rank, ens member, and tile number
    character(len=*), intent(in)    :: file_name  
    character(len=10), intent(in)   :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    integer, intent(in)             :: in3D_vdim ! 0 - input is 2D, 
                                                 ! >0, gives dim of 3rd dimension
    double precision, intent(out)   :: data_vec(len_land_vec) 

    double precision :: dummy2D(res, res) 
    double precision :: dummy3D(res, res, in3D_vdim)  
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//var_name//' id in '//trim(file_name) )
    if (in3D_vdim==0) then
        ierr=nf90_get_var(ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'reading '//var_name//' data in '//trim(file_name) )
    else  ! special case for reading in 3D variable, and retaining only 
          ! level 1
        ierr=nf90_get_var(ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'reading '//var_name//' data in '//trim(file_name) )
        dummy2D=dummy3D(:,:,1) 
    endif

    do nn=1,len_land_vec 
        data_vec(nn) = dummy2D(tile2vector(nn,1), tile2vector(nn,2))
    enddo

end subroutine read_nc_var2D

!--------------------------------------------------------
! Subroutine to read in a 3D variable from netcdf file, 
! and save to noahmp vector
!--------------------------------------------------------

subroutine read_nc_var3D(ncid, file_name, len_land_vec, res, vdim,  & 
                tile2vector, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res, vdim
    character(len=*), intent(in)    :: file_name
    character(len=10), intent(in)   :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    double precision, intent(out)   :: data_vec(len_land_vec, vdim)

    double precision :: dummy3D(res, res, vdim) 
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//var_name//' id in '//trim(file_name) )
    ierr=nf90_get_var(ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'reading '//var_name//' data in '//trim(file_name) )

    do nn=1,len_land_vec 
        data_vec(nn,:) = dummy3D(tile2vector(nn,1), tile2vector(nn,2), :) 
    enddo

end subroutine read_nc_var3D

!--------------------------------------------------------------
! write updated fields tofv3_restarts  open on ncid
!--------------------------------------------------------------
 subroutine write_fv3_restart(file_name, noahmp_state, grid_state, res, ncid, len_land_vec, &
                 frac_grid, tile2vector) 

 implicit none 

 integer, intent(in) :: ncid, res, len_land_vec
 character(len=*), intent(in)    :: file_name
 type(noahmp_type), intent(in) :: noahmp_state
 type(grid_type), intent(in) :: grid_state
 logical, intent(in) :: frac_grid
 integer, intent(in) :: tile2vector(len_land_vec,2)

 
   ! write swe over land (file name: sheleg, vert dim 1) 
    call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                        noahmp_state%name_swe, noahmp_state%swe)

    ! write snow_depth over land (file name: snwdph, vert dim 1)
    call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                        noahmp_state%name_snow_depth, noahmp_state%snow_depth)

    if (frac_grid) then
       ! write swe over grid (file name: sheleg, vert dim 1) 
        call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                            grid_state%name_swe, grid_state%swe)

        ! write snow_depth over grid (file name: snwdph, vert dim 1)
        call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                            grid_state%name_snow_depth, grid_state%snow_depth)
    endif 

    ! write active_snow_layers (file name: snowxy, vert dim: 1) 
    call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                        'snowxy    ', noahmp_state%active_snow_layers)

    ! write swe_previous (file name: sneqvoxy, vert dim: 1) 
    call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 0, & 
                        'sneqvoxy  ', noahmp_state%swe_previous)

    ! write snow_soil_interface (file name: zsnsoxy, vert dim: 7) 
    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, 7,  tile2vector, & 
                        'zsnsoxy   ', noahmp_state%snow_soil_interface)

    ! write temperature_snow (file name: tsnoxy, vert dim: 3) 
    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, 3, tile2vector, & 
                        'tsnoxy    ', noahmp_state%temperature_snow)

    ! write snow_ice_layer (file name:  snicexy, vert dim: 3) 
    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, 3, tile2vector, & 
                        'snicexy    ', noahmp_state%snow_ice_layer)

    ! write snow_liq_layer (file name: snliqxy, vert dim: 3) 
    call write_nc_var3D(ncid, trim(file_name), len_land_vec, res, 3, tile2vector, & 
                        'snliqxy   ', noahmp_state%snow_liq_layer)

    ! write temperature_soil (file name: stc, use layer 1 only, vert dim: 1) 
    call write_nc_var2D(ncid, trim(file_name), len_land_vec, res, tile2vector, 4, & 
                        'stc       ', noahmp_state%temperature_soil)


 end subroutine write_fv3_restart


!--------------------------------------------------------
! Subroutine to write a 2D variable to the netcdf file 
!--------------------------------------------------------

subroutine write_nc_var2D(ncid, file_name, len_land_vec, res, tile2vector,   & 
                in3D_vdim, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res
    character(len=*), intent(in)    :: file_name
    character(len=10), intent(in)   :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    integer, intent(in)             :: in3D_vdim ! 0 - input is 2D, 
                                                 ! >0, gives dim of 3rd dimension
    double precision, intent(in)    :: data_vec(len_land_vec)

    double precision :: dummy2D(res, res) 
    double precision :: dummy3D(res, res, in3D_vdim)
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    if (in3D_vdim==0) then 
        ierr=nf90_get_var(ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
    else  ! special case for reading in multi-level variable, and 
          ! retaining only first level.
        ierr=nf90_get_var(ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
        dummy2D = dummy3D(:,:,1)
    endif
    
    ! sub in updated locations (retain previous fields for non-land)  
    do nn=1,len_land_vec 
        dummy2D(tile2vector(nn,1), tile2vector(nn,2)) = data_vec(nn) 
    enddo

    ! overwrite
    if (in3D_vdim==0) then 
        ierr = nf90_put_var( ncid, id_var, dummy2D)
        call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    else 
        dummy3D(:,:,1) = dummy2D 
        ierr = nf90_put_var( ncid, id_var, dummy3D)
        call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    endif
    call remove_checksum(ncid, id_var)
 
end subroutine write_nc_var2D

!--------------------------------------------------------
! Subroutine to write a 3D variable to the netcdf file 
!--------------------------------------------------------

subroutine write_nc_var3D(ncid, file_name, len_land_vec, res, vdim, & 
                tile2vector, var_name, data_vec)

    integer, intent(in)             :: ncid, len_land_vec, res, vdim
    character(len=*), intent(in)    :: file_name
    character(len=10), intent(in)   :: var_name
    integer, intent(in)             :: tile2vector(len_land_vec,2)
    double precision, intent(in)    :: data_vec(len_land_vec, vdim)

    double precision :: dummy3D(res, res, vdim)
    integer          :: nn, ierr, id_var

    ierr=nf90_inq_varid(ncid, trim(var_name), id_var)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' id in '//trim(file_name) )
    ierr=nf90_get_var(ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'reading '//trim(var_name)//' data in '//trim(file_name) )
    
    ! sub in updated locations (retain previous fields for non-land)  
    do nn=1,len_land_vec 
        dummy3D(tile2vector(nn,1), tile2vector(nn,2),:) = data_vec(nn,:)
    enddo

    ! overwrite
    ierr = nf90_put_var( ncid, id_var, dummy3D)
    call netcdf_err(ierr, 'writing '//trim(var_name)//' in '//trim(file_name) )
    call remove_checksum(ncid, id_var)
 
end subroutine write_nc_var3D

!> Remove the checksum attribute from a netcdf record.
!!
!! @param[in] ncid netcdf file id
!! @param[in] id_var netcdf variable id.
!!
!! @author George Gayno NCEP/EMC
 subroutine remove_checksum(ncid, id_var)

 implicit none

 integer, intent(in)       :: ncid, id_var

 integer                   :: error

 error=nf90_inquire_attribute(ncid, id_var, 'checksum')

 if (error == 0) then ! attribute was found

   error = nf90_redef(ncid)
   call netcdf_err(error, 'entering define mode' )

   error=nf90_del_att(ncid, id_var, 'checksum')
   call netcdf_err(error, 'deleting checksum' )

   error= nf90_enddef(ncid)
   call netcdf_err(error, 'ending define mode' )

 endif

 end subroutine remove_checksum

 end program apply_incr_noahmp_snow
