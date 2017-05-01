! ============================================================================
! aux_module.f90
!
!  Module handling aux array fields that are based on provided files.
!  This is based on topography files that are used in GeoClaw and was
!  generalized by David George for work with D-Claw originally.
!
! ============================================================================
!      Copyright (C) 2012-04-05 David George <dgeorge@uw.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module aux_module

    implicit none

    integer, parameter :: NUM_DIM = 2

    ! Container for aux file information
    type aux_file_type
        ! File information
        character(len=150) :: path
        integer :: file_type, init_type
        real(kind=8) :: no_data_value

        ! Location of data
        real(kind=8) :: lower(NUM_DIM), upper(NUM_DIM), dx(NUM_DIM)

        ! Size of data
        integer :: num_cells(NUM_DIM)

        ! Control of use of data
        integer :: min_level, max_level

        ! data
        real(kind=8), allocatable :: data(NUM_DIM)

    end type aux_file_type

    logical, private :: moduel_setup = .false.

    ! Aux Files
    integer :: num_aux_files
    type(aux_file_type), allocatable :: aux_files(:)

contains

    ! ========================================================================
    ! Read auxinit files as specified in auxinit.data
    !
    ! Each file has a type stored in file_type(i).
    !   file_type = 1:  standard GIS format: 3 columns: x,y,z(m)
    !   file_type = 2:  Header as in DEM file, height(m) one value per line
    !   file_type = 3:  Header as in DEM file, height(m) one row per line
    !   file_type = 4:  NetCDf
    ! For other formats modify readauxinit routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Associated with each file is a initialization type, init_type:
    !     as follows:
    !     defines a perturbation (or definition of) aux(i,j,iauxinit)
    ! ========================================================================
    subroutine set_aux_files(fname)

        use geoclaw_module

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: fname

        ! Local storage
        integer, parameter :: unit = 7



      ! Locals
      ! integer, parameter :: iunit = 7
      ! integer :: i,j,iauxinitfile
      ! character*25 :: file_name
      ! logical :: found_file


        ! Open and begin parameter file output
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET_AUX_FILES:'
        write(GEO_PARM_UNIT,*) '---------'

        ! Open data file
        if (present(fname)) then
            call opendatafile(unit, fname)
        else
            call opendatafile(unit, 'aux_files.data')
        endif

        read(unit, *) num_aux_files

        if (num_aux_files == 0) then
            write(GEO_PARM_UNIT,*) '   num_aux_files = 0'
            return
        endif

        write(GEO_PARM_UNIT,*) '   num_aux_files = ', num_aux_files

        ! Read and allocate data parameters for each file
        allocate(aux_files(num_aux_files))

        do i = 1, num_aux_files
            read(unit, *) aux_files(i)%path
            read(unit, *) aux_files(i)%file_type, aux_files(i)%init_type
            read(unit, *) aux_files(i)%min_level, aux_files(i)%max_level

            write(GEO_PARM_UNIT,*) '   '
            write(GEO_PARM_UNIT,*) '   ', aux_files(i)%path
            write(GEO_PARM_UNIT,*) '  file_type = ', aux_files(i)%file_type
            write(GEO_PARM_UNIT,*) '  init_type = ', aux_files(i)%init_type
            write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                  aux_files(i)%min_level, aux_files(i)%max_level

            ! Read in the size of the data here
            call read_aux_file_header(aux_files(i))

            allocate(aux_files(i)%data(aux_files(i)%num_cells(1), &
                                       aux_files(i)%num_cells(2))

            ! Read in the data itself
            call read_aux_file_data(aux_files(i))
        end do
    end subroutine set_auxinit

    ! ========================================================================
    !  Read aux file header
    ! ========================================================================
    subroutine read_aux_file_header(aux_file)

#ifdef NETCDF
        use netcdf
#endif
        use utility_module, only: parse_values

        implicit none

        ! Input/Output
        type(aux_file_type), intent(in out) :: aux_file

        ! Locals
        integer, parameter :: unit = 8
        integer :: mx, my, total_size, status, n
        real(kind=8) :: xll, yll, xhi, yhi, x, y, z
        logical :: found_file
        character(len=80) :: str
        real(kind=8) :: values(10)

        logical, parameter :: verbose = .false.

        inquire(file=aux_file%path, exist=found_file)
        if (.not. found_file) then
            print *, 'Missing aux file:'
            print *, '   ', aux_file%path
            stop
        endif

        print *, ' '
        print *, 'Reading aux file  ', aux_file%path

        select case(abs(aux_file%file_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                open(unit=unit, file=aux_file%path, status='unknown', &
                     form='formatted')

                ! Initial size variables
                total_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(unit,*) xll, yhi
                total_size = total_size + 1
                mx = mx + 1

                ! Go through first row figuring out num_cells(0)
                y = yhi
                do while (yhi == y)
                    read(unit,*) x, y, z
                    total_size = total_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(unit, fmt=*, iostat=status) x, y, z
                    total_size = total_size + 1
                enddo
                if (status > 0) then
                    print *, "IO error occured in ", aux_file%path, ", aborting!"
                    stop
                endif

                ! Calculate remaining values
                my = total_size / mx
                xhi = x
                yll = y
                dx = (xhi - xll) / (mx - 1)
                dy = (yhi - yll) / (my - 1)

                ! Update the aux_file object
                aux_file%num_cells = (mx, my)
                aux_file%lower = (xll, yll)
                aux_file%upper = (xhi, yhi)
                aux_file%dx = (dx, dy)

                close(unit)

            ! ASCII file with header followed by z data
            case(2:3)
                open(unit=unit, file=aux_file%path, status='unknown', &
                     form='formatted')

                read(unit, '(a)') str
                call parse_values(str, n, values)
                aux_file%num_cells(1) = nint(values(1))

                read(unit, '(a)') str
                call parse_values(str, n, values)
                aux_file%num_cells(2) = nint(values(1))

                read(unit, '(a)') str
                call parse_values(str, n, values)
                aux_file%lower(1) = values(1)

                read(unit, '(a)') str
                call parse_values(str, n, values)
                aux_file%lower(2) = values(1)

                read(unit, '(a)') str
                call parse_values(str, n, values)
                aux_file%dx(1) = values(1)
                if (n == 2) then
                    aux_file%dx(2) = values(2)
                else
                    aux_file%dx(2) = aux_file%dx(1)
                endif

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                aux_file%no_data_value = values(1)

                aux_file%upper(:) = aux_file%lower(:) +         &
                                    (aux_file%num_cells(:) - 1) * aux_file%dx(:)

                close(unit)

            case(4)
#ifdef NETCDF

                ! Open file    
                call check_netcdf_error(nf90_open(aux_file%path, nf90_nowrite, nc_file))

                ! Get dimensions - Assume the lon-lat are dimensions 1 and 2
                call check_netcdf_error(nf90_inquire_dimension(nc_file, 1, x_dim_name, mx))
                call check_netcdf_error(nf90_inquire_dimension(nc_file, 2, y_dim_name, my))

                if (verbose) then
                    print *, "Names = (", x_dim_name, ", ", y_dim_name, ")"
                    print *, "Size = (", mx, ", ", my, ")"
                end if

                ! Read in variables
                call check_netcdf_error(nf90_inq_varids(nc_file, num_vars, var_ids))
                
                do n=1, num_vars
                    call check_netcdf_error(nf90_inquire_variable(nc_file, n, var_name, var_type, num_dims, dim_ids))
                    if (verbose) then
                        print *, n, ": ", var_name, "|", var_type, "|", num_dims, "|", dim_ids
                    end if

                    ! Assume dim names are same as var ids
                    if (var_name == x_dim_name) then
                        x_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, x_var_name, x_var_id))
                        ! x_var_id = n
                        ! x_var_name = var_name
                    else if (var_name == y_dim_name) then
                        y_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, y_var_name, y_var_id))
                        ! y_var_id = n
                        ! y_var_name = var_name
                    ! Assume that this is the aux file data
                    else if (num_dims == 2) then
                        z_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, z_var_name, z_var_id))
                        ! z_var_id = n
                        ! z_var_name = var_name
                    else
                        if (verbose) then
                            print *, "Not using var_id ", n
                        end if
                    end if

                end do

                if (verbose) then
                    print *, "x_var_name, x_var_id = ", x_var_name, x_var_id
                    print *, "x_var_name, x_var_id = ", y_var_name, y_var_id
                    print *, "x_var_name, x_var_id = ", z_var_name, z_var_id
                end if

                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xll, start=(/ 1 /)))
                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xhi, start=(/ mx /)))

                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, yll, start=(/ 1 /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, yhi, start=(/ my /)))

                call check_netcdf_error(nf90_close(nc_file))
                
                ! Update the aux_file object
                aux_file%num_cells = (mx, my)
                aux_file%lower = (xll, yll)
                aux_file%upper = (xhi, yhi)
                aux_file%dx(:) = (aux_file%upper(:) - aux_file%lower(:)) / (aux_file%num_cells(:) - 1)

#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of AMRClaw."
                stop
#endif                

            case default
                print *, 'ERROR:  Unrecognized aux file type'
                print *, '    file_type = ', aux_file%file_type
                print *, '    for aux file:', aux_file%path
                stop
        end select

        write(GEO_PARM_UNIT, *) '  mx = ',aux_file%num_cells(1),'  x = (',aux_file%lower(1),',',aux_file%upper(1),')'
        write(GEO_PARM_UNIT, *) '  my = ',aux_file%num_cells(2),'  y = (',aux_file%lower(2),',',aux_file%upper(2),')'
        write(GEO_PARM_UNIT, *) '  dx, dy (meters/degrees) = ', aux_file%dx

    end subroutine read_aux_file_header


    ! ========================================================================
    !
    !  Read auxinit file.
    !
    ! ========================================================================
    subroutine read_aux_file_data(aux_file)

#ifdef NETCDF
        use netcdf
#endif
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Input/Output
        type(aux_file_type), intent(in out) :: aux_file

        ! Locals
        ! integer, parameter :: iunit = 19, miss_unit = 17
        ! double precision, parameter :: auxinit_missing = -150.d0
        ! logical, parameter :: maketype2 = .false.
        ! integer :: i,j,num_points,missing,status,auxinit_start
        ! double precision :: no_data_value,x,y,z

        ! Locals
        integer, parameter :: unit = 8, missing_unit = 17
        integer :: i, j, status

        print *, ' '
        print *, 'Reading aux file  ', aux_file%path

        select case(abs(filetype))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                allocate(buffer(aux_file%num_cells(1) * aux_file%num_cells(2)))

                open(unit=unit, file=aux_file%path, status='unknown', &
                     form='formatted')
                i = 0
                status = 0
                do while (status == 0)
                    i = i + 1
                    read(unit, fmt=*, iostat=status) x, y, buffer(i)
                enddo

                aux_file%data = reshape(buffer, aux_file%num_cells)

                close(unit)

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if filetype=2 or
            ! mx values per line if filetype=3
            ! ================================================================
            case(2:3)
                open(unit=unit, file=aux_file%path, status='unknown', &
                     form='formatted')

                ! Read header
                do i = 1, 6
                    read(iunit,*)
                enddo

                ! Read in data
                missing = 0
                select case(abs(aux_file%file_type))
                    case(2)
                        do i=1, mx * my
                            read(iunit,*) buffer(i)
                            if (buffer(i) == aux_file%no_data_value) then
                                missing = missing + 1
                                buffer(i) = aux_file%missing_value
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(unit, *) (aux_file%data(i, j), i = 1, mx)
                            do i=1,mx
                                if (aux_file%data(i, j) == no_data_value) then
                                    missing = missing + 1
                                    aux_file%data(i, j) = aux_file%missing_value
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print *, '   WARNING...some missing data values this file'
                    print *, '       ',missing,' missing data values'
                    print *, '              (see fort.missing)'
                    print *, '   These values have arbitrarily been set to ',&
                        auxinit_missing
                endif

                close(unit)

            case(4)
#ifdef NETCDF
                stop "NetCDF reading not implemented"
#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of AMRClaw."
                stop
#endif
        end select

        close(unit=iunit)

   end subroutine read_auxinit

#ifdef NETCDF
    ! Check error return from NetCDF routine
    subroutine check_netcdf_error(ios)

        use netcdf

        implicit none

        integer, intent(in) :: ios

        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if

    end subroutine check_netcdf_error
#endif

end module aux_module