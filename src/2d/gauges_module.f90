! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine setbestsrc
!     Called each time regridding is done to determine which patch to 
!     use for interpolating to each gauge location.
!   subroutine print_gauges
!     Called each time step for each grid patch.
!     Refactored dumpgauge routine to interpolate for all gauges on patch.
!
!     Note: by default all components of q are printed at each gauge.
!     To print something different or a different precision, modify 
!     format statement 100 and/or the write statement that uses it.
!   
! Note: Updated for Clawpack 5.3.0:
!   - the dumpgauge and setbestsrc subroutines have been moved to this module 
!     and the dumpgauge subroutine has been refactored and renamed print_gauges.
!   - dumpgauge.f must be removed from Makefiles.
!   - setbestsrc uses quicksort to sort gauge numbers and
!     then figures out which gauges will be updated by grid, and stores this
!     information in new module variables mbestg1, mbestg2.
!   - print_gauges no longer uses binary search to locate first gauge handled
!     by a grid.  Instead loop over gauges specified by mbestg1, mbestg2.
!
! Note: Updated for Clawpack 5.4.0
!   - refactor so each gauge writes to its own file, and batches the writes instead of 
!     writing one at a time. This will remove the critical section and should speed up gauges a lot
!   - When array is filled, that gauge will write to file and start over. 
!   - Need to save index so know position in array where left off
!   - At checkpoint times, dump all gauges

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT = 89
    integer :: num_gauges

!     integer, parameter :: MAX_BUFFER = 1000
    integer, parameter :: MAX_BUFFER = 10

    ! Gauge data types
    type gauge_type
        ! Gauge number
        integer :: gauge_num

        character(len=14) :: file_name

        ! Location in time and space
        real(kind=8) :: x, y, t_start, t_end

        ! Output settings
        integer :: file_format
        character(len=10) :: display_format
        logical, allocatable :: q_out_vars(:)
        logical, allocatable :: aux_out_vars(:)
        integer :: num_out_vars

        ! Data buffers - data holds output and time
        real(kind=8), allocatable :: data(:, :)
        integer :: level(MAX_BUFFER)

        ! Where we are in the buffer
        integer :: buffer_index
    end type gauge_type

    ! Gague array
    type(gauge_type), allocatable :: gauges(:)

    ! real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable, dimension(:) ::  mbestsrc, mbestorder, &
                          igauge, mbestg1, mbestg2

contains

    subroutine set_gauges(restart, nvar, fname)

        use amr_module, only: maxgr
        use utility_module, only: get_value_count_alt, parse_values

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname
        logical, intent(in) :: restart
        integer, intent(in) :: nvar

        ! Locals
        integer :: i, ipos, idigit, inum, n, num_fields, index
        integer, parameter :: iunit = 7
        character(len=128) :: line
        character(len=4) temp_string
        character(len=20) :: q_column, aux_column
        real(kind=8) :: values(10)

        if (.not.module_setup) then
            ! Open file
            if (present(fname)) then
                call opendatafile(iunit,fname)
            else
                call opendatafile(iunit,'gauges.data')
            endif

            read(iunit,*) num_gauges

            allocate(gauges(num_gauges))
            allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
            allocate(mbestg1(maxgr), mbestg2(maxgr))
            
            mbestsrc = 0

            ! Original gauge information
            do i=1,num_gauges
                read(iunit,*) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, &
                              gauges(i)%t_start,gauges(i)%t_end
                gauges(i)%buffer_index = 1
            enddo

            ! Read in output formats
            read(iunit, *)
            read(iunit, *) (gauges(i)%file_format, i=1, num_gauges)
            read(iunit, *) (gauges(i)%display_format, i=1, num_gauges)

            ! Read in q fields
            read(iunit, *)
            do i=1, num_gauges
                read(iunit, "(a)") line
                num_fields = get_value_count_alt(line)
                allocate(gauges(i)%q_out_vars(num_fields))
                read(line, *) gauges(i)%q_out_vars
                
                gauges(i)%num_out_vars = 0
                do n=1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if  
                end do
            end do

            ! Read in aux fields
            read(iunit, *)
            do i=1, num_gauges
                read(iunit, "(a)") line
                num_fields = get_value_count_alt(line)
                allocate(gauges(i)%aux_out_vars(num_fields))
                read(line, *) gauges(i)%aux_out_vars

                do n=1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if  
                end do

                ! Allocate data buffer - time + num_out_vars + eta
                allocate(gauges(i)%data(gauges(i)%num_out_vars + 1, MAX_BUFFER))
                print *, "data size ", i, ": (", gauges(i)%num_out_vars + 1, ", ", MAX_BUFFER, ")"
            end do

            close(iunit)
            
            ! Create gauge output file names
            do i = 1, num_gauges
                gauges(i)%file_name = 'gaugexxxxx.txt'
                inum = gauges(i)%gauge_num
                do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
                    idigit = mod(inum,10)
                    gauges(i)%file_name(ipos:ipos) = char(ichar('0') + idigit)
                    inum = inum / 10
                end do

                ! status unknown since might be a restart run. maybe need to test and rewind?
                if (restart) then
                    open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name, status='old',        &
                         position='append', form='formatted')
                else
                    open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name, status='unknown',        &
                       position='append', form='formatted')
                    rewind OUTGAUGEUNIT
                    write(OUTGAUGEUNIT, 100) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, gauges(i)%num_out_vars + 1
 100                format("# gauge_id= ",i5," location=( ",1e15.7," ",1e15.7," ) num_eqn= ",i2)
                    ! write(OUTGAUGEUNIT, 101) 

                    ! Construct column labels
                    index = 0
                    q_column = "["
                    do n=1, size(gauges(i)%q_out_vars, 1)
                        if (gauges(i)%q_out_vars(n)) then
                            write(q_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if  
                    end do
                    q_column(3 * index + 2:4 + 3 * index) = "]"

                    aux_column = "["
                    index = 0
                    do n=1, size(gauges(i)%aux_out_vars, 1)
                        if (gauges(i)%aux_out_vars(n)) then
                            write(aux_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if  
                    end do
                    aux_column(3 * index + 2:4 + 3 * index) = "]"

                    write(OUTGAUGEUNIT, *) "# level, time, q",             &
                                           trim(q_column), ", aux",       &
                                           trim(aux_column), ", eta"
                endif

                close(OUTGAUGEUNIT)

            end do

            module_setup = .true.
        end if

    end subroutine set_gauges


!
! --------------------------------------------------------------------
!
    subroutine setbestsrc()
!
!     Called every time grids change, to set the best source grid patch
!     for each gauge, i.e. the finest level patch that includes the gauge.
!
!     lbase is grid level that didn't change, but since fine
!     grid may have disappeared, we still have to look starting
!     at coarsest level 1.
!
        use amr_module
        implicit none

        integer :: lev, mptr, i, k1, ki

!
! ##  set source grid for each loc from coarsest level to finest.
! ##  that way finest src grid left and old ones overwritten
! ##  this code uses fact that grids do not overlap

! # for debugging, initialize sources to 0 then check that all set
        do i = 1, num_gauges
            mbestsrc(i) = 0
        end do

        do lev = 1, lfine  
            mptr = lstart(lev)
            do
                do i = 1, num_gauges
                    if ((gauges(i)%x >= rnode(cornxlo,mptr)) .and. &
                        (gauges(i)%x <= rnode(cornxhi,mptr)) .and. &  
                        (gauges(i)%y >= rnode(cornylo,mptr)) .and. &
                        (gauges(i)%y <= rnode(cornyhi,mptr)) ) then
                        mbestsrc(i) = mptr
                    end if
                end do
                mptr = node(levelptr, mptr)
                if (mptr == 0) exit
            end do 
        end do

        do i = 1, num_gauges
            if (mbestsrc(i) .eq. 0) &
               print *, "ERROR in setting grid src for gauge data", i
        end do

        ! Sort the source arrays for easy testing during integration
        call qsorti(mbestorder, num_gauges, mbestsrc)

!     After sorting,  
!           mbestsrc(mbestorder(i)) = grid index to be used for gauge i
!     and mbestsrc(mbestorder(i)) is non-decreasing as i=1,2,..., num_gauges

!     write(6,*) '+++ mbestorder: ',mbestorder
!     write(6,*) '+++ mbestsrc: ',mbestsrc

!     Figure out the set of gauges that should be handled on each grid:  
!     after loop below, grid k should handle gauges numbered
!          mbestorder(i) for i = mbestg1(k), mbestg1(k)+1, ..., mbestg2(k)
!     This will be used for looping in print_gauges subroutine.

      ! initialize arrays to default indicating grids that contain no gauges:
        mbestg1 = 0
        mbestg2 = 0

        k1 = 0
        do i=1,num_gauges
            ki = mbestsrc(mbestorder(i))
            if (ki > k1) then
                ! new grid number seen for first time in list
                if (k1 > 0) then
                    ! mark end of gauges seen by previous grid
                    mbestg2(k1) = i-1
!                   write(6,*) '+++ k1, mbestg2(k1): ',k1,mbestg2(k1)
                endif
                mbestg1(ki) = i
!               write(6,*) '+++ ki, mbestg1(ki): ',ki,mbestg1(ki)
            endif
            k1 = ki
        enddo
        if (num_gauges > 0) then
            ! finalize 
            mbestg2(ki) = num_gauges
        endif

    end subroutine setbestsrc

!
! -------------------------------------------------------------------------
!
    subroutine update_gauges(q,aux,xlow,ylow,nvar,mitot,mjtot,naux,mptr)
!
!     This routine is called each time step for each grid patch, to output
!     gauge values for all gauges for which this patch is the best one to 
!     use (i.e. at the finest refinement level).  

!     It is called after ghost cells have been filled from adjacent grids
!     at the same level, so bilinear interpolation can be used to 
!     to compute values at any gauge location that is covered by this grid.  

!     The grid patch is designated by mptr.
!     We only want to set gauges i for which mbestsrc(i) == mptr.
!     The array mbestsrc is reset after each regridding to indicate which
!     grid patch is best to use for each gauge.

!     This is a refactoring of dumpgauge.f from Clawpack 5.2 
!     Loops over only the gauges to be handled by this grid, as specified
!     by indices from mbestg1(mptr) to mbestg2(mptr)

        use amr_module
  
        implicit none
  
        ! Input
        real(kind=8), intent(in) ::  q(nvar,mitot,mjtot)
        real(kind=8), intent(in) ::  aux(naux,mitot,mjtot)
        real(kind=8), intent(in) ::  xlow,ylow
        integer, intent(in) ::  nvar,mitot,mjtot,naux,mptr
  
        ! local variables:
        integer, parameter :: MAX_VARS = 20
        real(kind=8) :: var(MAX_VARS)
        real(kind=8) :: xcent,ycent,xoff,yoff,tgrid,hx,hy
        integer :: level,i,j,ioff,joff,iindex,jindex,ivar, ii,i1,i2
        integer :: icell,jcell, index
        integer :: var_index
  
!       write(*,*) '+++ in print_gauges with num_gauges, mptr = ',num_gauges,mptr
  
        if (num_gauges == 0) then
           return
        endif
  
        i1 = mbestg1(mptr)
        i2 = mbestg2(mptr)
  
        if (i1 == 0) then
           ! no gauges to be handled by this grid
           return
        endif
  
!       write(6,*) '+++ mbestg1(mptr) = ',mbestg1(mptr)
!       write(6,*) '+++ mbestg2(mptr) = ',mbestg2(mptr)
  
!       # this stuff the same for all gauges on this grid
        tgrid = rnode(timemult,mptr)
        level = node(nestlevel,mptr)
        hx    =  hxposs(level)
        hy    =  hyposs(level)

!     write(*,*) 'tgrid = ',tgrid

        do i = i1,i2
            ii = mbestorder(i)
            if (mptr /= mbestsrc(ii)) then
                print *, '*** should not happen... i, ii, mbestsrc(ii), mptr:'
                print *, i, ii, mbestsrc(ii), mptr
                stop
            endif
            if (tgrid < gauges(ii)%t_start .or. tgrid > gauges(ii)%t_end) then
                ! This gauge should not be output at this time
                cycle
            endif

            ! Bilinear interpolation at gauge location
            ! Note: changed 0.5 to  0.5d0 etc.
            iindex =  int(.5d0 + (gauges(ii)%x - xlow) / hx)
            jindex =  int(.5d0 + (gauges(ii)%y - ylow) / hy)
            if ((iindex < nghost .or. iindex > mitot-nghost) .or. &
                (jindex < nghost .or. jindex > mjtot-nghost)) then
                    print *, "ERROR in output of Gauge Data "
            end if
            xcent  = xlow + (iindex - 0.5d0) * hx
            ycent  = ylow + (jindex - 0.5d0) * hy
            xoff   = (gauges(ii)%x - xcent) / hx
            yoff   = (gauges(ii)%y - ycent) / hy
    !  IF WANT TO USE, MODIFY TO TEST FOR ROUNDOFF LEVEL DIFF
    !       if (xoff .lt. 0.d0 .or. xoff .gt. 1.d0 .or. &
    !           yoff .lt. 0.d0 .or. yoff .gt. 1.d0) then
    !          write(6,*)" BIG PROBLEM in DUMPGAUGE", i
    !       endif

    !       ## bilinear interpolation
            var_index = 0
            do ivar - 1, size(gauges(ii)%q_out_vars, 1)
                if (gauges(ii)%q_out_vars(ivar)) then
                    var_index = var_index + 1
                    var(var_index) = (1.d0 - xoff) * (1.d0 - yoff) * q(ivar, iindex, jindex) &
                                        + xoff * (1.d0 - yoff) * q(ivar, iindex + 1, jindex) &
                                        + (1.d0 - xoff) * yoff * q(ivar, iindex, jindex + 1) &
                                        + xoff * yoff * q(ivar, iindex + 1, jindex + 1)
                endif
                if (gauge(ii)%aux_out_vars(ivar)) then
                    var_index = var_index + 1
                    var(var_index) = (1.d0 - xoff) * (1.d0 - yoff) * aux(ivar, iindex, jindex) &
                                        + xoff * (1.d0 - yoff) * aux(ivar, iindex + 1, jindex) &
                                        + (1.d0 - xoff) * yoff * aux(ivar, iindex, jindex + 1) &
                                        + xoff * yoff * aux(ivar, iindex + 1, jindex + 1)
                endif
            end do

            ! Check to make sure we grabbed all the values
            if (gauges(ii)%num_out_vars /= var_index) then
                print *, gauges(ii)%num_out_vars, var_index
                print *, gauges(ii)%q_out_vars
                print *, gauges(ii)%aux_out_vars
                stop "Somehow we did not grab all the values we wanted..."
            end if

            ! Zero out tiny values to prevent underflow problems
            do j = 1, gauges(ii)%num_out_vars
                if (abs(var(j)) < 1d-90) var(j) = 0.d0
            end do

           ! save info for this time
           index = gauges(ii)%buffer_index
     
            gauges(ii)%level(index) = level
            gauges(ii)%data(1,index) = tgrid
            do j = 1, gauges(ii)%num_out_vars
                gauges(ii)%data(1 + j, index) = var(j)
            end do
            gauges(ii)%data(gauges(ii)%num_out_vars + 2, index) = eta
            
            gauges(ii)%buffer_index = index + 1
            if (gauges(ii)%buffer_index > MAX_BUFFER) then
                call print_gauges_and_reset_nextLoc(ii)  
            endif

        end do ! end of loop over gauges
 
      end subroutine update_gauges
!
! -------------------------------------------------------------------------
!
    subroutine print_gauges_and_reset_nextLoc(gauge_num)
        ! Write out gauge data for the gauge specified

        implicit none

        integer, intent(in) :: gauge_num

        ! Locals
        integer :: j, k, myunit
        integer :: omp_get_thread_num, mythread
        character(len=32) :: out_format

        ! Open unit dependent on thread number
        mythread = 0
!$      mythread = omp_get_thread_num()
        myunit = OUTGAUGEUNIT + mythread

        ! ASCII output
        if (gauges(gauge_num)%file_format == 1) then
            ! Construct output format based on number of output variables and
            ! request format
            write(out_format, "(A7, i2, A6, A1)") "(i5.2,",         &
               gauges(gauge_num)%num_out_vars + 1, gauges(gauge_num)%display_format, ")"

            open(unit=myunit, file=gauges(gauge_num)%file_name, status='old', &
                              position='append', form='formatted')
          
            ! Loop through gauge's buffer writing out all available data.  Also
            ! reset buffer_index back to beginning of buffer since we are emptying
            ! the buffer here
            do j = 1, gauges(gauge_num)%buffer_index - 1
                write(myunit, out_format) gauges(gauge_num)%level(j),    &
                    (gauges(gauge_num)%data(k, j), k=1, gauges(gauge_num)%num_out_vars + 1)
            end do
            gauges(gauge_num)%buffer_index = 1                        

            ! close file
            close(myunit)
        else
            print *, "Unhandled file format ", gauges(gauge_num)%file_format
            stop
        end if

    end subroutine print_gauges_and_reset_nextLoc

end module gauges_module
