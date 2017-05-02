subroutine setaux(mbc, mx, my, xlower, ylower, dx, dy, maux, aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This default version sets aux files based on files.  To maintain this
    ! feature derived setaux.f90 files should include the content here
 
    use aux_module, only: num_aux_files, aux_files

    implicit none

    ! Input/Output
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: n, i, j, i_start, i_end, j_start, j_end
    real(kind=8) :: x_hi, y_hi, x, x_im, x_ip, y, y_im, y_ip
    real(kind=8) :: x_ip_c, x_im_c, x_c, y_ip_c, y_im_c, y_c
    real(kind=8) :: daux

    ! Zero out aux fields that are going to be set below
    ! Is this right in general?  Is this how a generic value is set?  Kind of
    ! feel like this should be variable but this makes it harder to do this for
    ! generic aux files then (friction I am thinking of).  We would have to
    ! always compute offsets.  Maybe use negative values for having to use a 
    ! different value and replace rather than add aux fields?
    do n = 1, num_aux_file
        aux(aux_files(n)%init_type, :, :) = 0.d0
    end do

    ! Loop through each file to see if it intersects this grid
    x_hi = xlower + (mx - 0.5d0) * dx
    y_hi = ylower + (my - 0.5d0) * dy
    do n = 1, num_aux_files

        if ((xlower <= aux_files(n)%lower(1)) .and. &
            (ylower <= aux_files(n)%lower(2)) .and. &
            (x_hi >= aux_files(n)%upper(1)) .and. &
            (y_hi >= aux_files(n)%upper(2))) then

            ! Compute indexing
            ! TODO:  Implement this
            stop
            ! xintlow = dmax1(xlow,xlowauxinit(mf))
            ! xinthi  = dmin1(xhigher,xhiauxinit(mf))
            ! istart  = max(1-mbc,int(0.5 + (xintlow-xlow)/dx))
            ! iend    = min(mx+mbc,int(1.0 + (xinthi-xlow)/dx))

            ! yintlow = dmax1(ylow,ylowauxinit(mf))
            ! yinthi  = dmin1(yhigher,yhiauxinit(mf))
            ! jstart  = max(1-mbc,int(0.5 + (yintlow-ylow)/dy))
            ! jend    = min(my+mbc,int(1.0 + (yinthi-ylow)/dy))

            ! Compute integrals
            do i = i_start, i_end
                x = xlower + (i - 0.5d0) * dx
                x_im = x - 0.5d0 * dx
                x_ip = x + 0.5d0 * dx
                do j = j_start, j_end
                    y = ylower + (i - 0.5d0) * dy
                    y_im = y - 0.5d0 * dy
                    y_ip = y + 0.5d0 * dy

                    if ((x_ip > aux_file(n)%lower(1)) .and. &
                        () .and. &
                        () .and. &
                        () ) then

                        x_ip_c = min(x_ip, aux_file(n)%upper(1))
                        x_im_c = min(x_ip, aux_file(n)%lower(1))
                        x_c = 0.5d0 * (x_ip_c + x_im_c)
                        y_ip_c = min(y_ip, aux_file(n)%upper(2))
                        y_im_c = min(y_ip, aux_file(n)%lower(2))
                        y_c = 0.5d0 * (y_ip_c + y_im_c)

                        daux = integral()
     !                 daux =topointegral(ximc,xc,xipc,yjmc,yc,yjpc,
     ! &                        xlowauxinit(mf),ylowauxinit(mf),
     ! &                        dxauxinit(mf),dyauxinit(mf),
     ! &                        mxauxinit(mf),myauxinit(mf),
     ! &                        auxinitwork
     ! &                      (i0auxinit(mf):i0auxinit(mf)+mauxinit(mf)-1)
     ! &                        ,1)

                        ! TODO: What if kappa is not set?
                        daux = daux / ((x_ip_c - x_im_c) * (y_jp_c - y_jm_c) * kappa)
                        aux(aux_files(n)%astu, i, j) = aux(aux_files(n)%astu, i, j) + daux
                    end if
                end do
            end do
        end if
    end do



end subroutine setaux