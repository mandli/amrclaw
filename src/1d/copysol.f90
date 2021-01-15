! ----------------------------------------------------------
!
! copy solution into grid with different number ghost cells
!
! ----------------------------------------------------------
subroutine copysol(valbig, val, nvar, mitot, nghost, midub, ngbig)

    implicit none

    integer, intent(in) :: nvar, mitot, nghost, midub, ngbig
    real(kind=8), intent(in) :: val(nvar, mitot)
    real(kind=8), intent(inout) :: valbig(nvar, midub)

    integer :: i

    do i = nghost + 1, mitot - nghost
        valbig(:, i - nghost + ngbig) = val(:, i)
    end do

end subroutine copysol

