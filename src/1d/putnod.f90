! :::::::::::::::::::::::::::::: PUTNOD :::::::::::::::::::::;
!
!  return mptr node to the linked list kept in node array.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
subroutine putnod(mptr)

    use amr_module, only: ndfree, node, nextfree
    implicit none

    integer, intent(in) :: mptr

    node(nextfree, mptr) = ndfree
    ndfree = mptr

end subroutine putnod


! :::::::::::::::::::::::::::::: PUTNOD_BND :::::::::::::::::::::;
!
!  return bndry list node to the linked list 
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
subroutine putnod_bnd (mcell)

    use amr_module, only: ndfree_bnd, bndList, nextfree
    implicit none

    integer, intent(in) :: mcell

    bndList(mcell, nextfree) = ndfree_bnd
    ndfree_bnd = mcell

end subroutine putnod_bnd