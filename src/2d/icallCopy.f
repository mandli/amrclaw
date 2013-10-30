c
c --------------------------------------------------------------------
c
       subroutine icallCopy(val,aux,nrow,ncol,nvar,naux,
     .                  ilo,ihi,jlo,jhi,level,iputst,jputst,
     .                  auxflags,mgrid)

      use amr_module
       implicit double precision (a-h, o-z)

       dimension val(nvar,nrow,ncol)
       dimension aux(naux,nrow,ncol)
       integer (kind=1) auxflags(nrow,ncol)

       logical sticksout


c NEW INDEX ORDERING
       iadd   (ivar,i,j) = loc    + ivar-1 + nvar*((j-1)*mitot+i-1)
       iaddaux(ivar,i,j) = locaux + ivar-1 + naux*((j-1)*mitot+i-1)

c ::::::::::::::::::::::::::: ICALL :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    intersecting grids to both val and aux arrays.
c
c    use larger definition of grids here - boundary data already in.
c    aux arrays also enlarged size.
c
c    iputst, jputst: where to copy values into. may not be in
c                    location corresponding to ilo,ihi,etc. if
c                    the patch has been periodically wrapped.

c    
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       auxflags = 0
       
       mptr = lstart(level)

 10    if (mptr .eq. 0) go to 99
          iglo = node(ndilo,mptr) 
          ighi = node(ndihi,mptr) 
          jglo = node(ndjlo,mptr) 
          jghi = node(ndjhi,mptr) 

c         # does it intersect?
c  INCLUDING GHOST CELLS?? AT LEAST FOR COPYING AUX, IF NOT Q
          ixlo = max(iglo-nghost,ilo-nghost)
          ixhi = min(ighi+nghost,ihi+nghost)
          jxlo = max(jglo-nghost,jlo-nghost)
          jxhi = min(jghi+nghost,jhi+nghost)
c  how did ghost cells get in the allowable region? They are not filled
c  (since we may be interpolating from newly filled grids, not just grids
c  that have been primed with bcs to be advanced.
!--          ixlo = max(iglo,ilo)
!--          ixhi = min(ighi,ihi)
!--          jxlo = max(jglo,jlo)
!--          jxhi = min(jghi,jhi)


          if (ixlo .le. ixhi .and. jxlo .le. jxhi) then
              loc  = node(store1,mptr)
              locaux = node(storeaux,mptr)
              nx   = ighi - iglo + 1
              ny   = jghi - jglo + 1
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              do 30 j    = jxlo, jxhi
              do 30 i    = ixlo, ixhi
              do 20 ivar = 1, nvar
                  ialloc  =  iadd(ivar,i-iglo+nghost+1,j-jglo+nghost+1)
                  val(ivar,i-ilo+iputst,j-jlo+jputst)  =  alloc(ialloc)
 20           continue
!--              if (i-ilo+iputst .lt. 0 .or. i-ilo+iputst .gt.nrow .or.
!--     .           j-jlo+jputst .lt. 0 .or. j-jlo+jputst .gt.ncol) then
!--                    write(*,*)" overwriting memory ehrre"
!--              endif
              auxflags(i-ilo+iputst,j-jlo+jputst) = 1   ! set flag
              do 25 iaux = 1, naux
                  ialloc = iaddaux(iaux,i-iglo+nghost+1,j-jglo+nghost+1)
!--                  temp = alloc(ialloc)
!--                  if  ((aux(iaux,i-ilo+iputst,j-jlo+jputst) .ne. temp) 
!--     .               .and. iaux .le. 4)
!--     .                   write(*,*)" stop here"
                  aux(iaux,i-ilo+iputst,j-jlo+jputst)  =  alloc(ialloc)
 25           continue
 30           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   continue

c
c if enlarged patch sticks out but isnt periodic just stick values
c in there so qad and src1d doesnt complain. the values wont be used 
c dont fill using values where debuggers might complain about uninitialized values
      if (ilo .lt. 0 .or. ihi .ge. iregsz(level) .or.
     &    jlo .lt. 0 .or. jhi .ge. jregsz(level)) then
          sticksout = .true.
      else
          sticksout = .false.
      endif

      if (sticksout .and. .false.) then  !think about ghost cells . . . DEBUG
         if (xperdom .or. yperdom .or. spheredom) then
           write(*,*)" should not be in this code with periodic bc"
           write(*,*)"not writing into val correctly using mapping "
	         stop
         endif
         if (ilo .eq. -1) then
            do j = max(jlo,0), min(jhi,jregsz(level)-1)
              jj = jputst + j - jlo
              do ivar = 1, nvar
c               it has got to be the first row that sticks out
c               since called from filval with enlarged patch and
c               no ghost cells
                val(ivar,1,jj) = val(ivar,2,jj) 
              end do
              do iaux = 1, naux
c                do same for aux arrays
                 aux(iaux,1,jj) = aux(iaux,2,jj) 
              end do
            end do
         endif
         if (ihi .eq. iregsz(level)) then
            do j = max(jlo,0), min(jhi,jregsz(level)-1)
              jj = jputst + j - jlo
              do ivar = 1, nvar
                 val(ivar,nrow,jj) = val(ivar,nrow-1,jj) 
              end do
              do iaux = 1, naux
                 aux(iaux,nrow,jj) = aux(iaux,nrow-1,jj) 
              end do
            end do
         endif
         if (jlo .eq. -1) then
            do i = max(ilo,0), min(ihi,iregsz(level)-1)
              ii = iputst + i - ilo
              do ivar = 1, nvar
                val(ivar,ii,1) = val(ivar,ii,2) 
              end do
              do iaux = 1, naux
                aux(iaux,ii,1) = aux(iaux,ii,2) 
              end do
            end do
         endif
         if (jhi .eq. jregsz(level)) then
            do i = max(ilo,0), min(ihi,iregsz(level)-1)
              ii = iputst + i - ilo
              do ivar = 1, nvar
                 val(ivar,ii,ncol) = val(ivar,ii,ncol-1) 
              end do
              do iaux = 1, naux
                 aux(iaux,ii,ncol) = aux(iaux,ii,ncol-1) 
              end do
            end do
         endif
      endif

      return
      end
