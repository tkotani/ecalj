      subroutine bndconn_v2(alat,plat,qlat, 
     &      nbas, rmt, bas, ndimh, ldim2,
     &      evl,  ngp, ngp_p,  ngvecp,ngvecp_p,   geig,geig_p,   
     &      cphi, cphi_p,
     o iiyf,ovv)
c-Get the connenctivity of eigenfunctions along the q points.
co iiyf: connection. idim'th band is mapped to iiyf(idim) band.
co ovv : information of a deformed overlap matrix.
C     implicit none
      integer(4):: ndimh,ldim2,nbas,ngp,ngp_p,iix,iiy,iyiy,istart,j
      integer(4) :: iiyf(ndimh),ifdebug
     &  ,ngvecp(3,ngp),ngvecp_p(3,ngp_p),ifill(ndimh)
      real(8) ::  eigdeg, sumdeg,alat,plat(3,3),qlat(3,3)
     &           ,ovv(ndimh,ndimh),rmt(nbas),bas(3,nbas),evl(ndimh)
      complex(8):: cphi(ldim2,ndimh), cphi_p(ldim2,ndimh)
     & ,geig(ngp, ndimh),geig_p(ngp_p, ndimh)
      complex(8),allocatable::  gpg(:,:),ppovl2(:,:)
c----

      ifdebug =1198
      allocate(gpg(ndimh,ndimh), ppovl2(ngp,ngp_p))
c      print * ,' goto mkppovl2'
c     & ,sum(abs(geig)),sum(abs(geig_p))
c--- Get the connection. ovelap matrix of eigenfunctions
      call mkppovl2(alat,plat,qlat, 
     &      ngp,   ngvecp,   
     &      ngp_p, ngvecp_p, 
     &      nbas, rmt, bas,  
     o      ppovl2)
c      print * ,' end of mkppovl2'
c     & ,sum(abs(geig)),sum(abs(geig_p)),sum(abs(ppovl2))
      gpg = matmul(dconjg(transpose(geig)),matmul(ppovl2,geig_p))
      write(ifdebug,'("     ",255i4)') (j,j=1,ndimh)
c      print *, ' band mode xxx4'
      do iix =1,ndimh  !; print *, ' band mode xxx5 iix=',iix
        do iiy =1,ndimh
          ovv(iiy,iix) =abs( gpg(iiy,iix) 
     &  +sum( dconjg(cphi(1:ldim2,iiy))*cphi_p(1:ldim2,iix)) 
     &   )**2
        enddo
        iyiy = maxloc( ovv(1:ndimh,iix),dim=1 )
        write(ifdebug,'(2i3,255i4)') iix ,iyiy,
     &      (int(100*ovv(j,iix)),j=1,ndimh)
      enddo
      deallocate(ppovl2,gpg)
c      print *, ' band mode xxx5'

      ifill = 0
      do iix = 1,ndimh
!  Denenaracy treatment
        eigdeg = -1d99
        istart = 1
        sumdeg = 0d0
c      print *, ' band mode xxx6'
        do iiy=1,ndimh
          if(iiy==ndimh) then
            ovv(istart:ndimh,iix) = sumdeg + ovv(iiy,iix)
          elseif(abs(eigdeg-evl(iiy))>1d-6.and.iiy/=1) then
            ovv(istart:iiy-1,iix) = sumdeg
            eigdeg = evl(iiy)
            istart = iiy
            sumdeg = ovv(iiy,iix)
          else
            sumdeg = sumdeg + ovv(iiy,iix)
          endif
        enddo
c       print *, ' band mode xxx7'
c            write(ifdebug,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
c            write(ifdebug,'(2i3,255i4)') iix ,iix,
c     &      (int(100*ovv(j,iix)),j=1,ndimh)
! Get the connection ix to iy.
        do
          iyiy = maxloc( ovv(1:ndimh,iix),dim=1 )
          if( ifill(iyiy)==0 ) then
            iiyf(iix) = iyiy
            ifill(iyiy)=1
            ovv(1:iyiy-1    ,iix) = 0d0
            ovv(iyiy+1:ndimh,iix) = 0d0
            exit
          else
            ovv(iyiy,iix) = -1d0
          endif
        enddo
c            write(ifdebug,'(2i3,255i4)') iix ,iyiy,
c     &      (int(100*ovv(j,iix)),j=1,ndimh)
c      print *, ' band mode xxx8'
      enddo
c      print *, ' end of bndconn'
      end
