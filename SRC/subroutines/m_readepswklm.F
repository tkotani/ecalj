c$$$      module m_readepswklm
c$$$      use m_readq0p,only: nq0i
c$$$      implicit none
c$$$      
c$$$      real(8),allocatable,protected:: dmlx(:,:),epinvq0i(:,:),wklm(:)
c$$$      integer,protected:: lxklm
c$$$      
c$$$      contains
c$$$!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
c$$$      subroutine readepswklm()
c$$$      integer:: ifidmlx,nq0ix
c$$$!     ! Readin Vcoud and EPSwklm for newaniso()=T ===
c$$$      open(newunit=ifidmlx,file='EPSwklm',form='unformatted')
c$$$      read(ifidmlx) nq0ix,lxklm
c$$$      if(nq0i/=nq0ix) then
c$$$        write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
c$$$        call rx( 'nq0i from EPSwklm /= nq0i')
c$$$      endif
c$$$      allocate( dmlx(nq0i,9))
c$$$      allocate( epinvq0i(nq0i,nq0i) )
c$$$      allocate( wklm((lxklm+1)**2))
c$$$      read(ifidmlx) dmlx, epinvq0i
c$$$      read(ifidmlx) wklm
c$$$      close(ifidmlx)
c$$$      end subroutine
c$$$      end module
