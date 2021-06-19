      module m_readepswklm
      use m_readq0p,only: nq0i
      implicit none
      
      real(8),allocatable,protected:: dmlx(:,:),epinvq0i(:,:),wklm(:)
      integer,protected:: lxklm
      
      contains
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
      subroutine readepswklm()
      integer:: ifidmlx,nq0ix
!     ! Readin Vcoud and EPSwklm for newaniso()=T ===
      open(newunit=ifidmlx,file='EPSwklm',form='unformatted')
      read(ifidmlx) nq0ix,lxklm
      if(nq0i/=nq0ix) then
        write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
        call rx( 'nq0i from EPSwklm /= nq0i')
      endif
      allocate( dmlx(nq0i,9))
      allocate( epinvq0i(nq0i,nq0i) )
      allocate( wklm((lxklm+1)**2))
      read(ifidmlx) dmlx, epinvq0i
      read(ifidmlx) wklm
      close(ifidmlx)
      end subroutine
      end module
