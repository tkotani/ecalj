      module m_readVcoud
      implicit none
      real(8),allocatable,protected    :: vcousq(:),vcoud(:)    !,vcousq0(:),vcoudummy(:)
      complex(8),allocatable,protected :: zcousq(:,:) !,zcousqrsum(:,:,:),zcousqr(:,:)
      contains
      
      subroutine Readvcoud(q,iq,ngb)
      intent(in)::         q,iq
!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
      integer:: iq
      integer,optional ::ngb !dummy
      real(8):: q(3)
      integer:: ngb0,ifvcoud
      real(8):: qvv(3)
      character(128):: vcoudfile
      character(5) :: charnum5
      vcoudfile='Vcoud.'//charnum5(iq) ! iq was iqqv this is closed at the end of do 1001
      open(newunit=ifvcoud, file=trim(vcoudfile), action='read',form='unformatted')
      read(ifvcoud) ngb0
      read(ifvcoud) qvv
      if(sum(abs(qvv-q))>1d-10) then
         write(6,*)'qvv =',qvv
         call rx( 'readvcoud: qvv/=0 is not consistent')
      endif
c      if( ngb0/=ngb ) then      !sanity check
c         write(6,*)' qxx ngb0 ngb=',q,ngb0,ngb
c         call rx( 'readvcoud:ngb0/=ngb')
c      endif   
      if(allocated(zcousq)) deallocate( zcousq,vcousq,vcoud )
      allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0))
      read(ifvcoud) vcoud
      read(ifvcoud) zcousq
      close(ifvcoud)
      vcousq = sqrt(vcoud)
      end subroutine

      end module
