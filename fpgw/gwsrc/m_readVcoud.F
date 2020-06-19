      module m_readVcoud
      implicit none
      real(8),allocatable,protected    :: vcousq(:)   !,vcousq0(:),vcoudummy(:)
      complex(8),allocatable,protected :: zcousq(:,:) !,zcousqrsum(:,:,:),zcousqr(:,:)
      contains
      
      subroutine readvcoud(q,iq,ngb)
      intent(in)::         q,iq,ngb
!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
      integer:: iq,ngb
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
         call rx( 'hx0fp0_sc: qvv/=0 is not consistent')
      endif
      if(allocated(zcousq)) deallocate( zcousq,vcousq )
      allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
      read(ifvcoud) vcousq
      read(ifvcoud) zcousq
      close(ifvcoud)
      vcousq = sqrt(vcousq)
      end subroutine

      end module
