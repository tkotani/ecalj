      module m_readVcoud        ! Coulomb matrix for given q
      use m_readqg,only: Readqg0
      use m_rdpp,only: nbloch
      implicit none
!-----------------------------------      
      integer,protected :: ngb, ngc, iq
      real(8),protected :: q(3)
      real(8),allocatable,protected    :: vcousq(:), vcoud(:) !,vcousq0(:),vcoudummy(:)
      complex(8),allocatable,protected :: zcousq(:,:) !,zcousqrsum(:,:,:),zcousqr(:,:)
      contains
!-----------------------------------      
      subroutine Readvcoud(q,iq,NoVcou)
      intent(in)::         q,iq,NoVcou
!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
!! ===Readin diagonalized Coulomb interaction===
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!   A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!   Vould file is written in hvccfp0.m.F.
!! NOTE: vcoud and zcousq are used in m _zmelt.
      integer:: iq
      logical:: NoVcou
      real(8):: q(3),quu(3)
      integer:: ngb0,ifvcoud
      real(8):: qvv(3)
      character(128):: vcoudfile
      character(10) :: i2char
      call Readqg0('QGcou',q,   quu,ngc) ! ngc: the number of IPW for the interaction matrix (in QGcou),
                                         ! quu: equivalent to q generated in qg4g4.
      ngb = ngc+nbloch
      if(NoVcou) return
      vcoudfile='Vcoud.'//trim(i2char(iq)) ! iq was iqqv this is closed at the end of do 1001
      open(newunit=ifvcoud, file=trim(vcoudfile), action='read',form='unformatted')
      read(ifvcoud) ngb0
      read(ifvcoud) qvv
      if(sum(abs(qvv-q))>1d-10) then
         write(6,*)'qvv =',qvv
         call rx( 'readvcoud: qvv/=0 is not consistent')
      endif
      if( ngb0/=ngb ) then      !sanity check
         write(6,*)' qxx ngb0 ngb=',q,ngb0,ngb
         call rx( 'readvcoud:ngb0/=ngb')
      endif   
      if(allocated(zcousq)) deallocate( zcousq,vcousq,vcoud )
      allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0))
      read(ifvcoud) vcoud
      read(ifvcoud) zcousq
      close(ifvcoud)
      vcousq = sqrt(vcoud)
      end subroutine

      end module
