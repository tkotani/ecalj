      module m_ldau
      complex(8),allocatable,protected::  vorb(:,:,:,:)

      complex(8),allocatable,private::  dmato(:,:,:,:)
      logical,private:: init=.true.
      contains
!! LDA+U initialization
      subroutine m_ldau_init()
      use m_lmfinit,only: nlibu,nsp,lmaxu,ctrl_quit,lmaxu,nsp,nlibu
      integer:: i,idmatu
      idmatu=nsp*nlibu*(lmaxu*2+1)**2
      if(init) then
      allocate( dmato(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      allocate( vorb (-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      init=.false.
      endif
      call sudmtu(dmato, vorb)  !read dmatu from dmatu.ext and generate vorbdmat from dmatu.ext
      end 
!! vorb LDA+U mixing
      subroutine m_ldau_vorbset(dmatu)
      use m_chkdmu,only: Chkdmu
      use m_lmfinit,only: nlibu,nsp,lmaxu,ctrl_quit,lmaxu,nsp,nlibu
      complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      call chkdmu(dmatu, dmato, vorb)
      dmato=dmatu
      end
      end module
