  !$$$      module m_readvmem
  !$$$c------ a virtual memory system -----
  !$$$c all these are private.
  !$$$      implicit none
  !$$$      real(8),allocatable,private:: idat(:,:),qp(:,:)
  !$$$      logical,private:: init(2)=.true.
  !$$$      real(8),private:: QpGcut_cou, QpGcut_psi
  !$$$      integer(4),private::   nqnumc,nqnump,ngcmx,ngpmx
  !$$$      integer(4),allocatable,private:: ngvecp(:,:,:),ngp(:),ngvecc(:,:,:),ngc(:)
  !$$$c-----------------------------------------------------------------
  !$$$      contains
  !$$$c----------------------------------
  !$$$      subroutine read_vmem(key,ngmx)
  !$$$c- get ngcmx or mgpmx
  !$$$      implicit none
  !$$$      integer(4):: ngmx,ifiqg=4052
  !$$$      character*(*) key
  !$$$      if    (key=='QGpsi') then
  !$$$        open(ifiqg, file='QGpsi',form='unformatted')
  !$$$        read(ifiqg) nqnump, ngpmx, QpGcut_psi
  !$$$        ngmx=ngpmx
  !$$$      elseif(key=='QGcou') then
  !$$$        open(ifiqg, file='QGcou',form='unformatted')
  !$$$        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
  !$$$        ngmx=ngcmx
  !$$$      else
  !$$$        stop "readngmx: key is not QGpsi QGcou"
  !$$$      endif
  !$$$      close(ifiqg)
  !$$$      end subroutine
  !$$$
  !$$$c--------------------------
  !$$$      subroutine readqg(key,qin,ginv,  qu,ngv,ngvec)
  !$$$c- Get ngv and ngvec(3,ngv) for given qin(3)
  !$$$cr key=='QGcou' or 'QGpsi'
  !$$$      implicit none
  !$$$      integer(4):: ngv, ngvec(3,*), ifi, iq,verbose
  !$$$      real(8):: qin(3),qu(3),ginv(3,3)
  !$$$      character*(*) key
  !$$$      if    (key=='QGpsi') then
  !$$$        ifi=1
  !$$$        if(verbose()>=80) write (6,"('readqg psi: qin=',3f8.3,i5)") qin
  !$$$      elseif(key=='QGcou') then
  !$$$        ifi=2
  !$$$        if(verbose()>=80) write (6,"('readqg cou: qin=',3f8.3,i5)") qin
  !$$$      else
  !$$$        stop "readqg: wrongkey"
  !$$$      endif
  !$$$
  !$$$      if(init(ifi)) then
  !$$$        call init_readqg(ifi)
  !$$$        init(ifi)=.false.
  !$$$      endif
  !$$$
  !$$$
  !$$$      if(ifi==1) then
  !$$$        call iqindx2(qin,ginv, qp,nqnump, iq,qu)
  !$$$c        do iq=1, nqnump
  !$$$c          if(sum(abs(qp(:,iq)-qin))<1d-8 ) then
  !$$$        ngv  = ngp(iq)
  !$$$        ngvec(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
  !$$$        return
  !$$$c          endif
  !$$$c        enddo
  !$$$      elseif(ifi==2) then
  !$$$c        print *, "readqg cou: qin=",qin
  !$$$        call iqindx2(qin,ginv, qc,nqnumc, iq,qu)
  !$$$c        do iq=1, nqnumc
  !$$$c          if(sum(abs(qc(:,iq)-qin))<1d-8 ) then
  !$$$        ngv  = ngc(iq)
  !$$$        ngvec(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
  !$$$        return
  !$$$c          endif
  !$$$c        enddo
  !$$$      endif
  !$$$      stop "readqg: can not find QGpsi or QPcou for given q"
  !$$$      end subroutine
  !$$$
  !$$$
  !$$$c--------------------------
  !$$$      subroutine readqg0(key,qin,ginv,  qu,ngv)
  !$$$c- Get ngv
  !$$$cr key=='QGcou' or 'QGpsi'
  !$$$      implicit none
  !$$$      integer(4):: ngv, ifi, iq,verbose
  !$$$      real(8):: qin(3),qu(3),ginv(3,3)
  !$$$      character*(*) key
  !$$$      if    (key=='QGpsi') then
  !$$$        ifi=1
  !$$$        if(verbose()>=80) write (6,"('readqg0 psi: qin=',3f8.3,i5)") qin
  !$$$      elseif(key=='QGcou') then
  !$$$        ifi=2
  !$$$        if(verbose()>=80) write (6,"('readqg0 cou: qin=',3f8.3,i5)") qin
  !$$$      else
  !$$$        stop "readqg: wrongkey"
  !$$$      endif
  !$$$      if(init(ifi)) then
  !$$$        call init_readqg(ifi)
  !$$$        init(ifi)=.false.
  !$$$      endif
  !$$$      if(ifi==1) then
  !$$$        call iqindx2(qin,ginv, qp,nqnump, iq,qu)
  !$$$        ngv  = ngp(iq)
  !$$$      elseif(ifi==2) then
  !$$$        call iqindx2(qin,ginv, qc,nqnumc, iq,qu)
  !$$$        ngv  = ngc(iq)
  !$$$      endif
  !$$$      return
  !$$$      stop "readqg0: can not find QGpsi or QPcou for given q"
  !$$$      end subroutine
  !$$$
  !$$$c--------------------------
  !$$$      subroutine init_readqg(ifi)
  !$$$c- initialization. readin QGpsi or QGcou.
  !$$$      implicit none
  !$$$      integer(4):: ifi,ifiqg,iq,verbose
  !$$$      real(8)::qq(3)
  !$$$      print *,' init_readqg ifi=',ifi
  !$$$      ifiqg=4052
  !$$$      if(ifi==1) then
  !$$$        open(ifiqg, file='QGpsi',form='unformatted')
  !$$$        read(ifiqg) nqnump, ngpmx, QpGcut_psi
  !$$$        if(verbose()>49) write (6,"('init_readqg ngnumc ngcmx QpGcut_psi=',2i5,f8.3)")
  !$$$     &     nqnump, ngpmx, QpGcut_psi
  !$$$        allocate(ngvecp(3,ngpmx,nqnump),qp(3,nqnump),ngp(nqnump))
  !$$$        do iq=1, nqnump
  !$$$          read (ifiqg) qp(1:3,iq), ngp(iq)
  !$$$          read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
  !$$$          write (6,"('init_readqg psi qp ngp =',3f8.3,i5)") qp(1:3,iq),ngp(iq)
  !$$$        enddo
  !$$$      elseif(ifi==2) then
  !$$$        open(ifiqg, file='QGcou',form='unformatted')
  !$$$        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
  !$$$c         write (6,"('init_readqg ngnumc ngcmx QpGcut_cou=',2i5,f8.3)")
  !$$$c     &     nqnumc, ngcmx, QpGcut_cou
  !$$$        allocate(ngvecc(3,ngcmx,nqnumc),qc(3,nqnumc),ngc(nqnumc))
  !$$$        do iq=1, nqnumc
  !$$$          read(ifiqg) qc(1:3,iq), ngc(iq)
  !$$$c           write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
  !$$$          read (ifiqg) ngvecc(1:3,1:ngc(iq),iq)
  !$$$        enddo
  !$$$      endif
  !$$$      close(ifiqg)
  !$$$      end subroutine
  !$$$
  !$$$c--- release to save memory area.
  !$$$      subroutine releaseqg_notusednow(key)
  !$$$      implicit none
  !$$$      character*(*) key
  !$$$      integer(4):: ifi
  !$$$      if    (key=='QGpsi') then
  !$$$        ifi=1
  !$$$        deallocate(qp,ngvecp)
  !$$$      elseif(key=='QGcou') then
  !$$$        ifi=2
  !$$$        deallocate(qc,ngvecc)
  !$$$      else
  !$$$        stop "releaseqg: in readQGcou"
  !$$$      endif
  !$$$      init(ifi)=.false.
  !$$$      end subroutine
  !$$$
  !$$$      end module
  
  