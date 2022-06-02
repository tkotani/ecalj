c$$$      module m_readvmem
c$$$c------ a virtual memory system -----
c$$$c all these are private.
c$$$      implicit none
c$$$      real(8),allocatable,private:: idat(:,:),qp(:,:)
c$$$      logical,private:: init(2)=.true.
c$$$      real(8),private:: QpGcut_cou, QpGcut_psi
c$$$      integer(4),private::   nqnumc,nqnump,ngcmx,ngpmx
c$$$      integer(4),allocatable,private:: ngvecp(:,:,:),ngp(:),ngvecc(:,:,:),ngc(:)
c$$$c-----------------------------------------------------------------
c$$$      contains
c$$$c----------------------------------
c$$$      subroutine read_vmem(key,ngmx)
c$$$c- get ngcmx or mgpmx
c$$$      implicit none
c$$$      integer(4):: ngmx,ifiqg=4052
c$$$      character*(*) key
c$$$      if    (key=='QGpsi') then
c$$$        open(ifiqg, file='QGpsi',form='unformatted')
c$$$        read(ifiqg) nqnump, ngpmx, QpGcut_psi
c$$$        ngmx=ngpmx
c$$$      elseif(key=='QGcou') then
c$$$        open(ifiqg, file='QGcou',form='unformatted')
c$$$        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
c$$$        ngmx=ngcmx
c$$$      else
c$$$        stop "readngmx: key is not QGpsi QGcou"
c$$$      endif
c$$$      close(ifiqg)
c$$$      end subroutine
c$$$
c$$$c--------------------------
c$$$      subroutine readqg(key,qin,ginv,  qu,ngv,ngvec)
c$$$c- Get ngv and ngvec(3,ngv) for given qin(3)
c$$$cr key=='QGcou' or 'QGpsi'
c$$$      implicit none
c$$$      integer(4):: ngv, ngvec(3,*), ifi, iq,verbose
c$$$      real(8):: qin(3),qu(3),ginv(3,3)
c$$$      character*(*) key
c$$$      if    (key=='QGpsi') then
c$$$        ifi=1
c$$$        if(verbose()>=80) write (6,"('readqg psi: qin=',3f8.3,i5)") qin
c$$$      elseif(key=='QGcou') then
c$$$        ifi=2
c$$$        if(verbose()>=80) write (6,"('readqg cou: qin=',3f8.3,i5)") qin
c$$$      else
c$$$        stop "readqg: wrongkey"
c$$$      endif
c$$$
c$$$      if(init(ifi)) then
c$$$        call init_readqg(ifi)
c$$$        init(ifi)=.false.
c$$$      endif
c$$$
c$$$
c$$$      if(ifi==1) then
c$$$        call iqindx2(qin,ginv, qp,nqnump, iq,qu)
c$$$c        do iq=1, nqnump
c$$$c          if(sum(abs(qp(:,iq)-qin))<1d-8 ) then
c$$$        ngv  = ngp(iq)
c$$$        ngvec(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
c$$$        return
c$$$c          endif
c$$$c        enddo
c$$$      elseif(ifi==2) then
c$$$c        print *, "readqg cou: qin=",qin
c$$$        call iqindx2(qin,ginv, qc,nqnumc, iq,qu)
c$$$c        do iq=1, nqnumc
c$$$c          if(sum(abs(qc(:,iq)-qin))<1d-8 ) then
c$$$        ngv  = ngc(iq)
c$$$        ngvec(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
c$$$        return
c$$$c          endif
c$$$c        enddo
c$$$      endif
c$$$      stop "readqg: can not find QGpsi or QPcou for given q"
c$$$      end subroutine
c$$$
c$$$
c$$$c--------------------------
c$$$      subroutine readqg0(key,qin,ginv,  qu,ngv)
c$$$c- Get ngv
c$$$cr key=='QGcou' or 'QGpsi'
c$$$      implicit none
c$$$      integer(4):: ngv, ifi, iq,verbose
c$$$      real(8):: qin(3),qu(3),ginv(3,3)
c$$$      character*(*) key
c$$$      if    (key=='QGpsi') then
c$$$        ifi=1
c$$$        if(verbose()>=80) write (6,"('readqg0 psi: qin=',3f8.3,i5)") qin
c$$$      elseif(key=='QGcou') then
c$$$        ifi=2
c$$$        if(verbose()>=80) write (6,"('readqg0 cou: qin=',3f8.3,i5)") qin
c$$$      else
c$$$        stop "readqg: wrongkey"
c$$$      endif
c$$$      if(init(ifi)) then
c$$$        call init_readqg(ifi)
c$$$        init(ifi)=.false.
c$$$      endif
c$$$      if(ifi==1) then
c$$$        call iqindx2(qin,ginv, qp,nqnump, iq,qu)
c$$$        ngv  = ngp(iq)
c$$$      elseif(ifi==2) then
c$$$        call iqindx2(qin,ginv, qc,nqnumc, iq,qu)
c$$$        ngv  = ngc(iq)
c$$$      endif
c$$$      return
c$$$      stop "readqg0: can not find QGpsi or QPcou for given q"
c$$$      end subroutine
c$$$
c$$$c--------------------------
c$$$      subroutine init_readqg(ifi)
c$$$c- initialization. readin QGpsi or QGcou.
c$$$      implicit none
c$$$      integer(4):: ifi,ifiqg,iq,verbose
c$$$      real(8)::qq(3)
c$$$      print *,' init_readqg ifi=',ifi
c$$$      ifiqg=4052
c$$$      if(ifi==1) then
c$$$        open(ifiqg, file='QGpsi',form='unformatted')
c$$$        read(ifiqg) nqnump, ngpmx, QpGcut_psi
c$$$        if(verbose()>49) write (6,"('init_readqg ngnumc ngcmx QpGcut_psi=',2i5,f8.3)") 
c$$$     &     nqnump, ngpmx, QpGcut_psi
c$$$        allocate(ngvecp(3,ngpmx,nqnump),qp(3,nqnump),ngp(nqnump))
c$$$        do iq=1, nqnump
c$$$          read (ifiqg) qp(1:3,iq), ngp(iq)
c$$$          read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
c$$$          write (6,"('init_readqg psi qp ngp =',3f8.3,i5)") qp(1:3,iq),ngp(iq)
c$$$        enddo
c$$$      elseif(ifi==2) then
c$$$        open(ifiqg, file='QGcou',form='unformatted')
c$$$        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
c$$$c         write (6,"('init_readqg ngnumc ngcmx QpGcut_cou=',2i5,f8.3)")
c$$$c     &     nqnumc, ngcmx, QpGcut_cou
c$$$        allocate(ngvecc(3,ngcmx,nqnumc),qc(3,nqnumc),ngc(nqnumc))
c$$$        do iq=1, nqnumc
c$$$          read(ifiqg) qc(1:3,iq), ngc(iq)
c$$$c           write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
c$$$          read (ifiqg) ngvecc(1:3,1:ngc(iq),iq)
c$$$        enddo
c$$$      endif
c$$$      close(ifiqg)
c$$$      end subroutine
c$$$
c$$$c--- release to save memory area.
c$$$      subroutine releaseqg_notusednow(key)
c$$$      implicit none
c$$$      character*(*) key
c$$$      integer(4):: ifi
c$$$      if    (key=='QGpsi') then
c$$$        ifi=1
c$$$        deallocate(qp,ngvecp)
c$$$      elseif(key=='QGcou') then
c$$$        ifi=2
c$$$        deallocate(qc,ngvecc)
c$$$      else
c$$$        stop "releaseqg: in readQGcou"
c$$$      endif
c$$$      init(ifi)=.false.
c$$$      end subroutine
c$$$
c$$$      end module

