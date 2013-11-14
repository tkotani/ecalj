      module m_read_ppovl
!! == read PPOVLGG,PPOVLG,PPOVLI ==
!!  ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
      implicit none
      integer(4):: ngc2
      complex(8),allocatable :: ppx(:,:)
      integer(4),allocatable :: ngvecc2(:,:)
      complex(8),allocatable:: ggg(:),ppovlinv(:,:)
      integer,allocatable:: nvggg(:,:),nvgcgp2(:,:),ngvecc(:,:)
      integer:: nggg,ngcgp,ngcread
!! private
      integer(4),private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx, ngcgpmx,nqini,nqnumt
      logical,private:: ppovlclosed=.true.,init=.true.
      integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
      real(8),allocatable,private    :: qx_s(:,:)
      complex(8),allocatable,private :: ppx_s(:,:,:)
      logical,private:: debug=.false.
c-----------------------------------------------------
      contains
c-----------------------------------------------------
      subroutine getppx2(qbas,qi) 
!! This return nvggg,nvgcgp2,ngvecc,  nggg,ngcgp,ngcread, ggg,ppovlinv
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovlg = 9221,ippovli = 9222,idxk
      real(8)::qbas(3,3),qi(3),qx(3)
      complex(8),allocatable :: ppxs(:,:,:)
      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:)
      real(8),allocatable :: qxs(:,:)
      logical :: KeepPPOVL,debug=.false., ippovlggooo=.true.
      integer:: verbose,io,ippovlgg=9223,ngcread2
      character*3 charnum3
!! because of get_ngpgp iqix can be not consistent among PPOVLG and PPOVLI --> use iqix=0
      if(verbose()>=100) debug=.true.
      if(ippovlggooo) then
         open(ippovlgg,file= "PPOVLGG",form='unformatted')
         read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
         write(6,"('Readin getppx2: nggg ngcgp nqq=',3i5)") nggg, ngcgp, nqq
         allocate(nvggg(1:3,1:nggg),ggg(1:nggg),nvgcgp2(1:3,ngcgp))
         read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
         read(ippovlgg) nvggg(1:3,1:nggg)
         read(ippovlgg) ggg(1:nggg)
         close(ippovlgg)
         ippovlggooo=.false.
      endif
      do iqi = nqini,nqnumt
         open(ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
         read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
         if(debug) print *,'getppx2: iqi qi=',iqi,qi,qx
         ngc = ngcread
         if(sum(abs(qx-qi))>1d-10) then
           close(ippovlg)
         else
           open(ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
           read(ippovli) qx, ngcread2 !, ngcx_s(iqi),ngc2_s(iqi)
           if(ngcread/=ngcread2) call rx('rppovl.F: inconsistent PPOVLI PPOVLg')
           if(ngc==0) exit
           if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
           allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
           read(ippovlg) ngvecc(1:3,1:ngc)
           read(ippovli) ppovlinv(1:ngc,1:ngc)
           close(ippovlg)
           close(ippovli)
           goto 1012
         endif   
      enddo
      write(6,"('qi qx=',3f13.5,3x,3f13.5)") qi,qx
      call rx('getppx2: can not find given qi')
c$$$
c$$$      iqix=0
c$$$      do io=1,2 !The 1st path start from saved iqix. The second path(io=2) start from iqix=0.
c$$$        if(iqix<=0) then
c$$$           if(verbose()>=100) print *,'iqix <0 readin PPOVLG'
c$$$           if(iqix==0) close(ippovlg)
c$$$           if(iqix==0) close(ippovli)
c$$$           open(ippovlg,file= "PPOVLG",form='unformatted')
c$$$           open(ippovli,file= "PPOVLI",form='unformatted')
c$$$           iqix=0
c$$$        endif
c$$$        do iqi = iqix+1,nqq
c$$$           if(debug) print *,'getppx2: iqi qi=',iqi,qi
c$$$           read(ippovli) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
c$$$           read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
c$$$           if(debug) print *,'getppx2: iqi qx=',iqi,qx
c$$$           ngc=ngcread
c$$$           if(sum(abs(qx-qi))<1d-10) then !hit
c$$$              if(ngc/=0) then
c$$$               if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
c$$$               allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
c$$$               read(ippovlg) ngvecc(1:3,1:ngc)
c$$$               read(ippovli) ppovlinv(1:ngc,1:ngc)
c$$$              endif 
c$$$              iqix=iqi+1
c$$$              goto 1012
c$$$           else
c$$$              read(ippovlg)
c$$$              read(ippovli)
c$$$           endif  
c$$$        enddo
c$$$        if(iqix==0) exit
c$$$        iqix=0
c$$$      enddo
c$$$      call rx('getppx2: can not find given qi')
 1012 continue
c      print *,'eeeeeeee111111 iqix=',iqix
!! rotation of gvector rotgvec is moved to melpl2t nov2013
      end subroutine getppx2
c$$$c-----------------------------------------------------
c$$$      subroutine get_nvgcgp(qbas,qi,nadd) 
c$$$!! This return ngcgp,nvgcgp(1:3,ngcgp)
c$$$!! ngvcgp2 is for |qi+nadd+G|< QpGcutcou + QpGcutpsi (triangle inequality).
c$$$!!   Readin ngvgcp2 specify G'for |qi+G'|< QpGcutcou + QpGcutpsi
c$$$!!   Thus G= G'-nadd
c$$$      implicit none
c$$$      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovlg = 9221,ippovli = 9222,idxk
c$$$      real(8)::qbas(3,3),qi(3),qx(3)
c$$$      integer:: verbose,io, ndummy1,ndummy2,ndmmy3,ndummy4,nadd(3),igc
c$$$cccccccccccccccccccccccccccc
c$$$c      real(8):: rlatp(3,3),xmx2(3),ddd(3),dddx(3),qlati(3,3),vol
c$$$c      integer,parameter:: noutmx=48
c$$$c      integer:: nlatout(3,noutmx),nout,nimnx(3)
c$$$c      call shortn3_initialize(qbas,rlatp,xmx2)
c$$$      iqix=0
c$$$cccccccccccccccccccccccccccc
c$$$      if(allocated(nvgcgp2)) deallocate(nvgcgp2)
c$$$      print *,'xxxxxxxxxxxxxxxxxxxxxx init xxxxxxxxxxxx ppppppppppp'
c$$$      do io=1,2 !The 1st path start from saved iqix. The second path(io=2)start from iqix=0.
c$$$        if(iqix<=0) then
c$$$           if(iqix==0) close(ippovlg)
c$$$           open(ippovlg,file= "PPOVLG",form='unformatted')
c$$$           read(ippovlg) ndummy1,ndummy2,ndmmy3
c$$$           iqix=0
c$$$        endif
c$$$        do iqi = iqix+1,nqq
c$$$           read(ippovlg) qx, ndummy4, ngcgp
c$$$           print *,' qi=',qi
c$$$           print *,' qx=',qx
c$$$ccccccccccccccccccccccc
c$$$c           call dinv33(qbas,0,qlati,vol)
c$$$c           ddd = matmul(qlati,qi-qx)
c$$$c           call shortn3(rlatp,xmx2,ddd,noutmx, nout,nlatout) !need to call shortn3_initialize
c$$$c           dddx= matmul(qbas(:,:), ddd+nlatout(:,1))
c$$$cc           if(abs(sum(abs(nlatout(:,1))))==0) then
c$$$c           nimnx = nlatout(:,1)
c$$$c           write(6,"('pppppp ',f10.4,3x,3f13.6,3i3)")sum(abs(dddx**2)), dddx, nimnx
c$$$c           endif
c$$$ccccccccccccccccccccccc
c$$$           if(sum(abs(qx-qi))<1d-10) then !hit
c$$$               allocate(nvgcgp2(1:3,1:ngcgp))
c$$$               read(ippovlg) nvgcgp2(1:3,1:ngcgp)
c$$$               read(ippovlg) 
c$$$               read(ippovlg) 
c$$$               iqix=iqi+1
c$$$               goto 1112
c$$$           else
c$$$              read(ippovlg)
c$$$              read(ippovlg) 
c$$$              read(ippovlg) 
c$$$           endif  
c$$$        enddo
c$$$        if(iqix==0) exit
c$$$        iqix=0
c$$$        print *,'eeeeeee io=',io
c$$$      enddo
c$$$      call rx('get_nvgcgp: can not find given qi')
c$$$ 1112 continue 
c$$$      do igc=1,ngcgp
c$$$         nvgcgp2(1:3,igc)= nvgcgp2(1:3,igc) - nadd
c$$$      enddo   
c$$$      print *,'eeeeeeee2222 iqix=',iqix
c$$$!! rotation is moved to melpl2 nov2013
c$$$      end subroutine

c$$$!! ----------------------------------------------------------------------------
c$$$      subroutine getppx(qbas,  ngc,qi,symope) ! this return ppx and ngvecc2
c$$$      implicit none
c$$$      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovl = 9221,idxk
c$$$      real(8)::qbas(3,3),qi(3),symope(3,3),qx(3)
c$$$      complex(8),allocatable :: ppxs(:,:,:)
c$$$      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
c$$$     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:)
c$$$      real(8),allocatable :: qxs(:,:)
c$$$      logical :: KeepPPOVL
c$$$      integer:: verbose
c$$$      if(debug.and.allocated(ngvecc2_0_s)) then
c$$$      print *,'getppx :sumcheck ngvecc2_0= 3',ngc2_s(3),sum(abs(ngvecc2_0_s(1:3,1:ngc2_s(3),3)))
c$$$      endif
c$$$c--- KeepPPOVL=F to reduce memory usage.
c$$$      if(allocated(ppx)    ) deallocate(ppx)
c$$$      if(allocated(ngvecc2)) deallocate(ngvecc2)
c$$$      ngc2=999999
c$$$      if(init) then
c$$$        if(KeepPPOVL()      ) print *,' KeepPPOVL=T; keep PPOVL in m_read_ppovl'
c$$$        if(.not.KeepPPOVL())  print *,' KeepPPOVL=F; not keep PPOVL in m_read_ppovl'
c$$$        init=.false.
c$$$      endif
c$$$
c$$$      if(.not.KeepPPOVL()) then
c$$$        do
c$$$          if(ppovlclosed) then
c$$$            open(ippovl,file= "PPOVL",form='unformatted')
c$$$            read(ippovl) ngcmx, ngc2mx, nqq
c$$$            ppovlclosed =.false.
c$$$c            print *,' ppovl open =',ngcmx, ngc2mx, nqq
c$$$            iqix=0
c$$$          endif
c$$$c          print *,' getppx: ppovlclosed iqix nqq=',ppovlclosed, iqix,nqq
c$$$          do iqi = iqix+1,nqq
c$$$c            print *, ' iqi=',iqi
c$$$            read(ippovl) qx, ngcx,ngc2
c$$$c            print *, ' qx ngcx ngc2=', qx, ngcx,ngc2
c$$$            if(ngcx/=0) then !Aug2005
c$$$              allocate(ngvecc2_0(1:3,1:ngc2), ppx(ngcx,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
c$$$              read(ippovl) ngvecc2_0(1:3,1:ngc2),ppx(1:ngcx,1:ngc2)
c$$$            endif
c$$$            if(sum(abs(qx-qi))<1d-10) then
c$$$              ngc=ngcx
c$$$c              if(ngc/=ngcx) then
c$$$c                print *,iqi, qx, qi, ngcx, ngc
c$$$c                stop 'getppx: ngc/= ngc from PPOVL'
c$$$c              endif
c$$$              iqix=iqi
c$$$              goto 1110
c$$$            endif
c$$$            if(ngcx/=0) deallocate(ngvecc2_0,ppx)         !sf..5June02
c$$$          enddo
c$$$Cstop2rx 2013.08.09 kino          if(iqix==0) stop ' getppx: Wrong PPVOL?'
c$$$          if(iqix==0) call rx( ' getppx: Wrong PPVOL?')
c$$$          close(ippovl)
c$$$          ppovlclosed =.true.
c$$$        enddo
c$$$      else
c$$$        if(ppovlclosed) then
c$$$          open(ippovl,file= "PPOVL",form='unformatted')
c$$$          read(ippovl) ngcmx, ngc2mx, nqq
c$$$          if(verbose()>=100) print *,' ppovl open =',ngcmx, ngc2mx, nqq
c$$$          allocate(qx_s(3,nqq), ngcx_s(nqq),ngc2_s(nqq)
c$$$     &              ,ngvecc2_0_s(3,ngc2mx,nqq),ppx_s(ngcmx,ngc2mx,nqq))
c$$$          do iqi = 1,nqq
c$$$            if(verbose()>=100)  print *, ' iqi=',iqi
c$$$            read(ippovl) qx_s(:,iqi), ngcx_s(iqi),ngc2_s(iqi)
c$$$            ngcx = ngcx_s(iqi)
c$$$            ngc2 = ngc2_s(iqi)
c$$$            if(verbose()>=100) print *,' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
c$$$cccccccccccccccccccccccccccc
c$$$c           print *, ' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
c$$$cccccccccccccccccccccccccccccc
c$$$            if(ngcx/=0) then !Aug2005
c$$$              read(ippovl) ngvecc2_0_s(1:3,1:ngc2,iqi),ppx_s(1:ngcx,1:ngc2,iqi)
c$$$            endif
c$$$          enddo
c$$$          ppovlclosed =.false.
c$$$          if(verbose()>=100) print *,' end of readin PPOVL'
c$$$        endif
c$$$        if(verbose()>=100) print *,' goto iqi search'
c$$$ccccccccccccccccccccccccccccccccc
c$$$c           print *,' qi=',qi,nqq
c$$$ccccccccccccccccccccccccccccccccc
c$$$        iqi = idxk(qi,qx_s,nqq)
c$$$        ngc2= ngc2_s(iqi)
c$$$        if(verbose()>=100) print *,' iqi ngc2',iqi,ngc2
c$$$        ngc=ngcx_s(iqi)
c$$$        allocate(ngvecc2_0(3,ngc2), ppx(ngc,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
c$$$        ngvecc2_0(1:3,1:ngc2)= ngvecc2_0_s(1:3,1:ngc2,iqi)
c$$$ccccccccccccccccccccccccccccc
c$$$ccccccccccccccccccccccccccccc
c$$$        ppx(1:ngc,1:ngc2)   = ppx_s(1:ngc,1:ngc2,iqi)
c$$$c          if(ngc/=ngcx_s(iqi)) then
c$$$c            print *,iqi, qx, qi, ngcx, ngc
c$$$c            stop 'getppx: xxx ngc/= ngc from PPOVL'
c$$$c          endif
c$$$      endif
c$$$c---------------------------------------------
c$$$ 1110 continue
c$$$      allocate(ngvecc2(1:3,1:ngc2),ngc2s(1))
c$$$      ngc2s(1) = ngc2
c$$$      call rotgvec(symope, 1,
c$$$c     i    ngc2mx, ngc2s,qbas, ngvecc2_0, !alpha linux caused segmentation error
c$$$     i    ngc2, ngc2s,qbas, ngvecc2_0, 
c$$$     o    ngvecc2)
c$$$      if(allocated(ngc2s)    ) deallocate(ngc2s)
c$$$      if(allocated(ngvecc2_0)) deallocate(ngvecc2_0)
c$$$      end subroutine
      end module
