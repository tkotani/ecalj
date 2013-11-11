      module m_read_ppovl
!! == read PPOVLG ==
c   ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
      implicit none
      integer(4):: ngc2
      complex(8),allocatable :: ppx(:,:)
      integer(4),allocatable :: ngvecc2(:,:)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex(8),allocatable:: ggg(:),ppovlinv(:,:)
      integer,allocatable:: nvggg(:,:),nvgcgp2(:,:),ngvecc(:,:)
      integer:: nggg,ngcgp,ngcread
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! private
      integer(4),private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx, ngcgpmx
      logical,private:: ppovlclosed=.true.,init=.true.
      integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
      real(8),allocatable,private    :: qx_s(:,:)
      complex(8),allocatable,private :: ppx_s(:,:,:)
      logical,private:: debug=.false.
c-----------------------------------------------------
      contains
c-----------------------------------------------------
      subroutine get_nvgcgp(qbas,qi,nadd) 
!! This return ngcgp,nvgcgp(1:3,ngcgp)
!! ngvcgp2 is for |qi+nadd+G|< QpGcutcou + QpGcutpsi (triangle inequality).
!!   Readin ngvgcp2 specify G'for |qi+G'|< QpGcutcou + QpGcutpsi
!!   Thus G= G'-nadd
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovlg = 9221,ippovli = 9222,idxk
      real(8)::qbas(3,3),qi(3),qx(3)
      integer:: verbose,io, ndummy1,ndummy2,ndmmy3,ndummy4,nadd(3),igc
cccccccccccccccccccccccccccc
      iqix=0
cccccccccccccccccccccccccccc
      if(allocated(nvgcgp2)) deallocate(nvgcgp2)
      do io=1,2 !The 1st path start from saved iqix. The second path(io=2)start from iqix=0.
        if(iqix<=0) then
           if(iqix==0) close(ippovlg)
           open(ippovlg,file= "PPOVLG",form='unformatted')
           read(ippovlg) ndummy1,ndummy2,ndmmy3
           iqix=0
        endif
        do iqi = iqix+1,nqq
           read(ippovlg) qx, ndummy4, ngcgp
           if(sum(abs(qx-qi))<1d-10) then !hit
               allocate(nvgcgp2(1:3,1:ngcgp))
               read(ippovlg) nvgcgp2(1:3,1:ngcgp)
               read(ippovlg) 
               read(ippovlg) 
               iqix=iqi+1
               goto 1112
           else
              read(ippovlg)
              read(ippovlg) 
              read(ippovlg) 
           endif  
        enddo
        if(iqix==0) exit
        iqix=0
        print *,'eeeeeee io=',io
      enddo
      call rx('get_nvgcgp: can not find given qi')
 1112 continue 
      do igc=1,ngcgp
         nvgcgp2(1:3,igc)= nvgcgp2(1:3,igc) - nadd
      enddo   
      print *,'eeeeeeee2222 iqix=',iqix
!! rotation is moved to melpl2 oct2013
      end subroutine
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
      logical :: KeepPPOVL,debug=.false.
      integer:: verbose,io
!! because of get_ngpgp iqix can be not consistent among PPOVLG and PPOVLI --> use iqix=0
      iqix=0
      if(verbose()>=100) debug=.true.
      if(allocated(nvggg)) deallocate(nvggg,ngvecc,ggg,ppovlinv)
      do io=1,2 !The 1st path start from saved iqix. The second path(io=2)start from iqix=0.
        if(iqix<=0) then
           if(iqix==0) close(ippovlg)
           if(iqix==0) close(ippovli)
           open(ippovlg,file= "PPOVLG",form='unformatted')
           open(ippovli,file= "PPOVLI",form='unformatted')
           read(ippovlg) ngggmx, ngcgpmx, nqq
           iqix=0
        endif
        do iqi = iqix+1,nqq
           if(debug) print *,'getppx2: iqi qi=',iqi,qi
           read(ippovli) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
           read(ippovlg) qx, nggg, ngcgp, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
           if(debug) print *,'getppx2: iqi qx=',iqi,qx
           ngc=ngcread
           if(sum(abs(qx-qi))<1d-10) then !hit
              if(ngc/=0) then
               allocate(nvggg(1:3,1:nggg),ggg(1:nggg),ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
               read(ippovlg) !nvgcgp2(1:3,1:ngcgp) skip. This is by getngcgp2
               read(ippovlg) nvggg(1:3,1:nggg),ngvecc(1:3,1:ngc)
               read(ippovlg) ggg(1:nggg)
               read(ippovli) ppovlinv(1:ngc,1:ngc)
              endif 
              iqix=iqi+1
              goto 1012
           else
              read(ippovlg)
              read(ippovlg)
              read(ippovlg)
              read(ippovli)
           endif  
        enddo
        if(iqix==0) exit
        iqix=0
      enddo
      call rx('getppx2: can not find given qi')
 1012 continue
      print *,'eeeeeeee111111 iqix=',iqix
!! rotation is moved to melpl2 oct2013
      end subroutine getppx2


!! ----------------------------------------------------------------------------
      subroutine getppx(qbas,  ngc,qi,symope) ! this return ppx and ngvecc2
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovl = 9221,idxk
      real(8)::qbas(3,3),qi(3),symope(3,3),qx(3)
      complex(8),allocatable :: ppxs(:,:,:)
      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:)
      real(8),allocatable :: qxs(:,:)
      logical :: KeepPPOVL
      integer:: verbose
cccccccccccccccccccccccccccccccccccccccccc
      if(debug .and. allocated(ngvecc2_0_s)) then
        print *,' getppx :sumcheck ngvecc2_0= 3',ngc2_s(3),sum(abs(ngvecc2_0_s(1:3,1:ngc2_s(3),3)))
      endif
ccccccccccccccccccccccccccccccccccccccccc
c--- KeepPPOVL=F to reduce memory usage.
      if(allocated(ppx)    ) deallocate(ppx)
      if(allocated(ngvecc2)) deallocate(ngvecc2)
      ngc2=999999
      if(init) then
        if(KeepPPOVL()      ) print *,' KeepPPOVL=T; keep PPOVL in m_read_ppovl'
        if(.not.KeepPPOVL())  print *,' KeepPPOVL=F; not keep PPOVL in m_read_ppovl'
        init=.false.
      endif

      if(.not.KeepPPOVL()) then
        do
          if(ppovlclosed) then
            open(ippovl,file= "PPOVL",form='unformatted')
            read(ippovl) ngcmx, ngc2mx, nqq
            ppovlclosed =.false.
c            print *,' ppovl open =',ngcmx, ngc2mx, nqq
            iqix=0
          endif
c          print *,' getppx: ppovlclosed iqix nqq=',ppovlclosed, iqix,nqq
          do iqi = iqix+1,nqq
c            print *, ' iqi=',iqi
            read(ippovl) qx, ngcx,ngc2
c            print *, ' qx ngcx ngc2=', qx, ngcx,ngc2
            if(ngcx/=0) then !Aug2005
              allocate(ngvecc2_0(1:3,1:ngc2), ppx(ngcx,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
              read(ippovl) ngvecc2_0(1:3,1:ngc2),ppx(1:ngcx,1:ngc2)
            endif
            if(sum(abs(qx-qi))<1d-10) then
              ngc=ngcx
c              if(ngc/=ngcx) then
c                print *,iqi, qx, qi, ngcx, ngc
c                stop 'getppx: ngc/= ngc from PPOVL'
c              endif
              iqix=iqi
              goto 1110
            endif
            if(ngcx/=0) deallocate(ngvecc2_0,ppx)         !sf..5June02
          enddo
Cstop2rx 2013.08.09 kino          if(iqix==0) stop ' getppx: Wrong PPVOL?'
          if(iqix==0) call rx( ' getppx: Wrong PPVOL?')
          close(ippovl)
          ppovlclosed =.true.
        enddo
      else
        if(ppovlclosed) then
          open(ippovl,file= "PPOVL",form='unformatted')
          read(ippovl) ngcmx, ngc2mx, nqq
          if(verbose()>=100) print *,' ppovl open =',ngcmx, ngc2mx, nqq
          allocate(qx_s(3,nqq), ngcx_s(nqq),ngc2_s(nqq)
     &              ,ngvecc2_0_s(3,ngc2mx,nqq),ppx_s(ngcmx,ngc2mx,nqq))
          do iqi = 1,nqq
            if(verbose()>=100)  print *, ' iqi=',iqi
            read(ippovl) qx_s(:,iqi), ngcx_s(iqi),ngc2_s(iqi)
            ngcx = ngcx_s(iqi)
            ngc2 = ngc2_s(iqi)
            if(verbose()>=100) print *,' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
cccccccccccccccccccccccccccc
c           print *, ' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
cccccccccccccccccccccccccccccc
            if(ngcx/=0) then !Aug2005
              read(ippovl) ngvecc2_0_s(1:3,1:ngc2,iqi),ppx_s(1:ngcx,1:ngc2,iqi)
            endif
          enddo
          ppovlclosed =.false.
          if(verbose()>=100) print *,' end of readin PPOVL'
        endif
        if(verbose()>=100) print *,' goto iqi search'
ccccccccccccccccccccccccccccccccc
c           print *,' qi=',qi,nqq
ccccccccccccccccccccccccccccccccc
        iqi = idxk(qi,qx_s,nqq)
        ngc2= ngc2_s(iqi)
        if(verbose()>=100) print *,' iqi ngc2',iqi,ngc2
        ngc=ngcx_s(iqi)
        allocate(ngvecc2_0(3,ngc2), ppx(ngc,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
        ngvecc2_0(1:3,1:ngc2)= ngvecc2_0_s(1:3,1:ngc2,iqi)
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
        ppx(1:ngc,1:ngc2)   = ppx_s(1:ngc,1:ngc2,iqi)
c          if(ngc/=ngcx_s(iqi)) then
c            print *,iqi, qx, qi, ngcx, ngc
c            stop 'getppx: xxx ngc/= ngc from PPOVL'
c          endif
      endif
c---------------------------------------------
 1110 continue
      allocate(ngvecc2(1:3,1:ngc2),ngc2s(1))
      ngc2s(1) = ngc2
      call rotgvec(symope, 1,
c     i    ngc2mx, ngc2s,qbas, ngvecc2_0, !alpha linux caused segmentation error
     i    ngc2, ngc2s,qbas, ngvecc2_0, 
     o    ngvecc2)
      if(allocated(ngc2s)    ) deallocate(ngc2s)
      if(allocated(ngvecc2_0)) deallocate(ngvecc2_0)
      end subroutine


      end module
