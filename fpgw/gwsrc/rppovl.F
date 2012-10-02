C ... read PPOVL
      module m_read_ppovl
c--------readin dataset from PPOVL ------------
c   ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
      implicit none
      integer(4):: ngc2
      complex(8),allocatable :: ppx(:,:)
      integer(4),allocatable :: ngvecc2(:,:)

c---private
      integer(4),private:: iqix, ippovl=0, ngcmx, ngc2mx, nqq
      logical,private:: ppovlclosed=.true.,init=.true.
      integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
      real(8),allocatable,private    :: qx_s(:,:)
      complex(8),allocatable,private :: ppx_s(:,:,:)
      logical,private:: debug=.false.
c-----------------------------------------------------------------
      contains
c-----------------------------------------------------
      subroutine getppx(qbas,  ngc,qi,symope) ! this return ppx and ngvecc2
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovl = 9221,idxk
      real(8)::qbas(3,3),qi(3),symope(3,3),qx(3)
      complex(8),allocatable :: ppxs(:,:,:),ppxread(:,:)
      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:),nread(:,:)
      real(8),allocatable :: qxs(:,:)
      logical :: KeepPPOVL
      integer:: verbose
c      print *,'getppx:'
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
          if(iqix==0) stop ' getppx: Wrong PPVOL?'
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
              print *,' reading ppx_s'
              allocate(nread(1:3,ngc2), ppxread(ngcx,ngc2))
              read(ippovl) nread, ppxread
              ngvecc2_0_s(1:3,1:ngc2,iqi)=nread
              ppx_s(1:ngcx,1:ngc2,iqi)=ppxread
              deallocate(nread,ppxread)
c              read(ippovl) ngvecc2_0_s(1:3,1:ngc2,iqi),ppx_s(1:ngcx,1:ngc2,iqi)
              
              print *,' reading ppx_s finished'
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
