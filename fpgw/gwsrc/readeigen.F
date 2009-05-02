C ... read eigen
      module m_readeigen
c------Return eigenvalus and eigenfunctions for given q and isp.
c Routines are:
C   readeval: ev (eigenvalue)
c
      implicit none
      integer(4),private :: nqtt
      real(8),allocatable,private:: qtt(:,:),evud(:,:,:),ginv(:,:)
      logical,private:: init=.true.,init2=.true.,keepeig
      integer(4),allocatable,private:: ngp(:)
      integer(4),private:: nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   nqnum , ngpmx, nnnn, nsp,ldim2,mrecg,ifcphi,ifgeig
      complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
      real(8),private:: leval
ctm
      integer(4),private :: nwf,ifcphi_mlw,ifgeig_mlw
      complex(8),allocatable,private:: geig_mlw(:,:,:,:),
     &           cphi_mlw(:,:,:,:)
c-----------------------------------------------------------------
      contains

c------------------------------
      subroutine readeval(q,isp, ev)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx
      real(8):: q(3)
      real(8):: ev(nband)
      if(init) then
        stop 'readeigen: modele is not initialized yet'
      endif
c      iq = idxk(q, qtt,nqtt)
      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' iq=',iq
      ev(1:nband) = evud(1:nband,iq,isp)
      end subroutine

c------------------------------
      subroutine readgeig(q,ngp_in,isp, qu,geigen)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nband)
      if(init2) then
        stop 'readgeig: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgeig: q=',q
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.

      if(ngp_in < ngp(iq)) then
        print *,'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        stop 'readgeig: ngpmx<ngp(iq)'
      endif
! for ngp=0 Aug2005
      if(ngp(iq)==0) return

      if(keepeig) then
        geigen(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeig, rec=ikpisp) geigen(1:ngpmx,1:nband)
      endif
      end subroutine

c------------------------------
      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2_in,nband) !takao 1->*->nband
      if(init2) then
        stop 'readcphi: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgcphi: q=',q
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.

      if(keepeig) then
        cphif(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi, rec=ikpisp) cphif(1:ldim2,1:nband)
      endif
      end subroutine

c------------------------------
      subroutine readgeig_mlw(q,ngp_in,isp, qu,geigen)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) then
        stop 'readgeig_mlw: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgeig: q=',q
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.

      if(ngp_in < ngp(iq)) then
        print *,'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        stop 'readgeig_mlw: ngpmx<ngp(iq)'
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geig_mlw(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeig_mlw, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine

c------------------------------
      subroutine readcphi_mlw(q,ldim2_in,isp,  qu,cphif)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2_in,nwf)
      if(init2) then
        stop 'readcphi_mlw: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgcphi_mlw: q=',q
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.

      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphi(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine

c------------------------------
      subroutine init_readeigen(ginv_in,nsp_in,nband_in,mrece_in)
!initialization. Save QpGpsi EVU EVD to arrays.
      implicit none
      integer(4):: ifev(2),iq,is,ifiqg,nnnn,ikp,nsp_in, 
     & nband_in,mrece_in,iopen,iclose,isx,mrecb_in,ik,ib,verbose
      real(8):: ginv_in(3,3)
      integer(4),allocatable:: ngvecp(:,:)
      real(8):: QpGcut_psi
      write(6,*) 'init_readeigen:'
      allocate(ginv(3,3))
      ginv  = ginv_in
      nsp   = nsp_in
      if(nsp<0.or.nsp>2) stop 'init_reaeigen:nsp wrong'
      nband = nband_in
      mrece = mrece_in
c      write(6,*) 'init_readeigen:1'
c... Readin QGpsi
      ifiqg  = 101
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqtt , ngpmx, QpGcut_psi, nnnn
      write(6,*)'read(ifiqg)', nqtt , ngpmx, QpGcut_psi, nnnn
      allocate( ngvecp(3,ngpmx),qtt(3,nqtt),ngp(nqtt) )
      do ikp = 1,nqtt
        read (ifiqg) qtt(:,ikp), ngp(ikp)
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp))
      enddo
      close(ifiqg)
      deallocate(ngvecp)
c      write(6,*) 'init_readeigen:2'

c... Readin ev
      allocate(evud(nband,nqtt,nsp))
      ifev(1) = iopen('EVU', 0,0,mrece)
      if(nsp==2) ifev(2) = iopen('EVD', 0,0,mrece)
      do is= 1,nsp
        call rwdd (ifev(is),   !read eigenvalues
     d              nband,nqtt, 
     o              evud(1,1,is))
c        print *, 'evud is=',is,minval(evud(:,:,is))
        if(verbose()>100 ) then
          do ik=1,nqtt
            do ib=1,nband
              write(6,"('ib ik e=',2i5,f13.5)") ib,ik,evud(ib,ik,is)
            enddo
          enddo
        endif
      enddo
c      stop
      leval= minval(evud)
      isx = iclose('EVU')
      if(nsp==2) isx = iclose('EVD')
      init=.false.
      end subroutine

      real(8) function lowesteval()
      lowesteval=leval
      end function
c-------------------------------------------------------------
      subroutine init_readeigen2(mrecb_in,ldim2_in,mrecg_in)
c this should be called after init_readgeigen
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,ldim2_in, 
     & iopen,iclose,isx,mrecb_in,ikpisp,mrecg_in,verbose
      logical :: keepeigen
      mrecb = mrecb_in
      ldim2 = ldim2_in
      mrecg = mrecg_in
      keepeig = keepeigen()

      init2=.false.
c --- Readin CPHI and GEIG
      if(Keepeig     ) print *,' KeepEigen=T; readin geig and cphi into m_readeigen'
      if(.not.Keepeig) print *,' KeepEigen=F; not keep geig and cphi in m_readeigen'

      ifgeig = iopen('GEIG',0,0,mrecg)
      ifcphi = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction
      print *,' ifgeig ifcphi=',ifgeig,ifcphi
      if(.not.keepeig) return
c      ifgeig = iopen('GEIG',0,-1,0)
c      ifcphi = iopen('CPHI',0,-1,0) ! Coeeficients of eigenfunction
      allocate(geig(ngpmx,nband,nqtt,nsp))
      allocate(cphi(ldim2,nband,nqtt,nsp))
      do ikp= 1,nqtt
        do is= 1,nsp
c        read(ifgeig) geig(1:ngpmx,1:nband,ikp,is) !geig(1:ngp(ikp),1:nband,ikp,1:is)
c        read(ifcphi) cphi(1:ldim2,1:nband,ikp,is)
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
ccccccccccccccccccccccc
          if(verbose()>40) write(6,'("Readin geig sumcheck ",i5,3f9.4,"  ",i5,2d13.6)')
     &      ikp,qtt(:,ikp),ngp(ikp),sum(geig(1:ngp(ikp),1:nband,ikp,is))
cccccccccccccccccccccc
c        read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is)
          read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nband,ikp,is)
        enddo
      enddo
      isx= iclose('CPHI')
      isx= iclose('GEIG')

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw(nwf_o)
c allocate cphi_mlw and geig_mlw
c for hx0fp0_mlw (ixc=112 mode)
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,ib,iwf,iwf2,iko_ix,iko_fx,in,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c allocate geig_mlw and cphi_mlw
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(dnk,uum)


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig_mlw(ngpmx,nwf,nqtt,nsp))
        allocate(cphi_mlw(ldim2,nwf,nqtt,nsp))
        geig_mlw = 0d0
        cphi_mlw = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig_mlw(:,iwf,ikp,is) = geig_mlw(:,iwf,ikp,is) + 
     &           geig(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi_mlw(:,iwf,ikp,is) = cphi_mlw(:,iwf,ikp,is) + 
     &           cphi(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

          enddo
        enddo


      else
        ifcphi_mlw = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_mlw = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_mlw, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_mlw, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
      endif


      deallocate(cbwf)


c replace dimensions
      nwf_o = nwf

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_noeval1D(
     &                 nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2,
     & iset,nset,iset2,nset2,ini
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw_noeval1D'


c --- Readin MLWU/D and UU1dU/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifuu   = iopen('UU1dU',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifuu   = iopen('UU1dD',0,0,0)
        endif

c read UU1dU/D
        read(ifuu)nqbz,nq0i,nwf,nset

        iko_ix = 1
        iko_fx = nwf
        if (nqtt.ne.nset*nqbz) stop
     >      "init_readeigen_nowval1D: nqtt error"

c nqbz mesh-points

        if (is.eq.1) allocate(cbwf(nwf,nwf,nqtt,nsp))

        do iset = 1,nset
          ini = (iset-1)*nqbz
          read(ifmlw)nqbz2,nwf2,iset2,nset2
          if (nqbz.ne.nqbz2) stop
     >         "init_readeigen_nowval1D: nqbz error"
          if (nwf.ne.nwf2) stop
     >         "init_readeigen_nowval1D: nwf error"
          if (iset.ne.iset2) stop
     >         "init_readeigen_nowval1D: iset error"
          if (nset.ne.nset2) stop
     >         "init_readeigen_nowval1D: nset error"
          do iqbz = 1,nqbz
            read(ifmlw)iqbz2,q(1:3)
            if (iqbz2 .ne. iqbz)
     &           stop 'init_readeigen_mlw_noeval1D: iqbz error'
            read(ifmlw)cbwf(1:nwf,1:nwf,ini+iqbz,is)
          enddo
        enddo

c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifuu   = iclose('UU1dU')
        else
          ifmlw  = iclose('MLWD')
          ifuu   = iclose('UU1dD')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif


      deallocate(cbwf)


c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif

ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      deallocate(dnk,uum)


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif


      deallocate(cbwf)


c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_phi_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi and geig
c for hwmat_phi
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbz2,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,nsp2,nlmto2,ngpmx2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,ifdim
      real(8):: q(3),rnorm,cnorm
      complex(8),allocatable :: geig2(:,:),cphi2(:,:)

      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_phi_noeval'

      ifdim = iopen('PHIG.d',1,0,0)
      read(ifdim,*) nsp2,nqbz2,nwf2,nlmto2,ngpmx2
      isx= iclose('PHIG.d')
      if(nsp2.ne.nsp) then
        write(*,*)'nsp,nsp2',nsp,nsp2
        stop 'init_readeigen_phi: ns'
      endif
      if(nlmto2.ne.ldim2) then
        write(*,*)'nlmto,nlmto2',ldim2,nlmto2
        stop 'init_readeigen_phi: nlmto'
      endif
      if(ngpmx2.ne.ngpmx) then
        write(*,*)'ngpmx,ngpmx2',ngpmx,ngpmx2
        stop 'init_readeigen_phi: ngpmx'
      endif
      nwf = nwf2
      nqbz=nqbz2

c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0

c replace geig and cphi

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband

      isx= iclose('CPHI')
      isx= iclose('GEIG')
      ifgeig = iopen('GEIGg',0,0,mrecg_o)
      ifcphi = iopen('CPHIg',0,0,mrecb_o)

      if(keepeig) then
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        do ikp= 1,nqtt
          iqbz=ikp
c         iqbz=mod(ikp,nqbz)
c         if (iqbz.eq.0) iqbz=nqbz
          do is= 1,nsp
            ikpisp = is + nsp*(iqbz-1)
            read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nwf,ikp,is)
            read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nwf,ikp,is)
          enddo
        enddo
      else
c         allocate(geig2(ngpmx,nwf),cphi2(ldim2,nwf))
c         do ikp= 1,nqtt
c           if (ikp.le.nqbz) cycle
c           iqbz=mod(ikp,nqbz)
c           if (iqbz.eq.0) iqbz=nqbz
c         do is= 1,nsp
c           ikpisp = is + nsp*(iqbz-1)
c           read(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           read(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c           ikpisp = is + nsp*(ikp-1)
c           write(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           write(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c         enddo
c         enddo
c         deallocate(geig2,cphi2)
      endif

c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_diag(nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi, geig and evud
c for hx0fp0_mlw, ixc=11 mode
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:),
     &                          dnk2(:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),
     >                        dnk2(iko_ix:iko_fx,nwf))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c nqbze diag-points
        read(ifmlw)nqbze,nwf2
        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nqbze.ne.nqtt .and. nqbze.ne.nqtt-1)
     &      stop 'init_readeigen_mlw: nqbze error'
        read(ifmlwe)nqbze2,nwf2
        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nqbze2.ne.nqbze) stop 'init_readeigen_mlw: nqbze2 error'
        if (is.eq.1)allocate(evec(nwf,nwf,nqtt,nsp),eval(nwf,nqtt,nsp))
        evec = 0d0
        eval = 0d0
        do iqbz = 1,nqbze
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)evec(1:nwf,1:nwf,iqbz,is)

          read(ifmlwe)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen: iqbz error'
          read(ifmlwe)eval(1:nwf,iqbz,is)
        enddo

c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = eval
      deallocate(eval)

check write
c      do iq = 1,nqbze
c      write(6800,*)'***',iq
c      write(6800,*)evec(:,:,iq,1)
c      do iwf = 1,nwf
c         write(6800,*)evud(iwf,iq,1)
c      enddo
c      enddo

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqbze
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk(ib,in,iqbz,is)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
c
          else
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            dnk2 = 0
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                dnk2(ib,iwf) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
c
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk2(ib,in)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
          endif


ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      deallocate(dnk,evec,uum,dnk2)

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif


      deallocate(cbwf)


c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

      end module
