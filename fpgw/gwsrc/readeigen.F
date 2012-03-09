      module m_readeigen
!!-- Return eigenvalus and eigenfunctions for given q and isp.--
!! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
!! qtti(1:3,nqtti) :eivenvalues, eigenvectors are calculated only for them.
!!                  See lmfgw (q-vector with irr flag in QGpsi).
!!
!!* we have to call init_foobar to call readeval, readcphi, readgeig.
!!----------------
      implicit none
c      integer(4),private :: nqtt
c      real(8),allocatable,private:: qtt(:,:)
      real(8),allocatable,private:: evud(:,:,:), ginv(:,:)
      logical,private:: init=.true.,init2=.true.,keepeig
      integer(4),allocatable,private:: ngp(:)
      integer(4),private:: nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   nnnn, nsp,ldim2,mrecg,ifcphi,ifgeig
      complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
      real(8),private:: leval
ctm
      integer(4),private :: nwf,ifcphi_mlw,ifgeig_mlw
      complex(8),allocatable,private:: geig_mlw(:,:,:,:),
     &           cphi_mlw(:,:,:,:)
      logical,private:: debug=.true.
      integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
      integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:)
     & ,offset_tbl(:),offsetrev_tbl(:,:,:)
!!--------------------------------------------------------------
      contains

      subroutine readeval(q,isp, ev)
!!-- return ev(1:nband) for given q(1:3) and isp --
      use m_hamindex, only: qtt,nqtt,iqimap
      implicit none
      integer(4):: isp,iq,iqindx,i
      real(8):: q(3)
      real(8):: ev(nband)
      if(init) then
        stop 'readeigen: modele is not initialized yet'
      endif
      iq = iqindx(q, ginv,qtt,nqtt)
      print *,'iq iqimap(iq)=',iq,iqimap(iq)
      ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
      if(debug) then
        write(6,"('iq iqimap(iq)=',2i4,3f13.5)")iq,iqimap(iq),q
        write(6,"(9f9.4)")ev(1:9)
      endif
      end subroutine

!!-----------------------------------------------
      subroutine readgeig(q,ngp_in,isp, qu,geigen)
!!-- return ev(1:nband) for given q(1:3) and isp --
      use m_hamindex,only: iqmap,igmap,qtt,ngpmx,nqtt,shtvg,qlat,symops,plat,invgx,iqimap
      implicit none
      integer:: isp,iq,iqindx,ngp_in,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi
      real(8)   :: q(3),qu(3),ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3)
      complex(8):: geigen(ngp_in,nband),geigenr(ngp_in,nband),img=(0d0,1d0),img2pi
      logical ::debug=.true.
      img2pi=2d0*4d0*datan(1d0)*img
      platt=transpose(plat) !this is inverse of qlat
      if(init2) then
        stop 'readgeig: modele is not initialized yet'
      endif
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q. q-qu is a G vector.
      print *,' readgeig:xxx iq=',iq
c      do iq=1,nqtt
c         write(*,*)iq,qtt(:,iq)
c      enddo
      iqq=iqmap(iq)  
      iqi=iqimap(iq)  
      igg=igmap(iq)
      qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
      if(debug) then
        print *,' ddddd ginv=',ginv
        print *,' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        print *,' ddddd q      =',q
        print *,' ddddd qu      =',qu
        print *,' ddddd qtarget=',qtarget
        print *,' ddddd q (iqq)=',qtt(:,iqq)
      endif   
!!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
      if(ngp(iq)==0) return
      if(ngp(iq)/=ngp(iqq)) stop 'readgeig: ngp(iq)/=ngp(iqq)'
      if(debug) print *,' ngp=',ngp(iq)
      if(keepeig) then
        geigenr(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifgeig, rec=ikpisp) geigenr(1:ngpmx,1:nband)
      endif
ccccccccccccccccccccccccccccccccccccccccccc
c      geigen=0d0
c      return
ccccccccccccccccccccccccccccccccccccccccccc
!! rotation wave
      ddd = shtvg(:,invgx(igg))
      if(debug) then
        print *,' ddddd=',igg,ddd
        write(*,"('ddddd iq =',2i5,3f13.5)") iq,  ngp(iq), qtt(:,iqq)
        write(*,"('ddddd iqq=',2i5,3f13.5,2x,3f13.5)") iqq, ngp(iq), q,qu
        print *,' ddddd '
      endif
   
      do ig = 1,ngp(iqq)
        if(debug) print *,' ig=',ig
        qpg = qtt(:,iqq) + matmul( qlat(:,:),ngvecp(:,ig,iqq)) !q+G
        if(debug) print *,' qpg=',qpg
        qpgr = matmul(symops(:,:,igg),qpg)             !rotated q+G
        if(debug) print *,' qpgr=',qpgr
        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3),iq)
        if(debug) print *,' nnn=',nnn
        geigen(ig2,:)= geigenr(ig,:) * exp( img2pi*sum(qpg*ddd) )
      enddo
      if(ngp_in < ngp(iq)) then
        print *,'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        stop 'readgeig: ngpmx<ngp(iq)'
      endif
      end subroutine readgeig
c------------------------------
      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
      use m_hamindex,only: iqmap,igmap,qtt,ngpmx,nqtt,qlat,plat,invgx,iqimap,
     &  miat,tiat,dlmm,norbt,shtvg
!!-- return mto part of eigenfunction for given q(1:3) and isp
      implicit none
      integer:: isp,iq,iqindx,ldim2_in,ikpisp,iqq,iorb,ibaso,ibas,k,l,
     & ini1,ini2,iend1,iend2,igg,ig,iqi
      real(8)   :: q(3),qu(3),qtarget(3)
      complex(8):: img2pi,img=(0d0,1d0),phase,cphifr(ldim2_in,nband),cphif(ldim2_in,nband),phaseq !takao 1->*->nband
      img2pi=2d0*4d0*datan(1d0)*img
      if(init2) then
        stop 'readcphi: modele is not initialized yet'
      endif
      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
!! iqmap(iq): iq is rotated to igmap(iq) (qtt(1:nqtt)) by sympos(  ,igmap(iq)).
      iqq=iqmap(iq)
      iqi=iqimap(iq)
      qtarget=qtt(:,iq)
      igg=igmap(iq)
      if(keepeig) then
        cphifr(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifcphi, rec=ikpisp) cphifr(1:ldim2,1:nband)
      endif
cccccccccccccccccccccccccccccccccc
       write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ldim2_in,ldim2,norbt,
     &                sum(cphifr(1:ldim2,1:nband)),sum(abs(cphifr(1:ldim2,1:nband)))
cccccccccccccccccccccccccccccccccc

!! rotation  See mto part of rotwv xxxxxxxxxx
      phaseq= exp( img2pi*sum( qtt(:,iqq)*shtvg(:,invgx(igg)) ) )
      ibaso=-999
      do iorb=1,norbt !orbital-blocks are specified by ibas, l, and k.
        !Index of Hamiltonian is devided into these blocks.
        ibas = ibas_tbl(iorb)
        if(ibas/=ibaso) phase = exp( -img2pi*sum(qtarget*(tiat(:,ibas,igg)-shtvg(:,igg))) )
        ibaso=ibas
        l   = l_tbl(iorb)
        k   = k_tbl(iorb)
        ini1 = offset_tbl(iorb)+1
        iend1 = ini1+2*l
        ini2 = offsetrev_tbl(miat(ibas,igg),l,k)+1
        iend2 = ini2+2*l
        cphif(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),cphifr(ini1:iend1,:))*phase*phaseq
        write(6,*)ibaso,ibas,phase
        write(6,*)ini1,iend1,ini2,iend2
        write(6,*)sum(abs(dlmm(-l:l,-l:l,l,igg)))
        write(6,*)l,k,igg,sum(abs(cphifr(ini1:iend1,:)))
        write(6,*)l,k,igg,sum(abs(cphif(ini2:iend2,:)))
      enddo
ccccccccccccccccccccccccccccccccc
c       write(6,"('readcphi:: xxxx sum of cphifr=',2i4,4d23.16)")ldim2_in,ldim2,
c     &                sum(cphif(1:ldim2,1:nband)),sum(abs(cphif(1:ldim2,1:nband)))
c       stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxx takao test xxxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccccc
      end subroutine readcphi


c$$$c------------------------------
c$$$      subroutine readgeig(q,ngp_in,isp, qu,geigen)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: geigen(ngp_in,nband)
c$$$      if(init2) then
c$$$        stop 'readgeig: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(ngp_in < ngp(iq)) then
c$$$        print *,'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
c$$$        stop 'readgeig: ngpmx<ngp(iq)'
c$$$      endif
c$$$! for ngp=0 Aug2005
c$$$      if(ngp(iq)==0) return
c$$$      if(keepeig) then
c$$$        geigen(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifgeig, rec=ikpisp) geigen(1:ngpmx,1:nband)
c$$$      endif
c$$$      end subroutine readeig
c$$$c------------------------------
c$$$      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: cphif(ldim2_in,nband) !takao 1->*->nband
c$$$      if(init2) then
c$$$        stop 'readcphi: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(keepeig) then
c$$$        cphif(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifcphi, rec=ikpisp) cphif(1:ldim2,1:nband)
c$$$      endif
c$$$      end subroutine readcphi

c------------------------------
      subroutine readgeig_mlw(q,ngp_in,isp, qu,geigen)
      use m_hamindex
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
      use m_hamindex
!!-- return ev(1:nband) for given q(1:3) and isp
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
      use m_hamindex
!!-- initialization. Save QpGpsi EVU EVD to arrays.--
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,nsp_in, 
     & nband_in,mrece_in,iopen,iclose,isx,mrecb_in,ik,ib,verbose
      real(8):: ginv_in(3,3)
c      integer(4),allocatable:: ngvecp(:,:)
      real(8):: QpGcut_psi
      integer:: ifev,nband_ev, nqtti_, nsp_ev !,nqnum
      real(8),allocatable:: qtt_(:,:),qtti_(:,:)
      write(6,*) 'init_readeigen:'
      allocate(ginv(3,3))
      ginv  = ginv_in
      nsp   = nsp_in
      if(nsp<0.or.nsp>2) stop 'init_reaeigen:nsp wrong'
      nband = nband_in
      mrece = mrece_in
c      write(6,*) 'init_readeigen:1'

!!--- store requied index
      call readhamindex() 
c      write(*,*)ngrp,nbas,kxx,lxx
c      write(*,*)'nqtt=',nqtt,nqtti,nqnum
c      write(*,*)'imx=',imx,ngpmx
c      do iq=1,nqtt
c        write(*,*) qtt(:,iq)
c      enddo
c      stop 'xxxxxxxxxxxxxxxxxx test xxxxx'
!!--- Readin QGpsi ---
      ifiqg  = 101
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum , ngpmx, QpGcut_psi, nnnn,nqtti,imx
      write(6,*)'read(ifiqg)', nqnum , ngpmx, QpGcut_psi, nnnn
      allocate( qtt_(3,nqnum),ngp(nqnum) )
      allocate( ngvecp(3,ngpmx,nqnum))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
      do ikp = 1,nqnum
        read (ifiqg) qtt_(:,ikp), ngp(ikp)
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp) 
      enddo
      close(ifiqg)
      deallocate(qtt_)

!!--- Read eigenvalue ---
      ifev= iopen('EValue',0,-1,0)    !feb2012takao
      read(ifev) nband_ev, nqtti_, nsp_ev
      write(*,*)'read EValue: nband_ev,nqtti,nsp_ev', nband_ev, nqtti_, nsp_ev
      if(nband_ev/=nband) stop 'init_readeigen:nband_ev/=nband'
      if(nsp_ev  /=  nsp) stop 'init_readeigen:nsp_ev/=nsp'
      if(nqtti  /=  nqtti_) stop 'init_readeigen:nqtti/=nqtti_'
      allocate(evud(nband,nqtti_,nsp),qtti_(3,nqtti_))
      read(ifev) qtti_(1:3,1:nqtti_) !qtti
      read(ifev) evud(1:nband, 1:nqtti, 1:nsp )
      close(ifev)
      if(debug) then
        do is= 1,nsp
        do ik= 1,nqtt
        do ib= 1,nband
          write(6,"('ib ik e=',2i5,f13.5)") ib,ik,evud(ib,ik,is)
        enddo
        enddo
        enddo
      endif

c comment out mar2012takao
c      ifev(1) = iopen('EVU', 0,0,mrece)
c      if(nsp==2) ifev(2) = iopen('EVD', 0,0,mrece)
C$$$      do is= 1,nsp
C$$$        call rwdd (ifev(is),   !read eigenvalues
C$$$     d              nband,nqtti, 
C$$$     o              evud(1,1,is))
C$$$c        print *, 'evud is=',is,minval(evud(:,:,is))
C$$$        if(debug) then
C$$$          do ik=1,nqtt
C$$$            do ib=1,nband
C$$$              write(6,"('ib ik e=',2i5,f13.5)") ib,ik,evud(ib,ik,is)
C$$$            enddo
C$$$          enddo
C$$$        endif
C$$$      enddo
C$$$c      stop
C$$$      isx = iclose('EVU')
C$$$      if(nsp==2) isx = iclose('EVD')
      leval= minval(evud)
      init=.false.
      end subroutine

      real(8) function lowesteval()
      lowesteval=leval
      end function
c-------------------------------------------------------------
      subroutine init_readeigen2(mrecb_in,ldim2_in,mrecg_in)
      use m_hamindex
c this should be called after init_readgeigen
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,ldim2_in, 
     & iopen,iclose,isx,mrecb_in,ikpisp,mrecg_in,verbose,ifoc,
     & i1,i2,i3,i4,i5,iorb,iorbold
      logical :: keepeigen
      mrecb = mrecb_in
      ldim2 = ldim2_in
      mrecg = mrecg_in
      keepeig = keepeigen()

      init2=.false.

!! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      read(ifoc,*)
      norbt=0
      write(*,*) ' goto readin @MNLA_CPHI: norbt=',norbt
      do 
        read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
        write(6,*) i1,i2,i3,i4,i5,iorb
        if(iorb>norbt) norbt=iorb
      enddo
 106  continue
      write(*,*) ' end of readin @MNLA_CPHI: norbt=',norbt
      rewind ifoc
      read(ifoc,*)
      allocate(l_tbl(norbt),k_tbl(norbt),ibas_tbl(norbt),offset_tbl(norbt))
      iorbold=0
      do 
        read(ifoc,*,end=107)i1,i2,i3,i4,i5,iorb
        if(iorbold/=iorb) then
          k_tbl(iorb)=i2
          l_tbl(iorb)=i3
          ibas_tbl(iorb)=i4
          offset_tbl(iorb)=i5-1
          iorbold=iorb
        endif
      enddo
 107  continue
      allocate(offsetrev_tbl(maxval(ibas_tbl),0:maxval(l_tbl),maxval(k_tbl)))
      offsetrev_tbl=-99999
      do iorb=1,norbt
        offsetrev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))= offset_tbl(iorb)
      enddo
      write(*,*)'init_readeigen2: table for rotation of cphi' 
      do iorb=1,norbt
        write(*,*)ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb),
     &   offset_tbl(iorb),
     &  offsetrev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))
      enddo
      ifoc=iclose('@MNLA_CPHI')
c      stop 'init2 xxxxxxxxxxxxx'

c --- Readin CPHI and GEIG
      if(Keepeig     ) print *,' KeepEigen=T; readin geig and cphi into m_readeigen'
      if(.not.Keepeig) print *,' KeepEigen=F; not keep geig and cphi in m_readeigen'

      ifgeig = iopen('GEIG',0,0,mrecg)
      ifcphi = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction
      print *,' ifgeig ifcphi=',ifgeig,ifcphi
      if(.not.keepeig) return
c      ifgeig = iopen('GEIG',0,-1,0)
c      ifcphi = iopen('CPHI',0,-1,0) ! Coeeficients of eigenfunction
      allocate(geig(ngpmx,nband,nqtti,nsp))
      allocate(cphi(ldim2,nband,nqtti,nsp))
      do ikp= 1,nqtti
        do is= 1,nsp
c        read(ifgeig) geig(1:ngpmx,1:nband,ikp,is) !geig(1:ngp(ikp),1:nband,ikp,1:is)
c        read(ifcphi) cphi(1:ldim2,1:nband,ikp,is)
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
ccccccccccccccccccccccc
          if(verbose()>40) 
     &       write(6,'("Readin geig sumcheck ",i5,3f9.4,"  ",i5,2d13.6)')
     &      ikp,qtti(:,ikp),ngp(ikp),sum(geig(1:ngp(ikp),1:nband,ikp,is))
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
      use m_hamindex
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
      use m_hamindex

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
      use m_hamindex

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
      use m_hamindex

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
      use m_hamindex

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

c taken from lm7k/subs/rotwv.F
      subroutine conv2int(a,ii,n)
      implicit none
      integer::ii(n),i,n
      real(8):: a(n),eps=1d-8,aa
      do i=1,n
        if(a(i)>0) then
          aa = a(i)
          ii(i)=  int(aa+eps)
        else
          aa=-a(i)
          ii(i)= -int(aa+eps)
        endif
        if(abs(ii(i)-a(i))>eps)stop'conv2int: a is not close to integer'
      enddo
      end subroutine conv2int
c
      end module
