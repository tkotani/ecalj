      module m_rdpp
!> Read PPBRD_V2_*, radial integerals <p|p b> and rotated cg coefficients cgr.
      use NaNum,only: NaN
      use m_genallcf_v3,only: nl,nn,nclass,nspin !,qbas=>pos
      use m_readqg,only: ngcmx
      
      public:: Rdpp
      integer,protected,public:: mdimx=NaN,nbloch=NaN,nxx=NaN,nblochpmx=NaN, nprecx=NaN, mrecl=NaN
      integer,allocatable,protected,public:: nblocha(:) ,lx(:), nx(:,:)
      real(8),allocatable,protected,public:: ppbrd (:,:,:,:,:,:,:), cgr(:,:,:,:)
      logical,protected,public:: done_rdpp=.false.
      private
      integer:: ndble=8
      
      contains
c      subroutine rdpp( nl,ngrp, nn, nclass, nspin,symope,qbas) !nxx,
      subroutine Rdpp( ngrp, symope) !nxx,
!     o      nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr) ---> these are set in m_rdpp.
      implicit none
      integer:: ngrp !,nl,nn,nclass
      integer :: ngpmx,nqbz,nqibz, nband,! nl,ngrp, nxx, nn,
c     i      nblochpmx, !nspin,
     &      n1,n2,n3,iq0,
     &      ifppb(nclass)
      real(8) ::  symope(3,3,ngrp),  pi!,qbas(3,3)
      integer :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1,nwordr
      character*11 :: filename(nclass)
      write(6,*)" rdpp: nclass=",nclass
      if(done_rdpp) call rx('rdpp is already called')
!!  Radial integrals ppbrd
      allocate( nblocha(nclass) ,lx(nclass),
     &     nx(0:2*(nl-1),nclass))
      do ic = 1,nclass
        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        open(newunit=ifppb(ic),file=trim(filename(ic)),form='unformatted')
        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
      enddo
      nxx = maxval( nx )
      allocate( ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrp) )
      write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
      do ic = 1,nclass
        do isp= 1,nspin
          nxic = maxval( nx(0:2*(nl-1),ic) )
          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
        enddo
        close(ifppb(ic))
      enddo
c Belows overide the values given by genallc.
      mdimx  = maxval(nblocha)
      nbloch = sum(nblocha)
      nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB
!! ---------- WV.d
      nprecx = ndble          ! We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()!record size
c --- rotated CG setup
      write(6,*)' rdpp mdimx=',mdimx
      cgr=1d99
      call rotcg(nl-1,symope,ngrp,cgr)
      done_rdpp=.true.
      write(6,*)' rdpp:end '
      end subroutine rdpp
      end module m_rdpp
      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine rdpp_v3(nxx, nl,ngrp, nn, nclass, nspin,symope,!qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr)
      implicit none
      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
     i      nclass,nspin,nn,nblochpmx,nbloch,mdimx,
     &      n1,n2,n3,iq0,
     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
      real(8)    ::  symope(3,3,ngrp) !, pi !qbas(3,3)
      integer(4) :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
      character*11 :: filename(nclass)
      integer:: nx(0:2*(nl-1),nclass)
      real(8):: ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrp) 
      write(6,*)" rdpp_v3: "
!!  Radial integrals ppbrd
      do ic = 1,nclass
        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        open(newunit=ifppb(ic),file=trim(filename(ic)),form='unformatted')
        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
      enddo
      write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
      do ic = 1,nclass
        do isp= 1,nspin
          nxic = maxval( nx(0:2*(nl-1),ic) )
          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
        enddo
        close(ifppb(ic))
      enddo
      mdimx  = maxval(nblocha)
      nbloch = sum(nblocha)
      cgr=1d99
      call rotcg(nl-1,symope,ngrp,cgr)
      write(6,*)' rdpp_v3:end '
      end
!---------------------------------------------------------
      subroutine Getsrdpp2(nclass,nl,nxx) 
      integer(4),intent(in):: nclass,nl
      integer(4) :: nx(0:2*(nl-1),nclass),nxx,ifppb,ic,lxx,nblocha
      character*20 :: filename
      do ic = 1,nclass
        filename = 'PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        open(newunit=ifppb, file=trim(filename), action='read',form='unformatted')
        read(ifppb) nblocha,lxx, nx(0:2*(nl-1),ic)
        close(ifppb)
      enddo
      nxx   = maxval( nx )
      end

c$$$C-----------------------------------------------------------------------------
c$$$      subroutine rdpp_pln2( ngpmx,ngcmx,qibz,nqibz, qbz,nqbz,nband,nspin,
c$$$     o      ngpn,ngvecpB,  ngcni,ngveccB)
c$$$!      nblochpmx = nbloch + ngcmx
c$$$      implicit none
c$$$      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
c$$$     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
c$$$     &      ngpn(nqbz),
c$$$     &      ngvecpB(3,ngpmx,nqbz) ,
c$$$     &      ngcni(nqibz), ! IBZ !
c$$$     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
c$$$     &      iqib(nqbz)
c$$$!     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
c$$$!     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
c$$$c      complex(8) :: geigB(ngpmx,nband,nqbz,nspin) !,img=(0d0,1d0),phase2
c$$$      real(8)    :: qibz(3,nqibz), qbz(3,nqbz) !, symope(3,3,ngrp),
c$$$      complex(8),allocatable:: geig(:,:,:)
c$$$      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
c$$$      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
c$$$     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
c$$$c      logical:: ifgeigb
c$$$      write(6,*)" rdpp_pln2: "
c$$$c --- plane wave contributions 2000 May
c$$$      ifplane = iopen('PLN',0,-1,0)
c$$$      write(6,*)' readin ngp and geig xxxxxxxxxxx'
c$$$      read (ifplane) ngpmx_dum, ngcmx_dum
c$$$      ngveccB =0  !;ngveccBr =0
c$$$      write(6,*)' readin ngp and geig'
c$$$      iqib = 0
c$$$c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
c$$$      do iqi  = 1, nqibz
c$$$        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
c$$$        iqib(iqbzx) = iqi
c$$$cccccccccccccccccccccccc
c$$$c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
c$$$ccccccccccccccccccccccccc
c$$$      enddo
c$$$      do iq = 1,nqbz
c$$$        read(ifplane) ngp, ngc
c$$$        ngpn(iq) = ngp
c$$$        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
c$$$c        write(6,*)' xxx1=',iq
c$$$        read(ifplane) ngvecp, ngvecc, geig
c$$$c        write(6,*)' xxx2=',iq
c$$$c        if(.not.ifgeigb()) then
c$$$c          geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c$$$c        endif
c$$$        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
c$$$c        write(6,*)' xxx3=',iq
c$$$        iqi=iqib(iq)
c$$$        if(iqi/=0 ) then
c$$$ccccccccccccccccccccccccccc
c$$$c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c$$$c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
c$$$cccccccccccccccccccccccccccc
c$$$          ngcni(iqi) = ngc
c$$$          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
c$$$        endif
c$$$        deallocate( geig, ngvecp, ngvecc )
c$$$      enddo
c$$$      is = iclose('PLN')
c$$$      write(6,*)' end of PLN read'
c$$$      return
c$$$      end

c$$$c----------------------------------
c$$$      subroutine rdpp_pln_notused( ngpmx,ngcmx,qibz,nqibz, qbz,nqbz,nband,nspin,
c$$$     o      ngpn,geigB,ngvecpB,  ngcni,ngveccB)
c$$$!      nblochpmx = nbloch + ngcmx
c$$$      implicit none
c$$$      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
c$$$     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
c$$$     &      ngpn(nqbz),
c$$$     &      ngvecpB(3,ngpmx,nqbz) ,
c$$$     &      ngcni(nqibz), ! IBZ !
c$$$     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
c$$$     &      iqib(nqbz)
c$$$!     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
c$$$!     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
c$$$      complex(8) :: geigB(ngpmx,nband,nqbz,nspin) !,img=(0d0,1d0),phase2
c$$$      real(8)    :: qibz(3,nqibz), qbz(3,nqbz) !, symope(3,3,ngrp),
c$$$      complex(8),allocatable:: geig(:,:,:)
c$$$      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
c$$$      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
c$$$     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
c$$$c      logical:: ifgeigb
c$$$      write(6,*)" rdpp_pln: "
c$$$c --- plane wave contributions 2000 May
c$$$      ifplane = iopen('PLN',0,-1,0)
c$$$      read (ifplane) ngpmx_dum, ngcmx_dum
c$$$      ngveccB =0  !;ngveccBr =0
c$$$c      write(6,*)' readin ngp and geig'
c$$$      iqib = 0
c$$$c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
c$$$      do iqi  = 1, nqibz
c$$$        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
c$$$        iqib(iqbzx) = iqi
c$$$cccccccccccccccccccccccc
c$$$c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
c$$$ccccccccccccccccccccccccc
c$$$      enddo
c$$$      do iq = 1,nqbz
c$$$        read(ifplane) ngp, ngc
c$$$        ngpn(iq) = ngp
c$$$        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
c$$$        read(ifplane) ngvecp, ngvecc, geig
c$$$c        if(.not.ifgeigb()) then
c$$$        geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c$$$c        endif
c$$$        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
c$$$        iqi=iqib(iq)
c$$$        if(iqi/=0 ) then
c$$$ccccccccccccccccccccccccccc
c$$$c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c$$$c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
c$$$cccccccccccccccccccccccccccc
c$$$          ngcni(iqi) = ngc
c$$$          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
c$$$        endif
c$$$        deallocate( geig, ngvecp, ngvecc )
c$$$      enddo
c$$$      is = iclose('PLN')
c$$$      write(6,*)' end of PLN read'
c$$$      return
c$$$      end

c$$$c----------------------------------
c$$$      subroutine rdpp_v2( ngpmx,ngcmx,nxx,  qibz,nqibz, qbz,nqbz,
c$$$     i      nband,nl,ngrp, nn,
c$$$     i      nclass, nspin,
c$$$     i      symope,         !qbas,
c$$$     o      nblocha, lx, nx, 
c$$$     o      ppbrd ,
c$$$     o      mdimx,nbloch,  
c$$$     o      cgr,
c$$$     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB)
c$$$c simple ppbrd.
c$$$c-- read radial integerals <p|pb> and so on.
c$$$      implicit none
c$$$      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
c$$$     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
c$$$     &      ngpn(nqbz),
c$$$     &      ngvecpB(3,ngpmx,nqbz) ,
c$$$     &      ngcni(nqibz), ! IBZ !
c$$$     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
c$$$     &      iqib(nqbz),
c$$$     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
c$$$     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
c$$$      complex(8) :: geigB(ngpmx,nband,nqbz,nspin),img=(0d0,1d0),phase2
c$$$      real(8)    :: qibz(3,nqibz), qbz(3,nqbz), symope(3,3,ngrp),
c$$$     & cgr(nl**2,nl**2,(2*nl-1)**2,ngrp),    pi,!qbas(3,3),
c$$$     & ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass)
c$$$c
c$$$      complex(8),allocatable:: geig(:,:,:)
c$$$      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
c$$$      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
c$$$     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
c$$$      character*11 :: filename(nclass)
c$$$c
c$$$c      logical:: ifgeigb
c$$$
c$$$      write(6,*)" rdpp_v2: "
c$$$c --- Radial integrals ppbrd
c$$$      do ic = 1,nclass
c$$$        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
c$$$        ifppb(ic) = iopen(filename(ic),0,-1,0)
c$$$        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
c$$$      enddo
c$$$c      nxx   = maxval( nx )
c$$$c     nspin =ispin
c$$$      write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
c$$$      do ic = 1,nclass
c$$$        do isp= 1,nspin
c$$$c      do ip2= 1,2
c$$$c      do ip1= 1,2
c$$$c        ioff = 1 + (ic-1) + 4*nclass*(isp-1)
c$$$          nxic = maxval( nx(0:2*(nl-1),ic) )
c$$$          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
c$$$c      enddo;  enddo;
c$$$        enddo
c$$$        is= iclose(filename(ic))
c$$$      enddo
c$$$c Belows overide the values given by genallc.
c$$$      mdimx  = maxval(nblocha)
c$$$      nbloch = sum(nblocha)
c$$$c     write(6,*)' imdim',w(imdim),w(imdim+1),w(imdim+2),w(imdim+3)
c$$$c --- rotated CG setup
c$$$      write(6,*)' rdpp mdimx=',mdimx
c$$$c      lmxax = nl-1
c$$$      cgr=1d99
c$$$      call rotcg(nl-1,symope,ngrp,cgr)
c$$$      write(6,*)' end of rotcg'
c$$$c --- plane wave contributions 2000 May
c$$$      ifplane = iopen('PLN',0,-1,0)
c$$$      read (ifplane) ngpmx_dum, ngcmx_dum
c$$$      nblochpmx = nbloch + ngcmx
c$$$      ngveccB =0  !;ngveccBr =0
c$$$c      write(6,*)' readin ngp and geig'
c$$$      iqib = 0
c$$$c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
c$$$      do iqi  = 1, nqibz
c$$$        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
c$$$        iqib(iqbzx) = iqi
c$$$cccccccccccccccccccccccc
c$$$c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
c$$$ccccccccccccccccccccccccc
c$$$      enddo
c$$$      do iq = 1,nqbz
c$$$        read(ifplane) ngp, ngc
c$$$        ngpn(iq) = ngp
c$$$        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
c$$$        read(ifplane) ngvecp, ngvecc, geig
c$$$c        if(.not.ifgeigb()) then
c$$$        geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c$$$c        endif
c$$$        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
c$$$        iqi=iqib(iq)
c$$$        if(iqi/=0 ) then
c$$$ccccccccccccccccccccccccccc
c$$$c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c$$$c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
c$$$cccccccccccccccccccccccccccc
c$$$          ngcni(iqi) = ngc
c$$$          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
c$$$        endif
c$$$        deallocate( geig, ngvecp, ngvecc )
c$$$      enddo
c$$$      is = iclose('PLN')
c$$$      write(6,*)' end of PLN read'
c$$$      return
c$$$
c$$$c-----------------------------------------------------------------------
c$$$c      write(6,*)' goto 2111'
c$$$c      goto 2111
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cROTATION test for planewave part NO.1.
c$$$c
c$$$c        ngvecpB(1:3,1:ngp,  2)  --> ngvecpB(1:3, 1:ngp,  3or4)
c$$$c      write(6,*)' xxxgeig 2 and 6'
c$$$c      geigB(:,1:nband,1:1) =0d0
c$$$c      geigB(:,1:nband,3:7) =0d0
c$$$c      geigB(:,1:nband,5:8) =0d0
c$$$c      return
c$$$
c$$$c      pi = 4*atan(1d0)
c$$$c     iq0 = 4
c$$$c      do iq  = 4,4
c$$$c      do ig1 = 1,ngpn(iq)
c$$$c        phase2 = exp( img*2d0*pi *
c$$$c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c$$$c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c$$$c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c$$$c        n1 = ngvecpB(1,ig1, iq0)
c$$$c        n2 = ngvecpB(2,ig1, iq0)
c$$$c        n3 = ngvecpB(3,ig1, iq0)
c$$$c        ngvecpB(1,ig1, iq) = -n1
c$$$c        ngvecpB(2,ig1, iq) = -n2 -1
c$$$c        ngvecpB(3,ig1, iq) = n3
c$$$c      enddo
c$$$c      enddo
c$$$c
c$$$c      write(6,*)' end of pi/2 test 2'
c$$$c      return
c$$$c
c$$$ccccccccccccccccccccccccccccccccccccc
c$$$cROTATION test for planewave part No3.
c$$$c     iq0 = 2
c$$$c      do iq  = 2,2
c$$$c      do ig1 = 1,ngpn(iq)
c$$$c        phase2 = exp( img*2d0*pi *
c$$$c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c$$$c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c$$$c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c$$$c        n1 = ngvecpB(1,ig1, iq0)
c$$$c        n2 = ngvecpB(2,ig1, iq0)
c$$$c        n3 = ngvecpB(3,ig1, iq0)
c$$$c        ngvecpB(1,ig1, iq) = -n1
c$$$c        ngvecpB(2,ig1, iq) = -n2
c$$$c        ngvecpB(3,ig1, iq) = n3
c$$$c      enddo
c$$$c      enddo
c$$$c
c$$$c      write(6,*)' end of pi/2 test'
c$$$c      return
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c
c$$$c 2111 continue
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cROTATION test for planewave part No2.
c$$$c
c$$$c        ngvecpB(1:3,1:ngp,  2)  --> ngvecpB(1:3, 1:ngp,  3or4)
c$$$c      geigB(:,1:nband,1) =0d0
c$$$c      geigB(:,1:nband,3:5) =0d0
c$$$c      geigB(:,1:nband,7:8) =0d0
c$$$c
c$$$c      pi = 4*atan(1d0)
c$$$c     iq0=6
c$$$c      do iq  = 4,4
c$$$c      do ig1 = 1,ngpn(iq)
c$$$c        phase2 = exp( img*2d0*pi *
c$$$c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c$$$c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c$$$c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c$$$c        geigB(ig1,1:nband,iq) = 0d0
c$$$c        n1 = ngvecpB(1,ig1, iq0)
c$$$c        n2 = ngvecpB(2,ig1, iq0)
c$$$c        n3 = ngvecpB(3,ig1, iq0)
c$$$c        ngvecpB(1,ig1, iq) = -n2
c$$$c        ngvecpB(2,ig1, iq) = n1+n2
c$$$c        ngvecpB(3,ig1, iq) = n3
c$$$c      enddo
c$$$c      enddo
c$$$c      return
c$$$      end
