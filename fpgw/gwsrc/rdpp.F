      module m_rdpp
!> Read PPBRD_V2_*, radial integerals <p|p b> and rotated cg coefficients cgr.
      Use NaNum,only: NaN
c      integer,parameter:: NULL=-9999
      integer,protected:: mdimx=NaN,nbloch=NaN,nxx=NaN
      integer,allocatable,protected:: nblocha(:) ,lx(:), nx(:,:)
      real(8),allocatable,protected:: ppbrd (:,:,:,:,:,:,:), cgr(:,:,:,:)
      logical,protected:: done_rdpp=.false.
      contains
      subroutine rdpp( nl,ngrp, nn, nclass, nspin,symope,qbas) !nxx,
!     o      nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr) ---> these are set in m_rdpp.
      implicit none
      integer:: nl,ngrp,nn,nclass
      integer :: ngpmx,ngcmx, nqbz,nqibz, nband,! nl,ngrp, nxx, nn,
     i      nspin,nblochpmx,
     &      n1,n2,n3,iq0,
     &      ifppb(nclass)
      real(8)    ::  symope(3,3,ngrp),  pi,qbas(3,3)
      integer :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1,ifile_handle
      character*11 :: filename(nclass)
      write(6,*)" rdpp: "
      if(done_rdpp) call rx('rdpp is already called')
!!  Radial integrals ppbrd
      allocate( nblocha(nclass) ,lx(nclass),
     &     nx(0:2*(nl-1),nclass))
      do ic = 1,nclass
        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        open(newunit=ifppb(ic),file=filename(ic),form='unformatted')
        write(6,*)'ifppb=',ifppb(ic),nl
        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
        write(*,*)'srrrrr bbbbb',ic,nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
      enddo
      nxx = maxval( nx )
      allocate(  ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
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
c     call icopy(nclass, nblocha,w(imdim))
c     write(6,*)' imdim',w(imdim),w(imdim+1),w(imdim+2),w(imdim+3)
c --- rotated CG setup
      write(6,*)' rdpp mdimx=',mdimx
c      lmxax = nl-1
      cgr=1d99
      call rotcg(nl-1,symope,ngrp,cgr)
      done_rdpp=.true.
      write(6,*)' rdpp:end '
      end subroutine rdpp
      end module m_rdpp



      subroutine rdpp_v3(nxx, nl,ngrp, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr)
      implicit none
      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
     i      nclass,nspin,nn,nblochpmx,nbloch,mdimx,
     &      n1,n2,n3,iq0,
     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
      real(8)    ::  symope(3,3,ngrp),  pi,qbas(3,3)
      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
      character*11 :: filename(nclass)
      integer:: nx(0:2*(nl-1),nclass)
      real(8):: ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrp) 
      write(6,*)" rdpp_v3: "
!!  Radial integrals ppbrd
      do ic = 1,nclass
        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        ifppb(ic) = iopen(filename(ic),0,-1,0)
        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
      enddo
c      nxx   = maxval( nx )
c     nspin =ispin
      write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
      do ic = 1,nclass
        do isp= 1,nspin
c      do ip2= 1,2
c      do ip1= 1,2
c        ioff = 1 + (ic-1) + 4*nclass*(isp-1)
          nxic = maxval( nx(0:2*(nl-1),ic) )
          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
c      enddo;  enddo;
        enddo
        is= iclose(filename(ic))
      enddo
c Belows overide the values given by genallc.
      mdimx  = maxval(nblocha)
      nbloch = sum(nblocha)
c     call icopy(nclass, nblocha,w(imdim))
c     write(6,*)' imdim',w(imdim),w(imdim+1),w(imdim+2),w(imdim+3)
c --- rotated CG setup
      write(6,*)' rdpp mdimx=',mdimx
c      lmxax = nl-1
      cgr=1d99
      call rotcg(nl-1,symope,ngrp,cgr)
      write(6,*)' rdpp_v3:end '
      end

      subroutine Getsrdpp2(nclass,nl,nxx) 
!     ! Get dimension for PLN and PPBRD_V2_/atom
      integer(4),intent(in):: nclass,nl
      integer(4)  :: nx(0:2*(nl-1),nclass),nxx,ifppb
      character*11 :: filename
      do ic = 1,nclass
        filename = 'PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        open(newunit=ifppb, file=trim(filename), action='read',form='unformatted')
        read(ifppb) nblocha,lxx, nx(0:2*(nl-1),ic)
        close(ifppb)
      enddo
      nxx   = maxval( nx )
      end

      subroutine getsrdpp(nclass,nl, ngpmx,ngcmx,nxx)
c get dimensiton for PLN and PPBRD_V2_/atom
      integer(4)  :: nx(0:2*(nl-1),nclass)
      character*11 :: filename
      do ic = 1,nclass
        filename = 'PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        ifppb = iopen(filename,0,-1,0)
        read(ifppb) nblocha,lx, nx(0:2*(nl-1),ic)
        is = iclose(filename)
      enddo
      nxx   = maxval( nx )
      ifplane = iopen('PLN',0,-1,0)
      read (ifplane) ngpmx, ngcmx
      is = iclose('PLN')
      end


C-----------------------------------------------------------------------------
      subroutine rdpp_pln2( ngpmx,ngcmx,qibz,nqibz, qbz,nqbz,nband,nspin,
     o      ngpn,ngvecpB,  ngcni,ngveccB)
!      nblochpmx = nbloch + ngcmx
      implicit none
      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
     &      ngpn(nqbz),
     &      ngvecpB(3,ngpmx,nqbz) ,
     &      ngcni(nqibz), ! IBZ !
     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
     &      iqib(nqbz)
!     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
!     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
c      complex(8) :: geigB(ngpmx,nband,nqbz,nspin) !,img=(0d0,1d0),phase2
      real(8)    :: qibz(3,nqibz), qbz(3,nqbz) !, symope(3,3,ngrp),
      complex(8),allocatable:: geig(:,:,:)
      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
c      logical:: ifgeigb
      write(6,*)" rdpp_pln2: "
c --- plane wave contributions 2000 May
      ifplane = iopen('PLN',0,-1,0)
      write(6,*)' readin ngp and geig xxxxxxxxxxx'
      read (ifplane) ngpmx_dum, ngcmx_dum
      ngveccB =0  !;ngveccBr =0
      write(6,*)' readin ngp and geig'
      iqib = 0
c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
      do iqi  = 1, nqibz
        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
        iqib(iqbzx) = iqi
cccccccccccccccccccccccc
c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
ccccccccccccccccccccccccc
      enddo
      do iq = 1,nqbz
        read(ifplane) ngp, ngc
        ngpn(iq) = ngp
        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
c        write(6,*)' xxx1=',iq
        read(ifplane) ngvecp, ngvecc, geig
c        write(6,*)' xxx2=',iq
c        if(.not.ifgeigb()) then
c          geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c        endif
        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
c        write(6,*)' xxx3=',iq
        iqi=iqib(iq)
        if(iqi/=0 ) then
ccccccccccccccccccccccccccc
c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
cccccccccccccccccccccccccccc
          ngcni(iqi) = ngc
          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
        endif
        deallocate( geig, ngvecp, ngvecc )
      enddo
      is = iclose('PLN')
      write(6,*)' end of PLN read'
      return
      end

c----------------------------------
      subroutine rdpp_pln_notused( ngpmx,ngcmx,qibz,nqibz, qbz,nqbz,nband,nspin,
     o      ngpn,geigB,ngvecpB,  ngcni,ngveccB)
!      nblochpmx = nbloch + ngcmx
      implicit none
      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
     &      ngpn(nqbz),
     &      ngvecpB(3,ngpmx,nqbz) ,
     &      ngcni(nqibz), ! IBZ !
     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
     &      iqib(nqbz)
!     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
!     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
      complex(8) :: geigB(ngpmx,nband,nqbz,nspin) !,img=(0d0,1d0),phase2
      real(8)    :: qibz(3,nqibz), qbz(3,nqbz) !, symope(3,3,ngrp),
      complex(8),allocatable:: geig(:,:,:)
      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
c      logical:: ifgeigb
      write(6,*)" rdpp_pln: "
c --- plane wave contributions 2000 May
      ifplane = iopen('PLN',0,-1,0)
      read (ifplane) ngpmx_dum, ngcmx_dum
      ngveccB =0  !;ngveccBr =0
c      write(6,*)' readin ngp and geig'
      iqib = 0
c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
      do iqi  = 1, nqibz
        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
        iqib(iqbzx) = iqi
cccccccccccccccccccccccc
c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
ccccccccccccccccccccccccc
      enddo
      do iq = 1,nqbz
        read(ifplane) ngp, ngc
        ngpn(iq) = ngp
        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
        read(ifplane) ngvecp, ngvecc, geig
c        if(.not.ifgeigb()) then
        geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c        endif
        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
        iqi=iqib(iq)
        if(iqi/=0 ) then
ccccccccccccccccccccccccccc
c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
cccccccccccccccccccccccccccc
          ngcni(iqi) = ngc
          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
        endif
        deallocate( geig, ngvecp, ngvecc )
      enddo
      is = iclose('PLN')
      write(6,*)' end of PLN read'
      return
      end

c----------------------------------
      subroutine rdpp_v2( ngpmx,ngcmx,nxx,  qibz,nqibz, qbz,nqbz,
     i      nband,nl,ngrp, nn,
     i      nclass, nspin,
     i      symope,         qbas,
     o      nblocha, lx, nx, 
     o      ppbrd ,
     o      mdimx,nbloch,  
     o      cgr,
     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB)
c simple ppbrd.
c-- read radial integerals <p|pb> and so on.
      implicit none
      integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp,
     i      nclass,nspin,nn,mdimx,nbloch,nblochpmx,
     &      ngpn(nqbz),
     &      ngvecpB(3,ngpmx,nqbz) ,
     &      ngcni(nqibz), ! IBZ !
     &      ngveccB(3,ngcmx,nqibz), !,ngveccBr(3,ngcmx,nqibz),
     &      iqib(nqbz),
     &      nx(0:2*(nl-1),nclass), n1,n2,n3,iq0,
     &      nblocha(nclass) ,lx(nclass),ifppb(nclass)
      complex(8) :: geigB(ngpmx,nband,nqbz,nspin),img=(0d0,1d0),phase2
      real(8)    :: qibz(3,nqibz), qbz(3,nqbz), symope(3,3,ngrp),
     & cgr(nl**2,nl**2,(2*nl-1)**2,ngrp),    pi,qbas(3,3),
     & ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass)
c
      complex(8),allocatable:: geig(:,:,:)
      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:)
      integer(4) :: iclose,is,iopen,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,
     &  ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
      character*11 :: filename(nclass)
c
c      logical:: ifgeigb

      write(6,*)" rdpp_v2: "
c --- Radial integrals ppbrd
      do ic = 1,nclass
        filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
        ifppb(ic) = iopen(filename(ic),0,-1,0)
        read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
      enddo
c      nxx   = maxval( nx )
c     nspin =ispin
      write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
      do ic = 1,nclass
        do isp= 1,nspin
c      do ip2= 1,2
c      do ip1= 1,2
c        ioff = 1 + (ic-1) + 4*nclass*(isp-1)
          nxic = maxval( nx(0:2*(nl-1),ic) )
          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
c      enddo;  enddo;
        enddo
        is= iclose(filename(ic))
      enddo
c Belows overide the values given by genallc.
      mdimx  = maxval(nblocha)
      nbloch = sum(nblocha)
c     call icopy(nclass, nblocha,w(imdim))
c     write(6,*)' imdim',w(imdim),w(imdim+1),w(imdim+2),w(imdim+3)
c --- rotated CG setup
      write(6,*)' rdpp mdimx=',mdimx
c      lmxax = nl-1
      cgr=1d99
      call rotcg(nl-1,symope,ngrp,cgr)
      write(6,*)' end of rotcg'
c --- plane wave contributions 2000 May
      ifplane = iopen('PLN',0,-1,0)
      read (ifplane) ngpmx_dum, ngcmx_dum
      nblochpmx = nbloch + ngcmx
      ngveccB =0  !;ngveccBr =0
c      write(6,*)' readin ngp and geig'
      iqib = 0
c      call dcopy(3*nqibz, w(iqibz),1, qibz(:,:),1)
      do iqi  = 1, nqibz
        iqbzx =  idxk (qibz(1:3,iqi),qbz,nqbz)
        iqib(iqbzx) = iqi
cccccccccccccccccccccccc
c      write(6,"(' qibz=', i4,3f12.5)") iqi, qibz(1:3,iqi)
ccccccccccccccccccccccccc
      enddo
      do iq = 1,nqbz
        read(ifplane) ngp, ngc
        ngpn(iq) = ngp
        allocate( geig(ngp,nband,nspin), ngvecp(3,ngp), ngvecc(3,ngc) )
        read(ifplane) ngvecp, ngvecc, geig
c        if(.not.ifgeigb()) then
        geigB(1:ngp,1:nband,iq,1:nspin)= geig(1:ngp,1:nband,1:nspin)
c        endif
        ngvecpB(1:3,1:ngp,iq)  = ngvecp(1:3,1:ngp)
        iqi=iqib(iq)
        if(iqi/=0 ) then
ccccccccccccccccccccccccccc
c      write(6,"(' ## qibz=', i4,3f12.5)")iqi, qibz(1:3,iqi)
c      write(6,"(' ## qbz =', i4,3f12.5)")iq ,  qbz(1:3,iq)
cccccccccccccccccccccccccccc
          ngcni(iqi) = ngc
          ngveccB(1:3,1:ngc,iqi) = ngvecc(1:3,1:ngc)
        endif
        deallocate( geig, ngvecp, ngvecc )
      enddo
      is = iclose('PLN')
      write(6,*)' end of PLN read'
      return

c-----------------------------------------------------------------------
c      write(6,*)' goto 2111'
c      goto 2111
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cROTATION test for planewave part NO.1.
c
c        ngvecpB(1:3,1:ngp,  2)  --> ngvecpB(1:3, 1:ngp,  3or4)
c      write(6,*)' xxxgeig 2 and 6'
c      geigB(:,1:nband,1:1) =0d0
c      geigB(:,1:nband,3:7) =0d0
c      geigB(:,1:nband,5:8) =0d0
c      return

c      pi = 4*atan(1d0)
c     iq0 = 4
c      do iq  = 4,4
c      do ig1 = 1,ngpn(iq)
c        phase2 = exp( img*2d0*pi *
c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c        n1 = ngvecpB(1,ig1, iq0)
c        n2 = ngvecpB(2,ig1, iq0)
c        n3 = ngvecpB(3,ig1, iq0)
c        ngvecpB(1,ig1, iq) = -n1
c        ngvecpB(2,ig1, iq) = -n2 -1
c        ngvecpB(3,ig1, iq) = n3
c      enddo
c      enddo
c
c      write(6,*)' end of pi/2 test 2'
c      return
c
ccccccccccccccccccccccccccccccccccccc
cROTATION test for planewave part No3.
c     iq0 = 2
c      do iq  = 2,2
c      do ig1 = 1,ngpn(iq)
c        phase2 = exp( img*2d0*pi *
c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c        n1 = ngvecpB(1,ig1, iq0)
c        n2 = ngvecpB(2,ig1, iq0)
c        n3 = ngvecpB(3,ig1, iq0)
c        ngvecpB(1,ig1, iq) = -n1
c        ngvecpB(2,ig1, iq) = -n2
c        ngvecpB(3,ig1, iq) = n3
c      enddo
c      enddo
c
c      write(6,*)' end of pi/2 test'
c      return
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 2111 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cROTATION test for planewave part No2.
c
c        ngvecpB(1:3,1:ngp,  2)  --> ngvecpB(1:3, 1:ngp,  3or4)
c      geigB(:,1:nband,1) =0d0
c      geigB(:,1:nband,3:5) =0d0
c      geigB(:,1:nband,7:8) =0d0
c
c      pi = 4*atan(1d0)
c     iq0=6
c      do iq  = 4,4
c      do ig1 = 1,ngpn(iq)
c        phase2 = exp( img*2d0*pi *
c     &      sum((qbz(1:3,iq0)+ matmul(qbas, ngvecpB(1:3,ig1,iq0)) )
c     &          *(/.5773505d0,0d0,0.8135d0/) ) )
c        geigB(ig1,1:nband,iq) = geigB(ig1,1:nband,iq0)*phase2
c        geigB(ig1,1:nband,iq) = 0d0
c        n1 = ngvecpB(1,ig1, iq0)
c        n2 = ngvecpB(2,ig1, iq0)
c        n3 = ngvecpB(3,ig1, iq0)
c        ngvecpB(1,ig1, iq) = -n2
c        ngvecpB(2,ig1, iq) = n1+n2
c        ngvecpB(3,ig1, iq) = n3
c      enddo
c      enddo
c      return
      end
