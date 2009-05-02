      module m_genallcf_v3
c----------------------------------------------------------------
c> structure
c     o                   plat,alat,natom,nclass,pos,
c     o                   ngrp, symgg,
c     o                   invg, ef,
c>> l,n and dimensions
c     o                   clabl, nspin,nl,nn,nnv,nnc,
c     o                   nindx, nindxv, nindxc, iclass,
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c
c     o                   z,
c>> l,n,m indices for Phi
c     o                   il, in, im,  ilnm,  nlnm,
c     o                   ilv,inv,imv, ilnmv, nlnmv,
c     o                   ilc,inc,imc, ilnmc, nlnmc,
c>> core
c     o                   ncwf, ecore, konf, icore, ncore,nctot,
c> frequency
c     o                   niw,diw,nw,dw,delta,deltaw,esmr, freq)
c             symgrp
c             ,nocc, nunocc, occv, unoccv, occc, unoccc
c----------------------------------------------------------------
cccc      use m_rgwinf_v3
      implicit none
      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot, niw,nw
      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
      character(120):: symgrp
      character(6),allocatable :: clabl(:)
      integer(4),allocatable:: iclass(:)
     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
     &    occv(:,:,:),unoccv(:,:,:)   
     &   ,occc(:,:,:),unoccc(:,:,:),
     o    nocc(:,:,:),nunocc(:,:,:)
      real(8), allocatable:: 
     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:)   !w(igrp)
c-----------------------------------------------------------------------
      contains
      subroutine genallcf_v3(nwin,efin,incwfx)
C- Readin GWIN_V2 and LMTO(crystal) data and allocate all required.
C
Cr Return iclass=ibas.
Cr nwin,efin,incwfx, are used as switches.
Ci input: nwin,efin,incwfx,
Ci        GWIN_V2, LMTO
Co output: All the output are given in the declear section above.
C--------------------------------------------------------
      implicit none
      integer(4)::iflmto,ifinin,nwin,incwfx,ifec,i,j,
     & lmx, lmx2,nlmto2,nprodxc,nlnaxc,nlnaxv,nprodx,ifi,ig,is
     & ,iopen,iclose,nprodxv,nlnax
     & ,noflmto,maxnn
      integer(4):: infwfx
      integer(4):: n1,n2,n3,imagw,lcutmx,n,ic
      logical :: nocore
      real(8)::efin
      real(8),allocatable::tolbas(:)
      character(120):: symgrpt
      real(8),   allocatable:: ecoret(:)
      integer(4),allocatable::ncwf2(:,:,:)

c      logical ::readgwinput
c----------------------------------------------
      iflmto     = iopen('LMTO',1,0,0)
      if (iflmto < 0) stop 'unit file for GWIN_V2 < 0'
c      if(.not.readgwinput()) then
c        ifinin     = iopen('GWIN_V2',1,0,0)
c      endif

c--- readin these by rgwinf_v3
c      character*120 symgrp
c      integer(4)::nclass,natom,nspin,nl,nnv,nnc
c      real(8)::alat
c      integer(4),allocatable::
c     &   iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:)
c     &  ,occv(:,:,:),unoccv(:,:,:)
c     &  ,occc(:,:,:),unoccc(:,:,:)
c     &  ,ncwf(:,:,:)
c      real(8),allocatable:: plat(:,:),pos(:,:),z(:)
c      character*6,allocatable:: clabl(:)
c      print *,' goto rgwin'
c      call rgwinf_v3 (iflmto,ifinin,nwin,efin,incwfx) !these are inputs
c      print *,' end of rgwinf_v3'
c-------------------------------------------------------------------------
c--- rgwinf ---
      ifi = iflmto
      nw  = nwin
      ef  = efin
      read(ifi,*);  read(ifi,*)
      read(ifi,*)symgrpt        !SYMMETRY
      j = 0
      call rmvbl (symgrpt,120,j)
      symgrp = '  '//symgrpt(j+1:120)
      print *,' symgrp=', symgrp

      read(ifi,*)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)alat              !lattice constant
      allocate(plat(3,3))          !primitive lattice vectors
      read(ifi,*)
      read(ifi,*)plat(1:3,1)
      read(ifi,*)plat(1:3,2)
      read(ifi,*)plat(1:3,3)
      read(ifi,*)
      read(ifi,*) natom            !Number of atoms
c-------------------------------------
      nclass = natom  !We set nclass = natom through the GW calculations
c-------------------------------------
      print *,' alat natom=',alat,natom
      allocate(pos(3,natom))       !positions of atoms
      read(ifi,*)
      do n = 1,natom
        read(ifi,*) pos(1,n),pos(2,n),pos(3,n)
      end do
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)nspin           !spin (1=paramagnetic  2=ferromagnetic)
      read(ifi,*)
      read(ifi,*)nl              !max. no. valence and core l
      read(ifi,*)
      read(ifi,*)nnv,nnc   !max. no. valence and core n
      print *,' nspin nl nnv nnc =',nspin,nl,nnv,nnc
c--------------------------
      if(nnv==1) nnv=2 ! for backword compatibility!takao apr 2002
      ! nnv=2 corresponds to phi and phidot
      ! nnv=3 corresponds to
c--------------------------
      read(ifi,*)
      read(ifi,*)  !nrx is not readin
      read(ifi,*)
      allocate(clabl(nclass),z(nclass)) !class-label, z
      do ic = 1,nclass
        read(ifi,*) clabl(ic),z(ic) !,nrofi is not readin
      end do

      allocate(iclass(natom)) !atom and its class.
      do n = 1,natom          !!We set nclass = natom through the GW calculations
        iclass(n)=n
      end do

      allocate(nindxv(nl,nclass), nindxc(nl,nclass),
     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
     &        occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))

      allocate(tolbas(0:2*(nl-1)))
ctakao
c      incwfx = incwf
c      call defi    (incwf2,nl*nnc*nclass)
c      call defi    (incwf,nl*nnc*nclass)
c      print *,' goto rgwinaf use incfw2=',useincfw2
      call rgwinaf    (ifi,ifinin,nl,nnv,nnc,nclass,
c> BZ
     o                  n1,n2,n3,ef,
c> frequencies
     o                  niw,diw,nw,dw,delta,deltaw,esmr,imagw,
c> coulomb
c     o                  tolvc,alp,alptx,h,ng,
c> product basis
     o                  tolbas,lcutmx,nindxv,nindxc,
     o                  occv,unoccv, occc,unoccc,
c> core
     o                  ncwf,ncwf2 )
c----
      inquire(file='NoCore',exist=nocore)
      if(nocore) then
        occc=0    ! call iclear(nl*nnc*nclass, w(ioccc))
        unoccc=0    ! call iclear(nl*nnc*nclass, w(iunoccc))
        ncwf  =0    ! call iclear(nl*nnc*nclass, w(incwf))
      elseif( incwfx==-1 ) then
        print *,' ### incwf=-1 Use ForSxc for core'
        ncwf = ncwf2  !call icopy(nl*nnc*nclass,w(incwf2),w(incwf))
      elseif( incwfx==-2 ) then
        print *,' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
        call notbit(nl*nnc*nclass, ncwf2)
        ncwf  = ncwf2 ! call icopy (nl*nnc*nclass, w(incwf2),w(incwf))
        occc= ncwf  ! call icopy (nl*nnc*nclass, w(incwf),w(ioccc))
        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
      elseif( incwfx==-3 ) then
        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
        print *,' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
      elseif( incwfx==-4 ) then
        print *,' ### incwf=-4  occ=0 and unocc=0 for all core '
        occc=0  !call iclear(nl*nnc*nclass, w(ioccc))
        unoccc=0 !call iclear(nl*nnc*nclass, w(iunoccc))
        ncwf=0  !call iclear(nl*nnc*nclass, w(incwf))
      elseif(incwfx==0) then
        print *,' ### Use unocc occ ForX0 for core'
      else
        stop ' ### proper incwf is not given for genallcf2:rgwinf '
      endif
      deallocate(ncwf2)
C... End of rgwinf section ----------------------------


c> dimensions and constants
      lmx        = 2*(nl-1)
      lmx2       = (lmx+1)**2
      nlmto      = noflmto(nindxv,iclass,nl,nclass,natom)
      nlmto2     = nlmto*nlmto
      nn         = maxnn (nindxv,nindxc,nl,nclass)

c>> combine nocc,nunocc,nindx
      allocate(nindx(nl,nclass))
      allocate(nocc(nl,nn,nclass),nunocc(nl,nn,nclass))
      call reindx  (occv,unoccv,nindxv,
     i              occc,unoccc,nindxc,
     d              nl,nn,nnv,nnc,nclass,
     o              nocc,nunocc,nindx)
      call maxdim  (occc,unoccc,nindxc,nl,nnc,nclass,
     o              nprodxc,nlnxc,nlnmxc,nlnaxc)
      call maxdim  (occv,unoccv,nindxv,nl,nnv,nclass,
     o              nprodxv,nlnxv,nlnmxv,nlnaxv)
      call maxdim  (nocc,nunocc,nindx,nl,nn,nclass,
     o              nprodx,nlnx,nlnmx,nlnax)

c      nlnx4      = nlnx**4
c      nphi       = nrx*nl*nn*nclass
c      pi         = 4d0*datan(1d0)
c      tpia       = 2d0*pi/alat

c> frequency mesh
c      call defdr   (ifreq,nw)
      print *,' nw from rgwinaf=',nw
      if(nw>0) then
        allocate(freq(nw))
        call genfreq (nw,dw,0.d0,
     o              freq )
      endif

c> index for allowed core states
c      call defi    (iicore,nl*nl*nnc*nclass)
c      call defi    (incore,nclass)
      allocate(icore(nl**2*nnc,nclass),ncore(nclass))
      icore=9999999
      ncore=9999999
      call incor   (ncwf,nindxc,iclass,
     d              nl,nnc,nclass,natom,
     o              icore,ncore,nctot )
ccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' nnc=',nnc,nl,nclass,natom
c      print *,' ncwf ',ncwf
c      print *,' nindxc ',nindxc
c      print *,' iclass ',iclass
c      print *,' --- icore=',icore
c      print *,' --- ncore nctot=',ncore,nctot
ccccccccccccccccccccccccccccccccccccccccccccccccc

c> core energies
c      call defi    (ikonf,-nl*nclass)
c      call defdr   (iecore,2*nctot)
c      call defdr   (iecoret,-nl*nnc*2*nclass)
      ifec       = iopen ('ECORE',1,0,0)
      allocate(konf(nl,nclass),ecore(nctot,2))
      konf=0
      allocate(ecoret(nl*nnc*2*nclass))
      ecoret=0d0
      call rcore2   (ifec, ! minor changed version to rcore in rw.f
     d              nl,nnc,nclass,
     o              konf,ecoret)
      call gencor  (ecoret,ncwf,iclass,
     d              nl,nnc,natom,nclass,nctot,
     o              ecore)
c      call rlse    (iecoret)
      deallocate(ecoret)

c> index for core and LMTO basis
c      call defi    (iil,nlnmx*nclass)
c      call defi    (iin,nlnmx*nclass)
c      call defi    (iim,nlnmx*nclass)
c      call defi    (iilnm,nn*nl*nl*nclass)
c      call defi    (iilv,nlnmxv*nclass)
c      call defi    (iinv,nlnmxv*nclass)
c      call defi    (iimv,nlnmxv*nclass)
c      call defi    (iilnmv,nnv*nl*nl*nclass)
c      call defi    (iilc,nlnmxc*nclass)
c      call defi    (iinc,nlnmxc*nclass)
c      call defi    (iimc,nlnmxc*nclass)
c      call defi    (iilnmc,nnc*nl*nl*nclass)
      allocate(
     & il(nlnmx,nclass),
     & in(nlnmx,nclass),
     & im(nlnmx,nclass),
     & ilnm(nn*nl*nl*nclass),
     & ilv(nlnmxv*nclass),
     & inv(nlnmxv*nclass),
     & imv(nlnmxv*nclass),
     & ilnmv(nnv*nl*nl*nclass),
     & ilc(nlnmxc*nclass),
     & inc(nlnmxc*nclass),
     & imc(nlnmxc*nclass),
     & ilnmc(nnc*nl*nl*nclass)
     & )

      call idxlnmc ( nindxv,nindxc,
     d              nl,nn,nnv,nnc,nlnmx,nlnmxv,nlnmxc,nclass,
     o              il,in,im,ilnm,
     o              ilv,inv,imv,ilnmv,
     o              ilc,inc,imc,ilnmc)

c      call defi    (inlnmv,nclass)
c      call defi    (inlnmc,nclass)
c      call defi    (inlnm,nclass)
      allocate(nlnmv(nclass),nlnmc(nclass),nlnm(nclass))
      call nolnma  (nindxv,nl,nclass,
     o              nlnmv )
      call nolnma  (nindxc,nl,nclass,
     o              nlnmc )
      call nolnma  (nindx,nl,nclass,
     o              nlnm )
c takao
      i=0
      CALL SKIPBL(symgrp,120,I)
      if(symgrp(i+1:i+13)/= 'UseSYMOPSfile') then
        stop " Not: UseSYMOPSfile in LMTO file"
      endif
      write(6,*) ' symgrp==UseSYMOPSfile'
      ifi = 6661
      open (ifi, file='SYMOPS')
      read(ifi,*) ngrp
      allocate(symgg(3,3,ngrp))
      do ig = 1,ngrp
        read(ifi,*)
        do i=1,3
          read(ifi,"(3d24.16)") symgg(i,1:3,ig)
        enddo
      enddo
      close(ifi)
c      call defdr(igrp,9*ngrp)
c      call dcopy(9*ngrp,symgg,1,w(igrp),1)
c      deallocate(symgg)

c> index for inverse rotations
c      call defi    (iinvg,ngrp)
      allocate(invg(ngrp))
      call invgrp  (symgg,ngrp,
     o              invg)
      is = iclose ('LMTO')
c      is = iclose ('GWIN_V2')
      is = iclose ('ECORE')
      call cputid(0)
      write (*,*) 'genallcf_v3'
      end subroutine
      end module


c-------------------------------------------------------------------
c taken from bas.f
      subroutine gencor (eco,ncwf,iclass,
     d                   nl,nnc,natom,nclass,nctot,
     o                   eclnm)

c 92.03.25
c generate core energies for each l,n,m

c eco(l,n,cl)  = core energies
c ncwf(l,n,cl) = 1 ==> allowed core state
c                0 ==> not allowed
c nnc          = max. no. core states/l
c nctot        = total no. core states (all atoms)

c eclnm({l,n,m,atom}) = core energies

      implicit real*8 (a-h,o-z)
      dimension eco(0:nl-1,nnc,2,nclass),ncwf(0:nl-1,nnc,nclass),
     i          iclass(natom)
      dimension eclnm(nctot,2)

      i          = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        do       l = 0,nl-1
          do       n = 1,nnc
            do       m = -l,l
              if (ncwf(l,n,ic) .eq. 1) then
                i          = i + 1
                if (i .gt. nctot) stop 'gencor: wrong nctot'
cccccccccccccccc
c      write(6,"(' eco=',4i4,d13.5)")i,l,n,ic,eco(l,n,1,ic)
cccccccccccccccc
                eclnm(i,1) = eco(l,n,1,ic)
                eclnm(i,2) = eco(l,n,2,ic)
              endif
            end do
          end do
        end do
      end do
      end

c--------------------------------------------------------------------
      subroutine rcore2 (ifec, ! minor changed version to rcore in rw.f
     d                  nl,nnc,nclass,
     o                  konfig,eco)

c 92.03.17
c read core eigenvalues both up and down

      implicit real*8 (a-h,o-z)
      character*32 char32
      dimension konfig(0:nl-1,nclass),eco(0:nl-1,nnc,2,nclass)
c loop over class
      do      ic = 1,nclass
        write(6,*) ' rcore2 : ic=',ic
        read (ifec,6000) blank
        read (ifec,6000) blank
        read (ifec,6000) blank
        read (ifec,*) z,ic1,nr,a,b,nsp
        read (ifec,6000) blank
        read (ifec,*) (konfig(l,ic),l=0,nl-1)
        read (ifec,6000) blank
c
        write(6,*) ' rcore : goto loop over l n'

c loop over l,n
        do       l = 0,nl-1
          ncore      = konfig(l,ic)-l-1
          if (ncore .gt. nnc) stop 'rcore: wrong nnc'
          do       n = 1,ncore
            read (ifec,*) lt,nt,(eco(l,n,isp,ic),isp=1,nsp) !takao
            if(nsp==1) eco(l,n,2,ic) = eco(l,n,1,ic)        !
cccccccccccccccc
c      write(6,"(' eco=',3i4,d13.5)")l,n,ic,eco(l,n,isp,ic)
cccccccccccccccc
            if (lt .ne. l) stop 'rcore: wrong l'
            if (nt .ne. n) stop 'rcore: wrong n'
          end do
        end do

      end do
 6000 format(a)
 6300 format(1x,2i3,2d13.6)
 6500 format(1x,10i3)
      return
      end










