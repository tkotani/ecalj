      module m_genallcf_v3
!!-- Read basis data ----------------------------------
!!! here is old memo      
!! - structure
!!  -  o                   plat,alat,natom,nclass,pos,
!!  -  o                   ngrp, symgg,
!!  -  o                   invg, ef,
!! - l,n and dimensions
!!   - o                   clabl, nspin,nl,nn,nnv,nnc,
!!   - o                   nindx, nindxv, nindxc, iclass,
!!   - d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
!!   - o                   z,
!! - l,n,m indices for Phi (atomic basis)
!!   -  o                   il, in, im,  ilnm,  nlnm,
!!   -  o                   ilv,inv,imv, ilnmv, nlnmv,
!!   -  o                   ilc,inc,imc, ilnmc, nlnmc,
!! - core
!!   -  o                   ncwf, ecore, konf, icore, ncore,nctot,
!! - frequency
!!   -                   niw,diw,nw,dw,delta,deltaw,esmr, freq)
!!             symgrp
!!             ,nocc, nunocc, occv, unoccv, occc, unoccc
      implicit none
      public:: Setesmr,Genallcf_v3
      integer,protected,public:: nrx,lcutmx
      real(8),allocatable,public::cutbase(:)
      integer,allocatable,protected,public::
     &     iclass(:),
     &     nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
     o     invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
     o     ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
     o     ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
     o     nindx(:,:),konf(:,:),icore(:,:),ncore(:),
     &     occv(:,:,:),unoccv(:,:,:)   
     &     ,occc(:,:,:),unoccc(:,:,:),
     o     nocc(:,:,:),nunocc(:,:,:),   iantiferro(:),ooo(:,:,:)
      integer,protected,public::
     o     nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     o     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot, niw 
      real(8),protected,public::
     &     alat,deltaw,esmr,delta,tpioa
      real(8), allocatable,protected,public:: 
     &     plat(:,:),pos(:,:),z(:),symgg(:,:,:),ecore(:,:)
      logical,protected,public::
     &     done_genallcf_v3=.false.
      character(8),allocatable,protected,public::
     &     spid(:)
      character(120),protected,public:: symgrp
      character(6),allocatable,protected,public :: clabl(:)
      
      private
c-----------------------------------------------------------------------
      contains

      subroutine setesmr(esmr_in)
      intent(in)::       esmr_in
      real(8):: esmr_in
      esmr=esmr_in
      end
      
      subroutine genallcf_v3(incwfx)
      use m_keyvalue,only: getkeyvalue
      implicit none
      intent(in)::           incwfx
!! Readin GWIN_V2 and LMTO(crystal) data and allocate all required.
!! Return iclass=ibas.
!! efin,incwfx, are used as switches.
!! input: efin,incwfx,
!!        GWIN_V2, LMTO
!! output: All the output are given in the declear section above.
!! --------------------------------------------------------
      integer::incwfx,ifec,i,j,
     & lmx, lmx2,nlmto2,nprodxc,nlnaxc,nlnaxv,nprodx,ifi,ig,is
     & ,nprodxv,nlnax,ix,ixoff,lx
     & ,noflmto,maxnn
      integer:: infwfx,ret, n1,n2,n3,imagw,n,ic
      logical :: nocore,readon
      real(8)::efin
      character(120):: symgrpt
      character*1000 tolchar
      real(8),   allocatable:: ecoret(:,:,:,:)
      integer(4),allocatable::ncwf2(:,:,:)
      integer:: ia,l,m,ic1,isp,lt,nt,nsp,nr,ncorex,ifix
      real(8)::a,b,zz, efdummy,dw,diw,pi
      integer:: nwdummy,ifile_handle
      if(done_genallcf_v3) call rx('genallcf_v3 is already called')
      done_genallcf_v3=.true.
c--- rgwinf ---
      open(newunit=ifi,file='LMTO')
      read(ifi,*);  read(ifi,*)
      read(ifi,*)symgrpt        !SYMMETRY
      j = 0
      symgrp='  '//trim(adjustl(symgrpt))
      write(6,*)' symgrp=', symgrp
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)alat           !lattice constant
      pi=4d0*datan(1d0)
      tpioa=2d0*pi/alat
      
      allocate(plat(3,3))          !primitive lattice vectors
      read(ifi,*)
      read(ifi,*)plat(1:3,1)
      read(ifi,*)plat(1:3,2)
      read(ifi,*)plat(1:3,3)
      read(ifi,*)
      read(ifi,*) natom            !Number of atoms
!!
      nclass = natom  !We set nclass = natom through the GW calculations
      write(6,*)'genalloc: alat natom=',alat,natom
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
      write(6,*)' nspin nl nnv nnc =',nspin,nl,nnv,nnc
c--------------------------
      if(nnv==1) nnv=2 ! for backword compatibility!takao apr 2002
      ! nnv=2 corresponds to phi and phidot
      ! nnv=3 corresponds to
c--------------------------
      read(ifi,*)
      read(ifi,*)  nrx !is not readin
      read(ifi,*)
      allocate(clabl(nclass),z(nclass)) !class-label, z
      do ic = 1,nclass
        read(ifi,*) clabl(ic),z(ic) !,nrofi is not readin
      end do
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)n1,n2,n3
      write(6,"(a,3i5)")' divisions along reciprocal base vectors G1 G2 G3=', n1,n2,n3
      allocate(iantiferro(1:natom),spid(1:natom))
      read(ifi,*)
      read(ifi,*)iantiferro(1:natom) !may2015
      read(ifi,*)
      read(ifi,*)spid(1:natom)
      close(ifi)
      
c------------
c FREQUENCIES
c------------
      call getkeyvalue("GWinput","niw",   niw )
      call getkeyvalue("GWinput","delta", delta )
      call getkeyvalue("GWinput","deltaw",deltaw )
      call getkeyvalue("GWinput","esmr",  esmr )
      imagw=-999999 ! dummy
      write(6,*)' --- Freq ---'
      write(6,"(a,i6)")'    niw  =',niw
      write(6,"(a,f12.6)")'    delta=',delta
      write(6,"(a,f12.6)")'    esmr =',esmr
!
      allocate(iclass(natom)) !atom and its class.
      do n = 1,natom          !!We set nclass = natom through the GW calculations
        iclass(n)=n
      end do

      allocate(nindxv(nl,nclass), nindxc(nl,nclass),
     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
     &     occc(nl,nnc,nclass),unoccc(nl,nnc,nclass),
     &        ooo(nl,nnv,nclass))
      
      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))
      allocate(cutbase(0:2*(nl-1)))
      
c--------------
c PRODUCT BASIS
c--------------
      ncwf  =99 !This is for counting the number of nctot in gencor.
      ncwf2 =99
      write(6,*)' reading <PRODUCT_BASIS> section'
      call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifi,status=ret) !this open GWinput and locate file position.c> tolerance in percentage for optimal product basis
c  sum of lowest overlap eigenvalues / total eigenvalues
      read(ifi,*)
c... readin tlerance for each l jun2004
      read(ifi,"(a)") tolchar
      readon=.false.
      lx=0
      do ix=1,1000
        if(.not.readon.and.tolchar(ix:ix)/=' ') then
          readon=.true.
          ixoff=ix
        endif
        if(readon.and.tolchar(ix:ix)==' ') then
          read(tolchar(ixoff:ix),*,err=1097) cutbase(lx)
          if(lx==2*(nl-1)) goto 1098
          readon=.false.
          lx=lx+1
        endif
      enddo
 1097 continue
      cutbase(lx:)=cutbase(lx-1)
 1098 continue
      do lx=0,2*(nl-1)
        write(6,"(' lx=',i3,' readin tolerance=',d11.3)") lx, cutbase(lx)
      enddo
c> check completeness of optimal product basis
      read(ifi,*)
      read(ifi,*)lcutmx
      write(6,"(' --- prod section: lcutmx cutbase='i3,100d11.3)") lcutmx,cutbase
      
      call rgwinaf    (ifi,nl,nnv,nnc,nclass, !ifi can be changed.
c> BZ
     o                  n1,n2,n3,efdummy,
c> frequencies
     o                  niw,diw,nwdummy,dw,delta,deltaw,esmr,imagw,
c> product basis
     o                  cutbase,lcutmx,nindxv,nindxc,
     o                  occv,unoccv, occc,unoccc,
c> core
     o                  ncwf,ncwf2 )
c----
      inquire(file='NoCore',exist=nocore)
      if(nocore) then
        occc=0      
        unoccc=0    
        ncwf  =0    
      elseif( incwfx==-1 ) then
        write(6,*)' ### incwf=-1 Use ForSxc for core'
        ncwf = ncwf2 
      elseif( incwfx==-2 ) then
        write(6,*)' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
        call notbit(nl*nnc*nclass, ncwf2)
        ncwf  = ncwf2
        occc= ncwf  
        unoccc= 0   
cccccccccccccccccccccccccccccccccc
        ooo=0
        call ibitand(nl*nnv*nclass, unoccv,occv, ooo)
        unoccv = ooo
cccccccccccccccccccccccccccccccccc
      elseif( incwfx==-3 ) then
        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
        unoccc= 0   
        write(6,*)' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
      elseif( incwfx==-4 ) then
        write(6,*)' ### incwf=-4  occ=0 and unocc=0 for all core '
        occc=0  
        unoccc=0
        ncwf=0  
      elseif(incwfx==0) then
        write(6,*)' ### Use unocc occ ForX0 for core'
      else
        call rx( ' ### proper incwf is not given for genallcf_v3:rgwinf ')
      endif
      deallocate(ncwf2)
C...  End of rgwinf section ----------------------------
      
!! dimensions and constants
      lmx        = 2*(nl-1)
      lmx2       = (lmx+1)**2
      nlmto      = noflmto(nindxv,iclass,nl,nclass,natom)
      nlmto2     = nlmto*nlmto
      nn         = maxnn (nindxv,nindxc,nl,nclass)

!! combine nocc,nunocc,nindx
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

!! index for allowed core states
      allocate(icore(nl**2*nnc,nclass),ncore(nclass))
      icore=9999999
      ncore=9999999
      call incor   (ncwf,nindxc,iclass,
     d              nl,nnc,nclass,natom,
     o              icore,ncore,nctot )

!! core energies
      open(newunit=ifec,file='ECORE')
      allocate(konf(nl,nclass),ecore(nctot,2))
      konf=0
      allocate(ecoret(0:nl-1,nnc,2,nclass))
      ecoret=0d0
      do ic = 1,nclass
        write(6,*) ' read ECORE : ic=',ic
        read (ifec,*) 
        read (ifec,*) 
        read (ifec,*) 
        read (ifec,*) !zz,ic1,nr ,a,b,nsp
        read (ifec,*) 
        read (ifec,*) (konf(l+1,ic),l=0,nl-1)
        read (ifec,*) 
        do  l = 0,nl-1
          ncorex = konf(l+1,ic)-l-1
          if (ncorex .gt. nnc) call rx( 'ECORE: wrong nnc')
          do n = 1,ncorex
            read (ifec,*) lt,nt,(ecoret(l,n,isp,ic),isp=1,nspin) !takao
            if(nspin==1) ecoret(l,n,2,ic) = ecoret(l,n,1,ic)        !
c           write(6,"(' read ecore=',3i4,2d13.5)")l,n,ic,ecoret(l,n,1:nspin,ic)
            if (lt .ne. l) call rx( 'rcore: wrong l')
            if (nt .ne. n) call rx( 'rcore: wrong n')
          end do
        end do
      end do
      close(ifec)
      i = 0
      do ia = 1,nclass
        ic  = iclass(ia)
        do l = 0,nl-1
        do n = 1,nnc
        do m = -l,l
          if (ncwf(l+1,n,ic) .eq. 1) then
            i = i + 1
            if (i > nctot) call rx( 'genalloc_mod: wrong nctot')
            ecore(i,1:nspin) = ecoret(l,n,1:nspin,ic)
            write(6,"(' ecore=',4i4,2d13.5)")i, l,n,ic,ecore(i,1:nspin)
          endif
        enddo
        enddo
        enddo
      enddo
!! dummy to overlaid -check bounds sep2014 jun2020 moved from sxcf
      if(size(ecore)==0) then
         deallocate(ecore)
         allocate(ecore(1,2))
      endif
      deallocate(ecoret)
!! index for core and LMTO basis
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
      allocate(nlnmv(nclass),nlnmc(nclass),nlnm(nclass))
      call nolnma  (nindxv,nl,nclass,
     o              nlnmv )
      call nolnma  (nindxc,nl,nclass,
     o              nlnmc )
      call nolnma  (nindx,nl,nclass,
     o              nlnm )
      i=2 !see previous definition of symgrp
      if(symgrp(i+1:i+13)/= 'UseSYMOPSfile') then
        call rx( " Not: UseSYMOPSfile in LMTO file")
      endif
      write(6,*) ' symgrp==UseSYMOPSfile'
      open (newunit=ifi, file='SYMOPS')
      read(ifi,*) ngrp
      allocate(symgg(3,3,ngrp))
      do ig = 1,ngrp
        read(ifi,*)
        do i=1,3
          read(ifi,"(3d24.16)") symgg(i,1:3,ig)
        enddo
      enddo
      allocate(invg(ngrp))
      call invgrp  (symgg,ngrp,
     o              invg)
      close(ifi)
      call cputid(0)
      write(6,*) 'genallcf_v3'
      end subroutine genallcf_v3
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine bit99to0(n,idat)
      integer(4) :: idat(n),i,n
      do i=1,n
        if(idat(i)==99) idat(i)=0
      enddo
      end
c
      subroutine notbit(n,idat)
      integer(4) :: idat(n),i,n,ix
      do i=1,n
        ix =  idat(i)
        if(idat(i)==0) then
          idat(i)=1
        elseif(idat(i)==1) then
          idat(i)=0
        endif
c       write(6,*)'notbit=',i,ix,idat(i)
      enddo
      end
      subroutine iclear(n,idat)
      integer(4) :: idat(n),n
      idat=0
      end
      subroutine ibiton(nclass,nl,nnc,nindxc, noccc,ncwf)
      integer(4) ::nclass,nl,nnc,noccc(0:nl-1,nnc,nclass),nindxc(0:nl-1,nclass)
      integer(4) ::ncwf(0:nl-1,nnc,nclass),ic,l,n
      noccc=0
      ncwf=0
      do      ic = 1,nclass
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            noccc(l,n,ic)=1
            ncwf(l,n,ic) =1
          end do
        end do
      end do
      end
      
!ssssssssssssssssssssssssssssssssssssssssssss
      subroutine ibitand(n,a,b,c)
      integer(4) :: n, a(n),b(n),c(n),i
      do i = 1,n
        if( a(i)==1.or.b(i)==1 ) c(i)=1
        write (6,"('ibitand:: ',4i3)")i,a(i),b(i),c(i)
      enddo
      end

c$$$      subroutine writeemesh(ifi,freqi,niw,freqr,nnw,delta)
c$$$c Write energy mesh along imag axis and real axis. -----------------------
c$$$      implicit none
c$$$      integer(4):: iw,ifi,niw,nnw
c$$$      real(8) :: freqi(niw),freqr(0:nnw-1),delta
c$$$      complex(8):: fff,img =(0d0,1d0)
c$$$      write(6,*)" writeemesh: ifi=",ifi
c$$$      write(ifi,"(' iw   omega(Ry) on Imag-axis  niw=', i8)")niw
c$$$      do iw = 1,niw
c$$$        fff = img * freqi(iw)            ! along img axis
c$$$        write(ifi,"(i3,2d15.6)")iw,fff*2d0
c$$$      enddo
c$$$      write(ifi,"('  0   0.000000D+00   0.000000D+00  ! This line is a dummy not in niw.')")
c$$$      write(ifi,"(' iw   omega(Ry) on Real-axis  nw=', i8)") nnw
c$$$      do iw= 0,nnw-1
c$$$        fff = freqr(iw) + img*delta ! delta is in a.u.  ! along real axis
c$$$        write(ifi,"(i3,2d15.6)")iw,fff*2d0
c$$$      enddo
c$$$      end
      
      end module m_genallcf_v3
      
!!
      subroutine idxlnmc(nindxv,nindxc,
     d                   nl,nn,nnv,nnc,nlnmx,nlnmxv,nlnmxc,nclass,
     o                   il,in,im,ilnm,
     o                   ilv,inv,imv,ilnmv,
     o                   ilc,inc,imc,ilnmc)
c 92.jan.07
c 92.03.17 include core states
c indexing of core states and LMTO basis functions for all classes,
c follows that in TB-LMTO program
c il,in,im = l,n,m
c ilnm(n,lm) = index of n,l,m
c lm = l*l + l + m + 1
c NOTE: the indexing starts with core first and then valence on top
c       of core (not the same as index generated from nindx)
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension nindxv(0:nl-1,nclass),nindxc(0:nl-1,nclass)
      dimension ilnm(nn,nl*nl,nclass),
     o          ilnmv(nnv,nl*nl,nclass),
     o          ilnmc(nnc,nl*nl,nclass),
     o          in(nlnmx,nclass),il(nlnmx,nclass),im(nlnmx,nclass),
     o         inv(nlnmxv,nclass),ilv(nlnmxv,nclass),imv(nlnmxv,nclass),
     o         inc(nlnmxc,nclass),ilc(nlnmxc,nclass),imc(nlnmxc,nclass)
      do     ic = 1,nclass
        ind       = 0
c core
        do      l = 0,nl-1
          l2        = l*l
          do      n = 1,nindxc(l,ic)
            do      m = 1,2*l+1
              ind       = ind + 1
              if (ind .gt. nlnmx) call rx( 'idxlnmc: ind > nlnmx')
              lm        = l2 + m
              il(ind,ic)= l
              in(ind,ic)= n
              im(ind,ic)= m - l - 1
              ilnm(n,lm,ic) = ind
              ilc(ind,ic)= l
              inc(ind,ic)= n
              imc(ind,ic)= m - l - 1
              ilnmc(n,lm,ic)= ind
            end do
          end do
        end do
c valence
        indv      = 0
        do      l = 0,nl-1
          l2        = l*l
          ncore     = nindxc(l,ic)
          do      n = 1,nindxv(l,ic)
            if (ncore+n .gt. nn) call rx( 'idxlnmc: ncore+n > nn')
            do      m = 1,2*l+1
              ind       = ind + 1
              indv      = indv + 1
              if (ind .gt. nlnmx) call rx( 'idxlnmc: ind > nlnmx')
              lm        = l2 + m
              il(ind,ic)= l
              in(ind,ic)= ncore + n
              im(ind,ic)= m - l - 1
              ilnm(ncore+n,lm,ic) = ind
              ilv(indv,ic)= l
              inv(indv,ic)= n
              imv(indv,ic)= m - l - 1
              ilnmv(n,lm,ic) = indv
            end do
          end do
        end do
      end do
      return
      end

      integer function noflmto(nindx,iclass,nl,nclass,natom)
c total number of LMTO basis functions
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension nindx(0:nl-1,nclass),iclass(natom)
      noflmto   = 0
      do  i = 1,natom
        ic        = iclass(i)
      do  l = 0,nl-1
         noflmto   = noflmto + (2*l+1)*nindx(l,ic)
      enddo
      enddo
      return
      end

      integer function nalwln (nocc,nunocc,nindx,nl,nn)
c gives the number of allowed product radial phi
c nocc(l,n)   = 0,1 ==> unoccupied, occupied
c nunocc(l,n) = 1,0 ==> unoccupied,occupied
c nalwln    = number of allowed phi(l1,n1) phi(l2,n2)
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)
      if (nl-1 .gt. lmax) call rx( 'nalwln: increase lmax')
      if (nn .gt. nnx) call rx( 'nalwln: increase nnx')
      icheck=0
      nalwln     = 0
      do 101   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
        if(nocc(l1,n1) .eq. 0)goto 10
        do 201   l2 = 0,nl-1
        do 20   n2 = 1,nindx(l2)
          if(nunocc(l2,n2) .eq. 0)goto 20
          if((l1.ne.l2 .or. n1.ne.n2) .and. icheck(l2,n2,l1,n1).ne.0)
     . goto 20
          nalwln     = nalwln + 1
          icheck(l1,n1,l2,n2) = nalwln
   20   continue
 201  continue
   10 continue
 101  continue
      return
      end

      integer function nofln(nindx,nl)
c count the number of l,n
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension nindx(0:nl-1)
      nofln      = 0
      do       l = 0,nl-1
        nofln      = nofln + nindx(l)
      end do
      return
      end
c------------------------------------------------------------------
      integer function noflnm(nindx,nl)
c number of l,n,m
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension nindx(0:nl-1)
      noflnm    = 0
      do 1    l = 0,nl-1
        noflnm    = noflnm + nindx(l)*(2*l+1)
    1 continue
      return
      end

      integer function nallow (nocc,nunocc,nindx,nl,nn)
c gives the number of allowed product basis
c nocc(n,l) = 0,1 ==> unoccupied, occupied
c nallow    = number of allowed product basis
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)
      if(nl-1 .gt. lmax) call rx( 'nallow: increase lmax')
      if(nn .gt. nnx) call rx( 'nallow: increase nnx')
      icheck=0
      do      l1 = 0,nl-1
        do      n1 = 1,nindx(l1)
          do      l2 = 0,nl-1
            do      n2 = 1,nindx(l2)
              icheck(l1,n1,l2,n2) = nocc(l1,n1)*nunocc(l2,n2)
              if (l1 .ne. l2 .or. n1 .ne. n2) then
                if (icheck(l1,n1,l2,n2)*icheck(l2,n2,l1,n1) .ne. 0)
     . icheck(l1,n1,l2,n2) = 0
              endif
            end do
          end do
        end do
      end do
      nallow     = 0
      do    l1 = 0,nl-1
      do    n1 = 1,nindx(l1)
      do    m1 = 1,2*l1+1
      do    l2 = 0,nl-1
      do    n2 = 1,nindx(l2)
      do    m2 = 1,2*l2+1
c     if (nocc(l1,n1) .eq. 0)goto 10
c     if (nunocc(l2,n2) .eq. 0)goto 10
        if (icheck(l1,n1,l2,n2) .eq. 0) cycle
c temporary
        if (l1 .eq. l2 .and. n1.eq.n2 .and. m1.lt.m2) cycle
        nallow     = nallow + 1
      enddo  
      enddo  
      enddo  
      enddo  
      enddo  
      enddo  
      return
      end

      subroutine incor  (ncwf,nindxc,iclass,
     d                   nl,nnc,nclass,natom,
     o                   icore,ncore,nctot)
c 92.03.18
c sorts out allowed core states and count the number of core states
c ncwf(l,n,cl) = 1 ==> allowed, 0 ==> not allowed
c nindxc(l,cl)  = no. core states/l,class
c nl,nnc = max. no. l,n
c icore(i,cl) = index for allowed core states
c ncore(cl)   = no. allowed core states
c nctot       = total no. allowed core states
      implicit real*8 (a-h,o-z)
      implicit integer(i-n)
      dimension ncwf(0:nl-1,nnc,nclass),nindxc(0:nl-1,nclass),
     i          iclass(natom)
      dimension icore(nl*nl*nnc,nclass),ncore(nclass)
      ncx        = nl*nl*nnc
      do      ic = 1,nclass
        i          = 0
        j          = 0
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            do       m = -l,l
              j          = j + 1
              if (ncwf(l,n,ic) .eq. 1) then
                i          = i + 1
                if (i .gt. ncx) call rx( 'incore: wrong ncx')
                icore(i,ic)= j
              endif
            end do
          end do
        end do
        ncore(ic)  = i
      end do
c total no. allowed core states
      nctot      = 0
      do       i = 1,natom
        ic         = iclass(i)
        nctot      = nctot + ncore(ic)
      end do
      return
      end

      
!------------------------    
      module m_readhbe
      integer,protected:: nprecb,mrecb,mrece,nlmtot,nqbzt,mrecg,nband
c      integer:: nband !not yet protected!
      contains
      subroutine readhbe()
      integer:: ifhbe
      open(newunit=ifhbe, file='hbe.d', action='read')
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      close(ifhbe)
      end subroutine
      end module
c----------------------------------
 
!------------------------
      module m_ReadEfermi
      real(8),protected:: bandgap, ef, ef_kbt
      contains
      subroutine readefermi()
      implicit none
      integer:: ifief,ifile_handle
      ifief=ifile_handle()
      open(ifief,file='EFERMI')
      read(ifief,*) ef,bandgap
      close(ifief)
      write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
      end subroutine
!------------------------------------
      subroutine readefermi_kbt()
      implicit none
      integer:: ifief_kbt,ifile_handle
      ifief_kbt=ifile_handle()
      open(ifief_kbt,file='EFERMI_kbt')
      read(ifief_kbt,*) ef_kbt,bandgap
      close(ifief_kbt)
      write(6,"(a,f12.6)")' --- READIN ef from EFERMI_kbt. ef=',ef_kbt
      end subroutine
      end module m_ReadEfermi
   
      
c----------------------------------------------------------
!! we neede module reading m_genallfc_v3 should be here after it is dedined in this file
      module m_readclasst
      use m_genallcf_v3,only: natom
      integer,allocatable,protected :: iclasst(:)!, invgx(:)
      contains
      subroutine Readclasst()
      integer:: ibas,ibasx,ificlass
      allocate(iclasst(natom))
      write(6,*)'  --- Read true CLASS (crystalographyically equivalent sites) ---'
      open(newunit=ificlass, file='CLASS', action='read')
      do ibas = 1,natom
        read(ificlass,*)  ibasx, iclasst(ibas)
        write(6,"(2i10)") ibasx, iclasst(ibas)
      enddo
      close(ificlass)
      end subroutine
      end module m_readclasst

