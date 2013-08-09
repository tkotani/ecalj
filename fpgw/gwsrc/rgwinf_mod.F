      module m_rgwinf_v3
      integer(4)::nclass,natom,nspin,nl,nnv,nnc,lcutmx,nrx
      real(8)::alat
      integer(4),allocatable:: 
     &   iclass(:)
     &  ,nindxv(:,:),nindxc(:,:)
     &  ,occv(:,:,:),unoccv(:,:,:),ooo(:,:,:)
     &  ,occc(:,:,:),unoccc(:,:,:)
     &  ,ncwf(:,:,:)
      real(8),allocatable:: z(:),cutbase(:)
      contains
      subroutine rgwinf_v3 (iflmto,ifinin,incwfx) !these are inputs
C- readin GWIN_V2 and LMTO(crystal) data.
C all the output are given in the declear section
C!
Cr Return iclass=ibas.
Cr nwin,efin,incwfx, are used as switches.
C--------------------------------------------------------
      implicit real*8(a-h,o-z)
c      character(120):: symgrp
c      character(120):: symgrpt
      integer(4):: infwfx,nwin
      logical :: nocore
      integer(4),allocatable::ncwf2(:,:,:),nrofi(:)
      real(8)::efin
      character(6)::clablxxx
      ifi = iflmto
      ef  = -999d0 ! not readin efermi
      nw  = -999 !Not readin NW file

c---------
c SYMMETRY
c---------
c      print *,' goto sym'
c      read(ifi,*);  read(ifi,*)
c      read(ifi,*)symgrpt
c      j           = 0
c      call rmvbl   (symgrpt,120,j)
c      symgrp(1:2) = '  '
c      symgrp(3:120) = symgrpt(j+1:120)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)

c----------
c STRUCTURE
c----------
c      print *,' goto structure'
c> lattice constant
      read(ifi,*);  read(ifi,*); read(ifi,*)
      read(ifi,*)alat

c> primitive lattice vectors
c      allocate(plat(3,3))
      read(ifi,*)
      read(ifi,*)!plat(1:3,1)
      read(ifi,*)!plat(1:3,2)
      read(ifi,*)!plat(1:3,3)

c> no. atoms
      read(ifi,*)
      read(ifi,*)natom

c We assume nclass=natom
      nclass = natom

c> positions of atoms
c      allocate(pos(3,natom))
      read(ifi,*)
      do n = 1,natom
        read(ifi,*) !pos(1,n),pos(2,n),pos(3,n)
      end do

c-----
c LMTO
c-----
      print *,' goto lmto'
c> spin (1=paramagnetic  2=ferromagnetic)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)
      read(ifi,*)nspin
c      print *,' nspin=',nspin

c> max. no. valence and core l
      read(ifi,*)
      read(ifi,*)nl
c      print *,' ispin nl =',ispin,nl

c> max. no. valence and core n
      read(ifi,*)
      read(ifi,*)nnv,nnc

ccccccccccccccccccccccccccc
      if(nnv==1) nnv=2 ! for backword compatibility!takao apr 2002
ccccccccccccccccccccccccccc

c> max. no. radial mesh points
      read(ifi,*)
      read(ifi,*)nrx

c> class-label, z, no. radial points
      read(ifi,*)
      allocate(z(nclass),nrofi(nclass))
      do      ic = 1,nclass
        read(ifi,*)clablxxx,z(ic),nrofi(ic)
      end do

c> atom and its class
      allocate(iclass(natom))
      do n = 1,natom
        iclass(n)=n
      end do

      allocate(nindxv(nl,nclass),nindxc(nl,nclass),
     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
     &        ooo(nl,nnv,nclass),
     &        occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))
      allocate( cutbase(0:2*(nl-1)) )
ctakao
      call rgwinaf    (ifi,ifinin,nl,nnv,nnc,nclass,
c> BZ
     o                  n1,n2,n3,ef,
c> frequencies
     o                  niw,diw,nw,dw,delta,deltaw,esmr,imagw,
c> coulomb
c     o                  tolvc,alp,alptx,h,ng,
c> product basis
     o                  cutbase,lcutmx,nindxv,nindxc,
     o                  occv,unoccv, occc,unoccc,
c> core
     o                  ncwf,ncwf2 )
c----
c      cutbase = tolbas
      inquire(file='NoCore',exist=nocore)
      if(nocore) then
        occc=0    ! call iclear(nl*nnc*nclass, w(ioccc))
        unoccc=0    ! call iclear(nl*nnc*nclass, w(iunoccc))
        ncfw  =0    ! call iclear(nl*nnc*nclass, w(incwf))
      elseif( incwfx==-1 ) then
        print *,' ### incwf=-1 Use ForSxc for core'
        ncwf = ncwf2  !call icopy(nl*nnc*nclass,w(incwf2),w(incwf))
      elseif( incwfx==-2 ) then
        print *,' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
        call notbit(nl*nnc*nclass, ncwf2)
        ncwf  = ncwf2 ! call icopy (nl*nnc*nclass, w(incwf2),w(incwf))
        occc= ncwf  ! call icopy (nl*nnc*nclass, w(incwf),w(ioccc))
        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
cccccccccccccccccccccccccccccccccc
c 31May2006
        ooo=0
        call ibitand(nl*nnv*nclass, unoccv,occv, ooo)
        unoccv = ooo
cccccccccccccccccccccccccccccccccc
      elseif( incwfx==-3 ) then
        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
        print *,' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
      elseif( incwfx==-4 ) then
        print *,' ### incwf=-4  occ=0 and unocc=0 for all core '
        occc=0  !call iclear(nl*nnc*nclass, w(ioccc))
        unoccc=0  !call iclear(nl*nnc*nclass, w(iunoccc))
        ncwf  =0  !call iclear(nl*nnc*nclass, w(incwf))
      elseif(incwfx==0) then
        print *,' ### Use unocc occ ForX0 for core'
      else
Cstop2rx 2013.08.09 kino        stop ' ### proper incwf is not given for genallcf2:rgwinf '
        call rx( ' ### proper incwf is not given for genallcf2:rgwinf ')
      endif
      deallocate(ncwf2)
      end subroutine
      end module

c--------------------------------------------------------------------
      subroutine rgwinaf (ifix,ifinin,nl,nnv,nnc,nclass,
c> BZ
     o                  n1,n2,n3,ef,
c> frequencies
     o                  niw,diw,nw,dw,delta,deltaw,esmr, imagw,
c> product basis
     o                  tolbas,lcutmx,nindxv,nindxc,
     o                  noccv,nunoccv,noccc,nunoccc,
c> core
     o                  ncwf,ncwf2)
      use keyvalue
c read the rest of GW input data
      implicit real*8(a-h,o-z)
      dimension 
     o          nindxv(0:nl-1,nclass),
     o          nindxc(0:nl-1,nclass),
     o          noccv(0:nl-1,nnv,nclass),
     o          noccc(0:nl-1,nnc,nclass),
     o          nunoccv(0:nl-1,nnv,nclass),
     o          nunoccc(0:nl-1,nnc,nclass),
     o          ncwf(0:nl-1,nnc,nclass),
     o          ncwf2(0:nl-1,nnc,nclass),
     o  tolbas(0:2*(nl-1))
c      logical :: readgwinput
      integer(4):: ret ,ncc
      real(8)::tolbas0
      logical::readon
      character*1000 tolchar
c---------------
c BRILLOUIN ZONE
c---------------
      ifi=ifix
c> divisions along reciprocal base vectors G1 G2 G3
      read(ifi,6000)blank
      read(ifi,6000)blank
      read(ifi,6000)blank
      read(ifi,*)n1,n2,n3
      print *,' divisions along reciprocal base vectors G1 G2 G3=',
     & n1,n2,n3

c--- EFERMI
      if(ef/=-999d0) then
        open( 1101,file='EFERMI')
        read( 1101,*) ef
        close(1101)
        print *,' --- READIN ef from EFERMI. ef=',ef
      endif
c--- NW
      if(nw/=-999) then
        open( 1101,file='NW')
        read( 1101,*) nw
        close(1101)
        print *,' --- READIN nw from NW. nw=',nw
      endif

c------------
c FREQUENCIES
c------------
c      if(readgwinput()) then
      call getkeyvalue("GWinput","dw",    dw  )
      call getkeyvalue("GWinput","niw",   niw )
      call getkeyvalue("GWinput","delta", delta )
      call getkeyvalue("GWinput","deltaw",deltaw )
      call getkeyvalue("GWinput","esmr",  esmr )
c      else
c       ifi=ifinin
cc> no. frequencies along Re axis and mesh size (Hartree)
c       read(ifi,6000)blank
c       read(ifi,6000)blank
c       read(ifi,*) dw
cc> no. frequencies along Im axis and mesh size (Hartree)
c       read(ifi,6000)blank
c       read(ifi,*)niw !,diw
cc> broadening of the delta function in the response function (Hartree)
c       read(ifi,6000)blank
c       read(ifi,*) delta
c       read(ifi,6000)blank
c       read(ifi,*) deltaw
c       read(ifi,6000)blank
c       read(ifi,*) esmr
cc> Re (0), Im (1) or both (2) frequencies
c       read(ifi,6000)blank
c       read(ifi,6000)blank
c      endif
      diw=1d90 !dummy never used.
      imagw=-999999 ! dummy
      print *,' --- Freq ---'
      print *,'    dw   =',dw
      print *,'    niw  =',niw
      print *,'    delta=',delta
      print *,'    esmr =',esmr

c--------------
c PRODUCT BASIS
c--------------
ctakao
      ncwf  =99 !This is for counting the number of nctot in gencor.
      ncwf2 =99
      print *,' goto prod'
c      if(readgwinput()) then
      call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifi,status=ret)
      write(6,*)'<... > '
c      else
c       read(ifi,6000)blank
c       read(ifi,6000)blank
c      endif

c> tolerance in percentage for optimal product basis
c  sum of lowest overlap eigenvalues / total eigenvalues
      read(ifi,6000)blank

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
          read(tolchar(ixoff:ix),*,err=1097) tolbas(lx)
          if(lx==2*(nl-1)) goto 1098
          readon=.false.
          lx=lx+1
        endif
      enddo
 1097 continue
      tolbas(lx:)=tolbas(lx-1)
 1098 continue
      do lx=0,2*(nl-1)
        write(6,"(' lx=',i3,' readin tolerance=',d11.3)") lx, tolbas(lx)
      enddo
ccccccccccccccccccccccccccccc
c      read(ifi,*)tolbas0
c      tolbas(:)=tolbas0
ccccccccccccccccccccccccccccc

c> check completeness of optimal product basis
      read(ifi,6000)blank
      read(ifi,*)lcutmx
      write(6,"(' --- prod section: lcutmx tolbas='i3,100d11.3)") lcutmx,tolbas
c      print *,' nclass nl-1=',nclass,nl-1
c> n index for valence and core
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          read(ifi,*) ict,lt,nindxv(l,ic),nindxc(l,ic)
ccccccccc
c      write(6,*) ict,ic,lt,nindxv(l,ic),nindxc(l,ic)
cccccccccccccc
c      if(ict .ne. ic) then
c       stop 'rgwina: 1st wrong class '
c      endif
Cstop2rx 2013.08.09 kino          if(lt  .ne. l ) stop 'rgwina: 1st wrong l '
          if(lt  .ne. l ) call rx( 'rgwina: 1st wrong l ')
        end do
      end do
      print *,' --- valence product basis section'
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
c> criteria for product basis
c>> valence

ccccccccccccccccccccccccccc
c      noccv=0;nunoccv=0
c      goto 1099
ccccccccccccccccccccccccccccccc
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          do       n = 1,nindxv(l,ic)
            read(ifi,*)ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
            write(6,"(100i3)") ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
c      if(ict .ne. ic)stop 'rgwina: 2nd wrong class valence'
Cstop2rx 2013.08.09 kino            if(lt  .ne. l )stop 'rgwina: 2nd wrong l valence'
            if(lt  .ne. l )call rx( 'rgwina: 2nd wrong l valence')
Cstop2rx 2013.08.09 kino            if(nt  .ne. n )stop 'rgwina: wrong n valence'
            if(nt  .ne. n )call rx( 'rgwina: wrong n valence')
          end do
        end do
      end do
 1099 continue
c      print *,' goto prod 3'
c>> core
      print *,' --- core product basis section'
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            read(ifi,*)ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic)
     & ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
            write(6,"(100i3)") 
     & ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic)
     & ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
c      if(ict .ne. ic)stop 'rgwina: 2nd wrong class core'
Cstop2rx 2013.08.09 kino            if(lt  .ne. l )stop 'rgwina: 2nd wrong l core'
            if(lt  .ne. l )call rx( 'rgwina: 2nd wrong l core')
Cstop2rx 2013.08.09 kino            if(nt  .ne. n )stop 'rgwina: wrong n core'
            if(nt  .ne. n )call rx( 'rgwina: wrong n core')
          end do
        end do
      end do
c      print *,' goto prod 4'
 6000 format(a)
      close(ifi)
      return
      end
c

      subroutine bit99to0(n,idat)
      integer(4) :: idat(n)
      do i=1,n
        if(idat(i)==99) idat(i)=0
      enddo
      end
c
      subroutine notbit(n,idat)
      integer(4) :: idat(n)
      do i=1,n
        ix =  idat(i)
        if(idat(i)==0) then
          idat(i)=1
        elseif(idat(i)==1) then
          idat(i)=0
        endif
c       print *,'notbit=',i,ix,idat(i)
      enddo
      end
      subroutine iclear(n,idat)
      integer(4) :: idat(n)
      idat=0
      end
      subroutine ibiton(nclass,nl,nnc,nindxc, noccc,ncwf)
      integer(4) ::noccc(0:nl-1,nnc,nclass),nindxc(0:nl-1,nclass)
      integer(4) ::ncwf(0:nl-1,nnc,nclass)
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

      subroutine ibitand(n,a,b,c)
      integer(4) :: n, a(n),b(n),c(n),i
      do i = 1,n
        if( a(i)==1.or.b(i)==1 ) c(i)=1
        write (6,"('ibitand:: ',4i3)")i,a(i),b(i),c(i)
      enddo
      end

      subroutine writeemesh(ifi,freqi,niw,freqr,nnw,delta)
c Write energy mesh along imag axis and real axis. -----------------------
      implicit none
      integer(4):: iw,ifi,niw,nnw
      real(8) :: freqi(niw),freqr(0:nnw-1),delta
      complex(8):: fff,img =(0d0,1d0)
      write(6,*)" writeemesh: ifi=",ifi
      write(ifi,"(' iw   omega(Ry) on Imag-axis  niw=', i8)")niw
      do iw = 1,niw
        fff = img * freqi(iw)            ! along img axis
        write(ifi,"(i3,2d15.6)")iw,fff*2d0
      enddo
      write(ifi,"('  0   0.000000D+00   0.000000D+00  ! This line is a dummy not in niw.')")
      write(ifi,"(' iw   omega(Ry) on Real-axis  nw=', i8)") nnw
      do iw= 0,nnw-1
        fff = freqr(iw) + img*delta ! delta is in a.u.  ! along real axis
        write(ifi,"(i3,2d15.6)")iw,fff*2d0
      enddo
      end






