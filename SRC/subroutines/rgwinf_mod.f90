! sssssssssssssssssssssssssssssssssssssssssss
subroutine rgwinaf (ifi,nl,nnv,nnc,nclass, &
     !> BZ
     n1,n2,n3,ef, &
     !> frequencies
     niw,diw,nw,dw,delta,deltaw,esmr, imagw, &
     !> product basis
     cutbase,lcutmx,nindxv,nindxc, &
     noccv,nunoccv,noccc,nunoccc, &
     !> core
     ncwf,ncwf2)
  use m_keyvalue,only: getkeyvalue
  ! read the rest of GW input data
  implicit real*8(a-h,o-z)
  implicit integer(i-n)
  dimension &
       nindxv(0:nl-1,nclass), &
       nindxc(0:nl-1,nclass), &
       noccv(0:nl-1,nnv,nclass), &
       noccc(0:nl-1,nnc,nclass), &
       nunoccv(0:nl-1,nnv,nclass), &
       nunoccc(0:nl-1,nnc,nclass), &
       ncwf(0:nl-1,nnc,nclass), &
       ncwf2(0:nl-1,nnc,nclass), &
       cutbase(0:2*(nl-1))
  !      logical :: readgwinput
  integer(4):: ret ,ncc
  !      real(8)::cutbase0
  logical::readon
  !      write(6,*)' nclass nl-1=',nclass,nl-1
  !> n index for valence and core
  read(ifi,6000)blank
  do      ic = 1,nclass
     do       l = 0,nl-1
        read(ifi,*) ict,lt,nindxv(l,ic),nindxc(l,ic)
        if(lt  /= l ) call rx( 'rgwina: 1st wrong l ')
     end do
  end do
  write(6,*)' --- valence product basis section'
  !> criteria for product basis
  !>> valence
  read(ifi,6000)blank
  do      ic = 1,nclass
     do       l = 0,nl-1
        do       n = 1,nindxv(l,ic)
           read(ifi,*)ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
           write(6,"(100i3)") ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
           if(lt  /= l )call rx( 'rgwina: 2nd wrong l valence')
           if(nt  /= n )call rx( 'rgwina: wrong n valence')
        end do
     end do
  end do
1099 continue
  !      write(6,*)' goto prod 3'
  !>> core
  write(6,*)' --- core product basis section'
  read(ifi,6000)blank
  do      ic = 1,nclass
     do       l = 0,nl-1
        do       n = 1,nindxc(l,ic)
           read(ifi,*)ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic) &
                ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
           write(6,"(100i3)") &
                ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic) &
                ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
           if(lt  /= l )call rx( 'rgwina: 2nd wrong l core')
           if(nt  /= n )call rx( 'rgwina: wrong n core')
        end do
     end do
  end do
6000 format(a)
  close(ifi)
  return
end subroutine rgwinaf

!$$$      module m_rgwinf_v3
!$$$      integer,protected:: nclass,natom,nspin,nl,nnv,nnc,lcutmx,nrx
!$$$      real(8),protected:: alat
!$$$      integer,allocatable,protected::
!$$$     &   iclass(:)
!$$$     &  ,nindxv(:,:),nindxc(:,:)
!$$$     &  ,occv(:,:,:),unoccv(:,:,:),ooo(:,:,:)
!$$$     &  ,occc(:,:,:),unoccc(:,:,:)
!$$$     &  ,ncwf(:,:,:)
!$$$      real(8),allocatable,protected:: z(:),cutbase(:)
!$$$
!$$$      contains
!$$$
!$$$      subroutine rgwinf_v3 (iflmto,incwfx)
!$$$      implicit real*8(a-h,o-z)
!$$$      implicit integer(i-n)
!$$$      intent(in)::          iflmto,incwfx
!$$$C- readin GWIN_V2 and LMTO(crystal) data.
!$$$C all the output are given in the declear section
!$$$C!
!$$$Cr Return iclass=ibas.
!$$$Cr nwin,efin,incwfx, are used as switches.
!$$$C--------------------------------------------------------
!$$$c      character(120):: symgrp
!$$$c      character(120):: symgrpt
!$$$      integer(4):: infwfx,nwin,iflmto,ifiniin,incwfx
!$$$      logical :: nocore
!$$$      integer(4),allocatable::ncwf2(:,:,:),nrofi(:)
!$$$      real(8)::efin
!$$$      character(6)::clablxxx
!$$$      ifi = iflmto
!$$$      ef  = -999d0 ! not readin efermi
!$$$      nw  = -999 !Not readin NW file
!$$$
!$$$c---------
!$$$c SYMMETRY
!$$$c---------
!$$$c      write(6,*)' goto sym'
!$$$c      read(ifi,*);  read(ifi,*)
!$$$c      read(ifi,*)symgrpt
!$$$c      j           = 0
!$$$c      call rmvbl   (symgrpt,120,j)
!$$$c      symgrp(1:2) = '  '
!$$$c      symgrp(3:120) = symgrpt(j+1:120)
!$$$      read(ifi,*)
!$$$      read(ifi,*)
!$$$      read(ifi,*)
!$$$
!$$$c----------
!$$$c STRUCTURE
!$$$c----------
!$$$c      write(6,*)' goto structure'
!$$$c> lattice constant
!$$$      read(ifi,*);  read(ifi,*); read(ifi,*)
!$$$      read(ifi,*)alat
!$$$
!$$$c> primitive lattice vectors
!$$$c      allocate(plat(3,3))
!$$$      read(ifi,*)
!$$$      read(ifi,*)!plat(1:3,1)
!$$$      read(ifi,*)!plat(1:3,2)
!$$$      read(ifi,*)!plat(1:3,3)
!$$$
!$$$c> no. atoms
!$$$      read(ifi,*)
!$$$      read(ifi,*)natom
!$$$
!$$$c We assume nclass=natom
!$$$      nclass = natom
!$$$
!$$$c> positions of atoms
!$$$c      allocate(pos(3,natom))
!$$$      read(ifi,*)
!$$$      do n = 1,natom
!$$$        read(ifi,*) !pos(1,n),pos(2,n),pos(3,n)
!$$$      end do
!$$$
!$$$c-----
!$$$c LMTO
!$$$c-----
!$$$      write(6,*)' goto lmto'
!$$$c> spin (1=paramagnetic  2=ferromagnetic)
!$$$      read(ifi,*)
!$$$      read(ifi,*)
!$$$      read(ifi,*)
!$$$      read(ifi,*)nspin
!$$$c      write(6,*)' nspin=',nspin
!$$$
!$$$c> max. no. valence and core l
!$$$      read(ifi,*)
!$$$      read(ifi,*)nl
!$$$c      write(6,*)' ispin nl =',ispin,nl
!$$$
!$$$c> max. no. valence and core n
!$$$      read(ifi,*)
!$$$      read(ifi,*)nnv,nnc
!$$$
!$$$ccccccccccccccccccccccccccc
!$$$      if(nnv==1) nnv=2 ! for backword compatibility!takao apr 2002
!$$$ccccccccccccccccccccccccccc
!$$$
!$$$c> max. no. radial mesh points
!$$$      read(ifi,*)
!$$$      read(ifi,*)nrx
!$$$
!$$$c> class-label, z, no. radial points
!$$$      read(ifi,*)
!$$$      allocate(z(nclass),nrofi(nclass))
!$$$      do      ic = 1,nclass
!$$$        read(ifi,*)clablxxx,z(ic),nrofi(ic)
!$$$      end do
!$$$
!$$$c> atom and its class
!$$$      allocate(iclass(natom))
!$$$      do n = 1,natom
!$$$        iclass(n)=n
!$$$      end do
!$$$
!$$$      allocate(nindxv(nl,nclass),nindxc(nl,nclass),
!$$$     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
!$$$     &        ooo(nl,nnv,nclass),
!$$$     &        occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
!$$$      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))
!$$$      allocate( cutbase(0:2*(nl-1)) )
!$$$ctakao
!$$$      call rgwinaf    (ifi,nl,nnv,nnc,nclass,
!$$$c> BZ
!$$$     o                  n1,n2,n3,ef,
!$$$c> frequencies
!$$$     o                  niw,diw,nw,dw,delta,deltaw,esmr,imagw,
!$$$c> coulomb
!$$$c     o                  tolvc,alp,alptx,h,ng,
!$$$c> product basis
!$$$     o                  cutbase,lcutmx,nindxv,nindxc,
!$$$     o                  occv,unoccv, occc,unoccc,
!$$$c> core
!$$$     o                  ncwf,ncwf2 )
!$$$c----
!$$$c      cutbase = tolbas
!$$$      inquire(file='NoCore',exist=nocore)
!$$$      if(nocore) then
!$$$        occc=0
!$$$        unoccc=0
!$$$        ncfw  =0
!$$$      elseif( incwfx==-1 ) then
!$$$        write(6,*)' ### incwf=-1 Use ForSxc for core'
!$$$        ncwf = ncwf2
!$$$      elseif( incwfx==-2 ) then
!$$$        write(6,*)' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
!$$$        call notbit(nl*nnc*nclass, ncwf2)
!$$$        ncwf  = ncwf2
!$$$        occc= ncwf
!$$$        unoccc= 0
!$$$cccccccccccccccccccccccccccccccccc
!$$$c 31May2006
!$$$        ooo=0
!$$$        call ibitand(nl*nnv*nclass, unoccv,occv, ooo)
!$$$        unoccv = ooo
!$$$cccccccccccccccccccccccccccccccccc
!$$$      elseif( incwfx==-3 ) then
!$$$        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
!$$$        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
!$$$        write(6,*)' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
!$$$      elseif( incwfx==-4 ) then
!$$$        write(6,*)' ### incwf=-4  occ=0 and unocc=0 for all core '
!$$$        occc=0  !call iclear(nl*nnc*nclass, w(ioccc))
!$$$        unoccc=0  !call iclear(nl*nnc*nclass, w(iunoccc))
!$$$        ncwf  =0  !call iclear(nl*nnc*nclass, w(incwf))
!$$$      elseif(incwfx==0) then
!$$$        write(6,*)' ### Use unocc occ ForX0 for core'
!$$$      else
!$$$        call rx( ' ### proper incwf is not given for genallcf2:rgwinf ')
!$$$      endif
!$$$      deallocate(ncwf2)
!$$$      end subroutine
!$$$      end module
!$$$
!$$$c--------------------------------------------------------------------

!$$$
!$$$
!$$$
!$$$
!$$$
!$$$
!$$$
