!ssssssssssssssssssssssssssssssssssssssssssss
      subroutine rgwinaf (ifi,nl,nnv,nnc,nclass,
c> BZ
     o                  n1,n2,n3,ef,
c> frequencies
     o                  niw,diw,nw,dw,delta,deltaw,esmr, imagw,
c> product basis
     o                  cutbase,lcutmx,nindxv,nindxc,
     o                  noccv,nunoccv,noccc,nunoccc,
c> core
     o                  ncwf,ncwf2)
      use m_keyvalue,only: getkeyvalue
c read the rest of GW input data
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension 
     o          nindxv(0:nl-1,nclass),
     o          nindxc(0:nl-1,nclass),
     o          noccv(0:nl-1,nnv,nclass),
     o          noccc(0:nl-1,nnc,nclass),
     o          nunoccv(0:nl-1,nnv,nclass),
     o          nunoccc(0:nl-1,nnc,nclass),
     o          ncwf(0:nl-1,nnc,nclass),
     o          ncwf2(0:nl-1,nnc,nclass),
     o  cutbase(0:2*(nl-1))
c      logical :: readgwinput
      integer(4):: ret ,ncc
c      real(8)::cutbase0
      logical::readon
c      write(6,*)' nclass nl-1=',nclass,nl-1
c> n index for valence and core
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          read(ifi,*) ict,lt,nindxv(l,ic),nindxc(l,ic)
          if(lt  .ne. l ) call rx( 'rgwina: 1st wrong l ')
        end do
      end do
      write(6,*)' --- valence product basis section'
c> criteria for product basis
c>> valence
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          do       n = 1,nindxv(l,ic)
            read(ifi,*)ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
            write(6,"(100i3)") ict,lt,nt,noccv(l,n,ic),nunoccv(l,n,ic)
            if(lt  .ne. l )call rx( 'rgwina: 2nd wrong l valence')
            if(nt  .ne. n )call rx( 'rgwina: wrong n valence')
          end do
        end do
      end do
 1099 continue
c      write(6,*)' goto prod 3'
c>> core
      write(6,*)' --- core product basis section'
      read(ifi,6000)blank
      do      ic = 1,nclass
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            read(ifi,*)ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic)
     & ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
            write(6,"(100i3)") 
     & ict,lt,nt,noccc(l,n,ic),nunoccc(l,n,ic),ncwf(l,n,ic)
     & ,ncwf2(l,n,ic) !ncwf2 is for Sigma calcuation
            if(lt  .ne. l )call rx( 'rgwina: 2nd wrong l core')
            if(nt  .ne. n )call rx( 'rgwina: wrong n core')
          end do
        end do
      end do
 6000 format(a)
      close(ifi)
      return
      end
      
c$$$      module m_rgwinf_v3
c$$$      integer,protected:: nclass,natom,nspin,nl,nnv,nnc,lcutmx,nrx
c$$$      real(8),protected:: alat
c$$$      integer,allocatable,protected:: 
c$$$     &   iclass(:)
c$$$     &  ,nindxv(:,:),nindxc(:,:)
c$$$     &  ,occv(:,:,:),unoccv(:,:,:),ooo(:,:,:)
c$$$     &  ,occc(:,:,:),unoccc(:,:,:)
c$$$     &  ,ncwf(:,:,:)
c$$$      real(8),allocatable,protected:: z(:),cutbase(:)
c$$$
c$$$      contains
c$$$
c$$$      subroutine rgwinf_v3 (iflmto,incwfx) 
c$$$      implicit real*8(a-h,o-z)
c$$$      implicit integer(i-n)
c$$$      intent(in)::          iflmto,incwfx
c$$$C- readin GWIN_V2 and LMTO(crystal) data.
c$$$C all the output are given in the declear section
c$$$C!
c$$$Cr Return iclass=ibas.
c$$$Cr nwin,efin,incwfx, are used as switches.
c$$$C--------------------------------------------------------
c$$$c      character(120):: symgrp
c$$$c      character(120):: symgrpt
c$$$      integer(4):: infwfx,nwin,iflmto,ifiniin,incwfx
c$$$      logical :: nocore
c$$$      integer(4),allocatable::ncwf2(:,:,:),nrofi(:)
c$$$      real(8)::efin
c$$$      character(6)::clablxxx
c$$$      ifi = iflmto
c$$$      ef  = -999d0 ! not readin efermi
c$$$      nw  = -999 !Not readin NW file
c$$$
c$$$c---------
c$$$c SYMMETRY
c$$$c---------
c$$$c      write(6,*)' goto sym'
c$$$c      read(ifi,*);  read(ifi,*)
c$$$c      read(ifi,*)symgrpt
c$$$c      j           = 0
c$$$c      call rmvbl   (symgrpt,120,j)
c$$$c      symgrp(1:2) = '  '
c$$$c      symgrp(3:120) = symgrpt(j+1:120)
c$$$      read(ifi,*)
c$$$      read(ifi,*)
c$$$      read(ifi,*)
c$$$
c$$$c----------
c$$$c STRUCTURE
c$$$c----------
c$$$c      write(6,*)' goto structure'
c$$$c> lattice constant
c$$$      read(ifi,*);  read(ifi,*); read(ifi,*)
c$$$      read(ifi,*)alat
c$$$
c$$$c> primitive lattice vectors
c$$$c      allocate(plat(3,3))
c$$$      read(ifi,*)
c$$$      read(ifi,*)!plat(1:3,1)
c$$$      read(ifi,*)!plat(1:3,2)
c$$$      read(ifi,*)!plat(1:3,3)
c$$$
c$$$c> no. atoms
c$$$      read(ifi,*)
c$$$      read(ifi,*)natom
c$$$
c$$$c We assume nclass=natom
c$$$      nclass = natom
c$$$
c$$$c> positions of atoms
c$$$c      allocate(pos(3,natom))
c$$$      read(ifi,*)
c$$$      do n = 1,natom
c$$$        read(ifi,*) !pos(1,n),pos(2,n),pos(3,n)
c$$$      end do
c$$$
c$$$c-----
c$$$c LMTO
c$$$c-----
c$$$      write(6,*)' goto lmto'
c$$$c> spin (1=paramagnetic  2=ferromagnetic)
c$$$      read(ifi,*)
c$$$      read(ifi,*)
c$$$      read(ifi,*)
c$$$      read(ifi,*)nspin
c$$$c      write(6,*)' nspin=',nspin
c$$$
c$$$c> max. no. valence and core l
c$$$      read(ifi,*)
c$$$      read(ifi,*)nl
c$$$c      write(6,*)' ispin nl =',ispin,nl
c$$$
c$$$c> max. no. valence and core n
c$$$      read(ifi,*)
c$$$      read(ifi,*)nnv,nnc
c$$$
c$$$ccccccccccccccccccccccccccc
c$$$      if(nnv==1) nnv=2 ! for backword compatibility!takao apr 2002
c$$$ccccccccccccccccccccccccccc
c$$$
c$$$c> max. no. radial mesh points
c$$$      read(ifi,*)
c$$$      read(ifi,*)nrx
c$$$
c$$$c> class-label, z, no. radial points
c$$$      read(ifi,*)
c$$$      allocate(z(nclass),nrofi(nclass))
c$$$      do      ic = 1,nclass
c$$$        read(ifi,*)clablxxx,z(ic),nrofi(ic)
c$$$      end do
c$$$
c$$$c> atom and its class
c$$$      allocate(iclass(natom))
c$$$      do n = 1,natom
c$$$        iclass(n)=n
c$$$      end do
c$$$
c$$$      allocate(nindxv(nl,nclass),nindxc(nl,nclass),
c$$$     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
c$$$     &        ooo(nl,nnv,nclass),
c$$$     &        occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
c$$$      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))
c$$$      allocate( cutbase(0:2*(nl-1)) )
c$$$ctakao
c$$$      call rgwinaf    (ifi,nl,nnv,nnc,nclass,
c$$$c> BZ
c$$$     o                  n1,n2,n3,ef,
c$$$c> frequencies
c$$$     o                  niw,diw,nw,dw,delta,deltaw,esmr,imagw,
c$$$c> coulomb
c$$$c     o                  tolvc,alp,alptx,h,ng,
c$$$c> product basis
c$$$     o                  cutbase,lcutmx,nindxv,nindxc,
c$$$     o                  occv,unoccv, occc,unoccc,
c$$$c> core
c$$$     o                  ncwf,ncwf2 )
c$$$c----
c$$$c      cutbase = tolbas
c$$$      inquire(file='NoCore',exist=nocore)
c$$$      if(nocore) then
c$$$        occc=0    
c$$$        unoccc=0  
c$$$        ncfw  =0  
c$$$      elseif( incwfx==-1 ) then
c$$$        write(6,*)' ### incwf=-1 Use ForSxc for core'
c$$$        ncwf = ncwf2  
c$$$      elseif( incwfx==-2 ) then
c$$$        write(6,*)' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
c$$$        call notbit(nl*nnc*nclass, ncwf2)
c$$$        ncwf  = ncwf2 
c$$$        occc= ncwf  
c$$$        unoccc= 0     
c$$$cccccccccccccccccccccccccccccccccc
c$$$c 31May2006
c$$$        ooo=0
c$$$        call ibitand(nl*nnv*nclass, unoccv,occv, ooo)
c$$$        unoccv = ooo
c$$$cccccccccccccccccccccccccccccccccc
c$$$      elseif( incwfx==-3 ) then
c$$$        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
c$$$        unoccc= 0     ! call iclear(nl*nnc*nclass, w(iunoccc))
c$$$        write(6,*)' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
c$$$      elseif( incwfx==-4 ) then
c$$$        write(6,*)' ### incwf=-4  occ=0 and unocc=0 for all core '
c$$$        occc=0  !call iclear(nl*nnc*nclass, w(ioccc))
c$$$        unoccc=0  !call iclear(nl*nnc*nclass, w(iunoccc))
c$$$        ncwf  =0  !call iclear(nl*nnc*nclass, w(incwf))
c$$$      elseif(incwfx==0) then
c$$$        write(6,*)' ### Use unocc occ ForX0 for core'
c$$$      else
c$$$        call rx( ' ### proper incwf is not given for genallcf2:rgwinf ')
c$$$      endif
c$$$      deallocate(ncwf2)
c$$$      end subroutine
c$$$      end module
c$$$
c$$$c--------------------------------------------------------------------
c
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
