      subroutine genallcf_v3_dump(
c> structure
     o                   plat,alat,natom,nclass, pos,
c> symmetry
     o                   symgrp,gen,ngnmax,ngrp, symgg,  !igrp,
     o                   invg,ef,
     o                   clabl,nspin,nl,nn,nnv,nnc,
     o                   nindx,nindxv,nindxc,iclass,
     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
     o                   z, 
     o                   il, in, im, ilnm, nlnm,
     o                   ilv,inv,imv,ilnmv,nlnmv,
     o                   ilc,inc,imc,ilnmc,nlnmc,
     o                   ncwf,ecore,konf,icore,ncore,nctot,
     o                   imagw,niw,diw,nw,dw,delta,deltaw,esmr,freq )
c----------------------------------------------------------------------------
      implicit none
      integer(4):: nclass,natom,nspin,nl,nnv,nnc,nn
      character(6) :: clabl(nclass)
      integer(4)::ngrp,
     o        nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
     o        nctot,
     o        niw,nw
      real(8) :: alat,ef
      real(8) :: diw,dw,delta,deltaw,esmr
      integer(4):: 
     &   iclass(*)
     &  ,nindxv(*),nindxc(*)
     &  ,ncwf(*)
      integer(4)::
     o    invg(*),
     o    il(*), in(*), im(*),   ilnm(*),  nlnm(*),
     o           ilv(*),inv(*),imv(*),  ilnmv(*), nlnmv(*),
     o           ilc(*),inc(*),imc(*),  ilnmc(*), nlnmc(*),
     o   nindx(*),konf(*),icore(*),ncore(*)
      real(8) :: plat(9),pos(*),z(*)
      real(8) :: symgg(*), ecore(*), freq(*)
      character*120 symgrp
c      parameter (mxclass=100)
c      character*6 clabl(mxclass)
c      integer(4),parameter:: ngnmax=10
      integer(4):: ngnmax,imagw,ifi,ifi2
      real(8) :: gen(9,ngnmax)
c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=10)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)

c--------------------------------------------------------
      ifi =3012
      ifi2=3023
      write(6,*)'genallcf_v3_dump:'
      open (ifi, form='unformatted',file="genallcf3.dump")
      open (ifi2,file="genallcft3.dump")
      write(6,*)'genallcf_v3_dump:2'
      write(ifi)
     o                   plat,alat,natom,nclass,
     o                   symgrp,gen,ngnmax,ngrp,  ef,
     &     clabl, nl,nn,nnv,nnc,
     d     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
     &     nctot,
     &     niw,diw,nw,dw,delta,deltaw,esmr
      write(6,*)'genallcf_v3_dump:3'

      write(ifi2,*)
     o                   plat,alat,natom,nclass,
     o                   symgrp,gen,ngnmax,ngrp,  ef,
     &     clabl, nl,nn,nnv,nnc,
     d     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
     &     nctot,
     &     niw,diw,nw,dw,delta,deltaw,esmr
      write(6,*)'genallcf_v3_dump:4 nw=',nw
c      close(ifi)
c      close(ifi2)

      call ddump(ifi,ifi2,pos,3*natom)
      write(6,*)'genallcf_v3_dump:5'
      call ddump(ifi,ifi2,symgg,9*ngrp)

      call idump(ifi,ifi2,invg,ngrp)
      call idump(ifi,ifi2,nindx,nl*nclass)
      call idump(ifi,ifi2,nindxv,nl*nclass)
      call idump(ifi,ifi2,nindxc,nl*nclass)
      call idump(ifi,ifi2,iclass,natom)

      call ddump(ifi,ifi2,z,nclass)

      call idump(ifi,ifi2,il,nlnmx*nclass)
      call idump(ifi,ifi2,in,nlnmx*nclass)
      call idump(ifi,ifi2,im,nlnmx*nclass)
      call idump(ifi,ifi2,ilnm,nn*nl*nl*nclass)
      call idump(ifi,ifi2,ilv,nlnmxv*nclass)
      call idump(ifi,ifi2,inv,nlnmxv*nclass)
      call idump(ifi,ifi2,imv,nlnmxv*nclass)
      call idump(ifi,ifi2,ilnmv,nnv*nl*nl*nclass)
      call idump(ifi,ifi2,ilc,nlnmxc*nclass)
      call idump(ifi,ifi2,inc,nlnmxc*nclass)
      call idump(ifi,ifi2,imc,nlnmxc*nclass)
      call idump(ifi,ifi2,ilnmc,nnc*nl*nl*nclass)
      call idump(ifi,ifi2,nlnmv,nclass)
      call idump(ifi,ifi2,nlnmc,nclass)
      call idump(ifi,ifi2,nlnm,nclass)

c     o                   incwf,iecore,ikonf,iicore,incore,
      call idump(ifi,ifi2,ncwf,nl*nnc*nclass)
      call ddump(ifi,ifi2,ecore,2*nctot)
      call idump(ifi,ifi2,konf,nl*nclass)
      write(ifi2,*)"-------------------------------"
      call idump(ifi,ifi2,icore,nl*nl*nnc*nclass)
      call idump(ifi,ifi2,ncore,nclass)
      if(nw>0) call ddump(ifi,ifi2,freq,nw)
      end

      subroutine ddump(ifi,ifi2,a,n)
      real(8)::  a(*)
      write(ifi) a(1:n)
      write(ifi2,*) a(1:n)
      end
      subroutine idump(ifi,ifi2,a,n)
      integer(4)::a(*)
      write(ifi) a(1:n)
      write(ifi2,*) a(1:n)
      end
