!> This is in lm7K/subs/m_hamindex.F and in fpgw/gwsrc/m_hamindex.F 
!! We will need to unify make system and source code in fpgw and lmf.
!! norbtx is given in gwsrc/readeigen.F init_readeigen2
      module m_hamindex
      implicit none
      integer,parameter:: null=-999999
      integer:: ngrp=null, lxx=null, kxx=null,norbmto=null,norbtx=null
      integer:: imx=null,nbas,ndimham=null
      integer:: nqtt,nqi, nqnum,ngpmx
      integer,allocatable:: ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),offlrev(:,:,:),ibastab(:)
      integer,allocatable:: iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:) !,ngvecp(:,:,:),ngvecprev(:,:,:,:)
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:)
      real(8):: plat(3,3),qlat(3,3)
      real(8),allocatable:: qtt(:,:),qtti(:,:)
      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
      integer::  napwmx=null,lxxa=null
      logical,private:: debug=.false.
      contains

      integer function getikt(qin)
      integer::i
      real(8):: qin(3)
c      if(debug) print *,'nkt=',nkt
      do i=1, nqnum !*2 !nkt
        if(debug) print *,i,qin, qq(:,i)
        if(sum (abs(qin-qq(:,i)))<1d-8) then
          getikt=i
          return
        endif
      enddo
      print *,' getikt: xxx error nqnum qin=',nqnum,qin
      do i=1, nqnum !*2 !nkt
        write(*,"('i qq=',i3,3f11.5)")i, qq(:,i)
      enddo  
Cstop2rx 2013.08.09 kino      stop ' getikt can not find ikt for given q'
      call rx( ' getikt can not find ikt for given q')
      end function

      subroutine writehamindex()
      integer(4):: ifi
      logical::pmton
!! == write info for wave rotation. feb2012takao ==
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      write(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      write(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      write(ifi)lxxa
      write(ifi)dlmm
      write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
      write(ifi)qq !,ngvecp,ngvecprev
      write(ifi)plat,qlat,napwmx
      if(napwmx/=0) then !for APW rotation used in rotwvigg
         write(ifi) igv2,napwk,igv2rev
      endif   
      close(ifi)
      end subroutine writehamindex

      subroutine readhamindex()
!! == read info for wave rotation. feb2012takao ==
      integer(4):: ifi,nkt
      logical::pmton
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      allocate(symops(3,3,ngrp),ag(3,ngrp),qtt(3,nqtt),qtti(3,nqi))
      allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      read(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      allocate( offlrev(nbas,0:lxx,kxx))
      read(ifi) lxxa
      allocate( dlmm(-lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      read(ifi) dlmm
      read(ifi)ibastab,ltab,ktab,offl,offlrev
c      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
c      allocate( ngvecp(3,ngpmx,nqnum) )
      allocate( qq(3,nqnum)) !this was qq(3,nqnum*2) until Aug2012 when shorbz had been used.
      read(ifi)qq !,ngvecp,ngvecprev
      read(ifi)plat,qlat,napwmx
      if(napwmx/=0)then !for APW rotation used in rotwvigg
        nkt=nqtt
        allocate( igv2(3,napwmx,nkt) )
        allocate( napwk(nkt))
        allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
        read(ifi) igv2,napwk,igv2rev
      endif  
      close(ifi)
      end subroutine readhamindex
      end module


