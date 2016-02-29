!> This is in lm7K/subs/m_hamindex.F and in fpgw/gwsrc/m_hamindex.F 
!! We will need to unify make system and source code in fpgw and lmf.
!! norbtx is given in gwsrc/readeigen.F init_readeigen2
      module m_hamindex
      integer,parameter,private:: null=-999999
      logical,private:: debug=.false.

      integer,protected:: ngrp=null, lxx=null, kxx=null,norbmto=null
      integer,protected:: nbas,nqtt,ndimham=null
      integer,allocatable,protected:: ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),offlrev(:,:,:),ibastab(:)
      integer,allocatable,protected:: iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:) !,ngvecp(:,:,:),ngvecprev(:,:,:,:)
      real(8),allocatable,protected:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:)
      real(8),protected:: plat(3,3),qlat(3,3)
      real(8),allocatable,protected:: qtt(:,:),qtti(:,:)
      integer,allocatable,protected:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
      integer,protected::  napwmx=null,lxxa=null


      integer:: norbtx=null
      integer:: nqi, nqnum,ngpmx,imx=null
      contains

!> get index ikt such that for qin(:)=qq(:,ikt)
      integer function getikt(qin) !return
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
      call rx( ' getikt can not find ikt for given q')
      end function

!> write info for wave rotation. 
      subroutine writehamindex()
      integer(4):: ifi
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('writehamindex is already done')
      done=.true.
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

!> read info for wave rotation. 
      subroutine readhamindex()
      integer(4):: ifi,nkt
      logical::pmton
      logical,save:: done=.false.
      if(done) call rx('readhamindex is already done')
      done=.true.
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      allocate(symops(3,3,ngrp),ag(3,ngrp),qtt(3,nqtt),qtti(3,nqi))
      allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      write(6,*) 'ngrp=',ngrp
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
      done=.true.
      end subroutine readhamindex
      end module


