      module m_hamindex
c index for Hamiltonian
      implicit none
      integer,parameter:: null=-999999
      integer:: ngrp=null, lxx=null, kxx=null, norbt=null, napwmx=null, nkt=null,ndimham=null
      integer,allocatable:: ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),invgx(:),miat(:,:),offlrev(:,:,:),ibastab(:)
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:)
      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),ibasindex(:)
      real(8):: plat(3,3),qlat(3,3)
      integer:: imx=null,nbas
      contains

      integer function getikt(qin)
      integer::i,debugmode
      real(8):: qin(3)
      if(debugmode()>0) print *,'nkt=',nkt
      do i=1, nkt
        if(debugmode()>0) print *,i,qin, qq(:,i)
        if(sum (abs(qin-qq(:,i)))<1d-8) then
          getikt=i
          return
        endif
      enddo
      print *,' getikt: error qin=',qin
      call rx(' getikt can not find ikt for given q')
      end function

      subroutine writehamindex()
      integer(4):: ifi
c ... write info for rotation. feb2012takao
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      write(ifi)ngrp,nbas,kxx,lxx
      write(ifi)symops,ag,invgx,miat,tiat,shtvg
      write(ifi)ibastab,ltab,ktab,offl,offlrev,dlmm
      write(ifi)imx,napwmx,nkt
      write(ifi)qq,igv2,igv2rev
      write(ifi)plat,qlat
      close(ifi)
      end subroutine

      subroutine readhamindex()
      integer(4):: ifi
c ... read info for rotation. feb2012takao
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      read(ifi)ngrp,nbas,kxx,lxx
      allocate(symops(3,3,ngrp),ag(3,ngrp))
      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      read(ifi)symops,ag,invgx,miat,tiat,shtvg
      allocate( ltab(norbt),ktab(norbt),offl(norbt),ibastab(norbt) )
      allocate( offlrev(nbas,0:lxx,kxx))
      allocate( dlmm( -lxx:lxx, -lxx:lxx, 0:lxx, ngrp))
      read(ifi)ibastab,ltab,ktab,offl,offlrev,dlmm
      read(ifi)imx,napwmx,nkt
      allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
      allocate( igv2(3,napwmx,nkt) )
      allocate( qq(3,nkt) )
      read(ifi)qq,igv2,igv2rev
      read(ifi)plat,qlat
      close(ifi)
      end subroutine

      end module


