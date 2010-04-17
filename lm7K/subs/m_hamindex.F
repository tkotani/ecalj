      module m_hamindex
c index for Hamiltonian
      implicit none
      integer,parameter:: null=-999999
      integer:: ngrp=null, lxx=null, kxx=null, norbt=null, napwmx=null, nkt=null,ndimham=null
      integer,allocatable:: ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),invgx(:),miat(:,:),offlrev(:,:,:),ibastab(:)
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:)
      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),ibasindex(:)
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

      end module


