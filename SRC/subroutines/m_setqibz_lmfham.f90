!>qibz handling for lmfham1,2. a quick remedy.
module m_setqibz_lmfham
   real(8),allocatable,protected:: qibz(:,:),qbzii(:,:,:),wiqibz(:)
   integer,allocatable,protected:: irotq(:),irotg(:),ndiff(:,:),iqbzrep(:),iqii(:,:), ngx(:),igx(:,:)
   logical,allocatable,protected:: igiqibz(:,:)
   integer:: nqibz
   public:: set_qibz
contains
   subroutine set_qibz(plat,qbz,nqbz,symops,ngrp)
     use m_ftox
     use m_lgunit,only: stdo
     integer:: nqbz,ngrp,i,ig,ibz,iqibz,iqbz,nx
     real(8):: plat(3,3),qbz(3,nqbz),symops(3,3,ngrp),qp(3),qx(3)
     real(8),allocatable::wtemp(:)
     real(8),external::tolq !eps=1d-8
!      do ig=1,ngrp
!         write(6,ftox) 'gggggggggggggg',ig, ftof(reshape(symops(:,:,ig),[9]))
!      enddo   
      allocate(qibz(3,nqbz),irotq(nqbz),irotg(nqbz),ndiff(3,nqbz),iqbzrep(nqbz),wiqibz(nqbz))
      iqibz=0
      wiqibz=0d0
      do iqbz=1,nqbz
         qp = qbz(:,iqbz)
         do ig=1,ngrp
            do i=1,iqibz
               qx= matmul(transpose(plat),  qp-matmul(symops(:,:,ig),qibz(:,i)))
               qx=qx-nint(qx) !qx-ndiff !translation of qx
               if(sum(abs(qx))<tolq()) then
                  irotq(iqbz)=i
                  irotg(iqbz)=ig
                  ndiff(:,iqbz) = nint(qx)
                  wiqibz(i)=wiqibz(i)+1d0
                  goto 88
               endif
            enddo
         enddo
         iqibz=iqibz+1
         qibz(:,iqibz) = qp
         irotq(iqbz)=iqibz
!         write(6,ftox) 'iiiiiiiiiiii setqibz',iqbz,irotq(iqbz)
         irotg(iqbz)=1 !identical symops ig=1
         ndiff(:,iqbz)=0
         iqbzrep(iqibz)= iqbz !representative
         wiqibz(iqibz) = 1d0
88       continue
      enddo
      nqibz=iqibz
      allocate(ngx(nqibz),igx(ngrp,nqibz))
      do i=1,nqibz !ig to keep qibz
         nx=0
         do ig=1,ngrp
            qx= matmul(transpose(plat), qibz(:,i)-matmul(symops(:,:,ig),qibz(:,i)))
            qx=qx-nint(qx) !qx-ndiff !translation of qx
            if(sum(abs(qx))<tolq()) then
               nx=nx+1
               igx(nx,i) =ig
            endif
         enddo
         ngx(i)=nx
      enddo
      allocate(wtemp,source=wiqibz(1:nqibz))
      deallocate(wiqibz)
      call move_alloc(from=wtemp,to=wiqibz)
      allocate(igiqibz(ngrp,nqibz),qbzii(3,ngrp,nqibz),iqii(ngrp,nqibz))
      igiqibz=.false.
      do concurrent(iqbz=1:nqbz)
         iqibz= irotq(iqbz)
         ig   = irotg(iqbz)
         igiqibz(ig,iqibz) =.true.
         iqii(ig,iqibz)=iqbz
         qbzii(:,ig,iqibz) = qbz(:,iqbz)
      enddo
      wiqibz=wiqibz/nqbz
!      write(stdo,*)'sum of wiqibz=',sum(wiqibz)
!      do iqibz=1,nqibz
!         write(stdo,ftox)' qibz:  ',iqibz,' ',ftof(qibz(:,iqibz))
!      enddo
!      if(master_mpi) write(stdo,ftox) 'set_qibz: nqbz nqibz ngrp wiqibz=',nqbz,nqibz,ngrp,sum(wiqibz)
!    forall( iqibz=1:nqibz) nigiq(iqibz) = count(igiqibz(:,iqibz))
   endsubroutine set_qibz
end module m_setqibz_lmfham
