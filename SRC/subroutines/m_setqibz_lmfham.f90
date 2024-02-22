!>qibz handling for lmfham1,2. a quick remedy.
module m_setqibz_lmfham
   real(8),allocatable,protected:: qibz(:,:),qbzii(:,:,:)
   integer,allocatable,protected:: irotq(:),irotg(:),ndiff(:,:),iqbzrep(:),iqii(:,:)
   logical,allocatable,protected:: igiqibz(:,:)
   integer:: nqibz
   public:: set_qibz
contains
   subroutine set_qibz(plat,qbz,nqbz,symops,ngrp)
      integer:: nqbz,ngrp,i,ig,ibz,iqibz,iqbz
      real(8)::eps=1d-8
      real(8):: plat(3,3),qbz(3,nqbz),symops(3,3,ngrp),qp(3),qx(3)
      allocate(qibz(3,nqbz),irotq(nqbz),irotg(nqbz),ndiff(3,nqbz),iqbzrep(nqbz))
      iqibz=0
      do iqbz=1,nqbz
         qp = qbz(:,iqbz)
         do ig=1,ngrp
            do i=1,iqibz
               qx= matmul(transpose(plat),  qp-matmul(symops(:,:,ig),qibz(:,i)))
               qx=qx-nint(qx) !qx-ndiff !translation of qx
               if(sum(abs(qx))<eps) then
                  irotq(iqbz)=i
                  irotg(iqbz)=ig
                  ndiff(:,iqbz) = nint(qx)
                  goto 88
               endif
            enddo
         enddo
         iqibz=iqibz+1
         qibz(:,iqibz) = qp
         irotq(iqbz)=iqibz
         irotg(iqbz)=1 !identical symops ig=1
         ndiff(:,iqbz)=0
         iqbzrep(iqibz)= iqbz !representative
88       continue
      enddo
      nqibz=iqibz
      allocate(igiqibz(ngrp,nqibz),qbzii(3,ngrp,nqibz),iqii(ngrp,nqibz))
      igiqibz=.false.
      do concurrent(iqbz=1:nqbz)
         iqibz= irotq(iqbz)
         ig   = irotg(iqbz)
         igiqibz(ig,iqibz) =.true.
         iqii(ig,iqibz)=iqbz
         qbzii(:,ig,iqibz) = qbz(:,iqbz)
      enddo
!    forall( iqibz=1:nqibz) nigiq(iqibz) = count(igiqibz(:,iqibz))
      !write(6,*) 'nqbz nqibz ngrp=',nqbz,nqibz,ngrp
   endsubroutine set_qibz
end module m_setqibz_lmfham
