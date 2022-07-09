subroutine phispinsym_ssite_set()
  use m_lmfinit,only: nbas,nsp,sspec=>v_sspec,n0,ispec !ssite=>v_ssite
  use m_MPItk,only:master_mpi
  use m_struc_def
  use m_lgunit,only:stdo
  use m_density,only: pnuall,pnzall
  implicit none
  integer:: ib,is,lmxa,l
  real(8):: pmean
  real(8),pointer:: pnu(:,:),pnz(:,:)
  if(master_mpi)write(6,*)' --phispinsym use spin-averaged potential for phi and phidot'
  do ib = 1,nbas
!     pnu = ssite(ib)%pnu
!     pnz = ssite(ib)%pz
     pnu=>pnuall(:,1:nsp,ib)
     pnz=>pnzall(:,1:nsp,ib)
     is   = ispec(ib) !ssite(ib)%spec
     lmxa = sspec(is)%lmxa
     do l=0,lmxa
        pmean = sum(pnu(l+1,1:nsp))/nsp
        if(master_mpi .AND. nsp==2) write(stdo,"('  ibas l=',2i3,' pnu=',2f10.5,' -->',f10.5)") &
             ib,l,pnu(l+1,1:nsp),pmean
        pnu(l+1,1:nsp) = pmean
        if(pnz(l+1,1)/=0 ) then
           pmean = sum(pnz(l+1,1:nsp))/nsp
           if(master_mpi .AND. nsp==2)write(stdo,"('  ibas l=',2i3,' pnz=',2f10.5,' -->',f10.5)") &
                ib,l,pnz(l+1,1:nsp),pmean
           pnz(l+1,1:nsp) = pmean
        endif
     enddo
!     ssite(ib)%pnu=pnu
!     ssite(ib)%pz=pnz
     pnuall(:,1:nsp,ib)=pnu(:,1:nsp)
     pnzall(:,1:nsp,ib)=pnz(:,1:nsp)
  enddo
end subroutine phispinsym_ssite_set
