subroutine suldau(nbas,sspec,nlibu,lmaxu,lldau)
  use m_struc_def  
  use m_ftox
  use m_lgunit,only:stdo
  use m_lmfinit,only:idu,ispec
  !- Finds lda+U sites and counts number of blocks
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !o Outputs
  !i   lldau :lldau(ib)=0 => no U on this site otherwise
  !i         :U on site ib with dmat in dmats(*,lldau(ib))
  !o   nlibu :number of LDA+U blocks
  !o   lmaxu :highest l for which a U block found, used as
  !o         :dimensioning parameter for U matrix
  !r Remarks
  !r
  !u Updates
  !u   27 Apr 05 (Lambrecht) first created
  !------------------------------------
  implicit none
  integer :: nbas,nlibu,lmaxu,lldau(nbas),igetss,is,ib,l,lmxa !,idu(4),i_copy_size
  type(s_spec)::sspec(*)
  nlibu = 0
  lmaxu = 0
  do  ib = 1, nbas
     lldau(ib) = 0
     is  = ispec(ib) !ssite(ib)%spec
     lmxa= sspec(is)%lmxa
     !idu = sspec(is)%idu
     do  l = 0, min(lmxa,3)
        if (idu(l+1,is) /= 0) then
           if (lldau(ib) == 0) lldau(ib) = nlibu+1
           nlibu = nlibu+1
           lmaxu = max(lmaxu,l)
        endif
     enddo
  enddo
  if(nlibu/=0) write(stdo,ftox)'suldau:  ',nlibu,' U block(s)  lmaxu =',lmaxu
end subroutine suldau


