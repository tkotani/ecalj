subroutine suldau(nbas,sspec,ssite,nlibu,lmaxu,lldau)
  use m_struc_def  !Cgetarg
  !- Finds lda+U sites and counts number of blocks
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa idu
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
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
  ! ... Passed parameters
  integer :: nbas,nlibu,lmaxu,lldau(nbas),igetss,is,ib,l,lmxa,idu(4),i_copy_size
  type(s_spec)::sspec(*)
  type(s_site)::ssite(*)
  nlibu = 0
  lmaxu = 0
  do  ib = 1, nbas
     lldau(ib) = 0
     is = int(ssite(ib)%spec)


     lmxa=sspec(is)%lmxa
     i_copy_size=size(sspec(is)%idu)
     call icopy(i_copy_size,sspec(is)%idu,1,idu,1)

     do  l = 0, min(lmxa,3)
        if (idu(l+1) /= 0) then
           if (lldau(ib) == 0) lldau(ib) = nlibu+1
           nlibu = nlibu+1
           lmaxu = max(lmaxu,l)
        endif
     enddo
  enddo

  if (nlibu /= 0) then
     call info2(10,1,0, &
          ' suldau:  %i U block(s)  lmaxu = %i',nlibu,lmaxu)
  endif

end subroutine suldau

