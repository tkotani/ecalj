subroutine dfqkkl( sv_p_oqkkl )
  use m_lmfinit,only: nkaph,nspc,nsp,nbas,ispec,sspec=>v_sspec
  use m_struc_def  !Cgetarg
  !- Allocates arrays to accumulate local output density
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa lmxb kmxt
  !i     Stored:    *
  !i     Passed to: *
  !o Outputs
  !o   oqkkl :memory is allocated for qkkl
  !l Local variables
  !l   nkapi :number of envelope function types per l q.n. for spec is2
  !l   nkaph :number of orbital types for a given L quantum no. in basis
  !r Remarks
  !u Updates
  !u   01 Jul 05 handle lmxa=-1 -> no allocation
  !u   15 Jun 05 Allocation for noncollinear case
  !u   25 Aug 01 Extended to local orbitals
  !u   15 Jun 00 spin polarized
  !u   22 Apr 00 Adapted from nfp df_qkkl.f
  ! ----------------------------------------------------------------------
  implicit none
  type(s_rv1) :: sv_p_oqkkl(3,nbas)
  integer :: ib,igetss,is,kmax,lmxa,lmxh,nelt1,nlma,nlmh,nelt3,nelt2
  do  ib = 1, nbas
     is = ispec(ib) 
     lmxa=sspec(is)%lmxa
     lmxh=sspec(is)%lmxb
     kmax=sspec(is)%kmxt
     if (lmxa == -1) goto 10
     nlma = (lmxa+1)**2
     nlmh = (lmxh+1)**2
     !   ... Case Pkl*Pkl
     nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
     if(allocated(sv_p_oqkkl(1,ib)%v)) then
        deallocate(sv_p_oqkkl(1,ib)%v,sv_p_oqkkl(2,ib)%v,sv_p_oqkkl(3,ib)%v)
     endif
     allocate(sv_p_oqkkl(1,ib)%v(nelt1*nsp*nspc))
     sv_p_oqkkl(1,ib)%v=0d0
     !   ... Case Pkl*Hsm
     nelt2 = (kmax+1)*nkaph*nlma*nlmh
     allocate(sv_p_oqkkl(2,ib)%v(nelt2*nsp*nspc))
     sv_p_oqkkl(2,ib)%v=0d0
     !   ... Case Hsm*Hsm
     nelt3 = nkaph*nkaph*nlmh*nlmh
     allocate(sv_p_oqkkl(3,ib)%v(nelt3*nsp*nspc))
     sv_p_oqkkl(3,ib)%v=0d0
10   continue
  enddo
end subroutine dfqkkl


