subroutine dfqkkl( oqkkl ) !Allocates arrays to accumulate output site density
  use m_lmfinit,only: nkaph,nsp,nbas,ispec,sspec=>v_sspec
  use m_struc_def,only:s_rv5
  !o oqkkl : memory is allocated for qkkl
  !  nkaph :number of orbital types for a given L quantum no. in basis
  implicit none
  type(s_rv5) :: oqkkl(3,nbas)
  integer :: ib,is,kmax,lmxa,lmxh,nlma,nlmh
  do  ib = 1, nbas
     is = ispec(ib) 
     lmxa=sspec(is)%lmxa
     if (lmxa == -1) cycle
     nlma = (lmxa+1)**2
     nlmh = (sspec(is)%lmxb+1)**2
     kmax =  sspec(is)%kmxt
     if(allocated(oqkkl(1,ib)%v)) deallocate(oqkkl(1,ib)%v,oqkkl(2,ib)%v,oqkkl(3,ib)%v)
     allocate(oqkkl(1,ib)%v(0:kmax,0:kmax,nlma,nlma,nsp), source=0d0)  ! Pkl*Pkl
     allocate(oqkkl(2,ib)%v(nkaph,0:kmax,nlmh,nlma   ,nsp), source=0d0)! Pkl*Hsm
     allocate(oqkkl(3,ib)%v(nkaph,nkaph,nlmh,nlmh    ,nsp), source=0d0)! Hsm*Hsm
  enddo
end subroutine dfqkkl
