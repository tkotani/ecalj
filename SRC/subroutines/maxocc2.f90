integer function maxocc2 (nspin,ef, nband, qbz,nqbz) ! maximum no. occupied states
  use m_readeigen, only: readeval
  implicit none
  integer :: nspin,nqbz,nband,is,iq
  real(8) :: qbz(3,nqbz),ef,ekt(nband, nqbz,nspin )
  do is = 1,nspin
     do iq = 1,nqbz
        ekt(:,iq,is)= readeval(qbz(:,iq),is)
     enddo
  enddo
  maxocc2 = maxval(count(ekt(1:nband,1:nqbz,1:nspin)<ef,1))
end function maxocc2

