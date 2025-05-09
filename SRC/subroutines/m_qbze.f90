!>Get qbze qibze, whic are extented qbz qibz including q0p shifts
module m_qbze 
  use NaNum,only: NaN
  use m_read_bzdata,only: nqbz,nqibz, qbz,qibz, q0i,nq0i,nq0iadd
  use m_mpi,only:ipr
  implicit none
  integer,protected::    nqbze=NaN, nqibze=NaN
  real(8),allocatable,protected :: qbze(:,:),qibze(:,:)
contains
  subroutine setqbze()
    integer:: i,iq
    nqbze  = nqbz *(1 + nq0i+nq0iadd)
    nqibze = nqibz + nq0i+nq0iadd
    allocate( qbze(3, nqbze), qibze(3, nqibze))
    qibze(:,1:nqibz)= qibz(:,1:nqibz)
    qibze(:,nqibz+1:nqibz+nq0i+nq0iadd ) = q0i(:, 1:nq0i+nq0iadd)
    qbze (:,1:nqbz) = qbz (:,1:nqbz)
    do i = 1,nq0i+nq0iadd
       do iq = 1,nqbz
          qbze (:,nqbz*i + iq) = qbz(:,iq) + q0i(:,i)
       enddo
    enddo
    if(ipr) write(6,*) ' m_qbze: nqibz nqibze=',nqibz,nqibze
  end subroutine setqbze
end module m_qbze
