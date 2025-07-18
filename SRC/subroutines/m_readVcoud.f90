!> read Coulomb matrix for given q
module m_readVcoud        
  use m_readqg,only: Readqg0
  use m_rdpp,only: nbloch
  implicit none
  integer,protected :: ngb, ngc, iq
  real(8),protected :: q(3)
  real(8),allocatable,protected    :: vcousq(:), vcoud(:) 
  complex(8),allocatable,protected :: zcousq(:,:) 
contains
  subroutine Readvcoud(q,iq,NoVcou) !readin diagonalized Coulomb interaction
    intent(in)::       q,iq,NoVcou
    !! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
    !! ===Readin diagonalized Coulomb interaction===
    !!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
    !!   A possible choice for paralellization is "Vcoud.ID" files where ID=kx
    !!   Vould file is written in hvccfp0.m.F.
    !! NOTE: vcoud and zcousq are used in m _zmelt.
    integer:: iq
    logical:: NoVcou
    real(8):: q(3),quu(3)
    integer:: ngb0,ifvcoud
    real(8):: qvv(3)
    character(128):: vcoudfile
    character(10) :: i2char
    call Readqg0('QGcou',q,   quu,ngc) ! ngc: the number of IPW for the interaction matrix (in QGcou),
    ! quu: equivalent to q generated in qg4g4.
    ngb = ngc+nbloch
    if(NoVcou) return
    open(newunit=ifvcoud, file=trim('__Vcoud.'//trim(i2char(iq))),form='unformatted')
    read(ifvcoud) ngb0
    read(ifvcoud) qvv
    if(allocated(zcousq)) deallocate(zcousq)
    if(allocated(vcousq)) deallocate(vcousq)
    if(allocated(vcoud)) deallocate(vcoud)
    allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0))
    read(ifvcoud) vcoud
    read(ifvcoud) zcousq
    close(ifvcoud)
    vcousq = sqrt(vcoud)
!    write(6,*)'voud sumcheck1111=',sum(abs(vcousq)),sum(abs(zcousq))
!    write(6,*)'voud sumcheck1112=',sum(abs(qvv)),ngb
  end subroutine Readvcoud
  subroutine ReleaseZcousq
    if(allocated(zcousq)) deallocate(zcousq)
  end subroutine ReleaseZcousq
end module m_readVcoud
