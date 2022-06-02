module m_pomat
  implicit none

  complex(8),allocatable,protected:: pomat(:,:)
  integer,protected:: nn,no

  !-----------------------------------------------------------------
contains
  subroutine readpomat(q)
    intent(in)::         q
    ! This returns nn,no, and pomat(nn,no)
    real(8):: q_r(3),q(3)
    integer(4):: nn,iqx,isx
    integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,iopen,iclose
    !      complex(8),allocatable:: pomat(:,:)
    !      if( allocated(pomat) ) deallocate(pomat)
    open(newunit=ifpomat,file='POmat',form='unformatted')
    !... smoothed mixed basis !oct2005
    ! This replace original zmelt with new zmelt based on smoothed mixed basis.
    do
       read(ifpomat) q_r,nn,no,iqx !readin reduction matrix pomat
       !          write(6,"('ttt: q  =',3f12.5)") q
       !          write(6,"('ttt: q_r=',3f12.5)") q_r
       allocate( pomat(nn,no) )
       read(ifpomat) pomat
       if( sum(abs(q-q_r))<1d-10) then ! .AND. kx <= nqibz ) then
          write(6,*) 'ok find the section for give qibz_k'
          exit
          !         elseif (kx >nqibz ) then
          !           exit
       endif
       deallocate(pomat)
    enddo
    !       if( sum(abs(q-q_r))>1d-10 ) then
    !          write(6,"('q  =',3f12.5)") q
    !          write(6,"('q_r=',3f12.5)") q_r
    !          stop 'POmat reading err q/=q_r'
    !       endif
    close(ifpomat)
  end subroutine readpomat
end module m_pomat
