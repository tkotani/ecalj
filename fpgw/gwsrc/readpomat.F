C ... read pomat
      module m_pomat
c------Return eigenvalus and eigenfunctions for given q and isp.
c Routines are:
C   readeval: ev (eigenvalue)
c
      implicit none
      complex(8),allocatable:: pomat(:,:)
      integer(4):: nn,no
c-----------------------------------------------------------------
      contains

c------------------------------
      subroutine readpomat(q)
! This returns nn,no, and pomat(nn,no)
      implicit none
      real(8):: q_r(3),q(3)
      integer(4):: nn,iqx,isx
      integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,iopen,iclose
c      complex(8),allocatable:: pomat(:,:)
c      if( allocated(pomat) ) deallocate(pomat)
      ifpomat = iopen('POmat',0,-1,0) !oct2005
C... smoothed mixed basis !oct2005
C This replace original zmelt with new zmelt based on smoothed mixed basis.
      do
        read(ifpomat) q_r,nn,no,iqx !readin reduction matrix pomat
c          write(6,"('ttt: q  =',3f12.5)") q
c          write(6,"('ttt: q_r=',3f12.5)") q_r
        allocate( pomat(nn,no) )
        read(ifpomat) pomat
        if( sum(abs(q-q_r))<1d-10) then ! .and.kx <= nqibz ) then
          write(6,*) 'ok find the section for give qibz_k'
          exit
!         elseif (kx >nqibz ) then
!           exit
        endif
        deallocate(pomat)
      enddo
c       if( sum(abs(q-q_r))>1d-10 ) then
c          write(6,"('q  =',3f12.5)") q
c          write(6,"('q_r=',3f12.5)") q_r
c          stop 'POmat reading err q/=q_r'
c       endif
      isx = iclose('POmat')
      end subroutine
      end module
