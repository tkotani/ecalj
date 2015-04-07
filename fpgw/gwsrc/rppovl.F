!> read PPOVLGG,PPOVLG,PPOVLI 
!!  ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
      module m_read_ppovl
      implicit none
      integer(4):: ngc2
      complex(8),allocatable :: ppx(:,:)
      integer(4),allocatable :: ngvecc2(:,:)
      complex(8),allocatable:: ggg(:),ppovlinv(:,:)
      integer,allocatable:: nvggg(:,:),nvgcgp2(:,:),ngvecc(:,:)
      integer:: nggg,ngcgp,ngcread

      integer(4),private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx, ngcgpmx,nqini,nqnumt
      logical,private:: ppovlclosed=.true.,init=.true.
      integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
      real(8),allocatable,private    :: qx_s(:,:)
      complex(8),allocatable,private :: ppx_s(:,:,:)
      logical,private:: debug=.false.

      contains
      subroutine getppx2(qbas,qi) 
!! This return nvggg,nvgcgp2,ngvecc,  nggg,ngcgp,ngcread, ggg,ppovlinv
      implicit none
      real(8), intent(in)  ::qbas(3,3),qi(3)

      integer(4):: ngc, iqi,ippovlg = 9221,ippovli = 9222
      real(8)::qx(3)
      logical ::ippovlggooo=.true.
      integer:: ippovlgg=9223,ngcread2
      character*3 charnum3
      integer:: verbose
      if(verbose()>=100) debug=.true.
      if(ippovlggooo) then
         open(ippovlgg,file= "PPOVLGG",form='unformatted')
         read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
         write(6,"('Readin getppx2: nggg ngcgp nqq=',3i10)") nggg, ngcgp, nqq
         allocate(nvggg(1:3,1:nggg),ggg(1:nggg),nvgcgp2(1:3,ngcgp))
         read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
         read(ippovlgg) nvggg(1:3,1:nggg)
         read(ippovlgg) ggg(1:nggg)
         close(ippovlgg)
         ippovlggooo=.false.
      endif
      do iqi = nqini,nqnumt
         open(ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
         read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
         if(debug) print *,'getppx2: iqi qi=',iqi,qi,qx
         ngc = ngcread
         if(sum(abs(qx-qi))>1d-10) then
           close(ippovlg)
         else
           open(ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
           read(ippovli) qx, ngcread2 !, ngcx_s(iqi),ngc2_s(iqi)
           if(ngcread/=ngcread2) call rx('rppovl.F: inconsistent PPOVLI PPOVLg')
           if(ngc==0) exit
           if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
           allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
           read(ippovlg) ngvecc(1:3,1:ngc)
           read(ippovli) ppovlinv(1:ngc,1:ngc)
           close(ippovlg)
           close(ippovli)
           goto 1012
         endif   
      enddo
      write(6,"('qi qx=',3f13.5,3x,3f13.5)") qi,qx
      call rx('getppx2: can not find given qi')
 1012 continue
      end subroutine getppx2
      end module
