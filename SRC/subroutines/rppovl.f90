!> read PPOVLGG,PPOVLG,PPOVLI
!!  ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
module m_read_ppovl
  implicit none
  integer,protected:: ngc2
  complex(8),protected,allocatable :: ppx(:,:)
  integer,protected,allocatable :: ngvecc2(:,:)
  complex(8),protected,allocatable:: ggg(:),ppovlinv(:,:)
  integer,protected,allocatable:: nvggg(:,:),nvgcgp2(:,:),ngvecc(:,:)
  integer,protected:: nggg,ngcgp,ngcread
  integer,allocatable,protected:: igggi(:,:,:),igcgp2i(:,:,:)
  integer,protected:: nxi,nxe,nyi,nye,nzi,nze
  !     !
  integer(4),private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx, ngcgpmx,nqini,nqnumt
  logical,private:: ppovlclosed=.true.,init=.true.
  integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
  real(8),allocatable,private    :: qx_s(:,:)
  complex(8),allocatable,private :: ppx_s(:,:,:)
  logical,private:: debug=.false.
  real(8),allocatable,private:: qxtable(:,:)
  integer(4), private:: loopnum = 0, iex,gex
  integer,private:: nnxi,nnxe,nnyi,nnye,nnzi,nnze
contains
  subroutine getppx2(qbas,qi)
    !! This return nvggg,nvgcgp2,ngvecc,  nggg,ngcgp,ngcread, ggg,ppovlinv
    real(8), intent(in)  ::qbas(3,3),qi(3)
    integer(4):: ngc, iqi,ippovlg,ippovli, ippovlginit
    integer(4):: access,ifile_handle
    real(8)::qx(3)
    logical ::ippovlggooo=.true.
    integer:: ngcread2,ippovlgg
    character(3) :: charnum3
    integer:: verbose
    logical:: init=.true.
    integer:: iqi0,igcgp2,iggg
    if(verbose()>=100) debug= .TRUE. 
    if(ippovlggooo) then
       open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
       read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
       write(6,"('Readin getppx2: nggg ngcgp nqq=',3i10)") nggg, ngcgp, nqq
       allocate(nvggg(1:3,1:nggg),ggg(1:nggg),nvgcgp2(1:3,ngcgp))
       read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
       read(ippovlgg) nvggg(1:3,1:nggg)
       read(ippovlgg) ggg(1:nggg)
       close(ippovlgg)
       !!  Make igggi inversion table
       nxi =minval(nvggg(1,1:nggg))
       nxe =maxval(nvggg(1,1:nggg))
       nyi =minval(nvggg(2,1:nggg))
       nye =maxval(nvggg(2,1:nggg))
       nzi =minval(nvggg(3,1:nggg))
       nze =maxval(nvggg(3,1:nggg))
       allocate( igggi(nxi:nxe,nyi:nye,nzi:nze))
       igggi = -100000
       do iggg =1,nggg
          igggi(nvggg(1,iggg),nvggg(2,iggg),nvggg(3,iggg)) = iggg
       enddo
       !     ! inversion table for nvgcgp2
       nnxi = minval(nvgcgp2(1,1:ngcgp))
       nnxe = maxval(nvgcgp2(1,1:ngcgp))
       nnyi = minval(nvgcgp2(2,1:ngcgp))
       nnye = maxval(nvgcgp2(2,1:ngcgp))
       nnzi = minval(nvgcgp2(3,1:ngcgp))
       nnze = maxval(nvgcgp2(3,1:ngcgp))
       allocate(igcgp2i(nnxi:nnxe,nnyi:nnye,nnzi:nnze))
       igcgp2i = -100000
       do igcgp2 =1,ngcgp
          igcgp2i(nvgcgp2(1,igcgp2),nvgcgp2(2,igcgp2),nvgcgp2(3,igcgp2))=igcgp2
       enddo
       ippovlggooo=.false.
    endif
    !! cache qx for finding a file for given qi. dec2017
    if(init) then
       init=.false.
       allocate( qxtable(3,nqini:nqnumt) )
       ippovlginit=ifile_handle()
       do iqi = nqini,nqnumt
          open(ippovlginit,file="PPOVLG."//charnum3(iqi),form='unformatted')
          read(ippovlginit) qx
          qxtable(:,iqi) = qx
          close(ippovlginit) ! brought from outside of do iqi loop
       enddo
       loopnum=0
       write(6,"('init ok!:should be done only once')")
    endif
    !! find file name (=charnum3(iqi)) for given qi.
    do iqi0 = nqini,nqnumt
       qx = qxtable(:,iqi0)
       if(sum(abs(qx-qi))<1d-10) then
          iqi = iqi0
          goto 1011
       endif
    enddo
    write(6,"('nnnnnnq ',3f10.5)") qi
    call rx('rppovl.F: qi is not found. some bug.')
1011 continue
    !! read file of iqi. iqi is determined for given qx.
    loopnum=loopnum+1
    if(loopnum == 1) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum
    if(loopnum == 10) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum
    if(loopnum == 100) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum
    if(loopnum == 1000) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum
    if(loopnum == 5000) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum
    if(loopnum == 10000) write(6,"('iqi=,loop num=',i10,i10)") iqi,loopnum

    gex=access("PPOVLG."//charnum3(iqi),' ')
    iex=access("PPOVLI."//charnum3(iqi),' ')
    if(gex /= 0)then
       write(6,"('PPOVLG.00... does not exist! in iqi=)',i4,'( in loop ',i4)") iqi,loopnum
       call rx('some PPOLVG. file does not exist')
    endif
    if(iex /= 0)then
       write(6,"('PPOVLI.00... does not exist! in iqi=)',i4,'( in loop ',i4)") iqi,loopnum
       call rx('some PPOLVI. file does not exist')
    endif

    open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
    ippovli=ifile_handle()
    open(ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
    read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
    ngc = ngcread
    read(ippovli) qx, ngcread2 !, ngcx_s(iqi),ngc2_s(iqi)
    !! sanity checkcs
    if(ngc==0) then
       write(6,"('qi qx=',3f13.5,3x,3f13.5)") qi,qx
       call rx('getppx2: can not find given qi')
    endif
    if(sum(abs(qx-qi))>1d-10) then
       write(6,"('nnnnnnqiqx ',3f10.5,2x,3f10.5)") qi,qx
       write(6,"('nnnnfile=',a)")"PPOVLG."//charnum3(iqi)
       call rx('getppx2: qx\ne qi')
    endif
    if(ngcread/=ngcread2) call rx('rppovl.F: inconsistent PPOVLI PPOVLg')
    !! main do for ppovlg and ppovli
    if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
    allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
    read(ippovlg) ngvecc(1:3,1:ngc)     !main do 1st
    read(ippovli) ppovlinv(1:ngc,1:ngc) !main do 2nd
    close(ippovlg)
    close(ippovli)
  end subroutine getppx2
end module m_read_ppovl