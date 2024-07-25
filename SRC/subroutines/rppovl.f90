!> read PPOVLGG,PPOVLG,PPOVLI  ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
module m_read_ppovl
  use m_lgunit,only:stdo
  use m_nvfortran,only:findloc
  use m_ftox
  implicit none
  integer,public,protected:: nggg,ngcgp,ngcread,nxi,nxe,nyi,nye,nzi,nze,ngc2
  complex(8),public,protected,allocatable:: ppx(:,:),ggg(:),ppovlinv(:,:)
  integer,public,protected,allocatable:: ngvecc2(:,:),nvggg(:,:),nvgcgp2(:,:),ngvecc(:,:),igggi(:,:,:),igcgp2i(:,:,:)
  integer,public,protected:: nnxi,nnxe,nnyi,nnye,nnzi,nnze
  public:: getppx2
  private
  integer,private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx,nqini,nqnumt
  logical,private:: ppovlclosed=.true.,init=.true.
  integer,allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
  real(8),allocatable,private    :: qx_s(:,:)
  complex(8),allocatable,private :: ppx_s(:,:,:)
  logical,private:: debug=.false.
  real(8),allocatable,private:: qxtable(:,:)
  integer, private:: iex,gex
contains
  subroutine getppx2(qbas,qi,getngcgp) ! This return nvggg,nvgcgp2,ngvecc,  nggg,ngcgp,ngcread, ggg,ppovlinv
    real(8), intent(in)  ::qbas(3,3),qi(3)
    integer:: ngc, iqi,ippovlg,ippovli, ippovlginit
    integer:: access
    real(8)::qx(3)
    logical ::ippovlggooo=.true.
    integer:: ngcread2,ippovlgg
    character(3) :: charnum3
    integer:: verbose
!    logical:: init=.true.
    integer:: iqi0,igcgp2,iggg
    logical,optional:: getngcgp
    if(present(getngcgp)) then
       open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
       read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
       close(ippovlgg)
       return
    endif
    if(verbose()>=100) debug= .TRUE. 
    if(ippovlggooo) then !!  Make igggi inversion table
       open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
       read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
       if(debug)  write(stdo,"('Readin getppx2: nggg ngcgp nqq=',3i10)") nggg, ngcgp, nqq
       allocate(nvggg(1:3,1:nggg),ggg(1:nggg),nvgcgp2(1:3,ngcgp))
       read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
       read(ippovlgg) nvggg(1:3,1:nggg)
       read(ippovlgg) ggg(1:nggg)
       close(ippovlgg)
       nxi =minval(nvggg(1,1:nggg))
       nxe =maxval(nvggg(1,1:nggg))
       nyi =minval(nvggg(2,1:nggg))
       nye =maxval(nvggg(2,1:nggg))
       nzi =minval(nvggg(3,1:nggg))
       nze =maxval(nvggg(3,1:nggg))
       allocate( igggi(nxi:nxe,nyi:nye,nzi:nze),source= -100000)
       forall(iggg =1:nggg) igggi(nvggg(1,iggg),nvggg(2,iggg),nvggg(3,iggg)) = iggg
       nnxi = minval(nvgcgp2(1,1:ngcgp))      
       nnxe = maxval(nvgcgp2(1,1:ngcgp))
       nnyi = minval(nvgcgp2(2,1:ngcgp))
       nnye = maxval(nvgcgp2(2,1:ngcgp))
       nnzi = minval(nvgcgp2(3,1:ngcgp))
       nnze = maxval(nvgcgp2(3,1:ngcgp))
       allocate(igcgp2i(nnxi:nnxe,nnyi:nnye,nnzi:nnze),source= -100000)
       forall(igcgp2 =1:ngcgp) igcgp2i(nvgcgp2(1,igcgp2),nvgcgp2(2,igcgp2),nvgcgp2(3,igcgp2))=igcgp2 ! inversion table for nvgcgp2
       ippovlggooo=.false.
    endif
    if(init) then ! cache qx for finding a file for given qi. dec2017
       init=.false.
       allocate( qxtable(3,nqini:nqnumt) )
       do iqi = nqini,nqnumt
          open(newunit=ippovlginit,file="PPOVLG."//charnum3(iqi),form='unformatted')
          read(ippovlginit) qxtable(:,iqi) 
          close(ippovlginit) ! brought from outside of do iqi loop
       enddo
       if(debug) write(stdo,"('init ok!:should be done only once')")
    endif
    iqi=findloc([(sum(abs(qxtable(:,iqi0)-qi))<1d-10,iqi0=nqini,nqnumt)],value=.true.,dim=1)+nqini-1
    if(iqi<nqini) call rx('rppovl.F: qi is not found. some bug. qi='//ftof(qi))
    ! do iqi0 = nqini,nqnumt ! find file name (=charnum3(iqi)) for given qi.
    !    qx = qxtable(:,iqi0)
    !    if(sum(abs(qx-qi))<1d-10) then
    !       iqi = iqi0
    !       goto 1011
    !    endif
    ! enddo
    ! call rx('rppovl.F: qi is not found. some bug. qi='//ftof(qi))
    !1011 continue
    gex=access("PPOVLG."//charnum3(iqi),' ')
    iex=access("PPOVLI."//charnum3(iqi),' ')
    if(gex /= 0) call rx("PPOVLG."//charnum3(iqi)//" does not exist!") 
    if(iex /= 0) call rx("PPOVLI."//charnum3(iqi)//" does not exist!") 
    open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
    open(newunit=ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
    read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
    ngc = ngcread
    read(ippovli) qx, ngcread2 !, ngcx_s(iqi),ngc2_s(iqi)
    if(ngc==0) call rx('getppx2: can not find given qi ='//ftof(qi))
    if(sum(abs(qx-qi))>1d-10) call rx('getppx2: qx='//ftof(qx)//'.ne.'//ftof(qi))
    if(ngcread/=ngcread2)     call rx('rppovl.F: inconsistent PPOVLI PPOVLg')
    if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
    allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
    read(ippovlg) ngvecc(1:3,1:ngc)     !main do 1st
    read(ippovli) ppovlinv(1:ngc,1:ngc) !main do 2nd
    close(ippovlg)
    close(ippovli)
  end subroutine getppx2
end module m_read_ppovl
