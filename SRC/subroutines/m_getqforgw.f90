!> Set q and range of band index i for calculating the self-energy \Sigma(q,i). Range of band index is not q dependent
module m_getQforGW 
  use m_read_bzdata,only:  nqibz,qibz
  use m_keyvalue,only: getkeyvalue
  use m_ftox
  use m_nvfortran
  integer,public,protected :: nbmin,nbmax,nq
  real(8),public,allocatable,protected:: qx(:,:)
  public::getqforgw, getqonly
  private
contains
  subroutine getqforgw(lqall)
    use m_readefermi,only: readefermi,ef
    use m_genallcf_v3,only: nspin, nband
    use m_readeigen,only: readeval !init* is alreaday called.
    logical:: lqall
    integer:: ret,nnx,ip,is,k,ifqpnt,nnm
    real(8),allocatable::eqt(:)
    real(8)::rydberg,ecut,emin
    if(lqall) then !gwsc mode. Set all bands and nqibz
       nbmax=nband
       nbmin=1
       nq=nqibz
       allocate(qx(3,nq),source= qibz(1:3,1:nq))
       return
    endif
    call getqonly() !Set qx(1:3,nq)
    call readefermi()
    call getkeyvalue("GWinput","EMAXforGW",ecut,default= 99999d0)
    call getkeyvalue("GWinput","EMINforGW",emin,default=-99999d0)
    allocate(eqt(nband))
    nnx=-999999
    nnm= 999999
    do ip=1,nq
       do is=1,nspin
          eqt = readeval(qx(1:3,ip),is)-ef
          nnx= max(findloc(eqt>ecut/rydberg(),value=.true.,dim=1),nnx)
          nnm= min(findloc(eqt>emin/rydberg(),value=.true.,dim=1),nnm)
       enddo
    enddo
!    write(6,*)'vvvvvvvvvvvvvvvvv nq=',nq
    deallocate(eqt)
    nbmax=merge(nband,nnx-1,nnx==-999999)
    nbmin=nnm
  end subroutine getqforgw
  subroutine getqonly()
    integer:: ret,nnx,ip,is,k,ifqpnt
    logical:: ibzqq
    real(8),allocatable::eqt(:)
    real(8)::qq(3)
    call getkeyvalue("GWinput","<QforGW>",unit=ifqpnt,errstop='off',status=ret)
    nq=0
    if(ret>0) then !read <QforGW> section
       k=0
       do 
          k=k+1
          read(ifqpnt,*,end=1012,err=1012) qq(1:3)
          write(6,*)k,qq
       enddo
1012   continue
       nq=k-1
       close(ifqpnt)
    endif
    write(6,*)' Readin from QforGW :nq=',nq
    call getkeyvalue("GWinput","QforGWIBZ",ibzqq,default=.false.)
    if(nq==0.or.ibzqq) then
       nq=nqibz
       allocate(qx(3,nq),source= qibz(1:3,1:nq))
    else    
       call getkeyvalue("GWinput","<QforGW>",unit=ifqpnt,errstop='OFF',status=ret)
       allocate(qx(3,nq))
       do k=1,nq
          read(ifqpnt,*) qx(1:3,k)
       enddo
       if(ret>0) close(ifqpnt)
    endif
  end subroutine getqonly
endmodule m_getQforGW
