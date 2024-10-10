module m_getQforGW !reading q and band index for gw_lmfh !2024-oct
  use m_read_bzdata,only: read_bzdata, nqibz,qibz
  use m_keyvalue,only:getkeyvalue
  use m_readefermi,only: readefermi,ef
  use m_ftox
  integer,public:: nbmin,nbmax,nq
  real(8),public,allocatable:: qx(:,:)
  public::getqforgw,getqonly
  private
contains
  subroutine getqforgw(lqall)
    use m_genallcf_v3,only: nspin, nband
    use m_readeigen,only: readeval !init* is alreaday called.
    logical:: lqall
    integer:: ret,nnx,ip,is,k,ifqpnt
    real(8),allocatable::eqt(:)
    real(8)::rydberg,ecut
    call readefermi()
    call getkeyvalue("GWinput","EMAXforGW",ecut,default=5d0)
    if(lqall.or.ret==-1) ecut=1d6
    nbmin=1
    allocate(eqt(nband))
    nnx=-999999
    do ip=1,nqibz
       do is=1,nspin
          eqt = readeval(qibz(1:3,ip),is)-ef
          nnx= max(findloc(eqt>ecut/rydberg(),value=.true.,dim=1),nnx)
       enddo
    enddo
    deallocate(eqt)
    nbmax=nnx-1
    call getqonly()
  end subroutine getqforgw
  subroutine getqonly()
    integer:: ret,nnx,ip,is,k,ifqpnt
    real(8),allocatable::eqt(:)
    real(8)::qq(3)
    call getkeyvalue("GWinput","<QforGW>",unit=ifqpnt,errstop='OFF',status=ret)
    k=0
    do 
       k=k+1
       read(ifqpnt,*,end=1012,err=1012) qq(1:3)
       write(6,*)k,qq
    enddo
1012 continue
    nq=k-1
    write(6,*)' Readin from QforGW :nq=',nq
    close(ifqpnt)
    if(nq==0) then
       nq=nqibz
       allocate(qx(3,nq),source= qibz(1:3,1:nq))
    else    
       call getkeyvalue("GWinput","<QforGW>",unit=ifqpnt,errstop='OFF',status=ret)
       allocate(qx(3,nq))
       do k=1,nq
          read(ifqpnt,*) qx(1:3,k)
       enddo
       close(ifqpnt)
    endif
  end subroutine getqonly
endmodule m_getQforGW
