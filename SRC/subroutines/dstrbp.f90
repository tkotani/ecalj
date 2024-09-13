module m_dstrbp
  implicit none
  public dstrbp
  private
contains
  subroutine dstrbp(ndatx,np,nblk,index) !rank divider
    ! i=1,ndatx is divided into [index(procid),index(procid+1)-1] for each procid=0:np-1
    ! Note that index(procid)=ndatx+1,
    use m_lgunit,only:stdo
    integer:: ndatx,ndiv,nloop,np,nblk,index(0:np),ip,iqini
    if(nblk/=1) call rx('we allow only nblk=1 now')
    ndiv= ndatx/np !MPI division
    if(ndatx>ndiv*np) ndiv=ndiv+1
    !print *,'nnndiv=',ndiv
    do ip=0,np
       iqini = ndiv*ip+1
       if(iqini> ndatx) iqini=ndatx+1
       index(ip) = iqini
    enddo
  end subroutine dstrbp
end module m_dstrbp
