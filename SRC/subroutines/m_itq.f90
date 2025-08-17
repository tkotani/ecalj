!>Set itq for which we calculate self-energy. --ntqxx mechanism to memorize this setting. --ntqxx may stabilize QSGW iteration.
module m_itq
  use m_keyvalue,only: Getkeyvalue
  use m_readeigen,only: Readeval
  use m_genallcf_v3,only: nband
  use m_lgunit,only: stdo
  implicit none
  public itq,ntq,setitq_hsfp0sc,setitq,nbandmx,setitq_hsfp0
  integer,allocatable,protected :: itq(:),nbandmx(:,:)
  integer,protected  :: ntq
  logical,private:: rntq=.false.
  private
contains
  subroutine setitq()
    integer:: i
    ntq = nband
    allocate(itq,source=[(i,i=1,ntq)])
  end subroutine setitq
  subroutine setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx)
    use m_read_bzdata,only:qibz,nqibz
    use m_nvfortran,only: findloc
    use m_ftox
    use m_mpi,only: mpi__root
    intent(in)::            nbmx_sig,ebmx_sig,eftrue,nspinmx
    integer:: nbmx_sig,nspinmx,ifih,nspinmxin
    integer:: ntqxx,is,ip,iband,i,nqibzin,nspinin,iqibz,ierr
    real(8):: ebmx_sig,eftrue
    real(8),allocatable:: eqt(:)
    logical:: cmdopt0,readntqxx
    if(rntq) return
    allocate(nbandmx(nqibz,nspinmx),eqt(nband))
    readntqxx=.false.
    if(cmdopt0('--ntqxx')) then !NTQXX is to keep the same number of bands for Sigma for each k during iteration.
      open(newunit=ifih,file='NTQXX',status='old',iostat=ierr)
      if(ierr==0) then
        read(ifih,*) nqibzin,nspinmxin
        if(nqibzin==nqibz.or.nspinmxin==nspinmx) readntqxx=.true.
      endif
    endif   
    if(readntqxx) then
      do iqibz=1,nqibz
         read(ifih,*) nbandmx(iqibz,:)
      enddo
      close(ifih)
      rntq=.true.
    else  
      do is = 1,nspinmx
          do ip = 1,nqibz
             eqt= readeval(qibz(1,ip),is)
             nbandmx(ip,is) = min(findloc(eqt-eftrue>ebmx_sig,value=.true.,dim=1)-1, nbmx_sig)
          enddo
      enddo
      write(stdo,*)'vvvvvvvvvvv222 ntqxx=',cmdopt0('--ntqxx'),mpi__root
      if(mpi__root.and.cmdopt0('--ntqxx')) then
          open(newunit=ifih,file='NTQXX')
          write(ifih,ftox) nqibz,nspinmx,' !nqibz nspinmx. Note NTQXX is used when --ntqxx'
          do iqibz=1,nqibz
             write(ifih,ftox) nbandmx(iqibz,:),' !=nbandmx ', iqibz,' ! =iqibz'
          enddo
          close(ifih)
      endif
    endif
    ntq = maxval(nbandmx)+1
    allocate(itq,source=[(i,i=1,ntq)]) ! trivial case of itq itq(i)=i
  end subroutine setitq_hsfp0sc
  subroutine setitq_hsfp0 (ngcmx_in,ngpmx_in,tote,nbmin,nbmax,noccxv)
    intent(in)::           ngcmx_in,ngpmx_in,tote,nbmin,nbmax,noccxv
    integer:: nband_in,ngcmx_in,ngpmx_in,ifqpnt,noccxv,i,nss(2),nbmax,nbmin
    logical:: tote
    if (tote) then
       ntq = noccxv
       allocate( itq, source=[(i,i=1,ntq)])
    else
       ntq = nbmax-nbmin+1
       allocate( itq, source=[(i,i=nbmin,nbmax)])
    endif
  end subroutine setitq_hsfp0
end module m_itq

