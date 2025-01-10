!>Set itq for which we calculate self-energy. NTQXX mechanism to memorize this setting.
module m_itq
  use m_keyvalue,only: Getkeyvalue
  use m_readeigen,only: Readeval
  use m_genallcf_v3,only: nband
  implicit none
  public itq,ntq,setitq_hsfp0sc,setitq,nbandmx,setitq_hsfp0
  integer,allocatable,protected :: itq(:),nbandmx(:,:)
  integer,protected  :: ntq 
  private
contains
  subroutine setitq()
    integer:: i
    ntq = nband
    allocate(itq(ntq),source=[(i,i=1,ntq)])
  end subroutine setitq
  !!======================================================================
  subroutine setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx)
    use m_read_bzdata,only:qibz,nqibz
    use m_genallcf_v3,only: nspin
    use m_nvfortran,only: findloc
    intent(in)::            nbmx_sig,ebmx_sig,eftrue,nspinmx
    integer:: nbmx_sig,nspinmx
    real(8):: ebmx_sig,eftrue
    real(8),allocatable:: eqt(:)
    integer:: ntqxx,is,ip,iband,i
    allocate(nbandmx(nqibz,nspinmx))     !!   Get nbandmx(iq,isp)
    allocate(eqt(nband))
    do is = 1,nspin
      do ip = 1,nqibz
        eqt= readeval(qibz(1,ip),is)
        nbandmx(ip,is) = min(findloc(eqt-eftrue>ebmx_sig,value=.true.,dim=1)-1, nbmx_sig)
      enddo
    enddo
    ntq = maxval(nbandmx)+1
    allocate(itq(ntq),source=[(i,i=1,ntq)]) !itq is used also in hsfp0.m.F ! trivial case of itq itq(i)=i
  end subroutine setitq_hsfp0sc
  subroutine setitq_hsfp0 (ngcmx_in,ngpmx_in,tote,nbmin,nbmax,noccxv)
    intent(in)::           ngcmx_in,ngpmx_in,tote,nbmin,nbmax,noccxv
    integer:: nband_in,ngcmx_in,ngpmx_in,ifqpnt,noccxv,i,nss(2),nbmax,nbmin
    logical:: tote
    if (tote) then
       ntq = noccxv
       allocate( itq(ntq),source=[(i,i=1,ntq)])
    else
       ntq = nbmax-nbmin+1
       allocate( itq, source=[(i,i=nbmin,nbmax)])
    endif
  end subroutine setitq_hsfp0
end module m_itq

