!>Set itq for which we calculate self-energy. And NTQX mechanism
module m_itq
  use m_keyvalue,only: Getkeyvalue
  use m_readeigen,only: Readeval
  use m_readhbe,only: nband
  implicit none
  public itq,ntq,setitq_hsfp0sc,setitq,nbandmx,setitq_hsfp0
  integer,allocatable,protected :: itq(:),nbandmx(:,:)
  integer,protected  :: ntq 
  private
contains
  subroutine setitq()
    integer:: i
    ntq   = nband
    allocate(itq(ntq))
    do i=1,ntq
       itq(i)=i
    enddo
  end subroutine setitq
  !!======================================================================
  subroutine setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx)
    use m_read_bzdata,only:qibz,nqibz
    use m_genallcf_v3,only: nspin
    intent(in)::            nbmx_sig,ebmx_sig,eftrue,nspinmx
    integer:: nbmx_sig,nspinmx
    real(8):: ebmx_sig,eftrue
    real(8),allocatable:: eqt(:)
    integer:: ifih,ntqxx,is,ip,iband,i,nband_r,nq_r
    logical:: lntq
    inquire(file='NTQXX',EXIST=lntq) !NTQXX is to keep the same number of bands during iteration.
    open(newunit=ifih,file='NTQXX')
    !! Get ntq
    if(lntq) then
       read(ifih,*) nband_r,nq_r,ntq
       if(nband_r/=nband .OR. nq_r/=nqibz) then
          rewind ifih
          lntq=.false.
       endif
    else
       ntq=0
       allocate(eqt(nband))
       do is = 1,nspin
          do ip = 1,nqibz
             eqt= readeval(qibz(1,ip),is)
             do iband=1,nband
                ntq = max(iband,ntq)
                if(eqt(iband)-eftrue>ebmx_sig) exit
             enddo
          enddo
       enddo
       ntq = min(ntq, nbmx_sig)
       deallocate(eqt)
       write(ifih,"(3i10)") nband,nqibz,ntq
    endif
    !! Determine ntq.  See also in sxcf_fal.sc.F ntq should be common for all ixc modes.
    !! FIX NTQ during iteration by the file NTQ 15jun2015
    !!
    !! Determine nbandmx. Moved from sxcf_fal2.sc.F.
!!!! count number of band to calculate.
    !! I think it it better to determine nbandmx in a manner within LDA
    !! (need to care degeneracy...).
    allocate(nbandmx(nqibz,nspinmx))
    !!   Get nbandmx(iq,isp)
    allocate(eqt(nband))
    do is = 1,nspinmx
       do ip = 1,nqibz
          eqt= READEVAL(qibz(1,ip),is)
          if(lntq) then
             read(ifih,*) ntqxx ! ntqxx = ntq !jun2016
          else
             ntqxx = 0
             do i = 1,ntq
                if(eqt(i)-eftrue<ebmx_sig) ntqxx =ntqxx  + 1
             enddo
             ntqxx = min(ntqxx, nbmx_sig)
             write(ifih,"(i10)") ntqxx
          endif
          if(ntqxx<nband) then ! redudce ntqxx when band tops are degenerated.
             do i=ntqxx,1,-1
                if(eqt(i+1)-eqt(i)<1d-2) then !1d-2 is a tol to check degeneracy.
                   ntqxx=i-1
                else
                   exit
                endif
             enddo
          endif
          nbandmx(ip,is) = ntqxx !number of bands to be calculated
       enddo
    enddo
    deallocate(eqt)
    close(ifih)
    !! trivial case of itq itq(i)=i
    allocate (itq(ntq))
    do i = 1, ntq
       itq(i) = i !itq is used also in hsfp0.m.F
    enddo
  end subroutine setitq_hsfp0sc
  !!======================================================================
  subroutine setitq_hsfp0 (ngcmx_in,ngpmx_in,tote,ifqpnt,noccxv,nss)
    intent(in)::     ngcmx_in,ngpmx_in,tote,ifqpnt,noccxv,nss
    integer:: nband_in,ngcmx_in,ngpmx_in,ifqpnt,noccxv,i,nss(2)
    logical:: tote
    if (tote) then
       ntq = noccxv
       allocate( itq(ntq),source=[(i,i=1,ntq)])
    else
       read (ifqpnt,*) ntq
       allocate( itq(ntq)) 
       read (ifqpnt,*) (itq(i),i=1,ntq)
    endif
    if(nss(2)/=-99997) then
       if(allocated(itq)) deallocate(itq)
       ntq=nss(2)-nss(1)+1
       allocate( itq(ntq) )
       do i=max(1,nss(1)),min(nss(2),nband)
          itq(i-nss(1)+1) = i
       enddo
    endif
  end subroutine setitq_hsfp0
end module m_itq

!!------------------------------------------------------------------
! module m_selectqp
!   use m_keyvalue,only: GETKEYVALUE
!   implicit none
!   integer,protected:: nq
!   real(8),allocatable,protected ::qvec(:,:)
! contains
!   subroutine getqpoint(selectqp,nqibz,qibz)
!     real(8):: qibz(3,nqibz)
!     logical:: lqall,laff,selectqp
!     integer:: ifqpnt,ret,nqibz,i,iaf,iq,iqall,k
!     if(selectqp) then 
!        call GETKEYVALUE("GWinput","<QPNT>",unit=ifqpnt,status=ret)
!        lqall      = .false.
!        laff        = .false.
!        call readx   (ifqpnt,10)
!        read (ifqpnt,*) iqall,iaf
!        if (iqall == 1) lqall = .TRUE. 
!        if (iaf   == 1) laff = .TRUE. 
!        call readx   (ifqpnt,100)
!        if (lqall) then         !all q-points case
!           nq         = nqibz
!           allocate(qvec(3,nq))
!           call dcopy   (3*nqibz,qibz,1,qvec,1)
!        else
!           call readx   (ifqpnt,100)
!           read (ifqpnt,*) nq
!           allocate(qvec(3,nq))
!           do       k = 1,nq
!              read (ifqpnt,*) i,qvec(1,k),qvec(2,k),qvec(3,k)
!           enddo
!        endif
!        close(ifqpnt)
!     else
!        nq = nqibz
!        allocate(qvec(3,nq))
!        qvec(:,1:nq) = qibz(:,1:nq)
!     endif
!     do iq=1,nq
!        write(6,'(" Target iq q=",i6,3f9.4)')iq,qvec(:,iq)
!     enddo
!   end subroutine getqpoint
!end module m_selectqp
