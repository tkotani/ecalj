      module m_itq
      use m_keyvalue,only: Getkeyvalue
      use m_readeigen,only: Readeval
      integer,allocatable,protected :: itq(:),nbandmx(:,:)
      integer,protected  :: ntq,ngcmx,ngpmx, nband

!!======================================================================
      contains
      subroutine setnband(nband_in)
      integer:: nband_in
      nband=nband_in
      end

      subroutine setitq(ntq_in,nband_in,ngcmx_in,ngpmx_in)
      implicit none
      integer,intent(in):: nband_in,ngcmx_in,ngpmx_in,ntq_in
      integer:: i
      ntq   = ntq_in
      nband = nband_in
      ngcmx = ngcmx_in
      ngpmx = ngpmx_in
      allocate(itq(ntq))
      do i=1,ntq
        itq(i)=i
      enddo
      end subroutine setitq
!!======================================================================
      subroutine setitq_hsfp0sc(qibz,nqibz,nq,nspin,nbmx,ebmx,eftrue,nspinmx)
      implicit none
      integer,intent(in):: nqibz,nq,nspin,nbmx(2),nspinmx
      real(8),intent(in):: qibz(3,nqibz),ebmx(2),eftrue
      real(8),allocatable:: eqt(:)
      integer:: ifih,ntqxx,is,ip,iband,i,nband_r,nq_r
      logical:: lntq
      inquire(file='NTQXX',EXIST=lntq)
      open(newunit=ifih,file='NTQXX')
!! Get ntq
      if(lntq) then
         read(ifih,*) nband_r,nq_r,ntq
         if(nband_r/=nband.or.nq_r/=nq) then
            rewind ifih
            lntq=.false.
         endif  
      else
         ntq=0
         allocate(eqt(nband))
         do is = 1,nspin
            do ip = 1,nq
               call READEVAL(qibz(1,ip),is, eqt)
               do iband=1,nband
                  ntq = max(iband,ntq)
                  if(eqt(iband)-eftrue>ebmx(1)) exit
               enddo  
            enddo
         enddo
         ntq = min(ntq, nbmx(1))
         deallocate(eqt)
         write(ifih,"(3i10)") nband,nq,ntq
      endif
!! Determine ntq.  See also in sxcf_fal.sc.F ntq should be common for all ixc modes.
!! FIX NTQ during iteration by the file NTQ 15jun2015
!!
!! Determine nbandmx. Moved from sxcf_fal2.sc.F.
!!!! count number of band to calculate.
!! I think it it better to determine nbandmx in a manner within LDA 
!! (need to care degeneracy...).
      allocate(nbandmx(nq,nspinmx))
!!   Get nbandmx(iq,isp)
      allocate(eqt(nband))
      do is = 1,nspinmx
         do ip = 1,nq   
            call READEVAL(qibz(1,ip),is, eqt)
            if(lntq) then
               read(ifih,*) ntqxx ! ntqxx = ntq !jun2016
            else  
               ntqxx = 0
               do i = 1,ntq
                  if(eqt(i)-eftrue<ebmx(1)) ntqxx =ntqxx  + 1
               enddo
               ntqxx = min(ntqxx, nbmx(1))
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
      end subroutine
!!======================================================================
      subroutine setitq_hsfp0 (nband_in,ngcmx_in,ngpmx_in,tote,ifqpnt,noccxv,nss)
      integer,intent(in):: nband_in,ngcmx_in,ngpmx_in,ifqpnt,noccxv
      logical,intent(in):: tote
      integer:: i,nss(2)
      nband = nband_in
      ngcmx = ngcmx_in
      ngpmx = ngpmx_in
      if (tote) then
        ntq    = noccxv
        allocate( itq(ntq) )
        do i   = 1,ntq
          itq(i) = i
        enddo
      else
        read (ifqpnt,*) ntq
        allocate( itq(ntq) )
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

      end module 
