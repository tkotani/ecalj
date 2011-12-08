      module m_rsibl_ev
      contains
Ckino from rsibl write eigenvectors
CSFCPP#define F90 1
      subroutine rsibl_ev(ssite,sspec,slat,nbas,isp,q,iq,ndimh,nspc,
     .napw,igapw,iprmb,nevec,evec,k1,k2,k3, n_eiglist,eiglist)
      use m_struc_def  !Cgetarg
      use m_globalvariables
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
      use m_w_psir
C- Add smooth part of output density into smrho and forces.
C ----------------------------------------------------------------------
Ci Inputs
Ci  x lfrce :if nonzero, accumulate contribution to force
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: rsibl1
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: ngcut
Ci     Stored:    *
Ci     Passed to: tbhsi rsibl1 uspecb
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat gmax nabc ng ogv okv vol
Ci     Stored:    *
Ci     Passed to: *
Ci   lfrce :1 calculate contribution to forces
Ci   nbas  :size of basis
Ci   q     :Bloch vector
Ci   iq    :index to current k-point
Ci   ndimh :dimension of hamiltonian
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial Fermi levels
Ci   nevec :number of eigenvectors with nonzero weights
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   k1..3 :dimensions smpot,smrho
Ci   smpot :smooth potential on uniform mesh, needed for forces
Co Outputs
Co   smrho :smooth density accumulated for this qp
Co   f     :force contribution accumulated for this qp
Cr Remarks
Cm MPI
Cm   Parallelise over the eigenvector loop. The vector block size is
Cm   chosen (in the range 6-16, by dstrbp.f) so as to distribute the
Cm   work optimally across processes. Two work arrays of the size of
Cm   smrho are allocated from the heap as buffers. Only one will be
Cm   needed under MPI-2. See comments in hsibl.
Cb Bugs
Cb    replace call to gvgvcomp and pass ipv as input
Cb    The non-F90 version should work, but it is no longer tested
Cu Updates
Cu   29 Dec 08 Unsuccessful attempt to make work with openmp
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu   10 Sep 06 Added MPI parallelization in the spin-coupled case
Cu   23 Dec 04 Extended to spin-coupled case
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   12 Oct 00 Use q-dependent list of G vectors
Cu    6 Jul 00 attempt to vectorize by grouping eigenvectors in blocks
Cu   17 Jun 00 Spin polarized
Cu   23 May 00 Adapted from nfp rsif_q.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI
      include "mpif.h"
#if MPE
      include "mpef.h"
#endif
      integer ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1),i
      double precision starttime, endtime
      logical mlog,cmdopt
      integer lgunit
      character*120 strn
#endif
#if MPE
C Event numbers:
      include "events.ins"
#endif
      integer procid, master, nproc, mpipid
      integer nbas,isp,k1,k2,k3,ndimh,nevec,iprmb(*),iq,nspc
      integer,parameter:: numq=1
      integer napw,igapw(3,napw)
      real(8):: q(3) 
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex evec(ndimh,nspc,nevec)
      integer:: n_eiglist
      integer:: eiglist(n_eiglist)

C ... Local parameters
      integer n0,nkap0,nermx,npmx,nblk,nlmto
      parameter (n0=10,nkap0=3,nermx=100,npmx=128)
        integer:: nspec , ngabc(3) , n1 , n2 , n3 , nrt , net , ng , 
     .  nglob , ltop , nlmtop , ogq , og2 , ohe , ohr , oyl , oylw , 
     .  oiv , iprint 
       integer,allocatable :: iv_p_okv(:) 

ckino Dec.8.2011:         real(8),allocatable :: rv_p_ogv(:) 
       real(8),allocatable :: rv_a_ogv(:) 

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
      double precision alat,qlat(3,3),plat(3,3),q0(3),gmax,xx
C     Shared variables (openmp)
      double precision vol
      double precision etab(nermx),rtab(nermx)
      integer ivec,nvec
C     Local for parallel threads
      integer,allocatable:: ivp(:)
      complex(8),allocatable::psi(:,:,:),psir(:,:,:),vpsi(:,:,:),
     .wk(:,:,:)
      real(8),allocatable:: cosi(:),sini(:),wk2(:)

#if MPI
      integer, dimension(:), allocatable :: vproc
      integer nvsub
      integer osmbuf,ofrbuf
      complex(8),allocatable:: xsmrho(:,:,:,:,:)
      real(8),allocatable:: fr(:,:,:)
#endif

#if OPENMP
      complex(8),allocatable:: xsmrho(:,:,:,:,:)
      real(8),allocatable:: fr(:,:,:)
#endif

      integer:: ivecini,ivecend
      integer,allocatable:: w_oiv(:)
      real(8),allocatable:: w_ogq(:),w_oyl(:),w_oylw(:),w_og2(:),w_ohe(:),w_ohr(:)
      complex(8),allocatable:: w_osmbuf(:)
      real(8),allocatable:: w_ofrbuf(:)

ckino      complex(8),allocatable:: f0(:,:,:, :,:)

      real(8),allocatable:: ewgt(:,:)

      complex(8),allocatable:: smrho(:,:,:,:), smpot(:,:,:)


      logical:: cmdopt


C     call wkprnt(1)

      nproc  = mpipid(0)
      procid = mpipid(1)
#if MPI
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#endif

      if (nevec .le. 0) return
ckino      call tcn('rsibl_ev')

      nbas  = globalvariables%nbas
      nspec = globalvariables%nspec
      nlmto = ndimh-napw

C ... First setup

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      gmax=slat%gmax


      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ng=slat%ng

      vol=slat%vol


C ... Setup for q-dependent gv ... also makes kv, gv+q and iv
C     NB: gv generated by gvlst2 has q already added to it!
      call tcn('gvlist')
      call pshpr(iprint()-30)
C     call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,1000,0,ng,xx,xx,xx,xx)
      call gvlist(alat,plat,q,n1,n2,n3,gmax,500,0,ng,xx,xx,xx,xx)
ckino Dec.8.2011:         allocate(rv_p_ogv(abs(ng*3)))
       allocate(rv_a_ogv(abs(ng*3)))
ckino Dec.8.2011:         rv_p_ogv(:)=0.0d0
       rv_a_ogv(:)=0.0d0
       allocate(iv_p_okv(abs(ng*3)))
       iv_p_okv(:)=0


ckino dummy allocation
       allocate(smrho(k1,k2,k3,numq))
       allocate(smpot(k1,k2,3))

c
c      call defi (oiv, ng*3)
       allocate(w_oiv(ng*3))

        call gvlist ( alat , plat , q , n1 , n2 , n3 , gmax , 509 , ng 
ckino Dec.8.2011:       .  , ng , iv_p_okv , rv_p_ogv , w_oiv , w_oiv ) 
     .  , ng , iv_p_okv , rv_a_ogv , w_oiv , w_oiv ) 


      call poppr
      call tcx('gvlist')

C     For PW basis ... for now.
      if (napw .gt. 0) then
        allocate(ivp(napw))
        call gvgvcomp(ng,w_oiv,napw,igapw,ivp)
      else
        allocate(ivp(1))
      endif

C --- Tables of energies, rsm, indices to them ---
      call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)

C --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2

      allocate(w_ogq(ng*3), w_oyl(ng*nlmtop), w_oylw(ng*nlmtop), w_og2(ng), w_ohe(ng*net), w_ohr(ng*nrt))

Ckino H_L(G)= \frac{-4 pi}{e-G^2} {cal Y}_L(-iG) exp(gamma(e-G^2))
Ckino hsibl1 calculaets he=1/(e-G^2) and hr=exp(-gamma G^2)
Ckino the other parts are calculated in rsibl5.
      call dpzero(q0,3)
      if (nlmto .gt. 0) then
         call hsibl1 ( net , etab , nrt , rtab , ltop , alat , q0 , ng 
ckino Dec.8.2011:       .   , rv_p_ogv , w_ogq , w_og2 , w_oyl , w_ohe , 
     .   , rv_a_ogv , w_ogq , w_og2 , w_oyl , w_ohe , 
     .   w_ohr ) 

      endif
      deallocate(w_og2)

C ... Allocate local arrays; setup for fast phase factors
C     Setup already done in gvlst2
C     call defi(oiv,  ng*3)
C     call suphs0(plat,ng,w(ogv),w(oiv))
#if MPI
      nblk = 0
      nvsub = nevec/6
      if (nvsub .le. nproc) nblk = 1
      if (cmdopt('--noblk',7,0,strn)) nblk = 1
      allocate (vproc(0:nproc), stat=ierr)
      call dstrbp(nevec,nproc,nblk,vproc(0))
      if (mlog .and. iq .eq. 1) then
        call awrit3(' RSIBL: nevec=%i, np=%i, nblk=%i%N  proc   index',
     .  ' ',256,lgunit(1),nevec,nproc,nblk)
        write(lgunit(1),1) (i,vproc(i),i=0,nproc-1)
    1   format (3x,i3,3x,i5)
      endif
#elif OPENMP
      nblk = 16
      if (nproc .gt. 0) then
        nblk = max(min(16,nevec/nproc),4)
      endif
#else
      nblk = 16
#endif

C  --- Loop over eigenstates ---
C$OMP PARALLEL PRIVATE (psi,psir,vpsi,cosi,sini,wk,wk2,xsmrho,fr)
      allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
      allocate(psir(k1,k2,k3),cosi(ng),sini(ng),wk2(ng))
#if OPENMP
      allocate(xsmrho(k1,k2,k3,numq,nspc),fr(3,nbas,numq))
      call dpzero(xsmrho,k1*k2*k3*numq*nspc*2)
      call dpzero(fr,3*nbas*numq)
#endif

#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_RSIBL,procid,"rsibl")
#endif
      allocate(xsmrho(k1,k2,k3,numq,nspc))
      call dpzero(xsmrho,k1*k2*k3*numq*nspc*2)
ckino  Force is not necessary.
      lfrce=0 
ckino      if (lfrce .ne. 0) then
ckino        allocate(fr(3,nbas,numq))
ckino        call dpzero(fr,3*nbas*numq)
ckino      endif
c      do  ivec = vproc(procid), vproc(procid+1)-1, nblk
      ivecini= vproc(procid)
      ivecend= vproc(procid+1)-1
#else
C$OMP DO schedule(static,1) PRIVATE (ivec,nvec)
      ivecini= 1
      ivecend= nevec
#endif

        allocate(ewgt(numq,nevec))
        ewgt=1.0d0

ckino 
      write(*,'(a,9i5)') 'ivecinic,ivecend,nblk,nspc,nspec=',ivecini,ivecend,nblk,nspc,nspec

      ivecloop: do  ivec = ivecini,ivecend, nblk
#if MPI
        if (mlog .and. ivec .eq. vproc(procid)) then
          call gettime(datim)
          call awrit5(' rsibl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting evecs %i to %i nblk=%i',' ',256,lgunit(3),
     .    procid,nproc,vproc(procid),vproc(procid+1)-1,nblk)
        endif
        nvec = min(nblk, vproc(procid+1)-ivec)
        if (mlog)
     .  call awrit5(' procid=%i,vproc(procid)=%i,vproc(procid+1)-1=%i'//
     .  ',ivec=%i,nvec=%i',' ',256,lgunit(3),procid,vproc(procid),
     .  vproc(procid+1)-1,ivec,nvec)
#else
        nvec = min(nblk, nevec-ivec+1)
#endif


CSFCPP#if F90
C   ... Add together total smooth wavefunction
ckino   rsibl1 calculates fourier transformed smoothed Hankel, H(G), to psi
        call rsibl1(0,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2,
     .  n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,wk,
     .  wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,numq,
     .  ewgt,ivec,nvec,evec,w,psi,w)
Ckino    rsiblp adds PW(G) to psi 
        call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol),
     .  evec(1,1,ivec),psi)

Ckino now psi= H(G) + PW(G)

C   ... Add to real-space mesh, optionally make smpot*psi for forces
Ckino rsibl2 executes FFT to get psi(r), which is F0
Ckino and also calculates <psi|psi>(=F0F0) to get real space charge density.
Ckino psir is work area. 
ckino#if MPI | OPENMP
ckino         call rsibl2 ( ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2 
ckino     .   , k3 , iv_p_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
ckino     .   1 , 1 , 1  ) , psir , xsmrho , vpsi ) 
ckino
ckino#else
ckino         call rsibl2 ( ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2 
ckino     .   , k3 , iv_p_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
ckino     .   1 , 1 , 1  ) , psir , smrho ( 1 , 1 , 1 , 1  ) , vpsi ) 
ckino
ckino#endif

Ckino write psi=F0 part
         if (.true.) then
             call w_psir(ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2
     i   , k3 , iv_p_okv , isp, q ,iq 
     i   ,  n_eiglist,eiglist
     i   , plat,alat, nbas, slat%rv_p_opos, sspec(ssite(1:nbas)%spec)%z
     w  , psir )
Ckino  input      compelx(8):: psi(ng,nspc,nev)
Ckino  work area  complex(8):: psir(k1,k2,k3)
         endif



      enddo ivecloop

C$OMP END DO
      deallocate(psi,vpsi,wk,psir,cosi,sini,wk2)
#if OPENMP
CC kino's fix at feb27-2009
CC      deallocate(xsmrho,fr)
      deallocate(xsmrho)
      if (allocated(fr)) deallocate(fr)
#endif
C$OMP END PARALLEL

#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_RSIBL,procid,"rsibl")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
#endif

      allocate(w_osmbuf(k1*k2*k3*numq*nspc))
      call MPI_ALLREDUCE(xsmrho,w_osmbuf,2*k1*k2*k3*numq*nspc,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit2(' rsibl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce smrho',' ',256,lgunit(3),
     .  procid,nproc)
      endif
      call daxpy(2*k1*k2*k3*numq*nspc,1d0,w_osmbuf,1,
     .smrho(1,1,1,1,isp),1)
c      call rlse(osmbuf)
      deallocate(w_osmbuf)

      if (lfrce .ne. 0) then
        allocate(w_ofrbuf(3*nbas*numq))
        call MPI_ALLREDUCE(fr,w_ofrbuf,3*nbas*numq,
     .  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit2(' rsibl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' allreduce f',' ',256,lgunit(3),
     .    procid,nproc)
        endif
        call daxpy(3*nbas*numq,1d0,w_ofrbuf,1,f,1)
c        call rlse(ofrbuf)
        deallocate(w_ofrbuf)
      endif
      deallocate(xsmrho)
      if (allocated(fr)) deallocate(fr)
C
      deallocate(vproc, stat=ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
#endif
#endif
ckino add deallocations
       if (allocated(ewgt)) deallocate(ewgt)
ckino Dec.8.2011:         if (allocated(rv_p_ogv)) deallocate(rv_p_ogv)
       if (allocated(rv_a_ogv)) deallocate(rv_a_ogv)
       if (allocated(iv_p_okv)) deallocate(iv_p_okv) 

      deallocate(ivp)
C     call zprm3('smrho after rsibl',0,smrho,k1,k2,k3)
ckino      call tcx('rsibl_ev')
      end subroutine rsibl_ev
      end module m_rsibl_ev
