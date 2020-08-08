CSFCPP#define F90 1
      subroutine rsibl(ssite,sspec,slat,lfrce,nbas,isp,q,iq,ndimh,nspc,
     .napw,igapw,iprmb,numq,nevec,evec,ewgt,k1,k2,k3,smpot,smrho,f)
      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Add smooth part of output density into smrho and forces.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lfrce :if nonzero, accumulate contribution to force
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
      integer procid, master, nproc, mpipid
      integer lfrce,nbas,isp,k1,k2,k3,ndimh,nevec,numq,iprmb(1),iq,nspc
      integer napw,igapw(3,napw)
      real(8):: q(3) , ewgt(numq,nevec) , f(3,nbas,numq)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex evec(ndimh,nspc,nevec),smrho(k1,k2,k3,numq,isp),
     .smpot(k1,k2,k3,isp)
C ... Local parameters
      integer n0,nkap0,nermx,npmx,nblk,nlmto
      parameter (n0=10,nkap0=3,nermx=100,npmx=128)
      integer:: nspec , ngabc(3) , n1 , n2 , n3 , nrt , net , ng , 
     .  nglob , ltop , nlmtop , ogq , og2 , ohe , ohr , oyl , oylw , 
     .  oiv , iprint
ckino Dec.28.2011:  ckino 2011.05        integer,pointer :: iv_p_okv(:) =>NULL()
ckino 2011.05        integer,pointer :: iv_a_okv(:) =>NULL()
ckino Dec.28.2011:          integer,allocatable :: iv_p_okv(:)
      integer,allocatable :: iv_a_okv(:)

ckino 2011.05       real(8),pointer :: rv_p_ogv(:) =>NULL()
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
CSFCPP#if F90
      integer,allocatable:: ivp(:)
      complex(8),allocatable::psi(:,:,:),psir(:,:,:),vpsi(:,:,:),
     .wk(:,:,:)
      real(8),allocatable:: cosi(:),sini(:),wk2(:)
CSFCPP#else
CSFCPP      integer ivp(napw)
CSFCPP      integer ocos,of,opsi,osin,ovpsi,opsi0,owk2
CSFCPP#endif
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
      call tcn('rsibl')
      nbas  = globalvariables%nbas
      nspec = globalvariables%nspec
      nlmto = ndimh-napw
C ... First setup
      alat=slat%alat
      plat = slat%plat
      qlat = slat%qlat
      gmax=slat%gmax
      ngabc=slat%nabc
      ng=slat%ng
      vol=slat%vol
C ... Setup for q-dependent gv ... also makes kv, gv+q and iv
C     NB: gv generated by gvlst2 has q already added to it!
      call tcn('gvlst2')
      call pshpr(iprint()-30)
C     call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,1000,0,ng,xx,xx,xx,xx)
c      call gvlist(alat,plat,q,n1,n2,n3,gmax,500,0,ng,xx,xx,xx,xx)
      call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,500,0,ng,xx,xx,xx,xx)
      allocate(rv_a_ogv(abs(ng*3)))
      rv_a_ogv(:)=0.0d0
      allocate(iv_a_okv(abs(ng*3)))
      iv_a_okv(:)=0
      allocate(w_oiv(ng*3))
C      call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,1009,ng,ng,w(okv),
C     .  w(ogv),w(oiv),w(oiv))
c      call gvlist(alat , plat, q , n1 , n2 , n3 , gmax,509 , ng, ng, iv_a_okv, rv_a_ogv, w_oiv, w_oiv )
      call gvlst2(alat, plat, q, n1, n2, n3, 0d0,gmax,0,509, ng, ng, iv_a_okv, rv_a_ogv, w_oiv, w_oiv )
      call poppr
      call tcx('gvlst2')
C     For PW basis ... for now.
      if (napw .gt. 0) then
CSFCPP#if F90
        allocate(ivp(napw))
CSFCPP#endif
        call gvgvcomp(ng,w_oiv,napw,igapw,ivp)
      else
CSFCPP#if F90
        allocate(ivp(1))
CSFCPP#endif
      endif

C --- Tables of energies, rsm, indices to them ---
      call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)

C --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2
c      call defrr(ogq, ng*3)
c      call defrr(oyl, ng*nlmtop)
c      call defrr(oylw,ng*nlmtop)
c      call defrr(og2, ng)
c      call defrr(ohe, ng*net)
c      call defrr(ohr, ng*nrt)

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
CSFCPP#if F90
C$OMP PARALLEL PRIVATE (psi,psir,vpsi,cosi,sini,wk,wk2,xsmrho,fr)
      allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
      allocate(psir(k1,k2,k3),cosi(ng),sini(ng),wk2(ng))
#if OPENMP
      allocate(xsmrho(k1,k2,k3,numq,nspc),fr(3,nbas,numq))
      call dpzero(xsmrho,k1*k2*k3*numq*nspc*2)
      call dpzero(fr,3*nbas*numq)
#endif
CSFCPP#else
CSFCPP      call defcc(opsi,  ng*nspc*nblk)
CSFCPP      call defrr(ocos,  ng)
CSFCPP      call defrr(osin,  ng)
CSFCPP      call defcc(opsi0, ng*nspc*nblk)
CSFCPP      call defrr(owk2,  ng)
CSFCPP      call defcc(ovpsi, ng*nspc*nblk)
CSFCPP      call defcc(of,    k1*k2*k3)
CSFCPP#endif


#if MPI
      allocate(xsmrho(k1,k2,k3,numq,nspc))
      call dpzero(xsmrho,k1*k2*k3*numq*nspc*2)
      if (lfrce .ne. 0) then
        allocate(fr(3,nbas,numq))
        call dpzero(fr,3*nbas*numq)
      endif
c      do  ivec = vproc(procid), vproc(procid+1)-1, nblk
      ivecini= vproc(procid)
      ivecend= vproc(procid+1)-1
#else
C$OMP DO schedule(static,1) PRIVATE (ivec,nvec)
c      do  ivec = 1, nevec, nblk
c        nvec = min(nblk, nevec-ivec+1)
      ivecini= 1
      ivecend= nevec
#endif


      do  ivec = ivecini,ivecend, nblk
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
#if MPI | OPENMP
        call rsibl2 ( ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2
ckino Dec.28.2011:       .   , k3 , iv_p_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
     .   , k3 , iv_a_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
     .   1 , 1 , 1 , isp ) , psir , xsmrho , vpsi )

#else
        call rsibl2 ( ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2
ckino Dec.28.2011:       .   , k3 , iv_p_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
     .   , k3 , iv_a_okv , numq , ewgt ( 1 , ivec ) , lfrce , smpot ( 
     .   1 , 1 , 1 , isp ) , psir , smrho ( 1 , 1 , 1 , 1 , isp ) , vpsi 
     .   )

#endif

C    --- Add to forces ---
        if (lfrce .ne. 0) then
#if MPI | OPENMP
          call rsibl1(1,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2,
     .    n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,
     .    wk,wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,
     .    numq,ewgt,ivec,nvec,evec,vpsi,psi,fr)
#else
          call rsibl1(1,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2,
     .    n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,
     .    wk,wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,
     .    numq,ewgt,ivec,nvec,evec,vpsi,psi,f)
#endif
        endif

CSFCPP#else
CSFCPPC   ... Add together total smooth wavefunction
CSFCPP        call rsibl1(0,ssite,sspec,q,nbas,iprmb,ng,w(ogq),w(oiv),n1,n2,
CSFCPP     .  n3,qlat,w(ocos),w(osin),w(oyl),w_oylw,w_ohe,w_ohr,w(opsi0),
CSFCPP     .  w(owk2),vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,numq,
CSFCPP     .  ewgt,ivec,nvec,evec,w,w(opsi),w)
CSFCPP        call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol),
CSFCPP     .  evec(1,1,ivec),w(opsi))
CSFCPP
CSFCPPC   ... Add to real-space mesh, optionally make smpot*psi for forces
CSFCPP        call rsibl2(ng,nspc,nvec,w(opsi),n1,n2,n3,k1,k2,k3,w(okv),numq,
CSFCPP     .  ewgt(1,ivec),lfrce,smpot(1,1,1,isp),w(of),smrho(1,1,1,1,isp),
CSFCPP     .  w(ovpsi))
CSFCPP
CSFCPPC    --- Add to forces ---
CSFCPP        if (lfrce .ne. 0) then
CSFCPP          call rsibl1(1,ssite,sspec,q,nbas,iprmb,ng,w(ogq),w(oiv),n1,n2,
CSFCPP     .    n3,qlat,w(ocos),w(osin),w(oyl),w_oylw,w_ohe,w_ohr,
CSFCPP     .    w(opsi0),w(owk2),vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,
CSFCPP     .    numq,ewgt,ivec,nvec,evec,w(ovpsi),w(opsi),f)
CSFCPP        endif
CSFCPP#endif

      enddo
C$OMP END DO
CSFCPP#if F90
      deallocate(psi,vpsi,wk,psir,cosi,sini,wk2)
#if OPENMP
CC kino's fix at feb27-2009
CC      deallocate(xsmrho,fr)
      deallocate(xsmrho)
      if (allocated(fr)) deallocate(fr)
#endif
C$OMP END PARALLEL
CSFCPP#else
CSFCPP#endif

#if MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
c        call defrr(ofrbuf, 3*nbas*numq)
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
CC kino's fix at feb27-2009
CC     deallocate(xsmrho,fr)
      deallocate(xsmrho)
      if (allocated(fr)) deallocate(fr)
C
      deallocate(vproc, stat=ierr)
#endif
ckino 2011.05 if (associated(rv_p_ogv)) deallocate(rv_p_ogv)
ckino Dec.8.2011:         if (allocated(rv_p_ogv)) deallocate(rv_p_ogv)
      if (allocated(rv_a_ogv)) deallocate(rv_a_ogv)
ckino Dec.28.2011:         if (allocated(iv_p_okv)) deallocate(iv_p_okv)
      if (allocated(iv_a_okv)) deallocate(iv_a_okv)

CSFCPP#if F90
      deallocate(ivp)
CSFCPP#endif
C     call zprm3('smrho after rsibl',0,smrho,k1,k2,k3)
      call tcx('rsibl')
      end subroutine rsibl


      subroutine rsibl1(mode,ssite,sspec,q,nbas,iprmb,ng,gq,iv,n1,n2,n3,
     .qlat,cosgp,singp,yl,ylw,he,hr,psi0,wk2,vol,iprt,ipet,etab,rtab,
     .ndimh,nlmto,nspc,numq,ewgt,ivec,nvec,evec,vpsi,psi,f)

      use m_struc_def  !Cgetarg

C- Make wave function for a block of evecs, or add contr. to forces
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 make wave function
Ci         :1 Add 2*Re( (v psi+) grad(psi) ) to f
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   q     :Bloch wave number
Ci   nbas  :size of basis
Ci   ng    :number of G-vectors
Ci   gq    :2*pi/alat * (q+G) for all G-vectors
Ci   iv    :g-vectors as integer multiples of qlat (suphs0)
Ci   n1..3 :size uniform mesh for smooth density and potential
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   cosgp :cos(phase) for each g-vector
Ci   singp :sin(phase) for each g-vector
Ci   yl    :spherical harmonics for ng vectors
Ci   ylw   :work array of same dimension as yl
Ci   he    :table of energy factors
Ci   hr    :table of smoothing radius factors
Ci   psi0  :work array (dim ng*2*nspc*nev): psi sans phase factors
Ci   wk2   :work array of dimension ng
Ci   vol   :cell volume
Co   iprt  :index to which entry in rt a given orbital belongs
Ci   ipet  :index to which entry in etab a given orbital belongs
Ci   etab  :table of all inequivalent energies
Ci   rtab  :table of all inequivalent smoothing radii
Ci   ndimh :dimensions evec
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   numq  :number of trial fermi levels
Ci   ewgt  :weights for each of the trial fermi levels
Ci   ivec  :first of current block of eigenvectors
Ci   nvec  :number of eigenstates to generate
Ci   evec  :eigenvectors
Ci   vspi  :potential * wave function, needed only for mode=1
Co Outputs
Co   psi   :wave function (mode=0); work area (mode=1)
Co   f     :term added to forces (mode=1)
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      integer mode,nbas,ng,ndimh,nlmto,nspc,numq,ivec,nvec,iv(ng,3),
     .n1,n2,n3,n0,nkap0,iprmb(*)
      parameter (n0=10,nkap0=3)
      integer iprt(n0,nkap0,*),ipet(n0,nkap0,*)
      double precision vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*),
     .psi0(ng,2,nspc,nvec),wk2(ng),cosgp(ng),singp(ng),etab(*),
     .rtab(*),gq(ng,3),f(3,nbas,numq),ewgt(numq,nvec+ivec-1),qlat(3,3)
      double complex
     .psi(ng,nspc,nvec),evec(ndimh,nspc,ivec),vpsi(ng,nspc,nvec)
C ... Local parameters
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .blks(n0*nkap0),ntab(n0*nkap0),ncut(n0,nkap0),lh(nkap0),nkapi
      double precision e,rsm,eh(n0,nkap0),rsmh(n0,nkap0),f0(3)
      double precision xx(n0),wt,p(3)
      integer ib,is,io,jo,l2,kp,ie,ir,ioff,nlm1,nlm2,iq,kb,lt,i
c takao
      integer::ncutt
      call dpzero(psi, 2*ng*nspc*nvec)

      if (nlmto .eq. 0) return

      do  ib = 1, nbas

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,p,1)


        i_copy_size=size(sspec(is)%ngcut)
        call icopy(i_copy_size,sspec(is)%ngcut,1,ncut,1)

        call suphas(q,p,ng,iv,n1,n2,n3,qlat,cosgp,singp)
C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C       Block into groups with consecutive l and common (e,rsm)
        call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkapi)
        call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)

        call dpzero(psi0,ng*2*nspc*nvec)
        if (mode .eq. 1) call dpzero(psi, 2*ng*nspc*nvec)
        do  io = 1, norb
          if (blks(io) .ne. 0) then
            jo = ntab(io)
            l2 = ltab(io)
            lt = ltab(jo)
            kp = ktab(io)
            ie = ipet(l2+1,kp,is)
            ir = iprt(l2+1,kp,is)
            ioff = offl(io)
            nlm1 = l2**2+1
            nlm2 = nlm1 + blks(io)-1
            rsm = rtab(ir)
            e   = etab(ie)
c takao Apr2009
            ncutt=ncut(lt+1,kp)
            call ncutcorrect(ncutt,1,gq,ng)
c
Ckino    rsibl5 calculates 4 pi exp(e gamma) and {cal Y}_L(-iG)
Ckino    and make psi0(G), which is H(-iG), with hr and he.
            call rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,min(ng,ncutt)
     .      ,yl,ylw,he,hr,wk2,ioff,evec(1,1,ivec),ndimh,nspc,nvec,psi0)
          endif
        enddo
Ckino    multiply exp(i G R_i) * psi0 to make psi()
        call rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)

        if (mode .eq. 1) then
          do  i = 1, nvec
            call rsibl4(vol,ng,nspc,gq,vpsi(1,1,i),psi(1,1,i),f0)
            do  iq = 1, numq
              wt = ewgt(iq,i+ivec-1)
              f(1,ib,iq) = f(1,ib,iq) + wt*f0(1)
              f(2,ib,iq) = f(2,ib,iq) + wt*f0(2)
              f(3,ib,iq) = f(3,ib,iq) + wt*f0(3)
C             This shouldn't be necessary
              do  kb = 1, nbas
                f(1,kb,iq) = f(1,kb,iq) - wt*f0(1)/nbas
                f(2,kb,iq) = f(2,kb,iq) - wt*f0(2)/nbas
                f(3,kb,iq) = f(3,kb,iq) - wt*f0(3)/nbas
              enddo
            enddo
          enddo
        endif
      enddo
      end subroutine rsibl1


      subroutine rsibl2(ng,nspc,nev,psi,n1,n2,n3,k1,k2,k3,kv,numq,ewgt,
     .lfrce,smpot,f,smrho,vpsi)

C- FT wave function to real space and add square into mesh density
C  and optionally make smpot * psi
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nev   :number of wave functions
Ci   psi   :wave function in reciprocal space
Ci   n1..3 :size of FT mesh
Ci   k1..3 :dimensions smpot,smrho
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   numq  :number of trial fermi levels
Ci   ewgt  :weights for each of the trial fermi levels
Ci   lfrce :if nonzero, make vpsi
Ci   smpot :potential, needed if lfrce is nonzero
Co Outputs
Co   f     :psi in real space
Co   smrho :ewgt (f+)f added to smooth density
Co   vpsi  :FT (smpot * r.s. wave function) if lfrce is nonzero
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,ng,nspc,nev,numq,kv(ng,3),lfrce
      double precision ewgt(numq,nev)
      double complex psi(ng,nspc,nev),vpsi(ng,nspc,nev),f(k1,k2,k3)
      double complex smrho(k1,k2,k3,numq,nspc),smpot(k1,k2,k3,nspc)
C ... Local parameters
      integer i1,i2,i3,iq,i,ispc
      double precision wgt1

      call tcn('rsibl2')
      do  ispc = 1, nspc
        do  i = 1, nev
          call gvputf(ng,1,kv,k1,k2,k3,psi(1,ispc,i),f)
          call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)
          do  iq = 1, numq
            wgt1 = ewgt(iq,i)
            do  i3 = 1, n3
              do  i2 = 1, n2
                do  i1 = 1, n1
                  smrho(i1,i2,i3,iq,ispc) = smrho(i1,i2,i3,iq,ispc) +
     .            wgt1*dconjg(f(i1,i2,i3))*f(i1,i2,i3)
                enddo
              enddo
            enddo
          enddo

          if (lfrce .ne. 0) then
            do  i3 = 1, n3
              do  i2 = 1, n2
                do  i1 = 1, n1
                  f(i1,i2,i3) = f(i1,i2,i3)*smpot(i1,i2,i3,ispc)
                enddo
              enddo
            enddo
            call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
            call gvgetf(ng,1,kv,k1,k2,k3,f,vpsi(1,ispc,i))
          endif
        enddo
      enddo
      call tcx('rsibl2')
      end subroutine rsibl2


C      subroutine rsibl3(ng,n1,n2,n3,k1,k2,k3,kv,smpot,f,vpsi)
CC- Make f*smpot and transform to reciprocal space
C      implicit none
CC ... Passed parameters
C      integer k1,k2,k3,n1,n2,n3,ng,kv(ng,3)
C      double complex vpsi(ng),f(k1,k2,k3),smpot(k1,k2,k3)
CC ... Local parameters
C      integer i1,i2,i3
C
C      call tcn('rsibl3')
C      do  i3 = 1, n3
C        do  i2 = 1, n2
C          do  i1 = 1, n1
C            f(i1,i2,i3) = f(i1,i2,i3)*smpot(i1,i2,i3)
C          enddo
C        enddo
C      enddo
C      call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,kv,k1,k2,k3,f,vpsi)
C      call tcx('rsibl3')
C      end

      subroutine rsibl4(vol,ng,nspc,gq,vpsi,psi,f0)

C- Force term 2*Re( (psi_nu+) vsm grad(psi_nu) )
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   gq    :2*pi/alat * (q+G) for all G-vectors
Ci   vpsi  :(psi vsm) in reciprocal space
Ci   psi   :portion of wave function associated with one site ib
Co Outputs
Co   f0    :2*Re( (psi_nu+) vsm grad(psi_ib,nu) ) put in f0
Cr Remarks
Cr   gradient operator is i*G
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,nspc
      double precision vol,gq(ng,3),f0(3)
      double precision vpsi(2,nspc,ng),psi(2,nspc,ng)
C ... Local parameters
      integer ispc,ig
      double precision sum1,sum2,sum3,xx
C     double complex cc,ovl

      call tcn('rsibl4')
      sum1 = 0
      sum2 = 0
      sum3 = 0
C     ovl = 0
      do  ispc = 1, nspc
        do  ig = 1, ng
C        cc  = dcmplx(psi(1,ispc,ig),-psi(2,ispc,ig))
C     .       *dcmplx(vpsi(1,ispc,ig),vpsi(2,ispc,ig)) * vol
C        ovl = ovl + cc
C     .      * dcmplx(0d0,1d0)*gq(ig,1)
          xx = vpsi(2,ispc,ig)*psi(1,ispc,ig)
     .    - vpsi(1,ispc,ig)*psi(2,ispc,ig)
          sum1 = sum1 + gq(ig,1)*xx
          sum2 = sum2 + gq(ig,2)*xx
          sum3 = sum3 + gq(ig,3)*xx
        enddo
      enddo

      f0(1) = 2*vol*sum1
      f0(2) = 2*vol*sum2
      f0(3) = 2*vol*sum3

C     print *, ovl

      call tcx('rsibl4')

      end subroutine rsibl4


      subroutine rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,ncut,yl,ylw,he,hr,
     .wk,ioff,evec,ndimh,nspc,nvec,psi)

C- Add contribution to wave function from one block of orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci   ie    :index to appropriate entry in energy factor table he
Ci   ir    :index to appropriate entry in sm radius factor table hr
Ci   e     :hankel energy
Ci   rsm   :smoothing radius
Ci   vol   :cell volume
Ci   nlm1  :starting orbital L for which to accumulate wave function
Ci   nlm2  :final orbital L for which to accumulate wave function
Ci   ng    :number of G-vectors
Ci   ncut  :G-cutoff for wave function
Ci   yl    :spherical harmonics for ng vectors
Ci   ylw   :work array dimensioned same as yl
Ci   he    :table of energy factors
Ci   hr    :table of smoothing radius factors
Ci   wk    :work array of dimension at least ncut
Ci   ioff  :offset to hamiltonian (eigenvector) for this orbital block
Ci   ndimh :dimension of hamiltonian
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   evec  :eigenvectors
Ci   nvec  :number of eigenvectors
Co Outputs
Co   psi   :contribution to psi from this block accumulated
Cr Remarks
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ie,ioff,ir,ncut,ng,nlm1,nlm2,ndimh,nspc,nvec
      double precision e,rsm,vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*),
     .wk(ncut)
      double complex evec(ndimh,nspc,nvec),psi(ng,nspc,nvec)
C ... Local parameters
      integer i,ii,ilm,l,ll,lmax,m,iv,nlmx,ispc
      parameter (nlmx=100)
      double complex cfac(nlmx),cc,evp(nlmx,nspc,nvec)
      double precision pi,fac
      parameter (pi=3.1415926535897931d0)

      if (nlm2 .eq. 0) return
      call tcn('rsibl5')

C ... Phase and other factors
      lmax = ll(nlm2)
Ckino -4 pi exp(e gamma)
      fac = -4d0*pi*dexp(e*rsm*rsm*0.25d0)/vol
Ckino phase factor of {cal Y}_L(-iG)
      cc = (0d0,1d0)*fac
      ilm = 0
      do  l = 0, lmax
        cc = cc*(0d0,-1d0)
        do m = -l,l
          ilm = ilm+1
          cfac(ilm) = cc
        enddo
      enddo

C ... Combine G-dependent energy, rsm and YL factors
      do  i = 1, ncut
        wk(i) = he(i,ie)*hr(i,ir)
      enddo
      do  ilm = nlm1, nlm2
        do  i = 1, ncut
          ylw(i,ilm) = wk(i)*yl(i,ilm)
        enddo
      enddo

C ... Make vector evec*phase
      do  ispc = 1, nspc
        do  ilm = nlm1, nlm2
          ii = ilm-nlm1+ioff+1
          do  iv = 1, nvec
C          cc = evec(ii,ispc,iv)*cfac(ilm)
C          evpr(ilm,ispc,iv) = dble(cc)
C          evpi(ilm,ispc,iv) = dimag(cc)
            evp(ilm,ispc,iv) = evec(ii,ispc,iv)*cfac(ilm)
          enddo
        enddo

C ... For each orbital and evecs 1..nvec, accumulate psi
C     ii = ilm-nlm1+ioff+1
        do  ilm = nlm1, nlm2
          do  iv = 1, nvec
            do  i = 1, ncut
              psi(i,ispc,iv) = psi(i,ispc,iv)+ylw(i,ilm)*evp(ilm,ispc,iv)
            enddo
          enddo
        enddo
      enddo

      call tcx('rsibl5')
      end subroutine rsibl5


      subroutine rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)

C- Multiply by phase to make final FT of partial wave function
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nvec  :number of eigenvectors
Ci   cosgp :cos(phase)
Ci   singp :sin(phase)
Ci   wr    :real part of psi, unscaled by phase
Ci   wi    :imaginary part of psi, unscaled by phase
Co Outputs
Co   psi   :(wr,si)*(cosgp,singp) is added into psi
Cr Remarks
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,nspc,nvec
      double precision cosgp(ng),singp(ng)
      double complex psi0(ng,nspc,nvec),psi(ng,nspc,nvec)
C ... Local parameters
      integer i,iv,ispc

      call tcn('rsibl6')

Ckino phase factor exp(i G R_spc)
      do  iv = 1, nvec
        do  ispc = 1, nspc
          do  i = 1, ng
            psi(i,ispc,iv) = psi(i,ispc,iv)
     .      + psi0(i,ispc,iv)*dcmplx(cosgp(i),singp(i))
          enddo
        enddo
      enddo

      call tcx('rsibl6')
      end subroutine rsibl6


      subroutine rsibl7(xsmrho,k1,k2,k3,numq,smrho)

C- Combine rho from separate parallel threads
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,numq
      double complex smrho(k1,k2,k3,numq),xsmrho(k1,k2,k3,numq)
C ... Local parameters
      integer ik1,ik2,ik3,iq
      do  iq = 1, numq
        do  ik3 = 1, k3
          do  ik2 = 1, k2
            do  ik1 = 1, k1
              smrho(ik1,ik2,ik3,iq) = smrho(ik1,ik2,ik3,iq) +
     .        xsmrho(ik1,ik2,ik3,iq)
            enddo
          enddo
        enddo
      enddo

      end subroutine rsibl7


      subroutine gvgvcomp(ng,igv,napw,igapw,ivp)

C- Find pointer ivp that maps igv to igapw.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Co   igv   :list of reciprocal lattice vectors G (gvlist.f)
Ci   napw   :number of R.L.V for PW basis
Ci   igapw :reciprocal lattice vectors for PW basis.
Co Outputs
Co   ivp   :if ig = ivp(jg), igv(jg) and nvec(ig) are same vector
Cr Remarks
Cr  This routine should be be cleaned up and ivp
Cr  used by rest of program in place of igapw
Cu Updates
Cu   05 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,igv(ng,3),napw,igapw(3,napw),ivp(napw)
C ... Local parameters
      integer jg,ig

C     Redesign ... inefficient.
CSFCPP#if F90
      do  ig = 1, napw
        do  jg = 1, ng
          if(sum(abs(igv(jg,:)-igapw(:,ig))) == 0) then
            ivp(ig) = jg
            goto 333
          endif
        enddo
        ivp(ig) = -9999
        call rx('gvgvcomp wrong 111! maybe enlarge GMAX or so')
  333   continue
        if( sum(abs( igapw(:,ig)-igv(ivp(ig),:) )) /=0) 
     .  call rx('bug in gvgvcomp.  Cannot find ivp')
      enddo
CSFCPP#endif
      end subroutine gvgvcomp


      subroutine rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,sqv,evec,psi)

C- Plane wave part of evec
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nvec  :number of eigenstates to generate
Ci   evec  :eigenvectors
Ci   vspi  :potential * wave function, needed only for mode=1
Ci   sqv   :square root of volume
Co Outputs
Co   psi   :wave function
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,ndimh,nlmto,nspc,nvec
      integer napw,ivp(napw)
      double precision sqv
      double complex psi(ng,nspc,nvec),evec(ndimh,nspc,nvec)
C ... Local parameters
      integer i,ispc,igv

      if (napw .le. 0) return
C     sqvol = dsqrt()
      do  ispc = 1, nspc
        do  i = 1, nvec
          do  igv = 1, napw
            psi(ivp(igv),ispc,i) = psi(ivp(igv),ispc,i)
     .      + evec(nlmto+igv,ispc,i)/sqv
          enddo
        enddo
      enddo

      end subroutine rsiblp


