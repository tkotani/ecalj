Cdwdef1 have_common_w_in augmbl , #_of_w_access= 0
Cchp3 cmdvar= type(s_rv1),sv_p_,oppi
Cchp3 cmdvar= type(s_rv1),sv_p_,oppi
Cchp3 cmdvar= type(s_rv1),sv_p_,otau
Cchp3 cmdvar= type(s_rv1),sv_p_,osig
Cchp1 cmdvar= integer,iv_p_,ojcg
Cchp1 cmdvar= integer,iv_p_,oidxcg
Cchp1 cmdvar= real(8),rv_p_,ocy
Cchp1 cmdvar= real(8),rv_p_,ocg
c#define F90 1
Cchp3       subroutine augmbl(mode,ssite,sspec,slat,isp,lcplxp,q,
Cchp3      .osig,otau,oppi,iprmb,ndimh,
Cchp3      .napw,igapw,h,hso,s)
Cchp3        subroutine augmbl ( mode , ssite , sspec , slat , isp , lcplxp 
Cchp3      . , q , sv_p_osig , otau , oppi , iprmb , ndimh , napw , igapw 
Cchp3      . , h , hso , s ) 
Cchp3         subroutine augmbl ( mode , ssite , sspec , slat , isp , lcplxp 
Cchp3      .  , q , sv_p_osig , sv_p_otau , oppi , iprmb , ndimh , napw , igapw 
Cchp3      .  , h , hso , s ) 
         subroutine augmbl ( mode , ssite , sspec , slat , isp , lcplxp 
     .   , q , sv_p_osig , sv_p_otau , sv_p_oppi , iprmb , ndimh , napw 
     .   , igapw , h , hso , s ) 




      use m_struc_def  !Cgetarg
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
C- Adds augmentation part of H and S
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 compute both hamiltonian and overlap
Ci         :  otherwise, compute overlap only.
Ci         :  In this case, vavg is not used
Ci         :10s digit
Ci         :  0 do not compute hso
Ci         :  1 compute hso.
Ci         :    Note: only a portion of hso is computed for a
Ci         :    particular isp.  The total hso is assembled
Ci         :    after isp loops from 1..2.  hso should not be
Ci         :    initialized between isp=1 and isp=2 loops.
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa lmxb kmxt rsma
Ci     Stored:    *
Ci     Passed to: uspecb
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy
Ci     Stored:    *
Ci     Passed to: hxpbl
Ci   isp   :current spin channel
Ci   q     :Bloch wave number
Ci   osig  :overlap matrix of P_kL
Ci         :NB: also head-head, head-tail contributions; see augmat.f
Ci   otau  :kinetic energy matrix of P_kL
Ci         :NB: also head-head, head-tail contributions; see augmat.f
Ci         (otau is not needed because folded into ppi already)
Ci   oppi  :kinetic energy + potential matrix of P_kL
Ci         :NB: also head-head, head-tail contributions; see augmat.f
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndimh :dimension of h and s
Ci   napw  :number of PWs in APW part of basis
Ci   igapw :PWs in units of reciprocal lattice vectors
Co Outputs
Co   h     :augmentation part of hamiltonian matrix added to h
Co   hso   :spin off-diagonal block of spin-orbit hamiltonian
Co   s     :augmentation part of overlap matrix added to s
Cl Local variables
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cl         :at augmentation site ia, including local orbitals
Cl   nlmto :number of lmto basis functions
Cr Remarks
Cr   Some expressions labelled JMP refer to J.Math.Phys39, 3393 (1998)
Cb Bugs
Cb   Not really a bug, but an inefficiency:
Cb   Right now, strux are kept for all orbitals in the basis, including
Cb   expansions coffs for local orbitals (which are set to zero).
Cb   Better to condense strux to reduce computational effort for 2-
Cb   and 3-center terms.
Cm MPI
Cm   See remarks in hsibl. Buffers for h and s are taken from the heap.
Cm   In addition a buffer the same size as as h and s for ALLREDUCE.
Cu Updates
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   08 Sep 06 (WRL) updated MPI to work with SO coupling
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   03 Feb 05 (A. Chantis) calculate hso
Cu    1 Sep 04 Adapted to handle complex ppi.  S.O. folded into ppi
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   29 Jun 04 (A. Chantis) Include LzSz spin-orbit coupling
Cu   14 Aug 02 Added overlap-only option
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   11 Jan 02 Adapted for f90 constructs
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 Jun 00 spin polarized
Cu   18 May 00 Adapted from nfp augm_q.f
Cu   1998      (DLN) parallel version for SGI
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI
#if MPE
      include "mpef.h"
#endif
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog,cmdopt
      integer i,lgunit
      character*120 strn
#endif

Cchp3       integer mode,lcplxp,isp,ndimh,napw,igapw(3,napw),iprmb(ndimh),
Cchp3      .osig(3,1),otau(3,1),oppi(3,1)
Cchp3        integer:: mode , lcplxp , isp , ndimh , napw , igapw(3,napw) 
Cchp3      . , iprmb(ndimh) , otau(3,1) , oppi(3,1) 
Cchp3         integer:: mode , lcplxp , isp , ndimh , napw , igapw(3,napw) 
Cchp3      .  , iprmb(ndimh) , oppi(3,1) 
         integer:: mode , lcplxp , isp , ndimh , napw , igapw(3,napw) 
     .   , iprmb(ndimh) 
         type(s_rv1) :: sv_p_oppi(3,1)

        type(s_rv1) :: sv_p_otau(3,1)

       type(s_rv1) :: sv_p_osig(3,1)

      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex h(ndimh,ndimh),s(ndimh,ndimh),hso(ndimh,ndimh)
C ... Local parameters
      integer nlmbx,nlmax,ktop0,lofb,mode0
      parameter (ktop0=20, nlmbx=49, nlmax=49, lofb=(ktop0+1)*nlmax)
c#if F90
      complex(8),allocatable:: b(:)
c#else
c      double complex b(lofb*ndimh),b0(lofb*nlmbx),bos(lofb)
c#endif
      double complex g(lofb)
Cchp1       integer ia,isa,kmax,lmxa,lmxha,nbas,nglob,nlma,
Cchp1      .nlmha,nkaph,ocg,ocy,oidxcg,ojcg,mode1,nlmto
Cchp1        integer:: ia , isa , kmax , lmxa , lmxha , nbas , nglob , nlma 
Cchp1      . , nlmha , nkaph , ocy , oidxcg , ojcg , mode1 , nlmto 
Cchp1         integer:: ia , isa , kmax , lmxa , lmxha , nbas , nglob , nlma 
Cchp1      .  , nlmha , nkaph , oidxcg , ojcg , mode1 , nlmto 
Cchp1          integer:: ia , isa , kmax , lmxa , lmxha , nbas , nglob , nlma 
Cchp1      .   , nlmha , nkaph , ojcg , mode1 , nlmto 
          integer:: ia , isa , kmax , lmxa , lmxha , nbas , nglob , nlma 
     .    , nlmha , nkaph , mode1 , nlmto 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)

      double precision rsma,pa(3),xx,alat,qlat(3,3),vol
c#ifndef F90
c      integer ib,is,ik,nlmh
c      double precision p(3)
cC ... For the basis
c      integer nkap0,n0,nkapi
c      parameter (nkap0=3,n0=10)
c      integer lh(nkap0)
c      double precision eh(n0,nkap0),rsmh(n0,nkap0)
c      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
c#endif
#if MPI
      integer, dimension(:), allocatable :: bproc
      complex(8) ,allocatable :: h_zv(:)
      complex(8) ,allocatable :: hbuf_zv(:)
      complex(8) ,allocatable :: s_zv(:)
      complex(8) ,allocatable :: sbuf_zv(:)
      complex(8) ,allocatable :: hso_zv(:)

#endif
      integer:: initbas, endbas


C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w

      call tcn ('augmbl')

#if MPI
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#else
      numprocs=1
#endif

C --- Setup ---
      nbas  = nglob('nbas')
      nkaph = nglob('nkaph')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      nlmto = ndimh-napw

Cchp1       ocg=slat%ocg
       rv_p_ocg => slat%rv_p_ocg 

Cchp1       ojcg=slat%ojcg
       iv_p_ojcg => slat%iv_p_ojcg 

Cchp1       oidxcg=slat%oidxcg
       iv_p_oidxcg => slat%iv_p_oidxcg 

Cchp1       ocy=slat%ocy
       rv_p_ocy => slat%rv_p_ocy 



      alat=slat%alat
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      vol=slat%vol

C     tpiba = 2d0*4d0*datan(1d0)/alat

c#if F90
      allocate (b(lofb*ndimh))
c#endif

C --- Loop over augmentation sites --- ---
#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_AUGMBL,procid,"augmbl")
#endif
      allocate(h_zv(abs(-ndimh*ndimh)))
      if (-ndimh*ndimh<0) h_zv(:)=0.0d0

      allocate(s_zv(abs(-ndimh*ndimh)))
      if (-ndimh*ndimh<0) s_zv(:)=0.0d0

      if (lcplxp .ne. 0) then
        allocate(hso_zv(abs(-ndimh*ndimh)))
        if (-ndimh*ndimh<0) hso_zv(:)=0.0d0

      endif
      allocate (bproc(0:numprocs), stat=ierr)
      call dstrbp(nbas,numprocs,1,bproc(0))
      initbas=bproc(procid)
      endbas =bproc(procid+1)-1
#else
      initbas=1
      endbas=nbas
#endif

      do  ia = initbas, endbas
#if MPI
        if (mlog .and. ia .eq. bproc(procid)) then
          call gettime(datim)
          call awrit4(' augmbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting atoms %i to %i',' ',256,lgunit(3),
     .    procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
#endif


        isa=ssite(ia)%spec
        i_copy_size=size(ssite(ia)%pos)
        call dcopy(i_copy_size,ssite(ia)%pos,1,pa,1)


        lmxa=sspec(isa)%lmxa
        lmxha=sspec(isa)%lmxb
        kmax=sspec(isa)%kmxt
        rsma=sspec(isa)%rsma

        nlmha = (lmxha+1)**2
        nlma  = (lmxa+1)**2
        if (lmxa .eq. -1) goto 10

C   --- Make strux to expand all orbitals at site ia ---
        call rxx((kmax+1)*nlma .gt. lofb, 'augmbl: increase lofb')
c#if F90
C Check ...
C        allocate (b(0:kmax,nlma,ndimh))
Cchp1         call bstrux(0,slat,ssite,sspec,w(ocg),w(oidxcg),w(ojcg),w(ocy),
Cchp1      .  iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,
Cchp1      .  napw,igapw,b,xx)
Cchp1          call bstrux ( 0 , slat , ssite , sspec , rv_p_ocg , w ( oidxcg 
Cchp1      .   ) , w ( ojcg ) , w ( ocy ) , iprmb , nbas , ia , pa , rsma , 
Cchp1      .   q , kmax , nlma , ndimh , napw , igapw , b , xx ) 
Cchp1           call bstrux ( 0 , slat , ssite , sspec , rv_p_ocg , w ( oidxcg 
Cchp1      .    ) , w ( ojcg ) , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q 
Cchp1      .    , kmax , nlma , ndimh , napw , igapw , b , xx ) 
Cchp1            call bstrux ( 0 , slat , ssite , sspec , rv_p_ocg , iv_p_oidxcg 
Cchp1      .     , w ( ojcg ) , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
Cchp1      .     kmax , nlma , ndimh , napw , igapw , b , xx ) 
            call bstrux ( 0 , slat , ssite , sspec , rv_p_ocg , iv_p_oidxcg 
     .      , iv_p_ojcg , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
     .      kmax , nlma , ndimh , napw , igapw , b , xx ) 




C$$$#else
C$$$C       Zero out strux to eliminate contributions from local orbitals
C$$$        call dpzero(b,(kmax+1)*nlma*ndimh*2)
C$$$        do  ib = 1, nbas
C$$$Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
C$$$Cgetarg...info...           real(8) :: p(3) , real(8) :: ssite%pos 3
C$$$Cgetarg           call upack('site spec pos',ssite,ib,is,p,0,0)

C$$$           is=ssite(ib)%spec
C$$$           i_copy_size=size(ssite(ib)%pos)
C$$$           call dcopy(i_copy_size,ssite(ib)%pos,1,p,1)

C$$$C         This line augments no local orbitals
C$$$C         ik = 1
C$$$C         This line augments onsite extended local orbitals only
C$$$C         if (ia .eq. ib) ik = 2
C$$$C         This line augments extended local orbitals all sites
C$$$          ik = 2
C$$$          call uspecb(0,ik,sspec,is,is,lh,rsmh,eh,nkapi)
C$$$C         Position in h; l,k indices for orbitals connected w/ ib
C$$$          call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C$$$C         Loop over blocks of envelope functions
C$$$          do  ik = 1, nkapi
C$$$            nlmh = (lh(ik)+1)**2
C$$$            if (nlmh .gt. nlmbx) call rxi('augmbl: need nlmbx',nlmh)
C$$$            if (nlmh .gt. nlma .and. ia.eq.ib)
C$$$     .        call rx('augmbl: nlmh.gt.nlma')
C$$$            call hxpbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,nlma,
C$$$     .        kmax,nlma,w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,b0)
C$$$C       ... Subtract strux of on-site heads if ib=ia
C$$$#ifndef ALL3C
C$$$            if (ib .eq. ia) then
C$$$              call hxpos(rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,kmax,bos)
C$$$              call paugq2(kmax,nlmh,nlma,bos,b0)
C$$$            endif
C$$$#endif
C$$$            call paugq1(kmax,nlma,kmax,ik,norb,ltab,ktab,rsmh,offl,b0,b)
C$$$          enddo
C$$$        enddo
C$$$#endif

C   --- Add 1-center and 2-center terms ---
#if MPI
        if (lcplxp .eq. 0) then
Cchp3           call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
Cchp3      .    , kmax , nlma , w ( osig ( 3 , ia ) ) , w ( oppi ( 3 , ia ) ) 
Cchp3      .    , w ( osig ( 2 , ia ) ) , w ( oppi ( 2 , ia ) ) , b , ndimh , 
Cchp3      .    nlmto , s_zv , h_zv )
Cchp3            call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
Cchp3      .     , kmax , nlma , sv_p_osig( 3 , ia )%v , w ( oppi ( 3 , ia ) ) 
Cchp3      .     , sv_p_osig( 2 , ia )%v , w ( oppi ( 2 , ia ) ) , b , ndimh , 
Cchp3      .     nlmto , s_zv , h_zv ) 
            call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
     .      , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s_zv , h_zv ) 



        else
Cchp3           call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
Cchp3      .    , nlmha , kmax , nlma , w ( osig ( 3 , ia ) ) , w ( oppi ( 3 
Cchp3      .    , ia ) ) , w ( osig ( 2 , ia ) ) , w ( oppi ( 2 , ia ) ) , b 
Cchp3      .    , ndimh , nlmto , s_zv , h_zv , hso_zv )
Cchp3            call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
Cchp3      .     , nlmha , kmax , nlma , sv_p_osig( 3 , ia )%v , w ( oppi ( 3 
Cchp3      .     , ia ) ) , sv_p_osig( 2 , ia )%v , w ( oppi ( 2 , ia ) ) , b 
Cchp3      .     , ndimh , nlmto , s_zv , h_zv , hso_zv ) 
            call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
     .      , nlmha , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s_zv , h_zv , hso_zv ) 



        endif
#else
#ifndef ALL3C
        if (lcplxp .eq. 0) then
Cchp3           call augq12(mode0,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,nlma,
Cchp3      .    w(osig(3,ia)),w(oppi(3,ia)),w(osig(2,ia)),
Cchp3      .    w(oppi(2,ia)),b,ndimh,nlmto,s,h)
Cchp3            call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
Cchp3      .     , kmax , nlma , sv_p_osig( 3 , ia )%v , w ( oppi ( 3 , ia ) ) 
Cchp3      .     , sv_p_osig( 2 , ia )%v , w ( oppi ( 2 , ia ) ) , b , ndimh , 
Cchp3      .     nlmto , s , h ) 
            call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
     .      , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s , h ) 


        else
Cchp3           call augq2z(mode0,mode1,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,
Cchp3      .    nlma,w(osig(3,ia)),w(oppi(3,ia)),w(osig(2,ia)),
Cchp3      .    w(oppi(2,ia)),b,ndimh,nlmto,s,h,hso)
Cchp3            call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
Cchp3      .     , nlmha , kmax , nlma , sv_p_osig( 3 , ia )%v , w ( oppi ( 3 
Cchp3      .     , ia ) ) , sv_p_osig( 2 , ia )%v , w ( oppi ( 2 , ia ) ) , b 
Cchp3      .     , ndimh , nlmto , s , h , hso ) 
            call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
     .      , nlmha , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s , h , hso ) 


        endif
#endif
#endif

C   --- Add B+ sig B to S and B+ ppi B to H ---
#if MPI
Cchp3         call augqs3 ( kmax , lmxa , nlma , ndimh , isp , g , w ( osig 
Cchp3      .  ( 1 , ia ) ) , b , s_zv )
         call augqs3 ( kmax , lmxa , nlma , ndimh , isp , g , sv_p_osig( 1 , ia )%v 
     .   , b , s_zv ) 


        if (mode0 .eq. 0 .and. lcplxp .eq. 0) then
Cchp3           call augqp3 ( kmax , nlma , ndimh , isp , g , w ( oppi ( 1 , 
Cchp3      .    ia ) ) , b , h_zv )
           call augqp3 ( kmax , nlma , ndimh , isp , g , sv_p_oppi( 1 , ia )%v 
     .     , b , h_zv ) 


        elseif (mode0 .eq. 0 .and. lcplxp .ne. 0) then
Cchp3           call augq3z ( mode1 , kmax , nlma , ndimh , isp , g , w ( oppi 
Cchp3      .    ( 1 , ia ) ) , b , h_zv , hso_zv )
           call augq3z ( mode1 , kmax , nlma , ndimh , isp , g , sv_p_oppi( 1 , ia )%v 
     .     , b , h_zv , hso_zv ) 


        endif
#else
Cchp3         call augqs3(kmax,lmxa,nlma,ndimh,isp,g,w(osig(1,ia)),b,s)
         call augqs3 ( kmax , lmxa , nlma , ndimh , isp , g , sv_p_osig( 1 , ia )%v 
     .   , b , s ) 

        if (mode0 .eq. 0 .and. lcplxp .eq. 0) then
Cchp3           call augqp3(kmax,nlma,ndimh,isp,g,w(oppi(1,ia)),b,h)
           call augqp3 ( kmax , nlma , ndimh , isp , g , sv_p_oppi( 1 , ia )%v 
     .     , b , h ) 

        elseif (mode0 .eq. 0 .and. lcplxp .ne. 0) then
Cchp3           call augq3z(mode1,kmax,nlma,ndimh,isp,g,w(oppi(1,ia)),b,h,
Cchp3      .    hso)
           call augq3z ( mode1 , kmax , nlma , ndimh , isp , g , sv_p_oppi( 1 , ia )%v 
     .     , b , h , hso ) 

        endif
#endif

   10   continue
C ... end loop over ia
      enddo
#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_AUGMBL,procid,"augmbl")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
#endif
      allocate(hbuf_zv(ndimh*ndimh))
Ctakao_ZeroClear_NotRequiered        if (ndimh*ndimh<0) hbuf_zv(:)=0.0d0

      call mpi_allreduce ( h_zv , hbuf_zv , 2 * ndimh * ndimh 
     ., mpi_double_precision , mpi_sum , mpi_comm_world , ierr )

      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce h ndimh=%i',' ',256,lgunit(3),
     .  procid,numprocs,ndimh)
      endif
      call daxpy ( 2 * ndimh * ndimh , 1d0 , hbuf_zv , 1 , h , 1 
     .)

      if (allocated(hbuf_zv)) deallocate(hbuf_zv)

      allocate(sbuf_zv(ndimh*ndimh))
Ctakao_ZeroClear_NotRequiered        if (ndimh*ndimh<0) sbuf_zv(:)=0.0d0

      call mpi_allreduce ( s_zv , sbuf_zv , 2 * ndimh * ndimh 
     ., mpi_double_precision , mpi_sum , mpi_comm_world , ierr )

      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce s ndimh=%i',' ',256,lgunit(3),
     .  procid,numprocs,ndimh)
      endif
      call daxpy ( 2 * ndimh * ndimh , 1d0 , sbuf_zv , 1 , s , 1 
     .)

      if (allocated(sbuf_zv)) deallocate(sbuf_zv)

      if (lcplxp .ne. 0) then
        allocate(sbuf_zv(ndimh*ndimh))
Ctakao_ZeroClear_NotRequiered         if (ndimh*ndimh<0) sbuf_zv(:)=0.0d0

        call mpi_allreduce ( hso_zv , sbuf_zv , 2 * ndimh * ndimh 
     .  , mpi_double_precision , mpi_sum , mpi_comm_world , ierr )

        if (mlog) then
          call gettime(datim)
          call awrit3(' augmbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' allreduce hso ndimh=%i',' ',256,lgunit(3),
     .    procid,numprocs,ndimh)
        endif
        call daxpy ( 2 * ndimh * ndimh , 1d0 , sbuf_zv , 1 , hso , 
     .  1 )

        if (allocated(sbuf_zv)) deallocate(sbuf_zv)

        if (allocated(hso_zv)) deallocate(hso_zv)

      endif
      if (allocated(s_zv)) deallocate(s_zv)
      if (allocated(h_zv)) deallocate(h_zv)


Cki#error ERROR, try to release name= os ,but list=null at linenumber= 400 list= (None)

      deallocate(bproc, stat=ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
#endif
#endif

C      call z2herm('U',ndimh,ndimh,h)
C      call z2herm('U',ndimh,ndimh,s)
C      call zprm('h-aug',2,h,ndimh,ndimh,ndimh)
C      call zprm('s-aug',2,s,ndimh,ndimh,ndimh)

c#if F90
      deallocate (b)
c#endif

      call tcx ('augmbl')

      end subroutine augmbl


      subroutine augq12(mode,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,
     .nlma,sighh,ppihh,sighp,ppihp,b,ndimh,nlmto,s,h)

C- Add one and two-center terms to hamiltonian and overlap matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 compute both hamiltonian and overlap
Ci         :  otherwise, compute overlap only.
Ci         :  In this case, vavg is not used
Ci   ia    :augmentation site about which strux are expanded
Ci   isp   :current spin channel
Ci   nkaph :dimensions augmentation matrices
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions augmentation potential matrix at site a
Ci   lmxha :dimensions sighh at site a
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   sighh :augmentation head-head overlap matrix
Ci   ppihh :augmentation head-head potential matrix
Ci   sighp :augmentation head-Pkl overlap matrix
Ci   ppihp :augmentation head-Pkl potential matrix
Ci   b     :Bloch strux connecting site ia to all sites
Ci   ndimh :hamiltonian dimension
Ci   nlmto :Number of lmto basis functions
Co Outputs
Co   h     :1- and 2- center augmentation part of ham. added to h
Co   s     :1- and 2- center augmentation part of ovlp added to s
Cr Remarks
Cr  In this implementation, the augmentation matrices and the row
Cr  dimension of the structure constants b follow normal L order.
Cr  The column dimension of b is permuted in iprmb order.
Cu Updates
Cu   01 Sep 04 folded so into complex potential
Cu   29 Jun 04 (A. Chantis) added 1- and 2- center spherical so*Lz*Sz
Cu   14 Aug 02 Added overlap-only option
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ia,isp,kmax,nkaph,ndimh,nlma,lmxha,nlmha,iprmb(1)
      integer nlmto
      double precision
     .sighh(nkaph,nkaph,0:lmxha,1), ppihh(nkaph,nkaph,nlmha,nlmha,1),
     .sighp(nkaph,0:kmax,0:lmxha,1),ppihp(nkaph,0:kmax,nlmha,nlma,1)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),h(ndimh,ndimh)
C ... Local parameters
      integer iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      double precision xx
      double complex cadd

C     call zprm('strux',2,b,(kmax+1)*nlma,(kmax+1)*nlma,ndimh)

C --- Loop over basis functions at site ia (augentation index) ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
      do  iorb = 1, norb
C       l1,ik1 = l and kaph indices, needed for sigma
        l1  = ltab(iorb)
        ik1 = ktab(iorb)
C       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
        i1 = offl(iorb)
        do  ilm1 = l1**2+1, (l1+1)**2
          i1 = i1+1

C     ... Two-center terms
C         Loop over basis functions 1..ndimh from all sites
          if (mode .eq. 0) then
            do  j = 1, ndimh
              do  k = 0, kmax
                cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
                s(i1,j) = s(i1,j) + cadd
                s(j,i1) = s(j,i1) + dconjg(cadd)
                do  ilma = 1, nlma
                  cadd = ppihp(ik1,k,ilm1,ilma,isp)*b(k,ilma,j)
                  h(i1,j) = h(i1,j) + cadd
                  h(j,i1) = h(j,i1) + dconjg(cadd)
                enddo
              enddo
            enddo

C     ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
                if (ilm1 .eq. ilm2)
     .          s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
              enddo
            enddo
          else

            do  j = 1, ndimh
              do  k = 0, kmax
                cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
                s(i1,j) = s(i1,j) + cadd
                s(j,i1) = s(j,i1) + dconjg(cadd)
              enddo
            enddo

C     ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                if (ilm1 .eq. ilm2)
     .          s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
              enddo
            enddo
          endif

        enddo
      enddo

      end subroutine augq12


      subroutine augq2z(mode,mode1,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,
     .nlma,sighh,ppihh,sighp,ppihp,b,ndimh,nlmto,s,h,hso)

C- Add one and two-center terms to h,s for complex potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 compute both hamiltonian and overlap
Ci         :  otherwise, compute overlap only.
Ci         :  In this case, vavg is not used
Ci   mode1 :0 do not compute hso
Ci         :1 compute hso
Ci   ia    :augmentation site about which strux are expanded
Ci   isp   :current spin channel
Ci   nkaph :dimensions augmentation matrices
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions augmentation potential matrix at site a
Ci   lmxha :dimensions sighh at site a
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   sighh :augmentation head-head overlap matrix
Ci   ppihh :augmentation head-head potential matrix
Ci   sighp :augmentation head-Pkl overlap matrix
Ci   ppihp :augmentation head-Pkl potential matrix
Ci   b     :Bloch strux connecting site ia to all sites
Ci   ndimh :hamiltonian dimension
Co Outputs
Co   h     :1- and 2- center augmentation part of ham. added to h
Co   s     :1- and 2- center augmentation part of ovlp added to s
Co   hso   :1- and 2- center spin up-down spin orbit block
Cr Remarks
Cr  In this implementation, the augmentation matrices and the row
Cr  dimension of the structure constants b follow normal L order.
Cr  The column dimension of b is permuted in iprmb order.
Cr  The ppihh(i,i,i,i,3), ppihh(i,i,i,i,4) are the head-head matrix
Cr  elements of LxSx+LySy. The ppihp(i,i,i,i,3), ppihp(i,i,i,i,4) are
Cr  the corresponding head-tail elements.
Cr  The 2c term has the form h_{i,j} = Sum_kL(conjg(b_{i;k,L})*p_{j;k,L})+
Cr   Sum_kL(p_{i;k,L}*p_{j;k,L}); To get the second term for spin orbit
Cr   we rely on the hermicity of the ppi_{LxSx+LySy} block.
Cr   Symbolically:
Cr   hso_{i,j,u,d} =  Sum_kL[p_{i,j,u,d}*b_{j} + conjg(p_{j,i,d,u}*b_{i})]
Cr   where u = spin-up and d = spin-down.
Cr   If the structure constants become noncollinear, additional terms have
Cr   to be added in the matrix element above.
Cu Updates
Cu   03 Feb 05 (A. Chantis) added 1- and 2- center spherical so*(LxSx+LySy)
Cu   01 Sep 04 folded so into complex potential
Cu   29 Jun 04 (A. Chantis) added 1- and 2- center spherical so*Lz*Sz
Cu   14 Aug 02 Added overlap-only option
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,mode1,ia,isp,kmax,nkaph,nlma,lmxha,nlmha,iprmb(1),
     .ndimh,nlmto
      double precision
     .sighh(nkaph,nkaph,0:lmxha,1),sighp(nkaph,0:kmax,0:lmxha,1)
      double complex
     .ppihh(nkaph,nkaph,nlmha,nlmha,isp+2*mode1),
     .ppihp(nkaph,0:kmax,nlmha,nlma,isp+2*mode1)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),
     .h(ndimh,ndimh),hso(ndimh,ndimh)
C ... Local parameters
      integer iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2,
     .jsp,ksp
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      double precision xx
      double complex cadd,cadd1
c#if F90
      complex(8),allocatable:: tso(:,:,:,:)
c#else
c      double complex tso(ndimh,ndimh,2,2)
c#endif

      if (mode1 .eq. 1) then
c#if F90
        allocate (tso(ndimh,ndimh,2,2))
c#endif
        call dpzero(tso,ndimh*ndimh*4*2)
      endif

C     call zprm('strux',2,b,(kmax+1)*nlma,(kmax+1)*nlma,ndimh)
C --- Loop over basis functions at site ia (augentation index) ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
      do  iorb = 1, norb
C       l1,ik1 = l and kaph indices, needed for sigma
        l1  = ltab(iorb)
        ik1 = ktab(iorb)
C       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
        i1 = offl(iorb)
        do  ilm1 = l1**2+1, (l1+1)**2
          i1 = i1+1

C     ... Two-center terms
C         Loop over basis functions 1..ndimh from all sites
          if (mode .eq. 0) then
            do  j = 1, ndimh
              do  k = 0, kmax
                cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
                s(i1,j) = s(i1,j) + cadd
                s(j,i1) = s(j,i1) + dconjg(cadd)
                do  ilma = 1, nlma
                  cadd = ppihp(ik1,k,ilm1,ilma,isp)*b(k,ilma,j)
C           ... Make 2c term LxSx+LySy part of SO
                  if (mode1 .eq. 1) then
                    do  jsp = 1, 2
                      if (isp .ne. jsp) then
                        ksp = 2*isp + jsp - 1
                        cadd1 = ppihp(ik1,k,ilm1,ilma,ksp)*b(k,ilma,j)
                        tso(i1,j,isp,jsp) = tso(i1,j,isp,jsp)
     .                  + 0.5d0*cadd1
                        tso(j,i1,jsp,isp) = tso(j,i1,jsp,isp)
     .                  + 0.5d0*dconjg(cadd1)
                      endif
                    enddo
                  endif
                  h(i1,j) = h(i1,j) + cadd
                  h(j,i1) = h(j,i1) + dconjg(cadd)
                enddo
              enddo
            enddo

C     ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
C        ...  Make 1c LxSx+LySy part of SO
                if (mode1 .eq. 1 .and. isp .eq. 2) hso(i1,i2) = hso(i1,i2)
     .          + 0.5d0*ppihh(ik1,ik2,ilm1,ilm2,3)
                if (ilm1 .eq. ilm2) s(i1,i2) = s(i1,i2)
     .          + sighh(ik1,ik2,l1,isp)
              enddo
            enddo
          else

            do  j = 1, ndimh
              do  k = 0, kmax
                cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
                s(i1,j) = s(i1,j) + cadd
                s(j,i1) = s(j,i1) + dconjg(cadd)
              enddo
            enddo

C     ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                if (ilm1 .eq. ilm2)
     .          s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
              enddo
            enddo
          endif

        enddo
      enddo

      if (mode1 .eq. 1) then
        call dpadd(hso(1,1),tso(1,1,1,2),1,2*ndimh*ndimh,1d0)
c#if F90
        deallocate (tso)
c#endif
      endif

      end subroutine augq2z


      subroutine augqs3(kmax,lmxa,nlma,ndimh,isp,g,sig,b,s)

C- Add B+ sig B to s for L-diagonal sig
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff
Ci   lmxa  :dimensions sig at site a
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   sig   :augmentation Pkl-Pkl overlap matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   s     :overlap matrix
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,lmxa,nlma,ndimh,isp
      double precision sig(0:kmax,0:kmax,0:lmxa,isp)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),
     .g(0:kmax,nlma),csum
C ... Local parameters
      integer nlmax
      parameter (nlmax=49)
      integer kjtop,i1,i2,ilm,k1,k2,l,kjlm,lla(nlmax)
      data lla/0,3*1,5*2,7*3,9*4,11*5,13*6/

      if (nlma .gt. nlmax) call rxi('augqs3: increase nlmax to',nlma)

C     call tcn('augqs3')
      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh

C   ... Make sig*b in g
        do  ilm = 1, nlma
          l = lla(ilm)
          do  k1 = 0, kmax
            g(k1,ilm) = 0d0
            do  k2 = 0, kmax
              g(k1,ilm) = g(k1,ilm) + sig(k1,k2,l,isp)*b(k2,ilm,i2)
            enddo
          enddo
        enddo
C   ... Make dot products with vectors i1
        do  i1 = 1, i2
c         csum=zdotc(kjtop,b(0,1,i1),1,g,1)
          csum = (0d0,0d0)
          do  kjlm = 0, kjtop-1
            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
          enddo
          s(i1,i2) = s(i1,i2) + csum
        enddo
      enddo
C     call tcx('augqs3')

      end subroutine augqs3


      subroutine augqp3(kmax,nlma,ndimh,isp,g,ppi,b,h)

C- Add B+ ppi B to H for non-L-diagonal matrix ppi
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   ppi   :augmentation Pkl-Pkl potential matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   h     :3-center from this augmentation site added to h
Cr Remarks
Cu Updates
Cu 01 Sep 04 folded so into complex potential
Cu 29 Jun 04 (A. Chantis) added 3- center so*Sz*Lz (spherical part)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,nlma,ndimh,isp
      double precision ppi(0:kmax,0:kmax,nlma,nlma,isp)
      double complex b(0:kmax,nlma,ndimh),h(ndimh,ndimh),
     .g(0:kmax,nlma),csum
C ... Local parameters
      integer i1,i2,jlm1,jlm2,k1,k2,kjtop,kjlm

C     call tcn('augqp3')
      kjtop = nlma*(kmax+1)

      do  i2 = 1, ndimh
C   ... g <- ppi*b
        call dpzero(g,2*kjtop)
        do  jlm1 = 1, nlma
          do  jlm2 = 1, nlma
            do  k2 = 0, kmax
              do  k1 = 0, kmax
                g(k1,jlm1) = g(k1,jlm1) +
     .          ppi(k1,k2,jlm1,jlm2,isp)*b(k2,jlm2,i2)
              enddo
            enddo
          enddo
        enddo

C   ... Make dot products with vectors i1
        do  i1 = 1, i2
C         csum=zdotc(kjtop,b(0,1,i1),1,g,1)
          csum = (0d0,0d0)
          do  kjlm = 0, kjtop-1
            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
          enddo
          h(i1,i2) = h(i1,i2) + csum
        enddo
      enddo
C     call tcx('augqp3')
      end subroutine augqp3


      subroutine augq3z(mode1,kmax,nlma,ndimh,isp,g,ppi,b,h,hso)

C- Add B+ ppi B to H for non-L-diagonal, complex matrix ppi
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   ppi   :augmentation Pkl-Pkl potential matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   h     :3-center from this augmentation site added to h
Co   hso   :3-center from this augmentation site added to hso
Cr Remarks
Cu Updates
Cu 03 Feb 05 (A. Chantis) added 3- center spherical so*(LxSx+LySy)
Cu 01 Sep 04 folded so into complex potential
Cu 29 Jun 04 (A. Chantis) added 3- center so*Sz*Lz (spherical part)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,nlma,ndimh,isp,mode1
      double complex ppi(0:kmax,0:kmax,nlma,nlma,isp+2*mode1)
      double complex b(0:kmax,nlma,ndimh),h(ndimh,ndimh),
     .g(0:kmax,nlma),csum,gso(0:kmax,nlma),csum1,hso(ndimh,ndimh)
C ... Local parameters
      integer i1,i2,jlm1,jlm2,k1,k2,kjtop,kjlm

C     call tcn('augqp3')
      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh
C   ... g <- ppi*b
        call dpzero(g,2*kjtop)
        call dpzero(gso,2*kjtop)
        do  jlm1 = 1, nlma
          do  jlm2 = 1, nlma
            do  k2 = 0, kmax
              do  k1 = 0, kmax
                g(k1,jlm1) = g(k1,jlm1) +
     .          ppi(k1,k2,jlm1,jlm2,isp)*b(k2,jlm2,i2)
C     ...  Make LxSx+LySy part of SO
                if (mode1 .eq. 1 .and. isp .eq. 2)
     .          gso(k1,jlm1) = gso(k1,jlm1) + 
     .          ppi(k1,k2,jlm1,jlm2,3)*b(k2,jlm2,i2)
              enddo
            enddo
          enddo
        enddo

C   ... Make dot products with vectors i1
        do  i1 = 1, ndimh
C         csum=zdotc(kjtop,b(0,1,i1),1,g,1)
          csum = (0d0,0d0)
          if (mode1 .eq. 1 .and. isp .eq. 2) csum1 = (0d0,0d0)
          do  kjlm = 0, kjtop-1
            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
C     ...  Make LxSx+LySy part of SO
            if (mode1 .eq. 1 .and. isp .eq. 2)
     .      csum1 = csum1 + dconjg(b(kjlm,1,i1))*gso(kjlm,1)
          enddo
          h(i1,i2) = h(i1,i2) + csum
          if (mode1 .eq. 1 .and. isp .eq. 2) 
     .    hso(i1,i2) = hso(i1,i2) + 0.5d0*csum1
        enddo
      enddo
C     call tcx('augqp3')
      end subroutine augq3z


