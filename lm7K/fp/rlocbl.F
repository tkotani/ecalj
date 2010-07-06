Cchp2 cmdvar= real(8),rv_p_,oppipp
Cchp2 cmdvar= real(8),rv_p_,oppihh
Cchp2 cmdvar= real(8),rv_p_,oppihp
Cchp2 cmdvar= real(8),rv_p_,_oppihp
Cchp2 cmdvar= real(8),rv_p_,oppipp
Cchp3 cmdvar= type(s_rv1),sv_p_,oppi
Cchp3 cmdvar= type(s_rv1),sv_p_,otau
Cchp2 cmdvar= real(8),rv_p_,osighp
Cchp2 cmdvar= real(8),rv_p_,osigpp
Cchp3 cmdvar= type(s_rv1),sv_p_,osig
Cchp2 cmdvar= real(8),rv_p_,oeqhh
Cchp2 cmdvar= real(8),rv_p_,oeqhp
Cchp2 cmdvar= real(8),rv_p_,oeqpp
Cchp2 cmdvar= real(8),rv_p_,oqhh
Cchp2 cmdvar= real(8),rv_p_,oqhp
Cchp2 cmdvar= real(8),rv_p_,oqpp
Cchp3 cmdvar= type(s_rv1),sv_p_,oeqkkl
Cchp3 cmdvar= type(s_rv1),sv_p_,oqkkl
Cchp1 cmdvar= integer,iv_p_,ojcg
Cchp1 cmdvar= integer,iv_p_,oidxcg
Cchp1 cmdvar= real(8),rv_p_,ocy
Cchp1 cmdvar= real(8),rv_p_,ocg
CSFCPP#define F90 1
Cchp3       subroutine rlocbl(ssite,sspec,slat,lfrce,nbas,isp,q,ndham,ndimh,
Cchp3      .nspc,napw,igvapw,iprmb,numq,nevec,evec,ewgt,evl,osig,otau,oppi,
Cchp3      .lcplxp,lekkl,oqkkl,oeqkkl,f)
Cchp3        subroutine rlocbl ( ssite , sspec , slat , lfrce , nbas , isp 
Cchp3      . , q , ndham , ndimh , nspc , napw , igvapw , iprmb , numq , nevec 
Cchp3      . , evec , ewgt , evl , osig , otau , oppi , lcplxp , lekkl , sv_p_oqkkl 
Cchp3      . , oeqkkl , f ) 
Cchp3         subroutine rlocbl ( ssite , sspec , slat , lfrce , nbas , isp 
Cchp3      .  , q , ndham , ndimh , nspc , napw , igvapw , iprmb , numq , nevec 
Cchp3      .  , evec , ewgt , evl , osig , otau , oppi , lcplxp , lekkl , sv_p_oqkkl 
Cchp3      .  , sv_p_oeqkkl , f ) 
Cchp3          subroutine rlocbl ( ssite , sspec , slat , lfrce , nbas , isp 
Cchp3      .   , q , ndham , ndimh , nspc , napw , igvapw , iprmb , numq , nevec 
Cchp3      .   , evec , ewgt , evl , sv_p_osig , otau , oppi , lcplxp , lekkl 
Cchp3      .   , sv_p_oqkkl , sv_p_oeqkkl , f ) 
Cchp3           subroutine rlocbl ( ssite , sspec , slat , lfrce , nbas , isp 
Cchp3      .    , q , ndham , ndimh , nspc , napw , igvapw , iprmb , numq , nevec 
Cchp3      .    , evec , ewgt , evl , sv_p_osig , sv_p_otau , oppi , lcplxp , 
Cchp3      .    lekkl , sv_p_oqkkl , sv_p_oeqkkl , f ) 
           subroutine rlocbl ( ssite , sspec , slat , lfrce , nbas , isp 
     .     , q , ndham , ndimh , nspc , napw , igvapw , iprmb , numq , nevec 
     .     , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp 
     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f ) 





      use m_struc_def  !Cgetarg
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
C- Accumulates the local atomic densities.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa kmxt lmxb rsma
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy alat qlat
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   nbas  :size of basis
Ci   isp   :spin channel
Ci   q     :Bloch wave number
Ci   ndham :leanding dimension of evl
Ci   ndimh :dimension of evec
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial fermi levels
Ci   nevec :number of occupied eigenvectors
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   evl   :eigenvalues
Ci   osig  :overlap matrix of P_kL
Ci   otau  :kinetic energy matrix of P_kL (not used here)
Ci   oppi  :potential matrix of P_kL
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
Co Outputs
Co   oqkkl :local density-matrix; see Remarks
Co   oeqkkl:local part of energy-weighted density matrix
Co   f     :local contribution to forces is added
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cr Remarks
Cr   The qkkl are contractions of the proper density-matrix
Cr      Dij = {sum_n w_n evec*_in evec_jn}
Cr   and the coefficients to the one-center expansion of the wave
Cr   function inside the augmentation sphere
Cr     F~i = Fi + sum_kL C^i_kL (P~kL - PkL)
Cr   As usual, we neglect cross terms when making function products.
Cr   Thus function products are of the form
Cr     F~i F~j = Fi Fj +
Cr             = sum_kLk'L' C^i_kL (P~kL P~k'L' - PkL Pk'L') C^j_k'L'
Cr             = sum_kLk'L' C^i_kL (n1kLk'L' - n2kLk'L') C^j_k'L'
Cr   the qkkl are defined as, e.g.
Cr      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
Cr   so that the local part of the output density is
Cr      n1 - n2 = sum_kLk'L' qpp_kLk'L' (n1kLk'L' - n2kLk'L')
Cu Updates
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   16 Jun 05 Makes spin-off-diagonal density matrix, noncollinear case
Cu   23 Dec 04 Extended to spin-coupled case
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 Jun 00 spin polarized
Cu   25 May 00 Adapted from nfp rloc_q.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
#if MPI
      include "mpif.h"
#if MPE
      include "mpef.h"
#endif
      integer procid,master
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1), pid
      double precision starttime, endtime
      character*120 strn
      logical mlog,cmdopt
#endif
c#if MPE
cC Event numbers:
c      include "events.ins"
c#endif
      integer lfrce,nbas,isp,ndimh,nspc,numq,nevec,lcplxp,lekkl,
     .iprmb(1),ndham,napw,igvapw(3,napw)
Cchp3       integer oqkkl(3,1),oeqkkl(3,1),osig(3,1),otau(3,1),oppi(3,1)
Cchp3        integer:: oeqkkl(3,1) , osig(3,1) , otau(3,1) , oppi(3,1) 
Cchp3         integer:: osig(3,1) , otau(3,1) , oppi(3,1) 
Cchp3          integer:: otau(3,1) , oppi(3,1) 
Cchp3           integer:: oppi(3,1) 
           type(s_rv1) :: sv_p_oppi(3,1)

          type(s_rv1) :: sv_p_otau(3,1)

         type(s_rv1) :: sv_p_osig(3,1)

        type(s_rv1) :: sv_p_oeqkkl(3,1)

       type(s_rv1) :: sv_p_oqkkl(3,1)

      real(8):: q(3) , ewgt(numq,nevec) , evl(ndham,isp) , f(3,nbas,numq)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex evec(ndimh,nspc,nevec)
C ... Local parameters or process-shared variables
#if MPI
      double precision, dimension(:), allocatable :: buffer
      integer, dimension(:), allocatable :: bproc
      integer nelt(3),lgunit,lmxh,nlmh,nsp,i
#endif

      integer is,nlmbx,nlmx,ktop0,npmx,nkap0,n0
      parameter (nlmbx=25,  nkap0=3, n0=10) !npmx=32,
      integer kmaxx,nlmax,igetss,nglob,nlmto !mp,
Cchp1       integer ocg,ocy,oidxcg,ojcg
Cchp1        integer:: ocy , oidxcg , ojcg 
Cchp1         integer:: oidxcg , ojcg 
Cchp1          integer:: ojcg 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)

c      integer oib(npmx),oidb(npmx),oixf(npmx),o
      double precision alat,qlat(3,3)
C ... Local process-specific variables
      integer ia,isa,ivec,kmax,lmxa,nlma,ip,oxf,
     .lmxha,nlmha,nkaph,ispc,ksp
Cchp2       integer ob,odb,oppihp,oppipp,osighp,osigpp,
Cchp2      .oqhh,oqhp,oqpp,oeqhh,oeqhp,oeqpp
Cchp2        integer:: ob , odb , oppihp , oppipp , osighp , osigpp , oqhh 
Cchp2      . , oqhp , oeqhh , oeqhp , oeqpp 
Cchp2         integer:: ob , odb , oppihp , oppipp , osighp , osigpp , oqhh 
Cchp2      .  , oeqhh , oeqhp , oeqpp 
Cchp2          integer:: ob , odb , oppihp , oppipp , osighp , osigpp , oeqhh 
Cchp2      .   , oeqhp , oeqpp 
Cchp2           integer:: ob , odb , oppihp , oppipp , osighp , osigpp , oeqhh 
Cchp2      .    , oeqhp 
Cchp2            integer:: ob , odb , oppihp , oppipp , osighp , osigpp , oeqhh 
Cchp2             integer:: ob , odb , oppihp , oppipp , osighp , osigpp 
Cchp2              integer:: ob , odb , oppihp , oppipp , osighp 
Cchp2               integer:: ob , odb , oppihp , oppipp 
Cchp2                integer:: ob , odb , oppihp 
                integer:: ob , odb 
                real(8),pointer :: rv_p_oppihp(:)

               real(8),pointer :: rv_p_oppipp(:)

              real(8),pointer :: rv_p_osighp(:)

             real(8),pointer :: rv_p_osigpp(:)

            real(8),pointer :: rv_p_oeqhh(:)

           real(8),pointer :: rv_p_oeqhp(:)

          real(8),pointer :: rv_p_oeqpp(:)

         real(8),pointer :: rv_p_oqhh(:)

        real(8),pointer :: rv_p_oqhp(:)

       real(8),pointer :: rv_p_oqpp(:)

      double precision pa(3),rsma,pi,tpiba
CSFCPP#if F90
      complex(8),allocatable:: cPkL(:),da(:),wk(:)
CSFCPP#else
CSFCPP      integer ib,nlmhb,ik
CSFCPP      double precision p(3),xx
CSFCPP      parameter (nlmx=49,ktop0=10)
CSFCPP      double complex cPkL((ktop0+1)*nlmx),da((ktop0+1)*nlmx*3),
CSFCPP     .b0((ktop0+1)*nlmx*nlmbx),db0((ktop0+1)*nlmx*nlmbx*3),
CSFCPP     .wk((ktop0+1)*nlmx),bos((ktop0+1)*nlmbx)
CSFCPP      integer norb,nkapi
CSFCPP      integer lh(nkap0),ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
CSFCPP      double precision eh(n0,nkap0),rsmh(n0,nkap0)
CSFCPP#endif
C$    integer mp_numthreads,mp_my_threadnum
C ... Heap
      integer w(1)
      common /w/ w
      integer:: iaini,iaend
      complex(8),allocatable::w_ob(:),w_odb(:),w_oxf(:)

#if MPI
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      nsp = nglob('nsp')
      if (mlog) then
        do  pid = 0, numprocs-1
          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,
     .    MPI_COMM_WORLD,ierr)
          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
     .    MPI_COMM_WORLD,ierr)
        enddo
      endif
#endif

      if (nevec .le. 0) return
      call tcn('rlocbl')

C --- Setup ---
      nkaph = nglob('nkaph')

Cchp1       ocg=slat%ocg
       rv_p_ocg => slat%rv_p_ocg 

Cchp1       ojcg=slat%ojcg
       iv_p_ojcg => slat%iv_p_ojcg 

Cchp1       oidxcg=slat%oidxcg
       iv_p_oidxcg => slat%iv_p_oidxcg 

Cchp1       ocy=slat%ocy
       rv_p_ocy => slat%rv_p_ocy 


c      call defi(o,1)
C ... Find maximum sizes needed to allocate strux; allocate them
      nlmax = 0
      kmaxx = 0
      do  ia = 1, nbas
        isa = int(ssite(ia)%spec)


        lmxa=sspec(isa)%lmxa
        kmax=sspec(isa)%kmxt

        nlma = (lmxa+1)**2
        kmaxx = max(kmaxx,kmax)
        nlmax = max(nlmax,nlma)
      enddo
      nlmto = ndimh-napw
C     Needed for PW part

      alat=slat%alat
c      i_copy_size=size(slat%qlat)
c      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      qlat=slat%qlat

      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat

CSFCPP#if F90
C ... Allocate workspace for augmentation arrays
      nlmx  = nlmax
      ktop0 = kmaxx
      allocate(cPkL((ktop0+1)*nlmx),da((ktop0+1)*nlmx*3),
     .wk((ktop0+1)*nlmx))
C      allocate(cPkL((ktop0+1)*nlmx),da((ktop0+1)*nlmx*3),
C     .  b0((ktop0+1)*nlmx*nlmbx),db0((ktop0+1)*nlmx*nlmbx*3),
C     .  wk((ktop0+1)*nlmx),bos((ktop0+1)*nlmbx))
CSFCPP#else
CSFCPP      if (napw .gt. 0) call rx('rlocbl: PWs require F90 compiler')
CSFCPP#endif
      if (nlmax .gt. nlmx) call rxi('rlocbl: nlmx.lt.nlma=',nlmax)
      if (kmaxx .gt. ktop0) call rxi('rlocbl: ktop0.lt.kmax=',kmax)

C ... Allocate workspace for strux
c      mp = 1
C$    mp = mp_numthreads()

      allocate(w_ob(ndimh*nlmax*(kmaxx+1)))
      allocate(w_odb(ndimh*nlmax*(kmaxx+1)*3))
      if (lfrce .ne. 0) then
        allocate(w_oxf(3*nbas*numq))
        w_oxf=0d0
      endif  
c      do  ip = 1, mp
c        call defcc(oib(ip), ndimh*nlmax*(kmaxx+1))
c        call defcc(oidb(ip),ndimh*nlmax*(kmaxx+1)*3)
c        if (lfrce .ne. 0) call defrr(oixf(ip),-3*nbas*numq)
c      enddo

C --- Loop over augmentation sites ---
#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_RLOCBL,procid,"rlocbl")
#endif
      allocate (bproc(0:numprocs), stat=ierr)
      call dstrbp(nbas,numprocs,1,bproc(0))
c      do  ia = bproc(procid), bproc(procid+1)-1
      iaini=bproc(procid)
      iaend=bproc(procid+1)-1
#else
      iaini=1
      iaend=nbas
#endif
      do ia = iaini,iaend
#if MPI
        if (mlog .and. ia .eq. bproc(procid)) then
          call gettime(datim)
          call awrit4(' rlocbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting atoms %i to %i',' ',256,lgunit(3),
     .    procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
#endif
        ip=1
C$      ip=mp_my_threadnum()+1

        isa=ssite(ia)%spec
c        i_copy_size=size(ssite(ia)%pos)
c        call dcopy(i_copy_size,ssite(ia)%pos,1,pa,1)
        pa=ssite(ia)%pos(1:3)

        lmxa=sspec(isa)%lmxa
        lmxha=sspec(isa)%lmxb
        kmax=sspec(isa)%kmxt
        rsma=sspec(isa)%rsma

        nlmha = (lmxha+1)**2
        nlma  = (lmxa+1)**2
        if (lmxa .eq. -1) goto 10

        rv_p_oqpp => sv_p_oqkkl(1,ia)%v
        rv_p_oqhp => sv_p_oqkkl(2,ia)%v
        rv_p_oqhh => sv_p_oqkkl(3,ia)%v
        rv_p_oeqpp => sv_p_oeqkkl(1,ia)%v
        rv_p_oeqhp => sv_p_oeqkkl(2,ia)%v
        rv_p_oeqhh => sv_p_oeqkkl(3,ia)%v
        rv_p_oppipp => sv_p_oppi(1,ia)%v
        rv_p_osigpp => sv_p_osig(1,ia)%v
        rv_p_oppihp => sv_p_oppi(2,ia)%v
        rv_p_osighp => sv_p_osig(2,ia)%v

C   --- Strux to expand all orbitals and their gradients at site ia ---
c        ob = oib(ip)
c        odb = oidb(ip)
c        oxf = oixf(ip)

CSFCPP#if F90
Cchp1         call bstrux(1,slat,ssite,sspec,w(ocg),w(oidxcg),w(ojcg),w(ocy),
Cchp1      .  iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,napw,igvapw,w(ob),
Cchp1      .  w(odb))
Cchp1          call bstrux ( 1 , slat , ssite , sspec , rv_p_ocg , w ( oidxcg 
Cchp1      .   ) , w ( ojcg ) , w ( ocy ) , iprmb , nbas , ia , pa , rsma , 
Cchp1      .   q , kmax , nlma , ndimh , napw , igvapw , w ( ob ) , w ( odb 
Cchp1      .   ) ) 
Cchp1           call bstrux ( 1 , slat , ssite , sspec , rv_p_ocg , w ( oidxcg 
Cchp1      .    ) , w ( ojcg ) , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q 
Cchp1      .    , kmax , nlma , ndimh , napw , igvapw , w ( ob ) , w ( odb ) 
Cchp1      .    ) 
Cchp1            call bstrux ( 1 , slat , ssite , sspec , rv_p_ocg , iv_p_oidxcg 
Cchp1      .     , w ( ojcg ) , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
Cchp1      .     kmax , nlma , ndimh , napw , igvapw , w ( ob ) , w ( odb ) ) 
            call bstrux ( 1 , slat , ssite , sspec , rv_p_ocg , iv_p_oidxcg 
     .      , iv_p_ojcg , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
     .      kmax , nlma , ndimh , napw , igvapw , w_ob , w_odb ) 




CSFCPP#else
CSFCPPC       Zero out strux to eliminate contributions from local orbitals
CSFCPP        call dpzero(w(ob),(kmax+1)*nlma*ndimh*2)
CSFCPP        call dpzero(w(odb),(kmax+1)*nlma*ndimh*2*3)
CSFCPP        do  ib = 1, nbas
CSFCPP
CSFCPP          is=ssite(ib)%spec
CSFCPP          i_copy_size=size(ssite(ib)%pos)
CSFCPP          call dcopy(i_copy_size,ssite(ib)%pos,1,p,1)
CSFCPP
CSFCPPC         This line augments no local orbitals
CSFCPPC         ik = 1
CSFCPPC         This line augments onsite extended local orbitals only
CSFCPPC         if (ia .eq. ib) ik = 2
CSFCPPC         This line augments extended local orbitals all sites
CSFCPP          ik = 2
CSFCPP          call uspecb(0,ik,sspec,is,is,lh,rsmh,eh,nkapi)
CSFCPPC         Position in h; l,k indices for orbitals connected w/ ib
CSFCPP          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,xx)
CSFCPPC         Get list of blocks with envelope functions
CSFCPPC         call gtbsl1(8+16,norb,ltab,ktab,xx,xx,ntab,blks)
CSFCPPC         Loop over blocks of envelope functions
CSFCPP          do  ik = 1, nkapi
CSFCPP            nlmhb = (lh(ik)+1)**2
CSFCPP            if (nlmhb .gt. nlmbx) call rxi('rlocbl: need nlmbx',nlmhb)
CSFCPP            if (nlmhb .gt. nlma .and. ia.eq.ib)
CSFCPP     .      call rx('rlocbl: nlmh .gt. nlma')
CSFCPP
CSFCPP            call hxpgbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmhb,nlma,
CSFCPP     .      kmax,nlmbx,nlma,w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,b0,
CSFCPP     .      db0)
CSFCPP            if (ib .eq. ia) then
CSFCPP              call hxpos(rsmh(1,ik),rsma,eh(1,ik),kmax,nlmhb,kmax,bos)
CSFCPP              call paugq2(kmax,nlmhb,nlma,bos,b0)
CSFCPP            endif
CSFCPP            call prlcb1(1,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,nlma,
CSFCPP     .      kmax,b0,db0,w(ob),w(odb))
CSFCPP          enddo
CSFCPP        enddo
CSFCPP#endif

C   --- Loop over eigenstates ---
C       In noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       whereas ispc is spin index in the noncoll case, but 1 for coll.
        do  ivec = 1, nevec
          do  ispc = 1, nspc
            ksp = max(ispc,isp)

C     ... Pkl expansion of eigenvector
            call rlocb1(ndimh,nlma,kmax,evec(1,ispc,ivec),w_ob,cPkL)

C      ... needed for off-diagonal density matrix ... comment out now
C          if (nspc .eq. 2 .and. ispc .eq. 1) then
CC#ifdef F90
C            allocate(cPkL1((kmax+1)*nlma))
C            call zcopy((kmax+1)*nlma,cPkL,1,cPkL1,1)
CC#elseC
CC            call rx('rlocbl: noncoll mode requires f90')
CC#endif
C          endif

C     ... Add to local density coefficients for one state
Cchp2             call prlcb3(0,kmax,nlma,ksp,cPkL,numq,ewgt(1,ivec),w,w(oqpp))
             call prlcb3 ( 0 , kmax , nlma , ksp , cpkl , numq , ewgt ( 1 
     .       , ivec ) , w , rv_p_oqpp ) 

Cchp2             call prlcb2(0,ia,nkaph,iprmb,nlmha,kmax,nlma,ksp,cPkL,nlmto,
Cchp2      .      evec(1,ispc,ivec),ewgt(1,ivec),numq,w,w(oqhh),w(oqhp))
Cchp2              call prlcb2 ( 0 , ia , nkaph , iprmb , nlmha , kmax , nlma , 
Cchp2      .       ksp , cpkl , nlmto , evec ( 1 , ispc , ivec ) , ewgt ( 1 , ivec 
Cchp2      .       ) , numq , w , w ( oqhh ) , rv_p_oqhp ) 
              call prlcb2 ( 0 , ia , nkaph , iprmb , nlmha , kmax , nlma , 
     .        ksp , cpkl , nlmto , evec ( 1 , ispc , ivec ) , ewgt ( 1 , ivec 
     .        ) , numq , w , rv_p_oqhh , rv_p_oqhp ) 


            if (lekkl .eq. 1) then
Cchp2               call prlcb3(1,kmax,nlma,ksp,cPkL,numq,ewgt(1,ivec),
Cchp2      .        evl(ivec,isp),w(oeqpp))
               call prlcb3 ( 1 , kmax , nlma , ksp , cpkl , numq , ewgt ( 1 
     .         , ivec ) , evl ( ivec , isp ) , rv_p_oeqpp ) 

Cchp2               call prlcb2(1,ia,nkaph,iprmb,nlmha,kmax,nlma,ksp,cPkL,nlmto,
Cchp2      .        evec(1,ispc,ivec),ewgt(1,ivec),numq,evl(ivec,isp),
Cchp2      .        w(oeqhh),w(oeqhp))
Cchp2                call prlcb2 ( 1 , ia , nkaph , iprmb , nlmha , kmax , nlma , 
Cchp2      .         ksp , cpkl , nlmto , evec ( 1 , ispc , ivec ) , ewgt ( 1 , ivec 
Cchp2      .         ) , numq , evl ( ivec , isp ) , w ( oeqhh ) , rv_p_oeqhp ) 
                call prlcb2 ( 1 , ia , nkaph , iprmb , nlmha , kmax , nlma , 
     .          ksp , cpkl , nlmto , evec ( 1 , ispc , ivec ) , ewgt ( 1 , ivec 
     .          ) , numq , evl ( ivec , isp ) , rv_p_oeqhh , rv_p_oeqhp ) 


            endif

C          ... needed for off-diagonal density matrix ... comment out now
CC#ifdef F90
C          if (nspc .eq. 2 .and. ispc .eq. 2) then
C            call prlc23(kmax,nlma,2,cPkL1,cPkL,numq,ewgt(1,ivec),
C     .        w(oqpp))
C            call prlc22(ia,nkaph,iprmb,nlmha,kmax,nlma,2,cPkL,ndimh,
C     .        evec(1,1,ivec),evec(1,2,ivec),ewgt(1,ivec),numq,
C     .        w(oqhh),w(oqhp))
C            deallocate(cPkL1)
C          endif
CC#endif

C     ... Contribution to forces
            if (lfrce .ne. 0) then
              call rxx(nspc.ne.1,'forces not implemented in noncoll case')
Cchp2               call flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto,
Cchp2      .        ndimh,iprmb,ksp,evl(ivec,isp),evec(1,ispc,ivec),ewgt(1,
Cchp2      .        ivec),numq,cPkL,w_odb,da,wk,w(oppipp),w(oppipp),
Cchp2      .        w(osigpp),w(oppihp),w(oppihp),w(osighp),lcplxp,w_oxf)
Cchp2                call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma 
Cchp2      .         , lmxa , nlmto , ndimh , iprmb , ksp , evl ( ivec , isp ) , evec 
Cchp2      .         ( 1 , ispc , ivec ) , ewgt ( 1 , ivec ) , numq , cpkl , w_odb 
Cchp2      .         , da , wk , w ( oppipp ) , w ( oppipp ) , rv_p_osigpp , w ( oppihp 
Cchp2      .         ) , w ( oppihp ) , w ( osighp ) , lcplxp , w_oxf ) 
Cchp2                 call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma 
Cchp2      .          , lmxa , nlmto , ndimh , iprmb , ksp , evl ( ivec , isp ) , evec 
Cchp2      .          ( 1 , ispc , ivec ) , ewgt ( 1 , ivec ) , numq , cpkl , w_odb 
Cchp2      .          , da , wk , w ( oppipp ) , w ( oppipp ) , rv_p_osigpp , w ( oppihp 
Cchp2      .          ) , w ( oppihp ) , rv_p_osighp , lcplxp , w_oxf ) 
Cchp2                  call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma 
Cchp2      .           , lmxa , nlmto , ndimh , iprmb , ksp , evl ( ivec , isp ) , evec 
Cchp2      .           ( 1 , ispc , ivec ) , ewgt ( 1 , ivec ) , numq , cpkl , w_odb 
Cchp2      .           , da , wk , rv_p_oppipp , rv_p_oppipp , rv_p_osigpp , w ( oppihp 
Cchp2      .           ) , w ( oppihp ) , rv_p_osighp , lcplxp , w_oxf ) 
                  call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma 
     .            , lmxa , nlmto , ndimh , iprmb , ksp , evl ( ivec , isp ) , evec 
     .            ( 1 , ispc , ivec ) , ewgt ( 1 , ivec ) , numq , cpkl , w_odb 
     .            , da , wk , rv_p_oppipp , rv_p_oppipp , rv_p_osigpp , rv_p_oppihp 
     .            , rv_p_oppihp , rv_p_osighp , lcplxp , w_oxf ) 




            endif
          enddo
        enddo
   10   continue
C ... end loop over ia
      enddo
#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_RLOCBL,procid,"rlocbl")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_BCAST,procid,"broadcast")
#endif
      do  pid = 0, numprocs-1
        do  ia = bproc(pid), bproc(pid+1)-1
          is = int(ssite(ia)%spec)


          lmxa=sspec(is)%lmxa
          lmxh=sspec(is)%lmxb
          kmax=sspec(is)%kmxt

          nlma = (lmxa+1)**2
          nlmh = (lmxh+1)**2
          nelt(1) = (kmax+1)*(kmax+1)*nlma*nlma
          nelt(2) = (kmax+1)*nkaph*nlma*nlmh
          nelt(3) = nkaph*nkaph*nlmh*nlmh
          do  i = 1, 3
Cchp3             call MPI_BCAST(w(oqkkl(i,ia)),nelt(i)*numq*nsp*nspc,
Cchp3      .      MPI_DOUBLE_PRECISION,pid,MPI_COMM_WORLD,ierr)
             call mpi_bcast ( sv_p_oqkkl( i , ia )%v , nelt ( i ) * numq * 
     .       nsp * nspc , mpi_double_precision , pid , mpi_comm_world , ierr 
     .       ) 

            if (mlog) then
              call gettime(datim)
              call awrit6(' rlocbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' bcast w(oqkkl(%i,%i)) %i d.p. numbers'//
     .        ' from process %i on '
     .        //shortname(pid)(1:namelen(pid)),' ',
     .        256,lgunit(3),procid,numprocs,i,ia,
     .        nelt(i)*numq*nsp*nspc,pid)
            endif
            if (lekkl .eq. 1) then
Cchp3               call MPI_BCAST(w(oeqkkl(i,ia)),nelt(i)*numq*nsp*nspc,
Cchp3      .        MPI_DOUBLE_PRECISION,pid,MPI_COMM_WORLD,ierr)
               call mpi_bcast ( sv_p_oeqkkl( i , ia )%v , nelt ( i ) * numq 
     .         * nsp * nspc , mpi_double_precision , pid , mpi_comm_world , 
     .         ierr ) 

              if (mlog) then
                call gettime(datim)
                call awrit6(' rlocbl '//datim//' Process %i of %i on '
     .          //shortname(procid)(1:namelen(procid))//
     .          ' bcast w(oeqkkl(%i,%i)) %i d.p. numbers'//
     .          ' from process %i on '
     .          //shortname(pid)(1:namelen(pid)),' ',
     .          256,lgunit(3),procid,numprocs,i,ia,
     .          nelt(i)*numq*nsp*nspc,pid)
              endif
            endif
          enddo
        enddo
      enddo
      if (lfrce .ne. 0) then
        allocate(buffer(3*nbas*numq), stat=ierr)
c        call MPI_ALLREDUCE(w(oixf(1)),buffer,3*nbas*numq,
        call MPI_ALLREDUCE(w_oxf,buffer,3*nbas*numq,
     .  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit2(' rlocbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' allreduce forces',' ',256,lgunit(3),
     .    procid,numprocs)
        endif
        call daxpy(3*nbas*numq,1d0,buffer,1,f,1)
        deallocate(buffer, stat=ierr)
      endif
      deallocate(bproc, stat=ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BCAST,procid,"broadcast")
#endif
#else
      if (lfrce .ne. 0) then
ctakao think mp=1 is assumed in this routine 3rdJune2010
c        do  ip = 1, mp
          call daxpy(3*nbas*numq,1d0,w_oxf,1,f,1)
c        enddo
      endif
#endif

CSFCPP#if F90
      deallocate(cPkL,da,wk)
CSFCPP#endif
c      call rlse(o)
      deallocate(w_ob,w_odb)
      if(lfrce /=0) deallocate(w_oxf)
      call tcx('rlocbl')
      end subroutine rlocbl


      subroutine rlocb1(ndimh,nlma,kmax,evec,b,cPkL)

C- Pkl expansion of wave function at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimh :dimension of evec
Ci   nlma  :augmentation L-cutoff in PkL expansion
Ci   kmax  :k- cutoff in PkL expansion
Ci   evec  :eigenvector coefficients
Ci   b     :strux to expand of orbitals from other sites in PkL
Ci         :b = b(ndimh,nlma,0:kmax)
Co Outputs
Co   cPkL  :coefficients to PkL expansion of evec
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,ndimh,nlma
      double complex b(ndimh,nlma,0:kmax),cPkL(0:kmax,nlma),evec(ndimh)
C ... Local parameters
      integer i,k,ilma


      call tcn('rlocb1')
C     call zprm('b',2,b,ndimh,ndimh,nlma*(kmax+1))
      call dpzero(cPkL, 2*(kmax+1)*nlma)
      do  k = 0, kmax
        do  ilma = 1, nlma
          do  i = 1, ndimh
            cPkL(k,ilma) = cPkL(k,ilma) + evec(i)*b(i,ilma,k)
          enddo
        enddo
      enddo
C     call zprm('cPkL',2,cPkL,kmax+1,kmax+1,nlma)
      call tcx('rlocb1')
      end subroutine rlocb1


      subroutine prlcb2(job,ia,nkaph,iprmb,nlmha,kmax,nlma,isp,cPkL,
     .nlmto,evec,ewgt,numq,evl,qhh,qhp)

C- Add one and two-center terms to density coeffs
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   ia    :site of augmentation
Ci   nkaph :dimensions qhh,qhp
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions qhh,qhp
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   isp   :spin channel
Ci   cPkL  :PkL expansion eigenvector at site ia.
Ci   nlmto :dimension of lmto component of basis
Ci   evec  :eigenvector
Ci   numq  :number of trial fermi levels
Ci   ewgt  :eigenvector weight
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qhh   :one-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Co   qhp   :two-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
Cu   27 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,ia,kmax,nkaph,isp,nlmto,nlma,nlmha,numq,iprmb(1)
      double precision qhh(nkaph,nkaph,nlmha,nlmha,numq,isp),
     .qhp(nkaph,0:kmax,nlmha,nlma,numq,isp),ewgt(numq),evl
      double complex evec(nlmto),cPkL(0:kmax,nlma)
C ... Local parameters
      integer i1,i2,ilm1,ilm2,ilma,io1,io2,iq,k,ik1,ik2,
     .l1,l2,n0,nkap0,nlm11,nlm12,nlm21,nlm22
      parameter (n0=10,nkap0=3)
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .blks(n0*nkap0),ntab(n0*nkap0)
      double precision xx

      if (nlmto .eq. 0) return
      call tcn('prlcb2')

C --- Loop over all orbitals centered at site ia, incl. local orbs ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

      do  io1 = 1, norb

        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
C       i1 = hamiltonian offset for first orbital in block
        i1 = offl(io1)
        do  ilm1 = nlm11, nlm12
          i1 = i1+1
C     ... Accumulate products H*Pkl
          do  iq = 1, numq
            if (job .eq. 0) then
              do  k = 0, kmax
                do  ilma = 1, nlma
                  qhp(ik1,k,ilm1,ilma,iq,isp)= qhp(ik1,k,ilm1,ilma,iq,isp)
     .            + 2d0*dconjg(evec(i1))*cPkL(k,ilma)*ewgt(iq)
                enddo
              enddo
            else
              do  k = 0, kmax
                qhp(ik1,k,ilm1,ilm1,iq,isp)= qhp(ik1,k,ilm1,ilm1,iq,isp)
     .          + evl*2d0*dconjg(evec(i1))*cPkL(k,ilm1)*ewgt(iq)
              enddo
            endif
          enddo

C     ... Accumulate products H*H
          do  io2 = 1, norb

            l2  = ltab(io2)
            ik2 = ktab(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks(io2)-1
C           i2 = orbital index in iprmb order
            i2 = offl(io2)
            do  ilm2 = nlm21, nlm22
              i2 = i2+1
              if (job .eq. 0) then
                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,iq,isp) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,isp) +
     .            dconjg(evec(i1))*evec(i2)*ewgt(iq)
                enddo
              elseif (job .eq. 1 .and. ilm1 .eq. ilm2) then
                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,iq,isp) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,isp) +
     .            evl*dconjg(evec(i1))*evec(i2)*ewgt(iq)
                enddo
              endif
            enddo
          enddo

        enddo
      enddo

      call tcx('prlcb2')
      end subroutine prlcb2


      subroutine prlcb3(job,kmax,nlma,isp,cPkL,numq,ewgt,evl,qpp)

C- Add to local density coefficients for one state
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   kmax  :polynomial cutoff in PkL expansion
Ci   nlma  :L cutoff in PkL expansion
Ci   isp   :spin channel
Ci   cPkL  :coefficients to PkL expansion of evec
Ci   numq  :number of trial fermi levels
Ci   ewgt  :eigenvector weights
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qpp   :local density matrix for PkL expansion (job=0)
Co         :energy-weighted local density matrix (job=1)
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,kmax,nlma,isp,numq
      double complex cPkL(0:kmax,nlma)
      double precision qpp(0:kmax,0:kmax,nlma,nlma,numq,isp),ewgt(numq),
     .evl
C ... Local parameters
      double precision fac
      integer iq,ilm2,ilm1,k1,k2

      call tcn('prlcb3')
      do  iq = 1, numq
        fac = ewgt(iq)
        if (job .eq. 1) fac = evl*ewgt(iq)
        do  ilm2 = 1, nlma
          do  ilm1 = 1, nlma
            do  k1 = 0, kmax
              do  k2 = 0, kmax
                qpp(k1,k2,ilm1,ilm2,iq,isp)= qpp(k1,k2,ilm1,ilm2,iq,isp)
     .          + fac*dconjg(cPkL(k1,ilm1))*cPkL(k2,ilm2)
              enddo
            enddo
          enddo
        enddo
      enddo
      call tcx('prlcb3')

      end subroutine prlcb3


C  ... needed for off-diagonal density matrix ... comment out now
C      subroutine prlc22(ia,nkaph,iprmb,nlmha,kmax,nlma,ioff,cPkL,ndimh,
C     .  evec1,evec2,ewgt,numq,qhh,qhp)
CC- Add one and two-center terms to density coeffs
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ia    :site of augmentation
CCi   nkaph :dimensions qhh,qhp
CCi   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
CCi   nlmha :dimensions qhh,qhp
CCi   kmax  :polynomial cutoff
CCi   nlma  :augmentation L-cutoff
CCi   ioff  :spin channel
CCi   cPkL  :PkL expansion eigenvector at site ia.
CCi   ndimh :dimension of hamiltonian
CCi   evec1 :left eigenvector
CCi   evec2 :right eigenvector
CCi   numq  :number of trial fermi levels
CCi   ewgt  :eigenvector weight
CCo Outputs
CCo   qhh   :one-center density-matrix, spin 12 block
CCo   qhp   :two-center density-matrix, spin 12 block
CCr Remarks
CCu Updates
CCu   27 Aug 01 Extended to local orbitals.
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer ia,kmax,nkaph,ioff,ndimh,nlma,nlmha,numq,iprmb(1)
C      double precision ewgt(numq)
C      double complex qhh(nkaph,nkaph,nlmha,nlmha,numq,ioff),
C     .  qhp(nkaph,0:kmax,nlmha,nlma,numq,ioff)
C      double complex evec1(ndimh),evec2(ndimh),cPkL(0:kmax,nlma)
CC ... Local parameters
C      integer i1,i2,ilm1,ilm2,ilma,io1,io2,iq,k,ik1,ik2,
C     .  l1,l2,n0,nkap0,nlm11,nlm12,nlm21,nlm22
C      parameter (n0=10,nkap0=3)
C      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
C     .  blks(n0*nkap0),ntab(n0*nkap0)
C      double precision xx
C
C      call tcn('prlc22')
C
CC --- Loop over all orbitals centered at site ia, incl. local orbs ---
C      call orbl(ia,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,xx)
CC     Block into groups of consecutive l
C      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
C
C      do  io1 = 1, norb
C
C        l1  = ltab(io1)
C        ik1 = ktab(io1)
C        nlm11 = l1**2+1
C        nlm12 = nlm11 + blks(io1)-1
CC       i1 = hamiltonian offset for first orbital in block
C        i1 = offl(io1)
C        do  ilm1 = nlm11, nlm12
C          i1 = i1+1
CC     ... Accumulate products H*Pkl
C          do  iq = 1, numq
C            do  k = 0, kmax
C              do  ilma = 1, nlma
C                qhp(ik1,k,ilm1,ilma,iq,ioff) =
C     .          qhp(ik1,k,ilm1,ilma,iq,ioff)
C     .            + 2d0*dconjg(evec1(i1))*cPkL(k,ilma)*ewgt(iq)
C              enddo
C            enddo
C          enddo
C
CC     ... Accumulate products H*H
C          do  io2 = 1, norb
C
C            l2  = ltab(io2)
C            ik2 = ktab(io2)
C            nlm21 = l2**2+1
C            nlm22 = nlm21 + blks(io2)-1
CC           i2 = orbital index in iprmb order
C            i2 = offl(io2)
C            do  ilm2 = nlm21, nlm22
C              i2 = i2+1
C              do  iq = 1, numq
C                qhh(ik1,ik2,ilm1,ilm2,iq,ioff) =
C     .          qhh(ik1,ik2,ilm1,ilm2,iq,ioff) +
C     .          dconjg(evec1(i1))*evec2(i2)*ewgt(iq)
C              enddo
C            enddo
C          enddo
C
C        enddo
C      enddo
C
C      call tcx('prlc22')
C
C      end
C      subroutine prlc23(kmax,nlma,ioff,cPkL1,cPkL2,numq,ewgt,qpp)
CC- Add to local density coefficients for one state, off-diagonal spin block
CC ----------------------------------------------------------------------
CCi Inputs
CCi   kmax  :polynomial cutoff in PkL expansion
CCi   nlma  :L cutoff in PkL expansion
CCi   ioff  :offset in qpp array
CCi   cPkL1 :coefficients to PkL expansion of left evec
CCi   cPkL2 :coefficients to PkL expansion of right evec
CCi   numq  :number of trial fermi levels
CCi   ewgt  :eigenvector weights
CCo Outputs
CCo   qpp   :local density matrix for PkL expansion, spin 12 block
CCr Remarks
CCu Updates
CCu   16 Jun 05 Adapted from prlcb3
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer kmax,nlma,ioff,numq
C      double complex cPkL1(0:kmax,nlma),cPkL2(0:kmax,nlma)
C      double precision ewgt(numq)
C      double complex qpp(0:kmax,0:kmax,nlma,nlma,numq,ioff)
CC ... Local parameters
C      integer iq,ilm2,ilm1,k1,k2
C
C      call tcn('prlc23')
C      do  iq = 1, numq
C        do  ilm2 = 1, nlma
C          do  ilm1 = 1, nlma
C            do  k1 = 0, kmax
C              do  k2 = 0, kmax
C                qpp(k1,k2,ilm1,ilm2,iq,ioff) =
C     .          qpp(k1,k2,ilm1,ilm2,iq,ioff)
C     .            + ewgt(iq)*dconjg(cPkL1(k1,ilm1))*cPkL2(k2,ilm2)
C              enddo
C            enddo
C          enddo
C        enddo
C      enddo
C      call tcx('prlc23')
C      end

