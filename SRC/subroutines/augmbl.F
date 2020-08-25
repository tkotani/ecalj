      subroutine augmbl ( mode , ssite , sspec , slat , isp , lcplxp 
     .   , q , sv_p_osig , sv_p_otau , sv_p_oppi , iprmb , ndimh , napw 
     .   , igapw , h , hso , s )
      use m_work,only: rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy
      use m_struc_def           !Cgetarg
      use m_globalvariables
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
      integer:: mode , lcplxp , isp , ndimh , napw , igapw(3,napw), iprmb(ndimh)
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
      complex(8),allocatable:: b(:)
      double complex g(lofb)
      integer:: ia , isa , kmax , lmxa , lmxha , nbas , nglob , nlma 
     .    , nlmha , nkaph , mode1 , nlmto
      double precision rsma,pa(3),xx,alat,qlat(3,3),vol
      integer:: initbas, endbas
      logical:: debug=.false.
#if MPI
      integer, dimension(:), allocatable :: bproc
      complex(8) ,allocatable :: h_zv(:),
     & hbuf_zv(:),  s_zv(:), sbuf_zv(:),  hso_zv(:)
#endif
C--------------------------
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
      nbas  = globalvariables%nbas
      nkaph = globalvariables%nkaph
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      nlmto = ndimh-napw
      alat=slat%alat
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      vol=slat%vol
      allocate (b(lofb*ndimh))
C --- Loop over augmentation sites --- ---
#if MPI
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
        call bstrux ( 0 , slat , ssite , sspec , rv_a_ocg , iv_a_oidxcg
     .      , iv_a_ojcg , rv_a_ocy , iprmb , nbas , ia , pa , rsma , q , 
     .      kmax , nlma , ndimh , napw , igapw , b , xx )
        if(debug) print *,'end of bstrux'
C   --- Add 1-center and 2-center terms ---
#if MPI
        if (lcplxp .eq. 0) then
          call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
     .      , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s_zv , h_zv )
        else
          call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
     .      , nlmha , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s_zv , h_zv , hso_zv )
        endif
#else
        if (lcplxp .eq. 0) then
          call augq12 ( mode0 , ia , isp , nkaph , iprmb , lmxha , nlmha 
     .      , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s , h )
        else
          call augq2z ( mode0 , mode1 , ia , isp , nkaph , iprmb , lmxha 
     .      , nlmha , kmax , nlma , sv_p_osig ( 3 , ia ) %v , sv_p_oppi( 3 , ia )%v 
     .      , sv_p_osig ( 2 , ia ) %v , sv_p_oppi( 2 , ia )%v , b , ndimh 
     .      , nlmto , s , h , hso )
        endif
#endif
        if(debug) print *,'augg 1111111'
C   --- Add B+ sig B to S and B+ ppi B to H ---
#if MPI
        call augqs3 (kmax, lmxa, nlma, ndimh, isp, g, sv_p_osig( 1, ia )%v, b, s_zv )
        if (mode0 .eq. 0 .and. lcplxp .eq. 0) then
          call augqp3(kmax, nlma, ndimh, isp, g, sv_p_oppi( 1, ia )%v, b, h_zv )
        elseif (mode0 .eq. 0 .and. lcplxp .ne. 0) then
          call augq3z(mode1, kmax, nlma, ndimh, isp, g, sv_p_oppi( 1, ia )%v,b,h_zv,hso_zv)
        endif
#else
        call augqs3( kmax, lmxa, nlma, ndimh, isp, g, sv_p_osig( 1, ia )%v, b, s)
        if (mode0 .eq. 0 .and. lcplxp .eq. 0) then
          call augqp3( kmax, nlma, ndimh, isp, g, sv_p_oppi( 1, ia )%v, b, h )
        elseif (mode0 .eq. 0 .and. lcplxp .ne. 0) then
          call augq3z( mode1, kmax, nlma, ndimh, isp, g, sv_p_oppi( 1, ia )%v, b, h, hso)
        endif
#endif
        if(debug) print *,'augg 222222 ia=',ia
   10   continue
C ... end loop over ia
      enddo
#if MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      allocate(hbuf_zv(ndimh*ndimh))
      call mpi_allreduce ( h_zv , hbuf_zv , 2 * ndimh * ndimh 
     ., mpi_double_precision , mpi_sum , mpi_comm_world , ierr )
      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce h ndimh=%i',' ',256,lgunit(3),
     .  procid,numprocs,ndimh)
      endif
      call daxpy ( 2 * ndimh * ndimh , 1d0 , hbuf_zv , 1 , h , 1 )
      if (allocated(hbuf_zv)) deallocate(hbuf_zv)
      allocate(sbuf_zv(ndimh*ndimh))
      call mpi_allreduce ( s_zv , sbuf_zv , 2 * ndimh * ndimh 
     ., mpi_double_precision , mpi_sum , mpi_comm_world , ierr )
      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce s ndimh=%i',' ',256,lgunit(3),
     .  procid,numprocs,ndimh)
      endif
      call daxpy ( 2 * ndimh * ndimh , 1d0 , sbuf_zv , 1 , s , 1 )
      if (allocated(sbuf_zv)) deallocate(sbuf_zv)
      if (lcplxp .ne. 0) then
        allocate(sbuf_zv(ndimh*ndimh))
        call mpi_allreduce ( hso_zv , sbuf_zv , 2 * ndimh * ndimh 
     .  , mpi_double_precision , mpi_sum , mpi_comm_world , ierr )
        if (mlog) then
          call gettime(datim)
          call awrit3(' augmbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' allreduce hso ndimh=%i',' ',256,lgunit(3),
     .    procid,numprocs,ndimh)
        endif
        call daxpy ( 2 * ndimh * ndimh , 1d0 , sbuf_zv , 1 , hso , 1 )
        if (allocated(sbuf_zv)) deallocate(sbuf_zv)
        if (allocated(hso_zv)) deallocate(hso_zv)
      endif
      if (allocated(s_zv)) deallocate(s_zv)
      if (allocated(h_zv)) deallocate(h_zv)
      deallocate(bproc, stat=ierr)
#endif
      deallocate (b)
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
c
!!takao  NOTE: LzSz is alreay added h by locpot-augmat-gaugm
c
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
      complex(8),allocatable:: tso(:,:,:,:)
      if (mode1 .eq. 1) then
        allocate (tso(ndimh,ndimh,2,2))
        tso=0d0 !call dpzero(tso,ndimh*ndimh*4*2)
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
!! NOTE: LzSz is alreay added h by locpot-augmat-gaugm
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
!! NOTE: LzSz is alreay added h by locpot-augmat-gaugm
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
c        call dpadd(hso(1,1),tso(1,1,1,2),1,2*ndimh*ndimh,1d0)
        hso = hso + tso(:,:,1,2)
        deallocate (tso)
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
c      integer nlmax
c      parameter (nlmax=49)
      integer i1,i2,ilm,k1,k2,l,kjlm !,lla(nlmax) !kjtop,
c      data lla/0,3*1,5*2,7*3,9*4,11*5,13*6/
      integer :: ll
c      if (nlma .gt. nlmax) call rxi('augqs3: increase nlmax to',nlma)
C     call tcn('augqs3')
c      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh
C   ... Make sig*b in g
        do  ilm = 1, nlma
c          l = lla(ilm)
          l = ll(ilm)
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
ctakao
c          csum = (0d0,0d0)
c          do  kjlm = 0, kjtop-1
c            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
c          enddo
c          s(i1,i2) = s(i1,i2) + csum
          s(i1,i2) = s(i1,i2) + sum( dconjg(b(:,:,i1))*g(:,:) )
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
      implicit none
C ... Passed parameters
      integer kmax,nlma,ndimh,isp
      double precision ppi(0:kmax,0:kmax,nlma,nlma,isp)
      double complex b(0:kmax,nlma,ndimh),h(ndimh,ndimh),
     .g(0:kmax,nlma),csum
C ... Local parameters
      integer i1,i2,jlm1,jlm2,k1,k2,kjlm !,kjtop

C     call tcn('augqp3')
c      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh
C   ... g <- ppi*b
c        call dpzero(g,2*kjtop)
        g=0d0
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
cC   ... Make dot products with vectors i1
        do  i1 = 1, i2
cC         csum=zdotc(kjtop,b(0,1,i1),1,g,1)
c          csum = (0d0,0d0)
c          do  kjlm = 0, kjtop-1
c            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
c          enddo
c          h(i1,i2) = h(i1,i2) + csum
          h(i1,i2) = h(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
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
      implicit none
C ... Passed parameters
      integer kmax,nlma,ndimh,isp,mode1
      double complex ppi(0:kmax,0:kmax,nlma,nlma,isp+2*mode1)
      double complex b(0:kmax,nlma,ndimh),h(ndimh,ndimh),
     .g(0:kmax,nlma),csum,gso(0:kmax,nlma),csum1,hso(ndimh,ndimh)
C ... Local parameters
      integer i1,i2,jlm1,jlm2,k1,k2,kjlm !,kjtop

C     call tcn('augqp3')
c      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh
C   ... g <- ppi*b
        g=0d0   !call dpzero(g,2*kjtop)
        gso=0d0 !call dpzero(gso,2*kjtop)
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
cC         csum=zdotc(kjtop,b(0,1,i1),1,g,1)
c          csum = (0d0,0d0)
c          if (mode1 .eq. 1 .and. isp .eq. 2) csum1 = (0d0,0d0)
c          do  kjlm = 0, kjtop-1
c            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
cC     ...  Make LxSx+LySy part of SO
c            if (mode1 .eq. 1 .and. isp .eq. 2)
c     .      csum1 = csum1 + dconjg(b(kjlm,1,i1))*gso(kjlm,1)
c          enddo
c          h(i1,i2) = h(i1,i2) + csum
          h(i1,i2) = h(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
C     ...  Make LxSx+LySy part of SO
          if (mode1 .eq. 1 .and. isp .eq. 2) 
     .      hso(i1,i2) = hso(i1,i2) + .5d0*sum(dconjg(b(:,:,i1))*gso(:,:)) !+ 0.5d0*csum1
        enddo
      enddo
C     call tcx('augqp3')
      end subroutine augq3z


