c      subroutine pzhev(lov,n,oh,os,nb,nprow,npcol,emx,nmx,nev,e,ot)
      subroutine pzhev(lov,n,hamm,ovlm,nb,nprow,npcol,emx,nmx,nev,e,ot)

C- MPI parallel diagonaliser
C ----------------------------------------------------------------------
Ci Inputs: 
Ci   lov: true if overlap
Ci   n :  dimension
Ci   oh:  pointer to h allocated from heap in calling program
Ci   os:  pointer to s allocated from heap in calling program
Ci   nb,nprow,npcol: BLACS process configuration (defaults if nprow=-1)
Ci   emx,nmx,nev: as usual, see zhev
Co Outputs:
Co   e : eigenvalues
Co   ot: pointer to eigenvectors allocated here
Cr Remarks
Cr   pzhev needs to allocate local arrays from the heap which are
Cr   passed to PZHEGVX in place of h, o, and t. This can be done without
Cr   additional memory as follows. On entry oh and os are allocated but
Cr   not ot. An array oa is allocated and assigned after which it is
Cr   copied back into the heap at address oh. oa is then released and
Cr   ob allocated which is assigned and then copied back at the address
Cr   os. ob is released. Then a local array is allocated at oz; and oh,
Cr   os and oz passed to PZHEGVX. On exit the local arrays at oz have to
Cr   be assembled and returned at the address ot. However we don't need
Cr   oh or os anymore; so the local arrays at oz are copied back to oh
Cr   and oz is released and then ot is allocated. Finally the local
Cr   eigenvector arrays now at oh are distributed into the global
Cr   eigenvector array at ot. 
Cr
Cr   Process configuration: nb is a blocking factor; the processes can
Cr   be built into an nprow X npcol array (as long 
Cr   as nprow*npcol=numprocs). This can speed up PZHEGVX on some
Cr   architectures (see http://www.netlib.org/blacs/BLACS/QRef.html).
Cr   By default, if nprow=-1 on entry pzhev makes a linear array
Cr   (nprow=numprocs, npcol=1) this does no harm on a networked cluster
Cr   but it my be worth tuning for a high performance machine. 
Cu Updates
Cu   25 Apr 04 (K. Beleshenko) workaround scalapack bug
C ----------------------------------------------------------------------
#if MPI & BLACS
      implicit none
      include "mpif.h"
      integer procid, master, numprocs, ierr
C     integer status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 1028)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
C     character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
C     double precision starttime, endtime
      logical mlog,cmdopt
      character*120 strn
C Passed
      logical lov
      integer nmx,nev,n
      double precision emx
C BLACS process configuration
      integer nb,nprow,npcol
C E-vals (output)
      double precision e(n)
C Pointers to H, S (input) and Z (output)
c      integer oh,os,ot !oh --> hamm takao
      integer ot 
      complex(8):: hamm(*),ovlm(*)
C Pointers to local distributed arrays
      integer oa,ob,oz
C Work arrays
      integer lrwork,lwork,liwork,ifail
      integer owork,orwork,oiwork,oifail
C Work array sizes
      double complex swork
      double precision srwork(3)

C Local
      double precision zero,VL,VU
      parameter (zero = 0d0)
      integer BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     .        MB_, NB_, RSRC_, CSRC_, LLD_
      parameter ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     .          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     .          RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      integer context, i, iam, ibtype, info, m, mycol, myrow,
     .        nprocs, nz,
     .        NB_A, NB_B, NB_Z, CSRC_A, CSRC_B, CSRC_Z,
     .        lda, ldb, ldz, mda, mdb, mdz
      character jobz, range
      double precision abstol, d1mach
      integer desca(DLEN_), descb(DLEN_), descz(DLEN_),
     .        iclustr( MAX_PROCS*2 )
      double precision gap( MAX_PROCS )
      integer IU
      integer lgunit, numroc
      external blacs_exit, blacs_get, blacs_gridexit,
     .         blacs_gridinfo, blacs_gridinit, blacs_pinfo,
     .         blacs_setup, descinit, pzhegvx, pzlaprnt
C ... Heap
      integer w(1)
      common /w/ w

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
C MPI process configuration
      if (nprow .eq. -1) then
        nprow = numprocs
        npcol = 1
      endif
C     Initialize the BLACS
      call blacs_pinfo( iam, nprocs )
      if( ( nprocs.lt.1 ) ) then
         call blacs_setup( iam, nprow*npcol )
      end if
      if (mlog) then
        call gettime(datim)
        call awrit6(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' initialising BLACS; nprow=%i npcol=%i iam=%i nprocs=%i',
     .    ' ',256,lgunit(3),procid,numprocs,nprow,npcol,iam,nprocs)
        call ftflsh(-1)
      endif
C     Initialize a single BLACS context
      call blacs_get( -1, 0, context )
      call blacs_gridinit( context, 'r', nprow, npcol )
      call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
C     Bail out if this process is not a part of this context.
      if ( myrow .eq. -1 ) then
        if (mlog) then
          call gettime(datim)
          call awrit2(' pzhev '//datim//' Process %i of %i on '
     .       //shortname(procid)(1:namelen(procid))//
     .       ' is not in context, aborting ..',' ',256,lgunit(3),
     .        procid,numprocs)
        call ftflsh(-1)
        endif
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' is not in context, aborting ..',' ',256,lgunit(1),
     .    procid,numprocs)
        call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
        call fexit(0,0,' ',0)
      endif
C     These are basic array descriptors
      call descinit( desca, n, n, nb, nb, 0, 0, context, n, info )
      call descinit( descz, n, n, nb, nb, 0, 0, context, n, info )
      call descinit( descb, n, n, nb, nb, 0, 0, context, n, info )
C Get dimensions of local matrices a, b and z
      lda = desca( LLD_ )
      ldb = descb( LLD_ )
      ldz = descz( LLD_ )
      NB_A = desca ( NB_ )
      NB_B = descb ( NB_ )
      NB_Z = descz ( NB_ )
      CSRC_A = desca ( CSRC_ )
      CSRC_B = descb ( CSRC_ )
      CSRC_Z = descz ( CSRC_ )
      mda = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL )
      mdb = NUMROC( N, NB_B, MYCOL, CSRC_B, NPCOL )
      mdz = NUMROC( N, NB_Z, MYCOL, CSRC_Z, NPCOL )
      if (mlog) then
        call gettime(datim)
        if (lov) then
          call awrit7(' pzhev '//datim//
     .                ' getting local matrix dimensions:%N'//
     .                '   n=%i '//
     .                ' a:(%ix%i)'//
     .                ' b:(%ix%i)'//
     .                ' z:(%ix%i)'//
     .                ' allocating from heap ..',' ',256,lgunit(3),
     .                n,lda,mda,ldb,mdb,ldz,mdz)
        else
          call awrit5(' pzhev '//datim//
     .                ' getting local matrix dimensions:%N'//
     .                '   n=%i '//
     .                ' a:(%ix%i)'//
     .                ' z:(%ix%i)'//
     .                ' allocating from heap ..',' ',256,lgunit(3),
     .                n,lda,mda,ldz,mdz)
        endif
        call ftflsh(-1)
      endif
C Distribute h and s into local arrays
      call defcc(oa, lda*mda)
      call dstmt(desca,n,hamm,w(oa))
      call dcopy(2*lda*mda,w(oa),1,hamm,1)
      call rlse(oa)
      if (lov) then
        call defcc(ob, ldb*mdb)
        call dstmt(descb,n,ovlm,w(ob))
        call dcopy(2*ldb*mdb,w(ob),1,ovlm,1)
        call rlse(ob)
      endif
      call defcc(oz, ldz*mdz)
      ibtype = 1
      if (nmx .le. 0) then
        jobz = 'N'
        range = 'V'
        VL = -1d12
        VU = emx
        IU = n
      else
        jobz = 'V'
        range = 'I'
        VL = zero
        VU = zero
        IU = min(n,nmx)
      endif
      abstol = d1mach(3)
C Workspace query
      if (lov) then
        call PZHEGVX(ibtype,jobz,range,'U',n,hamm,1,1,desca,ovlm,1,1,
     .               descb,VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .               liwork,-1,ifail,iclustr,gap,info)
      else
        call PZHEEVX(jobz,range,'U',n,hamm,1,1,desca,
     .               VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .               liwork,-1,ifail,iclustr,gap,info)
      endif
      lwork = int(swork)
      lrwork = int(srwork(1))
      if (mlog) then
        call gettime(datim)
        call awrit3(' pzhev '//datim//' Optimal scalapack worksizes:'//
     .              '%N   lwork=%i lrwork=%i liwork=%i. '//
     .              ' Allocating from heap ..',' ',256,lgunit(3),
     .              lwork,lrwork,liwork)
        call ftflsh(-1)
      endif
      call defcc(owork,   lwork)
      call defrr(orwork,  lrwork)
      call defi (oiwork,  liwork)
      call defi (oifail,  n)
C Diagonalise
      if (lov) then
        call PZHEGVX(ibtype,jobz,range,'U',n,hamm,1,1,desca,ovlm,1,1,
     .               descb,VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
     .               lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,info)
      else
        call PZHEEVX(jobz,range,'U',n,hamm,1,1,desca,
     .               VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
     .               lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,info)
      endif
      if (info .ne. 0 .and. procid .eq. master) then
        if (lov) then
          call awrit1(' **** in pzhev, PZHEGVX returned info=%i',' ',128
     .      ,lgunit(1),info)
        else
          call awrit1(' **** in pzhev, PZHEEVX returned info=%i',' ',128
     .      ,lgunit(1),info)
        endif
      endif
      nev = nz
      call rlse(owork)
      call rlse(orwork)
      call rlse(oiwork)
      call rlse(oifail)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' is at the barrier',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C Poke distributed array into t (use heap location hamm for temp)
      call dcopy(2*ldz*mdz,w(oz),1,hamm,1)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' poked oz to hamm',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      call rlse(oz)
      call defcc(ot, n*n)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' ready to distribute ot',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      call udstmt(descz,n,hamm,w(ot))
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' done distribute ot',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
C Don't do this!
C      call blacs_gridexit(context)
C      call blacs_exit(1)
#if SUN
      call ieee_flags( 'clear', 'exception', 'underflow', '')
#endif
#endif
c##### MPI
Cgetarg       end
       end subroutine pzhev 

#if MPI & BLACS 
      subroutine dstmt(desc,n,ag,al)

C Distribute global matrix ag into local matrix al
      implicit none
      integer desc(1),n
      double complex ag(n,n),al(n,n)
      integer i,j
      do  i = 1, n
        do  j = 1, n
           CALL PZELSET( al, i, j, desc, ag(i,j))
         enddo
       enddo
Cgetarg        end
        end subroutine dstmt 

      subroutine udstmt(desc,n,al,ag)

C Undistribute local matrix al into global matrix ag
      implicit none
      integer desc(1),n
      double complex al(n,n),ag(n,n)

      integer i,j
      double complex alpha

      do  i = 1, n
        do j = 1, n
          call PZELGET( 'A', ' ', alpha, al, i, j, desc)
          ag(i,j) = alpha
        enddo
      enddo
Cgetarg       end
       end subroutine udstmt 

#endif
c###### MPI

