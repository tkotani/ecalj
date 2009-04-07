      subroutine pzhev(lov,n,oh,os,nb,nprow,npcol,emx,nmx,nev,e,ot)
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
C ----------------------------------------------------------------------
#if MPI|MPIK
      implicit none
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 1028)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
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
      integer oh,os,ot
C Pointers to local distributed arrays
      integer oa,ob,oz
C Work arrays
      integer lrwork,lwork,liwork,ifail
      integer owork,orwork,oiwork,oifail
C Work array sizes
      double complex swork(100)
      double precision srwork(100)

C Local
      double precision zero,VL,VU
      parameter (zero = 0d0)
      integer context, i, ibtype, info, m, mycol, myrow,
     .        nprocs, nz, lda, ldb, ldz, mda, mdb, mdz, llda, lldb, lldz
      character jobz, range
      double precision abstol, d1mach
      integer desca(9), descb(9), descz(9), iclustr( MAX_PROCS*2 )
      double precision gap( MAX_PROCS )
      integer IU,izero
      parameter (izero = 0)
      integer lgunit, numroc, iprint
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
C     Initialize the BLACS
      call blacs_pinfo( procid, nprocs )
      if (( nprocs .lt. 1 )) then
         call blacs_setup( procid, nprow*npcol )
      end if
      if (mlog) then
        call gettime(datim)
        call awrit6(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' initialising BLACS; nprow=%i npcol=%i procid=%i nprocs=%i',
     .    ' ',256,lgunit(3),procid,numprocs,nprow,npcol,procid,nprocs)
        call ftflsh(-1)
      endif
C     Initialize a single BLACS context
      call blacs_get( -1, 0, context )
      call blacs_gridinit( context, 'r', nprow, npcol )
      call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
      if (mlog) then
        call gettime(datim)
        call awrit4(' pzhev '//datim//
     .    ' BLACS_GRIDINFO returns nprow=%i npcol=%i myrow=%i mycol=%i',
     .    ' ',256,lgunit(3),nprow,npcol,myrow,mycol)
      endif
C     Bail out if this process is not a part of this context.
      if (myrow .eq. -1) then
        if (mlog) then
          call gettime(datim)
          call awrit2(' pzhev '//datim//' Process %i of %i on '
     .       //shortname(procid)(1:namelen(procid))//
     .       ' is not in context, aborting ..',' ',256,lgunit(3),
     .        procid,numprocs)
        call ftflsh(-1)
        endif
        call gettime(datim)
        if (procid .eq. master) then
          call awrit2(' pzhev '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' is not in context, aborting ..',' ',256,lgunit(1),
     .      procid,numprocs)
        endif
        call fexit(0,0,' ',0)
      endif
C Get dimensions of local matrices a, b and z
      llda = numroc(n,nb,myrow,0,nprow)
      lldb = numroc(n,nb,myrow,0,nprow)
      lldz = numroc(n,nb,myrow,0,nprow)
      mda  = numroc(n,nb,mycol,0,npcol)
      mdb  = numroc(n,nb,mycol,0,npcol)
      mdz  = numroc(n,nb,mycol,0,npcol)
      if (mlog) then
        call gettime(datim)
        call awrit7(' pzhev '//datim//
     .    ' NUMROC nb=%i mda=%i mdb=%i mdz=%i llda=%i lldb=%i lldz=%i',
     .    ' ',256,lgunit(3),nb,mda,mdb,mdz,llda,lldb,lldz)
      endif
C Make basic array descriptors
      call descinit(desca,n,n,nb,nb,izero,izero,context,llda,info)
      if (info .ne. 0) then
        if (iprint() .gt. 0) then
          call awrit1(' pzhev: illegal exit from descinit, info=%i',' ',
     .      256,lgunit(1),info)
        endif
        call fexit(0,0,' ',0)
      else
        if (mlog) then
          call gettime(datim)
          call awrit2(' pzhev '//datim//
     .      ' local array descriptor: %9:1i, context=%i',
     .      ' ',256,lgunit(3),desca,context)
        endif
      endif
      call descinit(descb,n,n,nb,nb,izero,izero,context,lldb,info)
      call descinit(descz,n,n,nb,nb,izero,izero,context,lldz,info)
      lda = desca(9)
      ldb = descb(9)
      ldz = descz(9)
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
      call dstmt(desca,n,w(oh),w(oa))
      call dcopy(2*lda*mda,w(oa),1,w(oh),1)
      call rlse(oa)
      if (lov) then
        call defcc(ob, ldb*mdb)
        call dstmt(descb,n,w(os),w(ob))
        call dcopy(2*ldb*mdb,w(ob),1,w(os),1)
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
        call PZHEGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),1,1,
     .               descb,VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .               liwork,-1,ifail,iclustr,gap,info)
      else
        call PZHEEVX(jobz,range,'U',n,w(oh),1,1,desca,
     .               VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .               liwork,-1,ifail,iclustr,gap,info)
      endif
      lwork = int(swork(1))
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
        call PZHEGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),1,1,
     .               descb,VL,VU,1,IU,abstol,m,nz,e,
     .               zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
     .               lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,info)
      else
        call PZHEEVX(jobz,range,'U',n,w(oh),1,1,desca,
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
C Poke distributed array into t (use heap location oh for temp)
      call dcopy(2*ldz*mdz,w(oz),1,w(oh),1)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' poked oz to oh',' ',256,lgunit(3),
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
      call udstmt(descz,n,w(oh),w(ot))
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
      end
#if MPI|MPIK
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
       end
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
      end
#endif

