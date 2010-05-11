Cdelw1 var ['oa', '0', 'oa', 'lda*mda', 'complex(8)', 'zv_w_', 'pzhev']
Cdelw1 var ['ob', '1', 'ob', 'ldb*mdb', 'complex(8)', 'zv_w_', 'pzhev']
Cdelw1 var ['oz', '0', 'oz', 'ldz*mdz', 'complex(8)', 'zv_w_', 'pzhev']
Cdelw1 var ['owork', '0', 'owork', 'lwork', 'complex(8)', 'zv_w_', 'pzhev']
Cdelw1 var ['orwork', '0', 'orwork', 'lrwork', 'real(8)', 'rv_w_', 'pzhev']
Cdelw1 var ['oiwork', '0', 'oiwork', 'liwork', 'integer', 'iv_w_', 'pzhev']
Cdelw1 var ['oifail', '0', 'oifail', 'n', 'integer', 'iv_w_', 'pzhev']
Cdelw1 var ['ot', '0', 'ot', 'n*n', 'complex(8)', 'zv_w_', 'pzhev']
Cdelw1 undel [['call', 'wref'], 'ot', 'pzhev']
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
Cdelw1 do not change  [['call', 'wref'], 'ot', 'pzhev']
      integer ot 
      complex(8):: hamm(*),ovlm(*)
C Pointers to local distributed arrays
Cdelw1       integer oa,ob,oz
       complex(8) ,allocatable :: zv_w_oa(:)
       complex(8) ,allocatable :: zv_w_ob(:)
       complex(8) ,allocatable :: zv_w_oz(:)

C Work arrays
      integer lrwork,lwork,liwork,ifail
Cdelw1       integer owork,orwork,oiwork,oifail
       complex(8) ,allocatable :: zv_w_owork(:)
       real(8) ,allocatable :: rv_w_orwork(:)
       integer ,allocatable :: iv_w_oiwork(:)
       integer ,allocatable :: iv_w_oifail(:)

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
Cdelw1       call defcc(oa, lda*mda)
       allocate(zv_w_oa(lda*mda))
Ctakao_ZeroClear_NotRequiered        if (lda*mda<0) zv_w_oa(:)=0.0d0

Cdelw1       call dstmt(desca,n,hamm,w(oa))
       call dstmt ( desca , n , hamm , zv_w_oa ) 

Cdelw1       call dcopy(2*lda*mda,w(oa),1,hamm,1)
       call dcopy ( 2 * lda * mda , zv_w_oa , 1 , hamm , 1 ) 

Cdelw1 rlse name= oa old_list= oa 
Cdelw1 rlse name= oa new_list= (None)
Cdelw1       call rlse(oa)
       if (allocated(zv_w_oa)) deallocate(zv_w_oa)

      if (lov) then
Cdelw1         call defcc(ob, ldb*mdb)
         allocate(zv_w_ob(ldb*mdb))
Ctakao_ZeroClear_NotRequiered          if (ldb*mdb<0) zv_w_ob(:)=0.0d0

Cdelw1         call dstmt(descb,n,ovlm,w(ob))
         call dstmt ( descb , n , ovlm , zv_w_ob ) 

Cdelw1         call dcopy(2*ldb*mdb,w(ob),1,ovlm,1)
         call dcopy ( 2 * ldb * mdb , zv_w_ob , 1 , ovlm , 1 ) 

Cdelw1 rlse name= ob old_list= ob 
Cdelw1 rlse name= ob new_list= (None)
Cdelw1         call rlse(ob)
         if (allocated(zv_w_ob)) deallocate(zv_w_ob)

      endif
Cdelw1       call defcc(oz, ldz*mdz)
       allocate(zv_w_oz(ldz*mdz))
Ctakao_ZeroClear_NotRequiered        if (ldz*mdz<0) zv_w_oz(:)=0.0d0

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
Cdelw1         call PZHEGVX(ibtype,jobz,range,'U',n,hamm,1,1,desca,ovlm,1,1,
Cdelw1      .               descb,VL,VU,1,IU,abstol,m,nz,e,
Cdelw1      .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
Cdelw1      .               liwork,-1,ifail,iclustr,gap,info)
         call pzhegvx ( ibtype , jobz , range , 'U' , n , hamm , 1 , 1 
     .   , desca , ovlm , 1 , 1 , descb , vl , vu , 1 , iu , abstol , 
     .   m , nz , e , zero , zv_w_oz , 1 , 1 , descz , swork , - 1 , srwork 
     .   , - 1 , liwork , - 1 , ifail , iclustr , gap , info ) 

      else
Cdelw1         call PZHEEVX(jobz,range,'U',n,hamm,1,1,desca,
Cdelw1      .               VL,VU,1,IU,abstol,m,nz,e,
Cdelw1      .               zero,w(oz),1,1,descz,swork,-1,srwork,-1,
Cdelw1      .               liwork,-1,ifail,iclustr,gap,info)
         call pzheevx ( jobz , range , 'U' , n , hamm , 1 , 1 , desca 
     .   , vl , vu , 1 , iu , abstol , m , nz , e , zero , zv_w_oz , 1 
     .   , 1 , descz , swork , - 1 , srwork , - 1 , liwork , - 1 , ifail 
     .   , iclustr , gap , info ) 

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
Cdelw1       call defcc(owork,   lwork)
       allocate(zv_w_owork(lwork))
Ctakao_ZeroClear_NotRequiered        if (lwork<0) zv_w_owork(:)=0.0d0

Cdelw1       call defrr(orwork,  lrwork)
       allocate(rv_w_orwork(lrwork))
Ctakao_ZeroClear_NotRequiered        if (lrwork<0) rv_w_orwork(:)=0.0d0

Cdelw1       call defi (oiwork,  liwork)
       allocate(iv_w_oiwork(liwork))
Ctakao_ZeroClear_NotRequiered        if (liwork<0) iv_w_oiwork(:)=0

Cdelw1       call defi (oifail,  n)
       allocate(iv_w_oifail(n))
Ctakao_ZeroClear_NotRequiered        if (n<0) iv_w_oifail(:)=0

C Diagonalise
      if (lov) then
Cdelw1         call PZHEGVX(ibtype,jobz,range,'U',n,hamm,1,1,desca,ovlm,1,1,
Cdelw1      .               descb,VL,VU,1,IU,abstol,m,nz,e,
Cdelw1      .               zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
Cdelw1      .               lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,info)
         call pzhegvx ( ibtype , jobz , range , 'U' , n , hamm , 1 , 1 
     .   , desca , ovlm , 1 , 1 , descb , vl , vu , 1 , iu , abstol , 
     .   m , nz , e , zero , zv_w_oz , 1 , 1 , descz , zv_w_owork , lwork 
     .   , rv_w_orwork , lrwork , iv_w_oiwork , liwork , iv_w_oifail , 
     .   iclustr , gap , info ) 

      else
Cdelw1         call PZHEEVX(jobz,range,'U',n,hamm,1,1,desca,
Cdelw1      .               VL,VU,1,IU,abstol,m,nz,e,
Cdelw1      .               zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
Cdelw1      .               lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,info)
         call pzheevx ( jobz , range , 'U' , n , hamm , 1 , 1 , desca 
     .   , vl , vu , 1 , iu , abstol , m , nz , e , zero , zv_w_oz , 1 
     .   , 1 , descz , zv_w_owork , lwork , rv_w_orwork , lrwork , iv_w_oiwork 
     .   , liwork , iv_w_oifail , iclustr , gap , info ) 

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
Cdelw1 rlse name= owork old_list= oz owork orwork oiwork oifail 
Cdelw1 rlse name= owork new_list= oz 
Cdelw1       call rlse(owork)
       if (allocated(iv_w_oifail)) deallocate(iv_w_oifail)
       if (allocated(iv_w_oiwork)) deallocate(iv_w_oiwork)
       if (allocated(rv_w_orwork)) deallocate(rv_w_orwork)
       if (allocated(zv_w_owork)) deallocate(zv_w_owork)

Cdelw1 rlse name= orwork old_list= oz 

#error ERROR, try to release name= orwork ,but list does not have  orwork at linenumber= 266 list= oz 

Cdelw1 rlse name= orwork new_list= oz 
Cdelw1       call rlse(orwork)

Cdelw1 rlse name= oiwork old_list= oz 

#error ERROR, try to release name= oiwork ,but list does not have  oiwork at linenumber= 267 list= oz 

Cdelw1 rlse name= oiwork new_list= oz 
Cdelw1       call rlse(oiwork)

Cdelw1 rlse name= oifail old_list= oz 

#error ERROR, try to release name= oifail ,but list does not have  oifail at linenumber= 268 list= oz 

Cdelw1 rlse name= oifail new_list= oz 
Cdelw1       call rlse(oifail)

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
Cdelw1       call dcopy(2*ldz*mdz,w(oz),1,hamm,1)
       call dcopy ( 2 * ldz * mdz , zv_w_oz , 1 , hamm , 1 ) 

      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' poked oz to hamm',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
Cdelw1 rlse name= oz old_list= oz 
Cdelw1 rlse name= oz new_list= (None)
Cdelw1       call rlse(oz)
       if (allocated(zv_w_oz)) deallocate(zv_w_oz)

Cdelw1 do not change ot because of ['call', 'wref']
      call defcc(ot, n*n)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' ready to distribute ot',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
Cdelw1 do not change ot because of ['call', 'wref']
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
Cdelw1 w_varlist remains: ot 
Cdelw1 w_varlistundel: [ot]
Cdelw1 w_varlist (undel), remains: [ot]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate ot because of [call wref]
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine udstmt 

#endif
c###### MPI

