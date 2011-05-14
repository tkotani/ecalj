CSFCPP#define F90 1
CSFCPP#ifndef F90
CSFCPP      subroutine yysbnv(mode,apr,api,ldap,ija,offs,nlev,
CSFCPP     .nn1,nn2,w,ldw,ar,ai,lda,ierr)
CSFCPP#else
      recursive subroutine yysbnv(mode,apr,api,ldap,ija,offs,nlev,
     .nn1,nn2,w,ldw,ar,ai,lda,ierr)
CSFCPP#endif
C- Inversion of a complex block sparse matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode   :1s digit
Ci           1 only the diagonal parts subblocks of matrix are complex
Ci             (only case tested so far)
Ci           2 full matrix is complex
Ci           3 matrix is hermitian (not tested)
Ci          10s digit fixes how array is partitioned
Ci             It is passed as first argument to psybnv
Ci           0 partitions array into approximately equal numbers
Ci             of nonzero subblocks
Ci           1 partitions array into approximately equal dimensions
Ci   nnmin  :fewest subblocks for which yysbnv calls itself recursively.
Ci   ar,ai  :real,imaginary parts of matrix inverse
Ci   lda    :leading dimension of ar,ai
Ci   apr,api:real,imaginary parts of matrix subblocks; see Remarks
Ci   ldap   :leading dimensions of apr,api
Ci   ija    :column index packed array pointer data; see Remarks
Ci   offs   :offsets to first entries in matrix subblocks. NB: offs(1)=0
Ci           Subblock dimension of row(or col) i = offs(i+1) - offs(i)
Ci   nn1,nn2:range of subblocks which comprise the matrix to be inverted
Ci           the matrix subblock to be inverted consist of the
Ci           rows and columns  offs(nn1)+1...offs(nn2+1)
Ci   nlev   :the maximum number of recursion levels allowed.
Ci           To avoid roundoff errors, nlev=2 is suggested.
Ci   w,ldw  :complex work array of dimension ldw*n
Co Outputs:
Co   ar,ai  :inverse of matrix, stored in a(i:j,i:j) with
Co           i = offs(nn1)+1 and j = offs(nn2+1).
Co   ierr is returned nonzero if matrix could not be inverted.
Cr Remarks:
Cr  *yysbnv uses the Strassen algorithm to invert a subblock of a matrix
Cr   stored in block packed form.  The inverse is not assumed to be
Cr   sparse, and is returned in conventional form in ar,ai.  Arrays
Cr   ija,offs,apr,api contain data for the matrix; see yysp2a for a
Cr   description of block matrix storage conventions, the use of these
Cr   arrays, and how a matrix subblock may be unpacked from them into
Cr   conventional form.
Cr
Cr  *The matrix to be inverted comprises the rows (and columns)
Cr   offs(nn1)+1...offs(nn2+1) of a.  ysqnv partitions these rows and
Cr   columns into four subblocks a11,a21,a12,a22.  (See below for how
Cr   the 1 and 2 partitions are apportioned.)  Let c be the inverse,
Cr   with subblocks c11,c21,c12,c22.  Then inversion proceeds by:
Cr      (1) invert   a22 (see below for recursive inversion)
Cr      (2) generate (c11)^-1 = (a11 - a12 a22^-1 a21)
Cr      (3) invert   (c11)^-1
Cr      (4) generate c21 = -(a22^-1 a21) c11 is
Cr      (5) generate c12 = -c11 a12 a22^-1 is generated
Cr      (6) generate c22 = a22^-1 + (a22^-1 a21 c11) (a12 a22^-1)
Cr   These steps require two inversions, three sparse and three normal
Cr   multiplications.
Cr
Cr  *Partitioning into subblocks.  Partitions 1 and 2 consist of
Cr   offs(nn1)+1...offs(nm1+1) and offs(nm1+1)+1...offs(nn2+1).
Cr   nm1 is set in psybnv according to the 10s digit of mode.
Cr
Cr  *Inversion of the subblocks.  It is more efficient to allow
Cr   inversion of a22 to proceed recursively, if your compiler allows
Cr   it.  Recursion proceeds provided nlev>0, the number of subblocks
Cr   exceeds nnmin and the dimension of the matrix to be inverted
Cr   exceeds nmin.  (c11)^-1 is inverted calling yyqinv, which can
Cr   also proceed recursively.
Cb Bugs
Cb   yysbnv has not been tested for the hermitian case, 
Cb   nor handle it efficiently.
C ----------------------------------------------------------------
C     implicit none
      integer mode,nn1,nn2,lda,ldw,ierr,nlev,ldap,ija(2,*),offs(nn2)
      double precision ar(lda,*),ai(lda,*),w(ldw,*),
     .apr(ldap,ldap,*),api(ldap,ldap,*)
      character*1 cs
      double precision xx
      integer nmin,n1,n2,nm1,nm2,m1,m2,i,j,nd,modl,mode0,mode1,nnmin
      logical:: isanrg, l_dummy_isanrg
#if DEBUG
      character*10 fmt
      data fmt /'(9f10.5)'/
#endif

      if (nn1 .gt. nn2) return
      cs = ' '
      if (mod(mode,10) .eq. 3) cs = 'h'
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
Ckino isanrg is logical function,       call isanrg(mode1,0,1,'yysbnv:','10s digit mode',.true.)
      l_dummy_isanrg=isanrg(mode1,0,1,'yysbnv:','10s digit mode',.true.)
      modl = 10*mod(mode0,3)
C     nmin = smallest rank below which standard inversion is used
      nmin = 32

C --- Partition a into four subblocks according to 10s digit mode ---
      call psybnv(mode1,nn1,nn2,ija,offs,nmin,nm1,nm2,n1,n2,m1,m2,nd,
     .nnmin)

#if DEBUG
      print 334, 'entering yysbnv ',mode,nlev,nn1,nm1,nn2,m1-n1,m2
  334 format(1x,a,': mode=',i2,': nlev=',i1,
     .'  nn1=',i2,'  nm1=',i2,'  nn2=',i2,
     .'  m1=',i3,'  m2=',i3)
#endif

C --- Straight yygefa,di if n2 lt nmin ---
      if (nd .le. nmin) then
        call yysp2a(mode0,nn1,nn2,nn1,nn2,apr,api,ldap,ija,offs,
     .  ar,ai,lda)
        if (cs .ne. 'h' .or. nd .eq. 1) then
          call yygefa(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,nd,w,ierr)
          if (ierr .ne. 0) return
          if (nd .eq. 1) then
            call yygedi(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,nd,w,xx,
     .      w(1,2),xx,1)
          else
            call yygedi(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,nd,w,xx,
     .      w(1,2),w(1,3),1)
          endif
        else
          call yyhifa(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,nd,w,ierr)
          if (ierr .ne. 0) return
          call yyhidi(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,nd,w,xx,i,
     .    w(1,2),w(1,3),1)
          do  10  i = n1+1, n2
            do  10  j = n1+1, i
              ar(i,j) =  ar(j,i)
   10     ai(i,j) = -ai(j,i)
        endif
#if DEBUG
        print *, 'exit yysbnv, used standard LU decomposition'
#endif
        return
      endif

      if (lda .lt. n2) call rx('yysbnv: lda lt n2')

C --- R1 = a22^-1 in a22 ---
      if (m2 .lt. nmin) then
        call yysp2a(mode0,1+nm1,nn2,1+nm1,nn2,apr,api,ldap,ija,offs,
     .  ar,ai,lda)
        if (cs .eq. 'h') then
          call yyhifa(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w,ierr)
          if (ierr .ne. 0) return
          call yyhidi(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w,
     .    xx,i,w(1,2),w(1,3),1)
          do  12  i = m1+1, n2
            do  12  j = m1+1, i
              ar(i,j) =  ar(j,i)
   12     ai(i,j) = -ai(j,i)
        else
          call yygefa(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w,ierr)
          if (ierr .ne. 0) return
          call yygedi(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w,xx,
     .    w(1,2),w(1,3),1)
        endif
CSFCPP#if RECUR | F90
      elseif (nm2 .gt. nnmin .and. nlev .gt. 0) then
        call yysbnv(mode,apr,api,ldap,ija,offs,nlev-1,
     .  nm1+1,(nn2),w,ldw,ar,ai,lda,ierr)
#if AIX | SUN_ULTRA 
        mode1 = mod(mode/10,10)
        call psybnv(mode1,nn1,nn2,ija,offs,nmin,nm1,nm2,n1,n2,m1,m2,nd,
     .  nnmin)
#endif
#if DEBUG
        print 333, 'recursive call exit',nlev,n2-n1,m1-n1,m2
#endif
CSFCPP#endif
      else
        call yysp2a(mode0,1+nm1,nn2,1+nm1,nn2,apr,api,ldap,ija,offs,
     .  ar,ai,lda)
        call yyqinv(cs,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,
     .  max(nlev-1,0),m2,w,ldw,ierr)
      endif
      if (ierr .ne. 0) return
*     call yprm(.false.,'R1',2,6,fmt,ar(1+m1,1+m1),lda,m2,lda,m2)

C --- R3 = a22^-1 a21 in w(2); uses w(1+m1:n2,1+n1+m1:n1+2*m1) ---
C     NB: R3 must be preserved to make c21
C      call yygemm('N','N',m2,m1-n1,m2,1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),
C     .  lda,ar(1+m1,1+n1),ai(1+m1,1+n1),lda,0d0,
C     .  w(1+m1,1+n1),w(1+m1,1+n1+m1),ldw)
      call yymsbm(modl,nm1+1,nn2,nn1,nm1,m2,apr,api,ldap,ija,offs,
     .ar(1+m1,1),ai(1+m1,1),lda,w(1+m1,1),w(1+m1,1+m1),ldw)
*     call yprm(.false.,'R3',2,6,fmt,w(1+m1,1),ldw,m2,m1,m1)

C --- R5 = (a11 - a12 a22^-1 a21) in a11 ---
C      call yygemm('N','N',m1-n1,m1-n1,m2,-1d0,ar(1+n1,1+m1),
C     .  ai(1+n1,1+m1),lda,w(1+m1,1+n1),w(1+m1,1+n1+m1),ldw,1d0,
C     .  ar(1+n1,1+n1),ai(1+n1,1+n1),lda)
*      call yprm(.false.,'R5',2,6,fmt,ar,lda,m1,lda,m1)
      call yysp2a(mode0,nn1,nm1,nn1,nm1,apr,api,ldap,ija,offs,ar,ai,lda)
      call yysbmm(modl+2,nn1,nm1,nm1+1,nn2,m1-n1,apr,api,ldap,ija,offs,
     .w(1,1+n1),w(1,1+n1+m1),ldw,ar(1,1+n1),ai(1,1+n1),lda)
*     call yprm(.false.,'R5',2,6,fmt,ar,lda,m1,lda,m1)

C --- c11 = -R6 = -R5^-1 = (a11 - a12 a22^-1 a21)^-1 in a11 ---
      if (m1-n1 .lt. nmin) then
        if (cs .eq. 'h') then
          call yyhifa(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,m1-n1,w,ierr)
          if (ierr .ne. 0) return
          call yyhidi(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,m1-n1,w,xx,i,
     .    w(1,2),w(1,3),1)
          do  16  i = 1, n2
            do  16  j = 1, i
              ar(i,j) =  ar(j,i)
   16     ai(i,j) = -ai(j,i)
        else
          call yygefa(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,m1-n1,w,ierr)
          if (ierr .ne. 0) return
          call yygedi(ar(1+n1,1+n1),ai(1+n1,1+n1),lda,m1-n1,w,xx,w(1,2),
     .    w(1,3),1)
        endif
      else
        call yyqinv(cs,ar(1+n1,1+n1),ai(1+n1,1+n1),lda,max(nlev-1,0),
     .  m1-n1,w,ldw,ierr)
      endif
      if (ierr .ne. 0) return
*     call yprm(.false.,'-R6',2,6,fmt,ar,lda,m1,lda,m1)

C --- c21 = R3 * R6 = -a22^-1 a21 c11 in a21 ---
      call yygemm('N','N',m2,m1-n1,m1-n1,-1d0,w(1+m1,1+n1),
     .w(1+m1,1+n1+m1),ldw,ar(1+n1,1+n1),ai(1+n1,1+n1),lda,0d0,
     .ar(1+m1,1+n1),ai(1+m1,1+n1),lda)
*     call yprm(.false.,'c21',2,6,fmt,ar(1+m1,1),lda,m2,lda,m1)

C ... Possibly add a special branch for the hermitian case
      if (.false.) then
      else

C   --- R2 = a12 * a22^-1 in w(1) ---
C       call yygemm('N','N',m1-n1,m2,m2,1d0,ar(1+n1,1+m1),
C    .    ai(1+n1,1+m1),lda,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,0d0,
C    .    w(1+n1,1),w(1+n1,1+m2),ldw)
        call yysbmm(modl,nn1,nm1,nm1+1,nn2,m2,apr,api,ldap,ija,offs,
     .  ar(1,1+m1),ai(1,1+m1),lda,w(1,1),w(1,1+m2),ldw)
*       call yprm(.false.,'R2',2,6,fmt,w,ldw,m1,m2,m2)

C   --- c12 = R6 * R2 = -c11 a12 a22^-1 ---
        call yygemm('N','N',m1-n1,m2,m1-n1,-1d0,
     .  ar(1+n1,1+n1),ai(1+n1,1+n1),lda,w(1+n1,1),w(1+n1,1+m2),ldw,
     .  0d0,ar(1+n1,1+m1),ai(1+n1,1+m1),lda)
*       call yprm(.false.,'c12',2,6,fmt,ar(1,1+m1),lda,m1,lda,m2)

C   --- c22 = R1 - c21 * R2 = a22^-1 + a22^-1 a21 c11 a12 a22^-1 ---
        call yygemm('N','N',m2,m2,m1-n1,-1d0,ar(1+m1,1+n1),ai(1+m1,1+n1)
     .  ,lda,w(1+n1,1),w(1+n1,1+m2),ldw,1d0,ar(1+m1,1+m1),ai(1+m1,1+
     .  m1),lda)
      endif

#if DEBUG
C     call yprm(.false.,'a^-1',2,6,fmt,ar(1+n1,1+n1),lda,nd,lda,nd)
      print 333, 'exiting  yysbnv',nlev,n2-n1,m1-n1,m2
  333 format(1x,a,': nlev=',i1,'  ndim=',i3,' partitioned into',2i3)
#endif
      end
      subroutine psybnv(mode,nn1,nn2,ija,offs,nmin,nm1,nm2,n1,n2,m1,m2,
     .nd,nnmin)
C- Return parameters that partition block sparse array into two parts
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode:fixes how the array is partitioned
Ci        mode=0 partitions array into approximately equal numbers
Ci               of nonzero subblocks
Ci        mode=1 partitions array into approximately equal dimensions
Ci        mode>1 uses input nm1
Ci   nn1 :range of subblocks which comprise the matrix to be inverted
Cr   nn2  the matrix subblock to be inverted consist of the
Cr        rows and columns  offs(nn1)+1...offs(nn2+1)
Ci   ija :column index packed array pointer data
Ci   offs:table of offsets to array subblocks
Ci   nmin:smallest rank below which standard inversion is used
Cio Inputs/Outputs
Cio  nm1 : number of blocks to lower partition of array; see Remarks
Cio        How nm1 is set depends on mode.
Co Outputs
Co   nm1 : number of blocks to lower partition of array
Co   nm2 : number of blocks to upper partition of array
Co   n1  : offset to first column in array - 1
Co   n2  : offset to last column in array
Co   m1  : dimension of lower partition + n1
Co   m2  : dimension of upper partition
Co   nd  : dimension of the matrix to be inverted
Co  nnmin:fewest subblocks for which yysbnv calls itself recursively.
Cr Remarks
Cr   psybnv supplies parameters partitioning a array in block sparse
Cr   format in two.  The array rows (or columns) may be described in
Cr   terms of number of subblocks, with nn1 and nn2 input:
Cr
Cr            (1  nn1     nm1     nn2 )
Cr                        <-- nm2 -->
Cr
Cr   nm1 may be calculated internally or is input, depending on mode.
Cr   In terms of true array dimensions, the parameters are:
Cr
Cr            (1  n1      m1      n2 )
Cr                <------ nd ------>
Cr                        <-- m2 -->
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nn1,nn2,ija(2,*),offs(nn2),nmin
      integer n1,n2,nm1,nm2,m1,m2,nd,nnmin
      integer iblk,nblk,ir,k,ic

      if (nn2 .gt. nn1) then
        if (mode .eq. 0) then
C   --- Partition into approx. equal nonzero subblocks ---
          nblk = 0
          do  10  ir = nn1, nn2
            nblk = nblk + 1
            do  11  k = ija(1,ir), ija(1,ir+1)-1
              ic  = ija(1,k)
              if (ic .lt. nn1 .or. ic .gt. nn2) goto 11
              nblk = nblk+1
   11       continue
   10     continue
C ...   Count to half the total number of subblocks
          iblk = 0
          do  20  ir = nn1, nn2
            iblk = iblk + 1
            do  21  k = ija(1,ir), ija(1,ir+1)-1
              ic  = ija(1,k)
              if (ic .lt. nn1 .or. ic .gt. nn2) goto 21
              iblk = iblk+1
   21       continue
            if (2*iblk .gt. nblk) then
              nm1 = ir-1
              goto 22
            endif
   20     continue
   22     continue

C --- Partition into approximately equal dimensions ---
        elseif (mode .eq. 1) then
          nd = (offs(nn2+1) + offs(nn1))/2
          nm1 = min((nn1+nn2)/2+1,nn2)
   31     continue
          if (offs(nm1) .lt. nd) then
            nm1 = nm1+1
            goto 31
          endif
   32     continue
          if (iabs(offs(nm1)-nd) .gt. iabs(offs(nm1-1)-nd)) then
            nm1 = nm1-1
            if (nm1 .gt. 1) goto 32
          endif
          nm1 = nm1-1
        endif
      endif
      nm1 = max(min(nm1,nn2-1),nn1)
      nm2 = nn2-nm1
      n1 = offs(nn1)
      n2 = offs(nn2+1)
      nd = n2-n1
      m1 = offs(nm1+1)
      m2 = n2-m1
C ... Find some reasonable nnmin
      nnmin = 0
   40 nnmin = nnmin+1
      if (nnmin .lt. nn2) then
        if (offs(nnmin+1) .lt. nmin) goto 40
      endif
      nnmin = nnmin-1
      end

