      subroutine yysbib(cs,mode,apr,api,ldap,ija,offs,nlev,
     .nn2,nb,w,ldw,w2,ar,ai,lda,br,bi,ldb,ierr)
C- Solve a x = b, with a = complex block sparse matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:  : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          'h'  a is assumed hermitian.
Ci          'b'  Assume partial inverse for a is already performed
Ci               and solve a x = b for a new b
Ci               ar,ai must be preserved between successive calls.
Ci          '1'  Only return b1 (solution of ax=b for lower partition)
Ci               Useful if only rhs for this lower block is needed.
Ci          '2'  Skip final matrix multiplication to generate x2
Ci               (see Remarks).  Useful if only some portions of
Ci               ax=b for upper partition are needed.
Ci          '4'  Do multiplications using standard four real operations
Ci               (slower, but avoids additions and subtractions that
Ci                can reduce machine precision)
Ci   mode   :1s digit
Ci           1  only the diagonal parts subblocks of matrix are complex
Ci              (only case tested so far)
Ci           2  full matrix is complex
Ci   ar,ai  :real,imaginary parts of matrix inverse
Ci   lda    :leading dimension of ar,ai
Ci   apr,api:real,imaginary parts of matrix subblocks; see Remarks
Ci   ldap   :leading dimensions of apr,api
Ci   ija    :column index packed array pointer data; see Remarks
Ci   offs   :offsets to first entries in matrix subblocks. NB: offs(1)=0
Ci           Subblock dimension of row(or col) i = offs(i+1) - offs(i)
Ci   nn1,nn2:range of subblocks which comprise the matrix to be inverted
Cr           the matrix subblock to be inverted consist of the
Cr           rows and columns  offs(nn1)+1...offs(nn2+1)
Ci   nlev   :the maximum number of recursion levels allowed in the
Ci           matrix inversion steps; see Remarks
Ci   w,ldw  :complex work array of dimension ldw*nd, nd= rank of a
Ci   w2     :a complex work array of dimension nb*nd
Ci           NB: w and w2 may use the same address space
Ci   br,bi  :real and imaginary parts of rhs of equation a x = b
Ci   ldb    :leading dimension of br,bi
Co Outputs:
Co   ar,ai  :partial matrix inverse; see Remarks
Co   ierr   :returned nonzero if matrix could not be inverted.
Co   br,bi  :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
Co           NB: yysbib may return before completing generation of
Co           a^-1 b, depending on input cs; see Remarks.
Co   w2     :b2 - a21 x1; see Remarks
Cr Remarks:
Cr  *yysbib uses a block decomposition to solve the linear system
Cr   a x = b.  a is partitioned into subblocks a11,a21,a12,a22 and
Cr   is stored in a block sparse format.  The partial inverse
Cr   is returned in in ar,ai, so that subsequent calls to yysbib with
Cr   new br,bi may be repeated without another decomposition.  By
Cr   setting switches via input variable cs, yysbib may return at
Cr   earlier stages than the full solution of a x = b.
Cr
Cr  *ija,offs,apr,api contain data for the matrix; see yysp2a for a
Cr   description of block matrix storage conventions, the use of these
Cr   arrays, and how a matrix subblock may be unpacked from them into
Cr   conventional form.
Cr
Cr  *The matrix to be inverted comprises the rows (and columns)
Cr   offs(nn1)+1...offs(nn2+1) of a.  ysqnvb partitions these rows and
Cr   columns into four subblocks a11,a21,a12,a22.  (See below for how
Cr   the 1 and 2 partitions are apportioned.)  Let c be the inverse,
Cr   with subblocks c11,c21,c12,c22.  yysbib follows these steps:
Cr
Cr     (1) sparse block inverter ysbnv generates a22^-1
Cr     (2) sparse multiplications make (c11)^-1 = (a11 - a12 a22^-1 a21)
Cr     (3) normal multiplication generates (b1 - a12 a22^-1 b2)
Cr     (4) (c11^-1) x1 = (b1 - a12 a22^-1 b2) is solved for x1
Cr         NB: ysqnvb returns at this point if cs contains '1'
Cr     (5) sparse multiplication generates (b2 - a21 x1)
Cr         NB: ysqnvb returns at this point if cs contains '2'
Cr     (6) normal multiplication generates x2 = a22^-1 (b2 - a21 x1)
Cr
Cr  *Partitioning into subblocks.  Partitions 1 and 2 consist of
Cr   offs(nn1)+1...offs(nm1+1) and offs(nm1+1)+1...offs(nn2+1).
Cr   nm1 is set in psybnv according to the 10s digit of mode.
Cr
Cr  *It is more efficient to allow inversion to proceed recursively,
Cr   if your compiler allows it.  (Steps (1) and (4) involve matrix
Cr   inversions.  Recursion proceeds provided nlev>0, the number of
Cr   subblocks exceeds nnmin and the dimension of the matrix to be
Cr   inverted exceeds nmin.  The latter should be chosen where the
Cr   tradoff between the extra overhead and fewer O(N^3) operations
Cr   takes place.  To avoid roundoff errors, nlev<=2 is suggested.
Cr
Cb Bugs
Cb   yysbib has not been tested for the hermitian case,
Cb   nor handle it efficiently.
C ----------------------------------------------------------------
C     implicit none
      integer mode,nn2,nb,lda,ldb,ldw,ierr,nlev,ldap,
     .ija(2,*),offs(nn2)
      double precision ar(lda,*),ai(lda,*),br(ldb,1),bi(ldb,1),
     .w(ldw,*),apr(ldap,ldap,*),api(ldap,ldap,*),w2(*)
      character cs*(*)
      integer nmin,n1,n2,nm1,nm2,m1,m2,i,j,nd,modl,mode0,mode1,jr,ji,nnm
      logical lsx,lsx1,lsx2,ltrns
      character cm*2
#if DEBUG
      character*10 fmt
      data fmt /'(9f12.6)'/
#endif

C --- Setup ---
      ierr = 0
      ltrns = .false.
      lsx = .false.
      lsx1 = .false.
      lsx2 = .false.

      cm = 'N '
      j = len(cs)
      do  2  i = 1, j
        if (cs(i:i) .eq. 't') then
          ltrns = .true.
        elseif (cs(i:i) .eq. 'b') then
          lsx = .true.
        elseif (cs(i:i) .eq. '4') then
          cm = 'N4'
        elseif (cs(i:i) .eq. '1') then
          lsx1 = .true.
        elseif (cs(i:i) .eq. '2') then
          lsx2 = .true.
        elseif (cs(i:i) .eq. 'h') then
        elseif (cs(i:i) .ne. ' ') then
          call rxs2('yysbib: unrecognized switch, cs=''',cs,'''')
        endif
    2 continue
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      modl = 10*mod(mode0,3)
C     Smallest rank below which standard inversion is used
      nmin = 32

C --- Partition a into four subblocks according to 10s digit mode ---
      call psybnv(mode1,1,nn2,ija,offs,nmin,nm1,nm2,n1,n2,m1,m2,nd,nnm)
#if DEBUG
      print 334, 'entering yysbib',mode,nlev,nm1,nn2,m1,m2
  334 format(1x,a,': mode=',i2,': nlev=',i1,
     .'  nm1=',i2,'  nn2=',i2,'  m1=',i3,'  m2=',i3)
#endif

C --- Straight inversion if n2 lt nmin ---
      if (nd .le. nmin) then
        if (.not. lsx) call yysp2a(mode0,1,nn2,1,nn2,apr,api,ldap,ija,
     .  offs,ar,ai,lda)
        if (nb .eq. 0) return
        call yyqnvb(cs,ar,ai,lda,nd,nb,w,ldw,w2,br,bi,ldb,ierr)
        return
      endif

      if (lda .lt. n2) call rx('yysbib: lda lt n')

C --- Generate a22^-1 and (a11 - a12 a22^-1 a21) ---
      if (.not. lsx) then

C --- R1 = a22^-1 in a22 ---
        call yysbnv(mode,apr,api,ldap,ija,offs,max(nlev-1,0),
     .  nm1+1,nn2,w,ldw,ar,ai,lda,ierr)
        if (ierr .ne. 0) return
C       call yprm(.false.,'R1',2,6,fmt,ar(1+m1,1+m1),lda,m2,lda,m2)

        if (.not. ltrns) then

C   --- a12 a22^-1 in w (uses w(1:m1+m2,1:m1)); copy back to a12 ---
C       call yygemm('N','N',m1,m2,m2,1d0,ar(1,1+m1),ai(1,1+m1),lda,
C    .    ar(1+m1,1+m1),ai(1+m1,1+m1),lda,0d0,w,w(1,1+m2),ldw)
          call yysbmm(modl,1,nm1,nm1+1,nn2,m2,apr,api,ldap,ija,offs,
     .    ar(1,1+m1),ai(1,1+m1),lda,w(1,1),w(1,1+m2),ldw)
          do  20  j = 1, m2
            do  20  i = 1, m1
              ar(i,j+m1) = w(i,j)
              ai(i,j+m1) = w(i,j+m2)
   20     continue

C   --- (a11 - a12 a22^-1 a21) in a11 ---
          call yysp2a(mode0,1,nm1,1,nm1,apr,api,ldap,ija,offs,ar,ai,lda)
C       call yygemm('N','N',m1,m1,m2,-1d0,ar(1,1+m1),ai(1,1+m1),lda,
C    .    ar(1+m1,1),ai(1+m1,1),lda,1d0,ar,ai,lda)
          call yymsbm(modl+2,nm1+1,nn2,1,nm1,m1,apr,api,ldap,ija,offs,
     .    ar,ai,lda,ar,ai,lda)

        else

C   --- a22^-1 a21 in w (uses w(1:m2,1:2m1)); copy back to a21
C        call yygemm('N','N',m2,m1,m2,1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),
C     .    lda,ar(1+m1,1),ai(1+m1,1),lda,0d0,w,w(1,1+m1),ldw)
          call yymsbm(modl,nm1+1,nn2,1,nm1,m2,apr,api,ldap,ija,offs,
     .    ar(1+m1,1),ai(1+m1,1),lda,w,w(1,1+m1),ldw)
          do  22  j = 1, m1
            do  22  i = 1, m2
              ar(i+m1,j) = w(i,j)
              ai(i+m1,j) = w(i,j+m1)
   22     continue

C   --- (a11 - a12 a22^-1 a21) in a11 ---
          call yysp2a(mode0,1,nm1,1,nm1,apr,api,ldap,ija,offs,ar,ai,lda)
          call yysbmm(modl+2,1,nm1,nm1+1,nn2,m1,apr,api,ldap,ija,offs,
     .    ar,ai,lda,ar,ai,lda)
        endif
      endif

      if (.not. ltrns) then

C   --- b1 <- b1 - a12 a22^-1 b2 ---
        call yygemm(cm,cm,m1,nb,m2,-1d0,ar(1,1+m1),ai(1,1+m1),lda,
     .  br(1+m1,1),bi(1+m1,1),ldb,1d0,br,bi,ldb)
C       call yprm(.false.,'b1 - a12 a22^-1 b2',2,6,fmt,br,ldb,m1,ldb,nb)

C   --- x1 = c11 (b1 - a12 a22^-1 b2) ---
        call yyqnvb(cs,ar,ai,lda,m1,nb,w,ldw,w2,br,bi,ldb,ierr)
        if (ierr .ne. 0 .or. nb .eq. 0 .or. lsx1) return

C   --- w2 <- b2 - a21 x1 ---
C       call yygemm(cm,cm,m2,nb,m1,1d0,ar(1+m1,1),ai(1+m1,1),lda,
C    .    br,bi,ldb,0d0,w2,w2(1+m2*nb),m2)
        call yysbmm(modl,nm1+1,nn2,1,nm1,nb,apr,api,ldap,ija,offs,
     .  br,bi,lda,w2(1-m1),w2(1+nb*m2-m1),m2)
        jr = -m2
        ji = m2*nb - m2
        do  40  j = 1, nb
          jr = jr+m2
          ji = ji+m2
          do  42  i = 1, m2
            w2(i+jr) = br(i+m1,j) - w2(i+jr)
   42     w2(i+ji) = bi(i+m1,j) - w2(i+ji)
   40   continue
        if (lsx2) return

C   --- x2 = a22^-1 (b2 - a21 x1) ---
        call yygemm(cm,cm,m2,nb,m2,1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,
     .  w2,w2(1+m2*nb),m2,0d0,br(1+m1,1),bi(1+m1,1),ldb)
      else

C   --- b1 <- b1 - b2 a22^-1 a21 ---
        call yygemm(cm,cm,nb,m1,m2,-1d0,br(1,1+m1),bi(1,1+m1),ldb,
     .  ar(1+m1,1),ai(1+m1,1),lda,1d0,br,bi,ldb)

C   --- x1 = c11 (b1 - b2 a22^-1 a21) ---
        call yyqnvb(cs,ar,ai,lda,m1,nb,w,ldw,w2,br,bi,ldb,ierr)
        if (ierr .ne. 0 .or. nb .eq. 0 .or. lsx1) return

C   --- w2 <- b2 - x1 a12 ---
C       call yygemm(cm,cm,nb,m2,m1,1d0,br,bi,ldb,ar(1,1+m1),ai(1,1+m1),
C    .    lda,0d0,w2,w2(1+nb*m2),nb)
        call yymsbm(modl,1,nm1,nm1+1,nn2,nb,apr,api,ldap,ija,offs,
     .  br,bi,lda,w2(1-nb*m1),w2(1+nb*(m2-m1)),nb)
        jr = -nb
        ji = nb*m2 - nb
        do  140  j = 1, m2
          jr = jr+nb
          ji = ji+nb
          do  142  i = 1, nb
            w2(i+jr) = br(i,j+m1) - w2(i+jr)
  142     w2(i+ji) = bi(i,j+m1) - w2(i+ji)
  140   continue
        if (lsx2) return

C   --- x2 = (b2 - a21 x1) a22^-1 ---
        call yygemm(cm,cm,nb,m2,m2,1d0,w2,w2(1+nb*m2),nb,
     .  ar(1+m1,1+m1),ai(1+m1,1+m1),lda,0d0,br(1,1+m1),bi(1,1+m1),ldb)

      endif

      end


