#define F90 1
#ifndef F90
      subroutine dqinvb(cs,a,lda,nlev,n,nb,w,ldw,w2,b,ldb,ierr)
#else
      recursive subroutine dqinvb(cs,a,lda,nlev,n,nb,w,ldw,w2,b,ldb,
     .  ierr)
#endif
C- Solution of a x = b by vectorizable multiplications and inversions
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          's'  a is assumed symmetric.
Ci          'b'  Assume partial inverse for a is already performed.
Ci               ar,ai must be preserved between successive calls.
Ci   a,lda :First matrix in eqation a x = b, and leading dimension
Ci   nlev  :the maximum number of recursion levels allowed.
Ci          To avoid roundoff errors, nlev=2 is suggested.
Ci   n     :solve a x = b for matrix a(1:n,1:n)
Ci   w,ldw :a double precision work array of dimension ldw*n:
Ci   w2    :a double precision work array of dimension nb*(n+1)
Ci          w and w2 may use the same address space
Ci   b,ldb :the rhs, and its leading dimension
Ci   nb    :the number of columns (rows, if cs contains 't') in x and b
Co Outputs:
Co   a     :is OVERWRITTEN, into a partially decomposed form
Co   ierr  :is returned nonzero if matrix was not successfully inverted
Co   a^-1 b (b a^-1) is returned in b(1..n,1..nb) (b(1..nb,1..n))
Cr Remarks
Cr   dqinvb uses a block decomposition to solve the linear system
Cr   a x = b.  a is partitioned into subblocks a11,a21,a12,a22 and
Cr   b is similarly partioned into b1,b2.  Let c be the inverse,
Cr   with subblocks c11,c21,c12,c22.  Solution proceeds by:
Cr     (1) invert   a22 -> a22^-1
Cr     (2) multiply a12 a22^-1
Cr     (3) multiply to obtain c11^-1 = (a11 - a12 a22^-1 a21)
Cr     ... At this point there are two options.  Standard option:
Cr     (4) invert   c11^-1 -> c11 = (a11 - a12 a22^-1 a21)^-1
Cr     (5) multiply to obtain x1 = c11 (b1 - (a12 a22^-1) b2)
Cr     (6) multiply to obtain x2 = a22^-1 (b2 - a21 x1)
Cr     ... the recursive option omits explicit inversion for c11:
Cr     (4) multiply to obtain (b1 - a12 a22^-1 b2)
Cr     (5) solve    c11^-1 x1 = (b1 - a12 a22^-1 b2) calling dqinvb
Cr     (6) multiply to obtain x2 = a22^-1 (b2 - a21 x1)
Cr
Cr   If the decomposition has already been effected, x may be calculated
Cr   for new vectors b without decomposing the matrix again (cs='b')
Cr   Debugging:
Cr   mc a -split a 1,4,6 1,4,6 a22 -i a12 -tog -x a21 -x -s-1 a11 -+ -i
Cr   -a P11 b -split b 1,4,6 1,6 a12 a22 -i -x b21 -x -s-1 b11 -+
Cr   P11 -tog -x -a x1 a21 x1 -x -s-1 b21 -+ a22 -i -tog -x -a x2
Cr   x1 x2 -rcat a -i b -x --
Cr
Cr   Alternative, swapping indices 1 and 2
Cr   mc a -split a 1,4,nr+1 1,4,nr+1 a21 a11 -i -x a12 -x -s-1 a22 -+ -i
Cr   -a P22 b -split b 1,4,nr+1 1,nc+1 a21 a11 -i -x b11 -x -s-1 b21 -+
Cr   P22 -tog -x -a x2 a12 x2 -x -s-1 b11 -+ a11 -i -tog -x x2 -rcat
Cr   a -i b -x -- -px
Cr
Cr   For the transpose case:
Cr   mc a -split a 1,4,6 1,4,6 a22 -i a12 -tog -x a21 -x -s-1 a11 -+ -i
Cr   -a P11 b -split b 1,5 1,4,6 a22 -i a21 -x b12 -tog -x -s-1 b11 -+
Cr   P11 -x -a x1 x1 a12 -x -s-1 b12 -+ a22 -i -x -a x2 x1 x2 -ccat
Cr   b a -i -x --
Cb Bugs
Cb   dqinvb fails if a22 is singular, even if a is not.
Cb   Similarly, dqinv may fail to invert a22 even if it is not singular.
C ----------------------------------------------------------------
C     implicit none
      integer n,nlev,lda,ldw,ldb,ierr,nb
      double precision a(lda,1),b(ldb,1),w(ldw,2),w2(nb*(n+1))
C Local variables
      logical lsx,ltrns,lrecur
      double precision ww(2)
      character cs*(*),css*1
      integer nmin,nmin0,m1,m2,i,j,jr,isw
C     parameter(nmin=70)
      save nmin
      data nmin /32/
#if DEBUG
c     character*10 fmt
C     data fmt /'(8f16.10)'/
C     call ywrm(0,'a',1,6,fmt,a,1,lda,n,n)
      print 334, 'entering dqinvb',cs,nlev,n,(n+1)/2,n-(n+1)/2
  334 format(1x,a,': cs=',a,': nlev=',i1,
     .       '  n=',i3,'  m1=',i3,'  m2=',i3)
#endif

      ierr = 0
      ltrns = .false.
      lrecur = .false.
      lsx = .false.
      css = ' '
      j = len(cs)
      do  2  i = 1, j
        if (cs(i:i) .eq. 't') then
          ltrns = .true.
        elseif (cs(i:i) .eq. 'b') then
          lsx = .true.
        elseif (cs(i:i) .eq. 's') then
          css = 's'
        elseif (cs(i:i) .eq. 'r') then
          lrecur = .true.
        elseif (cs(i:i) .ne. ' ') then
          call rxs('dqinvb: bad input cs',cs)
        endif
    2 continue
#if ! (RECUR | F90)
      lrecur = .false.
#endif

C --- For n<nmin, do directly by dgefa,dgesl ---
      if (n .lt. nmin) then
        if (.not. lsx) then
          if (css .eq. 's') call dsifa(a,lda,n,w,ierr)
          if (css .ne. 's') call dgefa(a,lda,n,w,ierr)
          if (ierr .ne. 0) return
        endif
        do  12  j = 1, nb
          if (.not. ltrns .and. css .eq. 's') then
            call dsisl(a,lda,n,w,b(1,j))
          elseif (.not. ltrns .and. css .ne. 's') then
            call dgesl(a,lda,n,w,b(1,j),0)
          else
C       ... Do for special case n<=2: w may not be large enough!
            if (n .le. 2) then
              do  116  i = 1, n
              ww(i) = b(j,i)
              if (css .eq. 's') call dsisl(a,lda,n,w,ww)
              if (css .ne. 's') call dgesl(a,lda,n,w,ww,1)
  116         continue
              do  118  i = 1, n
  118         b(j,i) = ww(i)
            else
              do  16  i = 1, n
   16         w(i,2) = b(j,i)
              if (css .eq. 's') call dsisl(a,lda,n,w,w(1,2))
              if (css .ne. 's') call dgesl(a,lda,n,w,w(1,2),1)
              do  18  i = 1, n
   18         b(j,i) = w(i,2)
            endif
          endif
   12   continue
        return
      endif

      if (ldw .lt. n) call rx('dqinvb: ldw lt n')
      m1 = (n+1)/2
      m2 = n-m1

C --- Decompose matrix ---
      if (.not. lsx) then

C --- a22^-1 in a22 ---
      call dqinv(css,a(1+m1,1+m1),lda,max(nlev-1,0),m2,w,ldw,ierr)
C     call ywrm(0,'a22^-1',1,6,fmt,a(1+m1,1+m1),1,lda,m2,m2)
      if (ierr .ne. 0) return

      if (.not. ltrns) then

C   --- a12 a22^-1 in w; copy back to a12 ---
        call dgemm('N','N',m1,m2,m2,1d0,a(1,1+m1),lda,a(1+m1,1+m1),lda,
     .    0d0,w,ldw)
        do  20  j = 1, m2
        do  20  i = 1, m1
   20   a(i,j+m1) = w(i,j)

      else

C   --- a22^-1 a21 in w; copy back to a21 ---
        call dgemm('N','N',m2,m1,m2,1d0,a(1+m1,1+m1),lda,
     .      a(1+m1,1),lda,0d0,w,ldw)
        do  22  j = 1, m1
        do  22  i = 1, m2
   22   a(i+m1,j) = w(i,j)

      endif

C --- c11^-1 = (a11 - a12 a22^-1 a21) in a11 ---
      if (cs .eq. 's') then
        call dsmpy(m1,m2,a(1,1+m1),lda,a(1+m1,1),lda,0d0,w,ldw)
        do  24  j = 1, m1
        do  24  i = 1, m1
   24   a(i,j) = a(i,j) - w(i,j)
      else
        call dgemm('N','N',m1,m1,m2,-1d0,a(1,1+m1),lda,a(1+m1,1),lda,
     .    1d0,a,lda)
      endif

C --- c11 = (a11 - a12 a22^-1 a21)^-1 in a11 ---
      if (.not. (lrecur .and. nlev .gt. 0 .and. m1 .ge. nmin))
     .  call dqinv(css,a,lda,max(nlev-1,0),m1,w,ldw,ierr)
      if (ierr .ne. 0) return
C     call yprm(.false.,'c11',1,6,fmt,a,lda,m1,lda,m1)
      endif
C ... End of matrix decomposition

      if (nb .eq. 0) return
      if (.not. ltrns) then

C --- Obtain x1 = c11 (b1 - a12 a22^-1 b2) ---
      if (lrecur .and. nlev .gt. 0 .and. m1 .ge. nmin) then
#if RECUR | F90
C   ... (b1 - a12 a22^-1 b2)
        call dgemm('N','N',m1,nb,m2,-1d0,a(1,1+m1),lda,b(1+m1,1),ldb,
     .    1d0,b,ldb)
C   ... Solve c11^-1 x1 - (b1 - a12 a22^-1 b2)
        call dqinvb(cs,a,lda,nlev-1,(m1),nb,w,ldw,w2,b,ldb,ierr)
C       AIX and gcc require this
        m1 = (n+1)/2
        m2 = n-m1
#endif
#if DEBUG
C       print 333, 'recursive call exit',nlev,n,m1,m2
#endif
      else
C   ... w2 <- b1 - a12 a22^-1 b2
        call dgemm('N','N',m1,nb,m2,1d0,a(1,1+m1),lda,b(1+m1,1),ldb,
     .    0d0,w2,m1)
        jr = -m1
        do  30  j = 1, nb
          jr = jr+m1
          do  32  i = 1, m1
   32     w2(i+jr) = b(i,j) - w2(i+jr)
   30   continue
C       call yprm(.false.,'b1 - a12 a22^-1 b2',1,6,fmt,w2,m1,m1,nb,nb)

C   ... x1 = c11 (b1 - a12 a22^-1 b2)
        call dgemm('N','N',m1,nb,m1,1d0,a,lda,w2,m1,0d0,b,ldb)
      endif

C --- w2 <- b2 - a21 x1 ---
      call dgemm('N','N',m2,nb,m1,1d0,a(1+m1,1),lda,b,ldb,0d0,w2,m2)
      jr = -m2
      do  40  j = 1, nb
        jr = jr+m2
        do  42  i = 1, m2
   42   w2(i+jr) = b(i+m1,j) - w2(i+jr)
   40 continue

C --- x2 = a22^-1 (b2 - a21 x1) ---
      call dgemm('N','N',m2,nb,m2,1d0,a(1+m1,1+m1),lda,w2,m2,
     .  0d0,b(1+m1,1),ldb)
C     call yprm(.false.,'a^-1 b',1,6,fmt,br,ldb,n,ldb,nb)

      else

C --- Obtain x1 = (b1 - b2 a22^-1 a21) c11 ---
      if (lrecur .and. nlev .gt. 0 .and. m1 .ge. nmin) then
#if RECUR | F90
C   ... (b1 - a12 a22^-1 b2)
        call dgemm('N','N',nb,m1,m2,-1d0,b(1,1+m1),ldb,
     .    a(1+m1,1),lda,1d0,b,ldb)
C   ... solve x1 c11^-1 = (b1 - b2 a22^-1 a21)
        call dqinvb(cs,a,lda,nlev-1,(m1),nb,w,ldw,w2,b,ldb,ierr)
C       AIX and gcc require the next two lines
        m1 = (n+1)/2
        m2 = n-m1
#endif
      else
C ... w2 <- b1 - b2 a22^-1 a21
      call dgemm('N','N',nb,m1,m2,1d0,b(1,1+m1),ldb,
     .  a(1+m1,1),lda,0d0,w2,nb)
      jr = -nb
      do  130  j = 1, m1
        jr = jr+nb
        do  132  i = 1, nb
  132   w2(i+jr) = b(i,j) - w2(i+jr)
  130 continue
C ... x1 = (b1 - b2 a22^-1 a21) c11
      call dgemm('N','N',nb,m1,m1,1d0,w2,nb,a,lda,0d0,b,ldb)
C     call yprm(.false.,'x1',1,6,fmt,br,ldb,nb,ldb,m1)
      endif

C --- w2 <- b2 - x1 a12 ---
      call dgemm('N','N',nb,m2,m1,1d0,b,ldb,a(1,1+m1),
     .  lda,0d0,w2,nb)
      jr = -nb
      do  140  j = 1, m2
        jr = jr+nb
        do  142  i = 1, nb
  142   w2(i+jr) = b(i,j+m1) - w2(i+jr)
  140 continue

C --- x2 = (b2 - a21 x1) a22^-1 ---
      call dgemm('N','N',nb,m2,m2,1d0,w2,nb,
     .  a(1+m1,1+m1),lda,0d0,b(1,1+m1),ldb)
C     call yprm(.false.,'b a^-1',1,6,fmt,br,ldb,nb,ldb,n)

      endif
#if DEBUG
      print 333, 'exiting dqinvb ',nlev,n,m1,m2
  333 format(1x,a,': nlev=',i1,'  ndim=',i4,' partitioned into',2i4)
#endif

      return


      entry dqnvb0(nmin0,isw)
      if (isw .gt. 0) then
        nmin = nmin0
      else
        nmin0 = nmin
      endif


      end


