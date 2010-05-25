      subroutine yyqnvb(cs,ar,ai,lda,n,nb,w,ldw,w2,br,bi,ldb,ierr)
C- Solution of a x = b by vectorizable multiplications and inversions
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          'h'  a is assumed hermitian.
Ci          'b'  Assume partial inverse for a is already performed.
Ci               ar,ai must be preserved between successive calls.
Ci          '1'  ignored by yyqnvb, for compatibility with ysbnvb
Ci          '2'  ignored by yyqnvb, for compatibility with ysbnvb
Ci          '4'  Do multiplications using standard four real operations
Ci               (slower, but avoids additions and subtractions that
Ci                can reduce machine precision)
Ci   ar,ai :real,imaginary parts of lhs of eqation a x = b
Ci   lda   :leading dimension of ar,ai
Ci   n     :solve a x = b for matrix a(1:n,1:n)
Ci   w,ldw :a double precision work array of dimension ldw*(n+1)
Ci   w2    :a double precision work array of dimension nb*(n+1)*2
Ci          NB: w and w2 may use the same address space
Ci   br,bi :real and imaginary parts of rhs of equation a x = b
Ci   ldb   :leading dimension of br,bi
Ci   nb    :the number of columns (rows, if cs contains 't') to solve.
Co Outputs:
Co   a     :is OVERWRITTEN, into a partially decomposed form
Co   ierr  :is returned nonzero if matrix was not successfully inverted
Co   br,bi :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
Cr Remarks
Cr   yyqnvb uses a block decomposition to solve the linear system
Cr   a x = b.  a is partitioned into subblocks a11,a21,a12,a22 and
Cr   b is similarly partioned into b1,b2. yyqnvb partially inverts by:
Cr     (1) inverting a22 using yyqinv
Cr     (2) generating (a^-1)_11 = (a11 - a12 a22^-1 a21)
Cr     (3) (a^-1)_11 = (a11 - a12 a22^-1 a21)^-1 generated using yyqinv
Cr   From this decomposition, x is solved as:
Cr     (3) x1 = (a^-1)_11 (b1 - a12 a22^-1 b2)
Cr     (4) x2 = a22^-1 (b2 - a21 x1)
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
Cr
Cr   See yqinvb for an equivalent version where real and imaginary 
Cr   parts are linked by an an offset.
Cb Bugs
Cb   yyqnvb fails if a22 is singular, even if a is not.
Cb   Similarly yyqinv may fail to invert a22 even if it is not singular.
C ----------------------------------------------------------------
C     implicit none
      integer n,lda,ldw,ldb,ierr,nb
      character cs*(*)
      double precision ar(lda,1),ai(lda,1),br(ldb,1),bi(ldb,1),
     .w(ldw,*),w2(nb*(n+1))
C Local variables
      logical lsx,ltrns
      integer nmin,m1,m2,i,j,jr,ji
      parameter(nmin=32)
      double precision ww(2,2)
      character css*1,cm*2
#if DEBUG
      character*10 fmt
      data fmt /'(8f16.10)'/
      print *, 'entering yyqnvb, n=', n
#endif

      ierr = 0
      ltrns = .false.
      lsx = .false.
      css = ' '
      cm = 'N '
      j = len(cs)
      do  2  i = 1, j
        if (cs(i:i) .eq. 't') then
          ltrns = .true.
        elseif (cs(i:i) .eq. 'b') then
          lsx = .true.
        elseif (cs(i:i) .eq. 'h') then
          css = 'h'
        elseif (cs(i:i) .eq. '4') then
          cm = 'N4'
        elseif (cs(i:i) .eq. '1' .or. cs(i:i) .eq. '2') then
        elseif (cs(i:i) .ne. ' ') then
          call rxs2('yyqnvb: unrecognized switch, cs=''',cs,'''')
        endif
    2 continue

C --- For n<nmin, do directly by yygefa,yygesl ---
      if (n .lt. nmin) then
        if (.not. lsx) then
          call yygefa(ar,ai,lda,n,w,ierr)
          if (ierr .ne. 0) return
        endif
        do  12  j = 1, nb
          if (cs(1:1) .ne. 't') then
            call yygesl(ar,ai,lda,n,w,br(1,j),bi(1,j),0)
          else
C       ... Do for special case n<=2: w may not be large enough!
            if (n .le. 2) then
              do  116  i = 1, n
                ww(i,1) = br(j,i)
  116         ww(i,2) = bi(j,i)
              call yygesl(ar,ai,lda,n,w,ww(1,1),ww(1,2),1)
              do  118  i = 1, n
                br(j,i) = ww(i,1)
                bi(j,i) = ww(i,2)
  118         continue
            else
              do  16  i = 1, n
                w(i,2) = br(j,i)
   16         w(i,3) = bi(j,i)
              call yygesl(ar,ai,lda,n,w,w(1,2),w(1,3),1)
              do  18  i = 1, n
                br(j,i) = w(i,2)
                bi(j,i) = w(i,3)
   18         continue
            endif
          endif
   12   continue
        return
      endif

      if (ldw .lt. n) call rx('yyqnvb: ldw lt n')
      m1 = (n+1)/2
      m2 = n-m1

C --- Decompose matrix ---
      if (.not. lsx) then

C --- a22^-1 in a22 ---
        call yyqinv(css,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,2,m2,
     .  w,ldw,ierr)
C     call yprm(.false.,'a22^-1',2,6,fmt,ar(1+m1,1+m1),lda,m2,lda,m2)
        if (ierr .ne. 0) return

        if (.not. ltrns) then

C   --- a12 a22^-1 in w; copy back to a12 ---
          call yygemm('N','N',m1,m2,m2,1d0,ar(1,1+m1),ai(1,1+m1),lda,
     .    ar(1+m1,1+m1),ai(1+m1,1+m1),lda,0d0,w,w(1,1+m1),ldw)
          do  20  j = 1, m2
            do  20  i = 1, m1
              ar(i,j+m1) = w(i,j)
   20     ai(i,j+m1) = w(i,j+m1)

        else

C   --- a22^-1 a21 in w; copy back to a21 ---
          call yygemm('N','N',m2,m1,m2,1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),
     .    lda,ar(1+m1,1),ai(1+m1,1),lda,0d0,w,w(1+m2,1),ldw)
          do  22  j = 1, m1
            do  22  i = 1, m2
              ar(i+m1,j) = w(i,j)
   22     ai(i+m1,j) = w(i+m2,j)

        endif

C --- c11 = (a11 - a12 a22^-1 a21)^-1 in a11 ---
        call yygemm('N','N',m1,m1,m2,-1d0,ar(1,1+m1),ai(1,1+m1),lda,
     .  ar(1+m1,1),ai(1+m1,1),lda,1d0,ar,ai,lda)
        call yyqinv(css,ar,ai,lda,2,m1,w,ldw,ierr)
        if (ierr .ne. 0) return
C     call yprm(.false.,'c11',2,6,fmt,ar,lda,m1,lda,m1)
      endif
C ... End of matrix decomposition

      if (nb .eq. 0) return
      if (.not. ltrns) then

C --- w2 <- b1 - a12 a22^-1 b2 ---
        call yygemm(cm,cm,m1,nb,m2,1d0,ar(1,1+m1),ai(1,1+m1),lda,
     .  br(1+m1,1),bi(1+m1,1),ldb,0d0,w2,w2(1+m1*nb),m1)
        jr = -m1
        ji = m1*nb - m1
        do  30  j = 1, nb
          jr = jr+m1
          ji = ji+m1
          do  32  i = 1, m1
            w2(i+jr) = br(i,j) - w2(i+jr)
   32     w2(i+ji) = bi(i,j) - w2(i+ji)
   30   continue
C     call yprm(.false.,'b1 - a12 a22^-1 b2',2,6,fmt,w2,m1,m1,nb,nb)

C --- x1 = c11 (b1 - a12 a22^-1 b2) ---
        call yygemm(cm,cm,m1,nb,m1,1d0,ar,ai,lda,w2,w2(1+m1*nb),m1,0d0,
     .  br,bi,ldb)

C --- w2 <- b2 - a21 x1 ---
        call yygemm(cm,cm,m2,nb,m1,1d0,ar(1+m1,1),ai(1+m1,1),lda,
     .  br,bi,ldb,0d0,w2,w2(1+m2*nb),m2)
        jr = -m2
        ji = m2*nb - m2
        do  40  j = 1, nb
          jr = jr+m2
          ji = ji+m2
          do  42  i = 1, m2
            w2(i+jr) = br(i+m1,j) - w2(i+jr)
   42     w2(i+ji) = bi(i+m1,j) - w2(i+ji)
   40   continue

C --- x2 = a22^-1 (b2 - a21 x1) ---
        call yygemm(cm,cm,m2,nb,m2,1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,
     .  w2,w2(1+m2*nb),m2,0d0,br(1+m1,1),bi(1+m1,1),ldb)
C     call yprm(.false.,'a^-1 b',2,6,fmt,br,ldb,n,ldb,nb)

      else

C --- w2 <- b1 - b2 a22^-1 a21 ---
        call yygemm(cm,cm,nb,m1,m2,1d0,br(1,1+m1),bi(1,1+m1),ldb,
     .  ar(1+m1,1),ai(1+m1,1),lda,0d0,w2,w2(1+nb*m1),nb)
        jr = -nb
        ji = nb*m1 - nb
        do  130  j = 1, m1
          jr = jr+nb
          ji = ji+nb
          do  132  i = 1, nb
            w2(i+jr) = br(i,j) - w2(i+jr)
  132     w2(i+ji) = bi(i,j) - w2(i+ji)
  130   continue

C --- x1 = R1 c11 ---
        call yygemm(cm,cm,nb,m1,m1,1d0,w2,w2(1+nb*m1),nb,
     .  ar,ai,lda,0d0,br,bi,ldb)
C     call yprm(.false.,'x1',2,6,fmt,br,ldb,nb,ldb,m1)

C --- w2 <- b2 - x1 a12 ---
        call yygemm(cm,cm,nb,m2,m1,1d0,br,bi,ldb,ar(1,1+m1),ai(1,1+m1),
     .  lda,0d0,w2,w2(1+nb*m2),nb)
        jr = -nb
        ji = nb*m2 - nb
        do  140  j = 1, m2
          jr = jr+nb
          ji = ji+nb
          do  142  i = 1, nb
            w2(i+jr) = br(i,j+m1) - w2(i+jr)
  142     w2(i+ji) = bi(i,j+m1) - w2(i+ji)
  140   continue

C --- x2 = (b2 - a21 x1) a22^-1 ---
        call yygemm(cm,cm,nb,m2,m2,1d0,w2,w2(1+nb*m2),nb,
     .  ar(1+m1,1+m1),ai(1+m1,1+m1),lda,0d0,br(1,1+m1),bi(1,1+m1),ldb)
C     call yprm(.false.,'b a^-1',2,6,fmt,br,ldb,nb,ldb,n)

      endif

      end

