CSFCPP#define F90 1
CSFCPP#ifndef F90
CSFCPP      subroutine yyqinv(cs,ar,ai,lda,nlev,n,w,ldw,ierr)
CSFCPP#else
      recursive subroutine yyqinv(cs,ar,ai,lda,nlev,n,w,ldw,ierr)
CSFCPP#endif
C- Inversion of a complex matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   :if 'h', a is assumed hermitian
Ci   a     :matrix to be inverted
Ci   lda   :leading dimension of a
Ci   n     :rank of the matrix to be inverted
Ci   nlev  :the maximum number of recursion levels allowed.
Ci          To avoid roundoff errors, nlev=2 is suggested.
Ci   w     :double precision work array of dimension ldw*(n+1)
Ci   ldw   :leading dimension of w
Co Outputs:
Co   a     :is overwritten by inverse of input a
Co   ierr  :returned nonzero if matrix was not fully inverted.
Cb Bugs
Cb   The algorithm fails if a22 is singular, even if a is not.
Cb   Similarly, if smaller subblocks are singular, yyqinv may fail
Cb   when called recursively.
Cr Remarks:
Cr   See Numerical Recipes, 2.11.
Cr   It is more efficient to allow inversion to proceed recursively,
Cr   if your compiler allows it.  Recursion proceeds provided nlev>0
Cr   and the dimension of the matrix to be inverted exceeds nmin.
Cr   The latter should be chosen where the tradoff between
Cr   the extra overhead and fewer O(N^3) operations takes place.
Cr   The AIX xlf90 compiler overwrites locally declared variables,
Cr   (NB: compile using -qrecur).  Compiling with AIX code uncommented
Cr   restores necessary variables after a recursive call.
Cr   (For debugging:)
Cr   mc -f4f16.10 a -split a 1,3,5 1,3,5  a11 -i -a r1  a21 r1 -x -a r2\
Cr   r2 a12 -x -a r4  r4 a22 -- -a r5  r5 -i -a r6  r1 a12 -x -a r3\
Cr   r3 r6 -x -a c12  r6 r2 -x -a c21  r3 c21 -x -a r7 r1 r7 -- -a c11\
Cr   r6 -s-1 -a c22 c11 c12 -ccat c21 c22 -ccat -rcat
Cr   Alternative, starting with a22
Cr   mc -f4f16.10 a -split a 1,3,5 1,3,5  a22 -i -a r1  a12 r1 -x -a r2\
Cr   r2 a21 -x -a r4  r4 a11 -- -a r5  r5 -i -a r6  r1 a21 -x -a r3\
Cr   r3 r6 -x -a c21  r6 r2 -x -a c12  r3 c12 -x -a r7 r1 r7 -- -a c22
Cr   r6 -s-1 -a c11 c11 c12 -ccat c21 c22 -ccat -rcat
Cr
Cr   See yqinv for an equivalent version where real and imaginary
Cr   parts are linked by an an offset.
C ----------------------------------------------------------------
C     implicit none
      integer n,lda,ldw,ierr,nlev
      double precision ar(lda,1),ai(lda,1),w(ldw,3)
      double precision xx
      character*1 cs
      integer nmin,m1,m2,i,j,isw,nmin0
      save nmin
#if DEBUG
      character*10 fmt
      data fmt /'(8f16.10)'/
#endif
      data nmin /32/

#if DEBUG
      print 334, 'entering yyqinv: cs=',cs,nlev,n
  334 format(1x,a,a,': nlev=',i1,'  n=',i4)
#endif
      m1 = n/2
      m2 = n-m1

C     call tcn('yyqinv')

C --- Straight yygefa,di if n lt nmin ---
      if (n .le. nmin) then
        if (cs .ne. 'h' .or. n .eq. 1) then
          call yygefa(ar,ai,lda,n,w,ierr)
          if (ierr .ne. 0) goto 99
          if (n .eq. 1) then
            call yygedi(ar,ai,lda,n,w,xx,w(1,2),xx,1)
          else
            call yygedi(ar,ai,lda,n,w,xx,w(1,2),w(1,3),1)
          endif
        else
          call yyhifa(ar,ai,lda,n,w,ierr)
          if (ierr .ne. 0) goto 99
          call yyhidi(ar,ai,lda,n,w,xx,i,w(1,2),w(1,3),1)
          do  10  i = 1, n
            do  10  j = 1, i
              ar(i,j) =  ar(j,i)
   10     ai(i,j) = -ai(j,i)
        endif
#if DEBUG
        print *, 'exit yyqinv, LU decomposition'
#endif
        goto 99
      endif

      if (lda .lt. n) call rx('yyqinv: lda lt n')
      if (ldw .lt. n) call rx('yyqinv: ldw lt n')

C --- R1 = a11^-1 in a11 ---
CSFCPP#if RECUR | F90
      if (m1 .lt. nmin .or. nlev .eq. 0) then
CSFCPP#endif
        if (cs .eq. 'h') then
          call yyhifa(ar,ai,lda,m1,w,ierr)
          if (ierr .ne. 0) goto 99
          call yyhidi(ar,ai,lda,m1,w,xx,i,w(1,2),w(1,3),1)
          do  12  i = 1, m1
            do  12  j = 1, i
              ar(i,j) =  ar(j,i)
   12     ai(i,j) = -ai(j,i)
        else
          call yygefa(ar,ai,lda,m1,w,ierr)
          if (ierr .ne. 0) goto 99
          call yygedi(ar,ai,lda,m1,w,xx,w(1,2),w(1,3),1)
        endif
CSFCPP#if RECUR | F90
      else
#if DEBUG
        print *, 'recursive call to yyqinv n,m1,m2=',n,m1,m2
#endif
#if AIX | SUN_ULTRA
        call yyqinv(cs,ar,ai,lda,nlev-1,(m1),w,ldw,ierr)
        m1 = n/2
        m2 = n-m1
#else
        call yyqinv(cs,ar,ai,lda,nlev-1,m1,w,ldw,ierr)
#endif
#if DEBUG
        print 333, 'recursive call exit',nlev,n,m1,m2
#endif
        if (ierr .ne. 0) goto 99
      endif
CSFCPP#endif
*     call yprm(.false.,'R1',2,6,fmt,ar,lda,m1,lda,m1)

C --- R3 = R1 a12 in w: uses w(1:m1,1:2*m2) ---
      call yygemm('N','N',m1,m2,m1,1d0,ar,ai,lda,ar(1,1+m1),ai(1,1+m1),
     .lda,0d0,w,w(1,1+m2),ldw)
*     call yprm(.false.,'R3',2,6,fmt,w,lda,m1,m2,m2)

C --- R4 = a21 R3 = a21 a11^-1 a12 in w2; uses w(1+m1:n,1:2*m2) ---
      call yygemm('N','N',m2,m2,m1,1d0,ar(1+m1,1),ai(1+m1,1),lda,
     .w,w(1,1+m2),ldw,0d0,w(1+m1,1),w(1+m1,1+m2),ldw)
*     call yprm(.false.,'R4',2,6,fmt,w(1+m1,1),ldw,m2,m2,m2)

C --- -R5 = a22 - R4 = a22 - a21 a11^-1 a12 in a22 ---
      do  20  j = 1, m2
        do  20  i = 1, m2
          ar(i+m1,j+m1) = ar(i+m1,j+m1) - w(i+m1,j)
   20 ai(i+m1,j+m1) = ai(i+m1,j+m1) - w(i+m1,j+m2)
C     R4 = w2 isn't needed beyond this point
*     call yprm(.false.,'-R5',2,6,fmt,ar(1+m1,1+m1),lda,m2,lda,m2)

C --- c22 = -R6 = -R5^-1 = (a22 - a21 a11^-1 a12)^-1 in a22 ---
CSFCPP#if RECUR | F90
      if (m2 .lt. nmin .or. nlev .eq. 0) then
CSFCPP#endif
        if (cs .eq. 'h') then
          call yyhifa(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w(1+m1,1),ierr)
          if (ierr .ne. 0) goto 99
          call yyhidi(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w(1+m1,1),
     .    xx,i,w(1+m1,2),w(1+m1,3),1)
          do  16  i = m1+1, n
            do  16  j = m1+1, i
              ar(i,j) =  ar(j,i)
   16     ai(i,j) = -ai(j,i)
        else
          call yygefa(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w(1+m1,1),ierr)
          if (ierr .ne. 0) goto 99
          call yygedi(ar(1+m1,1+m1),ai(1+m1,1+m1),lda,m2,w(1+m1,1),xx,
     .    w(1+m1,2),w(1+m1,3),1)
        endif
CSFCPP#if RECUR | F90
      else
#if DEBUG
        print *, 'recursive call to yyqinv n,m1,m2=',n,m1,m2
#endif
#if AIX | SUN_ULTRA
        call yyqinv(cs,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,
     .  nlev-1,(m2),w(1+m1,1),ldw,ierr)
        m1 = n/2
        m2 = n-m1
#else
        call yyqinv(cs,ar(1+m1,1+m1),ai(1+m1,1+m1),lda,
     .  nlev-1,m2,w(1+m1,1),ldw,ierr)
#endif
#if DEBUG
        print 333, 'recursive call exit',nlev,n,m1,m2
#endif
        if (ierr .ne. 0) goto 99
      endif
CSFCPP#endif
*     call yprm(.false.,'-R6',2,6,fmt,ar(1+m1,1+m1),lda,m2,lda,m2)

C --- c12 = R3 * R6 = -a11^-1 a12 c22 in a12 ---
      call yygemm('N','N',m1,m2,m2,-1d0,w,w(1,1+m2),ldw,ar(1+m1,1+m1),
     .ai(1+m1,1+m1),lda,0d0,ar(1,1+m1),ai(1,1+m1),lda)
*     call yprm(.false.,'c12',2,6,fmt,ar(1,1+m1),lda,m1,lda,m2)

C --- Hermitian case: R2 = R3+ and c21=c12+ ---
      if (cs .eq. 'h') then
        do  30  j = 1, m1
          do  30  i = 1, m2
            ar(i+m1,j) =  ar(j,i+m1)
   30   ai(i+m1,j) = -ai(j,i+m1)

C NB: do following in two steps to avoid beta=1 in yygemm call
C   ... R7 = R3 * R6 * R2 = R3 * c21
        call yygemm('N','N',m1,m1,m2,1d0,w,w(1,1+m2),ldw,
     .  ar(1+m1,1),ai(1+m1,1),lda,0d0,w(1+m1,1),w(1+m1,1+m2),ldw)
C   ... c11 = R1 - R7 = R1 - R3 * R6 * R2 = R1 - R3 * c21
        do  32  j = 1, m1
          do  32  i = 1, j
            ar(i,j) =  ar(i,j) - (w(i+m1,j)+w(j+m1,i))*.5d0
   32   ai(i,j) =  ai(i,j) - (w(i+m1,j+m2)-w(j+m1,i+m2))*.5d0

C       call yygemm('N','N',m1,m1,m2,-1d0,w,w(1,1+m2),ldw,
C    .    ar(1+m1,1),ai(1+m1,1),lda,1d0,ar,ai,lda)

C   ... Force hermitian explicitly
        do  34  i = 1, m1
          do  34  j = 1, i
            ar(i,j) =  ar(j,i)
   34   ai(i,j) = -ai(j,i)

C --- Nonhermitian case ---
      else

C   --- R2 = a21 * R1 in w(2); uses w(1+m1:n,1:2*m1) ---
        call yygemm('N','N',m2,m1,m1,1d0,ar(1+m1,1),ai(1+m1,1),lda,
     .  ar,ai,lda,0d0,w(1+m1,1),w(1+m1,1+m1),ldw)
*       call yprm(.false.,'R2',2,6,fmt,w(1+m1,1),ldw,m2,m2,m1)

C   --- c21 = R6 * R2 = -c22 a21 a11^-1 ---
        call yygemm('N','N',m2,m1,m2,-1d0,ar(1+m1,1+m1),ai(1+m1,1+m1),
     .  lda,w(1+m1,1),w(1+m1,1+m1),ldw,0d0,ar(1+m1,1),ai(1+m1,1),lda)
*       call yprm(.false.,'c21',2,6,fmt,ar(1+m1,1),lda,m2,lda,m1)

C   --- c11 = R1 - R3 * R6 * R2 = a11^-1 + a11^-1 a12 c22 a21 a11^-1 ---
        call yygemm('N','N',m1,m1,m2,-1d0,ar(1,1+m1),ai(1,1+m1),lda,
     .  w(1+m1,1),w(1+m1,1+m1),ldw,1d0,ar,ai,lda)
*       call yprm(.false.,'R7',2,6,fmt,w,ldw,m1,m2,m1)

      endif

#if DEBUG
C     call yprm(.false.,'a^-1',2,6,fmt,ar,lda,n,lda,n)
      print 333, 'exiting  yyqinv',nlev,n,m1,m2
  333 format(1x,a,': nlev=',i1,'  ndim=',i3,' partitioned into',2i3)
#endif

   99 continue
C     call tcx('yyqinv')      
      return
#if DEBUG
      call yprm(.false.,'a^-1',2,6,fmt,ar,lda,n,lda,n)
#endif

      entry yqinv0(nmin0,isw)
      if (isw .gt. 0) then
        nmin = nmin0
      else
        nmin0 = nmin
      endif

      end

