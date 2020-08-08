      subroutine zqinvb(cs,a,lda,n,nb, b,ldb,ierr)
C- Solution of a x = b by vectorizable multiplications and inversions. Lapack required.
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          'h'  a is assumed hermitian.
Ci          'b'  Assume partial inverse for a is already performed.
Ci               ar,ai must be preserved between successive calls.
Ci               If the Lapack branch is used, w must also be preserved
Ci          '1'  ignored by zqinvb, for compatibility with ysbnvb
Ci          '2'  ignored by zqinvb, for compatibility with ysbnvb
Ci          '4'  Do multiplications using standard four real operations
Ci               (slower, but avoids additions and subtractions that
Ci                can reduce machine precision)
!!
Ci          'l'  Always call Lapack routines  <---------- we assume this
C
Ci   a     :lhs of equation a x = b
Ci   lda   :leading dimension of ar,ai
Ci   n     :solve a x = b for matrix a(1:n,1:n)
Ci   b     :right hand side of equation a x = b
Ci   ldb   :leading dimension of b
Ci   nb    :the number of columns (rows, if cs contains 't') to solve.
Co Outputs:
Co   a     :is OVERWRITTEN, into a partially decomposed form
Co   ierr  :is returned nonzero if matrix was not successfully inverted
Co   b     :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
c
!! ---takao don't know why ldw=129 ???
C   w,ldw :double precision work array w(ldw,*), and leading dimension.
C         :The required size of w depends on which branch is executed.
C         : w = w(ldw,n); the req'd and optimal
C         :    sizes of ldw depend on the particular branch.
C         :   *If cs contains no 't' or 'h':
C              required ldw>=1; no optimal size
C         :   *If cs contains no 't':
C              required ldw>=1; optimal ldw=129
C         :   *If cs contains 't':
C              required ldw>=3; optimal ldw=129
C ----------------------------------------------------------------
      implicit none
      integer n,lda,ldw,ldb,ierr,nb
      character cs*(*)
      double precision a(lda),b(ldb)
      integer i,j,lwrk,nbi
      logical ltrns,lsx,lherm
      character*1 cc
      real(8),allocatable:: w(:)
      ldw=129
      allocate(w(ldw*n))
      ltrns = .false.
      lsx   = .false.
      lherm = .false.
C     cm = 'N'
      j = len(cs)
      do  2  i = 1, j
        if (cs(i:i) .eq. 't') then
          ltrns = .true.
        elseif (cs(i:i) .eq. 'b') then
          lsx = .true.
        elseif (cs(i:i) .eq. 'h') then
          lherm = .true.
        elseif (cs(i:i) .eq. '4') then
C         cm = 'N4'
        elseif (cs(i:i) .eq. '1' .or. cs(i:i) .eq. '2') then
c        elseif (cs(i:i) .eq. 'l') then
        elseif (cs(i:i) .ne. ' ') then
          call rxs2('zqinv: unrecognized switch, cs=''',cs,'''')
        endif
    2 continue
      ierr = 0
C     LU decomposition
      if (lherm) then
C     Use n elements for pivot; remainder is complex work
         lwrk = (ldw*n - n)/2
         ierr = 0
         if (.not. lsx) call zhetrf('U',n,a,lda,w,w(1+n),lwrk,ierr)
         cc = 'C'
      else
         if (.not. lsx) call zgetrf(n,n,a,lda,w,ierr)
         cc = 'N'
      endif
      if (ierr .ne. 0) return
C     Back substitution
      if (ltrns) then
C     Use n elements for pivot; remainder is complex work
         lwrk = (ldw - 1)/2
         do  i = 1, nb, lwrk
            nbi = min(nb-i+1,lwrk)
C     Copy subblock of b to transpose of w
            call zmcpy(cc,b(2*i-1),ldb,1,w(1+n),1,n,nbi,n)
C     call zprm('this is w',2,w(1+n),n,n,nbi)
C     Overwrite w with aT^-1 bT
            if (lherm) then
               call zhetrs('U',n,nbi,a,lda,w,w(1+n),n,ierr)
            else
               call zgetrs('T',n,nbi,a,lda,w,w(1+n),n,ierr)
            endif
            if (ierr .ne. 0) return
C     Copy transpose of w to subblock of b
            call zmcpy(cc,w(1+n),1,n,b(2*i-1),ldb,1,nbi,n)
         enddo
      elseif (lherm) then
         call zhetrs('U',n,nb,a,lda,w,b,ldb,ierr)
      else
         call zgetrs('N',n,nb,a,lda,w,b,ldb,ierr)
      endif
      deallocate(w)
      if (ierr .ne. 0) return
      end

