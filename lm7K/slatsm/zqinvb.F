      subroutine zqinvb(cs,a,lda,n,nb,w,ldw,w2,b,ldb,ierr)
C- Solution of a x = b by vectorizable multiplications and inversions
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
Ci          'l'  Always call Lapack routines
Ci   a     :lhs of equation a x = b
Ci   lda   :leading dimension of ar,ai
Ci   n     :solve a x = b for matrix a(1:n,1:n)
Ci   w,ldw :double precision work array w(ldw,*), and leading dimension.
Ci         :The required size of w depends on which branch is executed.
Ci         :*If Lapack not used: w = w(ldw,n+1) with ldw>=n.
Ci         :*If Lapack is used, w = w(ldw,n); the req'd and optimal
Ci         :    sizes of ldw depend on the particular branch.
Ci         :   *If cs contains no 't' or 'h':
Ci              required ldw>=1; no optimal size
Ci         :   *If cs contains no 't':
Ci              required ldw>=1; optimal ldw=129
Ci         :   *If cs contains 't':
Ci              required ldw>=3; optimal ldw=129
Ci   w2    :a double precision work array of dimension nb*(n+1)*2
Ci         :*If Lapack is used: w2 is not used
Ci         :*If Lapack is not used, w2 has dimensions nb*(n+1)*2
Ci         :    but w and w2 may use the same address space.
Ci   b     :right hand side of equation a x = b
Ci   ldb   :leading dimension of b
Ci   nb    :the number of columns (rows, if cs contains 't') to solve.
Co Outputs:
Co   a     :is OVERWRITTEN, into a partially decomposed form
Co   ierr  :is returned nonzero if matrix was not successfully inverted
Co   b     :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
Cb Bugs
Cb   yyqnvb operates recursively, and will fail if the smallest matrix
Cb   subblock of  a  is singular, even if  a  is not.
Cb   If Lapack is used, yyqnvb is not called.
Cr Remarks:
Cr   zqinvb uses one of two branches to solve the linear equation
Cr   One branch which makes lapack calls; the other calls yyqnvb.
Cr
Cr   The Lapack branch is recommended: it appears to be faster
Cr   and it is more accurate.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,lda,ldw,ldb,ierr,nb
      character cs*(*)
      double precision a(lda),b(ldb),w(ldw*n),w2(nb*(n+1))
C ... Local parameters
      integer i,j,lwrk,ncut,nbi
      logical lapack,ltrns,lsx,lherm
      character*1 cc
      parameter (ncut=32)

      ltrns = .false.
      lsx   = .false.
      lherm = .false.
C     cm = 'N'
      lapack = .false.
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
        elseif (cs(i:i) .eq. 'l') then
          lapack = .true.
        elseif (cs(i:i) .ne. ' ') then
          call rxs2('zqinv: unrecognized switch, cs=''',cs,'''')
        endif
    2 continue

      lapack = n .le. ncut .or. lapack

C --- Lapack branch ---
      if (lapack) then
        ierr = 0
C       LU decomposition
        if (lherm) then
C         Use n elements for pivot; remainder is complex work
          lwrk = (ldw*n - n)/2
          ierr = 0
          if (.not. lsx) call zhetrf('U',n,a,lda,w,w(1+n),lwrk,ierr)
          cc = 'C'
        else
          if (.not. lsx) call zgetrf(n,n,a,lda,w,ierr)
          cc = 'N'
        endif
        if (ierr .ne. 0) return
C       Back substitution
        if (ltrns) then
C         Use n elements for pivot; remainder is complex work
          lwrk = (ldw - 1)/2
          do  i = 1, nb, lwrk
            nbi = min(nb-i+1,lwrk)
C           Copy subblock of b to transpose of w
            call zmcpy(cc,b(2*i-1),ldb,1,w(1+n),1,n,nbi,n)
C           call zprm('this is w',2,w(1+n),n,n,nbi)
C           Overwrite w with aT^-1 bT
            if (lherm) then
              call zhetrs('U',n,nbi,a,lda,w,w(1+n),n,ierr)
            else
              call zgetrs('T',n,nbi,a,lda,w,w(1+n),n,ierr)
            endif
            if (ierr .ne. 0) return
C           Copy transpose of w to subblock of b
            call zmcpy(cc,w(1+n),1,n,b(2*i-1),ldb,1,nbi,n)
          enddo
C         call zprm('this is b',2,b,ldb,nb,n)
        elseif (lherm) then
          call zhetrs('U',n,nb,a,lda,w,b,ldb,ierr)
        else
          call zgetrs('N',n,nb,a,lda,w,b,ldb,ierr)
        endif
        if (ierr .ne. 0) return

C --- yyqnvb branch ---
      else
        call ztoy(a,lda,n,n,0)
        if (ltrns)       call ztoy(b,ldb,nb,n,0)
        if (.not. ltrns) call ztoy(b,ldb,n,nb,0)
C       call yyprm('a in zqinvb',2,a,lda,2*lda,n,n)
C       if (cs(1:1).ne.'t') call yyprm('b in zqinvb',2,b,ldb,2*ldb,n,nb)
C       if (cs(1:1).eq.'t') call yyprm('t in zqinvb',2,b,ldb,2*ldb,nb,n)
        call yyqnvb(cs,a,a(1+lda),2*lda,n,nb,w,ldw,w2,b,b(1+ldb),2*ldb,
     .  ierr)
C       if (cs(1:1).ne.'t') call yyprm('a^-1 b',2,b,ldb,2*ldb,n,nb)
C       if (cs(1:1).eq.'t') call yyprm('b a^-1',2,b,ldb,2*ldb,nb,n)
        if (ierr .ne. 0) return
        call ztoy(a,lda,n,n,1)
        if (ltrns)       call ztoy(b,ldb,nb,n,1)
        if (.not. ltrns) call ztoy(b,ldb,n,nb,1)
      endif

      end

