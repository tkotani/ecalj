      subroutine dlnref(lqinv,mode,nu,nv,nlma,Ai,ldi,A,lda,b,nb,cs,wk,
     .x,dx,dxmx)
C- Refine of a subblock in a solution of linear equations
C ----------------------------------------------------------------------
Ci Inputs  (see Remarks)
Ci   lqinv :specifies how the inverse Ai of A is represented.
Ci          0 Ai holds uu subblock of A, LU-decomposed by dsytrf
Ci          1 Ai holds uu subblock of A or its partially decomposed
Ci            inverse as obtained by an initial call to dqinvb.
Ci            In the latter case, caller should ensure that 'cs'
Ci            passed here has an 'b' appended (cf 'cs' in dqinvb.f).
Ci          2 Ai holds inverse of uu subblock of A
Ci   mode  :0 calculates u subblock x_u, given outer subblock x_v
Ci            or estimates x_u if x_v is approximate.  This mode is
Ci            suitable if the inverse A_uu^-1 is already accurate
Ci        <>0 refines subblock x_u given estimate x_u, x_v
Ci            This mode provides an iterative refinement to x if
Ci            Ai is approximate.
Ci   nu    :size of inner block for which beta_u is to be updated
Ci   nv    :size of inner+outer block which beta_u is to be updated
Ci   nlma  :column dimension of x
Ci   Ai,ldi:inverse of uu subblock of A, and leading dimension
Ci   A,lda :matrix for linear system of equations, and leading dimension
Ci          NB: for mode=1, uu subblock of A is not used,
Ci          in which case Ai and A may point to the same address space.
Ci   b,nb  :right-hand side for A x = b, and dimension.
Ci          NB: nb may be less then nu; dlnref assumes b(i>nb,*) = 0
Ci   cs    :string passed to dqinvb; used for lqinv=1 (cf dqinvb.f)
Ci   wk    :for lqinv=0, pivot array passed to dsytrf
Ci          for lqinv=1, work array passed to dqinvb, 
Ci   x     :estimate for solution of equations
Ci   niter :number of iterations
Co Outputs
Co   x     :updated estimate in subblock (nu+1..ndimW,1..nlma)
Co   dx    :work array holding changes in x_u
Co   dxmx  :maximum change in x (lqinv=0,1 only)
Cr Remarks
Cr   A and x are in the block form:          ( A_uu  A_uv )  ( x_u )
Cr   nu and nv are the subblock dimensions.  ( A_vu  A_vv )  ( x_v )
Cr   NB: nv may be zero.
Cr
Cr   mode = 0: dlnref solves the linear system
Cr       A x = b   =>  A_uu x_u = b - A_uv x_v
Cr
Cr   mode = 1:  solves the u subblock of the linear system A (x+dx) = 1:
Cr       A dx = A (x+dx) - b =>  A_uu dx_u = b - A_uu x_u - A_uv x_v
Cr
Cr   for dx and and adds it to x.  This mode requires an extra
Cr   multiplication, namely A_uu x_u.  But it is useful in that,
Cr   owing to numerical truncation errors that render inaccuracies
Cr   in the solution of A_uu x_u = b, it refines the estimate for x_u.
Cr Debugging: mode=0
Cr    mc -vndimW=171 -vnu=135 -vnlma=9 -qr s0a \
Cr   -split a 1,nu+1,ndimW+1 1,nu+1,ndimW+1 a11 -i -a a11i\
Cr   -i -sub 1,nr,1,nlma -split x 1,nu+1,ndimW+1 1,nc+1 -pop \
Cr   a12 x21 -x -s-1 -1:nlma -sub 1,nu,1,nlma -+ a11i -tog -x
Cr Debugging: mode=1
Cr   mc -vndimW=171 -vnu=135 -vnlma=9 -qr s0a \
Cr   -split a 1,nu+1,ndimW+1 1,nu+1,ndimW+1 a11 -i -a a11i \
Cr   -i -sub 1,nr,1,nlma -split x 1,nu+1,ndimW+1 1,nc+1 -pop \
Cr   a11i -sub 1,nu,1,nlma -a x11x\
Cr   a11 x11x -x a12 x21 -x -+\
Cr   -1:nlma -sub 1,nu,1,nlma -- a11i -tog -x
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) cs
      integer lqinv,lda,nu,nv,nb,nlma,mode,ldi
      double precision Ai(ldi,ldi),A(lda,lda),wk(*),b(nb,nlma),dxmx
      double precision x(lda,nlma),dx(lda,nlma)
C Local variables
      integer i1,i2,ii,i,j,ndef
      double precision add,omadd

C --- Set up rhs = A x - 1; store in dx ---
      i1 = 1
      if (mode .eq. 0) i1 = nu+1
      i2 = nv
      if (i2 .ge. i1) then
        call dgemm('N','N',nu,nlma,i2-i1+1,1d0,A(1,i1),lda,
     .  x(i1,1),lda,0d0,dx,lda)
      else
        do  10  j = 1, nlma
          do  10  i = 1, nu
            dx(i,j) = 0
   10   continue
      endif
      ii = min(nb,nu)
      do  12  j = 1, nlma
        do  12  i = 1, ii
          dx(i,j) = dx(i,j) - b(i,j)
   12 continue

C --- Solve A dx = rhs ---
      dxmx = 0
      add = 1d0
      if (mode .eq. 0) add = 0d0
      if (lqinv .eq. 0 .or. lqinv .eq. 1) then
        if (lqinv .eq. 0) then
          call dsytrs('L',nu,nlma,ai,ldi,wk,dx,lda,ndef)
        else
          call word(cs,1,i,j)
          if (cs(j:j) .ne. 'b')
     .    call rx('dlnref: bad input, cs missing ''b''')
          call dqinvb(cs,Ai,ldi,0,nu,nlma,wk,nu,wk,dx,lda,ndef)
        endif
        omadd = 1-add
        do  20  j = 1, nlma
          do  20  i = 1, nu
            dxmx = max(dxmx,abs(omadd*x(i,j)+dx(i,j)))
            x(i,j) = add*x(i,j) - dx(i,j)
   20   continue
      elseif (lqinv .eq. 2) then
        call dgemm('N','N',nu,nlma,nu,-1d0,Ai,ldi,dx,lda,add,x,lda)
      else
        call rxi('dlnref: bad input, lqinv=',lqinv)
      endif

C      call awrit4(' dlnref:  mode=%i with %?;n==1;partial;full;'
C     .  //' inverse.%?#n==1#   max dx = %,3;3g##',' ',80,6,
C     .  mode,lqinv,lqinv,dxmx)

c     call prmx('x',x,lda,nv,nlma)
      end

