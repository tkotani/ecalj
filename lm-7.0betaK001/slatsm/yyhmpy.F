#define BLAS3 1
      subroutine yyhmpy(ta,tb,n,l,a,b,lc,c)
C- Complex matrix multiply, result assumed hermitian
C ----------------------------------------------------------------
Ci Inputs:
Ci   ta,tb: follow BLAS3 conventions
Ci   n,l: dimension of c and length of vector product
Ci   lc:  true: lower triangle copied into upper triangle
Co Outputs:
Co   c
Cr Remarks:
Cr   Adapted from zhmpy, alternatively can call zampy in block form
Cm Memory
Cm   zampy call requires 2*nrow*l double precision words.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*1 ta,tb
      logical lc
      integer n,l
      double precision a(1), b(1), c(n,1)
C Local variables
      integer nrow,owk,ir,ic,nr,nc,ns,n2
      logical ca,cb,tra,trb
#if BLAS3
      parameter (nrow=48)
      integer w(1)
      common /w/ w
#else
      integer la,lb,nra,nrb
#endif

      ca = ta.eq.'c' .or. ta.eq.'C'
      cb = tb.eq.'c' .or. tb.eq.'C'
      tra = ca .or. ta .eq. 't' .or. ta .eq. 'T'
      trb = cb .or. tb .eq. 't' .or. tb .eq. 'T'

#if BLAS3
      n2 = n**2
      ns = n*l
      call defdr(owk,nrow*l*2)
      do  10  ir = 1, n, nrow
        nr = min(n-ir+1,nrow)
        nc = nr+ir-1
        if (tra) then
          call xyhmpy(nr,l,ns,a(1+(ir-1)*l),ca,w(owk))
          call zampy(w(owk),nr,1,l*nr,b,l,1,ns,c(ir,1),n,1,n2,nr,nc,l)
        else
          call zampy(a(ir),n,1,ns,b,l,1,ns,c(ir,1),n,1,n2,nr,nc,l)
        endif
   10 continue

      if (.not. lc) return
      do  12  ir = 1, n
        do  12  ic = ir+1, n
          c(ir,ic)    =  c(ic,ir)
          c(ir+n2,ic) = -c(ic+n2,ir)
   12   continue

#else
      la = n
      nra = 1
      if (tra) then
        nra = l
        la = 1
      endif
      lb = l
      nrb = 1
      if (trb) then
        nrb = n
        lb = 1
      endif
      if (ca) call dscal(n*l,-1d0,a(1+n*l),1)
      if (cb) call dscal(n*l,-1d0,b(1+n*l),1)
      call zhmpy(a,la,nra,n*l,b,lb,nrb,n*l,c,n,lc,l)
      if (ca) call dscal(n*l,-1d0,a(1+n*l),1)
      if (cb) call dscal(n*l,-1d0,b(1+n*l),1)
#endif

      end
#if BLAS3
      subroutine xyhmpy(nr,l,ns,a,ca,w)
C Kernel called by yyhmpy
C     implicit none
      integer nr,l,ns
      double precision w(nr,l,2),a(l,1)
      integer j
      logical ca

      call dpzero(w(1,1,2),nr*l)
      do  10  j = 1, nr
        call dcopy(l,a(1,j),1,   w(j,1,1),nr)
        call dcopy(l,a(1+ns,j),1,w(j,1,2),nr)
   10 continue
      if (ca) call dscal(nr*l,-1d0,w(1,1,2),1)

C      call zprm('(5f15.10)',w,nr,l)
      end
#endif

