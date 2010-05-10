c takao modified this based on delw_2.2.
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
      integer:: nrow , ir , ic , nr , nc , ns , n2 
      real(8) ,allocatable :: rv_w_owk(:)

      logical ca,cb,tra,trb
      parameter (nrow=48)

      ca = ta.eq.'c' .or. ta.eq.'C'
      cb = tb.eq.'c' .or. tb.eq.'C'
      tra = ca .or. ta .eq. 't' .or. ta .eq. 'T'
      trb = cb .or. tb .eq. 't' .or. tb .eq. 'T'

      n2 = n**2
      ns = n*l
      allocate(rv_w_owk(nrow*l*2))

      do  10  ir = 1, n, nrow
        nr = min(n-ir+1,nrow)
        nc = nr+ir-1
        if (tra) then
           call xyhmpy( nr , l , ns , a ( 1 + ( ir - 1 ) * l ) , ca , rv_w_owk ) 
           call zampy ( rv_w_owk , nr , 1 , l * nr , b , l , 1 , ns , c 
     .     ( ir , 1 ) , n , 1 , n2 , nr , nc , l ) 
        else
          call zampy(a(ir),n,1,ns,b,l,1,ns,c(ir,1),n,1,n2,nr,nc,l)
        endif
   10 continue
      if ( .not. lc ) then 
         if (allocated(rv_w_owk)) deallocate(rv_w_owk)
         return
      endif 

      do  12  ir = 1, n
        do  12  ic = ir+1, n
          c(ir,ic)    =  c(ic,ir)
          c(ir+n2,ic) = -c(ic+n2,ir)
   12   continue

      if (allocated(rv_w_owk)) deallocate(rv_w_owk)
      end

      subroutine xyhmpy(nr,l,ns,a,ca,w)
C Kernel called by yyhmpy
C     implicit none
      integer nr,l,ns
      double precision w(nr,l,2),a(l,1)
      integer j
      logical ca
c      call dpzero(w(1,1,2),nr*l)
      w(:,:,2)=0d0
      do  10  j = 1, nr
        call dcopy(l,a(1,j),1,   w(j,1,1),nr)
        call dcopy(l,a(1+ns,j),1,w(j,1,2),nr)
   10 continue
      if (ca) call dscal(nr*l,-1d0,w(1,1,2),1)
      end

