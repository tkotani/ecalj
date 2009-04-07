      integer function falsi(ax,bx,fa,fb,f,cof,ir,tol,itmax,x,fx)
C- Find a root by a modified regula falsi
C A modified version of the regula falsi method (secant method) to
C converge to a root, already assumed bracketed.  Falsi uses a
C simple scaling to avoid the usual convergence problems a highly
C nonlinear function near a root.  This scaling artifically divides
C one of the dependent variables by a scale factor.  If this
C point is discarded in the next iteration, the scaling has no
C effect, so the scaling only takes effect when one point is
C repeatedly thrown away.
C Inputs:
C AX and BX are two values of the independent variable;
C FA and FB are the associated function values (precomputed);
C F is the function name;
C TOL is the largest approximate relative error (can be zero); 
C ITMAX is the largest number of function iterations before FALSI quits
C Outputs
C X and FX are the best estimates of the root and associated value;
C ITER the number of iterations required to reach the root.

C     implicit none
      double precision ax,bx,fa,fb,f,fx,tol,cof
      double precision x,scale,eps,d1mach
      integer iter,itmax,ir,iprint
      external f

C EPS is approximately the machine precision times a number that is
C representative of the scale of the X coordinate, estimated by BX-AX
      eps = 2*d1mach(3)*(bx-ax)
      if (tol .ne. 0) eps = tol
      scale = 1
      do  10  iter = 1, itmax

C ...   debugging
        if (iprint() .ge. 100) then
          print 333, ax,bx,fa,fb
  333     format('falsi: ax,bx=',2f12.6,' fa,fb=',2f12.6)
        endif

C ...   New iterate for function
        x = (ax*fb - bx*fa)/(fb-fa)
        fx = f(x,cof,ir)
        falsi = iter
        if (fx .eq. 0d0 .or. dabs(bx-ax) .lt. eps) return

C ...   Decide which bound to update
        if (dsign(fx,fx) .eq. dsign(fx,fa)) then
          scale=scale*2
          fb = fb/scale
        else
          bx = ax
          fb = fa
          scale = 1
        endif
        ax = x
        fa = fx
   10 continue
      falsi = -iter

      end

