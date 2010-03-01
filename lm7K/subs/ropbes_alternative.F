C ----------------------------------------------------------------
C   Bessel function at r(i)*sqrt(e) divided by r**l 
c
Ci Inputs
Ci   r    list of points
Ci   e    energy
Ci   y,h  work vectors of length n each (not used)
Co Outputs
Co   xi   J(r,l)/r**l, according to standard definition
C
c origianl: Feb. 15, 2010, Hiori Kino
      subroutine ropbes(r,e,lmax,y,h,xi,n)
      implicit none
      integer lmax,n
      double precision e,r(n),xi(n,0:lmax),h(n),y(n)

      double precision:: e2,f(0:lmax),eps=1.0d-10,x
      integer:: i,l 

      e2=sqrt(e)
      do i=1,n
         x=r(i)*e2
        if (x< eps)  x=x+eps 
        call bessel(lmax,x,f,y)
        do l=0,lmax
          f(l)=f(l)/(x**l)
        enddo
        xi(i,0:lmax)= f(0:lmax)
      enddo   
 
      end subroutine ropbes
c----------------------------------------------------
#if 1
      subroutine bessel_yl(lmax,x0,f,d)
      implicit none
      integer:: lmax
      real(8):: f(0:lmax), d(0:lmax), x0
      real(8):: x,eps=1.0d-10
      integer:: l
      x=x0
      if ( x0 < eps ) x=x0+eps 
      call bessel(lmax,x,f,d)
      do l=0,lmax
        f(l)=f(l)/(x**l)
      enddo 
      end subroutine bessel_yl
#endif 
c----------------------------------------------------
C---- BESSEL
C     SPHERICAL BESSEL FUNCTION AND ITS DERIVATIVE
C     CODED BY T.OGUCHI  10/01/87
C
C     LMAX : MAXIMUN ANGULAR MOMENTUM
C     X    : ARGUMENT OF BESSEL FUNCTION
C     F    : SPHERICAL BESSEL FUNCTION
C     D    : DERIVATIVE OF SPHERICAL BESSEL FUNCTION
C
      subroutine bessel(lmax,x,f,d)
c                           @(#)bessel.f 9.1 97/05/08 14:47:47 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(0:LMAX),D(0:LMAX)
      DIMENSION Q(0:50)
      IF (LMAX.GT.50) STOP '=== STOP AT BESSEL ==='
      IF(X.LT.1.0D-3) THEN
        XN=1.0
        F1=1.0
        F2=3.0D0
        F(0)=XN/F1-0.5D0*XN*X*X/F2
        DO 4000 L=1,LMAX
        F1=F2
        F2=F2*(2.0D0*L+3.0D0)
        Q(L)=XN/F1-0.5D0*XN*X*X/F2
        F(L)=Q(L)*X
        XN=XN*X
 4000   CONTINUE
      ELSEIF(X.LT.10.0D0) THEN
        FVAL=1.0D-30
        FPLUS=0.0
        FSUM=0.0
        DO 2000 L=LMAX+20,LMAX+2,-1
        FMIN=FVAL*(2.0*L+1.0)/X-FPLUS
        FSUM=FSUM+(2.0*L-1.0)*FMIN*FMIN
        FPLUS=FVAL
        FVAL=FMIN
        IF(ABS(FMIN).GT.1.0D25) THEN
          FSUM=1.0D0/SQRT(FSUM)
          FPLUS=FPLUS*FSUM
          FVAL=FVAL*FSUM
          FSUM=1.0D0
        ENDIF
 2000   CONTINUE
        IF(X.GT.1.0D0) THEN
          XL=1.0D0
        ELSE
          XL=X**LMAX
        ENDIF
        FSUM=XL/SQRT(FSUM)
        FPLUS=FPLUS*FSUM
        FVAL=FVAL*FSUM
        FSUM=XL*XL
        DO 2100 L=LMAX+1,1,-1
        FMIN=FVAL*(2.0*L+1.0)/X-FPLUS
        FSUM=FSUM+(2.0*L-1.0)*FMIN*FMIN
        FPLUS=FVAL
        FVAL=FMIN
        F(L-1)=FVAL
 2100   CONTINUE
        FSUM=1.0/SQRT(FSUM)
        DO 3000 L=0,LMAX
        F(L)=F(L)*FSUM
        Q(L)=F(L)/X
 3000   CONTINUE
      ELSE
        F(0)=SIN(X)/X
        F(1)=SIN(X)/X/X-COS(X)/X
        Q(1)=F(1)/X
        DO 5000 L=1,LMAX-1
        F(L+1)=(2.0D0*L+1.0D0)*F(L)/X-F(L-1)
        Q(L+1)=F(L+1)/X
 5000   CONTINUE
      ENDIF
      D(0)=-F(1)
      DO 4100 L=1,LMAX
      D(L)=F(L-1)-(L+1.0D0)*Q(L)
 4100 CONTINUE
      RETURN
      END
