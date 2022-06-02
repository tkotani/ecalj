subroutine polint(xa,ya,n,npoly,x,dymx,iopt,jx,y,dy)
  !- Polynomial interpolation of tabulated data at a specified point
  ! ----------------------------------------------------------------
  !i Inputs
  !i    xa,ya,n: table of x-y pairs and number
  !i    npoly: maximum order of polynomial
  !i    x: point at which to evaluate function
  !i    dymx: return when estimated value of dy is less than dymx
  !i    jx:  initial guess for xa(jx) that brackets x for
  !i    iopt: not used; should be zero
  !o Outputs
  !o    y: interpolated value of ya at x
  !o    dy:error estimate
  !o    jx:  calculated value for xa(jx) that brackets x
  !r Remarks
  !     Adapted from Numerical Recipes
  !     Bracketing needs work if points irregularly spaced.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,npoly,iopt,jx
  double precision :: dy,dymx,x,y,xa(1),ya(1)
  ! Local variables
  integer :: np,j1

  np = min(npoly,n)
  call hunty(xa,n,np,x,jx,j1)
  jx = min(max(1,jx-np/2),n-np+1)
  !      if (iprint() .ge. 100) then
  !        print 333, n,np,jx,x,xa(n),ya(n)
  !  333   format(' n,np,jx,x,xa(n),ya(n)=',3i4,3g12.6)
  !        pause
  !      endif
  call polinx(xa(j1),ya(j1),min(npoly,n-j1+1),x,dymx,y,dy)
end subroutine polint
subroutine polinx(xa,ya,n,x,dymx,y,dy)
  !- Kernel called by polint
  ! ----------------------------------------------------------------
  !i Inputs
  !i    xa,ya,n: table of x-y pairs and number
  !i    x: point at which to evaluate function
  !i    dymx: return when estimated value of dy is less than dymx
  !o Outputs
  !o    y: interpolated value of ya at x
  !o    dy:error estimate
  !r Remarks
  !     Uses Neville algorithm (adapted from Numerical Recipes)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  double precision :: dy,dymx,x,y,xa(1),ya(1)
  integer :: n
  ! Local variables
  double precision :: den,dif,dift,ho,hp,w
  integer :: i,m,nmax,ns
  parameter (nmax=50)
  double precision :: c(nmax),d(nmax)

  ns = 1
  dif = dabs(x-xa(1))
  ! --- Find the index ns of the closest table entry ---
  do  11  i = 1, n
     dift = dabs(x-xa(i))
     if (dift < dif) then
        ns = i
        dif = dift
     endif
     c(i) = ya(i)
     d(i) = ya(i)
11 enddo

  y = ya(ns)
  ns = ns-1
  ! --- For each column of the Neville recursion update C and D ---
  do  13  m = 1, n-1
     do  12  i = 1, n-m
        ho = xa(i)-x
        hp = xa(i+m)-x
        w  = c(i+1)-d(i)
        den = ho-hp
        if (den == 0d0) call rx('polinx: den=0')
        den = w/den
        d(i) = hp*den
        c(i) = ho*den
12   enddo
     ! ---   Decide which whether to accumulate correction C or D to y ---
     if (2*ns < n-m) then
        dy = c(ns+1)
     else
        dy = d(ns)
        ns = ns-1
     endif
     y = y+dy
     if (dabs(dy) < dymx) return
13 enddo
end subroutine polinx
SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
  implicit double precision (a-h,p-z)
  integer:: nmax,i,m,ns,n
  PARAMETER (NMAX=50,TINY=1.d-25)
  double precision :: XA(N),YA(N),C(NMAX),D(NMAX),xx
  NS=1
  HH=dABS(X-XA(1))
  DO 11 I=1,N
     H=dABS(X-XA(I))
     IF (H == 0.d0)THEN
        Y=YA(I)
        DY=0.0d0
        RETURN
     ELSE IF (H < HH) THEN
        NS=I
        HH=H
     ENDIF
     C(I)=YA(I)
     D(I)=YA(I)+TINY
11 enddo
  Y=YA(NS)
  NS=NS-1
  DO 13 M=1,N-1
     DO 12 I=1,N-M
        W=C(I+1)-D(I)
        H=XA(I+M)-X
        T=(XA(I)-X)*D(I)/H
        DD=T-C(I+1)
        ! akao
        !          IF(DD.EQ.0.d0)PAUSE
        IF(DD == 0.d0) read(5,*) xx
        DD=W/DD
        D(I)=C(I+1)*DD
        C(I)=T*DD
12   enddo
     IF (2*NS < N-M)THEN
        DY=C(NS+1)
     ELSE
        DY=D(NS)
        NS=NS-1
     ENDIF
     Y=Y+DY
13 enddo
END SUBROUTINE RATINT
subroutine hunty(xa,n,m,x,ix,low)
  !- Finds points nearest a value within an ordered array of points
  ! ----------------------------------------------------------------
  !i    xa,n: array of points and number
  !i    x: value to bracket
  !i    ix: initial guess for output ix
  !i    m: number of points to find closest to x
  !o    ix: xa(ix) < x < xa(ix+1); low: xa(low) ... xa(low+m) nearest to x
  !r    Adapted from Numerical Recipes
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,m,ix,low
  double precision :: xa(1),x
  ! Local variables
  integer :: inc,jhi,jm
  logical :: ascnd

  ascnd = xa(n) .gt. xa(1)
  if (ix <= 0 .OR. ix > n) then
     ix = 0
     jhi = n+1
     goto 3
  endif
  inc = 1
  if (x >= xa(ix) .eqv. ascnd) then
1    jhi = ix+inc
     if (jhi > n) then
        jhi = n+1
     else if (x >= xa(jhi) .eqv. ascnd) then
        ix = jhi
        inc = inc+inc
        goto 1
     endif
  else
     jhi = ix
2    ix = jhi-inc
     if (ix < 1) then
        ix = 0
     else if (x < xa(ix) .eqv. ascnd) then
        jhi = ix
        inc = inc+inc
        goto 2
     endif
  endif
3 if (jhi-ix == 1) goto 10
  jm = (jhi+ix)/2
  if (x > xa(jm) .eqv. ascnd) then
     ix = jm
  else
     jhi = jm
  endif
  goto 3

  ! --- Given ix, ix+1 that bracket root, find low ---
10 continue
  jhi = (min(ix+1,n))
  if (dabs(xa(max(ix,1))-x) < dabs(xa(min(ix+1,n))-x)) jhi  = ix
  low = jhi
  do  12  jm = 2, m
     if (low > 1) then
        if (dabs(xa(max(low-1,1))-x) < dabs(xa(min(jhi+1,n))-x) &
             .OR. jhi == n) then
           low = low-1
        else
           jhi = jhi+1
        endif
     else
        jhi = min(jhi+1,n)
     endif
12 enddo
end subroutine hunty
! #if TESTP
! subroutine fmain
!   implicit double precision (a-h,p-z)
!   PARAMETER(NP=10)
!   double precision :: XA(NP),YA(NP)
!   pi = datan(1d0)
!   WRITE(*,*) 'Generation of interpolation tables'
!   WRITE(*,*) ' ... sin(x)    0<x<pi'
!   WRITE(*,*) ' ... exp(x)    0<x<1 '
!   WRITE(*,*) 'dymx, max order of interpolation (<=10)?'
!   read(*,*) dymx,nn
!   n = 10
!   DO 14 NFUNC=1,2
!      jx = 1
!      IF (NFUNC == 1) THEN
!         WRITE(*,*) 'sine function from 0 to pi'
!         DO 11 I=1,N
!            XA(I)=I*PI/N
!            YA(I)=dSIN(XA(I))
! 11      enddo
!      ELSE IF (NFUNC == 2) THEN
!         WRITE(*,*) 'exponential function from 0 to 1'
!         DO 12 I=1,N
!            XA(I)=I*1.0d0/N
!            YA(I)=dEXP(XA(I))
! 12      enddo
!      ELSE
!         STOP
!      ENDIF
!      WRITE(*,'(T10,A1,T20,A4,T28,A12,T46,A5,T58,A15)') &
!           'x','f(x)','interpolated','error','error est    jx'
!      jx = 0
!      DO 13 I=1,10
!         IF (NFUNC == 1) THEN
!            X=(-0.05d0+I/10.0d0)*PI
!            F=dSIN(X)
!         ELSE IF (NFUNC == 2) THEN
!            X=(-0.05d0+I/10.0d0)
!            F=dEXP(X)
!         ENDIF
!         call polint(xa,ya,n,nn,x,dymx,0,jx,y,dy)
!         WRITE(*,'(1x,3f12.6,2e15.4,i5)') x,f,y,f-y,dy,jx
! 13   enddo
!      pause
! 14 enddo
! end subroutine fmain
! #endif
! #if TESTR
! subroutine fmain
!   implicit double precision (a-h,p-z)
!   PARAMETER(NPT=6,EPSSQ=1.0d0)
!   double precision :: X(NPT),Y(NPT)
!   F(X)=X*dEXP(-X)/((X-1.0d0)**2+EPSSQ)
!   DO 11 I=1,NPT
!      X(I)=I*2.0d0/NPT
!      Y(I)=F(X(I))
! 11 enddo
!   WRITE(*,'(/1X,A/)') 'Diagonal rational function interpolation'
!   WRITE(*,'(1X,T6,A,T13,A,T26,A,T40,A,T53,A)') &
!        'x','interp.','est err ','actual','true err'
!   DO 12 I=1,10
!      XX=0.2d0*I
!      CALL RATINT(X,Y,NPT,XX,YY,DYY)
!      YEXP=F(XX)
!      WRITE(*,'(1X,F6.2,F12.6,E15.4,F12.6,E15.4,)') &
!           XX,YY,DYY,YEXP,yexp-yy
! 12 enddo
! end subroutine fmain
! #endif

