      subroutine ropyln(n,x,y,z,lmax,nd,yl,rsq)
C- Normalized spheric harmonic polynomials (vectorizes).
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of points for which to calculate yl
Ci   x     :x component of Cartesian coordinate
Ci   y     :y component of Cartesian coordinate
Ci   z     :z component of Cartesian coordinate
Ci   lmax  :compute yl for l=0..lmax
Ci   nd    :leading dimension of yl; nd must be >= n
Co Outputs
Co   yl    :Ylm(i,ilm) are the (real) spherical harmonics
Co         :for l=0..lmax and point i.
Co   rsq   :rsq(i) square of length of point i
Cr Remarks
Cr   yl = real harmonics (see Takao's GW note) * r^l
Cu Updates
Cu  25 Jun 03 (Kino) initialize cx to zero
C ----------------------------------------------------------------------
C     implicit none
      integer:: nd , n , i , m , lmax , l , kk
      real(8) ,allocatable :: cm_rv(:)
      real(8) ,allocatable :: sm_rv(:)
      real(8) ,allocatable :: q_rv(:)
      real(8) ,allocatable :: h_rv(:)

      double precision x(1),y(1),z(1),yl(nd,1),rsq(1),cx(3)
      double precision fpi,f2m
      integer w(1)
      common /w/ w

C     call tcn('ropyln')
      if (n .gt. nd) call rx('ropyln: nd too small')
      fpi = 16*datan(1d0)
      allocate(cm_rv(n))
Ctakao_ZeroClear_NotRequiered        if (n<0) cm_rv(:)=0.0d0

      allocate(sm_rv(n))
Ctakao_ZeroClear_NotRequiered        if (n<0) sm_rv(:)=0.0d0

      allocate(q_rv(n*2))
Ctakao_ZeroClear_NotRequiered        if (n*2<0) q_rv(:)=0.0d0

      allocate(h_rv(n))
Ctakao_ZeroClear_NotRequiered        if (n<0) h_rv(:)=0.0d0

      do  2  i = 1, n
        rsq(i) = x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
    2 continue
      cx(1) = 0
      cx(2) = 0
      cx(3) = 0

C --- Loop over m: cx are normalizations for l, l-1, l-2 ---
      f2m = 1d0
      do  10  m = 0, lmax
        call ropcsm ( m , n , x , y , h_rv , cm_rv , sm_rv )

        if (m .eq. 0) then
          cx(1) = dsqrt(1/fpi)
        else
          f2m = f2m*2*m*(2*m-1)
          cx(1) = dsqrt((2*m+1)*2/fpi/f2m)
        endif
        do  11  l = m, lmax
          call ropqln ( m , l , n , rsq , z , cx , q_rv , kk )

          call ropynx ( m , l , kk , n , q_rv , cm_rv , sm_rv 
     .   , nd , yl )

          cx(3) = cx(2)
          cx(2) = cx(1)
          cx(1) = cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
   11   continue
   10 continue
      if (allocated(h_rv)) deallocate(h_rv)
      if (allocated(q_rv)) deallocate(q_rv)
      if (allocated(sm_rv)) deallocate(sm_rv)
      if (allocated(cm_rv)) deallocate(cm_rv)

C     call tcx('ropyln')
      end
C     Separate subroutines below to avoid problems with some
C     optimizing compilers.
C      subroutine ropqln(m,l,n,r2,z,cx,q,kk)
CC- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
Cc  Returns kk, which points to the current component of q.
C      implicit none
C      integer mm,n,i,l,m,kk,k2,k1
C      double precision q(n,2),r2(n),z(n),cx(3)
C      double precision a,b,xx,yy
C
CC --- Case l=m ---
C      if (l .eq. m) then
C        a = 1d0
C        do  1  mm = 0, m-1
C    1   a = a*(2*mm+1)
C        kk = 1
C        a = a*cx(1)
C        do  2  i = 1, n
C    2   q(i,kk) = a
C        return
C      endif
C
CC --- Case l=m+1 ---
C      if (l .eq. m+1) then
C        b = 1d0
C        do  3  mm = 0, m
C    3   b = b*(2*mm+1)
C        b = b*cx(1)
C        kk = 2
C        do  4  i = 1, n
C    4   q(i,kk) = b*z(i)
C        return
C      endif
C
CC --- Case l=m+2 and higher by recursion ---
C      if (l .ge. m+2) then
C        k2 = kk
C        k1 = kk+1
C        if (k1 .eq. 3) k1 = 1
C        xx = -(l+m-1d0)/(l-m)*cx(1)/cx(3)
C        yy = (2*l-1d0)/(l-m)*cx(1)/cx(2)
C        do  6  i = 1, n
C    6   q(i,k1) = xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
C        kk = k1
C        return
C      endif
C      end
C      subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
C      implicit none
C      integer lav,n,nd,l,i,m,kk
C      double precision q(n,2),cm(n),sm(n),yl(nd,1)
C      lav = l*(l+1)+1
C      do  1  i = 1, n
C    1 yl(i,lav+m) = cm(i)*q(i,kk)
C      if (m .eq. 0) return
C      do  2  i = 1, n
C    2 yl(i,lav-m) = sm(i)*q(i,kk)
C      end
#if TEST
C Test program to check ropyln
      subroutine fmain
      implicit none
      integer nrx,nlmx,nr,lmax,nlm1,ir,ii,l,ilm,i1,i2,nsize
      parameter (nrx=20,nlmx=49,nsize=100000)
      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),
     .  ylv(nrx,nlmx),ylok(nrx,nlmx),dr(3),tops,ylm(nlmx)
      integer w(1)
      common /w/ w

      call wkinit(nsize)
      call wkfast(.true.)
      call sylmnc(cy,15)

      lmax = 2
   99 print *, 'lmax='
      read(*,*) lmax

      call makr(0d0,nr,x,y,z)

C ... Make nonvectorized ylm's
      nlm1 = (lmax+1)**2
      do  ir = 1, nr
        dr(1) = x(ir)
        dr(2) = y(ir)
        dr(3) = z(ir)
        call sylm(dr,ylm,lmax,r2)
        do  ilm = 1, nlm1
          ylok(ir,ilm) = cy(ilm)*ylm(ilm)
        enddo
      enddo
C     Test: Y_1-1 = sqrt(3/4/pi) y
C     print *, y(1) * dsqrt(0.75/4/atan(1d0))
C     print *, ylok(1,2)

      call ropyln(nr,x,y,z,lmax,nrx,ylv,r2)

      tops = 0
      do  10  ir = 1, nr
        do  12  l = 0, lmax
          i1 = l*l+1
          i2 = (l+1)**2
          print 333, (ylok(ir,ii),ii=i1,i2)
   12   continue
  333   format(9f8.5)
        print *
        do  14  l = 0, lmax
          i1 = l*l+1
          i2 = (l+1)**2
          do  16  ii = i1, i2
            tops = max(tops,dabs(ylok(ir,ii)-ylv(ir,ii)))
   16     continue
          print 333, (ylok(ir,ii)-ylv(ir,ii),ii=i1,i2)
   14   continue
        print *, '----------------'
   10 continue

      print 335, tops
  335 format(' max error for ylm:',1pe12.2)
      end
      subroutine makr(rsm,nr,x,y,z)
      implicit none
      integer nr,i,ir
      double precision rs,rsm,x(1),y(1),z(1)
      real ran1
      rs = rsm
      if (rsm .lt. 1d-9) rs = .5d0
      call ran1in(1)
      nr = 5
      do  10  i = 1, nr
        ir = i+1
        x(i) = abs((ran1()-.5d0)*5*rs)
        y(i) = (ran1()-.5d0)*5*rs
        z(i) = (ran1()-.5d0)*5*rs
   10 continue

      end
#endif

