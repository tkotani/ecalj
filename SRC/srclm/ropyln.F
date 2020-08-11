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

C     call tcn('ropyln')
      if (n .gt. nd) call rx('ropyln: nd too small')
      fpi = 16*datan(1d0)
      allocate(cm_rv(n),sm_rv(n),q_rv(n*2),h_rv(n))

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
     .    , nd , yl )
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
c$$$#if TEST
c$$$C Test program to check ropyln
c$$$      subroutine fmain
c$$$      implicit none
c$$$      integer nrx,nlmx,nr,lmax,nlm1,ir,ii,l,ilm,i1,i2,nsize
c$$$      parameter (nrx=20,nlmx=49,nsize=100000)
c$$$      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),
c$$$     .ylv(nrx,nlmx),ylok(nrx,nlmx),dr(3),tops,ylm(nlmx)
c$$$
c$$$      call wkinit(nsize)
c$$$      call wkfast(.true.)
c$$$      call sylmnc(cy,15)
c$$$
c$$$      lmax = 2
c$$$   99 print *, 'lmax='
c$$$      read(*,*) lmax
c$$$
c$$$      call makr(0d0,nr,x,y,z)
c$$$
c$$$C ... Make nonvectorized ylm's
c$$$      nlm1 = (lmax+1)**2
c$$$      do  ir = 1, nr
c$$$        dr(1) = x(ir)
c$$$        dr(2) = y(ir)
c$$$        dr(3) = z(ir)
c$$$        call sylm(dr,ylm,lmax,r2)
c$$$        do  ilm = 1, nlm1
c$$$          ylok(ir,ilm) = cy(ilm)*ylm(ilm)
c$$$        enddo
c$$$      enddo
c$$$C     Test: Y_1-1 = sqrt(3/4/pi) y
c$$$C     print *, y(1) * dsqrt(0.75/4/atan(1d0))
c$$$C     print *, ylok(1,2)
c$$$
c$$$      call ropyln(nr,x,y,z,lmax,nrx,ylv,r2)
c$$$
c$$$      tops = 0
c$$$      do  10  ir = 1, nr
c$$$        do  12  l = 0, lmax
c$$$          i1 = l*l+1
c$$$          i2 = (l+1)**2
c$$$          print 333, (ylok(ir,ii),ii=i1,i2)
c$$$   12   continue
c$$$  333   format(9f8.5)
c$$$        print *
c$$$        do  14  l = 0, lmax
c$$$          i1 = l*l+1
c$$$          i2 = (l+1)**2
c$$$          do  16  ii = i1, i2
c$$$            tops = max(tops,dabs(ylok(ir,ii)-ylv(ir,ii)))
c$$$   16     continue
c$$$          print 333, (ylok(ir,ii)-ylv(ir,ii),ii=i1,i2)
c$$$   14   continue
c$$$        print *, '----------------'
c$$$   10 continue
c$$$
c$$$      print 335, tops
c$$$  335 format(' max error for ylm:',1pe12.2)
c$$$      end
c$$$      subroutine makr(rsm,nr,x,y,z)
c$$$      implicit none
c$$$      integer nr,i,ir
c$$$      double precision rs,rsm,x(1),y(1),z(1)
c$$$      real ran1
c$$$      rs = rsm
c$$$      if (rsm .lt. 1d-9) rs = .5d0
c$$$      call ran1in(1)
c$$$      nr = 5
c$$$      do  10  i = 1, nr
c$$$        ir = i+1
c$$$        x(i) = abs((ran1()-.5d0)*5*rs)
c$$$        y(i) = (ran1()-.5d0)*5*rs
c$$$        z(i) = (ran1()-.5d0)*5*rs
c$$$   10 continue
c$$$
c$$$      end
c$$$#endif

