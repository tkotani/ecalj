      subroutine ropyln(n,x,y,z,lmax,nd,yl,rsq)
C- Normalized spheric harmonic polynomials (vectorizes).
C- Also returns squares of points in rsq.
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(1),y(1),z(1),yl(nd,1),rsq(1)
c      real w(1)
      real(8),allocatable:: cx(:),wocm(:),wosm(:),woq(:),woh(:)

c      common /w/ w

      allocate(cx(3))
      cx=0d0 !kino
      if (n .le. 0) return
      fpi = 16*datan(1d0)
ccccccccccccccccccccccccc
c      call defrr(ocm,  n)
c      call defrr(osm,  n)
c      call defrr(oq,   n*2)
c      call defrr(oh,   n)
cccccccccccccccccccccc
      allocate(wocm(n),wosm(n),woq(n*2),woh(n))

      do 2 i=1,n
        rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
    2 continue
      if (n .gt. nd) call rx('ropyln: nd too small')

C --- Loop over m: cx are normalizations for l, l-1, l-2 ---
      f2m = 1d0
      do 10 m=0,lmax
        call ropcsm(m,n,x,y,woh,wocm,wosm)
        if (m .eq. 0) then
          cx(1)=dsqrt(1/fpi)
        else
          f2m=f2m*2*m*(2*m-1)
          cx(1)=dsqrt((2*m+1)*2/fpi/f2m)
        endif
        do 11 l=m,lmax
ccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) 'roplyn 4 m l= ',m,l
c      write(6,*) ' 1 cx=',cx(1:3)
ccccccccccccccccccccccccccccccccccccccccccc
          call ropqln(m,l,n,rsq,z,cx,woq,kk)
          call ropynx(m,l,kk,n,woq,wocm,wosm,nd,yl)
ccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) ' 2 cx=',cx(1:3)
ccccccccccccccccccccccccccccccccccccccccccc

c This gives strange results due to the compilar bug.
c So I reserve cx(3) by allocate.
c It makes things OK, but I don't know why.
          cx(3)=cx(2)
          cx(2)=cx(1)
          cx(1)=cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
ccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) ' 3 cx=',cx(1:3)
ccccccccccccccccccccccccccccccccccccccccccc
   11   continue
   10 continue
c      call rlse(ocm)
      deallocate(wocm,wosm,woq,woh,cx)
      end

c --- ropqln
      subroutine ropqln(m,l,n,r2,z,cx,q,kk)
C- makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
c  Returns kk, which points to the current component of q.
      implicit real*8 (a-h,p-z), integer (o)
      dimension q(n,2),r2(n),z(n),cx(3)
ccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) 'ropqln l m cx ',l,m,cx
ccccccccccccccccccccccccccccccccccccccccccc

c --- case l=m ---
      if(l.eq.m) then
        a=1d0
        do 1 mm=0,m-1
          a=a*(2*mm+1)
    1   continue
        kk=1
        a=a*cx(1)
        do 2 i=1,n
          q(i,kk)=a
    2   continue
        return
      endif
c --- case l=m+1 ---
      if(l.eq.m+1) then
        b=1d0
        do 3 mm=0,m
          b=b*(2*mm+1)
    3   continue
        b=b*cx(1)
        kk=2
        do 4 i=1,n
          q(i,kk)=b*z(i)
    4   continue
        return
      endif
c --- case l=m+2 and higher by recursion ---
      if(l.ge.m+2) then
        k2=kk
        k1=kk+1
        if(k1.eq.3) k1=1
        xx=-(l+m-1d0)/(l-m)*cx(1)/cx(3)
        yy=(2*l-1d0)/(l-m)*cx(1)/cx(2)
ccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) 'ropqln 2'
ccccccccccccccccccccccccccccccccccccccccccc
        do 6 i=1,n
          q(i,k1)=xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
    6   continue
        kk=k1
        return
      endif
      end

c --- ropynx
      subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
      implicit real*8 (a-h,p-z), integer (o)
      dimension q(n,2),cm(n),sm(n),yl(nd,1)
      lav=l*(l+1)+1
      do 1 i=1,n
        yl(i,lav+m)=cm(i)*q(i,kk)
    1 continue
      if(m.eq.0) return
      do 2 i=1,n
        yl(i,lav-m)=sm(i)*q(i,kk)
    2 continue
      end

C#ifdefC TEST
CC Test program to check ropyln
C      subroutine fmain
C      implicit none
C      integer nrx,nlmx,nr,lmax,nlm1,ir,ii,i,l,ilm,i1,i2,nsize
C      parameter (nrx=20,nlmx=49,nsize=10000)
C      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),
C     .  ylv(nrx,nlmx),ylok(nrx,nlmx),dr(3),tops
C      common /static/ cy
C
C      call wkinit(nsize)
C      call sylmnc(cy,15)
C
C      lmax = 2
C   99 print *, 'lmax='
C      read(*,*) lmax
C
C      call makr(0d0,nr,x,y,z)
C
CC ... Make ok ylm's
C      nlm1 = (lmax+1)**2
C      call ropylm(nr,x,y,z,lmax,nrx,ylok,r2)
C      do  15  ilm = 1, nlm1
C      do  15  ir = 1, nr
C   15 ylok(ir,ilm) = cy(ilm)*ylok(ir,ilm)
C
C      call ropyln(nr,x,y,z,lmax,nrx,ylv,r2)
C
C      tops = 0
C      do  10  ir = 1, nr
C        do  12  l = 0, lmax
C        i1 = l*l+1
C        i2 = (l+1)**2
C   12   print 333, (ylok(ir,ii),ii=i1,i2)
C  333   format(9f8.5)
C        print *
C        do  14  l = 0, lmax
C        i1 = l*l+1
C        i2 = (l+1)**2
C        do  16  ii = i1, i2
C   16   tops = max(tops,dabs(ylok(ir,ii)-ylv(ir,ii)))
C   14   print 333, (ylok(ir,ii)-ylv(ir,ii),ii=i1,i2)
C        print *, '----------------'
C   10 continue
C
C      print 335, tops
C  335 format(' max errors for grad h:',f12.6)
C      end
C      subroutine makr(rsm,nr,x,y,z)
C      implicit none
C      integer nr,i,ir
C      double precision rs,rsm,x(1),y(1),z(1)
C      real ran1
C      rs = rsm
C      if (rsm .lt. 1d-9) rs = .5d0
C      call ran1in(1)
C      nr = 5
C      do  10  i = 1, nr
C        ir = i+1
C        x(i) = abs((ran1()-.5d0)*5*rs)
C        y(i) = (ran1()-.5d0)*5*rs
C        z(i) = (ran1()-.5d0)*5*rs
C   10 continue
C
C      end
C#endif















