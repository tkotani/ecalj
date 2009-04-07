      subroutine platsl(r,npt)
C- Generates points corresponding platonic and related solids.
C  npt (number of points) must be 4, 6, 8, 12, 20, 30 or 60
C  30 generates points midway between the sides of the icosahedron.
C  60 generates points for the buckeyball.
C     implicit none
      integer npt
      double precision r(3,20)
      double precision d,cp,cx,pi,sx,sp,ddot,x(3,12),r2(3),wb,wi,rm(3,3)
      integer i

C --- npt 4 or 8 ---
      if (npt .eq. 4 .or. npt .eq. 8) then
        cx = 1/dsqrt(3d0)
        do  410  i = 1, 12
  410   r(i,1) = cx
        do  411  i = 1, 3
        r(i,i) = -cx
  411   r(i,4) = -cx
        if (npt .eq. 8) then
          do  810  i = 1, 12
  810     r(i,5) = -r(i,1)
        endif
        return
      endif

C --- npt 6 ---
      if (npt .eq. 6) then
        do  610  i = 1, 18
  610   r(i,1) = 0
        do  611  i = 1, 3
        r(i,i) = 1
  611   r(i,i+3) = -1
        return
      endif

      if (npt .ne. 12 .and. npt .ne. 20 .and. npt .ne. 30 .and.
     .    npt .ne. 60) call rx('platsl: bad npt')
      pi = 4*datan(1d0)
      cp = dcos(2*pi/5)
      sp = dsqrt(1-cp**2)
      cx = 1/dsqrt(5d0)
      sx = dsqrt(1-cx**2)

C --- Icosahedron ---
C ... Top 6 points
      r(1,6) = 0d0
      r(2,6) = 0d0
      r(3,6) = 1d0
      r(1,1) = sx*r(3,6)
      r(2,1) = 0d0
      r(3,1) = cx*r(3,6)
      do  1210  i = 2, 5
        r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
        r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
        r(3,i) = r(3,1)
 1210 continue

C ... Bottom 6 points of icosahedron
      r(1,12) = 0d0
      r(2,12) = 0d0
      r(3,12) =-1d0
      r(1,7) = sx*r(3,12)
      r(2,7) = 0d0
      r(3,7) = cx*r(3,12)
      do  1220  i = 8,11
        r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
        r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
        r(3,i) = r(3,7)
 1220 continue

C --- Points on the twenty faces of the icosahedron ---
      if (npt .eq. 20) then
        do  2030  i = 1, 3
 2030   r2(i) = r(i,1)+r(i,2)+r(i,6)
        d = dsqrt(ddot(3,r2,1,r2,1))
        call dpscop(r,x,36,1,1,1/d)
        do  2040  i = 1, 5
          r(1,i) = x(1,i) + x(1,1+mod(i,5)) + x(1,6)
          r(2,i) = x(2,i) + x(2,1+mod(i,5)) + x(2,6)
          r(3,i) = x(3,i) + x(3,1+mod(i,5)) + x(3,6)

          r(1,5+i) = x(1,6+i) + x(1,7+mod(i,5)) + x(1,12)
          r(2,5+i) = x(2,6+i) + x(2,7+mod(i,5)) + x(2,12)
          r(3,5+i) = x(3,6+i) + x(3,7+mod(i,5)) + x(3,12)

          r(1,10+i) = x(1,i) + x(1,1+mod(i,5)) + x(1,7+mod(i+2,5))
          r(2,10+i) = x(2,i) + x(2,1+mod(i,5)) + x(2,7+mod(i+2,5))
          r(3,10+i) = x(3,i) + x(3,1+mod(i,5)) + x(3,7+mod(i+2,5))

          r(1,15+i) = x(1,6+i) + x(1,7+mod(i,5)) + x(1,1+mod(i+2,5))
          r(2,15+i) = x(2,6+i) + x(2,7+mod(i,5)) + x(2,1+mod(i+2,5))
          r(3,15+i) = x(3,6+i) + x(3,7+mod(i,5)) + x(3,1+mod(i+2,5))

 2040   continue
      endif

C --- Points on the thirty edges of the icosahedron ---
      if (npt .eq. 30) then
        do  3030  i = 1, 3
 3030   r2(i) = r(i,1)+r(i,2)
        d = dsqrt(ddot(3,r2,1,r2,1))
        call dpscop(r,x,36,1,1,1/d)
        do  3040  i = 1, 5
          r(1,i) = x(1,i) + x(1,6)
          r(2,i) = x(2,i) + x(2,6)
          r(3,i) = x(3,i) + x(3,6)

          r(1,5+i) = x(1,6+i) + x(1,12)
          r(2,5+i) = x(2,6+i) + x(2,12)
          r(3,5+i) = x(3,6+i) + x(3,12)

          r(1,10+i) = x(1,i) + x(1,7+mod(i+2,5))
          r(2,10+i) = x(2,i) + x(2,7+mod(i+2,5))
          r(3,10+i) = x(3,i) + x(3,7+mod(i+2,5))

          r(1,15+i) = x(1,6+i) + x(1,1+mod(i+2,5))
          r(2,15+i) = x(2,6+i) + x(2,1+mod(i+2,5))
          r(3,15+i) = x(3,6+i) + x(3,1+mod(i+2,5))

          r(1,20+i) = x(1,i) + x(1,1+mod(i,5))
          r(2,20+i) = x(2,i) + x(2,1+mod(i,5))
          r(3,20+i) = x(3,i) + x(3,1+mod(i,5))

          r(1,25+i) = x(1,6+i) + x(1,7+mod(i,5))
          r(2,25+i) = x(2,6+i) + x(2,7+mod(i,5))
          r(3,25+i) = x(3,6+i) + x(3,7+mod(i,5))

 3040   continue
      endif

C --- Buckeyball ---
C The 30 edges of the icosahedron lie at the midpoints of the arcs
C connecting the 60 buckeyball points.  Determine arc length wb from
C sol'n of sin^2(wi-wb)cos(2pi/5) + cos^2(wi-wb) = cos(2*wb)
C where wi is half the icosahedral angle.  Orient making a
C pentagon at the top, with symmetry r5z.  This is replicated 5 times
C by rotation by euler angles z:i*2pi/5,y:2*wi,z:pi for i=0..4 to
C generate the 30 points for the top half of the buckeyball.
      if (npt .eq. 60) then
        wi = datan(2/(1+sqrt(5d0)))
        wb = 0.2031689460497291d0
C ...   Five points of top pentagon
        r(1,1) = dsin(wi-wb)
        r(2,1) = 0
        r(3,1) = dcos(wi-wb)
        do  6010  i = 2, 5
          r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
          r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
          r(3,i) = r(3,1)
 6010   continue
C ...   Replicate these 5 points by rotations to other pentagons
        do  6020  i = 0, 4
          call eulerm(i*2*pi/5,2*wi,pi,rm)
          call dgemm('T','N',3,5,3,1d0,rm,3,r,3,0d0,r(1,6+5*i),3)
 6020   continue
C ...   Bottom 60 points are inversions of top 60
        do  6030  i = 1, 30
          r(1,61-i) = -r(1,i)
          r(2,61-i) = -r(2,i)
          r(3,61-i) = -r(3,i)
 6030   continue
      endif

C      print '(/'' platsl:'',i5,'' points generated:'')', npt
C      do  90  i = 1, npt
C   90 print 333, i, r(1,i), r(2,i), r(3,i)
C  333 format(i3,3f20.15)

C  999 continue
C      do  91  i = 1, npt
C      do  91  j = 1, npt
C   91 print 334, i,j,ddot(3,r(1,i),1,r(1,j),1)
C  334 format(2i3,f20.10)

      end
      subroutine eulerm(alpha,beta,gamma,r)
C- Generate the rotation matrix corresponding to Euler angles
C     implicit none
      double precision r(3,3),alpha,beta,gamma
      double precision ca,cb,cg,sa,sb,sg

      ca = dcos(alpha)
      sa = dsin(alpha)
      cb = dcos(beta)
      sb = dsin(beta)
      cg = dcos(gamma)
      sg = dsin(gamma)

      r(1,1) =  ca*cb*cg - sa*sg
      r(2,1) = -ca*cb*sg - sa*cg
      r(3,1) =  ca*sb
      r(1,2) =  sa*cb*cg + ca*sg
      r(2,2) = -sa*cb*sg + ca*cg
      r(3,2) =  sa*sb
      r(1,3) = -sb*cg
      r(2,3) =  sb*sg
      r(3,3) =  cb

C      print 335, ((r(i,j),j=1,3),i=1,3)
C  335 format((3f15.9))
      end

