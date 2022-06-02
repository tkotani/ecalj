double precision function pythag(a,b)
  !     implicit none
  double precision :: a,b
  !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
  double precision :: p,r,s,t,u
  p = dmax1(dabs(a),dabs(b))
  if (p == 0.0d0) go to 20
  r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
  t = 4.0d0 + r
  if (t == 4.0d0) go to 20
  s = r/t
  u = 1.0d0 + 2.0d0*s
  p = u*p
  r = (s/u)**2 * r
  go to 10
20 pythag = p
  return
END function pythag

