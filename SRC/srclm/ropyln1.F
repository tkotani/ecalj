      subroutine ropqln(m,l,n,r2,z,cx,q,kk)
C- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
c  Returns kk, which points to the current component of q.
C  These subroutines are utility routines called by ropyln.f.
C  Kept separate from ropyln because some optimizing compilers have bugs
C  (e.g. intel ifort version 11).
C  These routines are the time-critical steps.
C     implicit none
      integer mm,n,i,l,m,kk,k2,k1
      double precision q(n,2),r2(n),z(n),cx(3)
      double precision a,b,xx,yy

C --- Case l=m ---
      if (l .eq. m) then
        a = 1d0
        do  1  mm = 0, m-1
          a = a*(2*mm+1)
    1   continue
        kk = 1
        a = a*cx(1)
        do  2  i = 1, n
          q(i,kk) = a
    2   continue
        return
      endif

C --- Case l=m+1 ---
      if (l .eq. m+1) then
        b = 1d0
        do  3  mm = 0, m
          b = b*(2*mm+1)
    3   continue
        b = b*cx(1)
        kk = 2
        do  4  i = 1, n
          q(i,kk) = b*z(i)
    4   continue
        return
      endif

C --- Case l=m+2 and higher by recursion ---
      if (l .ge. m+2) then
        k2 = kk
        k1 = kk+1
        if (k1 .eq. 3) k1 = 1
        xx = -(l+m-1d0)/(l-m)*cx(1)/cx(3)
        yy = (2*l-1d0)/(l-m)*cx(1)/cx(2)
        do  6  i = 1, n
          q(i,k1) = xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
    6   continue
        kk = k1
        return
      endif
      end

      subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
C     implicit none
      integer lav,n,nd,l,i,m,kk
      double precision q(n,2),cm(n),sm(n),yl(nd,1)
      lav = l*(l+1)+1
      do  1  i = 1, n
        yl(i,lav+m) = cm(i)*q(i,kk)
    1 continue
      if (m .eq. 0) return
      do  2  i = 1, n
        yl(i,lav-m) = sm(i)*q(i,kk)
    2 continue
      end

