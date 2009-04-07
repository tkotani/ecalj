      double precision function cksumf(a,n)
C- Makes a checksumf from an array
C     implicit none
      integer n,i
      double precision a(n)
      double precision pi,fac
      pi = 4*datan(1d0)

      fac = pi

      cksumf = 0d0
      do  10  i = 1, n
        fac = fac+pi
        fac = fac - int(fac)
        cksumf = cksumf + fac*a(i)
   10 continue
      cksumf = cksumf - int(cksumf)

      end

