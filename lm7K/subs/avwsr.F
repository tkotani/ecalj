      double precision function avwsr(plat,alat,vol,nbas)
C- Calculate the average ws radius
C     implicit none
      integer nbas
      double precision plat(3,3),alat,vol
      vol = alat**3*dabs(
     .      plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .      plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .      plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))
      avwsr = (3/(16*datan(1.d0)*nbas)*vol)**(1.d0/3.d0)
      end

