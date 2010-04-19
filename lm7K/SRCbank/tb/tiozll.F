      subroutine tiozll(lgamma,nfilet,ldim,nev,z)
C- Put e-vecs on disc
C     implicit none
      integer nfilet, ldim,nev
      logical lgamma
      double precision z(ldim,ldim,2)
      integer kmax,i,j,k
      kmax = 2
      if (lgamma) kmax = 1
      if (nfilet .lt. 0) then
         write (-nfilet) nev
         write (-nfilet) (((z(i,j,k),i=1,ldim),j=1,nev),k=1,kmax)
      else
         call dcopy(kmax*ldim*ldim,0d0,0,z,1)
         read   (nfilet) nev
         read   (nfilet) (((z(i,j,k),i=1,ldim),j=1,nev),k=1,kmax)
      endif
      end
