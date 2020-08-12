      subroutine xxxdif(emin,emax,npts,ndos,mode,dos)
C- Differentiate Number of States to make DOS
C     implicit none
      integer npts,ndos,mode
      double precision dos(npts,ndos+1),emin,emax
      integer ip,idos,k
      double precision bin
      bin = 2*(emax-emin)/dble(npts-1)
      k = ndos+1
      do  10  idos = 1, ndos
        do  20  ip = 2, npts-1
          dos(ip,k) = (dos(ip+1,idos) - dos(ip-1,idos)) / bin
   20   continue
        dos(1,k) = dos(2,k)
        dos(npts,k) = dos(npts-1,k)
        call dcopy(npts,dos(1,k),1,dos(1,idos),1)
   10 continue
      end

