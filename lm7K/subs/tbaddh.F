      subroutine tbaddh(nl,nbas,nsp,it,nitmax,npr,delta,hrs)
C- 
C ----------------------------------------------------------------------
Ci Inputs:
Ci   
Co Outputs:
Co   
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nbas,nsp,it,nitmax,npr(0:1,nbas)
      double precision delta(nl**2,nl**2,nbas,nitmax+1,2),
     .                 hrs(nl**2,nl**2,1)

C Local Variables
      integer ib,j,ilm,ilmp,iprint

      if (iprint() .gt. 50) print 100
      do  3  ib = 1, nbas
        j = npr(1,ib) + 1
        do  2  ilm = 1, nl**2
          do  1  ilmp = 1, nl**2
            if (ilm .eq. ilmp .and. iprint() .gt. 50) write (*,200)
     .        ib,ilm,hrs(ilm,ilmp,j),delta(ilm,ilmp,ib,it,1),
     .        delta(ilm,ilmp,ib,it+1,1),hrs(ilm,ilmp,j)-
     .        delta(ilm,ilmp,ib,it,1)+delta(ilm,ilmp,ib,it+1,1)
            hrs(ilm,ilmp,j) = hrs(ilm,ilmp,j)
     .                      - delta(ilm,ilmp,ib,it,1)
     .                      + delta(ilm,ilmp,ib,it+1,1)
    1     continue
    2   continue
    3 continue
  100 format ('TBADDH: site  ilm  H_in  Delta_last  Delta_now  H_out')
  200 format (7x,2i4,4f10.4)
      end

