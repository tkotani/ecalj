      subroutine hoffs(nbas,lmx,ips,ioffb)
C- Make offset tables ioffb for basis
C     implicit none
      integer nbas,lmx(2),ioffb(nbas),ips(1),ib,is,ioff

C --- Offsets table for hamiltonian ---
      ioffb(1) = 0
      ioff = 0
      do  10  ib = 1, nbas
        is = ips(ib)
        ioff = ioff + (lmx(is)+1)**2
        ioffb(ib+1) = ioff
   10 continue

      end

