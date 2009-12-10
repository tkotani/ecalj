      subroutine spbndw(nbas,ntab,iax,offH,kl,ku)
C- Return width of band subdiagonal and superdiagonal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Co Outputs
Co   kl    :size of band subdiagonal
Co   ku    :size of band superdiagonal
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
      integer nbas,niax,kl,ku,nkap0,n0H
      parameter (nkap0=3,n0H=5,niax=10)
      integer ntab(nbas+1),iax(niax,1),offh(n0H,nkap0,nbas+1)
C ... Local parameters
      integer ib,jb,ibmn,ibmx,ic,kli,kui

      kl = 0
      do  10  ib = 1, nbas
        ibmn = nbas
        ibmx = 0
        do  12  ic = ntab(ib)+1, ntab(ib+1)
          jb = iax(2,ic)
          if (jb .le. nbas) then
            ibmn = min(ibmn,jb)
            ibmx = max(ibmx,jb)
          endif
          kli = offH(1,1,ibmx+1) - (offH(1,1,ib)+1)
          kui = offH(1,1,ib+1) -(offH(1,1,ibmn)+1)
          kl = max(kl,kli)
          ku = max(ku,kui)
   12   continue
C        print 333, ib,ibmx,offH(1,ibmx)+1,offH(1,ibmx+1),kli,
C     .                ibmn,offH(1,ibmn)+1,offH(1,ibmn+1),kui
C  333   format(i4,2x,4i4,2x,4i4)
   10 continue

      end


