      subroutine npr2tb(mode,nbas,npr,ntab)
C- Convert npr array to ntab array
C  mode=0: npr -> ntab
C  mode=1: ntab -> npr
C     implicit none
      integer nbas,mode,npr(2,nbas),ntab(nbas+1)
      integer ib,npi,nttab

      if (mode .eq. 0) then
        nttab = 0
        do  10  ib = 1, nbas
          npi = npr(1,ib)
          ntab(ib) = nttab
          nttab = nttab+npi
   10   continue
        ntab(nbas+1) = nttab

C        print *, 'convert to ntab'
C        print 333, (ntab(ib+1)-ntab(ib), ib=1,nbas)
C        print 333, (ntab(ib), ib=1,nbas)
C  333   format(12i3)
        return
      else
        do  20  ib = nbas, 1, -1
          npi   = ntab(ib+1)-ntab(ib)
          nttab = ntab(ib)
          npr(1,ib) = npi
          npr(2,ib) = nttab
   20   continue
C        print *, 'convert to npr'
C        print 333, (npr(1,ib), ib=1,nbas)
C        print 333, (npr(2,ib), ib=1,nbas)
      endif
      end
C      implicit none
C
C      integer nbas
C      parameter (nbas=12)
C      integer ntab(2*nbas+1),npr(2,1),i
C      equivalence (ntab,npr)
C      data ntab / 0,3,7,7,9,12,20,33,35,40,60,61,71/
C      
C      print 333, (ntab(i+1)-ntab(i), i=1,nbas)
C      print 333, (ntab(i), i=1,nbas)
C  333 format(12i3)
C
C      call npr2tb(1,nbas,ntab,npr)
C
C      print 333, (npr(1,i), i=1,nbas)
C      print 333, (npr(2,i), i=1,nbas)
C
C      call npr2tb(0,nbas,ntab,npr)
C      print 333, (ntab(i+1)-ntab(i), i=1,nbas)
C      print 333, (ntab(i), i=1,nbas)
C
C      end

