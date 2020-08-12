      subroutine nghbor(nbas,plat,bas,rmax1,rmax,iat1,mxnbr,
     .npr,iax,wk)
C- Create a table of all neighbors within spec'd range of a spec'd atom
C ----------------------------------------------------------------
Ci Inputs:
Ci   nbas,plat,bas
Ci   rmax:  restrict neighbors to within r < rmax.
Ci   rmax1: see npr(1) below
Ci   iat1: index to spec'd atom in bas
Ci   mxnbr:maximum number of neighbors for which memory is allocated
Co Outputs:
Co   iax
Co   npr(0) returns number of neighbors within rmax
Co   npr(1) returns number of neighbors within rmax1 < rmax
Co   wk: a work array of dimension mxnbr
Cr Remarks:
Cr   Creates a neighour list for a specified atom, generating iax
Cr   Table is sorted according to increasing length from central atom.
Cr   Sort of iax is according to following priority:
Cr   (1) distance, (2) iat2, (3), i, (4) j, (5) k.
Cr   may be used.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer iat1,mxnbr,nbas,niax
      parameter (niax=10)
      integer iax(niax,1),npr(0:1)
      double precision rmax1,rmax
      double precision plat(3,3),bas(3,nbas),wk(mxnbr)
C Local parameters
      integer i1,i2,i3,iprint,ix,i,j,k,iat2
      double precision dsqr,drr2
      double precision dr(3),rmxsqr
      external iprint,drr2

      rmxsqr = rmax**2
      call latlim(plat,rmax,i1,i2,i3)
C
C Sweep these lattice translations to find all neighbors within rmax
C Note: range i +/-1 is not sufficient; test case of LaCuO
      npr(0) = 0
      npr(1) = 0
      do  100  iat2 = 1, nbas
        do  110  i = -i1-2, i1+2
        do  110  j = -i2-2, i2+2
        do  110  k = -i3-2, i3+2
          dsqr = drr2(plat,bas(1,iat1),bas(1,iat2),i,j,k,dr)
          if (dsqr .le. rmxsqr) then
            npr(0) = npr(0) + 1
            if (npr(0) .gt. mxnbr)
     .          call fexit(-1,111,
     .          ' Exit -1 NGHBOR: too many pairs: npair=%i',npr)
            iax(1,npr(0)) = 10000*dsqrt(dsqr)
            iax(2,npr(0)) = iat2
            iax(3,npr(0)) = i
            iax(4,npr(0)) = j
            iax(5,npr(0)) = k
            if (dsqr .le. rmax1**2) npr(1) = npr(1) + 1
          endif
  110   continue
  100 continue

      call ivshel(10,npr,iax,wk,.false.)

C for debugging ...
c      npr(0) = 2
c      npr(1) = 2

      if (iprint() .ge. 50) print 300, iat1
  300 format(/' NGHBOR: neighbors around atom ',i3,
     ./,12x,'D',4x,'ATOM',3x,'--x PLAT---',14x,'----TAU----')
      do  200  i = 1, npr(0)
        iax(1,i) = iat1
        dsqr = dsqrt(drr2(plat,bas(1,iat1),bas(1,iax(2,i)),
     .  iax(3,i),iax(4,i),iax(5,i),dr))
        wk(i) = dsqr
        if (iprint() .ge. 50) then
          print 301, i, dsqr, (iax(j,i), j=2,5), (dr(ix), ix=1,3)
  301     format(i4,f11.4,4i5,3f11.4)
        endif
  200 continue

      if (iprint() .ge. 50) print *
      end

