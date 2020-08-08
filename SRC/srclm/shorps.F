      subroutine shorps(nbas,plat,mode,pin,pout)
C- Shift each basis vector by multiples of plat according to mode
C ----------------------------------------------------------------
Ci Inputs:  nbas,plat
Ci   nbas<0: sign used as a switch to make shorps return in
Ci         array pout the change in pin, in units of plat
Ci   pin:  position (basis) vectors
Ci   mode  vector of length 3 governing shifts along selected axes.
Ci         0 suppresses shifts along plat(j)
Ci         1 shifts to unit cell at origin (pos in 1st quadrant)
Ci         2 shifts to minimize length of pos
Co Outputs:
Co   pout  (may point to the same address space as pin).
Co   iat:  multiples of plat added
Cr Remarks
Cr   pos = f . plat, with f = integer + fraction for each plat.
Cr   Integer part according to mode.
Cu Updates
Cu   09 Jan 09 Do not change basis vectors if they exactly fall
Cu             on the cell boundary
Cu   10 Apr 02 Patch to handle mixed boundary conditions
C ----------------------------------------------------------------
C     implicit none
      integer nbas,mode(3)
      double precision plat(3,3),pin(3,nbas),pout(3,nbas)
C Local
      double precision qlat(3,3),xx,x0,x(3),a2,ap,p0(3),xmin(3),amin
      integer ib,i,m,j1,j2,j3,nbas1,j1max,j2max,j3max
      double precision tol
      parameter (tol = 1d-12)

C     call tcn('shorps')
      nbas1 = iabs(nbas)
C     call prmx('starting pos',pin,3,3,nbas1)

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,xx)

      do  10  ib = 1, nbas1

        call dpcopy(pin(1,ib),p0,1,3,1d0)

C   --- Reduce to unit cell centered at or near origin ---
        do  12  i = 1, 3
C   ... x0 is projection of pin along plat(i)
          x0 = pin(1,ib)*qlat(1,i)+pin(2,ib)*qlat(2,i)+pin(3,ib)*qlat(3,i)
          if (mode(i) .le. 0) then
            x(i) = x0
          else
c   ... leave basis vectors intact if |pin(i,ib)| = 0.5
            if (dabs(x0) .le. 0.5d0+tol) then
              xx = 0d0
            else
              xx = idnint(x0)
            endif
C     ... first octant for mode=1
            if (mode(i) .eq. 1 .and. x0-xx .lt. -tol) xx = xx-1
            x(i) = x0-xx
          endif
   12   continue
        do  14  m = 1, 3
          pout(m,ib) = x(1)*plat(m,1) + x(2)*plat(m,2) + x(3)*plat(m,3)
   14   continue

C   --- Try shortening by adding +/- lattice vectors ---
        j1max = 1
        if (mode(1) .le. 1) j1max = 0
        j2max = 1
        if (mode(2) .le. 1) j2max = 0
        j3max = 1
        if (mode(3) .le. 1) j3max = 0
   15   continue
        amin = 0
        do  16  j1 = -j1max, j1max
        do  16  j2 = -j2max, j2max
        do  16  j3 = -j3max, j3max

C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
          do  17  i = 1, 3
            x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
   17     continue
          a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          ap = pout(1,ib)*x(1) + pout(2,ib)*x(2) + pout(3,ib)*x(3)
          if (a2+2*ap .lt. amin) then
            xmin(1) = x(1)
            xmin(2) = x(2)
            xmin(3) = x(3)
            amin = a2+2*ap
          endif
   16   continue
c       if (amin .lt. 0) then
        if (amin .lt. -tol) then
          pout(1,ib) = pout(1,ib) + xmin(1)
          pout(2,ib) = pout(2,ib) + xmin(2)
          pout(3,ib) = pout(3,ib) + xmin(3)
C         In cases w/ mixed boundary conditions, (-1,0,1) may not be enough
C         Patched 10 Apr 02
          if (mode(1).eq.0 .or. mode(2).eq.0 .or. mode(3).eq.0) goto 15
        endif

C   --- pout <- pout - pin, units of plat ---
        if (nbas .lt. 0) then
          call dpadd(p0,pout(1,ib),1,3,-1d0)
          do  20  i = 1, 3
            xx = -p0(1)*qlat(1,i) - p0(2)*qlat(2,i) - p0(3)*qlat(3,i)
            if (dabs(xx-nint(xx)) .gt. 1d-10) call rx('bug in shorps')
            pout(i,ib) = xx
   20     continue
        endif
   10 continue

C     call prmx('plat',plat,3,3,3)
C     call prmx('qlat',qlat,3,3,3)
C     call prmx('shortened pos',pout,3,3,nbas1)
C     mc out.dat -t qlat -x gives pos as multiples of plat

C     call tcx('shorps')
      end

#if TEST
      subroutine fmain
C to see that the change in position vectors are multiples of
C the lattice vector, copy input,output pos to 'pos','posf'; invoke
C mc posf pos -- -t plat -t -i -x
      implicit none
      integer ix(3),i
      double precision pos(3,48),pos2(3,48),plat(3,3),qlat(9),xx

      call wkinit(10000)

      data plat /
     .0.5d0,          .5d0, 0d0,
     .0.0d0,          0.d0, 1d0,
     .2.570990255d0, -2.570990255d0, 0d0/
      data pos /
     .-0.697107d0,  1.197107d0,  0.250000d0,
     .-0.697107d0,  1.197107d0,  0.750000d0,
     .-0.770330d0,  0.770330d0,  0.000000d0,
     .-0.770330d0,  0.770330d0,  0.500000d0,
     .-0.343553d0,  0.843553d0,  0.250000d0,
     .-0.343553d0,  0.843553d0,  0.750000d0,
     .-0.416777d0,  0.416777d0,  0.000000d0,
     .-0.416777d0,  0.416777d0,  0.500000d0,
     .0.010000d0,  0.490000d0,  0.250000d0,
     .0.010000d0,  0.490000d0,  0.750000d0,
     .0.250000d0,  0.250000d0,  0.500000d0,
     .0.500000d0,  0.500000d0,  0.750000d0,
     .0.750000d0,  0.750000d0,  1.000000d0,
     .1.000000d0,  1.000000d0,  1.250000d0,
     .0.250000d0, -0.250000d0,  0.000000d0,
     .0.500000d0,  0.000000d0,  0.250000d0,
     .0.750000d0,  0.250000d0,  0.500000d0,
     .1.000000d0,  0.500000d0,  0.750000d0,
     .0.750000d0, -0.250000d0,  0.500000d0,
     .1.000000d0,  0.000000d0,  0.750000d0,
     .1.250000d0,  0.250000d0,  1.000000d0,
     .1.500000d0,  0.500000d0,  1.250000d0,
     .0.740000d0, -0.740000d0,  0.000000d0,
     .0.740000d0, -0.740000d0,  0.500000d0,
     .1.166777d0, -0.666777d0,  0.250000d0,
     .1.166777d0, -0.666777d0,  0.750000d0,
     .1.093553d0, -1.093553d0,  0.000000d0,
     .1.093553d0, -1.093553d0,  0.500000d0,
     .1.520330d0, -1.020330d0,  0.250000d0,
     .1.520330d0, -1.020330d0,  0.750000d0,
     .1.447107d0, -1.447107d0,  0.000000d0,
     .1.447107d0, -1.447107d0,  0.500000d0,
     .-1.050660d0,  1.550660d0,  0.250000d0,
     .-1.050660d0,  1.550660d0,  0.750000d0,
     .-1.123883d0,  1.123883d0,  0.000000d0,
     .-1.123883d0,  1.123883d0,  0.500000d0,
     .1.873883d0, -1.373883d0,  0.250000d0,
     .1.873883d0, -1.373883d0,  0.750000d0,
     .1.800660d0, -1.800660d0,  0.000000d0,
     .1.800660d0, -1.800660d0,  0.500000d0,
     .-1.404214d0,  1.904214d0,  0.250000d0,
     .-1.404214d0,  1.904214d0,  0.750000d0,
     .-1.477437d0,  1.477437d0,  0.000000d0,
     .-1.477437d0,  1.477437d0,  0.500000d0,
     .2.227437d0, -1.727437d0,  0.250000d0,
     .2.227437d0, -1.727437d0,  0.750000d0,
     .2.154214d0, -2.154214d0,  0.000000d0,
     .2.154214d0, -2.154214d0,  0.500000d0/


      call prmx('plat',plat,3,3,3)
      call prmx('starting pos',pos,3,3,48)
      ix(1) = 2
      ix(2) = 2
      ix(3) = 2
      call shorps(48,plat,ix,pos,pos2)
      call prmx('final pos',pos2,3,3,48)

      call mkqlat(plat,qlat,xx)
      do  10  i = 1, 48
        call shorbz(pos(1,i),pos2(1,i),plat,qlat)
   10 continue

      call prmx('from shorbz',pos2,3,3,48)
      end
#endif
#if TEST2
      subroutine fmain
C Check special case in which a bug fix, mixed boundary conditions
      implicit none
      integer ix(3),i
      double precision pos(3,1),pos2(3,1),plat(3,3),qlat(9),xx
      double precision dd1,dd2

      call wkinit(10000)

      data plat /-0.5d0,0.5d0,0.0d0,
     .0.0d0,0.0d0,1.0d0,
     .7.0d0,7.0d0,4.0d0/
      data ix /2,2,0/

      pos(1,1) = 2.5d0
      pos(2,1) = 3.0d0
      pos(3,1) = 0.0d0

      dd1 = dsqrt(pos(1,1)**2 + pos(2,1)**2 + pos(3,1)**2)
      call shorps(1,plat,ix,pos,pos2)
      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
      print *, dd1, dd2

      call mkqlat(plat,qlat,xx)
      do  10  i = 1, 1
        call shorbz(pos(1,i),pos2(1,i),plat,qlat)
   10 continue
      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
      print *, dd1, dd2

      end
#endif

