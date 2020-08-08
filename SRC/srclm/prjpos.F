      subroutine prjpos(mode,ix,plat1,plat2,nbas,pos,pos2)
C- Transform position vectors as multiples of a specified set of plat
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode     1's digit modifies plat1 and plat2:
Ci              0  unit vectors are used instead of plat1 and plat2
Ci                 Thus, the passed plat1 and plat2 are irrelevant.
Ci              1  Input pos are in units of the given plat1.
Ci              2  Output pos are in units of the given plat2.
Ci             >2  both 1+2
Ci           10's shifts the positions by multiples of plat to
Ci                achieve one of the following results:
Ci              0 no shifts
Ci              1  take fractional part of output pos to
Ci                 place in first octant of plat2.
Ci              2  shorten output pos2 by adding multiples
Ci                 of plat1
Ci          Add 4  For each component of pos, add a single
Ci                 integer to all sites that makes the smallest
Ci                 value of that component between 0 and 1.
Ci   ix       used to shorten pos (unused if 10s digit of mode ne 2)
Ci   plat1    lattice vectors in which input pos are represented
Ci            (depending on mode, unit vector my substitute for plat1)
Ci   plat2    lattice vectors in which output pos2 are represented
Ci            (depending on mode, unit vector my substitute for plat2)
Ci   pos,nbas position vectors as multiples of plat1, and number;
Co Outputs
Co   pos2     position vectors as multiples of plat2
Cr Remarks
Cr   If there are two triplets of lattice vectors P and S, then
Cr   the k_th component of any vector t may be expressed as
Cr   a linear combination of either P or S:
Cr      t_k = sum_m p_m P_km  = sum_m s_m S_km
Cr            sum_m p_m P+_mk = sum_m s_m (Q^-1)_km
Cr   where Q is the reciprocal lattice of S: (Q+ S) = 1
Cr   Coefficients p_m and s_m are related:
Cr      s_k = sum_m p_m (P+ Q)_mk
Cr   If the P's are cartesian coordinates, i.e.  P_mk = delta_mk
Cr      s_k = sum_m p_m Q_mk
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nbas,ix(3)
      double precision plat1(3,3),plat2(3,3),pos(3,nbas),pos2(3,nbas)
C Local
      logical lshrink,lcnst
      integer ib,i,nd,j1max,j2max,j3max,j1,j2,j3,mode1,n
      double precision xpos(3),qlat2(3,3),vol,xx(3),add,dlat(3,27),
     .dl2(3,27),a2,ap,tol,P(3,3),S(3,3)
      parameter (tol=-1d-10)

      mode1 = mod(mode/10,10)
      if (mod(mode,10) .eq. 1 .or. mod(mode,10) .gt. 2) then
        call dcopy(9,plat1,1,P,1)
      else
        do  2  j1 = 1, 3
          do  3  j2 = 1, 3
            P(j1,j2) = 0
    3     continue
          P(j1,j1) = 1
    2   continue
      endif
      if (mod(mode,10) .ge. 2) then
        call dcopy(9,plat2,1,S,1)
      else
        do  4  j1 = 1, 3
          do  5  j2 = 1, 3
            S(j1,j2) = 0
    5     continue
          S(j1,j1) = 1
    4   continue
      endif

C --- Create multiples of P ---
      if (mod(mode1,4) .eq. 2) then
        j1max = 1
        if (ix(1) .eq. 0) j1max = 0
        j2max = 1
        if (ix(2) .eq. 0) j2max = 0
        j3max = 1
        if (ix(3) .eq. 0) j3max = 0
   15   continue
        nd = 0
        do  16  j1 = -j1max, j1max
        do  16  j2 = -j2max, j2max
        do  16  j3 = -j3max, j3max
          nd = nd+1
          do  17  i = 1, 3
            dlat(i,nd) = P(i,1)*j1 + P(i,2)*j2 + P(i,3)*j3
   17     continue
   16   continue
      endif

C ... Basis vectors in units of S
      call dinv33(S,1,qlat2,vol)

C ... Multiples of P in units of S
      if (mod(mode1,4) .eq. 2) then
C       call prmx('dlat',dlat,3,3,nd)
        call dgemm('T','N',3,nd,3,1d0,qlat2,3,dlat,3,0d0,dl2,3)
C       call prmx('dlat (S)',dl2,3,3,nd)
      endif

      do  10  ib = 1, nbas
        call dgemm('N','N',3,1,3,1d0,P,3,pos(1,ib),3,0d0,xpos,3)
        call dgemm('N','N',1,3,3,1d0,xpos,1,qlat2,3,0d0,pos2(1,ib),1)
C   ... Shift to first octant
        if (mod(mode1,4) .eq. 1) then
          do  20  i = 1, 3
   22       continue
            if (pos2(i,ib) .lt. 0) then
              pos2(i,ib) = pos2(i,ib) + 1
              goto 22
            endif
            pos2(i,ib) = pos2(i,ib) - int(pos2(i,ib))
   20     continue
C   ... Shorten by adding multiples of P
        elseif (mod(mode1,4) .eq. 2) then
   28     continue
          do  27  n = 1, nd
            a2 = dl2(1,n)**2 + dl2(2,n)**2 + dl2(3,n)**2
            ap = pos2(1,ib)*dl2(1,n) +
     .      pos2(2,ib)*dl2(2,n) +
     .      pos2(3,ib)*dl2(3,n)
            if (a2 + 2*ap .lt. tol) then
              pos2(1,ib) = pos2(1,ib) + dl2(1,n)
              pos2(2,ib) = pos2(2,ib) + dl2(2,n)
              pos2(3,ib) = pos2(3,ib) + dl2(3,n)
              goto 28
            endif
   27     continue
        endif
        if (mode1 .ge. 4) then
          do  24  i = 1, 3
            if (ib .eq. 1) xx(i) = pos2(i,ib)
            xx(i) = min(xx(i),pos2(i,ib))
   24     continue
        endif
   10 continue

C ... Case 10s digit mode is >=4
      if (mode1 .ge. 4) then
        do  30  i = 1, 3
          if (xx(i) .ge. 0) then
            add = -int(xx(i))
          else
            add = int(-xx(i)) + 1
          endif
        do  30  ib = 1, nbas
   32     pos2(i,ib) = pos2(i,ib) + add
   30   continue
      endif

      if (mode1 .ge. 4) then
        do  50  ib = 1, nbas
   58     continue
          do  57  n = 1, nd
            a2 = dl2(1,n)**2 + dl2(2,n)**2 + dl2(3,n)**2
            ap = pos2(1,ib)*dl2(1,n) +
     .      pos2(2,ib)*dl2(2,n) +
     .      pos2(3,ib)*dl2(3,n)
            lshrink = a2 + 2*ap .lt. tol
            lcnst = pos2(1,ib) + dl2(1,n) .ge. 0d0
     .      .and. pos2(2,ib) + dl2(2,n) .ge. 0d0
     .      .and. pos2(3,ib) + dl2(3,n) .ge. 0d0
            if (lshrink .and. lcnst) then
              pos2(1,ib) = pos2(1,ib) + dl2(1,n)
              pos2(2,ib) = pos2(2,ib) + dl2(2,n)
              pos2(3,ib) = pos2(3,ib) + dl2(3,n)
              goto 58
            endif
   57     continue
   50   continue
      endif
      end

C$$$#if TEST
C$$$      subroutine fmain
C$$$      implicit none
C$$$      integer ix(3),i,j,m
C$$$      double precision pos(3,48),pos2(3,48),posp2(3,48),plat(3,3)
C$$$      double precision posp(3,48),plat2(3,3),pinv(3,3),p2inv(3,3),
C$$$     .  det,errmx
C$$$      real ran1
C$$$      logical ltmp,dcmp

C$$$      data plat /
C$$$     .  0.5d0,          .5d0, 0d0,
C$$$     .  0.0d0,          0.d0, 1d0,
C$$$     .  2.570990255d0, -2.570990255d0, 0d0/
C$$$      data plat2 / 1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/

C$$$      data pos /
C$$$     .  -0.697107d0,  1.197107d0,  0.250000d0,
C$$$     .  -0.697107d0,  1.197107d0,  0.750000d0,
C$$$     .  -0.770330d0,  0.770330d0,  0.000000d0,
C$$$     .  -0.770330d0,  0.770330d0,  0.500000d0,
C$$$     .  -0.343553d0,  0.843553d0,  0.250000d0,
C$$$     .  -0.343553d0,  0.843553d0,  0.750000d0,
C$$$     .  -0.416777d0,  0.416777d0,  0.000000d0,
C$$$     .  -0.416777d0,  0.416777d0,  0.500000d0,
C$$$     .   0.010000d0,  0.490000d0,  0.250000d0,
C$$$     .   0.010000d0,  0.490000d0,  0.750000d0,
C$$$     .   0.250000d0,  0.250000d0,  0.500000d0,
C$$$     .   0.500000d0,  0.500000d0,  0.750000d0,
C$$$     .   0.750000d0,  0.750000d0,  1.000000d0,
C$$$     .   1.000000d0,  1.000000d0,  1.250000d0,
C$$$     .   0.250000d0, -0.250000d0,  0.000000d0,
C$$$     .   0.500000d0,  0.000000d0,  0.250000d0,
C$$$     .   0.750000d0,  0.250000d0,  0.500000d0,
C$$$     .   1.000000d0,  0.500000d0,  0.750000d0,
C$$$     .   0.750000d0, -0.250000d0,  0.500000d0,
C$$$     .   1.000000d0,  0.000000d0,  0.750000d0,
C$$$     .   1.250000d0,  0.250000d0,  1.000000d0,
C$$$     .   1.500000d0,  0.500000d0,  1.250000d0,
C$$$     .   0.740000d0, -0.740000d0,  0.000000d0,
C$$$     .   0.740000d0, -0.740000d0,  0.500000d0,
C$$$     .   1.166777d0, -0.666777d0,  0.250000d0,
C$$$     .   1.166777d0, -0.666777d0,  0.750000d0,
C$$$     .   1.093553d0, -1.093553d0,  0.000000d0,
C$$$     .   1.093553d0, -1.093553d0,  0.500000d0,
C$$$     .   1.520330d0, -1.020330d0,  0.250000d0,
C$$$     .   1.520330d0, -1.020330d0,  0.750000d0,
C$$$     .   1.447107d0, -1.447107d0,  0.000000d0,
C$$$     .   1.447107d0, -1.447107d0,  0.500000d0,
C$$$     .  -1.050660d0,  1.550660d0,  0.250000d0,
C$$$     .  -1.050660d0,  1.550660d0,  0.750000d0,
C$$$     .  -1.123883d0,  1.123883d0,  0.000000d0,
C$$$     .  -1.123883d0,  1.123883d0,  0.500000d0,
C$$$     .   1.873883d0, -1.373883d0,  0.250000d0,
C$$$     .   1.873883d0, -1.373883d0,  0.750000d0,
C$$$     .   1.800660d0, -1.800660d0,  0.000000d0,
C$$$     .   1.800660d0, -1.800660d0,  0.500000d0,
C$$$     .  -1.404214d0,  1.904214d0,  0.250000d0,
C$$$     .  -1.404214d0,  1.904214d0,  0.750000d0,
C$$$     .  -1.477437d0,  1.477437d0,  0.000000d0,
C$$$     .  -1.477437d0,  1.477437d0,  0.500000d0,
C$$$     .   2.227437d0, -1.727437d0,  0.250000d0,
C$$$     .   2.227437d0, -1.727437d0,  0.750000d0,
C$$$     .   2.154214d0, -2.154214d0,  0.000000d0,
C$$$     .   2.154214d0, -2.154214d0,  0.500000d0/

C$$$      ix(1) = 2
C$$$      ix(2) = 2
C$$$      ix(3) = 2

C$$$      call ran1in(12)
C$$$      do  10  i = 1, 3
C$$$      do  10  j = 1, 3
C$$$      plat(i,j) = 2*ran1()-1
C$$$      plat2(i,j) = 2*ran1()-1
C$$$   10 continue

C$$$      call prmx('plat',plat,3,3,3)
C$$$      call dinv33(plat,0,pinv,det)
C$$$      call prmx('plat2',plat2,3,3,3)
C$$$      call dinv33(plat2,0,p2inv,det)

C$$$      call dgemm('N','N',3,48,3,1d0,pinv,3,pos,3,0d0,posp,3)
C$$$      call dgemm('N','N',3,48,3,1d0,p2inv,3,pos,3,0d0,posp2,3)
C$$$      call prmx('starting pos, cartesian coord',pos,3,3,48)
C$$$      call prmx('starting pos, units of plat',posp,3,3,48)
C$$$      call prmx('starting pos, units of plat2',posp2,3,3,48)

C$$$      print *, 'test 1s digit = 0: test pos -> pos'
C$$$      call prjpos(0,ix,plat,plat2,48,pos,pos2)
C$$$      if (.not. dcmp(pos2,pos,3*38,1d-6,m,errmx)) then
C$$$        call prmx('test failed: write pos2',pos2,3,3,48)
C$$$      else
C$$$        print *, 'test passed'
C$$$      endif

C$$$      print *, 'test 1s digit = 1: test pos(plat1) -> pos'
C$$$      call prjpos(1,ix,plat,plat2,48,posp,pos2)
C$$$      if (.not. dcmp(pos2,pos,3*38,1d-6,m,errmx)) then
C$$$        call prmx('test failed: write pos2',pos2,3,3,48)
C$$$      else
C$$$        print *, 'test passed'
C$$$      endif

C$$$      print *, 'test 1s digit = 2: test pos -> pos(plat2)'
C$$$      call prjpos(2,ix,plat,plat2,48,pos,pos2)
C$$$      if (.not. dcmp(pos2,posp2,3*38,1d-6,m,errmx)) then
C$$$        call prmx('test failed: write pos2',pos2,3,3,48)
C$$$      else
C$$$        print *, 'test passed'
C$$$      endif

C$$$      print *, 'test 1s digit = 3: test pos(plat1) -> pos(plat2)'
C$$$      call prjpos(3,ix,plat,plat2,48,posp,pos2)
C$$$      if (.not. dcmp(pos2,posp2,3*38,1d-6,m,errmx)) then
C$$$        call prmx('test failed: write pos2',pos2,3,3,48)
C$$$      else
C$$$        print *, 'test passed'
C$$$      endif

C$$$      print *, 'test 10s digit = 1: test pos(plat1) -> pos(plat2)'
C$$$      call prjpos(13,ix,plat,plat2,48,posp,pos2)
C$$$      call prmx('calc pos mode 13, units of plat2',pos2,3,3,48)
C$$$C ... compare with: mc out.dat posplat2 -- -t
C$$$      print *, 'test 10s digit = 2: test pos(plat1) -> pos(plat2)'
C$$$      call prjpos(23,ix,plat,plat2,48,posp,pos2)
C$$$      call prmx('calc pos mode 23, units of plat2',pos2,3,3,48)
C$$$C ... compare with: mc plat -i plat2 out.dat -x pos -- -x -t
C$$$      print *, 'test 10s digit = 6: test pos(plat1) -> pos(plat2)'
C$$$      call prjpos(63,ix,plat,plat2,48,posp,pos2)
C$$$      call prmx('calc pos mode 63, units of plat2',pos2,3,3,48)
C$$$C ... check: for above, rename out.dat as out and then:
C$$$C     compare with: mc out.dat out -- -t

C$$$      call prjpos(0,ix,plat,plat2,48,pos,pos2)
C$$$      call prmx('calc pos mode 0, units of plat2',pos2,3,3,48)

C$$$      end
C$$$#endif

