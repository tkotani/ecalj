subroutine prjpos(mode,ix,plat1,plat2,nbas,pos,pos2)
  !- Transform position vectors as multiples of a specified set of plat
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode     1's digit modifies plat1 and plat2:
  !i              0  unit vectors are used instead of plat1 and plat2
  !i                 Thus, the passed plat1 and plat2 are irrelevant.
  !i              1  Input pos are in units of the given plat1.
  !i              2  Output pos are in units of the given plat2.
  !i             >2  both 1+2
  !i           10's shifts the positions by multiples of plat to
  !i                achieve one of the following results:
  !i              0 no shifts
  !i              1  take fractional part of output pos to
  !i                 place in first octant of plat2.
  !i              2  shorten output pos2 by adding multiples
  !i                 of plat1
  !i          Add 4  For each component of pos, add a single
  !i                 integer to all sites that makes the smallest
  !i                 value of that component between 0 and 1.
  !i   ix       used to shorten pos (unused if 10s digit of mode ne 2)
  !i   plat1    lattice vectors in which input pos are represented
  !i            (depending on mode, unit vector my substitute for plat1)
  !i   plat2    lattice vectors in which output pos2 are represented
  !i            (depending on mode, unit vector my substitute for plat2)
  !i   pos,nbas position vectors as multiples of plat1, and number;
  !o Outputs
  !o   pos2     position vectors as multiples of plat2
  !r Remarks
  !r   If there are two triplets of lattice vectors P and S, then
  !r   the k_th component of any vector t may be expressed as
  !r   a linear combination of either P or S:
  !r      t_k = sum_m p_m P_km  = sum_m s_m S_km
  !r            sum_m p_m P+_mk = sum_m s_m (Q^-1)_km
  !r   where Q is the reciprocal lattice of S: (Q+ S) = 1
  !r   Coefficients p_m and s_m are related:
  !r      s_k = sum_m p_m (P+ Q)_mk
  !r   If the P's are cartesian coordinates, i.e.  P_mk = delta_mk
  !r      s_k = sum_m p_m Q_mk
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: mode,nbas,ix(3)
  double precision :: plat1(3,3),plat2(3,3),pos(3,nbas),pos2(3,nbas)
  ! Local
  logical :: lshrink,lcnst
  integer :: ib,i,nd,j1max,j2max,j3max,j1,j2,j3,mode1,n
  double precision :: xpos(3),qlat2(3,3),vol,xx(3),add,dlat(3,27), &
       dl2(3,27),a2,ap,tol,P(3,3),S(3,3)
  parameter (tol=-1d-10)

  mode1 = mod(mode/10,10)
  if (mod(mode,10) == 1 .OR. mod(mode,10) > 2) then
     call dcopy(9,plat1,1,P,1)
  else
     do  2  j1 = 1, 3
        do  3  j2 = 1, 3
           P(j1,j2) = 0
3       enddo
        P(j1,j1) = 1
2    enddo
  endif
  if (mod(mode,10) >= 2) then
     call dcopy(9,plat2,1,S,1)
  else
     do  4  j1 = 1, 3
        do  5  j2 = 1, 3
           S(j1,j2) = 0
5       enddo
        S(j1,j1) = 1
4    enddo
  endif

  ! --- Create multiples of P ---
  if (mod(mode1,4) == 2) then
     j1max = 1
     if (ix(1) == 0) j1max = 0
     j2max = 1
     if (ix(2) == 0) j2max = 0
     j3max = 1
     if (ix(3) == 0) j3max = 0
15   continue
     nd = 0
     do   j1 = -j1max, j1max
        do   j2 = -j2max, j2max
           do   j3 = -j3max, j3max
              nd = nd+1
              do   i = 1, 3
                 dlat(i,nd) = P(i,1)*j1 + P(i,2)*j2 + P(i,3)*j3
              enddo
           enddo
        enddo
     enddo

  endif

  ! ... Basis vectors in units of S
  call dinv33(S,1,qlat2,vol)

  ! ... Multiples of P in units of S
  if (mod(mode1,4) == 2) then
     !       call prmx('dlat',dlat,3,3,nd)
     call dgemm('T','N',3,nd,3,1d0,qlat2,3,dlat,3,0d0,dl2,3)
     !       call prmx('dlat (S)',dl2,3,3,nd)
  endif

  do  10  ib = 1, nbas
     call dgemm('N','N',3,1,3,1d0,P,3,pos(1,ib),3,0d0,xpos,3)
     call dgemm('N','N',1,3,3,1d0,xpos,1,qlat2,3,0d0,pos2(1,ib),1)
     !   ... Shift to first octant
     if (mod(mode1,4) == 1) then
        do  20  i = 1, 3
22         continue
           if (pos2(i,ib) < 0) then
              pos2(i,ib) = pos2(i,ib) + 1
              goto 22
           endif
           pos2(i,ib) = pos2(i,ib) - int(pos2(i,ib))
20      enddo
        !   ... Shorten by adding multiples of P
     elseif (mod(mode1,4) == 2) then
28      continue
        do  27  n = 1, nd
           a2 = dl2(1,n)**2 + dl2(2,n)**2 + dl2(3,n)**2
           ap = pos2(1,ib)*dl2(1,n) + &
                pos2(2,ib)*dl2(2,n) + &
                pos2(3,ib)*dl2(3,n)
           if (a2 + 2*ap < tol) then
              pos2(1,ib) = pos2(1,ib) + dl2(1,n)
              pos2(2,ib) = pos2(2,ib) + dl2(2,n)
              pos2(3,ib) = pos2(3,ib) + dl2(3,n)
              goto 28
           endif
27      enddo
     endif
     if (mode1 >= 4) then
        do  24  i = 1, 3
           if (ib == 1) xx(i) = pos2(i,ib)
           xx(i) = min(xx(i),pos2(i,ib))
24      enddo
     endif
10 enddo

  ! ... Case 10s digit mode is >=4
  if (mode1 >= 4) then
     do  301  i = 1, 3
        if (xx(i) >= 0) then
           add = -int(xx(i))
        else
           add = int(-xx(i)) + 1
        endif
        do  30  ib = 1, nbas
32         pos2(i,ib) = pos2(i,ib) + add
30      enddo
301  enddo
  endif

  if (mode1 >= 4) then
     do  50  ib = 1, nbas
58      continue
        do  57  n = 1, nd
           a2 = dl2(1,n)**2 + dl2(2,n)**2 + dl2(3,n)**2
           ap = pos2(1,ib)*dl2(1,n) + &
                pos2(2,ib)*dl2(2,n) + &
                pos2(3,ib)*dl2(3,n)
           lshrink = a2 + 2*ap .lt. tol
           lcnst = pos2(1,ib) + dl2(1,n) .ge. 0d0 &
                .and. pos2(2,ib) + dl2(2,n) .ge. 0d0 &
                .and. pos2(3,ib) + dl2(3,n) .ge. 0d0
           if (lshrink .AND. lcnst) then
              pos2(1,ib) = pos2(1,ib) + dl2(1,n)
              pos2(2,ib) = pos2(2,ib) + dl2(2,n)
              pos2(3,ib) = pos2(3,ib) + dl2(3,n)
              goto 58
           endif
57      enddo
50   enddo
  endif
end subroutine prjpos

!$$$#if TEST
!$$$      subroutine fmain
!$$$      implicit none
!$$$      integer ix(3),i,j,m
!$$$      double precision pos(3,48),pos2(3,48),posp2(3,48),plat(3,3)
!$$$      double precision posp(3,48),plat2(3,3),pinv(3,3),p2inv(3,3),
!$$$     .  det,errmx
!$$$      real ran1
!$$$      logical ltmp,dcmp

!$$$      data plat /
!$$$     .  0.5d0,          .5d0, 0d0,
!$$$     .  0.0d0,          0.d0, 1d0,
!$$$     .  2.570990255d0, -2.570990255d0, 0d0/
!$$$      data plat2 / 1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/

!$$$      data pos /
!$$$     .  -0.697107d0,  1.197107d0,  0.250000d0,
!$$$     .  -0.697107d0,  1.197107d0,  0.750000d0,
!$$$     .  -0.770330d0,  0.770330d0,  0.000000d0,
!$$$     .  -0.770330d0,  0.770330d0,  0.500000d0,
!$$$     .  -0.343553d0,  0.843553d0,  0.250000d0,
!$$$     .  -0.343553d0,  0.843553d0,  0.750000d0,
!$$$     .  -0.416777d0,  0.416777d0,  0.000000d0,
!$$$     .  -0.416777d0,  0.416777d0,  0.500000d0,
!$$$     .   0.010000d0,  0.490000d0,  0.250000d0,
!$$$     .   0.010000d0,  0.490000d0,  0.750000d0,
!$$$     .   0.250000d0,  0.250000d0,  0.500000d0,
!$$$     .   0.500000d0,  0.500000d0,  0.750000d0,
!$$$     .   0.750000d0,  0.750000d0,  1.000000d0,
!$$$     .   1.000000d0,  1.000000d0,  1.250000d0,
!$$$     .   0.250000d0, -0.250000d0,  0.000000d0,
!$$$     .   0.500000d0,  0.000000d0,  0.250000d0,
!$$$     .   0.750000d0,  0.250000d0,  0.500000d0,
!$$$     .   1.000000d0,  0.500000d0,  0.750000d0,
!$$$     .   0.750000d0, -0.250000d0,  0.500000d0,
!$$$     .   1.000000d0,  0.000000d0,  0.750000d0,
!$$$     .   1.250000d0,  0.250000d0,  1.000000d0,
!$$$     .   1.500000d0,  0.500000d0,  1.250000d0,
!$$$     .   0.740000d0, -0.740000d0,  0.000000d0,
!$$$     .   0.740000d0, -0.740000d0,  0.500000d0,
!$$$     .   1.166777d0, -0.666777d0,  0.250000d0,
!$$$     .   1.166777d0, -0.666777d0,  0.750000d0,
!$$$     .   1.093553d0, -1.093553d0,  0.000000d0,
!$$$     .   1.093553d0, -1.093553d0,  0.500000d0,
!$$$     .   1.520330d0, -1.020330d0,  0.250000d0,
!$$$     .   1.520330d0, -1.020330d0,  0.750000d0,
!$$$     .   1.447107d0, -1.447107d0,  0.000000d0,
!$$$     .   1.447107d0, -1.447107d0,  0.500000d0,
!$$$     .  -1.050660d0,  1.550660d0,  0.250000d0,
!$$$     .  -1.050660d0,  1.550660d0,  0.750000d0,
!$$$     .  -1.123883d0,  1.123883d0,  0.000000d0,
!$$$     .  -1.123883d0,  1.123883d0,  0.500000d0,
!$$$     .   1.873883d0, -1.373883d0,  0.250000d0,
!$$$     .   1.873883d0, -1.373883d0,  0.750000d0,
!$$$     .   1.800660d0, -1.800660d0,  0.000000d0,
!$$$     .   1.800660d0, -1.800660d0,  0.500000d0,
!$$$     .  -1.404214d0,  1.904214d0,  0.250000d0,
!$$$     .  -1.404214d0,  1.904214d0,  0.750000d0,
!$$$     .  -1.477437d0,  1.477437d0,  0.000000d0,
!$$$     .  -1.477437d0,  1.477437d0,  0.500000d0,
!$$$     .   2.227437d0, -1.727437d0,  0.250000d0,
!$$$     .   2.227437d0, -1.727437d0,  0.750000d0,
!$$$     .   2.154214d0, -2.154214d0,  0.000000d0,
!$$$     .   2.154214d0, -2.154214d0,  0.500000d0/

!$$$      ix(1) = 2
!$$$      ix(2) = 2
!$$$      ix(3) = 2

!$$$      call ran1in(12)
!$$$      do  10  i = 1, 3
!$$$      do  10  j = 1, 3
!$$$      plat(i,j) = 2*ran1()-1
!$$$      plat2(i,j) = 2*ran1()-1
!$$$   10 continue

!$$$      call prmx('plat',plat,3,3,3)
!$$$      call dinv33(plat,0,pinv,det)
!$$$      call prmx('plat2',plat2,3,3,3)
!$$$      call dinv33(plat2,0,p2inv,det)

!$$$      call dgemm('N','N',3,48,3,1d0,pinv,3,pos,3,0d0,posp,3)
!$$$      call dgemm('N','N',3,48,3,1d0,p2inv,3,pos,3,0d0,posp2,3)
!$$$      call prmx('starting pos, cartesian coord',pos,3,3,48)
!$$$      call prmx('starting pos, units of plat',posp,3,3,48)
!$$$      call prmx('starting pos, units of plat2',posp2,3,3,48)

!$$$      print *, 'test 1s digit = 0: test pos -> pos'
!$$$      call prjpos(0,ix,plat,plat2,48,pos,pos2)
!$$$      if (.not. dcmp(pos2,pos,3*38,1d-6,m,errmx)) then
!$$$        call prmx('test failed: write pos2',pos2,3,3,48)
!$$$      else
!$$$        print *, 'test passed'
!$$$      endif

!$$$      print *, 'test 1s digit = 1: test pos(plat1) -> pos'
!$$$      call prjpos(1,ix,plat,plat2,48,posp,pos2)
!$$$      if (.not. dcmp(pos2,pos,3*38,1d-6,m,errmx)) then
!$$$        call prmx('test failed: write pos2',pos2,3,3,48)
!$$$      else
!$$$        print *, 'test passed'
!$$$      endif

!$$$      print *, 'test 1s digit = 2: test pos -> pos(plat2)'
!$$$      call prjpos(2,ix,plat,plat2,48,pos,pos2)
!$$$      if (.not. dcmp(pos2,posp2,3*38,1d-6,m,errmx)) then
!$$$        call prmx('test failed: write pos2',pos2,3,3,48)
!$$$      else
!$$$        print *, 'test passed'
!$$$      endif

!$$$      print *, 'test 1s digit = 3: test pos(plat1) -> pos(plat2)'
!$$$      call prjpos(3,ix,plat,plat2,48,posp,pos2)
!$$$      if (.not. dcmp(pos2,posp2,3*38,1d-6,m,errmx)) then
!$$$        call prmx('test failed: write pos2',pos2,3,3,48)
!$$$      else
!$$$        print *, 'test passed'
!$$$      endif

!$$$      print *, 'test 10s digit = 1: test pos(plat1) -> pos(plat2)'
!$$$      call prjpos(13,ix,plat,plat2,48,posp,pos2)
!$$$      call prmx('calc pos mode 13, units of plat2',pos2,3,3,48)
!$$$C ... compare with: mc out.dat posplat2 -- -t
!$$$      print *, 'test 10s digit = 2: test pos(plat1) -> pos(plat2)'
!$$$      call prjpos(23,ix,plat,plat2,48,posp,pos2)
!$$$      call prmx('calc pos mode 23, units of plat2',pos2,3,3,48)
!$$$C ... compare with: mc plat -i plat2 out.dat -x pos -- -x -t
!$$$      print *, 'test 10s digit = 6: test pos(plat1) -> pos(plat2)'
!$$$      call prjpos(63,ix,plat,plat2,48,posp,pos2)
!$$$      call prmx('calc pos mode 63, units of plat2',pos2,3,3,48)
!$$$C ... check: for above, rename out.dat as out and then:
!$$$C     compare with: mc out.dat out -- -t

!$$$      call prjpos(0,ix,plat,plat2,48,pos,pos2)
!$$$      call prmx('calc pos mode 0, units of plat2',pos2,3,3,48)

!$$$      end
!$$$#endif

