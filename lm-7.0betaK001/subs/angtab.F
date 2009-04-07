      subroutine angtab(nbas,bas,alat,rmax,qss,qlat,dlat,nkd,
     .  ipc,neul,eula)
C- Print angles table
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   rmax  :augmentation radius, in a.u.,
Ci   qss   :Parameters defining spin-spiral
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   dlat  :direct lattice vectors
Ci   nkd   :number of dlat
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   eula  :Euler angles for noncollinear spins
Co Outputs
Cl Local variables
Cl         :
C  BUG: SS addition to Euler angle assumes along z, and no beta!
Cr Remarks
Cr
Cu Updates
Cu   26 Jan 03 got working again
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nkd,ipc(1),neul
      double precision dlat(3,2),bas(3,1),rmax(1),alat,
     .  eula(nbas,neul,3),qss(4),qlat(3,3)
C ... Local parameters
      integer n1,i1,i2,j1,j2,irep,jrep,nrep,ifi,i,j,i1mach,ii,ic,jc,m
      double precision z1(3),z2(3),angle(50),ddot,pi,
     .  dd0(3),dd(3),xd(3),dd1,sumrs,ovlprs,dmax,alfa,beta
      character outs*80
      logical cmdopt,a2bin
C     double precision rotm(3,3)

      if (neul .ne. 1) call rx('angtab not ready for l-dependent eula')

C --- Print out Euler angles ---
      if (ddot(3,qss,1,qss,1) .gt. 0) call
     .  awrit2('%x  Qss =%3:1;6d angle '//'%;6d',outs,80,
     .  -i1mach(2),qss,qss(4))
      print 344, (i, (eula(i,1,j), j=1,3), i=1,nbas)
  344 format('  ib      alpha       beta       gamma'/(i4,3f12.6))
      print *, ' '

C --- Print out Euler angles between each pair ---
      dmax = 0
      if (cmdopt('-dmax=',6,0,outs)) then
        j = 6
        call rxx(.not.a2bin(outs,dmax,4,0,' ',j,len(outs)),
     .           'angtab could not parse -dmax=..')
      endif
      pi = 4*datan(1d0)
      ifi = i1mach(2)
      n1 = 8
      nrep = (nbas-1)/n1 + 1
      do  10  irep = 1, nrep
      i1 = (irep-1)*n1+1
      i2 = min0(irep*n1,nbas)
      do  10  jrep = irep, nrep
        j1 = (jrep-1)*n1+1
        j2 = min0(jrep*n1,nbas)
        write(ifi,121) (j,j=j1,j2)
  121   format(3x,8i8)
        write(ifi,122) ('--------',j=j1,j2)
  122   format(7x,'---',8a8)
        do  11  i = i1, i2
        ic = ipc(i)
C        call eua2rm(eula(i,1,1),eula(i,1,2),eula(i,1,3),rotm)
CC       print 335, ((rotm(ii,jj),jj=1,3),ii=1,3)
CC  335  format((3f15.9))
C        z1(1) = rotm(3,1)
C        z1(2) = rotm(3,2)
C        z1(3) = rotm(3,3)

        alfa = eula(i,1,1)
        beta = eula(i,1,2)
        z1(1) = dcos(alfa)*dsin(beta)
        z1(2) = dsin(alfa)*dsin(beta)
        z1(3) = dcos(beta)

        do  12  j = j1, j2

C     ... Add any rotation from qss(z).
          do  16  m = 1, 3
   16     dd0(m) = bas(m,j) - bas(m,i)
          call shortn(dd0,dd,dlat,nkd)
C     ... This gives multiples of plat shortn took away
          do  17  ii = 1, 3
            xd(ii) =  (dd(1)-dd0(1))*qlat(1,ii)
     .              + (dd(2)-dd0(2))*qlat(2,ii)
     .              + (dd(3)-dd0(3))*qlat(3,ii)
   17     continue
C     ... assume here qss is along z, qlat(3) also
          alfa = eula(j,1,1) - 2*pi*xd(3)*qss(3)/qlat(3,3)
          beta = eula(j,1,2)
          z2(1) = dcos(alfa)*dsin(beta)
          z2(2) = dsin(alfa)*dsin(beta)
          z2(3) = dcos(beta)
Cc          call eua2rm(alfa,eula(j,1,2),eula(j,1,3),rotm)
C          z2(1) = rotm(3,1)
C          z2(2) = rotm(3,2)
C          z2(3) = rotm(3,3)
          angle(j) = dacos(max(-1d0,min(1d0,ddot(3,z1,1,z2,1))))
          if (angle(j) .gt. pi) angle(j) = angle(j) - 2*pi
   12   continue
        write(outs,130) i,(180/pi*angle(j),j=j1,j2)
  130   format(1x,i4,8f8.3)
        do  14  j = j1, j2
          jc = ipc(j)
          do  33  m = 1, 3
   33     dd0(m) = bas(m,j) - bas(m,i)
          call shortn(dd0,dd,dlat,nkd)
          do  34  m = 1, 3
   34     dd(m) = dd(m)*alat
          dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          sumrs = rmax(ic) + rmax(jc)
          ovlprs = sumrs - dd1
          if (i .eq. j) then
            outs(7+8*(j-j1):5+8+8*(j-j1)) = '  --'
          elseif (dmax .ne. 0) then
            if (dd1 .gt. dmax) outs(7+8*(j-j1):5+8+8*(j-j1)) = ' '
          elseif (100*ovlprs/dd1 .lt. -1) then
            outs(7+8*(j-j1):5+8+8*(j-j1)) = ' '
          endif
   14   continue
        call awrit0('%a',outs,-80,-i1mach(2))
   11 continue
   10 continue
      end

