      subroutine bzmesh(plat,qb,n1,n2,n3,lshft,g,ng,ipq,qp,wgt,nq,nqmx,
     .igstar,lpbc)
C- Divides the reciprocal lattice into microcells
C-----------------------------------------------------------------------
Ci Inputs:
Ci  plat     :primitive lattice vectors
Ci  n1,n2,n3 :no. divisions for the 3 recip. latt. vecs; (see Remarks)
Ci  g,ng     :symmetry point group operations, and number
Ci  wgt(1)   :if nonzero, nsgrp=input wgt(1) holds number of space group
Ci           :operations.  Setting wgt(1)>0 will cause bzmesh to
Ci           :flag which irreducible points have a star that consists
Ci           :of symmetry operations between nsgrp+1 and ng.  (Normally
Ci           :these ae symmetry operations which obtain not from the
Ci           :space group but from time-reversal symmetry.)  Thus the
Ci           :calling program can distinguish which irreducible points
Ci           :use this special symmetry.  bzmesh flags a special i
Ci           :points by returning wgt(i) as a negative number.
Ci  nqmx     :abort if number of k-points exceeds this maximum
Ci  igstar(0):if nonzero, make igstar.   See remarks for what is made
Ci  lshft    :logical switch, for each recip. latt. vec. :
Ci           :center the mesh through the origin (k=0)
Ci           :center the mesh straddling the origin --- no point at 0.
Ci  lpbc     :1s + 10s digit:
Ci           :  0 make usual 3D mesh
Ci           :  1 or 11 make 2D mesh
Ci           :100s digit
Ci           :  0 standard order of generation of q-points:
Ci                outer loop : i3=1..n3 inner loop :i1=1..n1
Ci           :  1 reverse order of generation of q-points:
Ci                outer loop : i1=1..n1 inner loop :i3=1..n3
Co Outputs:
Co   ipq   :ipq(i1,i2,i3) points to irreducible qp corresponding to
Co         :qp labelled (i1,i2,i3), where i1=1..n1, i2=1..n2, i3=1..n3
Co   wgt   :weight assigned to this qp (see Remarks)
Co         :Note the sign of wgt may be used to flag whether the star
Co         :of this irr qp comes at least in part from time-reversal
Co         :symmetry (see wgt(1) above).  Thus the true weight assigned
Co         :to this qp is abs(wgt(i)).  The sign of wgt(i) flags whether
Co         :this point has extra weighting from time-reversal symmetry.
Co   nq    :number of irreducible qp
Co   igstar:contains contains needed for the  mapping of qp
Co         :to another set.
Co         :*For a qp specified by the triplet (i1,i2,i3), let
Co         : let i = some irreducible qp with triplet (i1,i2,i3).
Co         : Define 'compressed index' j = i1 + ndmx*i2 + ndmx**2*i3
Co         : where ndmx = is defined in mxxyz()
Co         : Thus j contains triplet (i1,i2,i3) in compressed form
Co         :*If input igstar(0)=0, igstar is not touched.
Co         :*If input igstar(0)=2, igstar(i) = compressed index j
Co         : Thus igstar contains data for the 'inverse' of the ipq:
Co         : j = igstar(ipq(i1,i2,i3)) = compressed form of (i1,i2,i3)
Co         : In this mode igstar(1..nq) is generated.
Co         : Also igstar(0) is overwritten with n1+ndmx*n2+ndmx**2*n3
Co         :*If input igstar(0)=-2, igstar(i) contains the group
Co         : operation ig that maps points in the full BZ to the one of
Co         : the nq irreducible points.  Thus igstar contains
Co         : information needed to rotate a hamiltonian or wave function
Co         : from the irreducible qp any symmetry- equivalent point.
Co         : igstar() contains the group operation that rotates
Co         : irreducible qp ipq(i1,i2,i3) into qp(i1,i2,i3).
Co         : In this mode igstar(1..n1*n2*n3) is generated
Co         : with igstar(i1+(i2-1)*n1+(i3-1)*n1*n2) = ig.
Co         : Also igstar(0) is overwritten with -(n1+ndmx*n2+ndmx**2*n3)
Co   qb    :vectors of first microcell for input to BZINTS (see bzmsh0)
Cl Local variables
Cl  lsgrp  :if true, a symmetry op used to map this qp to another qp
Cl         :exceeded input wgt(1)=>start of q contains symmetry from
Cl         :not from space group but from time-reversal symmetry
Cr Remarks:
Cr  The reciprocal lattice is divided into n1*n2*n3 microcells which
Cr  are parallelipipeds with 8 corners.  The corners are nodes of the
Cr  k-space mesh in the whole reciprocal lattice unit cell.
Cr  Thus, for i1=1..n1, i2=1..n2, i3=1..n3 the qp(i1,i2,i3) are
Cr    q_k = (i1*ifac(1)-1)*qb(k,1) +
Cr          (i2*ifac(2)-1)*qb(k,2) +
Cr          (i3*ifac(3)-1)*qb(k,3)
Cr  where ifac is 1 or 2; see bzmsh0.
Cr
Cr  Some of the qp will be symmetry-related, leaving nq irreducible
Cr  k-points, which are returned in qp(1..3,j), j=1,nq.
Cr
Cr  wgt(j) contains the sampling weight associated with qp(j), i.e.
Cr    2/(n1*n2*n3) * no. points equivalent to qp(j); factor 2 for spin.
Cr
Cr  ipq(i1,i2,i3) marks which irreducible qp to which each point in the
Cr  full BZ belongs: point (i1,i2,i3) is equivalent to irreducible
Cr  point ipq(i1,i2,i3).
Cr
Cu Updates
Cu   09 Jan 09 Package calculation of ndmx into mxxyz, to circumvent
Cu             compiler bugs
Cu   09 Jan 03 Can pass ng=0
Cu   15 Sep 02 Use sign of wgt to flag which irr points contain
Cu             equivalent points from time-reversal symmetry
Cu   21 Jul 02 Further changes to 18 Jun revisions.
Cu   18 Jun 02 New 100s digit in lpbc option to change order in
Cu             generation of qpoints.  Stripped some options in igstar.
Cu   23 Nov 01 Bug fix for long unit cells
Cr   19 Nov 97 (WRL) added lpbc option, projecting qp to 2D
C-----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lshft(3)
      integer n1,n2,n3,nqmx,ng,nq,igstar(0:*),ipq(n1,n2,n3),lpbc
      double precision qb(3,3),wgt(nqmx),plat(3,3),qp(3,nqmx),g(9,1)
C ... Local parameters
      logical lsgrp
      integer i1,i2,i3,ifac(3),ig,igcnt,ii,ii1,ii2,ii3,ipr,iq,
     .is(3),iwgt,j1,j2,j3,lgunit,lpbc01,lpbc2,lstar,m1,m2,m3,ndmx,nn1,
     .nn2,nn3,stdo,nsgrp,mxxyz
      double precision w0,x1,x2,x3,swgt,v(3),v1(3),rb(3,3)
C ... External calls
      external asymop,awrit5,awrit6,bzmsh0,dinv33,dpzero,fexit,getpr,
     .grpop,iinit,projql,rx,rxx
      character*1 chr(0:2)
#if DEBUG
      character*30 sg
      double precision ag(3)
#endif
#if QP0
C     --- use only points centered within dq from q0
      integer ifi,fopn
      double precision q0(3),dq,ddot
C      data q0 /0.370925d0,-0.401855d0,0.370925d0/,dq/.15d0/
      data q0 /0.339995d0,0.339995d0,0.339995d0/,dq/.15d0/
      integer k,m
#endif

      call getpr(ipr)
      stdo = lgunit(1)
      ndmx = mxxyz()
      lstar = igstar(0)
      print *,'lstar xxx=',lstar
      lpbc01 = mod(lpbc,100)
      lpbc2  = lpbc/100
      nsgrp = wgt(1)
      chr(2) = ' '
      chr(0) = '*'
      if (nsgrp .ne. wgt(1)) call rx('bzmesh: invalid input wgt(1)')
C     lpbc2 = 1
      if (lstar .ne. 0 .and. max(n1,n2,n3) .gt. ndmx)
     .call rx('bzmesh: too many divisions to accomodate ndmx')
      if (min(n1,n2,n3) .lt. 1) 
     .call rx('bzmesh: improper specification of k-mesh')

      call bzmsh0(plat,lshft,lpbc01,n1,n2,n3,is,ifac,rb,qb)
#if DEBUG
      if (ipr .ge. 80) print 356
  356 format(/' BZMESH: qp mapping'/' i1..i3',25x,'qp',16x,'iq   ig g')
      call dpzero(ag,3)
#endif
      m1 = n1*ifac(1)
      m2 = n2*ifac(2)
      m3 = n3*ifac(3)
      call iinit(ipq,n1*n2*n3)
      w0 = 2d0/(n1*n2*n3)
      nq = 0
      swgt = 0d0
      igcnt = 0
      nn1 = 6*m1
      nn2 = 6*m2
      nn3 = 6*m3

C --- For each of (n1*n2*n3) qp, find irreducible set ---
C     Loops written without 'do' construct for freedom in choosing
C     loop order
C     do  20  i3 = 0, n3-1
C     do  20  i2 = 0, n2-1
C     do  20  i1 = 0, n1-1
      i1 = 0
      i3 = 0
   23 continue
      i2 = 0
   22 continue
      if (lpbc2 .eq. 0) then
        i1 = 0
      else
        i3 = 0
      endif
   21 continue

C     print 973, i1,i2,i3,igcnt
C 973 format(5i4)

C   ... Add qp to list if not flagged as symmetry-related to a prior
      if (ipq(i1+1,i2+1,i3+1) .eq. 0) then
        ii1 = i1*ifac(1)+is(1)
        ii2 = i2*ifac(2)+is(2)
        ii3 = i3*ifac(3)+is(3)
        v(1) = ii1*qb(1,1) + ii2*qb(1,2) + ii3*qb(1,3)
        v(2) = ii1*qb(2,1) + ii2*qb(2,2) + ii3*qb(2,3)
        v(3) = ii1*qb(3,1) + ii2*qb(3,2) + ii3*qb(3,3)
C         call prmx('kpt',v,3,3,1)
#if QP0
        call dcopy(3,v,1,v1,1)
        call daxpy(3,-1d0,q0,1,v1,1)
C          print 156, ddot(3,v1,1,v1,1), v1
C  156     format(4f15.6)
        if (ddot(3,v1,1,v1,1) .gt. dq*dq) goto 20
#endif

C     --- Mark each qp in the star of q as equivalent to this qp ---
        iwgt = 0
        lsgrp = .false.
        call dcopy(3,v,1,v1,1)
        do  25  ig = 1, max(ng,1)
          if (ng .gt. 0) call grpop(v,v1,g,ig)
          x1 = v1(1)*rb(1,1) + v1(2)*rb(2,1) + v1(3)*rb(3,1) - is(1)
          x2 = v1(1)*rb(1,2) + v1(2)*rb(2,2) + v1(3)*rb(3,2) - is(2)
          if (lpbc01 .eq. 0) then
            x3 = v1(1)*rb(1,3) + v1(2)*rb(2,3) + v1(3)*rb(3,3) - is(3)
          else
            x3 = 0
          endif
          j1 = idnint(x1)
          j2 = idnint(x2)
          j3 = idnint(x3)
          if (max(dabs(x1-j1),dabs(x2-j2),dabs(x3-j3)) .gt. 1d-4) then
C              print *, 'ii1,ii2,ii3=',ii1,ii2,ii3
C              print *, 'x1,x2,x3=',sngl(x1),sngl(x2),sngl(x3)
C              print *, 'ig,qb,rb,v,g,v1=',ig
C              call prm('(3f12.6)',qb,3,3)
C              call prm('(3f12.6)',rb,3,3)
C              call prm('(3f12.6)',v,3,1)
C              call prm('(3f12.6)',g(9*ig-8),3,3)
C              call prm('(3f12.6)',v1,3,1)
            call awrit2(' qp%3:1,3;3d -> %3:1,3;3d is not'//
     .      ' on k-mesh',' ',80,stdo,v,v1)
            call rx('BZMESH: symops incompatible with this mesh')
          endif
C           call prmx('rotated kpt',v1,3,3,1)
C            if (ig .ne. 1 .and. ipr .ge.
C            call awrit2(' qp%3:1,3;3d rotated from %3:1,3;3d',
C     .        ' ',80,stdo,v1,v)
C ..        scale shifted point or discard if shifted off mesh
          if (lshft(1) .and. mod(abs(j1),2).eq.1) goto 25
          if (lshft(2) .and. mod(abs(j2),2).eq.1) goto 25
          if (lshft(3) .and. mod(abs(j3),2).eq.1) goto 25
          if (lshft(1)) j1 = j1/2
          if (lshft(2)) j2 = j2/2
          if (lshft(3)) j3 = j3/2
C ...       Ensure (j1,j2,j3) in first quadrant of Q
          j1 = mod(j1+2*nn1,n1) + 1
          j2 = mod(j2+2*nn2,n2) + 1
          j3 = mod(j3+2*nn3,n3) + 1
          call rxx(j1.le.0.or.j2.le.0.or.j3.le.0,'neg j in bzmesh')

          if (ipq(j1,j2,j3) .eq. 0) then
            ipq(j1,j2,j3) = nq+1
            iwgt = iwgt+1
            if (ig .gt. nsgrp .and. nsgrp .gt. 0) lsgrp = .true.
            igcnt = igcnt+1
            if (lstar .eq. 2) then
              igstar(nq+1) = j1 + ndmx*j2 + ndmx**2*j3
            elseif (lstar .eq. -2) then
              ii = j1 + (j2-1)*n1 + (j3-1)*n1*n2
              igstar(ii) = ig
            elseif (lstar .ne. 0) then
              call rx('bzmesh: bad igstar(0)')
            endif
C              OLD conventions
C              if (lstar .eq. 1) then
C                igstar(igcnt) = j1 + ndmx*j2 + ndmx**2*j3
C              elseif (lstar .eq. 2) then
C                igstar(nq+1) = j1 + ndmx*j2 + ndmx**2*j3
C              elseif (lstar .eq. -1) then
C                igstar(igcnt) = ig
C              elseif (lstar .eq. -2) then
C                ii = j1 + (j2-1)*n1 + (j3-1)*n1*n2
C                igstar(ii) = ig
C              endif
#if DEBUG
            if (ipr .ge. 80 .and. ng .gt. 0) then
              call dpzero(ag,3)
              call asymop(g(1,ig),ag,' ',sg)
              call awrit6(' (%i,%i,%i)%15p%3:1,6;11D'//
     .        ' %,5i  %,3i '//sg//'%a',' ',80,stdo,j1,j2,j3,
     .        v1,ipq(j1,j2,j3),ig)
            endif
#endif
          endif
          call rxx(j1.lt.0.or.j2.lt.0.or.j3.lt.0,'neg j in bzmesh')
   25   continue
        nq = nq+1
        qp(1,nq) = v(1)
        qp(2,nq) = v(2)
        qp(3,nq) = v(3)
        wgt(nq) = iwgt*w0
        if (lsgrp) wgt(nq) = -iwgt*w0
        swgt = swgt + abs(wgt(nq))
      endif

C     End-of-loop for i1,i2,i3
C  20 continue
      if (lpbc2 .eq. 0) then
        i1 = i1+1
        if (i1 .lt. n1) goto 21
      else
        i3 = i3+1
        if (i3 .lt. n3) goto 21
      endif
      i2 = i2+1
      if (i2 .lt. n2) goto 22
      if (lpbc2 .eq. 0) then
        i3 = i3+1
        if (i3 .lt. n3) goto 23
      else
        i1 = i1+1
        if (i1 .lt. n1) goto 23
      endif
C ... Done accumulating inequivalent qp

      if (ipr .ge. 20) then
        call awrit6(' BZMESH:  %i irreducible QP from %i'//
     .  ' ( %i %i %i )  shift=%3:1l',' ',80,
     .  stdo,nq,n1*n2*n3,n1,n2,n3,lshft)
        call awrit5('bz  nq %i ( %i %i %i )  shift %3:1l',' ',80,
     .  lgunit(2),nq,n1,n2,n3,lshft)
      endif
#if  QP0
      call awrit3('          %i points within dq=%1;6d  swgt=%1;6d',
     .' ',80,stdo,igcnt,dq,swgt)
      ifi = fopn('QPTS')
      write(ifi,333) nq, n1,n2,n3,0
  333 format(4i5,i8)
      do  70  iq = 1, nq
        write(ifi,334) iq,(qp(m,iq),m=1,3), wgt(iq)
   70 continue
  334 format(i5,1p,4d18.10)
#else
      if (lstar .ne. 0) igstar(0) = n1 + ndmx*n2 + ndmx**2*n3
      if (lstar .lt. 0) igstar(0) = -igstar(0)
      if (igcnt .ne. n1*n2*n3) call rx('bug in bzmesh')
      if (dabs(swgt-2) .gt. 1.d-9) call fexit(-1,111,
     .' Exit -1 BZMESH: QP weights sum to %1;6d but must sum to 2',swgt)
#endif
      if (ipr .ge. 50) then
        write(stdo,663)
  663   format(14x,'Qx',10x,'Qy',10x,'Qz',6x,'Multiplicity    Weight')
        do  51  iq = 1, nq
          ii = 1+dsign(1d0,wgt(iq))
          iwgt = abs(wgt(iq)/w0) + .1d0
          write(stdo,661)
     .    iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,chr(ii),abs(wgt(iq))
   51   continue
  661   format(i5,2x,3f12.6,i10,1x,a,f14.6)
      endif
      end

      subroutine bzmsh0(plat,lshft,lpbc,n1,n2,n3,is,ifac,rb,qb)
C- Setup for a uniform mesh in the Brillouin zone
C-----------------------------------------------------------------------
Ci   plat:      primitive lattice vectors
Ci   n1,n2,n3:  number of divisions along the three R.L.V
Ci   lshft:     F center mesh points on the origin for i'th R.L.V
Ci              T set mesh points off center from the origin
Ci   lpbc:      0 make 3D mesh
Ci              1 or 11 make 2D mesh
Co Outputs:
Co   is(i)      0 mesh points centered at the origin for i'th axis
Co              1 mesh points off-centered from the origin for i'th axis
Co   ifac(i)    1 if not shifted, 2 if shifted for i'th axis; see qb.
Co   rb:        a microcell in the Wigner Seitz cell.
Co   qb:        a microcell in the Brillouin zone
Co              This, together with ifac provides information how to
Co              generate the actual q-point from a triplet of integers
Co              specifying a point on the mesh.  Given a
Co              triplet (j_1,j_2,j_3) of ipq, the i'th component q_i is
Co                 q_i(j_1,j_2,j_3) = sum_n (j'_n)*qb(i,n),  with
Co                 j'_n = j_n*ifac(n)-1
Cu Updates
Cu   19 Nov 97 added lpbc option
C-----------------------------------------------------------------------
C     implicit none
      logical lshft(3)
      integer n1,n2,n3,is(3),ifac(3),lpbc
      double precision plat(3,3),rb(3,3),qb(3,3),qlat(3,3),vol
      integer k,m,m1,m2,m3,iprint,lgunit,stdo

      stdo = lgunit(1)
      is(1)   = 0
      is(2)   = 0
      is(3)   = 0
      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      if (lshft(1)) then
        is(1) = 1
        ifac(1) = 2
      endif
      if (lshft(2)) then
        is(2) = 1
        ifac(2) = 2
      endif
      if (lshft(3)) then
        is(3) = 1
        ifac(3) = 2
      endif
      m1 = n1*ifac(1)
      m2 = n2*ifac(2)
      m3 = n3*ifac(3)
      call dinv33(plat,1,qlat,vol)
      if (lpbc .eq. 1 .or. lpbc .eq. 11) call projql(qlat)
      if (iprint() .gt. 80) then
        write(stdo,351)
        do  35  k = 1, 3
          write(stdo,350) (plat(m,k),m=1,3),(qlat(m,k),m=1,3)
   35   continue
      endif
  350 format(3f10.5,5x,3f10.5)
  351 format(' BZMESH : ',5X,'Plat',31X,'Qlat')
      do  8  m = 1, 3
        qb(m,1) = qlat(m,1)/m1
        qb(m,2) = qlat(m,2)/m2
        qb(m,3) = qlat(m,3)/m3
        rb(m,1) = plat(m,1)*m1
        rb(m,2) = plat(m,2)*m2
        rb(m,3) = plat(m,3)*m3
    8 continue

      end
      integer function mxxyz()
C- Return maximum integer whose cube fits into integer word
C  Package as a subroutine to avoid compiler bugs, e.g. result changing
C  as compiler switches change.
C  e.g. intel ifort v 11, for x86_64, these should produce same results:
C  ifort -g -cm -axW -WB -c bzmesh.f
C  ifort -g -cm -c bzmesh.f
C     implicit none
      integer i1mach

      mxxyz = (dble((i1mach(9)))/2.01d0) ** (1d0/3d0)
      end

