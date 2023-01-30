subroutine bzmesh(plat,qb,n1,n2,n3,lshft,g,ng,ipq,qp,wgt,nq,nqmx,igstar,lpbc)
  use m_ext,only: sname   !file extension. Open a file like file='ctrl.'//trim(sname)
  use m_lgunit,only:stdo
  !- Divides the reciprocal lattice into microcells
  !-----------------------------------------------------------------------
  !i Inputs:
  !i  plat     :primitive lattice vectors
  !i  n1,n2,n3 :no. divisions for the 3 recip. latt. vecs; (see Remarks)
  !i  g,ng     :symmetry point group operations, and number
  !i  wgt(1)   :if nonzero, nsgrp=input wgt(1) holds number of space group
  !i           :operations.  Setting wgt(1)>0 will cause bzmesh to
  !i           :flag which irreducible points have a star that consists
  !i           :of symmetry operations between nsgrp+1 and ng.  (Normally
  !i           :these ae symmetry operations which obtain not from the
  !i           :space group but from time-reversal symmetry.)  Thus the
  !i           :calling program can distinguish which irreducible points
  !i           :use this special symmetry.  bzmesh flags a special i
  !i           :points by returning wgt(i) as a negative number.
  !i  nqmx     :abort if number of k-points exceeds this maximum
  !i  igstar(0):if nonzero, make igstar.   See remarks for what is made
  !i  lshft    :logical switch, for each recip. latt. vec. :
  !i           :center the mesh through the origin (k=0)
  !i           :center the mesh straddling the origin --- no point at 0.
  !i  lpbc     :1s + 10s digit:
  !i           :  0 make usual 3D mesh
  !i           :  1 or 11 make 2D mesh
  !i           :100s digit
  !i           :  0 standard order of generation of q-points:
  !i                outer loop : i3=1..n3 inner loop :i1=1..n1
  !i           :  1 reverse order of generation of q-points:
  !i                outer loop : i1=1..n1 inner loop :i3=1..n3
  !o Outputs:
  !o   ipq   :ipq(i1,i2,i3) points to irreducible qp corresponding to
  !o         :qp labelled (i1,i2,i3), where i1=1..n1, i2=1..n2, i3=1..n3
  !o   wgt   :weight assigned to this qp (see Remarks)
  !o         :Note the sign of wgt may be used to flag whether the star
  !o         :of this irr qp comes at least in part from time-reversal
  !o         :symmetry (see wgt(1) above).  Thus the true weight assigned
  !o         :to this qp is abs(wgt(i)).  The sign of wgt(i) flags whether
  !o         :this point has extra weighting from time-reversal symmetry.
  !o   nq    :number of irreducible qp
  !o   igstar:contains contains needed for the  mapping of qp
  !o         :to another set.
  !o         :*For a qp specified by the triplet (i1,i2,i3), let
  !o         : let i = some irreducible qp with triplet (i1,i2,i3).
  !o         : Define 'compressed index' j = i1 + ndmx*i2 + ndmx**2*i3
  !o         : where ndmx = is defined in mxxyz()
  !o         : Thus j contains triplet (i1,i2,i3) in compressed form
  !o         :*If input igstar(0)=0, igstar is not touched.
  !o         :*If input igstar(0)=2, igstar(i) = compressed index j
  !o         : Thus igstar contains data for the 'inverse' of the ipq:
  !o         : j = igstar(ipq(i1,i2,i3)) = compressed form of (i1,i2,i3)
  !o         : In this mode igstar(1..nq) is generated.
  !o         : Also igstar(0) is overwritten with n1+ndmx*n2+ndmx**2*n3
  !o         :*If input igstar(0)=-2, igstar(i) contains the group
  !o         : operation ig that maps points in the full BZ to the one of
  !o         : the nq irreducible points.  Thus igstar contains
  !o         : information needed to rotate a hamiltonian or wave function
  !o         : from the irreducible qp any symmetry- equivalent point.
  !o         : igstar() contains the group operation that rotates
  !o         : irreducible qp ipq(i1,i2,i3) into qp(i1,i2,i3).
  !o         : In this mode igstar(1..n1*n2*n3) is generated
  !o         : with igstar(i1+(i2-1)*n1+(i3-1)*n1*n2) = ig.
  !o         : Also igstar(0) is overwritten with -(n1+ndmx*n2+ndmx**2*n3)
  !o   qb    :vectors of first microcell for input to BZINTS (see bzmsh0)
  !l Local variables
  !l  lsgrp  :if true, a symmetry op used to map this qp to another qp
  !l         :exceeded input wgt(1)=>start of q contains symmetry from
  !l         :not from space group but from time-reversal symmetry
  !r Remarks:
  !r  The reciprocal lattice is divided into n1*n2*n3 microcells which
  !r  are parallelipipeds with 8 corners.  The corners are nodes of the
  !r  k-space mesh in the whole reciprocal lattice unit cell.
  !r  Thus, for i1=1..n1, i2=1..n2, i3=1..n3 the qp(i1,i2,i3) are
  !r    q_k = (i1*ifac(1)-1)*qb(k,1) +
  !r          (i2*ifac(2)-1)*qb(k,2) +
  !r          (i3*ifac(3)-1)*qb(k,3)
  !r  where ifac is 1 or 2; see bzmsh0.
  !r
  !r  Some of the qp will be symmetry-related, leaving nq irreducible
  !r  k-points, which are returned in qp(1..3,j), j=1,nq.
  !r
  !r  wgt(j) contains the sampling weight associated with qp(j), i.e.
  !r    2/(n1*n2*n3) * no. points equivalent to qp(j); factor 2 for spin.
  !r
  !r  ipq(i1,i2,i3) marks which irreducible qp to which each point in the
  !r  full BZ belongs: point (i1,i2,i3) is equivalent to irreducible
  !r  point ipq(i1,i2,i3).
  !r
  !u Updates
  !u   09 Jan 09 Package calculation of ndmx into mxxyz, to circumvent
  !u             compiler bugs
  !u   09 Jan 03 Can pass ng=0
  !u   15 Sep 02 Use sign of wgt to flag which irr points contain
  !u             equivalent points from time-reversal symmetry
  !u   21 Jul 02 Further changes to 18 Jun revisions.
  !u   18 Jun 02 New 100s digit in lpbc option to change order in
  !u             generation of qpoints.  Stripped some options in igstar.
  !u   23 Nov 01 Bug fix for long unit cells
  !r   19 Nov 97 (WRL) added lpbc option, projecting qp to 2D
  !-----------------------------------------------------------------------
  implicit none
  logical :: lshft(3)
  integer :: n1,n2,n3,nqmx,ng,nq,igstar(0:*),ipq(n1,n2,n3),lpbc
  double precision :: qb(3,3),wgt(nqmx),plat(3,3),qp(3,nqmx),g(3,3,*)
  logical :: lsgrp
  integer :: i1,i2,i3,ifac(3),ig,igcnt,ii,ii1,ii2,ii3,ipr,iq, &
       is(3),iwgt,j1,j2,j3,lgunit,lpbc01,lpbc2,lstar,m1,m2,m3,ndmx,nn1, &
       nn2,nn3,nsgrp,mxxyz
  double precision :: w0,x1,x2,x3,swgt,v(3),v1(3),rb(3,3)
  character(1) :: chr(0:2)
  real(8),parameter:: tolq=1d-3
  call getpr(ipr)
  !      stdo = lgunit(1)
  ndmx = mxxyz()
  lstar = igstar(0)
  lpbc01 = mod(lpbc,100)
  lpbc2  = lpbc/100
  nsgrp = wgt(1)
  chr(2) = ' '
  chr(0) = '*'
  if (nsgrp /= wgt(1)) call rx('bzmesh: invalid input wgt(1)')
  if (lstar /= 0 .AND. max(n1,n2,n3) > ndmx) &
       call rx('bzmesh: too many divisions to accomodate ndmx')
  if (min(n1,n2,n3) < 1) call rx('bzmesh: improper specification of k-mesh')
  call bzmsh0(plat,lshft,lpbc01,n1,n2,n3,is,ifac,rb,qb)
  m1 = n1*ifac(1)
  m2 = n2*ifac(2)
  m3 = n3*ifac(3)
  ipq=0 !call iinit(ipq,n1*n2*n3)
  w0 = 2d0/(n1*n2*n3)
  nq = 0
  swgt = 0d0
  igcnt = 0
  nn1 = 6*m1
  nn2 = 6*m2
  nn3 = 6*m3

  ! --- For each of (n1*n2*n3) qp, find irreducible set ---
  i1 = 0
  i3 = 0
  !   23 continue
  do 23  !main loop here
     i2 = 0
22   continue
     if (lpbc2 == 0) then
        i1 = 0
     else
        i3 = 0
     endif
21   continue
     !   ... Add qp to list if not flagged as symmetry-related to a prior
     if (ipq(i1+1,i2+1,i3+1) == 0) then
        ii1 = i1*ifac(1)+is(1)
        ii2 = i2*ifac(2)+is(2)
        ii3 = i3*ifac(3)+is(3)
        v(1) = ii1*qb(1,1) + ii2*qb(1,2) + ii3*qb(1,3)
        v(2) = ii1*qb(2,1) + ii2*qb(2,2) + ii3*qb(2,3)
        v(3) = ii1*qb(3,1) + ii2*qb(3,2) + ii3*qb(3,3)
        !     --- Mark each qp in the star of q as equivalent to this qp ---
        iwgt = 0
        lsgrp = .false.
        call dcopy(3,v,1,v1,1)
        do  25  ig = 1, max(ng,1)
           if (ng > 0) v1=matmul(g(:,:,ig),v) !call grpop(v,v1,g,ig)
           x1 = v1(1)*rb(1,1) + v1(2)*rb(2,1) + v1(3)*rb(3,1) - is(1)
           x2 = v1(1)*rb(1,2) + v1(2)*rb(2,2) + v1(3)*rb(3,2) - is(2)
           if (lpbc01 == 0) then
              x3 = v1(1)*rb(1,3) + v1(2)*rb(2,3) + v1(3)*rb(3,3) - is(3)
           else
              x3 = 0
           endif
           j1 = idnint(x1)
           j2 = idnint(x2)
           j3 = idnint(x3)
           if (max(dabs(x1-j1),dabs(x2-j2),dabs(x3-j3)) > tolq) then
              !            call awrit2(' qp%3:1,3;3d -> %3:1,3;3d is not on k-mesh',' ',80,stdo,v,v1)
              write(stdo,"(a,3f9.4,' ',3f9.4)") ' qp mapped to is not on k-mesh',v,v1
              write(stdo,"(a,3f9.4,' ',3i5)")   '             x j=',x1,x2,x3,j1,j2,j3
              call rx('BZMESH: symops incompatible with this mesh')
           endif
           ! ..        scale shifted point or discard if shifted off mesh
           if (lshft(1) .AND. mod(abs(j1),2) == 1) goto 25
           if (lshft(2) .AND. mod(abs(j2),2) == 1) goto 25
           if (lshft(3) .AND. mod(abs(j3),2) == 1) goto 25
           if (lshft(1)) j1 = j1/2
           if (lshft(2)) j2 = j2/2
           if (lshft(3)) j3 = j3/2
           ! ...       Ensure (j1,j2,j3) in first quadrant of Q
           j1 = mod(j1+2*nn1,n1) + 1
           j2 = mod(j2+2*nn2,n2) + 1
           j3 = mod(j3+2*nn3,n3) + 1
           call rxx(j1.le.0.or.j2.le.0.or.j3.le.0,'neg j in bzmesh')
           if (ipq(j1,j2,j3) == 0) then
              ipq(j1,j2,j3) = nq+1
              iwgt = iwgt+1
              if (ig > nsgrp .AND. nsgrp > 0) lsgrp = .TRUE. 
              igcnt = igcnt+1
              if (lstar == 2) then
                 igstar(nq+1) = j1 + ndmx*j2 + ndmx**2*j3
              elseif (lstar == -2) then
                 ii = j1 + (j2-1)*n1 + (j3-1)*n1*n2
                 igstar(ii) = ig
              elseif (lstar /= 0) then
                 call rx('bzmesh: bad igstar(0)')
              endif
           endif
           call rxx(j1.lt.0.or.j2.lt.0.or.j3.lt.0,'neg j in bzmesh')
25      enddo
        nq = nq+1
        qp(1,nq) = v(1)
        qp(2,nq) = v(2)
        qp(3,nq) = v(3)
        wgt(nq) = iwgt*w0
        if (lsgrp) wgt(nq) = -iwgt*w0
        swgt = swgt + abs(wgt(nq))
     endif
     !     End-of-loop for i1,i2,i3
     !  20 continue
     if (lpbc2 == 0) then
        i1 = i1+1
        if (i1 < n1) goto 21
     else
        i3 = i3+1
        if (i3 < n3) goto 21
     endif
     i2 = i2+1
     if (i2 < n2) goto 22
     if (lpbc2 == 0) then
        i3 = i3+1
        if (i3 < n3) cycle !goto 23
     else
        i1 = i1+1
        if (i1 < n1) cycle !goto 23
     endif
     exit
23 enddo
  ! ... Done accumulating inequivalent qp
  if (ipr>= 20) then
     write(stdo,"(a,i5,a,3i4,a,3l)") " BZMESH:  ",nq," irreducible QP from ",n1,n2,n3," shift=",lshft
  endif
  if (lstar /= 0) igstar(0) = n1 + ndmx*n2 + ndmx**2*n3
  if (lstar < 0) igstar(0) = -igstar(0)
  if (igcnt /= n1*n2*n3) call rx('bug in bzmesh')
  if (dabs(swgt-2) > 1.d-9) call rx1('BZMESH: QP weights sum to ',swgt)
  if (ipr >= 50) then
     write(stdo,663)
663  format(14x,'Qx',10x,'Qy',10x,'Qz',6x,'Multiplicity    Weight')
     do  51  iq = 1, nq
        ii = 1+dsign(1d0,wgt(iq))
        iwgt = abs(wgt(iq)/w0) + .1d0
        write(stdo,661) &
             iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,chr(ii),abs(wgt(iq))
51   enddo
661  format(i5,2x,3f12.6,i10,1x,a,f14.6)
  endif
end subroutine bzmesh

subroutine bzmsh0(plat,lshft,lpbc,n1,n2,n3,is,ifac,rb,qb)
  use m_lgunit,only:stdo
  !- Setup for a uniform mesh in the Brillouin zone
  !-----------------------------------------------------------------------
  !i   plat:      primitive lattice vectors
  !i   n1,n2,n3:  number of divisions along the three R.L.V
  !i   lshft:     F center mesh points on the origin for i'th R.L.V
  !i              T set mesh points off center from the origin
  !i   lpbc:      0 make 3D mesh
  !i              1 or 11 make 2D mesh
  !o Outputs:
  !o   is(i)      0 mesh points centered at the origin for i'th axis
  !o              1 mesh points off-centered from the origin for i'th axis
  !o   ifac(i)    1 if not shifted, 2 if shifted for i'th axis; see qb.
  !o   rb:        a microcell in the Wigner Seitz cell.
  !o   qb:        a microcell in the Brillouin zone
  !o              This, together with ifac provides information how to
  !o              generate the actual q-point from a triplet of integers
  !o              specifying a point on the mesh.  Given a
  !o              triplet (j_1,j_2,j_3) of ipq, the i'th component q_i is
  !o                 q_i(j_1,j_2,j_3) = sum_n (j'_n)*qb(i,n),  with
  !o                 j'_n = j_n*ifac(n)-1
  !u Updates
  !u   19 Nov 97 added lpbc option
  !-----------------------------------------------------------------------
  !     implicit none
  logical :: lshft(3)
  integer :: n1,n2,n3,is(3),ifac(3),lpbc
  double precision :: plat(3,3),rb(3,3),qb(3,3),qlat(3,3),vol
  integer :: k,m,m1,m2,m3,iprint
  !      print *,'bzmsh0:',plat,n1,n2,n3,lshft,lpbc
  !      stdo = lgunit(1)
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
  if (lpbc == 1 .OR. lpbc == 11) call projql(qlat)
  if (iprint() > 80) then
     write(stdo,351)
     do  35  k = 1, 3
        write(stdo,350) (plat(m,k),m=1,3),(qlat(m,k),m=1,3)
35   enddo
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
8 enddo
  !      print *,'end of bzmesh0'
end subroutine bzmsh0
integer function mxxyz()
  !- Return maximum integer whose cube fits into integer word
  !  Package as a subroutine to avoid compiler bugs, e.g. result changing
  !  as compiler switches change.
  !  e.g. intel ifort v 11, for x86_64, these should produce same results:
  !  ifort -g -cm -axW -WB -c bzmesh.f
  !  ifort -g -cm -c bzmesh.f
  !     implicit none
  integer :: i1mach
  mxxyz = (dble((i1mach(9)))/2.01d0) ** (1d0/3d0)
end function mxxyz

subroutine projql(qlat)
  use m_lmfinit,only: stdo
  !- Project 3D reciprocal lattice vecs on surface (2D recip lattice vecs)
  ! ----------------------------------------------------------------------
  !i Input qlat(3,3)  3D reciprocal lattice vectors
  !o Output qlat(3,3) 2D reciprocal lattice vectors
  !r if bi i=1,3 are 3D vectors
  !r    Bi i=1,2 are 2D vectors
  !r    Bi(u)=bi(u)-(Bi.B3) B3(u)/(B3.B3) for u=x,y,z
  ! ----------------------------------------------------------------------
  implicit none
  double precision :: qlat(3,3),b3b3,bib3
  integer :: i,j,iprint
  b3b3 = qlat(1,3)**2 + qlat(2,3)**2 + qlat(3,3)**2
  do  20  i = 1, 2
     bib3 = qlat(1,i)*qlat(1,3) + qlat(2,i)*qlat(2,3) &
          + qlat(3,i)*qlat(3,3)
     do  21  j = 1, 3
        qlat(j,i) = qlat(j,i) - qlat(j,3)*bib3/b3b3
21   enddo
20 enddo
  ! ... Ignore qlat(3,.) by making it zero
  !      do  22  j = 1, 3
  !   22 qlat(j,3)=0d0
  if (iprint() >= 20) then
     ! hangenglob        stdo = nglob('stdo')
     !        stdo = globalvariables%stdo
     write(stdo,*) ' 2D-reciprocal vectors:'
     do  23  i = 1, 2
        write(stdo,222) (qlat(j,i),j=1,3)
23   enddo
  endif
222 format(3f12.6)
end subroutine projql

