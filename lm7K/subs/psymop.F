      subroutine psymop(t,plat,g,ag,ng)
C- Parse symbolic representation of symmetry group operations
C ----------------------------------------------------------------------
Ci Inputs:
Ci   t,nt  string of symmetry operations, separated by spaces
Ci   plat  lattice vectors that scale translation part ag
Ci         (if needed, i.e. if translation specified by '::')
Co Outputs:
Co   g,ng  group op (3x3 matrix) for each input, and number
Cr Remarks:
Cr   Symbols have two parts, first the point group part, followed
Cr   By an optional translation.  The point group part has the form
Cr   O(nx,ny,nz) where O is one of M, I or Rj for mirror, inversion
Cr   and j-fold rotations, respectively, and nx,ny,nz are a triplet
Cr   of indices specifying the axis of operation.
Cr   (nx,ny,nz) is one of (1,0,0), (0,1,0), (0,0,1) and (1,1,1),
Cr   it can be abbreviated as x,X, y,Y, z,Z and d,D, respectively.
Cr   Also permissible are products, eg I*R4X.
Cr   The translation is also of the form (n1,n2,n3)
Cr   Example: the following input
Cr     R3D(0,0,0) Mx R2(1/2,sqrt(3)/2,0)(pi,0,0) my*i'
Cr   is nonsensical, but permissible and generates four group ops.
Cr   10 Jan 1997 now generates g=transpose of prior versions.
C ----------------------------------------------------------------------
C     implicit none
      character*(*) t
      double precision plat(3,3),g(9,1),h(9),hh(9),ag(3,1),vec(3)
      integer nt,ng,i
      logical parsvc,flgp
      character*1:: leftp='('
C --- Do until no more symbolic representation, do ---
      nt = len(t)
      ng = 0
      i = 0
   90 call skipbl(t,nt,i)
      if (i .ge. nt) return
      ng = ng+1
      call parsop(t,i,g(1,ng))
      if (t(i+1:i+1) .eq. '*') then
        i = i+1
        call parsop(t,i,h)
        call grpprd(g(1,ng),h,hh)
C       call dmpy(g(1,ng),3,1,h,3,1,hh,3,1,3,3,3)
C       call dvcpy(hh,1,g(1,ng),1,9)
        call dcopy(9,hh,1,g(1,ng),1)
      endif
      call dpzero(ag(1,ng),3)
C ... Compatibility with old :T(x,y,z)
      if (t(i+1:i+2) .eq. ':T' .or. t(i+1:i+2) .eq. ':t') i=i+2
C ... Compatibility with ::(x,y,z)
      flgp = .false.
      if (t(i+1:i+2) .eq. '::') then
        flgp = .true.
        i=i+2
      elseif (t(i+1:i+1) .eq. ':') then
        i=i+1
      endif
      if (t(i+1:i+1) .eq. leftp) then
        if (.not. parsvc(-1,t,i,ag(1,ng)))
     .  call fexit(-1,111,' Exit -1 PSYMOP: '//
     .  'failed to parse translation, ig=%i',ng)
        if (flgp) then
          call dcopy(3,ag(1,ng),1,vec,1)
          call grpop(vec,ag(1,ng),plat,1)
C         call dgemm('N','N',3,1,3,1d0,plat,3,vec,3,0d0,ag(1,ng),3)
        endif
      endif
      goto 90
      end
      subroutine parsop(t,i,a)
C- Parse string for a point group operator
      double precision v(3),sp,c,s,pi2,a(3,3),ddot
      character*1 t(0:*)
      logical parsvc
      integer i,j,k,nrot,iii
      pi2 = 8*datan(1d0)
      if (t(i) .eq. 'r' .or. t(i) .eq. 'R') then
        i = i+1
        read(t(i),'(i1)',err=99) nrot
        i = i+1
        if (.not. parsvc(-1,t,i,v)) goto 99
        sp = ddot(3,v,1,v,1)
        sp = 1d0/dsqrt(sp)
        do  14  k = 1, 3
          v(k) = v(k)*sp
   14   continue
        c = dcos(pi2/nrot)
        s = dsin(pi2/nrot)
        do  16  k = 1, 3
          do  15  j = 1, 3
            a(k,j) = (1-c)*v(j)*v(k)
   15     continue
          a(k,k) = a(k,k) + c
   16   continue
        a(2,1) = a(2,1) + s*v(3)
        a(1,3) = a(1,3) + s*v(2)
        a(3,2) = a(3,2) + s*v(1)
        a(1,2) = a(1,2) - s*v(3)
        a(3,1) = a(3,1) - s*v(2)
        a(2,3) = a(2,3) - s*v(1)
      else if (t(i) .eq. 'm' .or. t(i) .eq. 'M') then
        i = i+1
        if (.not. parsvc(-1,t,i,v)) goto 99
        sp = ddot(3,v,1,v,1)
        do  11  j = 1, 3
          do  12  k = 1, 3
            a(j,k) = -2.d0*v(k)*v(j)/sp
   12     continue
          a(j,j) = a(j,j) + 1d0
   11   continue
      else if (t(i) .eq. 'i' .or. t(i) .eq. 'I') then
        i = i+1
C       call dvcpy(0d0,0,a,1,9)
C       call dvcpy(-1d0,0,a,4,3)
        call dpzero(a,9)
        a(1,1) = -1
        a(2,2) = -1
        a(3,3) = -1
      else if (t(i) .eq. 'e' .or. t(i) .eq. 'E') then
        i = i+1
C       call dvcpy(0d0,0,a,1,9)
C       call dvcpy(-1d0,0,a,4,3)
        call dpzero(a,9)
        a(1,1) = 1
        a(2,2) = 1
        a(3,3) = 1
      else
        goto 99
      endif
      return
   99 print *, 'PARSOP: parse error at ',(t(iii),iii = 0,i),'  ...'
      call fexit(-1,119,' ',0)
      end

C      subroutine fmain
C      implicit none
C      double precision g(9,9),ag(3,9)
C      integer ng,ig
C      character *50 t
C
C      t = 'R3D(.1,.2,.3) Mx R2(1/2,pi,1)(pi,0,-1) my*i'
C      call psymop(t,0d0,g,ag,ng)
C
C      print *, ng
C      do  10  ig = 1, ng
C        call awrit1('%9:1d',' ',80,6,g(1,ig))
C        call awrit1('%3:1d',' ',80,6,ag(1,ig))
C   10 continue
C
CC should produce:
CC            4
CC 0 1 0 0 0 1 1 0 0
CC 0 0 1 1 0 0 0 1 0
CC 0.1 0.2 0.3
CC -1 0 0 0 1 0 0 0 1
CC 0 0 0
CC -0.955034 0.282527 0.089931 0.282527 0.775172 0.565055 0.089931 ...
CC 3.141593 0 -1
CC -1 0 0 0 1 0 0 0 -1
CC 0 0 0
C      end

