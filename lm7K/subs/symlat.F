      subroutine symlat(platcp,ngrp,grp,isym)
C- Generates the (point) symmetry operations of the lattice
C ----------------------------------------------------------------------
Ci Inputs:
Ci   platcp:lattice vectors of most compact primitive unit cell
Co Outputs:
Co   ngrp  :number of allowed symmetry operations
Co   grp   :symmetry operation matrix
Co   isym  :index to lattice type, calculated from ngrp:
Co          ngrp   isym    name
Co                  0     shouldn't happen
Co            2     1     triclinic
Co            4     2     monoclinic
Co            8     3     orthorhombic
Co           16     4     tetragonal
Co           12     5     rhombohedral
Co           24     6     hexagonal
Co           48     7     cubic
Cr Remarks:
Cr   symlat analyzes the primitive translations of the bravais
Cr   lattice in order to supply the symmetry operations of the lattice.
Cr   It gives the number ngrp of allowed operations as well as
Cr   these operations themselves.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer ngrp,isym
      double precision platcp(3,3),grp(9,*)
C Local parameters:
      integer i,iprint,ltmax,ll1,lgunit,m,m1,m2,m3,mm,nrot(4)
      parameter(ltmax=3,ll1=ltmax*2+1)
      double precision platt(9),qlatcp(3,3),mat(9),vecg(3),vol
      logical latvec,lirr
      character*12 csym1(0:7)

C External calls:
      external dcopy,dinv33,dmpy,errmsg,iprint,latvec,lgunit
      data nrot /2,3,4,6/
      data csym1 /'indefinite','triclinic','monoclinic','orthorhombic',
     .          'tetragonal','rhombohedral','hexagonal','cubic'/

      mm(i,m) = ltmax-(mod(i,ll1**m)-mod(i,ll1**(m-1)))/ll1**(m-1)

      call dinv33(platcp,1,qlatcp,vol)
C --- Start out with E and I ---
      ngrp = 2
      call csymop(-1,grp(1,1),.false.,1,0d0)
      call csymop(-1,grp(1,2),.true.,1,0d0)

C --- Find all possible rotation axes ---
      do  10  i = 0, (ll1**3-1)/2-1
        m1 = mm(i,1)
        m2 = mm(i,2)
        m3 = mm(i,3)
        lirr = .true.
        do  12  m = 2, ll1
          lirr = lirr.and.(mod(m1,m).ne.0.or.mod(m2,m).ne.0.or.
     .                   mod(m3,m).ne.0)
   12   continue
        if (lirr) then
          do  14  m = 1, 3
            vecg(m) = m1*platcp(m,1) + m2*platcp(m,2) + m3*platcp(m,3)
   14     continue

          do  16  m = 1, 4
C       ... Matrix for this symmetry operation
            call csymop(-1,mat,.false.,nrot(m),vecg)
            call grpprd(mat,platcp,platt)
C           call dmpy(mat,3,1,platcp,3,1,platt,3,1,3,3,3)
C       ... Add it and i*symop, if allowed
            if (latvec(3,1d-5,qlatcp,platt)) then
              call csymop(-1,grp(1,ngrp+1),.false.,nrot(m),vecg)
              call csymop(-1,grp(1,ngrp+2),.true. ,nrot(m),vecg)
              ngrp = ngrp+2
              if (m .ne. 1) then
                call csymop(-1,grp(1,ngrp+1),.false.,-nrot(m),vecg)
                call csymop(-1,grp(1,ngrp+2),.true. ,-nrot(m),vecg)
                ngrp = ngrp+2
              endif
            endif
   16     continue
        endif
   10 continue

      isym = 0
      if (ngrp .eq. 2) isym=1
      if (ngrp .eq. 4) isym=2
      if (ngrp .eq. 8) isym=3
      if (ngrp .eq.16) isym=4
      if (ngrp .eq.12) isym=5
      if (ngrp .eq.24) isym=6
      if (ngrp .eq.48) isym=7

      if (iprint() .ge. 30) call awrit1(' SYMLAT: Bravais system is '
     .  //csym1(isym)//'%a with %i symmetry operations.',' ',80,
     .  lgunit(1),ngrp)

      end

