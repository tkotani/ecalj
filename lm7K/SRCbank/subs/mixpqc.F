      subroutine mixpqc(nclass,nrclas,nl,nsp,igroup,pnu,qnu)
C- Averages together classes of equivalent groups
C ----------------------------------------------------------------
Ci Inputs
Ci   nclass,nrclas,nl,nsp,igroup
Ci   jgroup,pwk,qwk: work arrays
Co Outputs
Co   pnu,qnu are averaged by group
Cr Remarks
Cl Variables
Cl   cgrp is current group; grps is size of current group;
Cl   ngrp is number of groups; ic is class
C ----------------------------------------------------------------
C     implicit none
      integer nclass,nrclas(0:1),nl,nsp,igroup(0:nclass)
      double precision pnu(nl,nsp,0:nclass-1),
     .                 qnu(3,nl,nsp,0:nclass-1)
C Local variables:
      integer ncmx
      parameter (ncmx=1024)
      integer ngrp,ic,i,grps,grpw,cgrp,iprint,jgroup(0:2*ncmx),isp,ispc,
     .  i0
      double precision pwk(10,2), qwk(3,10,2)

C --- Copy iabs(igroup) to jgroup ---
      if (nclass .gt. ncmx) call rx('mixpqc: increase ncmx')
      call icopy(nclass,igroup,1,jgroup,1)
      do  70  ic = 0, nclass-1
   70 jgroup(ic) = iabs(jgroup(ic))

C --- Sort groups and find first with group>0 ---
      call ivshel(1,nclass,jgroup,jgroup(nclass),.true.)
      do  80  i = 0, nclass-1
        i0 = i
        if (jgroup(jgroup(nclass+i)) .ne. 0) goto 82
   80 continue
   82 if (i0 .eq. nclass-1) return

C --- Accumulate P, Q by group ---
      if (iprint() .ge. 30) print '(/'' Mixpqc: ic  ngrp  group'')'
      ngrp = 0
      grps = 0
      grpw = 0
      cgrp = 0
      call dpzero(pwk,nl*nsp)
      call dpzero(qwk,nl*nsp*3)
C --- Loop through classes in order of permutation table ---
      do  40  i = i0, nclass-1
C ---   New group if group ID changes ... ---
        if (cgrp .ne. jgroup(jgroup(nclass+i))) then
          cgrp = jgroup(jgroup(nclass+i))
          if (grps .ne. 0)
     .    call mixpqx(i,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .      jgroup,pwk,qwk,pnu,qnu)
          grps = 0
          grpw = 0
          ngrp = ngrp+1
          call dpzero(pwk,20)
          call dpzero(qwk,60)
        endif
        grps = grps + 1
        ic = jgroup(nclass+i)
        grpw = grpw + nrclas(ic)
C ---   Sum P,Q of this class into its group ---
        do  50  isp = 1, nsp
          ispc = isp
          if (igroup(ic) .lt. 0) ispc = nsp+1-isp
          call daxpy(nl,dble(nrclas(ic)),pnu(1,isp,ic),1,pwk(1,ispc),1)
          call daxpy(3*nl,dble(nrclas(ic)),qnu(1,1,isp,ic),1,
     .      qwk(1,1,ispc),1)
   50   continue
   40 continue
      call mixpqx(nclass,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .  jgroup,pwk,qwk,pnu,qnu)

      end
      subroutine mixpqx(i,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .  jgroup,pwk,qwk,pnu,qnu)
C- Kernel called by mixpqc to redistribute group avg into each class
C     implicit none
      integer i,ngrp,grpw,grps,nl,nsp,igroup(0:3),
     .  jgroup(0:3),cgrp,nclass
      double precision pnu(nl,nsp,0:1),qnu(3,nl,nsp,0:1),
     .  pwk(10,2),qwk(3,10,2),q
      integer i1,iprint,ic,isp,ispc
      double precision dsum


      call dscal(10*nsp,1/dble(grpw),pwk,1)
      call dscal(30*nsp,1/dble(grpw),qwk,1)
      q = dsum(nl,qwk,3)
      if (nsp .eq. 2) q = q + dsum(nl,qwk(1,1,2),3)

C --- Redistribute P,Q of this group into each class ---
      do  10  i1 = i-grps, i-1
        ic = jgroup(nclass+i1)
        if (iprint() .ge. 40) print 333, ic+1,cgrp,ngrp
  333   format(i11,i5,i7)
        do  20  isp = 1, nsp
          ispc = isp
          if (igroup(ic) .lt. 0) ispc = nsp+1-isp
          call dcopy(nl,pwk(1,ispc),1,pnu(1,isp,ic),1)
          call dcopy(3*nl,qwk(1,1,ispc),1,qnu(1,1,isp,ic),1)
   20   continue
   10 continue
      if (iprint().ge.30) print 334, grps, ngrp, q
  334 format(i6,' elts in group',i3,'  Q=',f9.6)

      end

