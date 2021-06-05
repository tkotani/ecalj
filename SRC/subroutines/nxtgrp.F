      subroutine nxtgrp(nlst,match,grpid,iprm,i1,i2)
C- Return range in list of next group of equivalent elements
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlst : size of list
Ci   grpid: a list of group identification numbers.  Groups are
Ci          equivalent when they have the same identification number
Ci   iprm : permutation table that sorts grpid
Ci          NB iprm(1) <= 0 => assumes iprm(i) = i; iprm not referenced
Ci   i1   : first element to check
Co Outputs
Co   i1   : first element for which grpid(iprm(i)) = match
Co   i2   : last contiguous element for which grpid(iprm(i)) = match
Co   mxsiz: largest of largest group
Cr Remarks
Cr   Finds all elements in grpid(k) equivalent to 'match', for
Cr   k = iprm(i) and i= i1,i1+1,...nlst, stopping when i=nlst or when
Cr   grpid(i) differs from 'match'.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlst,match,grpid(nlst),iprm(nlst),i1,i2
      integer i,j,ii,jj
      logical liprm

      liprm = iprm(1) .gt. 0
      i2 = i1-1
      do  10  i = i1, nlst
        if (liprm) then
          ii = iprm(i)
        else
          ii = i
        endif
        if (grpid(ii) .ne. match) goto 10
C   ... Found a match ; now find last contiguous member
        i1 = i
        do  20  j = i1+1, nlst+1
          i2 = j-1
          if (i2 .eq. nlst) return
          if (liprm) then
            jj = iprm(j)
          else
            jj = j
          endif
          if (grpid(jj) .ne. match) return
   20   continue
   10 continue
      end
      subroutine ngroup(nlst,grpid,iprm,ngrp,mnsiz,mxsiz)
C- Sort a list of group identification numbers, count number of groups
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlst  :size of list
Ci   grpid :a list of group identification numbers.  Groups are
Ci          equivalent when they have the same identification number.
Co Outputs
Co   iprm  :permutation table that sorts grpid
Co   ngrp  :number of groups
Co   mnsiz :size of smallest group
Co   mxsiz :largest of largest group
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlst,grpid(nlst),iprm(nlst),ngrp,mxsiz,mnsiz
C ... Local parameters
      integer i,now,ib,gsize

      call ivshel(1,nlst,grpid,iprm,.true.)
      do   i = 1, nlst
         iprm(i) = iprm(i)+1
      enddo
      ngrp = 0
      gsize = 0
      mxsiz = 0
      mnsiz = nlst
      now = grpid(iprm(1))
      do  20  ib = 1, nlst+1
        if (ib .eq. nlst+1 .or. now.ne.grpid(iprm(min(ib,nlst)))) then
          ngrp = ngrp+1
C     ... Sort elements within this group
          call ishell(gsize,iprm(ib-gsize))
          now = grpid(iprm(min(ib,nlst)))
          mxsiz = max(mxsiz,gsize)
          mnsiz = min(mnsiz,gsize)
          gsize = 0
        endif
        gsize = gsize+1
   20 continue

      end

