      subroutine sympad(nlsqr,nbas,npadl,npadr,nsite0,nsites,iax,nkap,
     .  sflg,s)
C- Symmetrize padded part (bulk) structure constants (pgf)
C ----------------------------------------------------------------
Co Outputs: a portion of s is symmetrized (see Remarks).
Cr Remarks
Cr Symmetrize s(r1,l1,r2,l2) = s(r2,l2,r1,l1) by averaging the two.
Cr   This routine symmetrizes only the bulk (padded) part connected
Cr   to the next padded layer.  Symmetrization is accomplished by
Cr   finding the equivalent strux connecting a PL to the padded layer,
Cr   since strux originating in the second padded layer are missing.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nlsqr,nsite0,nsites,nkap,nbas,npadl,npadr,niax
      parameter (niax=10)
      integer iax(niax,nsites),sflg(nsites)
      double precision s(nlsqr,nlsqr,nkap,nkap,nsites)
C Local parameters
      integer i,j,lm1,lm2,ikap,jkap,nbasp,iat1,iat2,imap,ii,jj
      double precision temp

C Cludge for now ...
      if (nkap .ne. 1) return

C --- For each iat1>nbas, iat2>nbas+npad, find the corresponding pair --
      nbasp = nbas + npadl + npadr
      do  10  i = nsite0+1, nsites
        iat1 = iax(1,i)
        iat2 = iax(2,i)
        if (iat2 .gt. nbasp+npadl) then
C ...   Use Rpad -> Rppad  equivalent to R -> Rpad
          iat1 = iat1 - npadl - npadr
          iat2 = iat2 - npadl - npadr
        elseif (iat2 .gt. nbasp) then
C ...   Use Lpad -> Lppad  equivalent to L -> Lpad
          iat1 = iat1 - nbas
          iat2 = iat2 - npadl - npadr
        else
C ...   Not one of (iat1>nbas, iat2>nbas+npad)
          goto 10
        endif
C  ---  For this pair connecting pad->ppad, find the equivalent PL->pad
        do  30  j = 1, nsite0
          if (iax(1,j) .eq. iat1 .and. iax(2,j) .eq. iat2 .and. 
     .        iax(3,i) .eq. iax(3,j) .and.
     .        iax(4,i) .eq. iax(4,j) .and.
     .        iax(5,i) .eq. iax(5,j))  then
            imap = j
            goto 32
          endif
   30   continue
        call fexit3(-1,111,' Exit -1 sympad: missing pair%2:1i '//
     .    '(mapped to %i %i)',iax(1,i),iat1,iat2)
   32   continue
C       call awrit5(' sympad site %i->%i: (iax=%2:1i mapped to %i %i)',
C    .    ' ',255,lgunit(1),i,imap,iax(1,i),iat1,iat2)

        do  40  j = nsite0+1, nsites

          if (iax(1,j) .eq. iat2 .and. iax(2,j) .eq. iat1 .and.
     .        iax(3,imap) .eq. -iax(3,j) .and.
     .        iax(4,imap) .eq. -iax(4,j) .and.
     .        iax(5,imap) .eq. -iax(5,j))  then

            ii = iax(8,i)
            if (ii .eq. 0) ii = i
            jj = iax(8,j)
            if (jj .eq. 0) jj = j

C       ... Symmetrize ii as (ii+jj)/2 if neither was symmetrized
            if (sflg(ii) .eq. 0 .and. sflg(jj) .eq. 0) then
              do  42  ikap = 1, nkap
              do  42  jkap = 1, nkap
              do  42  lm1 = 1, nlsqr
              do  42  lm2 = 1, nlsqr
              temp =(s(lm1,lm2,ikap,jkap,ii)+s(lm2,lm1,jkap,ikap,jj))/2
              s(lm1,lm2,ikap,jkap,ii) = temp
C             s(lm2,lm1,jkap,ikap,jj) = temp
   42         continue
C         ... Flag this s(isite) as being symmetrized
              sflg(ii) = 1
C       ... Use jj if it has been symmetrized
            elseif (sflg(ii) .eq. 0) then
              do  44  ikap = 1, nkap
              do  44  jkap = 1, nkap
              do  44  lm1 = 1, nlsqr
              do  44  lm2 = 1, nlsqr
   44         s(lm2,lm1,jkap,ikap,ii) = s(lm2,lm1,jkap,ikap,jj)
C         ... Flag this s(isite) as being symmetrized
              sflg(ii) = 1
            endif
            goto 10
          endif
   40   continue
        call fexit3(-1,111,'SYMPAD: no pair corresponding to site'//
     .    ' %i (mapped to pair %i, iat =%2:1i)',i,imap,iax(1,imap))
   10 continue
      end

