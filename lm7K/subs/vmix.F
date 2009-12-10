      subroutine vmix(nclass,nl,nsp,nbas,ves,vesold,pp,
     .  iclass,nrclas,vtil,votil,uv,evmad,modcst,igroup,
     .  rmax,mad,nvmix,lmixpq)
C- Mixing for constrained potential variation
C ----------------------------------------------------------------
Ci Inputs
Ci   nclass,nl,nsp
Ci   a:    Workspace of four vectors
Ci   beta: mixing parameter: mixes beta * new  +  (1-beta) * old
Ci         (see function amix, below)
Ci   mmix: number of previous iterations for mixing.
Ci   nvmix: number of elements to mix
Ci   ves,vesold: input for last iteration
Ci   vtil, votil, uv: modcst:  mix potentials according to modcst
Ci   wmix: weight potentials by 1/evmad for mixing
Co Outputs
Co   ves is mixed
Co   rmsdel is the rms change in parameters
Cr Remarks
Cr   sign of beta is passed as sign of npmix<0 in function amix
Cr   Improvments: suggestion of PYB, take lc of moms 0 and 2.
Cr   In the case modst=-1, vtil must be rotated by M^-1 (see rotmad).
Cr   Also vtil(2) is constrained to be votil(2) (before rotation)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lmixpq
      integer nbas,nclass,nl,nsp,nmix,mmix,nvmix,nrclas(0:1),iclass(1),
     .  modcst,igroup(0:7)
      double precision pp(0:*),rmax(1),mad(1),ves(0:3),vesold(0:3),
     .  vtil(2),votil(2),uv(nbas,1),evmad(1)
C Local variables:
      integer na,ic,i,j,k,n1,imix,nelts,onorm,okpvt,
     .        fopn,iprint,amix,jmix,lgunit,ido,ib
      integer nlspi6,nlspic
      logical lpr
      external fopn,iprint
      double precision ddot,dsqrt,newrms,newbet,oldrms,oldbet,d1mach,xx
      integer ngrp,now,nvmix0

C heap:
      integer w(1)
      common /w/ w

      print *, 'vmix: not implemented, revert to normal mixing' 
      lmixpq = .true.
      return

C --- OLD ---
CC --- Copy P's and Q's into mixing array ---
C      do  110  ib = 1, nvmix
C        a(ib-1,0,1) = vtil(ib)
C        a(ib-1,0,2) = votil(ib)
C  110 continue
C      na = nvmix
C      if (na .gt. nvmix) stop 'VMIX: too many elements'
C
CC --- Prepare for evaluation of new beta ----
C      oldrms = rmsdel(1)
C      oldbet = rmsdel(2)
C
CC --- Calculate this rmsdel and estimate new beta  ---
C      newrms =  dsqrt(dabs(ddot(nvmix,a,1,a,1)
C     .       -2*ddot(nvmix,a,1,a(0,0,2),1)
C     .        + ddot(nvmix,a(0,0,2),1,a(0,0,2),1))/nvmix)
C      newbet = beta
C
C      if (iprint() .ge. 30) then
C      do  15  j = 1, 2
C   15 write(lgunit(j),334) oldrms, oldbet, newrms, newbet
C  334   format(/
C     .    ' VMIX: Old rms del   Old beta    New rms del   New beta'/
C     .    4x,4f13.6)
C        call query('beta=',4,beta)
C      endif
C      rmsdel(1) = newrms
C      rmsdel(2) = beta
C
CC --- Do the mixing ---
C      imix = 0
C      if (mmix .gt. 0) then
C        n1 = fopn('MIXV')
C        read(n1,103,err=8,end=8) imix, nelts
C  103   format(2i5)
C        if (nelts .eq. nvmix .and. imix .le. mmix) goto 12
C    8   imix = 0
C        goto 40
CC When reached this point, have found previous iter w/ correct # elts
C   12   continue
C        read(n1,104) (((a(i,j,k), i=0,nvmix-1), j=1, imix), k=1,2)
C  104   format(1p,4d18.11)
C      endif
CC Now have read in imix previous iterations, if there were any
C   40 continue
C
C      if (wmix) then
C        do  50  ib = 2, nbas
C   50   call dscal((mmix+2)*2,1/(.5d0+evmad(ib)),a(ib-1,0,1),nvmix)
C      endif
C
C      call defdr(onorm,mmix**2)
C      call defi(okpvt,mmix)
C      jmix = imix
C      if (beta .lt. 0) jmix = -imix
C      nmix = amix(nvmix,jmix,mmix,0,dabs(beta),iprint()+1,tjmax,
C     .            w(onorm),w(okpvt),a,tj,rms2)
C
C      if (wmix) then
C        do  51  ib = 2, nbas
C   51   call dscal((mmix+2)*2,(.5d0+evmad(ib)),a(ib-1,0,1),nvmix)
C      endif
C
CC --- Save this iteration into mixing file ---
C      imix = min(imix+1,mmix)
C      if (mmix .gt. 0) then
C        rewind n1
C        write(n1,103) imix, nvmix
C        write(n1,104) (((a(i,j,k), i=0,nvmix-1), j=1, imix), k=1,2)
C        call fclose(n1)
C      endif
C
CC --- Backtransform vtil into ves ---
C      call dcopy(nvmix,a(0,0,2),1,vtil,1)
C      if (modcst .eq. -1) then
C        vtil(2) = votil(2)
C        xx =     -evmad(1)*vtil(1) + evmad(2)*vtil(2)
C        vtil(1) = evmad(2)*vtil(1) + evmad(1)*vtil(2)
C        vtil(2) = xx
C        xx =      -evmad(1)*votil(1) + evmad(2)*votil(2)
C        votil(1) = evmad(2)*votil(1) + evmad(1)*votil(2)
C        votil(2) = xx
C      endif
C      if (iabs(modcst) .eq. 1) then
C        call ivshel(1,nclass,igroup,igroup(nclass),.true.)
C
CC --- Disperse change in V back into classes ---
C        nvmix0 = 0
C        now = igroup(igroup(nclass))-1
C        ngrp = 0
C        do  60  ic = 0, nclass
C          if (ic.eq.nclass .or. now .ne. igroup(igroup(nclass+ic))) then
C            do  62  i = ic-ngrp, ic-1
C              j = igroup(nclass+i)
C              ves(j) = vesold(j) + (vtil(nvmix0)-votil(nvmix0))
C   62       continue
C            now = igroup(igroup(nclass+ic))
C            ngrp = 0
C            nvmix0 = nvmix0+1
C          endif
C          ngrp = ngrp+1
C   60   continue
C        if (modcst .ne. -1 .and. nvmix .ne. nvmix0-1)
C     .    stop 'mismatch in nvmix'
C      else if (iabs(modcst) .eq. 2) then
C        call dmpy(uv,nbas,1,a(0,0,2),nbas,1,a,nbas,1,nbas,1,nbas)
C        call dpzero(ves,nclass)
C        do  25  ib = 1, nbas
C          ic = iclass(ib)
C          ves(ic) = ves(ic) + a(ib-1,0,1)/nrclas(ic)
C   25   continue
C      endif
C
CC --- Copy back to ves, Shift enu and c by ves ---
C      lpr = iprint() .ge. 50 .or.
C     .      iprint() .ge. 40 .and. modcst .ne. 0
C      if (lpr)
C     .  write(*,'(/'' VMIX:  ic       V          Vold        Del'')')
C      na = 0
C      do  11  ic = 0, nclass-1
C        if (lpr) print 440, ic,ves(ic),vesold(ic),ves(ic)-vesold(ic)
C  440     format(i9, 3f12.5)
C        nlspic = nl*nsp*ic
C        nlspi6 = 6*nlspic
C        call daxpy(nl*nsp,-1d0,vesold(ic),0,pp(nlspi6),6)
C        call daxpy(nl*nsp,-1d0,vesold(ic),0,pp(nlspi6+1),6)
C        call daxpy(nl*nsp,1d0,ves(ic),0,pp(nlspi6),6)
C        call daxpy(nl*nsp,1d0,ves(ic),0,pp(nlspi6+1),6)
C   11 continue
C      rms2 =  dsqrt(dabs(ddot(nclass,ves,1,ves,1)
C     .       -2*ddot(nclass,ves,1,vesold,1)
C     .        + ddot(nclass,vesold,1,vesold,1))/nclass)/beta
C      print 441, rms2
C  441 format(' rms delta v/beta',16x,f12.5)
C      if (iprint() .ge. 40) call query('abort?',-1,0)
C
      end

