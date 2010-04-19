      subroutine moment(mode,nl,nlo,nrhos,isp,nsp,nspc,nkp,ldim,nevx,
     .  nev,ikp,wgt,nclass,doswt,dosw2,qnu,rhos,orbtm,swgt)
C- Calculate moments from band and k-point weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :See outputs.
Ci         :0 accumulate qnu from doswt(0..2)
Ci         :1 accumulate qnu from doswt(0..2), rhos from dosw2
Ci            and orbtm from doswt(3)
Ci         :2 accumulate orbtm from doswt(1)
Ci         :100s digit: accumlate orbtm
Ci   nl    :(global maximum l) + 1
Ci   nlo   :number of channels for which to accumulate ASA moments.
Ci          nl0=nl when moments contracted by m, nl**2 when distinct.
Ci   nrhos :number of channels which spin density-matrix is stored
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nevx  :dimensions doswt, dosw2
Ci   nev   :actual number of eigenvectors generated
Ci   ikp   :k-point label
Ci   wgt   :weights for accumulating moments (written to disc in makwts)
Ci   nclass:number of inequivalent classes
Ci   doswt :weights for moments
Ci   dosw2 :weights for off-diagonal part of spin-density matrix
Co Outputs:
Co  data for this k-point is added into the following:
Co   qnu:  Charge and energy moments (1s digit mode nonzero)
Co   rhos: Spin density matrix (10s digit mode nonzero)
Co   orbtm:orbital moments (100s digit mode nonzero)
Co   swgt  :sum of the weight factors added to swgt
Cr Remarks
Cr   Accumulates moments as defined in CECAM notes, 6/9/87, eqn 2.67
Cu Updates
Cu   27 May 08 allow orbtm to be accumulated inpependently
Cu   21 Apr 04 rhos,dosw2 can be m-dependent.  New argument list
Cu   19 Sep 01 Returns swgt.  Altered argument list.
Cu   08 Dec 00 Orbital moments computed in spin-coupled case
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nlo,nrhos,isp,nsp,nspc,nkp,ldim,nclass,nevx,nev,ikp,
     .  mode
      double precision qnu(0:2,nl,nsp,nclass),
     .  doswt(nlo,nclass,nevx,nspc,0:3),
     .  rhos(2,0:2,nrhos,2,2,nclass),orbtm(nlo,2,nclass),
     .  dosw2(2,nrhos,nclass,nevx,0:2),wgt(ldim,nsp,nkp)
C Local variables
      logical lallm,lrhol,lrhos,lqnu,lorbtm
      integer ib,ic,l,m,lm,k,jsp,ksp,i,ll,k2,ilm
      double precision wt,swgt,fac

C --- Setup: read moments, spin density matrix from disk ---
      lallm = nlo .eq. nl**2
      lrhol = nrhos .eq. nl
      if (nsp .lt. nspc) call rx('MOMENT not set up for nsp lt nspc')
      if (isp .eq. 2 .and. nspc .eq. 2)
     .  call rx('MOMENT called with isp eq 2 and nspc eq 2')
      call isanrg(mode,0,2,'moment','mode',.true.)
      lqnu = mode .lt. 2
      lrhos = mode .eq. 1
      lorbtm = mode .gt. 0 .and. nspc .eq. 2

C --- Loop over energy levels, accumulating with appropriate weights ---
      do  5  ib = 1, nev
      do  5  jsp = 1, nspc
C ...   ksp is isp in collinear case, and jsp in noncollinear
        ksp = max(jsp,isp)
        wt = wgt(ib,isp,ikp)
        swgt = swgt + wt

C   ... Accumulate moments of the density, and spin density matrix
        do  4  ic = 1, nclass
          do  3  lm = 1, nlo
            l = lm
            if (lallm) l = ll(lm)+1
            k2 = l
            do  2  i = 0, 2
              if (lqnu) then
                qnu(i,l,ksp,ic) = qnu(i,l,ksp,ic) +
     .            doswt(lm,ic,ib,jsp,i)*wt
              endif
              if (lrhos .and. lrhol) then
                rhos(1,i,k2,ksp,ksp,ic)   = qnu(i,l,ksp,ic)
                rhos(1,i,k2,ksp,3-ksp,ic) = rhos(1,i,k2,ksp,3-ksp,ic) +
     .                                      dosw2(1,k2,ic,ib,i)*wt
                rhos(2,i,k2,ksp,3-ksp,ic) = rhos(2,i,k2,ksp,3-ksp,ic) +
     .                                dosw2(2,k2,ic,ib,i)*wt*(3-2*ksp)

              endif
    2     continue
          if (lorbtm) then
            k = 3
            if (mode .eq. 2) k = 1
            orbtm(lm,jsp,ic) = orbtm(lm,jsp,ic) +
     .                         doswt(lm,ic,ib,jsp,k)*wt
C            print 333,lm,ib,ic,wt,doswt(lm,ic,ib,jsp,k)
C  333       format(3i4,2f15.10)
          endif
    3     continue
          if (lrhos .and. .not. lrhol) then
            ilm = 0
            do  13  l = 0, nl-1
              do  15  m = -l, l
                ilm = ilm+1
                if (lallm) then
                  k = ilm
                  fac = 1
                else
                  k = l+1
                  fac = 2*l+1
                endif
                k2 = ilm
                do  12  i = 0, 2
                rhos(1,i,k2,ksp,ksp,ic)   = rhos(1,i,k2,ksp,ksp,ic) +
     .                                    doswt(k,ic,ib,jsp,i)*wt/fac
                rhos(1,i,k2,ksp,3-ksp,ic) = rhos(1,i,k2,ksp,3-ksp,ic) +
     .                                      dosw2(1,k2,ic,ib,i)*wt
                rhos(2,i,k2,ksp,3-ksp,ic) = rhos(2,i,k2,ksp,3-ksp,ic) +
     .                                dosw2(2,k2,ic,ib,i)*wt*(3-2*ksp)
   12           continue
   15         continue
   13       continue
          endif
    4   continue
    5 continue

      end

