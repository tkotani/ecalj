      subroutine paug1(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,
     .  lmxa,nlml,cg,jcg,indxcg,vum,nlx1,nlx2,ppi)
C- Add to ppi constribution from true pot, true wave functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of functions of first kind for each l
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there is a smooth part to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   v1    :values of f1 at rofi(nr)
Ci   d1    :slopes of f1 at rofi(nr)
Ci   nf2   :number of functions of second kind for each l
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :a smooth part is to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   v1    :values of f1 at rofi(nr)
Ci   d1    :slopes of f1 at rofi(nr)
Ci   lmxa  :augmentation l-cutoff
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
Ci         :decomposed by M where v1 = sum_M v1_M Y_M
Ci         :vum(l1,l2,M,1) = (u_l1 v1_M u_l2)
Ci         :vum(l1,l2,M,2) = (u_l1 v1_M s_l2)
Ci         :vum(l1,l2,M,3) = (s_l1 v1_M s_l2)
Ci         :vum(l1,l2,M,4) = (u_l1 v1_M g_l2)
Ci         :vum(l1,l2,M,5) = (s_l1 v1_M g_l2)
Ci         :vum(l1,l2,M,6) = (g_l1 v1_M g_l2)
Ci         :Note that
Ci         :vum(l2,l1,M,2) = (s_l1 v1_M u_l2)
Ci         :vum(l2,l1,M,4) = (g_l1 v1_M u_l2)
Ci         :vum(l2,l1,M,5) = (g_l1 v1_M s_l2)
Ci   nlx1  :dimensions ppi
Ci   nlx2  :dimensions ppi
Co Outputs
Co   ppi   :partial matrix element of potential; see Remarks
Cr Remarks
Cr    Makes the first half of the first term in Eq. 29, Springer book.
Cr    But see Remarks in augmat.f:  there are three flavors of this
Cr    contribution to pi:
Cr         P~ V1 P~      H~ V1 P~      H~ V1 H~
Cr    where V1 is true potential.
Cu Updates
Cu   27 Aug 01 Extended to local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmx1,lmx2,lmxa,nf1,nf1s,nf2,nf2s,nlml,nlx1,nlx2
      integer lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
      double precision vum(0:lmxa,0:lmxa,nlml,6),
     .  v1(0:lmx1,nf1),d1(0:lmx1,nf1),
     .  v2(0:lmx2,nf2),d2(0:lmx2,nf2),
     .  ppi(nf1,nf2,nlx1,nlx2),cg(1)
C ... Local parameters
      integer i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2
      double precision add

C ... Combine with CG coefficents
      do  i1 = 1, nf1s
        do  i2 = 1, nf2s
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2
          do  ilm1 = 1, nlm1
            l1 = ll(ilm1)
            do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              ix = max0(ilm1,ilm2)
              ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
              do  icg = indxcg(ix),indxcg(ix+1)-1
                mlm = jcg(icg)
                if (mlm.gt.1 .and. mlm.le.nlml) then

                  add = v1(l1,i1) * v2(l2,i2) * vum(l1,l2,mlm,1)
     .                + v1(l1,i1) * d2(l2,i2) * vum(l1,l2,mlm,2)
     .                + d1(l1,i1) * v2(l2,i2) * vum(l2,l1,mlm,2)
     .                + d1(l1,i1) * d2(l2,i2) * vum(l1,l2,mlm,3)
                  ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2)
     .                                 + cg(icg)*add

                endif
              enddo
            enddo
          enddo
        enddo
      enddo


C --- Matrix elements involving local orbitals ---
      if (nf1s .ge. nf1 .and. nf2s .ge. nf2) return
      do  i1 = 1, nf1
        do  i2 = 1, nf2
          if (i1 .gt. nf1s .or. i2 .gt. nf2s) then
            nlm1 = (lx1(i1)+1)**2
            nlm2 = (lx2(i2)+1)**2
            do  ilm1 = 1, nlm1
              l1 = ll(ilm1)
              do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                ix = max0(ilm1,ilm2)
                ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                do  icg = indxcg(ix),indxcg(ix+1)-1
                  mlm = jcg(icg)
                  if (mlm.gt.1 .and. mlm.le.nlml) then
C                 <g | V | g>
                    if (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                      add = vum(l1,l2,mlm,6)
C                 <g | V | (u,s)>
                    elseif (i1 .gt. nf1s) then
C                   Bug fix 27 Jun 05
                      add = vum(l2,l1,mlm,4) * v2(l2,i2)
     .                  + vum(l2,l1,mlm,5) * d2(l2,i2)
C                   add = v2(l2,i2) * vum(l1,l2,mlm,4)
C     .                 + d2(l2,i2) * vum(l1,l2,mlm,5)
C                 <(u,s) | V | g>
                    elseif (i2 .gt. nf2s) then
C                   Bug fix 27 Jun 05
                      add = v1(l1,i1) * vum(l1,l2,mlm,4)
     .                  + d1(l1,i1) * vum(l1,l2,mlm,5)
C                   add = v1(l1,i1) * vum(l1,l2,mlm,4)
C     .                 + d1(l1,i1) * vum(l2,l1,mlm,5)
                    endif
                    ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2)
     .                                 + cg(icg)*add
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      end

      subroutine paug2(nr,nlml,v2,rwgt,cg,jcg,indxcg,
     .  nf1,nf1s,lmx1,lx1,f1,nf2,nf2s,lmx2,lx2,f2,sum,nlx1,nlx2,ppi)
C- Put in ppi constribution from smooth pot, smooth wave functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   v2    :smooth potential, seen by unaugmented functions
Ci   rwgt  :radial mesh weights
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   nf1   :number of functions of first kind for each l
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there is a smooth part to be subtracted
Ci   lmx1  :dimensions f1,sum
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   f1    :`bra' radial functions tabulated numerically; see Outputs
Ci   nf2   :number of functions of second kind for each l
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :a smooth part is to be subtracted.
Ci   lmx1  :dimensions f2,sum
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   f2    :`ket' radial functions tabulated numerically; see Outputs
Ci   sum   :work array holding integrations
Ci   nlx1  :dimensions ppi
Ci   nlx2  :dimensions ppi
Co Outputs
Co   ppi   :<f1^ | V2~ | f2^> subtracted for each f1,f2 pair
Cr Remarks
Cr    Makes the 2nd half of the first term in Eq. 29, Springer book.
Cr    But see Remarks in augmat.f:  there are three flavors of this
Cr    contribution to pi:
Cr         P V2~ P      H V2~ P      H V2~ H
Cr    where V2~ is the one-center repsn'f of the smooth potential.
Cu Updates
Cu   24 Aug 01 Extended to local orbitals, which have no smooth part
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmx1,lmx2,nf1,nf2,nf1s,nf2s,nlml,nlx1,nlx2,nr
      integer lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
      double precision v2(nr,nlml),rwgt(nr),
     .  ppi(nf1,nf2,nlx1,nlx2),f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s),
     .  cg(1),sum(0:lmx1,0:lmx2,nlml)
C ... Local parameters
      integer i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2
      double precision sam

C ... Zero out ppi
      call dpzero(ppi, nf1*nf2*nlx1*nlx2)

C ... Sum over CG coefficients, make radial integrals as needed
      do  i1 = 1, nf1s
        do  i2 = 1, nf2s
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2

C         Set flag for all integrals --- integral not yet calculated
          do  l1 = 0, lx1(i1)
            do  l2 = 0, lx2(i2)
              do  mlm = 1, nlml
                sum(l1,l2,mlm) = 2d10
              enddo
            enddo
          enddo

          do  ilm1 = 1, nlm1
            l1 = ll(ilm1)
            do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              ix = max0(ilm1,ilm2)
              ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
              do icg = indxcg(ix),indxcg(ix+1)-1
                mlm = jcg(icg)
                if (mlm.gt.1 .and. mlm.le.nlml) then
                  if (sum(l1,l2,mlm) .gt. 1d10) then
                    call paug4(lmx1,f1,lmx2,f2,i1,i2,l1,l2,mlm,
     .                         nr,rwgt,v2,sam)
                    sum(l1,l2,mlm) = sam
C                   if (sam .ne. 0) write(*,500) i1,i2,ilm1,ilm2,sam
C 500               format(2i4,2x,2i4,f14.8)
                  endif
                  ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2)
     .                                 - cg(icg)*sum(l1,l2,mlm)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

      end

      subroutine paug3(nf1,lmx1,lx1,nf2,lmx2,lx2,lmxl,nlml,
     .  cg,jcg,indxcg,qm,gpotb,gpot0,lmux,ppi0,nlx1,nlx2,ppi)
C- Assemble the final potential augmentation matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of functions of first kind for each l
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   nf2   :number of functions of second kind for each l
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   lmxl  :l-cutoff for density, potential on the radial mesh
Ci   nlml  :(lmxl+1)*(lmxl+1)
Ci   cg    :Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   qm    :Moments of the products f1~*f2~ - f1*f2
Ci         :For local orbitals, the term f1*f2 is absent
Ci   gpotb :integrals of compensating gaussians * local smooth estat
Ci         :pot calculated from the compensated smoothed local density
Ci   gpot0 :integrals of local gaussians * phi0~ (smves.f)
Ci         :phi0~ is the estatic potential of the interstitial
Ci         :smooth density including compensating local charges.
Ci   lmux  :l-cutoff for sigma,tau, and spherical part of ppi
Ci   ppi0  :contribution to ppi from spherical part of potential
Ci   nlx1  :dimensions ppi
Ci   nlx2  :dimensions ppi
Ci   ppi   :nonspherical parts only of potential integrals (paug1,paug2)
Co Outputs
Co   ppi   :augmentation potential integrals assembled
Cr Remarks
Cr  The terms proportional to the spherical part of V are added to
Cr  ppi (generated in pvagm1, pvagm1c, pvagm2), and also the term
Cr     qm (gpot0-gpotb)
Cr  which are the terms proportional to Qkk'LL'M of Eqs. 28,29 in
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr   However, see comments in augmat.f about indices k,k'.
Cu Updates
Cu   14 Sep 01 Extended to local orbitals, which have no smooth part
Cu   17 May 00 Adapted from nfp paug1.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmux,lmx1,lmx2,lmxl,nf1,nf2,nlml,nlx1,nlx2
      integer lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
      double precision ppi0(nf1,nf2,0:lmux),ppi(nf1,nf2,nlx1,nlx2),
     .  cg(1),gpotb(nlml),gpot0(nlml),qm(nf1,nf2,0:lmx1,0:lmx2,0:lmxl)
C ... Local parameters
      integer i1,i2,nlm1,nlm2,ilm1,ilm2,l1,l2,ll,ix,icg,mlm,lm

C     print *, 'paug3 zero'
C     ppi = 0
C     ppi0 = 0
C     gpot0 = 0
C     gpotb = 0
C     gpot0(1) = 1

C ... Add terms from moments of f~g~ - fg and ppi0 to ppi
      do  i1 = 1, nf1
        do  i2 = 1, nf2
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2
          do  ilm1 = 1, nlm1
            l1 = ll(ilm1)
            if (ilm1 .le. nlm2)
     .      ppi(i1,i2,ilm1,ilm1) = ppi(i1,i2,ilm1,ilm1) + ppi0(i1,i2,l1)
C           if (i1 .le. nf1s .and. i2 .le. nf2s) then
            do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              ix = max0(ilm1,ilm2)
              ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
              do  icg = indxcg(ix),indxcg(ix+1)-1
                mlm = jcg(icg)
                if (mlm .le. nlml) then
                  lm = ll(mlm)
                  ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) +
     .               cg(icg)*qm(i1,i2,l1,l2,lm)*(gpot0(mlm)-gpotb(mlm))
                endif
              enddo
            enddo
C           endif
          enddo
        enddo
      enddo

C      print *, nf1,nf2
C      write(66) ppi
      end

      subroutine paug4(lmx1,f1,lmx2,f2,i1,i2,l1,l2,mlm,nr,rwgt,v2,sam)
C- Make an integral f1*V_mlm*f2
C     implicit none
C ... Passed parameters
      integer i1,i2,l1,l2,lmx1,lmx2,mlm,nr
      double precision sam,v2(nr,1),rwgt(nr),
     .  f1(nr,0:lmx1,1),f2(nr,0:lmx2,1)
C ... Local parameters
      integer i

      sam = 0d0
      do  i = 2, nr
        sam = sam + rwgt(i)*v2(i,mlm)*f1(i,l1,i1)*f2(i,l2,i2)
      enddo

      end
      subroutine paugnl(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,
     .                  lmaxu,vumm,nlx1,nlx2,ppiz,isp,idu)
C- Add to ppi contribution from non local part of pot (LDA+U)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of functions of first kind for each l
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there is a smooth part to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   v1    :values of f1 at rofi(nr)
Ci   d1    :slopes of f1 at rofi(nr)
Ci   nf2   :number of functions of second kind for each l
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :a smooth part is to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   v1    :values of f1 at rofi(nr)
Ci   d1    :slopes of f1 at rofi(nr)
Ci   lmaxu : used to dimension vumm
Ci   vumm  :matrix elements of non-local potential
Ci         :vumm(m1,m2,1) = <u| vorb(m1,m2) |u>
Ci         :vumm(m1,m2,2) = <u| vorb(m1,m2) |s>
Ci         :vumm(m1,m2,3) = <s| vorb(m1,m2) |u>
Ci         :vumm(m1,m2,4) = <s| vorb(m1,m2) |s>
Ci         :vumm(m1,m2,5) = <u| vorb(m1,m2) |z>
Ci         :vumm(m1,m2,6) = <s| vorb(m1,m2) |z>
Ci         :vumm(m1,m2,7) = <z| vorb(m1,m2) |z>
Ci         :vumm(m1,m2,8) = <z| vorb(m1,m2) |u>
Ci         :vumm(m1,m2,9) = <z| vorb(m1,m2) |s>
Ci   nlx1  :dimensions ppiz
Ci   nlx2  :dimensions ppiz
Ci   isp   : spin we're working on
Ci   idu   : idu(l)=1 => this l has a U
Co Outputs
Co   ppiz  :partial matrix element of potential; see Remarks
Cr Remarks
Cu Updates
Cu   09 Nov 05 Convert ppi to complex form
Cu   08 Jun 05 (MvS) extended to local orbitals
Cu   27 Apr 05 (Lambrecht) LDA+U
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmx1,lmx2,nf1,nf1s,nf2,nf2s,nlx1,nlx2,lmaxu
      integer lx1(nf1),lx2(nf2),isp,idu(4),nab
      parameter (nab=9)
      double precision
     .  v1(0:lmx1,nf1),d1(0:lmx1,nf1),
     .  v2(0:lmx2,nf2),d2(0:lmx2,nf2)
      double complex ppiz(nf1,nf2,nlx1,nlx2)
      double complex vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
C ... Local parameters
      integer i1,i2,ilm1,ilm2,l1,l2,m1,m2
      double complex add

C     call zprm('paugnl, init ppiz',2,ppiz,nf1*nf2,nf1*nf2,nlx1*nlx2)
C     ilm1 = 2*lmaxu+1
C     call zprm('paugnl, vumm',2,vumm,ilm1,ilm1,ilm1)

      do  i1 = 1, nf1
        do  i2 = 1, nf2

C     ... Matrix elements of vumm constructed from (u,s)
          ilm1 = 0
          do  l1 = 0, min(lx1(i1),lmaxu)
            do  m1 = -l1, l1
              ilm1 = ilm1+1
              ilm2 = 0
              do  l2 = 0, min(lx2(i2),lmaxu)
                do  m2 = -l2, l2
                  ilm2 = ilm2+1
                  if (idu(l1+1) .ne. 0 .and. l1 .eq. l2) then

C           ... (u,s)V(u,s)
                    if (i1 .le. nf1s .and. i2 .le. nf2s) then
                      add = v1(l1,i1) * v2(l2,i2) * vumm(m1,m2,1,isp,l1)
     .                + v1(l1,i1) * d2(l2,i2) * vumm(m1,m2,2,isp,l1)
     .                + d1(l1,i1) * v2(l2,i2) * vumm(m1,m2,3,isp,l1)
     .                + d1(l1,i1) * d2(l2,i2) * vumm(m1,m2,4,isp,l1)

C           ... zVz
                    elseif (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                      add = vumm(m1,m2,7,isp,l1)
C           ... zV(u,s)
                    elseif (i1 .gt. nf1s) then
                      add = vumm(m1,m2,8,isp,l1) * v2(l2,i2)
     .                + vumm(m1,m2,9,isp,l1) * d2(l2,i2)
C           ... (u,s)Vz
                    elseif (i2 .gt. nf2s) then
                      add = v1(l1,i1) * vumm(m1,m2,5,isp,l1)
     .                + d1(l1,i1) * vumm(m1,m2,6,isp,l1)
                    endif

                    ppiz(i1,i2,ilm1,ilm2) = ppiz(i1,i2,ilm1,ilm2) + add
                  endif
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo

      end

