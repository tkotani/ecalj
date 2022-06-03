subroutine paug1(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2, &
     lmxa,nlml,cg,jcg,indxcg,vum,nlx1,nlx2,ppi)
  !- Add to ppi constribution from true pot, true wave functions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nf1   :number of functions of first kind for each l
  !i   nf1s  :number of functions of first kind for each l, for which
  !i         :there is a smooth part to be subtracted (which also
  !i         :corresponds to the functions which connect to envelope
  !i         :functions)
  !i   lmx1  :dimensions f1
  !i   lx1   :l-cutoffs for each of the nf1 functions
  !i   v1    :values of f1 at rofi(nr)
  !i   d1    :slopes of f1 at rofi(nr)
  !i   nf2   :number of functions of second kind for each l
  !i   nf2s  :number of functions of second kind for each l, for which
  !i         :a smooth part is to be subtracted (which also
  !i         :corresponds to the functions which connect to envelope
  !i         :functions)
  !i   lmx2  :dimensions f2
  !i   lx2   :l-cutoffs for each of the nf2 functions
  !i   v1    :values of f1 at rofi(nr)
  !i   d1    :slopes of f1 at rofi(nr)
  !i   lmxa  :augmentation l-cutoff
  !i   nlml  :L-cutoff for charge density on radial mesh
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
  !i         :decomposed by M where v1 = sum_M v1_M Y_M
  !i         :vum(l1,l2,M,1) = (u_l1 v1_M u_l2)
  !i         :vum(l1,l2,M,2) = (u_l1 v1_M s_l2)
  !i         :vum(l1,l2,M,3) = (s_l1 v1_M s_l2)
  !i         :vum(l1,l2,M,4) = (u_l1 v1_M g_l2)
  !i         :vum(l1,l2,M,5) = (s_l1 v1_M g_l2)
  !i         :vum(l1,l2,M,6) = (g_l1 v1_M g_l2)
  !i         :Note that
  !i         :vum(l2,l1,M,2) = (s_l1 v1_M u_l2)
  !i         :vum(l2,l1,M,4) = (g_l1 v1_M u_l2)
  !i         :vum(l2,l1,M,5) = (g_l1 v1_M s_l2)
  !i   nlx1  :dimensions ppi
  !i   nlx2  :dimensions ppi
  !o Outputs
  !o   ppi   :partial matrix element of potential; see Remarks
  !r Remarks
  !r    Makes the first half of the first term in Eq. 29, Springer book.
  !r    But see Remarks in augmat.f:  there are three flavors of this
  !r    contribution to pi:
  !r         P~ V1 P~      H~ V1 P~      H~ V1 H~
  !r    where V1 is true potential.
  !u Updates
  !u   27 Aug 01 Extended to local orbitals.  Altered argument list.
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmx1,lmx2,lmxa,nf1,nf1s,nf2,nf2s,nlml,nlx1,nlx2
  integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
  double precision :: vum(0:lmxa,0:lmxa,nlml,6), &
       v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
       v2(0:lmx2,nf2),d2(0:lmx2,nf2), &
       ppi(nf1,nf2,nlx1,nlx2),cg(1)
  ! ... Local parameters
  integer :: i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2
  double precision :: add=1d99

  ! ... Combine with CG coefficents
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
                 if (mlm > 1 .AND. mlm <= nlml) then

                    add = v1(l1,i1) * v2(l2,i2) * vum(l1,l2,mlm,1) &
                         + v1(l1,i1) * d2(l2,i2) * vum(l1,l2,mlm,2) &
                         + d1(l1,i1) * v2(l2,i2) * vum(l2,l1,mlm,2) &
                         + d1(l1,i1) * d2(l2,i2) * vum(l1,l2,mlm,3)
                    ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) &
                         + cg(icg)*add

                 endif
              enddo
           enddo
        enddo
     enddo
  enddo


  ! --- Matrix elements involving local orbitals ---
  if (nf1s >= nf1 .AND. nf2s >= nf2) return
  do  i1 = 1, nf1
     do  i2 = 1, nf2
        if (i1 > nf1s .OR. i2 > nf2s) then
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
                    if (mlm > 1 .AND. mlm <= nlml) then
                       !                 <g | V | g>
                       if (i1 > nf1s .AND. i2 > nf2s) then
                          add = vum(l1,l2,mlm,6)
                          !                 <g | V | (u,s)>
                       elseif (i1 > nf1s) then
                          !                   Bug fix 27 Jun 05
                          add = vum(l2,l1,mlm,4) * v2(l2,i2) &
                               + vum(l2,l1,mlm,5) * d2(l2,i2)
                          !                   add = v2(l2,i2) * vum(l1,l2,mlm,4)
                          !     .                 + d2(l2,i2) * vum(l1,l2,mlm,5)
                          !                 <(u,s) | V | g>
                       elseif (i2 > nf2s) then
                          !                   Bug fix 27 Jun 05
                          add = v1(l1,i1) * vum(l1,l2,mlm,4) &
                               + d1(l1,i1) * vum(l1,l2,mlm,5)
                          !                   add = v1(l1,i1) * vum(l1,l2,mlm,4)
                          !     .                 + d1(l1,i1) * vum(l2,l1,mlm,5)
                       endif
                       ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) &
                            + cg(icg)*add
                    endif
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo

end subroutine paug1

subroutine paug2(nr,nlml,v2,rwgt,cg,jcg,indxcg, &
     nf1,nf1s,lmx1,lx1,f1,nf2,nf2s,lmx2,lx2,f2,sum,nlx1,nlx2,ppi)
  !- Put in ppi constribution from smooth pot, smooth wave functions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nr    :number of radial mesh points
  !i   nlml  :L-cutoff for charge density on radial mesh
  !i   v2    :smooth potential, seen by unaugmented functions
  !i   rwgt  :radial mesh weights
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   nf1   :number of functions of first kind for each l
  !i   nf1s  :number of functions of first kind for each l, for which
  !i         :there is a smooth part to be subtracted
  !i   lmx1  :dimensions f1,sum
  !i   lx1   :l-cutoffs for each of the nf1 functions
  !i   f1    :`bra' radial functions tabulated numerically; see Outputs
  !i   nf2   :number of functions of second kind for each l
  !i   nf2s  :number of functions of second kind for each l, for which
  !i         :a smooth part is to be subtracted.
  !i   lmx1  :dimensions f2,sum
  !i   lx2   :l-cutoffs for each of the nf2 functions
  !i   f2    :`ket' radial functions tabulated numerically; see Outputs
  !i   sum   :work array holding integrations
  !i   nlx1  :dimensions ppi
  !i   nlx2  :dimensions ppi
  !o Outputs
  !o   ppi   :<f1^ | V2~ | f2^> subtracted for each f1,f2 pair
  !r Remarks
  !r    Makes the 2nd half of the first term in Eq. 29, Springer book.
  !r    But see Remarks in augmat.f:  there are three flavors of this
  !r    contribution to pi:
  !r         P V2~ P      H V2~ P      H V2~ H
  !r    where V2~ is the one-center repsn'f of the smooth potential.
  !u Updates
  !u   24 Aug 01 Extended to local orbitals, which have no smooth part
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmx1,lmx2,nf1,nf2,nf1s,nf2s,nlml,nlx1,nlx2,nr
  integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
  double precision :: v2(nr,nlml),rwgt(nr), &
       ppi(nf1,nf2,nlx1,nlx2),f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s), &
       cg(1),sum(0:lmx1,0:lmx2,nlml)
  ! ... Local parameters
  integer :: i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2
  double precision :: sam

  ! ... Zero out ppi
  call dpzero(ppi, nf1*nf2*nlx1*nlx2)

  ! ... Sum over CG coefficients, make radial integrals as needed
  do  i1 = 1, nf1s
     do  i2 = 1, nf2s
        nlm1 = (lx1(i1)+1)**2
        nlm2 = (lx2(i2)+1)**2

        !         Set flag for all integrals --- integral not yet calculated
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
                 if (mlm > 1 .AND. mlm <= nlml) then
                    if (sum(l1,l2,mlm) > 1d10) then
                       call paug4(lmx1,f1,lmx2,f2,i1,i2,l1,l2,mlm, &
                            nr,rwgt,v2,sam)
                       sum(l1,l2,mlm) = sam
                       !                   if (sam .ne. 0) write(*,500) i1,i2,ilm1,ilm2,sam
                       ! 500               format(2i4,2x,2i4,f14.8)
                    endif
                    ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) &
                         - cg(icg)*sum(l1,l2,mlm)
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo

end subroutine paug2

subroutine paug3(nf1,lmx1,lx1,nf2,lmx2,lx2,lmxl,nlml, &
     cg,jcg,indxcg,qm,gpotb,gpot0,lmux,ppi0,nlx1,nlx2,ppi)
  !- Assemble the final potential augmentation matrix
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nf1   :number of functions of first kind for each l
  !i   lmx1  :dimensions f1
  !i   lx1   :l-cutoffs for each of the nf1 functions
  !i   nf2   :number of functions of second kind for each l
  !i   lmx2  :dimensions f2
  !i   lx2   :l-cutoffs for each of the nf2 functions
  !i   lmxl  :l-cutoff for density, potential on the radial mesh
  !i   nlml  :(lmxl+1)*(lmxl+1)
  !i   cg    :Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   qm    :Moments of the products f1~*f2~ - f1*f2
  !i         :For local orbitals, the term f1*f2 is absent
  !i   gpotb :integrals of compensating gaussians * local smooth estat
  !i         :pot calculated from the compensated smoothed local density
  !i   gpot0 :integrals of local gaussians * phi0~ (smves.f)
  !i         :phi0~ is the estatic potential of the interstitial
  !i         :smooth density including compensating local charges.
  !i   lmux  :l-cutoff for sigma,tau, and spherical part of ppi
  !i   ppi0  :contribution to ppi from spherical part of potential
  !i   nlx1  :dimensions ppi
  !i   nlx2  :dimensions ppi
  !i   ppi   :nonspherical parts only of potential integrals (paug1,paug2)
  !o Outputs
  !o   ppi   :augmentation potential integrals assembled
  !r Remarks
  !r  The terms proportional to the spherical part of V are added to
  !r  ppi (generated in pvagm1, pvagm1c, pvagm2), and also the term
  !r     qm (gpot0-gpotb)
  !r  which are the terms proportional to Qkk'LL'M of Eqs. 28,29 in
  !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
  !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
  !r      ed. (Springer-Verlag, Berlin) 2000.
  !r   However, see comments in augmat.f about indices k,k'.
  !u Updates
  !u   14 Sep 01 Extended to local orbitals, which have no smooth part
  !u   17 May 00 Adapted from nfp paug1.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmux,lmx1,lmx2,lmxl,nf1,nf2,nlml,nlx1,nlx2
  integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
  double precision :: ppi0(nf1,nf2,0:lmux),ppi(nf1,nf2,nlx1,nlx2), &
       cg(1),gpotb(nlml),gpot0(nlml),qm(nf1,nf2,0:lmx1,0:lmx2,0:lmxl)
  ! ... Local parameters
  integer :: i1,i2,nlm1,nlm2,ilm1,ilm2,l1,l2,ll,ix,icg,mlm,lm

  !     print *, 'paug3 zero'
  !     ppi = 0
  !     ppi0 = 0
  !     gpot0 = 0
  !     gpotb = 0
  !     gpot0(1) = 1

  ! ... Add terms from moments of f~g~ - fg and ppi0 to ppi
  do  i1 = 1, nf1
     do  i2 = 1, nf2
        nlm1 = (lx1(i1)+1)**2
        nlm2 = (lx2(i2)+1)**2
        do  ilm1 = 1, nlm1
           l1 = ll(ilm1)
           if (ilm1 <= nlm2) &
                ppi(i1,i2,ilm1,ilm1) = ppi(i1,i2,ilm1,ilm1) + ppi0(i1,i2,l1)
           !           if (i1 .le. nf1s .and. i2 .le. nf2s) then
           do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              ix = max0(ilm1,ilm2)
              ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
              do  icg = indxcg(ix),indxcg(ix+1)-1
                 mlm = jcg(icg)
                 if (mlm <= nlml) then
                    lm = ll(mlm)
                    ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) + &
                         cg(icg)*qm(i1,i2,l1,l2,lm)*(gpot0(mlm)-gpotb(mlm))
                 endif
              enddo
           enddo
           !           endif
        enddo
     enddo
  enddo

  !      print *, nf1,nf2
  !      write(66) ppi
end subroutine paug3

subroutine paug4(lmx1,f1,lmx2,f2,i1,i2,l1,l2,mlm,nr,rwgt,v2,sam)
  !- Make an integral f1*V_mlm*f2
  !     implicit none
  ! ... Passed parameters
  integer :: i1,i2,l1,l2,lmx1,lmx2,mlm,nr
  double precision :: sam,v2(nr,1),rwgt(nr), &
       f1(nr,0:lmx1,1),f2(nr,0:lmx2,1)
  ! ... Local parameters
  integer :: i

  sam = 0d0
  do  i = 2, nr
     sam = sam + rwgt(i)*v2(i,mlm)*f1(i,l1,i1)*f2(i,l2,i2)
  enddo

end subroutine paug4
subroutine paugnl(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2, &
     lmaxu,vumm,nlx1,nlx2,ppiz,isp,idu)
  !- Add to ppi contribution from non local part of pot (LDA+U)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nf1   :number of functions of first kind for each l
  !i   nf1s  :number of functions of first kind for each l, for which
  !i         :there is a smooth part to be subtracted (which also
  !i         :corresponds to the functions which connect to envelope
  !i         :functions)
  !i   lmx1  :dimensions f1
  !i   lx1   :l-cutoffs for each of the nf1 functions
  !i   v1    :values of f1 at rofi(nr)
  !i   d1    :slopes of f1 at rofi(nr)
  !i   nf2   :number of functions of second kind for each l
  !i   nf2s  :number of functions of second kind for each l, for which
  !i         :a smooth part is to be subtracted (which also
  !i         :corresponds to the functions which connect to envelope
  !i         :functions)
  !i   lmx2  :dimensions f2
  !i   lx2   :l-cutoffs for each of the nf2 functions
  !i   v1    :values of f1 at rofi(nr)
  !i   d1    :slopes of f1 at rofi(nr)
  !i   lmaxu : used to dimension vumm
  !i   vumm  :matrix elements of non-local potential
  !i         :vumm(m1,m2,1) = <u| vorb(m1,m2) |u>
  !i         :vumm(m1,m2,2) = <u| vorb(m1,m2) |s>
  !i         :vumm(m1,m2,3) = <s| vorb(m1,m2) |u>
  !i         :vumm(m1,m2,4) = <s| vorb(m1,m2) |s>
  !i         :vumm(m1,m2,5) = <u| vorb(m1,m2) |z>
  !i         :vumm(m1,m2,6) = <s| vorb(m1,m2) |z>
  !i         :vumm(m1,m2,7) = <z| vorb(m1,m2) |z>
  !i         :vumm(m1,m2,8) = <z| vorb(m1,m2) |u>
  !i         :vumm(m1,m2,9) = <z| vorb(m1,m2) |s>
  !i   nlx1  :dimensions ppiz
  !i   nlx2  :dimensions ppiz
  !i   isp   : spin we're working on
  !i   idu   : idu(l)=1 => this l has a U
  !o Outputs
  !o   ppiz  :partial matrix element of potential; see Remarks
  !r Remarks
  !u Updates
  !u   09 Nov 05 Convert ppi to complex form
  !u   08 Jun 05 (MvS) extended to local orbitals
  !u   27 Apr 05 (Lambrecht) LDA+U
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmx1,lmx2,nf1,nf1s,nf2,nf2s,nlx1,nlx2,lmaxu
  integer :: lx1(nf1),lx2(nf2),isp,idu(4),nab
  parameter (nab=9)
  double precision :: &
       v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
       v2(0:lmx2,nf2),d2(0:lmx2,nf2)
  double complex ppiz(nf1,nf2,nlx1,nlx2)
  double complex vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
  ! ... Local parameters
  integer :: i1,i2,ilm1,ilm2,l1,l2,m1,m2
  double complex add

  !     call zprm('paugnl, init ppiz',2,ppiz,nf1*nf2,nf1*nf2,nlx1*nlx2)
  !     ilm1 = 2*lmaxu+1
  !     call zprm('paugnl, vumm',2,vumm,ilm1,ilm1,ilm1)

  do  i1 = 1, nf1
     do  i2 = 1, nf2

        !     ... Matrix elements of vumm constructed from (u,s)
        ilm1 = 0
        do  l1 = 0, min(lx1(i1),lmaxu)
           do  m1 = -l1, l1
              ilm1 = ilm1+1
              ilm2 = 0
              do  l2 = 0, min(lx2(i2),lmaxu)
                 do  m2 = -l2, l2
                    ilm2 = ilm2+1
                    if (idu(l1+1) /= 0 .AND. l1 == l2) then

                       !           ... (u,s)V(u,s)
                       if (i1 <= nf1s .AND. i2 <= nf2s) then
                          add = v1(l1,i1) * v2(l2,i2) * vumm(m1,m2,1,isp,l1) &
                               + v1(l1,i1) * d2(l2,i2) * vumm(m1,m2,2,isp,l1) &
                               + d1(l1,i1) * v2(l2,i2) * vumm(m1,m2,3,isp,l1) &
                               + d1(l1,i1) * d2(l2,i2) * vumm(m1,m2,4,isp,l1)

                          !           ... zVz
                       elseif (i1 > nf1s .AND. i2 > nf2s) then
                          add = vumm(m1,m2,7,isp,l1)
                          !           ... zV(u,s)
                       elseif (i1 > nf1s) then
                          add = vumm(m1,m2,8,isp,l1) * v2(l2,i2) &
                               + vumm(m1,m2,9,isp,l1) * d2(l2,i2)
                          !           ... (u,s)Vz
                       elseif (i2 > nf2s) then
                          add = v1(l1,i1) * vumm(m1,m2,5,isp,l1) &
                               + d1(l1,i1) * vumm(m1,m2,6,isp,l1)
                       endif

                       ppiz(i1,i2,ilm1,ilm2) = ppiz(i1,i2,ilm1,ilm2) + add
                    endif
                 enddo
              enddo
           enddo
        enddo

     enddo
  enddo

end subroutine paugnl

