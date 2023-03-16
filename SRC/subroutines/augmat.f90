    !o Outputs
    !o   osig  :augmentation overlap integrals; see Remarks.
    !o   otau  :augmentation kinetic energy integrals; see Remarks.
    !o   oppi  :augmentation kinetic + potential integrals; see Remarks.
    !o   phzdphz  :phz dphz
    !o   hab   :matrix elements of the ham. with true w.f.  See Remarks.
    !o   vab   :matrix elements of the pot. with true w.f.  See Remarks.
    !o   sab   :matrix elements of    unity with true w.f.  See Remarks.
    !o   hsozz,hsopm !spin orbit integrals
    !i Inputs
    !i   lso   :if nonzero, calculate radial spin-orbit integrals
    !i   z     :nuclear charge
    !i   rmt   :augmentation sphere radius
    !i   rsma  :augmentation smoothing radius
    !i   lmxa  :augmentation L-cutoff
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pnz   :boundary conditions for second p.q.n (local orbital).
    !i          10s digit controls how local orbital included in hamiltonian
    !i   kmax  :polynomial cutoff
    !i   nlml  :L-cutoff for density, potential on the radial mesh
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rwgt  :radial mesh weights
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   v0    :spherical potential that defines phi and phidot
    !i   v1    :true nonspherical potential, seen by augmented functions, excluding
    !i         :nuclear contribution 2*Z/r
    !i   v2    :smooth nonspherical potential, seen by unaugmented functions
    !i         :estat part from n2~ = n2 + compensating gaussians + pseudocore charge
    !i   gpotb :integrals of local gaussians times local smooth ves;
    !i         :see Remarks
    !i   gpot0 :integrals of local gaussians * phi0~
    !i         :phi0~ is the estatic potential of the interstitial;
    !i         :see Remarks
    !i   nkaph :number of orbital types for a given L quantum no. in basis
    !i         :dimensions sig,tau,ppi
    !i   nkapi :number of valence envelope function types per l q.n.
    !i   lmxh  :largest l in basis; must be >= max(lh)
    !i   lh    :list of nkaph l-cutoffs in basis
    !i   eh    :energy of smoothed Hankel
    !i   rsmh  :smoothing radii of smoothed hankel for each l, energy.
    !i         :rsmh(l+1,ik) for must be zero for all ik=1..nkaph unless
    !i         :it channel is a valence or extended local orbital.
    !i     ... The following four parameters are used to extrapolate
    !i         quantities outside the MT radius.  Extrapolation is used to
    !i         correct matrix elements for local orbitals.
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    !i   rs3   :smoothing radius for extrapolation of MT potential
    !i   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
    !i         :The asymptotic form of V-vmtz is taken to be zero.
    !i   ...   The following are LDA+U-related
    !i   lmaxu :dimensioning parameter for U matrix
    !i   vorb  :orbital dependent potential matrices
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i          U on site ib with dmat beginning at dmats(*,lldau(ib))
    !i   idu   :l-dependent switch flagging which l's have U
    ! o Inputs/Outputs
    ! o  iblu  :index to current LDA+U block
    ! o        :On input, index to last LDA+U block that was accessed
    ! o        :iblu will be incremented to from blocks at this site
    !
    ! l  rofi  :radial mesh points.  rofi(1..nr) are made at first by radmsh.
    ! l        :if V is to be extrapolated outside its MT sphere, to
    ! l        :V(1..nrbig), rofi(nr+1,nrbig) are also generated
    ! l        :Thus MUST be dimensioned at least rofi(1..nrbig)
    ! l        :nrbig is internally generated, but will not
    ! l        :exceed parameter nrx defined vxtrap.
    !r Remarks
    !r   This subroutine implements the computation of matrices
    !r   sigma, tau, ppi that comprise the local (augmented) part of
    !r   the hamiltonian in the full-potential method described in
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r   See discussion in Eqns 20-29.
    !r
    !r  In this method, the integral of a product of two augmented functions
    !r  F~i and F~j, with corresponding (smooth) envelopes Fi and Fj, and
    !r  their one-center expansions F^i and F^j, is calculated as
    !r
    !r     int_all F~i F~j  =  int_all Fi Fj  + int_MT (F~i F~j - F^i F^j)
    !r
    !r  Augmented F~i matches continuously and differentiably to Fi.
    !r  The one-center expansion F^i of Fi is identical except the former
    !r  is truncated to a finite L-cutoff in its on-center expansion.
    !r
    !r     Fi  = sum_L=0..infinity Fi_L
    !r     F^i = sum_L=0..lmxa     Fi_L
    !r
    !r  In the following description, we assume one atom per cell and
    !r  one orbital per L channel to simplify notation: i->L and j->L'.
    !r  If FL is further approximated by a polynomial expansion inside
    !r  one augmentation sphere, viz
    !r
    !r     F~L = sum_kL CkL P~kL                                 (cf Eq 16)
    !r     F^L = sum_kL CkL PkL
    !r
    !r  the integral becomes
    !r
    !r     int_all F~L F~L' = int_all FL FL' + sum CkL sig_kk'l Ck'L' (20)
    !r
    !r  where
    !r
    !r     sig_kk'l = int_MT (P~kL P~k'L' - PkL Pk'L')              (21)
    !r
    !r  is independent of the shape of any particular augmentation
    !r  function --- thus according to (20) one can compute sig once and
    !r  for all for a set of (kL) pairs, and evaluate the augmentation `on
    !r  the fly' when the coefficients C are known.  There is a
    !r  corresponding expression for kinetic energy and a (slightly more
    !r  complicated) form for pot. matrix elements pi:
    !r
    !r     int F~L V F~L'  = int FL V0~ FL' + ppi_LL'
    !r
    !r  where
    !r
    !r     ppi_LL' = int (F~L V1 F~L' - F^L V2~ F^L')
    !r             + sum_M Q_LL'M * int [(V0~ - V2~) G_M]
    !r
    !r     Q_LL'M  = int ( F~L F~L' - F^L F^L') r^m Y_M
    !r
    !r  As before, we can further approximate Fi~ and Fi^ as in Eq. 16; then
    !r  each of these terms can be expressed as a linear combination of
    !r  matrix elements ppi_kk'LL' independent of the shape of any
    !r  particular augmentation function:
    !r
    !r     ppi_LL' = sum_kk' C*_kL C_k'L'  ppi_kk'LL'
    !r
    !r  where
    !r
    !r     ppi_kk'LL' = int (P~kL V1 P~k'L' - PkL V2~ Pk'L')
    !r                + sum_M Q_kk'LL'M * int [(V0~ - V2~) G_M]
    !r
    !r     Q_kk'LL'M  = int ( P~kL P~k'L' -  PkL Pk'L') r^m Y_M
    !r
    !r  In the Springer book the augmentation for all spheres is presumed
    !r  to proceed along the lines of (20).  For numerical reasons it has
    !r  been found necessary to separate out heads from tails and compute
    !r  integrals of involving heads directly.  This is because for the
    !r  heads the P_kL expansion must be expanded to very high k, so
    !r  high that the expansion sometimes has numerical problems with
    !r  convergence.  Only tail-tail matrix elements are computed according
    !r  to (20).  For 1-center (head-head) augmentation integrals, the
    !r  integral (Fi~ Fj~ - Fi Fj) is explicitly computed; head-tail
    !r  integrals are mixed explicit repsn in one index and polynomial
    !r  expansions in the other; the (original) poly-poly integrals are
    !r  kept and used for the 3-center expansions.  Consequently, three
    !r  kinds of integrals are kept.
    !r           sig1             sig2             sig3
    !r         P~P~ - PP        H~P~ - HP        H~H~ - HH
    !r  with corresponding integrals tau and ppi.
    !r
    !r  In this version, augmented functions are linear combinations
    !r  of radial wave functions u and s that match continuously and
    !r  differentiably onto the envelope functions at rmt.
    !r  u and s are linear combinations of and phi,phidot defined as:
    !r  u has val=1, slo=1 at rmax, s has val=0, slo=1 at rmax.
    !r
    !r  Local orbitals : local orbitals occur in one of these types:
    !r  1. val,slo=0 at rmt                         (10s digit pnz=0)
    !rxxx  2. val,slo matches sm hankel                (10s digit pnz=1)
    !rxxx  3. val,slo matches sm hankel, perturbative  (10s digit pnz=2)
    !r  In all cases tails are not expanded about other sites.  Thus local
    !r  orbitals involve only one- and two-center terms.
    !r
    !r  Consider the overlap matrix where Fi is a local orbital:
    !r    int_all F~i F~j  =  int_all Fi Fj  + int_MT (F~i F~j - F^i F^j)
    !r
    !r  Case 1: Fi = F~i = F^i.  Then no error occurs in the L truncation
    !r  of F~j because
    !r     int_MT Fi Fj - F^i F^j = int_MT Fi (Fj-F^j) = 0
    !r  This is because the second factor is zero for l<=lmxa, and the
    !r  first factor is zero for l>lmxb.  This is similarly true for the
    !r  kinetic energy matrix elements.
    !r
    !r  Cases 2 and 3: The true envelope function and its one-center
    !r  expansion are identical. Also we have Fi=F~i for r<rmt, and Fi=F^i
    !r  for r>rmt; and the value and slope of F~i matches F^i at rmt.
    !r  There is a also exact cancellation in the overlap
    !r     int Fi Fj - F^i F^j = int_r<rmt F~i (Fj-F^j) +
    !r                           int_r>rmt F^i (Fj-F^j)
    !r  In each term, the second factor is zero for l<=lmxa, and the first
    !r  factor is zero for l>lmxb.
    !r
    !r  For the potential matrix element there is an additional
    !r  complication and an approximation.  Let F~i be a local orbital
    !r  with angular momentum L and consider a one-center expansion of the
    !r  matrix element a one-center expansion of the matrix element
    !r
    !r     int F~i V F~j  = F~i sum_L'' V_L'' sum_L' F~L'
    !r
    !r  If F~j has a finite Y_L expansion at this site, then the integral
    !r  is again exact provided that L'' is expanded at least up to the
    !r  difference in L's between F~i F~j (the l-cutoff for V is lmxl),
    !r  because F~i as only a single L in its Ylm expansion; therefore
    !r  only the L projection of the product (sum_L'' V_L'' sum_L' F~L')
    !r  makes any contribution to the integral.
    !r
    !r  But when F~j has higher L components (as it does when it is
    !r  centered at another site) but is truncated in its L expansion,
    !r  there is an error because there is a missing projection onto L of
    !r  the following:
    !r
    !r    sum_L'' sum_L'>La  V_L'' F~L' ~ sum_L'>La V_(L-L') F~L'
    !r
    !r  Since both factors are small (V_(L-L') is small for L-L'>0, and
    !r  F~L' is supposed to be small for L'>La), and moreover the product
    !r  V_(L-L') F~L is largest near the MT radius, at which point the
    !r  local orbital's value and slope are both zero, the integral of a
    !r  local orbital with this product clearly is very rapidly convergent
    !r  in its L-cutoff.
    !r
    !r  The complication arises that, if no smoothed analog is to be
    !r  computed, it is no longer true that the local representation of
    !r  the true potential can be shifted by an arbitrary harmonic
    !r  function r^m Y_M, as in the usual 'three-fold' representation of
    !r  the potential.  Therefore, in considering these matrix elements,
    !r  the local expansion of V must be written as
    !r
    !r    V(r) = sum_L V_L(r) Y_L,   V_L(r) = V1_l(r) + VVAL_L/R^l r^l
    !r
    !r  where the VVAL_L is the L projection of the e.s. potential at the
    !r  MT radius R.
    !r
    !r  *Documentation of hab,vab,sab and sig,tau,ppi :
    !r   hab,vab,sab are generated in potpus and are matrix elements of the
    !r   triplet of wave functions (u,s,gz) in the spherical part of V(r)
    !r   u and s are linear combinations of and phi,phidot defined as:
    !r     u has val=1, slo=1 at rmax;   s has val=0, slo=1
    !r   There may additionally be local orbitals gz, specified by nonzero
    !r   pnz, of one of two types, as described in potpus.f
    !r   Orbitals of the first type are confined to the augmentation sphere;
    !r   they have no corresponding smooth part.  Orbitals of the second type
    !r   extend into the interstitial by attaching a Hankel tail.
    !r
    !RRRR important note
    !r  *Structure of sig,tau,ppi, and hab,vab,sab.
    !r   All these arrays contain matrix elements between combinations
    !r   of valence states and local orbitals and are linear combinations of
    !r      int (ul,sl,gz,r*P_k)_l {h,v,1} (ul,sl,gz,r*P_k)_l'
    !r
    !r   There is some duplication in the two sets of matrix elements,
    !r   which can be confusing.  The following points help to explain
    !r   their relationship.
    !r
    !r  *{h,v,s}ab are matrix elements of products of true wave functions
    !r   in the spherical V only, in the (ul,sl,gz) form discussed above.
    !r
    !r  *sig,tau,ppi are matrix elements of products of augmented wave
    !r   functions, with the smooth part subtracted off.  ppi contains
    !r   matrix elements of the full potential, not merely the spherical
    !r   part. {h,v,s}ab are used to help assemble sig,tau,ppi.
    !r
    !r  *With the exception of ppi, all matrix elements are diagonal in l,
    !r   and independent of m; require only one l index.
    !r   In the case of ppi, the full LL' matrix is required.
    !r   Contributions from the nonspherical part of V are calculated with
    !r   the large component of the radial wave functions only.
    !r
    !r  *Matrix elements with a local orbital of the first type do not have
    !r   a smooth contribution subtracted.
    !r
    !r  *sig,tau,ppi come in three flavors as discussed above, e.g.
    !r   H~H~-HH; H~P~-HP; P~P~-PP, there is a separate array for each:
    !r     sig1(Pk,Pk) is dimensioned (1+kmax,1+kmax)
    !r     sig2(H, Pk) is dimensioned (nkaph,1+kmax)
    !r     sig3(H, H)  is dimensioned (nkaph,nkaph)
    !r   Here nkapi is 1 or 2 depending on whether there are one
    !r   or two kinds of envelope functions in the basis.
    !r
    !r   Thus the sig matrix elements are dimensioned
    !r     sig{1,2,3}(nf1,nf2,0..lmax,1..nsp)
    !r   with nf{1,2} = nkaph or 1+kmax, depending on which sig
    !r   and lmax = lmxa for those involving P's and lmxh for those not.
    !r
    !r * Documentation of LDA+U
    !r   the orbital dependent potential is of the form Vmm' for each LDA+U
    !r   block (fixed l and site), which makes matrix elements
    !r   <phi_m|Vnonlocal|phi_m'>
    !r   To add to ppi, these need to be first rotated to (u,s) basis which is
    !r   done in potpusnl whose output is vumm array for this site (possibly
    !r   for different l's). They are passed on in each gaugm call

    ! --- Augmentation matrices for cases P*P, H*H, H*P ---
    ! NOTE for SO
    !     We calculate limited parts of SO interaction. H(head),P(tail)
    !
    !     ohsozz(3,ib)%sdiag(:,isp) = <H isp | Lz(diag)|H isp>
    !     ohsozz(2,ib)%sdiag(:,isp) = <H isp | Lz(diag)|P isp>
    !     ohsozz(1,ib)%sdiag(:,isp) = <P isp | Lz(diag)|P isp>
    !
    !     ohsopm is desined for (1,2) element of LdotS matrix.
    !     ohsopm(1,ib)%soffd(:,isp) = <P isp | L-|P ispo> ;ispo=2 for isp=1, ispo=1 for isp=2
    !                                  (ispo means opposite spin to isp), (up:isp=1,dn:isp=2)
    !     ohsopm(2,ib)%soffd(:,isp) = <H isp | L-|P ispo>, this means
    !        > ohsopm(2,ib)%soffd(:,1)= <H isp=1|L-|P isp=2>
    !        > ohsopm(2,ib)%soffd(:,2)= <H isp=2|L+|P isp=1> => dagger(soffd(:,2))=<P isp=1|L-|H isp=2>
    !     ohsopm(3,ib)%soffd(:,isp) = <H isp | L-(diag)|H ispo> ;ispo=2 for isp=1, ispo=1 for isp=2
    !
    !   For SO=1 with spin=001 axis (liqinPRB2019 001 case),
    !   we neither use ohsopm(1,ib)%soffd(:,isp=2) and ohsopm(3,ib)%soffd(:,isp=2).
    !
    !     Generally speaking, we may need to calculate all blocks ohsozz%soffd, ohsopm%sdiag, ohsopp
module m_augmat !- Make augmentation matrices sig,tau,pi for one site
  use m_ll,only: ll
  use m_lmfinit,only: n0
  public vlm2us,momusl
  private
contains
  subroutine momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp,rofi,rwgt,v0,v1, qum,vum)!Moments of ul*ul,ul*sl,sl*sl and their integrals with true pot.
    !i   z     :nuclear charge
    !i   rmt   :augmentation radius, in a.u.
    !i   lmxa  :augmentation l-cutoff
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   lmxl  :l-cutoff for density, potential on the radial mesh
    !i   nlml  :(lmxl+1)*(lmxl+1)
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nsp   :number of spin channels
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights: int f(r) dr = sum_ir f(ir) * rwgt(ir)
    !i   v0    :spherical potential defining wave functions, without 2Z/r
    !i   v1    :true potential.  Spherical part need not be v0.
    !o Outputs
    !o   qum   :moments (u,s,gz)_l1 * (u,s,gz)_l2 * r**l
    !o         :qum(l1,l2,l,i,j) = (ui_l1 uj_l2) * r**l
    !o   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
    !o         :decomposed by M where v1 = sum_M v1_M Y_M
    !o         :vum(l1,l2,M,i,j) = (ui_l1 v1_M uj_l2)
    !l Local variables
    !l   ul    :radial wave functions; see makusp for definition
    !l   sl    :radial wave functions; see makusp for definition
    !l   ruu   :diagonal w.f. products; see makusp for definition
    !l   rus   :diagonal w.f. products; see makusp for definition
    !l   rss   :diagonal w.f. products; see makusp for definition
    !r Remarks
    !r   Diagonal integrals (l1=l2, l=0 potential) are treated specially:
    !r   small component of w.f. is explicitly taken into account.
    !r   This routine does not rely on any special properties of (u,s,gz)
    implicit none
    integer :: lmxa,lmxl,nlml,nr,nsp
    real(8) :: rmt,a,z,rofi(nr),rwgt(nr),rsml(n0),ehl(n0), &
         v0(nr,nsp),v1(nr,nlml,nsp),pnu(n0,nsp),pnz(n0,nsp), &
         qum(0:lmxa,0:lmxa,0:lmxl,3,3,nsp), &
         vum(0:lmxa,0:lmxa,nlml,3,3,nsp),&
         pi,srfpi,xx, &
         quu,qus,qss,quz,qsz,qzz,vuu,vus,vss,vv,vuz,vsz,vzz,vvv(nr), &
         ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp), &
         sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp), &
         gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)
    logical :: lpz1,lpz2
    integer :: l1,l2,lm,i,mlm,isp
    call tcn('momusl')
    pi = 4d0*datan(1d0)
    srfpi = dsqrt(4d0*pi)
    gz=0d0
    call makusp(n0,z,nsp,rmt,lmxa,v0,a,nr,pnu,pnz,rsml,ehl, ul,sl,gz,ruu,rus,rss)
    ! --- Moments ---
    do  isp = 1, nsp
       do  l1 = 0, lmxa
          do  l2 = 0, lmxa
             lpz1 = pnz(l1+1,1) .ne. 0
             lpz2 = pnz(l2+1,1) .ne. 0
             quu = 0d0
             qus = 0d0
             qss = 0d0
             quz = 0d0
             qsz = 0d0
             qzz = 0d0
             do  lm = 0, lmxl
                if (l1 == l2 .AND. lm == 0) then
                   quu= sum(rwgt*ruu(:,l1,1,isp))
                   qus= sum(rwgt*rus(:,l1,1,isp))
                   qss= sum(rwgt*rss(:,l1,1,isp))
                   if (lpz1) then
                      quz= sum(rwgt*ruu(:,l1,2,isp))
                      qsz= sum(rwgt*rus(:,l1,2,isp))
                      qzz= sum(rwgt*rss(:,l1,2,isp))
                   endif
                else
                   quu= sum(rwgt*rofi**lm*ul(:,l1,isp)*ul(:,l2,isp))
                   qus= sum(rwgt*rofi**lm*ul(:,l1,isp)*sl(:,l2,isp))
                   qss= sum(rwgt*rofi**lm*sl(:,l1,isp)*sl(:,l2,isp))
                   if (lpz2) then
                      quz= sum(rwgt*rofi**lm*ul(:,l1,isp)*gz(:,l2,isp))
                      qsz= sum(rwgt*rofi**lm*sl(:,l1,isp)*gz(:,l2,isp))
                      qzz= sum(rwgt*rofi**lm*gz(:,l1,isp)*gz(:,l2,isp))
                   endif
                endif
                qum(l1,l2,lm,1,1,isp) = quu
                qum(l1,l2,lm,1,2,isp) = qus
                qum(l1,l2,lm,2,2,isp) = qss
                qum(l1,l2,lm,1,3,isp) = quz
                qum(l1,l2,lm,2,3,isp) = qsz
                qum(l1,l2,lm,3,3,isp) = qzz
                qum(l2,l1,lm,2,1,isp) = qus !transposed symmetric
                qum(l2,l1,lm,3,1,isp) = quz !
                qum(l2,l1,lm,3,2,isp) = qsz !
             enddo
          enddo
       enddo
    enddo
    ! --- Integrals ul*ul*V_M, ul*sl*V_M, sl*sl*V_M ---
    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  l1 = 0, lmxa
             do  l2 = 0, lmxa
                lpz1 = pnz(l1+1,1) .ne. 0
                lpz2 = pnz(l2+1,1) .ne. 0
                vuu = 0d0
                vus = 0d0
                vss = 0d0
                vuz = 0d0
                vsz = 0d0
                vzz = 0d0
                vvv = v1(:,mlm,isp)
                if (mlm == 1) vvv(2:) = vvv(2:) - [(srfpi*2d0*z/rofi(i),i=2,nr)]
                vvv(1)=0d0
                if (l1 == l2 .AND. mlm == 1) then
                   vuu= sum(rwgt*vvv* ruu(:,l1,1,isp))
                   vus= sum(rwgt*vvv* rus(:,l1,1,isp))
                   vss= sum(rwgt*vvv* rss(:,l1,1,isp))
                   if (lpz1) then
                      vuz= sum(rwgt*vvv* ruu(:,l1,2,isp))
                      vsz= sum(rwgt*vvv* rus(:,l1,2,isp))
                      vzz= sum(rwgt*vvv* rss(:,l1,2,isp))
                   endif
                else
                    vuu= sum(rwgt*vvv* ul(:,l1,isp)*ul(:,l2,isp))
                    vus= sum(rwgt*vvv* ul(:,l1,isp)*sl(:,l2,isp))
                    vss= sum(rwgt*vvv* sl(:,l1,isp)*sl(:,l2,isp))
                    if (lpz2) then
                       vuz= sum(rwgt*vvv* ul(:,l1,isp)*gz(:,l2,isp))
                       vsz= sum(rwgt*vvv* sl(:,l1,isp)*gz(:,l2,isp))
                       vzz= sum(rwgt*vvv* gz(:,l1,isp)*gz(:,l2,isp))
                    endif
                endif
                vum(l1,l2,mlm,1,1,isp) = vuu
                vum(l1,l2,mlm,1,2,isp) = vus
                vum(l1,l2,mlm,2,2,isp) = vss
                vum(l1,l2,mlm,1,3,isp) = vuz
                vum(l1,l2,mlm,2,3,isp) = vsz
                vum(l1,l2,mlm,3,3,isp) = vzz
                vum(l2,l1,mlm,2,1,isp) = vus  !transposed symmetric
                vum(l2,l1,mlm,3,1,isp) = vuz !
                vum(l2,l1,mlm,3,2,isp) = vsz !
             enddo
          enddo
       enddo
    enddo
    call tcx('momusl')
  end subroutine momusl
  subroutine vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,phzdphz,rotp,vumm)!- Rotate vorb from (phi,phidot) to (u,s) and store in vumm
    use m_lmfinit,only: nppn,stdo
    use m_ftox
    !i   lmaxu :dimensioning parameter for U matrix
    !i   lmxa  :augmentation l-cutoff
    !i   vorb  :orbital-dependent potential matrices
    !i   phzdphz  : phz dphz
    !o Inputs/Outputs
    !o  iblu  :index to current LDA+U block
    !          :on input, index to last LDA+U block that was accessed
    ! WARN     :iblu will be incremented to from blocks at this site
    !o Outputs
    !o   vumm  :vorb for this site in (us) representation u_i=(u,s,gz)
    !o         :vumm(m1,m2,i,j) = <u_i| vorb(m1,m2) |u_j>
    implicit none
    integer :: lmaxu,lmxa,iblu,idu(4),m1,m2,l,i
    real(8):: rmt,phi,dlphi,phip,dlphip,dphi,dphip, r12,r21,r11,r22,det, phz,dphz,&
         phzdphz(nppn,n0,2),rotpp(2,2),rotppt(2,2),rotp(0:lmxa,2,2,2) !nsp=2 expected
    complex(8):: vzz,vuz,vsz,vzu,vzs, Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,*), vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,3,3,2,0:lmaxu)
    if(lmaxu>lmxa) call rx('vlm2ux:lmaxu>lmxa')
    vumm=0d0 !corrected at 2023feb3 
    do  l = 0, lmaxu !min(lmxa,3) ! ... Rotate Vorb from phi,phidot basis to u,s basis
       if (idu(l+1) == 0) cycle
       iblu = iblu+1
       do  i = 1, 2
          rotpp  = rotp(l,i,:,:)
          rotppt = transpose(rotpp) !
          ! Vorb is for E_vorb= a_phi *Vorb *a_phi.
          !   [a_phi,a_phidot]= matmul([au,as],rotpp)
          !  Thus we have E_vorb = [a_phi,a_phidot] * vumm* [a_phi,a_phidot]^t
          ! NOTE Vorb is for <phi|Vorb|phi>. vumm is for (u,s,pz), where u,s pz contains phi components.
          vumm(:,:,1,1,i,l) = rotpp(1,1)*Vorb(:,:,i,iblu)*rotppt(1,1)
          vumm(:,:,1,2,i,l) = rotpp(1,1)*Vorb(:,:,i,iblu)*rotppt(1,2)
          vumm(:,:,2,1,i,l) = rotpp(2,1)*Vorb(:,:,i,iblu)*rotppt(1,1)
          vumm(:,:,2,2,i,l) = rotpp(2,1)*Vorb(:,:,i,iblu)*rotppt(1,2)
          phz  = phzdphz(1,l+1,i)
          dphz = phzdphz(2,l+1,i)
          if (phz /= 0) then ! (au-phz, as-dphz)*vumm*(au-phz, as-dphz)^t for gz is expanded to be
             vumm(:,:,1:2,3,i,l)= - phz*vumm(:,:,1:2,1,i,l) - dphz*vumm(:,:,1:2,2,i,l) !vuz,vsz
             vumm(:,:,3,1:2,i,l)= - phz*vumm(:,:,1,1:2,i,l) - dphz*vumm(:,:,2,1:2,i,l) !vzu,vzs
             vumm(:,:,3,3,i,l)  =   phz**2*vumm(:,:,1,1,i,l) + &
                  phz*dphz*(vumm(:,:,1,2,i,l)+vumm(:,:,2,1,i,l)) + dphz**2*vumm(:,:,2,2,i,l) !vzz
          endif
       enddo
    enddo
  end subroutine vlm2us
end module m_augmat
