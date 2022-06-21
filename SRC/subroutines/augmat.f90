subroutine augmat ( z , rmt , rsma , lmxa , pnu , pnz , kmax &
     , nlml , a , nr , nsp , lso , rofi , rwgt , cg , jcg , indxcg &
     , v0 , v1 , v2 , gpotb , gpot0 , nkaph , nkapi , lmxh , lh , &
     eh , rsmh , ehl , rsml , rs3 , vmtz ,  lmaxu , vorb ,  &
     lldau , iblu , idu , sv_p_osig , sv_p_otau , sv_p_oppi,ohsozz,ohsopm, ppnl &
     , hab , vab , sab )
  use m_lmfinit,only: n0,nkap0,nppn,nab
  use m_struc_def, only: s_rv1,s_cv1,s_sblock
  !- Make augmentation matrices sig,tau,pi for one site
  ! ----------------------------------------------------------------------
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
  !i  lcplxp :0 if ppi is real; 1 if ppi is complex
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
  ! o  rofi  :radial mesh points.  On input, rofi(1..nr) are made.
  ! o        :if V is to be extrapolated outside its MT sphere, to
  ! o        :V(1..nrbig), rofi(nr+1,nrbig) are also generated
  ! o        :Thus MUST be dimensioned at least rofi(1..nrbig)
  ! o        :nrbig is internally generated, but will not
  ! o        :exceed parameter nrx defined vxtrap.
  !o Outputs
  !o   osig  :augmentation overlap integrals; see Remarks.
  !o   otau  :augmentation kinetic energy integrals; see Remarks.
  !o   oppi  :augmentation kinetic + potential integrals; see Remarks.
  !o   ppnl  :NMTO potential parameters
  !o   hab   :matrix elements of the ham. with true w.f.  See Remarks.
  !o   vab   :matrix elements of the pot. with true w.f.  See Remarks.
  !o   sab   :matrix elements of    unity with true w.f.  See Remarks.
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
  !r  2. val,slo matches sm hankel                (10s digit pnz=1)
  !r  3. val,slo matches sm hankel, perturbative  (10s digit pnz=2)
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
  !r
  !u Updates
  !u   09 Nov 05 (wrl) Convert dmat to complex form
  !u   27 Apr 05  LDA+U (Lambrecht)
  !u   24 Dec 04 (A. Chantis) ppi matrix elements for full L.S
  !u    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
  !u   12 Aug 04 First implementation of extended local orbitals
  !u   15 Jul 04 (Chantis) radial integrals for spin-orbit coupling
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   27 Aug 01 Extended to local orbitals.  Altered argument list.
  !u   20 Feb 01 Added ppnl to potential parameters generated
  !u   13 Jun 00 spin polarized
  !u   17 May 00 Adapted from nfp augmats.f
  ! ----------------------------------------------------------------------

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
  
  
  implicit none
  integer :: lmxa,kmax,nlml,nr,nsp,nkaph,nkapi,lmxh,lso!,lcplxp
  integer :: lmaxu,lldau,iblu,idu(4)
  complex(8):: vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
  double precision :: z,rmt,rsma,a
  integer:: jcg(1) , indxcg(1) , lh(nkap0)
  type(s_cv1) :: sv_p_oppi(3)
  type(s_sblock):: ohsozz(3),ohsopm(3)
  type(s_rv1) :: sv_p_otau(3)
  type(s_rv1) :: sv_p_osig(3)
  double precision :: rofi(nr),rwgt(nr),v0(nr,nsp),cg(1), &
       pnu(n0,nsp),pnz(n0,nsp),ppnl(nppn,n0,2), &
       v1(nr,nlml,nsp),v2(nr,nlml,nsp),gpot0(nlml),gpotb(nlml), &
       hab(nab,n0,nsp),vab(nab,n0,nsp),sab(nab,n0,nsp), &
       eh(n0,nkaph),rsmh(n0,*),ehl(n0),rsml(n0), &
       rs3,vmtz
  ! ... Local parameters
  integer :: k,ll,lmxl,nlma,nlmh,i
  double precision :: pp(n0,2,5)
  integer :: lxa(0:kmax)
  double precision :: vdif(nr*nsp),sodb(nab,n0,nsp,2), &
       vum((lmxa+1)**2*nlml*6*nsp), &
       fh(nr*(lmxh+1)*nkap0),xh(nr*(lmxh+1)*nkap0), &
       vh((lmxh+1)*nkap0),fp(nr*(lmxa+1)*(kmax+1)), &
       dh((lmxh+1)*nkap0),xp(nr*(lmxa+1)*(kmax+1)), &
       vp((lmxa+1)*(kmax+1)),dp((lmxa+1)*(kmax+1))
  complex(8):: vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
  real(8),allocatable:: qum(:)
  real(8),parameter:: pi   = 4d0*datan(1d0),  y0   = 1d0/dsqrt(4d0*pi)
  call tcn('augmat')
  do  k = 1, lmxh+1 ! check; see description of rsmh above
     if (pnz(k,1) /= 0 .AND. pnz(k,1) < 10 .AND. rsmh(k,nkapi+1) /= 0) &
          call rx1('augmat: illegal value for rsmh',rsmh(k,nkapi+1))
  enddo
  nlma = (lmxa+1)**2
  lmxl = ll(nlml)
  lxa=lmxa
  allocate(qum((lmxa+1)**2*(lmxl+1)*6*nsp))
  do  i = 1, nsp
     vdif(1+nr*(i-1):1+nr*i)= y0*v1(1:nr,1,i) - v0(1:nr,i)
  enddo
  ! --- Make hab,vab,sab and potential parameters pp ---
  call potpus(z,rmt,lmxa,v0,vdif,a,nr,nsp,lso,rofi,pnu,pnz,ehl,rsml, &
       rs3,vmtz,nab,n0,pp,ppnl,hab,vab,sab,sodb)
  ! --- Moments and potential integrals of ul*ul, ul*sl, sl*sl ---
  call momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp,rofi, &
       rwgt,v0,v1,qum,vum)
  ! --- Set up all radial head and tail functions, and their BC's ---
  nlmh = (lmxh+1)**2
  fh=0d0
  call fradhd(nkaph,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh)
  call fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)
  ! ... LDA+U: rotate vorb from (phi,phidot) to (u,s) for all l with U at this site and store in vumm
  if (lldau > 0) call vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,ppnl,vumm)
  ! ... Pkl*Pkl !tail x tail
  call gaugm ( nr , nsp , lso , rofi , rwgt , lmxa , lmxl &
  , nlml , v2 , gpotb , gpot0 , hab , vab , sab , sodb , qum , &
       vum , cg , jcg , indxcg , kmax + 1 , kmax + 1 , lmxa , lxa , &
       fp , xp , vp , dp , kmax + 1 , kmax + 1 , lmxa , lxa , fp , xp, vp , dp , lmxa , &
       sv_p_osig(1)%v , sv_p_otau (1)%v , nlma , nlma , &
       sv_p_oppi(1)%cv,ohsozz(1)%sdiag,ohsopm(1)%soffd , lmaxu , vumm , lldau , idu )
  ! ... Hsm*Pkl !head x tail
  call gaugm ( nr , nsp , lso ,  rofi , rwgt , lmxa , lmxl &
  , nlml , v2 , gpotb , gpot0 , hab , vab , sab , sodb , qum , &
       vum , cg , jcg , indxcg , nkaph , nkapi , lmxh , lh , fh , xh &
       , vh , dh , kmax + 1 , kmax + 1 , lmxa , lxa , fp , xp , vp , dp , lmxh , &
       sv_p_osig(2)%v , sv_p_otau(2) %v , nlmh , nlma, &
       sv_p_oppi(2)%cv,ohsozz(2)%sdiag,ohsopm(2)%soffd, lmaxu , vumm , lldau , idu )
  ! ... Hsm*Hsm !head x head
  call gaugm ( nr , nsp , lso ,  rofi , rwgt , lmxa , lmxl &
  , nlml , v2 , gpotb , gpot0 , hab , vab , sab , sodb , qum , &
       vum , cg , jcg , indxcg , nkaph , nkapi , lmxh , lh , fh , xh &
       , vh , dh , nkaph , nkapi , lmxh , lh , fh , xh , vh , dh , lmxh, &
       sv_p_osig(3)%v, sv_p_otau(3)%v, nlmh,nlmh, &
       sv_p_oppi(3)%cv,ohsozz(3)%sdiag,ohsopm(3)%soffd, lmaxu , vumm , lldau , idu )
  call tcx('augmat')
end subroutine augmat
