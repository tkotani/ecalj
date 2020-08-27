      subroutine surho ( nbas , ssite , sspec , lmet , ldos , lrout 
     .  , lekkl , numq , sv_p_oqkkl , sv_p_oeqkkl , k1 , k2 , k3 , smrout 
     .  , ndos , dos , sumev , sumqv )



      use m_struc_def  !Cgetarg
      use m_lmfinit,only: globalvariables
c      use m_globalvariables
C- Initialize output density coeffs, dos, and eigenvalue sum
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa lmxb kmxt
Ci     Stored:    *
Ci     Passed to: *
Ci   lmet  :nonzero => assume metal and init dos
Ci   ldos  :nonzero => init dos
Ci   lrout :nonzero => init qkkl,smrout
C    lekkl :nonzero => init eqkkl
Ci   numq  :number of trial Fermi levels
Ci   oqkkl :offsets to the local (site) density matrices
Ci   k1..3 :dimensions smrout
Ci   smrout:output density
Ci   ndos  :number of energy mesh points
Ci   dos   :density of states
Ci   sumev :band sum of eigenvalues
Ci   sumqv :band sum of charges
Co Outputs
Co   See Remarks
Cr Remarks
Cr   The following steps are taken:
Cr     1.  Local arrays qkkl are zeroed
Cr     2.  smrout is zeroed
Cr     3.  DOS and eigenvalue sums are zeroed.
Cu Updates
Cu   27 Jul 08 (T. Kotani) added eqkkl
Cu   02 Jan 06 sumqv doubled for two spin channels
Cu   01 Jul 05 handle lmxa=-1 -> no local arrays
Cu   15 Jun 05 qkkl zeroed for noncollinear case
Cu   23 Jan 01 Added lrout switch
Cu   19 Jun 00 spin polarized
Cu    5 May 00 Adapted from nfp init_density
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,ldos,lmet,lrout,lekkl,nbas,ndos,numq
      type(s_rv1) :: sv_p_oeqkkl(3,nbas)

      type(s_rv1) :: sv_p_oqkkl(3,nbas)

      real(8):: sumev(2,3) , sumqv(3,2) , dos(ndos,2)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double complex smrout(k1,k2,k3,numq)
C ... Local parameters
      integer ib,i,is,kmax,lmxa,lmxh,nelt1,nelt2,nelt3,nlma,nlmh,igetss
      integer nglob,nsp,nspc,nkaph
C ... Heap

Changenglob      nsp = nglob('nsp')
      nsp = globalvariables%nsp
Changenglob      nspc = nglob('nspc')
      nspc = globalvariables%nspc
Changenglob      nkaph = nglob('nkaph')
      nkaph = globalvariables%nkaph

C ... Reset coeffs of local density
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)


        lmxa=sspec(is)%lmxa
        lmxh=sspec(is)%lmxb
        kmax=sspec(is)%kmxt

        if (lmxa .eq. -1) goto 10

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2

        nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
        nelt3 = nkaph*nkaph*nlmh*nlmh
        nelt2 = (kmax+1)*nkaph*nlma*nlmh

        if (lrout .ne. 0) then
          call dpzero ( sv_p_oqkkl( 1 , ib )%v , nelt1 * numq * nsp * nspc 
     .     )

          call dpzero ( sv_p_oqkkl( 2 , ib )%v , nelt2 * numq * nsp * nspc 
     .     )

          call dpzero ( sv_p_oqkkl( 3 , ib )%v , nelt3 * numq * nsp * nspc 
     .     )

        endif
        if (lrout .ne. 0 .and. lekkl .ne. 0) then
          call dpzero ( sv_p_oeqkkl( 1 , ib )%v , nelt1 * numq * nsp * 
     .     nspc )

          call dpzero ( sv_p_oeqkkl( 2 , ib )%v , nelt2 * numq * nsp * 
     .     nspc )

          call dpzero ( sv_p_oeqkkl( 3 , ib )%v , nelt3 * numq * nsp * 
     .     nspc )

        endif

   10   continue
      enddo

C ... Reset array for smooth density
      if (lrout .ne. 0) call dpzero(smrout,2*k1*k2*k3*numq*nsp)

C ... Dos and eigenvalue sum
      do  i = 1, numq
        sumev(1,i) = 0d0
        sumev(2,i) = 0d0
        sumqv(i,1) = 0d0
        sumqv(i,2) = 0d0
      enddo
      if (ldos.eq.1 .or. lmet.gt.0) call dpzero(dos,2*ndos*nsp)

      end subroutine surho


