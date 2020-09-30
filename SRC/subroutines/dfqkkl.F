      subroutine dfqkkl ( nbas , ssite , sspec , numq , sv_p_oqkkl )
      use m_lmfinit,only: nkaph,nspc,nsp
      use m_struc_def  !Cgetarg
C- Allocates arrays to accumulate local output density
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
Ci   numq  :number of Fermi levels for which to accumulate c.d.
Co Outputs
Co   oqkkl :memory is allocated for qkkl
Cl Local variables
Cl   nkapi :number of envelope function types per l q.n. for spec is2
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle lmxa=-1 -> no allocation
Cu   15 Jun 05 Allocation for noncollinear case
Cu   25 Aug 01 Extended to local orbitals
Cu   15 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp df_qkkl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer:: nbas , numq
      type(s_rv1) :: sv_p_oqkkl(3,nbas)

      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer ib,igetss,is,kmax,lmxa,lmxh,nelt1,nlma,nlmh,nelt3,nelt2
c      nsp = globalvariables%nsp
c      nspc = globalvariables%nspc
c      nkaph = globalvariables%nkaph
C --- Loop over sites, allocating qkkl for each site ---
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        lmxh=sspec(is)%lmxb
        kmax=sspec(is)%kmxt
        if (lmxa .eq. -1) goto 10
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
C   ... Case Pkl*Pkl
        nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
        allocate(sv_p_oqkkl(1,ib)%v(abs(nelt1*numq*nsp*nspc)))
C   ... Case Pkl*Hsm
        nelt2 = (kmax+1)*nkaph*nlma*nlmh
        allocate(sv_p_oqkkl(2,ib)%v(abs(nelt2*numq*nsp*nspc)))
C   ... Case Hsm*Hsm
        nelt3 = nkaph*nkaph*nlmh*nlmh
        allocate(sv_p_oqkkl(3,ib)%v(abs(nelt3*numq*nsp*nspc)))
c|        write(6,836) nelt1,nelt3,nelt2
c|  836   format('   nelt=',3i6)
   10   continue
      enddo

      end subroutine dfqkkl


