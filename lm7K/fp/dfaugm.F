Cchp3 cmdvar= type(s_rv1),sv_p_,oppi
Cchp3 cmdvar= type(s_rv1),sv_p_,otau
Cchp3 cmdvar= type(s_rv1),sv_p_,osig
Cchp3       subroutine dfaugm(nbas,lcplxp,lso,ssite,sspec,osig,otau,oppi)
Cchp3        subroutine dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig 
Cchp3      . , otau , oppi ) 
Cchp3         subroutine dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig 
Cchp3      .  , sv_p_otau , oppi ) 
         subroutine dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig 
     .   , sv_p_otau , sv_p_oppi ) 




      use m_struc_def  !Cgetarg

C- Allocate augmentation matrices sigma,tau,pi for all atoms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   lcplxp:if nonzero allocate space for complex ppi
Ci   lso   :flag for SO coupling
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa lmxb kmxt
Ci     Stored:
Ci     Passed to:
Co Outputs
Co   osig  :memory allocated
Co   otau  :memory allocated
Co   oppi  :memory allocated
Cr Remarks
Cr   Pointers are specified as osig(itype,ibas) where
Cr     itype=1: case Pkl*Pkl
Cr     itype=2: case Pkl*Hsm
Cr     itype=3: case Hsm*Hsm
Cr   sig and tau are l diagonal, ppi is full matrix
Cr   Thus integral (P~_kL P~_k'L' - P_kL P_k'L') is diagonal in LL',
Cr       sig(nf1,nf2,0..lmax) with lmax the l-cutoff
Cr   For sig(Pkl,Pkl), nf1=nf2==1+kmax; lmax=lmxa
Cr   For sig(Hsm,Pkl), nf1=nkaph and nf2=1+kmax; lmax=lmxh
Cr   For sig(Hsm,Hsm), nf1=nf2=nkaph; lmax = lmxh
Cl Local variables
Cl   nkapi :number of envelope function types per l q.n. for spec is2
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cl   nso   :1 if lso = 0, otherwise 2 (for allocating pi matrix)
Cu Updates
Cu   01 Jul 05 handle lmxa=-1 -> no allocation
Cu   29 Jun 05 Adapted to store SO in ppi separately from ppi
Cu    1 Sep 04 Adapted to handle complex ppi.  so folded into ppi
Cu   29 Jun 04 (A. Chantis) memory allocation for LzSz matrix elements.
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   25 Aug 01 Extended to local orbitals
Cu   11 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp df_augm.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
Cchp3       integer nbas,lcplxp,lso,osig(3,nbas),otau(3,nbas),oppi(3,nbas)
Cchp3        integer:: nbas , lcplxp , lso , otau(3,nbas) , oppi(3,nbas) 
Cchp3         integer:: nbas , lcplxp , lso , oppi(3,nbas) 
         integer:: nbas , lcplxp , lso 
         type(s_rv1) :: sv_p_oppi(3,nbas)

        type(s_rv1) :: sv_p_otau(3,nbas)

       type(s_rv1) :: sv_p_osig(3,nbas)

      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer ib,igetss,is,kmax,lmxa,lmxh,nelt1,nelt2,nglob,nlma,nlmh,
     .nsp,nspc,nkaph,nso
C     double precision z,rmt

      nsp = nglob('nsp')
      nspc = nglob('nspc')
      nkaph = nglob('nkaph')
c      print *,'lso nspc=',lso,nspc
      call rxx(lso.eq.1.and.nspc.ne.2,'dfaugm: inconsistent parameters')
      nso = 1
      if (lso .ne. 0) nso = 2

C --- Loop over sites, allocating sig,tau,pi for each site ---
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)


        lmxa=sspec(is)%lmxa
        lmxh=sspec(is)%lmxb
        kmax=sspec(is)%kmxt

C       call upack('spec z rmt',sspec,is,z,rmt,0,0)
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        if (lmxa .eq. -1) goto 10

C   ... Case Pkl*Pkl
        nelt1 = (kmax+1)*(kmax+1)*(lmxa+1)*nsp
        nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp*nspc*nso
Cchp3         call defrr(osig(1,ib), nelt1)
         allocate(sv_p_osig(1,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_osig(1,ib)%v(:)=0.0d0

Cchp3         call defrr(otau(1,ib), nelt1)
         allocate(sv_p_otau(1,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_otau(1,ib)%v(:)=0.0d0

        if (lcplxp .eq. 0) then
Cchp3           call defrr(oppi(1,ib), nelt2)
           allocate(sv_p_oppi(1,ib)%v(abs(nelt2)))
           if (nelt2<0) sv_p_oppi(1,ib)%v(:)=0.0d0

        else
Cchp3           call defrr(oppi(1,ib), nelt2*2)
           allocate(sv_p_oppi(1,ib)%v(abs(nelt2*2)))
           if (nelt2*2<0) sv_p_oppi(1,ib)%v(:)=0.0d0

        endif

C   ... Case Hsm*Hsm
        nelt1 = nkaph*nkaph*(lmxh+1)*nsp
        nelt2 = nkaph*nkaph*nlmh*nlmh*nsp*nspc*nso
Cchp3         call defrr(osig(3,ib), nelt1)
         allocate(sv_p_osig(3,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_osig(3,ib)%v(:)=0.0d0

Cchp3         call defrr(otau(3,ib), nelt1)
         allocate(sv_p_otau(3,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_otau(3,ib)%v(:)=0.0d0

        if (lcplxp .eq. 0) then
Cchp3           call defrr(oppi(3,ib), nelt2)
           allocate(sv_p_oppi(3,ib)%v(abs(nelt2)))
           if (nelt2<0) sv_p_oppi(3,ib)%v(:)=0.0d0

        else
Cchp3           call defrr(oppi(3,ib), nelt2*2)
           allocate(sv_p_oppi(3,ib)%v(abs(nelt2*2)))
           if (nelt2*2<0) sv_p_oppi(3,ib)%v(:)=0.0d0

        endif

C   ... Case Hsm*Pkl
        if (lmxh.gt.lmxa) call rx('dfaugm: lmxh.gt.lmxa unexpected')
        nelt1 = nkaph*(kmax+1)*(lmxh+1)*nsp
        nelt2 = nkaph*(kmax+1)*nlmh*nlma*nsp*nspc*nso
Cchp3         call defrr(osig(2,ib), nelt1)
         allocate(sv_p_osig(2,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_osig(2,ib)%v(:)=0.0d0

Cchp3         call defrr(otau(2,ib), nelt1)
         allocate(sv_p_otau(2,ib)%v(abs(nelt1)))
         if (nelt1<0) sv_p_otau(2,ib)%v(:)=0.0d0

        if (lcplxp .eq. 0) then
Cchp3           call defrr(oppi(2,ib), nelt2)
           allocate(sv_p_oppi(2,ib)%v(abs(nelt2)))
           if (nelt2<0) sv_p_oppi(2,ib)%v(:)=0.0d0

        else
Cchp3           call defrr(oppi(2,ib), nelt2*2)
           allocate(sv_p_oppi(2,ib)%v(abs(nelt2*2)))
           if (nelt2*2<0) sv_p_oppi(2,ib)%v(:)=0.0d0

        endif

   10   continue
      enddo

      end subroutine dfaugm


