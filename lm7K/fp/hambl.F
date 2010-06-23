Cchp1 cmdvar= integer,iv_p_,oiprmb
Cchp1 cmdvar= integer,iv_p_,oindxo
      subroutine hambl(mode,nbas,ssite,sspec,slat,sham,isp,q,k1,k2,k3,
     .smpot,vconst,osig,otau,oppi,lcplxp,alfa,ndimh,napw,igvapw,
     .h,s,hso)

      use m_struc_def  !Cgetarg

C- Make the LDA hamiltonian and overlap for one k-point.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit specifies LDA matrix elements
Ci         :  0 compute LDA hamiltonian and overlap
Ci         :  1 compute the overlap only
Ci         :10s digit
Ci         :  0 do not compute hso
Ci         :  1 compute hs0
Ci         :    Note: only a portion of hso is computed for a
Ci         :    particular isp.  The total hso is assembled
Ci         :    after isp loops from 1..2.  hso should not be
Ci         :    initialized between isp=1 and isp=2 loops.
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oindxo
Ci     Stored:    *
Ci     Passed to: *
Ci   isp   :spin index
Ci   q     :Bloch vector (k-point)
Ci   k1,k2,k3 dimensions of smpot
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst:additional constant potential
Ci   osig,otau,oppi  augmentation matrices
Ci   alfa  :add alfa * overlap to hamiltonian
Ci         :This is for stability in evals.  Preferably alfa=0
Ci   ndimh :dimension of hamiltonian and ovrlap h,s
Co Outputs
Co   h     :Hamiltonian matrix
Co   s     :overlap matrix
Co   hso   :spin off-diagonal block of spin-orbit hamiltonian
Cr Remarks
Cu Updates
Cu   04 Jul 08 (T. Kotani) New PW addition to basis
Cu   03 Feb 05 (A. Chantis) calculate hso
Cu    1 Sep 04 Adapted to handle complex ppi; S.O. put into ppi
Cu   25 Aug 04 modifications for extended local orbitals
Cu   15 Jul 04 (Chantis) Add Lz.Sz spin-orbit coupling
Cu   10 Jan 03 Remove addition from hambl.  See hambls.f
Cu   10 Jan 03 put sigma back into Bloch transform only
Cu   14 Aug 02 Added overlap-only option and option for orthog sigm
Cu   20 Jul 02 Can add Bloch transform of sigma matrix to ham
Cu   18 May 00 Adapted from nfp mk_hamiltonian.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,isp,k1,k2,k3,ndimh,lcplxp,napw,igvapw(3,napw)
      integer osig(3,nbas),otau(3,nbas),oppi(3,nbas)
      double precision alfa
      real(8):: q(3) , vconst
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_ham)::sham
      type(s_lat)::slat

      double complex smpot(k1,k2,k3),h(ndimh,ndimh),s(ndimh,ndimh)
      double complex hso(ndimh,ndimh)
C ... Local parameters
      integer mode0,mode1,ipr
Cchp1       integer oiprmb
       integer,pointer :: iv_p_oiprmb(:)

      double precision vavg
C     double precision plat(3,3)
C ... for sigma
C      integer ndhrs,ontabs,oiaxs,ohrs
C      integer i,hreal,nl,nglob,mxorb,nttabs,ival,nsp,lrsig
C      integer os,oh,oz
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('hambl')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

Cchp1       oiprmb=sham%oindxo
       iv_p_oiprmb => sham%iv_p_oindxo 



C ... Initialization of quantities to be computed
      if (mode0 .eq. 0) call dpzero(h,ndimh*ndimh*2)
      if (mode0 .le. 1) call dpzero(s,ndimh*ndimh*2)
      if (mode1 .eq. 1 .and. isp .eq. 1) call dpzero(hso,ndimh*ndimh*2)
C ... Skip LDA H, S
      if (mode0 .gt. 1) goto 100

C --- Augmentation part of hamiltonian and overlap ---
Cchp1       call augmbl(mode0+10*mode1,ssite,sspec,slat,isp,lcplxp,q,
Cchp1      .osig,otau,oppi,w(oiprmb),ndimh,napw,igvapw,h,hso,s)
       call augmbl ( mode0 + 10 * mode1 , ssite , sspec , slat , isp 
     . , lcplxp , q , osig , otau , oppi , iv_p_oiprmb , ndimh , napw 
     . , igvapw , h , hso , s ) 


C --- Add alfa S+ S ---
      if (dabs(alfa) .gt. 1d-10 .and. mode0 .eq. 0)  then
        call z2herm('U',ndimh,ndimh,s)
        call zgemm('C','N',ndimh,ndimh,ndimh,dcmplx(alfa,0d0),
     .  s,ndimh,s,ndimh,dcmplx(1d0,0d0),h,ndimh)
      endif

C ... Optionally set average istl potential to zero
      vavg = 0
C     call pvhmb1(0,slat,k1,k2,k3,smpot,vavg)
C     call pvhmb1(1,slat,k1,k2,k3,smpot,-vavg)
C --- Smoothed hamiltonian and overlap, S and T done analytically ---
C     print *, '!!'; s = 0

Cchp1       call smhsbl(mode0,ssite,sspec,slat,vavg+vconst,q,ndimh,w(oiprmb),
Cchp1      .napw,igvapw,h,s)
       call smhsbl ( mode0 , ssite , sspec , slat , vavg + vconst , 
     . q , ndimh , iv_p_oiprmb , napw , igvapw , h , s ) 

C     call zprm('s',12,s,ndimh,ndimh,ndimh)

      if (mode0 .eq. 1) goto 100
Cchp1       call hsibl(ssite,sspec,slat,k1,k2,k3,smpot,isp,q,ndimh,w(oiprmb),
Cchp1      .napw,igvapw,h)
       call hsibl ( ssite , sspec , slat , k1 , k2 , k3 , smpot , isp 
     . , q , ndimh , iv_p_oiprmb , napw , igvapw , h ) 


C ... Restore istl potential
C      call pvhmb1(1,slat,k1,k2,k3,smpot,vavg)

C --- Alternatively, smoothed H and S done completely numerically ---
C      call upack1('lat nabc',slat,ngabc)
C      call hsmi_q(ssite,sspec,slat,n1,n2,n3,k1,k2,k3,smpot,
C     .  vconst,q,ndimh,h,s)

C ... End of assembling LDA H, S
  100 continue

C     if (mode0 .eq. 1) goto 200
C --- Bloch transform sigma matrix and add to h ---
C      call upack('ham ndhrs ontabs oiaxs ohrs',sham,
C     .  ndhrs,ontabs,oiaxs,ohrs,0)
C      lrsig = mod(mode/10,10)
C      call isanrg(lrsig,0,2,'hambl:','10s digit mode',.true.)
C      if (lrsig .ne. 0 .and. ohrs .ne. 0) then
C        hreal = 2-mod(lrsig,4)
C        nl = nglob('nl')
C        mxorb = nglob('mxorb')
C        nttabs = ival(w(ontabs),nbas+1)
C        nsp  = nglob('nsp')
C        call isanrg(isp,1,nsp,'hambl:','isp',.true.)
C        call upack('lat plat',slat,plat,0,0,0,0)
C
C        call defcc(oh,ndimh**2)
C        call dcopy(2*ndimh**2,h,1,w(oh),1)
C        call zprm('h',2,h,ndimh,ndimh,ndimh)
C
CC       Call bloch : 104010=perm orb, transpose, c*16
C        i = 100000 + 4000 + 40*(1-hreal) + 10
C        call bloch(i,q,nl,plat,mxorb,w(oiprmb),1,nttabs,w(oiaxs),
C     .    w(ohrs),ndhrs,isp,nsp,ndimh,ndimh,0,ndimh,0,ndimh,0,h,w,w)
C        call zprm('sigm-orth',2,h,ndimh,ndimh,ndimh)
C
CC       try to work with orthogonalized basis.
CC       Rotation to orthogonalized basis doesn't work well.
CC        call defcc(os,ndimh**2)
CC        call dcopy(2*ndimh**2,s,1,w(os),1)
CCC       call zhev2(5,ndimh,ndimh,ndimh,.true.,h,w(os))
CCC       call zprm('L sigm-orth L+',2,h,ndimh,ndimh,ndimh)
CC        call defcc(oz,ndimh**2)
CC        call zhev3(5,ndimh,h,w(os),w(oz))
CC        call rlse(os)
CC
CC        call daxpy(2*ndimh**2,1d0,w(oh),1,h,1)
CCC       call zprm('h(lda)+sig',2,h,ndimh,ndimh,ndimh)
C
C      endif

C --- Fill second half of matrix ---
C 200 continue
      if (mode0 .eq. 0) call z2herm('U',ndimh,ndimh,h)
      if (mode0 .le. 1) call z2herm('U',ndimh,ndimh,s)

C     Test routine that doctors hamiltonian
C     call pvhmb2(ndimh,h,s)

      call getpr(ipr)
      if (ipr .ge. 90) then
        call info(0,0,0,' h and s for q=%3:1;6,6d',q,0)
        if (mode0 .eq. 0) call zprm('h',12,h,ndimh,ndimh,ndimh)
        if (mode0 .le. 1) call zprm('s',12,s,ndimh,ndimh,ndimh)
        if (mode1 .eq. 1) call zprm('hso',2,hso,ndimh,ndimh,ndimh)
      endif
      call tcx('hambl')
      end subroutine hambl


      subroutine pvhmb1(mode,slat,k1,k2,k3,smpot,vavg)

      use m_struc_def  !Cgetarg

C- Add subtract a constant from smooth potential
C ----------------------------------------------------------------------
Ci Inputs
Ci    mode :0 do not adjust smpot; just return avg smpot
Ci         :1 add vavg to smpot
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc
Ci   k1,k2,k3 dimensions of smpot
Cio Inputs/Outputs
Cio  smpot :(input for mode=0, altered for mode=1)
Cio        :smooth potential on uniform mesh (mkpot.f)
Cio  vavg  :(output for mode=0, input for mode=1) average potential
Cr Remarks
Cr
Cu Updates
Cu   16 Aug 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,k1,k2,k3
      type(s_lat)::slat

      double complex smpot(k1,k2,k3)
C ... Local parameters
      double precision vavg
      double complex xx
      integer ngabc(3),i

      if (mode .eq. 0) then

        i_copy_size=size(slat%nabc)
        call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)

        call mshint(1d0,1,ngabc(1),ngabc(2),ngabc(3),
     .  k1,k2,k3,smpot,vavg,xx)
      elseif (mode .eq. 1) then
        do i = 1, k1*k2*k3
          smpot(i,1,1) = smpot(i,1,1) + vavg
        enddo
      else
        call rx('pvhmb1: bad mode')
      endif
      end subroutine pvhmb1


C      subroutine pvhmb2(ndimh,h,s)
CC- test routine to doctor hamiltonian
C      implicit none
C      integer ndimh
C      double complex h(ndimh,ndimh),s(ndimh,ndimh)
C
C      integer low,high,i1,i2,k1,k2
C      return
C
CC ... Zeros out off-diagonal blocks between low+1, high
CC      low = 20
CC      high = 23
CC      print *, 'hambl zero out off-diagonal blocks', low+1,high
CC      do  i1 = 1, ndimh
CC        k1 = -1
CC        if (i1 .gt. low) k1 = 0
CC        if (i1 .gt. high) k1 = 1
CC        do  i2 = 1, ndimh
CC          k2 = -1
CC          if (i2 .gt. low) k2 = 0
CC          if (i2 .gt. high) k2 = 1
CC
CC          if (k1 .eq. 0  .neqv.  k2 .eq. 0) then
CC            h(i1,i2) = 0
CC            s(i1,i2) = 0
CC          endif
CC        enddo
CC      enddo
C
C      end

