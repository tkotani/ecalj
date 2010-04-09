Cgetarg...info...           structure ['sstrxq', 'lat', 'slat']
Cgetarg...info...           structure ['sstrxq', 'ham', 'sham']
Cgetarg...info...           use_to_add ['sstrxq', 'm_struc_def']
      subroutine sstrxq(job,slat,sham,q,ne,ei,tral,trad,sll,sdll,shl,
     .  sdhl)
       
       use m_struc_def  !Cgetarg

C- Screened strux S and optionally Sdot, by Ewald summation
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :1s digit
Ci         :number of energy derivatives (0 or 1)
Ci   slat  :struct containing information about the lattice
Ci   sham  :struct containing information about the hamiltonian
Ci   q     :Bloch vector (k-point)
Ci   ne    :number of energies
Ci   ei    :table of (real) energies
Ci   tral  :transformation matrix for head and tail functions (mktral.f)
Ci         :tral = tral(4,nl,nl,nbas,ne)
Ci   trad  :energy derivative of tral
Co Outputs:
Co   sll   :l+intermediate blocks of S in tral-representation
Co         :sll = sll(lidim,lidim,ne),
Co   sdll  :energy derivative of sll
Co   shl   :hl-block of S in b-representation  (not implemented)
Co         :shl = shl(hdim,lidim,ne)
Co   sdhl  :energy derivative of shl
Cb Bugs
Cb   If the energy derivative is sought, strux for ei=0 are actually
Cb   evaluated for a slightly positive energy, since sdot is singular.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer job,ne
Cgetarg       double precision slat(1),sham(1),q(3),ei(ne),tral(*),trad(*),
Cgetarg      .  sll(*),sdll(*),shl(*),sdhl(*)
       real(8):: q(3) , ei(ne) , tral(*) , trad(*) , sll(*) , sdll(*) 
     . , shl(*) , sdhl(*) 
       type(s_lat)::slat
       type(s_ham)::sham

C Local variables:
      integer hdim,i,ldham(16),ldim,lidim,lihdim,nbas,nderiv,nglob,nkd,
     .  nkq,nl,obas,ocg,ocy,odlv,offt,ofhl,ofll,oidxcg,oiprmb,ojcg,oqlv
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      double precision avw,awald,kap2i,vol,alat,plat(3,3)
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('sstrxq')
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: obas , integer(8) :: slat%opos 1
Cgetarg       call upack('lat alat plat opos',slat,alat,plat,obas,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       obas=slat%opos

Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
       
       ocg=slat%ocg
       ojcg=slat%ojcg
       oidxcg=slat%oidxcg
       ocy=slat%ocy

Cgetarg...info...           integer :: nkd , integer(8) :: slat%nkd 1
Cgetarg...info...           integer :: nkq , integer(8) :: slat%nkq 1
Cgetarg...info...           integer :: odlv , integer(8) :: slat%odlv 1
Cgetarg...info...           integer :: oqlv , integer(8) :: slat%oqlv 1
Cgetarg       call upack('lat nkd nkq odlv oqlv',slat,nkd,nkq,odlv,oqlv,0)
       
       nkd=slat%nkd
       nkq=slat%nkq
       odlv=slat%odlv
       oqlv=slat%oqlv

Cgetarg...info...           real(8) :: awald , real(8) :: slat%awald 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg       call upack('lat awald vol avw',slat,awald,vol,avw,0,0)
       
       awald=slat%awald
       vol=slat%vol
       avw=slat%avw

Cgetarg...info...           integer :: oiprmb , integer(8) :: sham%oindxo 1
Cgetarg...info...           integer :: ldham(16) , integer(8) :: sham%ldham 16
Cgetarg       call upack('ham oindxo ldham',sham,oiprmb,ldham,0,0,0)
       
       oiprmb=sham%oindxo
       i_copy_size=size(sham%ldham) 
       call i8icopy(i_copy_size,sham%ldham,1,ldham,1) 

      nl = nglob('nl')
      nbas = nglob('nbas')
      hdim = lihdim - lidim
      nderiv = mod(job,10)
      call isanrg(nderiv,0,1,'sstrxq:','nderiv',.true.)

C --- For each energy, make S, Sdot ---
      do  i = 1, ne
        kap2i = ei(i)
        if (nderiv .ne. 0 .and. dabs(kap2i) .lt. .01d0) kap2i = .01d0
        ofhl = hdim*lidim*(i-1)
        ofll = lidim*lidim*(i-1)
        offt = 4*nl*nl*nbas*(i-1)
        call sstrq1(alat,avw,w(obas),q,hdim,w(oiprmb),kap2i,nbas,nderiv,
     .    lidim,awald,vol,nkd,nkq,w(odlv),w(oqlv),w(ocg),w(ojcg),
     .    w(oidxcg),nl,tral(1+offt),trad(1+offt),sll(1+2*ofll),
     .    sdll(1+2*ofll),shl(1+2*ofhl),sdhl(1+2*ofhl))
      enddo

      call tcx('sstrxq')
Cgetarg       end
       end subroutine sstrxq 


      subroutine sstrq1(alat,avw,bas,q,hdim,iprmb,kap2,nbas,nderiv,
     .  lidim,awald,vol,nkd,nkq,dlv,qlv,cg,jcg,indxcg,nl,tral,trad,sll,
     .  sdll,shl,sdhl)

C- Screened strux S and Sdot, by Ewald summation, for one kinetic energy
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   bas   :basis vectors, in units of alat
Ci   q     :k-point
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   hdim  :dimension of higher block
Ci   iprmb :iprmb(i) points to row corresponding to i of shuffled matrix
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   kap2  :kinetic energy
Ci   nbas  :size of basis
Ci   nderiv:0 to calculate S only, 1 to calculate both S and Sdot
Ci   lidim :number of lower+intermediate orbitals
Ci   awald :Ewald parameter, scales with the lattice as: as/(vol)**(1/3)
Ci   nkd   :number of direct-space lattice vectors
Ci   nkq   :number of reciprocal-space lattice vectors
Ci   dlv   :direct-space lattice vectors, units of alat
Ci   qlv   :reciprocal lattice vectors, units of 2pi/alat
Ci   nl    :(global maximum l) + 1
Ci   trad  :(kappa*avw)^2-derivative of tral
Ci   tral  :transformation matrix for head and tail functions (mktral.f)
Ci   vol   :cell volume
Co Outputs:
Co   sll   :ll-block of Bloch sum of real-space matrix S-beta
Co   sdll  :energy derivative of sll
Co   shl   :not implemented
Co   sdhl  :not implemented
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nbas,nl,hdim,iprmb(72),indxcg(*),jcg(*),
     .  nderiv,lidim,nkd,nkq
C      integer nbas,nl,hdim,iprmb(nl*nl*nbas),indxcg(*),jcg(*),
C     .  nderiv,lidim,nkd,nkq
      double precision alat,avw,awald,bas(3,*),q(3),cg(*),dlv(3,*),
     .  qlv(3,*),kap2,shl(hdim,lidim,2),sdhl(hdim,lidim,2),
     .  tral(4,nl,nl,nbas),trad(4,nl,nl,nbas),vol
      double complex sll(lidim,lidim),sdll(lidim,lidim)
C Local variables:
      integer i,ia,ib,ib0,iprint,k,lmxa,lmxb,lmxx,mxint,n0,nl2,nlma,
     .  nlmb,nlmx,nkap0,offa,offb,olrfac,owk
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      parameter (lmxx=8,nlmx=(lmxx+1)**2)
      double precision dum,phma(3)
      double complex s(nlmx,nlmx),sdot(nlmx,nlmx)
      character*20 outs
C ... Heap
      integer w(1)
      common /w/ w

C --- Unscreened structure constants for all pairs ---
      nl2 = nl**2
      call dpzero(sll,2*lidim**2)
      if (nderiv .ne. 0) call dpzero(sdll,2*lidim**2)
      do  ia = 1, nbas
        call orbl(ia,0,lidim,iprmb,norb,ltab,ktab,dum,offl,dum)
        lmxa = mxint(norb,ltab)
        nlma = (lmxa+1)**2
        ib0 = ia
        do  ib = ib0, nbas
          call orbl(ib,0,lidim,iprmb,norb,ltab,ktab,dum,offl,dum)
          lmxb = mxint(norb,ltab)
          nlmb = (lmxb+1)**2
          do  k = 1, 3
            phma(k) = bas(k,ib)-bas(k,ia)
          enddo
          call strxq(3,kap2,q,phma,nlma,nlmb,nlmx,alat,vol,awald,nkd,
     .      nkq,dlv,qlv,cg,indxcg,jcg,s,sdot)
          call s2oka(111,nlmx,nlmb,nlma,nlmb,avw,s,s)
          offa = nl2*(ia-1)
          offb = nl2*(ib-1)
          call ztoy(s,nlmx,nlma,nlmb,0)
          call pblch1(1010,nl2,offa,offb,2,lidim,0,iprmb,
     .      0,lidim,0,lidim,s,s,nlmx,1d0,0d0,sll)
          if (nderiv .eq. 1) then
            call s2oka(111,nlmx,nlmb,nlma,nlmb,avw,sdot,sdot)
            offa = nl2*(ia-1)
            offb = nl2*(ib-1)
            call ztoy(sdot,nlmx,nlma,nlmb,0)
            call pblch1(1010,nl2,offa,offb,2,lidim,0,iprmb,
     .        0,lidim,0,lidim,sdot,sdot,nlmx,1d0,0d0,sdll)
          endif
        enddo
      enddo
      call dscal(2*lidim**2,1/avw**2,sdll,1)
C      call zprm('bare sll',2,sll,lidim,lidim,lidim)
C      call zprm('bare sdll',2,sdll,lidim,lidim,lidim)

C --- Scale sll according to tral matrix ---
      call defdr(olrfac,lidim*(nderiv+1))
      call sstrq2(nbas,nderiv,lidim,nl,tral,trad,iprmb,w(olrfac),sll,
     .  sdll)
      call rlse(olrfac)
C     call zprm('sll after scaling',2,sll,lidim,lidim,lidim)
C     call zprm('sdll after scaling',2,sdll,lidim,lidim,lidim)

C --- B = [cd +cS0c] = [S^a + D]/(-det) => make A = B^-1 ---
      call tcn('zqinv')
      call defrr(owk,lidim*(lidim+1))
      call zqinv('h',sll,lidim,0,lidim,w(owk),lidim,i)
      call rlse(owk)
      call tcx('zqinv')

C     call zprm('A',2,sll,lidim,lidim,lidim)

C --- Adot = A Bdot A ---
      if (nderiv .eq. 1) then
        call defcc(owk,lidim*lidim)
        call zgemm('N','N',lidim,lidim,lidim,dcmplx(1d0,0d0),
     .    sdll,lidim,sll,lidim,(0d0,0d0),w(owk),lidim)
        call zgemm('N','N',lidim,lidim,lidim,dcmplx(1d0,0d0),
     .    sll,lidim,w(owk),lidim,(0d0,0d0),sdll,lidim)

C      call zprm('Adot = A Bdot A',2,sdll,lidim,lidim,lidim)

      endif

C --- Convert B to S^a ---
      call defdr(olrfac,lidim*3)
      call sstrq3(nbas,nderiv,lidim,nl,tral,trad,iprmb,w(olrfac),sll,
     .  sdll)
      call rlse(olrfac)
      if (iprint() .ge. 110) then
        call awrit1('%xS^a, e=%d',outs,20,0,kap2)
        call zprm(outs,2,sll,lidim,lidim,lidim)
        if (nderiv .eq. 1) call zprm('Sd^a',2,sdll,lidim,lidim,lidim)
      endif

Cgetarg       end
       end subroutine sstrq1 


      subroutine sstrq2(nbas,nder,lidim,nl,tral,trad,indxsh,lrfac,sll,
     .  slld)

C- Scales and adds diagonal to Ewald-summed -S0 to make C matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nder  :1 if to scale slld, otherwise 0
Ci   lidim :number of lower+intermediate orbitals
Ci   nl    :(global maximum l) + 1
Ci   tral  :transformation matrix for head and tail functions (mktral.f)
Ci         :tral is in RL order
Ci   trad  :(kappa*avw)^2-derivative of tral
Ci         :trad is in RL order
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lrfac :work array of dimension 2*lidim, holding tral and trad in
Ci         :downfolding order
Co Outputs
Cio   sll  :on input are the bare strux, in downfolding order
Cio        :sll is scaled, i.e. overwritten by t3 t4 - t3 sll t3
Cio   sdll :on input are the energy derivative of sll
Cio        :sdll is overwritten by the energy derivative of scaled sll
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer lidim,nder,nbas,nl,indxsh(*)
      double complex sll(lidim,lidim),slld(lidim,lidim)
      double precision tral(4,nl*nl*nbas),trad(4,nl*nl*nbas),
     .  lrfac(lidim,0:*)
C Local variables:
      integer i,j,lmi,ib,il,im,lmj,jb,jl,jm
      double precision cd(0:1)
      double complex sc(0:1)

C --- c and its energy derivative, in downfolding order ---
      lmi = 0
      do  ib = 1, nbas
      do  il = 0, nl-1
      do  im = -il, il
        lmi = lmi+1
        i = indxsh(lmi)
        if (i .gt. 0 .and. i .le. lidim) then
          lrfac(i,0) = tral(3,lmi)
          if (nder .eq. 1) lrfac(i,1) = trad(3,lmi)
        endif
      enddo
      enddo
      enddo

C --- Overwrite sll with scaled -sll ---
      lmj = 0
      do  jb = 1, nbas
      do  jl = 0, nl-1
      do  jm = -jl, jl
        lmj = lmj+1
        j = indxsh(lmj)
        if (j .gt. 0 .and. j .le. lidim) then

          lmi = 0
          do  ib = 1, jb
          do  il = 0, nl-1
          do  im = -il, il
            lmi = lmi+1
            i = indxsh(lmi)
            if (i .gt. 0 .and. i .le. lidim) then

              if (nder .eq. 1) then
                sc(0) = sll(i,j)*lrfac(j,0)
                sc(1) = sll(i,j)*lrfac(j,1) + slld(i,j)*lrfac(j,0)
                slld(i,j) = -(lrfac(i,1)*sc(0) + lrfac(i,0)*sc(1))
                if (ib .ne. jb) slld(j,i) =  dconjg(slld(i,j))

              endif
              sll(i,j) = -lrfac(i,0)*sll(i,j)*lrfac(j,0)
              if (ib .ne. jb) sll(j,i) =  dconjg(sll(i,j))

            endif
          enddo
          enddo
          enddo
        endif

      enddo
      enddo
      enddo

C --- Add diagonal ---
      lmi = 0
      do  ib = 1, nbas
      do  il = 0, nl-1
      do  im = -il, il
        lmi = lmi+1
        i = indxsh(lmi)
        if (i .gt. 0 .and. i .le. lidim) then
          cd(0) = tral(3,lmi)*tral(4,lmi)
          sll(i,i) = sll(i,i) + cd(0)
          if (nder .eq. 1) then
            cd(1) = trad(3,lmi)*tral(4,lmi) + tral(3,lmi)*trad(4,lmi)
            slld(i,i) = slld(i,i) + cd(1)
          endif
        endif
      enddo
      enddo
      enddo

Cgetarg       end
       end subroutine sstrq2 


      subroutine sstrq3(nbas,nder,lidim,nl,tral,trad,indxsh,fac,sll,
     .  sdll)

C- Scales and adds diagonal to C matrix to make S^a
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nder  :1 if to scale sdll, otherwise 0
Ci   lidim :number of lower+intermediate orbitals
Ci   nl    :(global maximum l) + 1
Ci   tral  :transformation matrix for head and tail functions (mktral.f)
Ci         :tral is in RL order
Ci   trad  :(kappa*avw)^2-derivative of tral
Ci         :trad is in RL order
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   fac :work array of dimension 2*lidim, holding tral and trad in
Ci         :downfolding order
Co Outputs
Cio   sll  :on input are the bare strux, in downfolding order
Cio        :sll is scaled, i.e. overwritten by t3 t4 - t3 sll t3
Cio   sdll :on input are the energy derivative of sll
Cio        :sdll is overwritten by the energy derivative of scaled sll
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C- Passed variables:
      integer lidim,nder,nbas,nl,indxsh(*)
      double complex sll(lidim,lidim),sdll(lidim,lidim)
      double precision tral(4,nl*nl*nbas),trad(4,nl*nl*nbas),
     .  fac(lidim,3)
C- Local variables:
      integer i,j,lmi,ib,il,im,lmj,jb,jl,jm

C --- c and its energy derivative, in downfolding order ---
      lmi = 0
      do  ib = 1, nbas
      do  il = 0, nl-1
      do  im = -il, il
        lmi = lmi+1
        i = indxsh(lmi)
        if (i .gt. 0 .and. i .le. lidim) then
          fac(i,1) = tral(1,lmi)*tral(4,lmi)-tral(2,lmi)*tral(3,lmi)
          fac(i,2) = tral(1,lmi)/tral(3,lmi)
          if (nder .eq. 1) fac(i,3) = trad(1,lmi)/tral(3,lmi) -
     .      tral(1,lmi)*trad(3,lmi)/tral(3,lmi)**2
        endif
      enddo
      enddo
      enddo

C --- Overwrite sll with scaled -sll*dt ---
      lmj = 0
      do  jb = 1, nbas
      do  jl = 0, nl-1
      do  jm = -jl, jl
        lmj = lmj+1
        j = indxsh(lmj)
        if (j .gt. 0 .and. j .le. lidim) then

          lmi = 0
          do  ib = 1, nbas
          do  il = 0, nl-1
          do  im = -il, il
            lmi = lmi+1
            i = indxsh(lmi)
C           Scale by determinant
            if (i .gt. 0 .and. i .le. lidim) then
              sll(i,j) = -sll(i,j)*fac(j,1)
              if (nder .eq. 1) sdll(i,j) = sdll(i,j)*fac(j,1)
            endif
          enddo
          enddo
          enddo
        endif

      enddo
      enddo
      enddo

C --- Add diagonal ---
      call daxpy(lidim,1d0,fac(1,2),1,sll,2*(1+lidim))
      call daxpy(lidim,1d0,fac(1,3),1,sdll,2*(1+lidim))

Cgetarg       end
       end subroutine sstrq3 



