      subroutine flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto,
     .ndimh,iprmb,isp,evl,evec,ewgt,cPkL,db,da,wk,ppipp,ppippz,
     .sigpp,ppihp,ppihpz,sighp,lcplxp,f)
C- Force contribution from augmentation at site ia.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ia    :site of augmentation
Ci   kmax  :polynomial cutoff
Ci   nkaph :leading dimension of sigph,ppihp at site ia
Ci   nlmha :dimensions ppihp
Ci   nlma  :augmentation L-cutoff
Ci   lmxa  :augmentation L-cutoff
Ci   ndimh :dimension of hamiltonian
Ci   evl   :eigenvalue
Ci   evec  :eigenvector
Ci   ewgt  :eigenvector weight
Ci   cPkL  :PkL expansion eigenvector at site ia.
Ci   b     :structure constants, needed for PW contribution
Ci   db    :gradient of structure constants, needed to make grad psi
Ci   da    :work array holding grad psi
Ci   wk    :work array holding (ppi-evl*sig)*evec
Ci   ppipp :local tail-tail potential matrix
Ci   ppippz:local tail-tail potential matrix, complex form
Ci         :NB: only ppipp or ppippz is used, depending on lcplxp
Ci   sigpp :local tail-tail overlap matrix
Ci   ppihp :local head-tail potential matrix
Ci   ppihpz:local head-tail potential matrix, complex form
Ci         :NB: only ppihp or ppihpz is used, depending on lcplxp
Ci   sighp :local head-tail overlap matrix
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Co Outputs
Co   f     :local contribution to forces is added
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani) Contribution from PW part of basis
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   17 Jun 00 spin polarized
Cu   25 May 00 Adapted from nfp floc_q.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ia,kmax,lmxa,nkaph,nbas,nlmto,ndimh,nlma,lmxha,lcplxp,
     .nlmha,isp,iprmb(1)
      double precision evl,f(3,nbas),ewgt,
     .sigpp(0:kmax,0:kmax,0:lmxa,isp),sighp(nkaph,0:kmax,0:lmxha,isp),
     .ppihp(nkaph,0:kmax,nlmha,nlma,isp),
     .ppipp(0:kmax,0:kmax,nlma,nlma,isp)
      double complex db(ndimh,nlma,0:kmax,3),da(0:kmax,nlma,3),
     .evec(ndimh),cPkL(0:kmax,nlma),wk(0:kmax,nlma)
      double complex ppihpz(nkaph,0:kmax,nlmha,nlma,isp)
      double complex ppippz(0:kmax,0:kmax,nlma,nlma,isp)
C ... Local parameters
      integer i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,nlm2,norb,n0,nkap0
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      double precision sum,wt,xx

      if (nlmto .eq. 0) return
      call tcn('flocbl')

C ... Make (ppi-evl*sig)*psi in wk
      call flocb2(ia,nlmto,iprmb,kmax,nkaph,nlmha,nlma,evl,evec,
     .ppipp(0,0,1,1,isp),ppippz(0,0,1,1,isp),sigpp(0,0,0,isp),
     .ppihp(1,0,1,1,isp),ppihpz(1,0,1,1,isp),sighp(1,0,0,isp),
     .lcplxp,cPkL,wk)

C ... Loop over ib, virtual shift of wavefct part centered there
      do  ib = 1, nbas
        if (ib .eq. ia) goto 10
C       Block orbitals into groups with consecutive l
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

C   ... Grad of psi expansion coeffs from a virtual shift at site ib
        call dpzero(da, 2*(kmax+1)*nlma*3)
        do  io = 1, norb
          if (blks(io) .ne. 0) then
C         l2 = l index, needed to locate block in da,db
            l2   = ltab(io)
            nlm1 = l2**2+1
            nlm2 = nlm1 + blks(io)-1
            i    = offl(io)
C         i = index to head term, (ilm,k) to augmentation
            do  ilmb = nlm1, nlm2
              i = i+1
              do  m = 1, 3
                do  ilm = 1, nlma
                  do  k = 0, kmax !da can be written as {d cPkL}/{d R}. See rlocb1 in rlocbl.F for generating cPkL. b is used instead of db
                    da(k,ilm,m) = da(k,ilm,m) + evec(i)*db(i,ilm,k,m) 
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo

C   --- Force term is (grad psi_kL) * (ppi-evl*sig)*evec ---
        do  m = 1, 3
          sum = 0d0
          do  ilm = 1, nlma
            do  k = 0, kmax
              sum = sum + 2d0*dconjg(da(k,ilm,m))*wk(k,ilm)
            enddo
          enddo
c|          write (6,300) ia,ib,m,sum
c|  300     format(' ia=',i3,'  ib=',i3,'   m=',i2,'   sum=',2f12.6)
          wt = ewgt
          f(m,ib) = f(m,ib) - wt*sum
          f(m,ia) = f(m,ia) + wt*sum
        enddo

   10   continue
      enddo

C --- Force at site ia from PWs ---
      call dpzero(da, 2*(kmax+1)*nlma*3)
      do  m = 1, 3
        do  k = 0, kmax
          do  ilm = 1, nlma
            do  i = nlmto+1, ndimh
              da(k,ilm,m) = da(k,ilm,m) + evec(i)*db(i,ilm,k,m)
            enddo
          enddo
        enddo
      enddo
C ... Force term is (grad psi_kL) * (ppi-evl*sig)*evec
      do  m = 1, 3
        sum = 0d0
        do  ilm = 1, nlma
          do  k = 0, kmax
            sum = sum + 2d0*dconjg(da(k,ilm,m))*wk(k,ilm)
          enddo
        enddo
        wt = ewgt
        f(m,ia) = f(m,ia) + wt*sum
      enddo

      call tcx('flocbl')
      end

      subroutine flocb2(ia,nlmto,iprmb,kmax,nkaph,nlmha,nlma,evl,evec,
     .ppipp,ppippz,sigpp,ppihp,ppihpz,sighp,lcplxp,cPkL,wk)
C- Make (ppi-evl*sig)*evec
C ----------------------------------------------------------------------
Ci Inputs
Ci   ia    :site of augmentation
Cl   nlmto :dimension of LMTO part of hamiltonian, for hh and ht blocks
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions ppihp
Ci   nlma  :augmentation L-cutoff
Ci   kmax  :polynomial cutoff
Ci   evl   :eigenvalue
Ci   evec  :eigenvector
Ci   nkaph :leading dimension of sigph,ppihp
Ci   ppipp :local tail-tail potential matrix
Ci   ppippz:local tail-tail potential matrix, complex form
Ci         :NB: only ppipp or ppippz is used, depending on lcplxp
Ci   sigpp :local tail-tail overlap matrix
Ci   ppihp :local head-tail potential matrix
Ci   ppihpz:local head-tail potential matrix, complex form
Ci         :NB: only ppihp or ppihpz is used, depending on lcplxp
Ci   sighp :local head-tail overlap matrix
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   cPkL  :PkL expansion eigenvector at site ia.
Co Outputs
Co   wk    :(ppi-evl*sig)*evec
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ia,kmax,lcplxp,nkaph,nlma,nlmha,nlmto,iprmb(*)
      double precision
     .evl,sighp(nkaph,0:kmax,0:*),ppihp(nkaph,0:kmax,nlmha,nlma),
     .sigpp(0:kmax,0:kmax,0:*),ppipp(0:kmax,0:kmax,nlma,nlma)
      double complex wk(0:kmax,nlma),evec(1),cPkL(0:kmax,nlma)
      double complex ppihpz(nkaph,0:kmax,nlmha,nlma)
      double complex ppippz(0:kmax,0:kmax,nlma,nlma)
C ... Local parameters
      integer i,ilm,ilma,ilmb,io,jlm,k,k1,k2,l,ll,nlm1,nlm2,n0,nkap0,ik,
     .norb
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .blks(n0*nkap0),ntab(n0*nkap0)
      double precision xx

      call dpzero(wk, 2*(kmax+1)*nlma)

C ... Add Hsm*Pkl block of ppi-evl*sig times evec
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block evl*ppi contribution in groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
      do  io = 1, norb
C       l,ik = l and kaph indices, needed for sigma
        l  = ltab(io)
        ik = ktab(io)
C       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
        nlm1 = l**2+1
        nlm2 = (l+1)**2
        i  = offl(io)
C   ... evl*sig contribution requires explicit knowledge of l
        do  ilmb = nlm1, nlm2
          i = i+1
          do  k = 0, kmax
            wk(k,ilmb) = wk(k,ilmb) - evl*sighp(ik,k,l)*evec(i)
          enddo
        enddo
C   ... ppi contribution: loop over largest blocks possible
        if (blks(io) .ne. 0 .and. lcplxp .eq. 0) then
          nlm2 = nlm1 + blks(io)-1
          i  = offl(io)
          do  ilmb = nlm1, nlm2
            i = i+1
            do  ilma = 1, nlma
              do  k = 0, kmax
                wk(k,ilma) = wk(k,ilma) + ppihp(ik,k,ilmb,ilma)*evec(i)
              enddo
            enddo
          enddo
        endif
        if (blks(io) .ne. 0 .and. lcplxp .eq. 1) then
          nlm2 = nlm1 + blks(io)-1
          i  = offl(io)
          do  ilmb = nlm1, nlm2
            i = i+1
            do  ilma = 1, nlma
              do  k = 0, kmax
                wk(k,ilma) = wk(k,ilma) + ppihpz(ik,k,ilmb,ilma)*evec(i)
              enddo
            enddo
          enddo
        endif
      enddo

C ... Add Pkl*Pkl block of ppi-evl*sig times cPkL
      if (lcplxp .eq. 0) then
        do  ilm = 1, nlma
          l = ll(ilm)
          do  k1 = 0, kmax
            do  jlm = 1, nlma
              do  k2 = 0, kmax
                wk(k1,ilm) = wk(k1,ilm) +ppipp(k1,k2,ilm,jlm)*cPkL(k2,jlm)
              enddo
            enddo
            do  k2 = 0, kmax
              wk(k1,ilm) = wk(k1,ilm) - evl*sigpp(k1,k2,l)*cPkL(k2,ilm)
            enddo
          enddo
        enddo

      else
        do  ilm = 1, nlma
          l = ll(ilm)
          do  k1 = 0, kmax
            do  jlm = 1, nlma
              do  k2 = 0, kmax
                wk(k1,ilm) = wk(k1,ilm)+ppippz(k1,k2,ilm,jlm)*cPkL(k2,jlm)
              enddo
            enddo
            do  k2 = 0, kmax
              wk(k1,ilm) = wk(k1,ilm) - evl*sigpp(k1,k2,l)*cPkL(k2,ilm)
            enddo
          enddo
        enddo
      endif


      end

