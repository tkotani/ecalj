subroutine flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto, &
     ndimh,iprmb,isp,evl,evec,ewgt,cPkL,db, sigpp,sighp,ppippz,ppihpz,f)
  use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !- Force contribution from augmentation at site ia.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ia    :site of augmentation
  !i   kmax  :polynomial cutoff
  !i   nkaph :leading dimension of sigph,ppihp at site ia
  !i   nlmha :dimensions ppihp
  !i   nlma  :augmentation L-cutoff
  !i   lmxa  :augmentation L-cutoff
  !i   ndimh :dimension of hamiltonian
  !i   evl   :eigenvalue
  !i   evec  :eigenvector
  !i   ewgt  :eigenvector weight
  !i   cPkL  :PkL expansion eigenvector at site ia.
  !i   b     :structure constants, needed for PW contribution
  !i   db    :gradient of structure constants, needed to make grad psi
  !i   da    :work array holding grad psi
  !i   wk    :work array holding (ppi-evl*sig)*evec
  !i   ppipp :local tail-tail potential matrix
  !i   ppippz:local tail-tail potential matrix, complex form
  !i         :NB: only ppipp or ppippz is used, depending on lcplxp
  !i   sigpp :local tail-tail overlap matrix
  !i   ppihp :local head-tail potential matrix
  !i   ppihpz:local head-tail potential matrix, complex form
  !i         :NB: only ppihp or ppihpz is used, depending on lcplxp
  !i   sighp :local head-tail overlap matrix
  !i   lcplxp=1 only now: ppi is complex
  !o Outputs
  !o   f     :local contribution to forces is added
  !r Remarks
  !u Updates
  !u   05 Jul 08 (T. Kotani) Contribution from PW part of basis
  !u    1 Sep 04 Adapted to handle complex ppi
  !u   17 Jun 00 spin polarized
  !u   25 May 00 Adapted from nfp floc_q.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ia,kmax,lmxa,nkaph,nbas,nlmto,ndimh,nlma,lmxha,nlmha,isp,iprmb(1)
  double precision :: evl,f(3,nbas),ewgt, &
       sigpp(0:kmax,0:kmax,0:lmxa,isp),sighp(nkaph,0:kmax,0:lmxha,isp)
  !     .ppihp(nkaph,0:kmax,nlmha,nlma,isp)
  !     .ppipp(0:kmax,0:kmax,nlma,nlma,isp)
  double complex db(ndimh,nlma,0:kmax,3),da(0:kmax,nlma,3), &
       evec(ndimh),cPkL(0:kmax,nlma),wk(0:kmax,nlma)
  double complex ppihpz(nkaph,0:kmax,nlmha,nlma,isp)
  double complex ppippz(0:kmax,0:kmax,nlma,nlma,isp)
  ! ... Local parameters
  integer :: i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,nlm2,n0,nkap0
  parameter (n0=10,nkap0=3)
  !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
  integer :: blks(n0*nkap0),ntab(n0*nkap0)
  double precision :: sum,wt,xx
  if (nlmto == 0) return

  call tcn('flocbl')
  ! ... Make (ppi-evl*sig)*psi in wk
  call flocb2(ia,nlmto,iprmb,kmax,nkaph,nlmha,nlma,evl,evec, &
       ppippz(0,0,1,1,isp),sigpp(0,0,0,isp), &
       ppihpz(1,0,1,1,isp),sighp(1,0,0,isp), &
       cPkL,wk) !lcplxp,

  ! ... Loop over ib, virtual shift of wavefct part centered there
  do  ib = 1, nbas
     if (ib == ia) goto 10
     !       Block orbitals into groups with consecutive l
     !        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
     call orblib(ib)
     call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

     !   ... Grad of psi expansion coeffs from a virtual shift at site ib
     call dpzero(da, 2*(kmax+1)*nlma*3)
     do  io = 1, norb
        if (blks(io) /= 0) then
           !         l2 = l index, needed to locate block in da,db
           l2   = ltab(io)
           nlm1 = l2**2+1
           nlm2 = nlm1 + blks(io)-1
           i    = offl(io)
           !         i = index to head term, (ilm,k) to augmentation
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
     !   --- Force term is (grad psi_kL) * (ppi-evl*sig)*evec ---
     do  m = 1, 3
        sum = 0d0
        do  ilm = 1, nlma
           do  k = 0, kmax
              sum = sum + 2d0*dconjg(da(k,ilm,m))*wk(k,ilm)
           enddo
        enddo
        wt = ewgt
        f(m,ib) = f(m,ib) - wt*sum
        f(m,ia) = f(m,ia) + wt*sum
     enddo
10   continue
  enddo

  ! --- Force at site ia from PWs ---
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
  ! ... Force term is (grad psi_kL) * (ppi-evl*sig)*evec
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
end subroutine flocbl

subroutine flocb2(ia,nlmto,iprmb,kmax,nkaph,nlmha,nlma,evl,evec, &
     !     .ppipp,ppippz,sigpp,ppihp,ppihpz,sighp,cPkL,wk)!lcplxp,
     ppippz,sigpp,ppihpz,sighp,cPkL,wk)!lcplxp,
  use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !- Make (ppi-evl*sig)*evec
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ia    :site of augmentation
  !l   nlmto :dimension of LMTO part of hamiltonian, for hh and ht blocks
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !i   nlmha :dimensions ppihp
  !i   nlma  :augmentation L-cutoff
  !i   kmax  :polynomial cutoff
  !i   evl   :eigenvalue
  !i   evec  :eigenvector
  !i   nkaph :leading dimension of sigph,ppihp
  !i   ppipp :local tail-tail potential matrix
  !i   ppippz:local tail-tail potential matrix, complex form
  !i         :NB: only ppipp or ppippz is used, depending on lcplxp
  !i   sigpp :local tail-tail overlap matrix
  !i   ppihp :local head-tail potential matrix
  !i   ppihpz:local head-tail potential matrix, complex form
  !i         :NB: only ppihp or ppihpz is used, depending on lcplxp
  !i   sighp :local head-tail overlap matrix
  !i   lcplxp=1 only now: ppi is complex
  !i   cPkL  :PkL expansion eigenvector at site ia.
  !o Outputs
  !o   wk    :(ppi-evl*sig)*evec
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ia,kmax,nkaph,nlma,nlmha,nlmto,iprmb(*) !,lcplxp
  double precision :: &
       evl,sighp(nkaph,0:kmax,0:*),& ! & ppihp(nkaph,0:kmax,nlmha,nlma),
  sigpp(0:kmax,0:kmax,0:*)!,ppipp(0:kmax,0:kmax,nlma,nlma)
  double complex wk(0:kmax,nlma),evec(1),cPkL(0:kmax,nlma)
  double complex ppihpz(nkaph,0:kmax,nlmha,nlma)
  double complex ppippz(0:kmax,0:kmax,nlma,nlma)
  ! ... Local parameters
  integer :: i,ilm,ilma,ilmb,io,jlm,k,k1,k2,l,ll,nlm1,nlm2,n0,nkap0,ik
  !     .norb
  parameter (n0=10,nkap0=3)
  integer :: blks(n0*nkap0),ntab(n0*nkap0)!ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
  double precision :: xx
  wk=0d0
  ! ... Add Hsm*Pkl block of ppi-evl*sig times evec
  !      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
  call orblib(ia)
  !     Block evl*ppi contribution in groups of consecutive l
  call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
  do  io = 1, norb
     !       l,ik = l and kaph indices, needed for sigma
     l  = ltab(io)
     ik = ktab(io)
     !       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
     nlm1 = l**2+1
     nlm2 = (l+1)**2
     i  = offl(io)
     !   ... evl*sig contribution requires explicit knowledge of l
     do  ilmb = nlm1, nlm2
        i = i+1
        do  k = 0, kmax
           wk(k,ilmb) = wk(k,ilmb) - evl*sighp(ik,k,l)*evec(i)
        enddo
     enddo
     !   ... ppi contribution: loop over largest blocks possible
     if (blks(io) /= 0) then ! .AND. lcplxp == 1) then
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
  ! ... Add Pkl*Pkl block of ppi-evl*sig times cPkL
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
end subroutine flocb2

