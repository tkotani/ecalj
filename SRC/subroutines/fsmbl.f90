subroutine fsmbl(nbas,ssite,sspec,vavg,q,ndimh,nlmto,iprmb, &! & ,slat
  nevec,evl,evec,ewgt,f)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,lhh,nkaphh
  use m_uspecb,only:uspecb
  use m_struc_def
  use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
!  use m_shankel,only: hhigbl
  !- Force from smoothed hamiltonian (constant potential) and overlap
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: *
  !i     Stored:    *
  !i     Passed to: uspecb
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: ocg ojcg oidxcg ocy
  !i     Stored:    *
  !i     Passed to: hhigbl
  !i   vavg  :constant potential (MT zero) to be added to h
  !i   q     :Bloch wave vector
  !i   ndimh :dimension of hamiltonian
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !i   nevec :number of occupied eigenvectors
  !i   evl   :eigenvalues
  !i   evec  :eigenvectors
  !i   ewgt  :eigenvector weights
  !o Outputs
  !o   f
  !r Remarks
  !u Updates
  !u   05 Jul 08 Decouple ndimh from nlmto, for PW basis
  !u   10 Apr 02 Redimensioned eh,rsmh to accommodate larger lmax
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   10 Sep 01 Extended to local orbitals.
  !u   23 May 00 Adapted from nfp fsm_q.f
  ! ----------------------------------------------------------------------
  implicit none
  integer:: iloop
  integer :: nbas,ndimh,nlmto,nevec,iprmb(nlmto)
  real(8):: q(3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double precision :: evl(ndimh),f(3,nbas),ewgt(nevec),vavg
  double complex evec(ndimh,ndimh)
  ! ... Local parameters
  integer :: nlms,k0,n0,nkap0
  parameter (nlms=25, k0=1, n0=10, nkap0=3)
  integer:: i1 , i2 , ib1 , ib2 , ilm1 , ilm2 , io1 , io2 , iq &
       , is1 , is2 , l1 , l2 , ik1 , ik2 , ivec , m , nglob , nlm1 , &
       nlm2
  integer :: lh1(nkap0),lh2(nkap0),nkap1,nkap2,nlm21,nlm22,nlm11,nlm12
  !      integer norb1,ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),
  integer:: blks1(n0*nkap0),ntab1(n0*nkap0)
  !      integer norb2,ltab2(n0*nkap0),ktab2(n0*nkap0),offl2(n0*nkap0),
  integer:: blks2(n0*nkap0),ntab2(n0*nkap0)
  double precision :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3),wt, &
       e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx(n0)
  double complex  s(nlms,nlms,0:k0,nkap0,nkap0)
  double complex ds(nlms,nlms,0:k0,3,nkap0,nkap0),ccc,sum
  ! ... Heap
  integer:: iloopend
  if (nevec <= 0) return
  call tcn ('fsmbl')
  iloopend=nbas
  do  iloop = 1, iloopend
     ib1=iloop
     is1=ssite(ib1)%spec
     p1 =ssite(ib1)%pos
     call uspecb(is1,rsm1,e1)
     !       Row info telling fsmbl where to poke s0 made by hhibl
     !        call orbl(ib1,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
     call orblib1(ib1) !,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
     call gtbsl1(4+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
     do  ib2 = ib1+1, nbas
        is2=ssite(ib2)%spec
        p2=ssite(ib2)%pos
        call uspecb(is2,rsm2,e2)
        !         Column info telling fsmbl where to poke s0 made by hhibl
        !          call orbl(ib2,0,nlmto,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
        call orblib2(ib2)
        call gtbsl1(4+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)
        !     ... M.E. <1> and <KE> between all envelopes connecting ib1 and ib2
        do  i1 = 1, nkaphh(is1) !nkap1
           do  i2 = 1, nkaphh(is2) !nkap2
              nlm1 = (lhh(i1,is1)+1)**2
              nlm2 = (lhh(i2,is2)+1)**2
              if (nlm1 > nlms .OR. nlm2 > nlms) &
                   call rx('fsmbl: increase nlms')
              call hhigbl ( 11 , p1 , p2 , q , rsm1 ( 1 , i1 ) , rsm2 ( 1 , &
                   i2 ) , e1 ( 1 , i1 ) , e2 ( 1 , i2 ) , nlm1 , nlm2 , 1 , nlms &
                   , nlms , k0 , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy , &
                   s ( 1 , 1 , 0 , i1 , i2 ) , ds ( 1 , 1 , 0 , 1 , i1 ,  i2 ) )
           enddo
        enddo

        !     ... Loop over pairs of orbital groups, multiply bloc of gradients
        do  io2 = 1, norb2
           if (blks2(io2) /= 0) then
              !           l2,ik2 = l and kaph indices, needed to locate block in s0
              l2  = ltab2(io2)
              ik2 = ktab2(io2)
              nlm21 = l2**2+1
              nlm22 = nlm21 + blks2(io2)-1
              do  io1 = 1, norb1
                 if (blks1(io1) /= 0) then
                    !             l1,ik1 = l and kaph indices, needed to locate block in s0
                    l1  = ltab1(io1)
                    ik1 = ktab1(io1)
                    nlm11 = l1**2+1
                    nlm12 = nlm11 + blks1(io1)-1
                    !         ... Loop over occupied eigenstates and x,y,z
                    do  ivec = 1, nevec
                       do  m = 1, 3
                          !           ... Loop over orbital pairs within the groups
                          sum = 0d0
                          !               i2 = hamiltonian offset
                          i2 = offl2(io2)
                          do  ilm2 = nlm21, nlm22
                             i2 = i2+1
                             !                 i1 = orbital index in iprmb order
                             i1 = offl1(io1)
                             do  ilm1 = nlm11, nlm12
                                i1 = i1+1
                                ccc = vavg*ds(ilm1,ilm2,0,m,ik1,ik2) &
                                     -      ds(ilm1,ilm2,1,m,ik1,ik2) &
                                     - evl(ivec)*ds(ilm1,ilm2,0,m,ik1,ik2)
                                sum = sum + dconjg(evec(i1,ivec))*ccc*evec(i2,ivec)
                             enddo
                          enddo
                          wt = ewgt(ivec)
                          f(m,ib1) = f(m,ib1) - 2*wt*sum
                          f(m,ib2) = f(m,ib2) + 2*wt*sum
                       enddo
                    enddo

                 endif
              enddo
           endif
        enddo

     enddo
  enddo
  call tcx ('fsmbl')
end subroutine fsmbl

