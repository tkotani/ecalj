subroutine fsmbl(vavg,q,ndimh,nlmto, &
  nevec,evl,evec,ewgt,f)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,lhh,nkaphh,ispec,nbas,sspec=>v_sspec
  use m_uspecb,only:uspecb
  use m_struc_def
  use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
  use m_smhankel,only: hhigbl
  use m_lattic,only: rv_a_opos
  !- Force from smoothed hamiltonian (constant potential) and overlap
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
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
  integer :: ndimh,nlmto,nevec
  real(8):: q(3)
  double precision :: evl(ndimh),f(3,nbas),ewgt(nevec),vavg
  double complex evec(ndimh,ndimh)
  integer :: nlms,k0,n0,nkap0
  parameter (nlms=25, k0=1, n0=10, nkap0=3)
  integer:: i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,iq,is1,is2,l1,l2,ik1,ik2,ivec,m,nlm1,nlm2
  integer :: lh1(nkap0),lh2(nkap0),nkap1,nkap2,nlm21,nlm22,nlm11,nlm12
  integer:: blks1(n0*nkap0),ntab1(n0*nkap0)
  integer:: blks2(n0*nkap0),ntab2(n0*nkap0)
  real(8) :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3),wt,e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx(n0)
  complex(8):: s(nlms,nlms,0:k0,nkap0,nkap0), ds(nlms,nlms,0:k0,3,nkap0,nkap0),ccc(3),ssum(3)
  integer:: iloopend,ibl1,ibl2,ol1,oi1,ol2,oi2
  if (nevec <= 0) return
  call tcn ('fsmbl')
  do ib1=1,nbas
     is1=ispec(ib1) !ssite(ib1)%spec
     p1 =rv_a_opos(:,ib1) !ssite(ib1)%pos
     call uspecb(is1,rsm1,e1)
     call orblib1(ib1) !norb1,ltab1,ktab1,offl1
     call gtbsl1(4+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
     do ib2=ib1+1,nbas
        is2=ispec(ib2) !ssite(ib2)%spec
        p2 =rv_a_opos(:,ib2) !ssite(ib2)%pos
        call uspecb(is2,rsm2,e2)
        call orblib2(ib2)
        call gtbsl1(4+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)
        !     ... M.E. <1> and <KE> between all envelopes connecting ib1 and ib2
        do i1 = 1, nkaphh(is1) !nkap1
           do i2 = 1, nkaphh(is2) !nkap2
              nlm1 = (lhh(i1,is1)+1)**2
              nlm2 = (lhh(i2,is2)+1)**2
              if (nlm1 > nlms .OR. nlm2 > nlms) call rx('fsmbl: increase nlms')
              call hhigbl(11,p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),e2(1,i2),nlm1,nlm2,1,nlms, &
                   nlms , k0 ,& !rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy , &
                   s(1,1,0,i1,i2), ds(1,1,0,1,i1,i2) )
           enddo
        enddo
        do io2= 1, norb2
           ik2= ktab2(io2)
           ol2= ltab2(io2)**2
           oi2= offl2(io2)
           do io1= 1, norb1
              ik1= ktab1(io1)
              ol1= ltab1(io1)**2
              oi1= offl1(io1)
              do  ivec = 1, nevec
                 ssum = 0d0
                 do ibl2=1,blks2(io2) 
                    do ibl1=1,blks1(io1) 
                       ccc = [(vavg*ds(ol1+ibl1,ol2+ibl2,0,m,ik1,ik2) &
                            -       ds(ol1+ibl1,ol2+ibl2,1,m,ik1,ik2) &
                            - evl(ivec)*ds(ol1+ibl1,ol2+ibl2,0,m,ik1,ik2),m=1,3)]
                       ssum = ssum + dconjg(evec(oi1+ibl1,ivec))*ccc*evec(oi2+ibl2,ivec)
                    enddo
                 enddo
                 wt = ewgt(ivec)
                 f(:,ib1) = f(:,ib1) - 2*wt*ssum
                 f(:,ib2) = f(:,ib2) + 2*wt*ssum
              enddo
           enddo
        enddo
     enddo
  enddo
  call tcx ('fsmbl')
end subroutine fsmbl
