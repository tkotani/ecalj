subroutine fsmbpw(nbas,ssite,sspec,vavg,ndimh,nlmto, &
     nevec,evl,evec,ewgt,napw,qpgv,qpg2v,ylv,nlmax,lmxax,alat,sqv,f)
  use m_struc_def
  use m_uspecb,only:uspecb
  use m_orbl,only: Orblib1,ktab1,ltab1,offl1,norb1
  use m_lattic,only: rv_a_opos
  !- Force from smoothed hamiltonian (constant potential), PW contribution
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
  !i   vavg  :constant potential (MT zero) to be added to h
  !i   ndimh :dimension of hamiltonian
  !i   nlmto :dimension of lmto part of hamiltonian
  !i   nevec :number of occupied eigenvectors
  !i   evl   :eigenvalues
  !i   evec  :eigenvectors
  !i   ewgt  :eigenvector weights
  !i   napw  :number of augmented PWs in basis
  !i   qpgv
  !i   qpg2v
  !i   ylv
  !i   nlmax
  !i   lmxax
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   sqv   :square root of volume
  !o Outputs
  !o   f     :PW contribution to force is added to f
  !r Remarks
  !u Updates
  !u   04 Jul 08 (T. Kotani) first created
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nbas,ndimh,napw,nlmax,nlmto,nevec,lmxax
  integer,parameter:: n0=10, nkap0=3
  real(8),parameter:: pi = 4d0*datan(1d0), fpi = 4*pi
  real(8):: evl(ndimh),f(3,nbas),ewgt(nevec),vavg,qpgv(3,napw),qpg2v(napw),qpg2,alat,sqv
  real(8):: gam,denom, e1(n0,nkap0),rsm1(n0,nkap0),p1(3),xx(n0),wt,ylv(napw,nlmax),ssum(3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  integer :: i1,i2,ib1,ilm1,io1,iq,is1,l1,ik1,ig,ivec,nglob, nlm11,nlm12,m
  integer:: blks1(n0*nkap0),ntab1(n0*nkap0),lh1(nkap0),nkap1
  complex(8):: phase,fach,ovl,ccc(3),sum, srm1l(0:n0),evec(ndimh,ndimh),img=(0d0,1d0)
  integer:: ibl1,oi1,ol1
  ! --- Setup ---
  if (nevec <= 0) return
  call tcn ('fsmbpw')
  srm1l(0)=1d0
  do  l1 = 1, lmxax
     srm1l(l1) = img**l1
  enddo
  ! --- Loop over first and second site indices ---
  do 1000 ib1=1,nbas
     is1=ssite(ib1)%spec
     p1=rv_a_opos(:,ib1) !ssite(ib1)%pos
     call uspecb(is1,rsm1,e1)
     call orblib1(ib1) !norb1,ltab1,ktab1,xx,offl1,xx)
     call gtbsl1(8+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
     !   ... Hsm (i1) \times i(q+G)[(q+G)**2+const] PW (i2) Takao. Taken from smhsbl.f
     !       i1--> Hsm, i2--> PW
     do 2000 ig = 1, napw
        i2 = ig + nlmto
        qpg2 = qpg2v(ig)
        phase = exp(img*alat*sum(qpgv(:,ig)*p1))
        do 3000 io1 = 1, norb1
           l1  = ltab1(io1)
           ik1 = ktab1(io1)
           ol1 = ltab1(io1)**2
           oi1 = offl1(io1)
           denom = e1(l1+1,ik1) - qpg2
           gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
           fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
           do 3010 ibl1 = 1,blks1(io1) 
              !             s(i1,i2) = ovl
              !             h(i1,i2) = (qpg2 + vavg) * ovl
              ovl = fach * ylv(ig,ol1+ibl1)/sqv ! Eq. 9.4 in JMP39 3393
              do ivec = 1, nevec  !        gradient PW * (H - E S)
                 ccc = [(ovl * img*qpgv(m,ig) * (qpg2 + vavg - evl(ivec)),m=1,3)]
                 ssum = ewgt(ivec)*[(dconjg(evec(oi1+ibl1,ivec))*ccc(m)*evec(i2,ivec),m=1,3)]
                 f(:,ib1) = f(:,ib1) - 2d0*ssum
              enddo
3010       enddo
3000    enddo
2000 enddo
1000 enddo
  call tcx ('fsmbpw')
end subroutine fsmbpw


