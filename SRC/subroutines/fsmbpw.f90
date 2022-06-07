subroutine fsmbpw(nbas,ssite,sspec,vavg,ndimh,nlmto, &
     nevec,evl,evec,ewgt,napw,qpgv,qpg2v,ylv,nlmax,lmxax,alat,sqv,f)
  use m_struc_def
  use m_uspecb,only:uspecb
  use m_orbl,only: Orblib1,ktab1,ltab1,offl1,norb1
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
  integer:: i_copy_size,iloop
  ! ... Passed parameters
  integer :: nbas,ndimh,napw,nlmax,nlmto,nevec,lmxax
  real(8):: ylv(napw,nlmax)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  double precision :: evl(ndimh),f(3,nbas),ewgt(nevec),vavg
  double complex evec(ndimh,ndimh)
  double precision :: qpgv(3,napw),qpg2v(napw),qpg2,alat,sqv
  ! ... Local parameters
  integer :: n0,nkap0
  parameter (n0=10, nkap0=3)
  integer :: i1,i2,ib1,ilm1,io1,iq,is1,l1,ik1,ig,ivec,nglob, &
       nlm11,nlm12,m
  !      integer norb1,ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),
  integer:: blks1(n0*nkap0),ntab1(n0*nkap0),lh1(nkap0),nkap1
  double precision :: gam,denom,pi,fpi,ddot
  double precision :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3),xx(n0),wt
  double complex phase,srm1,fach,ovl,ccc(3),sum,srm1l(0:n0)
  parameter (srm1=(0d0,1d0))
  ! ... Heap
  integer:: ibini,ibend

  ! --- Setup ---
  if (nevec <= 0) return
  call tcn ('fsmbpw')

  pi = 4d0*datan(1d0)
  fpi = 4*pi
  srm1l(0) = 1d0
  do  l1 = 1, lmxax
     srm1l(l1) = (srm1)**l1
  enddo

  ! --- Loop over first and second site indices ---
  ibini=1
  ibend=nbas

  do iloop = ibini,ibend
     ib1=iloop
     is1=ssite(ib1)%spec
     i_copy_size=size(ssite(ib1)%pos)
     call dcopy(i_copy_size,ssite(ib1)%pos,1,p1,1)

     call uspecb(is1,rsm1,e1)
     !       Row info telling fsmbpw where to poke s0 made by hhibl
     call orblib1(ib1) !norb1,ltab1,ktab1,xx,offl1,xx)

     !       For now, do not allow l blocks to be grouped.
     !       To do so will require rewriting loops below.
     call gtbsl1(8+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)

     !   ... Hsm (i1) \times i(q+G)[(q+G)**2+const] PW (i2) Takao. Taken from smhsbl.f
     !       i1--> Hsm, i2--> PW
     do  ig = 1, napw
        i2 = ig + nlmto
        qpg2 = qpg2v(ig)
        phase = exp(srm1*alat*ddot(3,qpgv(1,ig),1,p1,1))
        do  io1 = 1, norb1
           if (blks1(io1) /= 0) then
              !           l1,ik1 = l and kaph indices, needed to locate block in s0
              l1  = ltab1(io1)
              ik1 = ktab1(io1)
              nlm11 = l1**2+1
              nlm12 = nlm11 + blks1(io1)-1
              !           Note:  using srm1l => l must be fixed in ilm loop below
              denom = e1(l1+1,ik1) - qpg2
              gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
              fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
              i1 = offl1(io1)
              do  ilm1 = nlm11, nlm12
                 i1 = i1+1
                 !             s(i1,i2) = ovl
                 !             h(i1,i2) = (qpg2 + vavg) * ovl
                 ovl = fach * ylv(ig,ilm1)/sqv ! Eq. 9.4 in JMP39 3393

                 !         ... Loop over occupied eigenstates and x,y,z
                 do  m = 1, 3
                    sum = 0
                    do  ivec = 1, nevec
                       !             gradient PW * (H - E S)
                       ccc(m) = ovl * srm1*qpgv(m,ig) * (qpg2 + vavg - evl(ivec))
                       sum = dconjg(evec(i1,ivec))*ccc(m)*evec(i2,ivec)
                       wt = ewgt(ivec)
                       f(m,ib1) = f(m,ib1) - 2*wt*sum

                    enddo               !ivec
                 enddo               !m
              enddo                 !ilm1
           endif
        enddo                 !io1
     enddo                   !ig
  enddo                     !ib1
  call tcx ('fsmbpw')
end subroutine fsmbpw


