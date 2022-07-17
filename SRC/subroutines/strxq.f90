subroutine strxq(mode,e,q,p,nlma,nlmh,ndim,alat,vol,awald,nkd,nkq, &
     dlv,qlv,cg,indxcg,jcg,s,sd)
  use m_lldata,only: ll
  !use m_shortn4,only: shortn4,shortn4_initialize,nout,nlatout
  use m_hamindex,only:   plat,qlat
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  use m_hsmq,only: hsmq,hsmqe0
  !- One-center expansion coefficents to j of Bloch summed h (strux)
  ! ----------------------------------------------------------------
  !r  Onsite contribution is not contained in the bloch sum in the case of p=0.
  !r  See job=10 for hsmq.
  !i Inputs:
  !i   mode  :1's digit (not implemented)
  !i         :1: calculate s only
  !i         :2: calculate sd only
  !i         :any other number: calculate both s and sdot
  !i   e     :energy of Hankel function.  e must be <=0
  !i   q     :Bloch wave number
  !i   p     :position of Hankel function center;
  !i         :structure constants are for expansion about the origin
  !i   nlma  :Generate coefficients S_R'L',RL for L' < nlma
  !i   nlmh  :Generate coefficients S_R'L',RL for L  < nlmh
  !i   ndim  :leading dimension of s,sdot
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   vol   :cell volume
  !i   awald :Ewald smoothing parameter
  !i   nkq   :number of direct-space lattice vectors
  !i   nkq   :number of reciprocal-space lattice vectors
  !i   dlv   :direct-space lattice vectors, units of alat
  !i   qlv   :reciprocal lattice vectors, units of 2pi/alat
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L quantum number for the C.G. coefficients (scg.f)
  !o Outputs
  !o   s      :structure constant matrix S_R'L',RL
  !o   sd     :Energy derivative of s
  ! ----------------------------------------------------------------
  !1.  Bloch phase.  For translation vectors T, it's sum_T exp(+i q T)

  !2.  Methfessel's definitions of Hankels and Bessel functions:

  !    h_0 = Re e^(ikr)/r and j = sin(kr)/kr, k = sqrt(e), Im k >=0.
  !    H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)

  !    They are related to the usual n_l and j_l by factors (I think)
  !       H_L  =  (i k)^(l+1) n_l (kr) Y_L   (E < 0)
  !       J_L  =  (i k)^(-l)  j_l (kr) Y_L   (E < 0)

  !   which explains how the energy-dependence is extracted out.
  !   Also cases e .ne. 0 and e .eq. 0 have the same equations.

  !r Expansion Theorem: H_{RL}(r) = H_L(r-R)
  !r   H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
  !r   S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1)^l (-E)^(l+l'-l")/2 H_L"(E,R-R')
  ! ---
  implicit none
  integer :: mode,ndim,nlma,nlmh
  integer :: indxcg(*),jcg(*),nkd,nkq
  double precision :: p(3),q(3),alat,awald,vol,e,cg(*), dlv(*),qlv(*)
  double complex s(ndim,nlmh),sd(ndim,nlmh)
  integer :: lmxx,nrxmx,nlm0
  double precision :: fpi,p1(3),sp
  real(8),allocatable :: wk(:),yl(:),efac(:)
  complex(8),allocatable :: dl(:),dlp(:)
  integer(4),allocatable :: sig(:)
  double complex phase,sumx,sud !dl(nlm0),dlp(nlm0)
  integer :: icg,icg1,icg2,ii,indx,ipow,l,lmax,nrx,nlm, ilm,ilma,la,ilmb,lh !sig(0:lmxx),
  logical :: ldot
  integer(4) :: job
  integer ::lmax_(1)
  real(8):: e_(1),rsm_(1),pp(3)
  ldot = .false.
  lmax = ll(nlma)+ll(nlmh)
  nlm = (lmax+1)**2
  nrx  = max(nkd,nkq)
  fpi  = 16d0*datan(1d0)
  if (nlma > ndim) call rxi('strxq: increase ndim: need',nlma)
  lmxx = lmax
  nlm0 =(lmxx+1)**2
  nrxmx= nrx
  allocate( wk((lmxx*2+10)*nrxmx),yl(nrxmx*(lmxx+1)**2), &
       efac(0:lmxx),sig(0:lmxx),dl(nlm0),dlp(nlm0))
  pp= matmul(transpose(qlat),p)
  call shortn3_plat(pp)
  p1 = matmul(plat,pp+nlatout(:,1))
  sp = fpi/2*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  job = 0
  if( sum(abs(p))<1d-10 ) job = 10
  lmax_(1)=lmax
  e_(1)=e
  rsm_(1)=0d0
  if (e < 0) then
     call hsmq(1,0,lmax_,e_,rsm_,job,q,p1,nrx,nlm0,yl,awald,alat,qlv,nkq,dlv,nkd,vol,dl,dlp)
  else
     call hsmqe0(lmax,0d0,job,q,p1,nrx,nlm0,wk,yl, awald,alat,qlv,nkq,dlv,nkd,vol,dl)
     ldot = .false.
  endif
  if (sp /= 0d0) then
     do  20  ilm = 1, nlm
        dl(ilm) = phase*dl(ilm) ! ... Put in phase to undo shortening
        if (ldot) dlp(ilm) = phase*dl(ilm)
20   enddo
  endif
  ! --- Combine with Clebsch-Gordan coefficients ---
  ! ... efac(l)=(-e)**l; sig(l)=(-)**l
  efac(0) = 1
  sig(0) = 1
  do  l = 1, lmax
     efac(l) = -e*efac(l-1)
     sig(l) = -sig(l-1)
  enddo
  do  11  ilma = 1, nlma
     la = ll(ilma)
     do  14  ilmb = 1, nlmh
        lh = ll(ilmb)
        ii = max0(ilma,ilmb)
        indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        sumx = 0d0
        sud = 0d0
        if (ldot) then
           do  16  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx = sumx + cg(icg)*efac(ipow)*dl(ilm)
              sud = sud + cg(icg)*efac(ipow)*(dlp(ilm)+ipow*dl(ilm)/e)
16         enddo
        else
           do  15  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx  = sumx + cg(icg)*efac(ipow)*dl(ilm)
15         enddo
        endif
        s(ilma,ilmb) = fpi*sig(lh)*dconjg(sumx)
        if (ldot) sd(ilma,ilmb) = fpi*dconjg(sud)*sig(lh)
14   enddo
11 enddo
  if (allocated(wk))deallocate(wk)
  if (allocated(yl))deallocate(yl)
  if (allocated(efac))deallocate(efac)
  if (allocated(sig))deallocate(sig)
  if (allocated(dl))deallocate(dl)
  if (allocated(dlp))deallocate(dlp)
end subroutine strxq
