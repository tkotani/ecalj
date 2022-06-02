subroutine mkorbm(isp,nev,iq,qp,evec, orbtm)
  use m_lmfinit,only: ssite=>v_ssite,sspec=>v_sspec,nbas,nlmax,nsp,nspc,nl,n0,nppn,nab
  use m_igv2x,only: napw,ndimh,ndimhx,igvapw=>igv2x
  use m_mkpot,only: ppnl=>ppnl_rv
  use m_subzi, only: wtkb=>rv_a_owtkb
  use m_qplist,only: nkp
  use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx
  !- Decomposition of norm from projection of w.f into MT sphere
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored:
  !i     Passed to:
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa rmt
  !i     Stored:
  !i     Passed to:
  !i   isp   :current spin channel (1 or 2)
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   nlmax :leading dimension of aus
  !i   ndham :dimensions aus,wtkp
  !i   nev   :number of eigenvectors to accumulate orbital moment
  !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
  !i   iq    :current k-point
  !i   nbas  :size of basis
  !i   ppnl  :NMTO potential parameters; see eg potpus.f
  !i   aus   :values of (phi,phidot) MT sphere boundary; see makusq
  !i   nl    :(global maximum l) + 1
  !i   nkp   :number of irreducible k-points (bzmesh.f)
  !o Outputs
  !o   orbtm :orbital moments accumulated for this qp
  !l Local variables
  !l   ispc  :the current spin index in the coupled spins case.
  !l         :Some quantities have no separate address space for each
  !l         :spin in the indepedent-spins case (evec,evl,ewgt) but do
  !l         :in the coupled-spins case.  A separate loop ispc=1..nspc
  !l         :must be added for the latter case
  !l         :ispc is the appropriate index for objects which distinguish
  !l         :spins in the spin-coupled case only
  !l   isp   :isp  is the appropriate index for objects which distinguish
  !l         :spins in the spin-uncoupled case only
  !l   ksp   :the current spin index in both independent and coupled
  !l         :spins cases.
  !l         :ksp is appropriate spin index for quantities that have
  !l         :separate address space for each spin in every case
  !l         :(potential- and density-like objects).
  !u Updates
  !u   25 Apr 05 (A. Chantis) extended to local orbitals
  !u   24 Dec 04 Extended to spin-coupled case
  !u   30 Aug 04 (A. Chantis) first written, adapted from mkpdos
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: isp,nev
  integer :: iq
  !      parameter (n0=10,nppn=12,nab=9)
  !      real(8):: wtkp(ndham,nsp,nkp)
  !      type(s_site)::ssite(*)
  !      type(s_spec)::sspec(*)

  double precision :: orbtm(nl,nsp,*) !ppnl(nppn,n0,nsp,*),
  !      double complex aus(nlmax,ndham*nspc,3,nsp,*)
  ! ... Local parameters
  integer :: lmxa,lmxax,lmdim,ichan,ib,is,igetss,iv,ilm,l,m,nlma, &
       ll,lc,em,ispc,ksp
  double precision :: suml(11),s11,s22,s12,s33,s31,s32,s13,s23, &
       suma,rmt,sab(nab,n0,2)
  double complex au,as,az,iot
  complex(8),allocatable ::aus(:,:,:,:,:)
  double complex evec(ndimh,nsp,nev)
  real(8):: qp(3)
  allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))
  aus=0d0
  call makusq(0 , nbas,0, nev, isp,1,qp,evec, aus )
  lmxax = ll(nlmax)
  iot = dcmplx(0d0,1d0)
  ichan = 0
  do  ib = 1, nbas
     is = int(ssite(ib)%spec)


     lmxa=sspec(is)%lmxa
     rmt=sspec(is)%rmt

     lmxa = min(lmxa,lmxax)
     if (lmxa == -1) goto 10

     nlma = (lmxa+1)**2
     lmdim = nlma

     call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)

     !       In noncollinear case, isp=1 always => need internal ispc=1..2
     !       ksp is the current spin index in both cases:
     !       ksp = isp  in the collinear case
     !           = ispc in the noncollinear case
     !       ispc=1 for independent spins, and spin index when nspc=2
     do  ispc = 1, nspc
        ksp = max(ispc,isp)
        do  iv = 1, nev
           call dpzero(suml,1+lmxa)
           suma = 0
           ilm = 0

           !  ....  Rotate from real to spherical harmonics (order assumed: m,...,-m).
           !        |Psi>_l = \Sum_{m}(A_l,m * u_l + B_l,m * s_l)*R_l,m --->
           !        |Psi>_l = \Sum_{m}(C_l,m * u_l + D_l,m * s_l)*Y_l,m
           !        R_l,m and Y_l,m are the real and spherical harmonics respectively.

           !              | (-1)^m/sqrt(2)*A_l,-m + i*(-1)^m/sqrt(2)*A_l,m , m>0
           !        C_l,m=|  A_l,m                                         , m=0
           !              |  1/sqrt(2)*A_l,-m -  i*1/sqrt(2)*A_l,m         , m<0

           !       Same relationships are valid between D and B.

           do  l = 0, lmxa
              lc = (l+1)**2 - l
              do  m = -l, l
                 em = abs(m)
                 ilm = ilm+1
                 !    ...    m,...,-m order
                 !            if (m .lt. 0) then
                 !            au = 1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) -
                 !     .           iot/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
                 !            as = 1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) -
                 !     .           iot/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
                 !            az = 1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) -
                 !     .           iot/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
                 !            else if (m .gt. 0) then
                 !            au = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) +
                 !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
                 !            as = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) +
                 !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
                 !            az = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) +
                 !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
                 !            else
                 !    ...   -m,...,m order
                 if (m < 0) then
                    au = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) + &
                         1d0/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
                    as = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) + &
                         1d0/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
                    az = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) + &
                         1d0/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
                 else if (m > 0) then
                    au = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) + &
                         (-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
                    as = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) + &
                         (-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
                    az = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) + &
                         (-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
                 else
                    au = aus(ilm,iv,1,ksp,ib)
                    as = aus(ilm,iv,2,ksp,ib)
                    az = aus(ilm,iv,3,ksp,ib)
                 endif

                 !           If (au,as) are coefficients to (u,s), use this
                 s11 = dconjg(au)*au*sab(1,l+1,ksp)
                 s12 = 2*dconjg(au)*as*sab(2,l+1,ksp)
                 s22 = dconjg(as)*as*sab(4,l+1,ksp)
                 s33 = dconjg(az)*az*sab(5,l+1,ksp)
                 s31 = dconjg(az)*au*sab(6,l+1,ksp)
                 s32 = dconjg(az)*as*sab(7,l+1,ksp)
                 s13 = dconjg(au)*az*sab(6,l+1,ksp)
                 s23 = dconjg(as)*az*sab(7,l+1,ksp)

                 orbtm(l+1,ksp,ib)=orbtm(l+1,ksp,ib)+m*(s11+s12+s22+ &
                      s33+s32+s23+s31+s13)*wtkb(iv,isp,iq) !corrected. it should be wtkb

              enddo
           enddo

        enddo
     enddo
     !        print*, isp, ib, 'ORB.MOMNT=',orbtm(isp,ib)
10   continue
  enddo
  deallocate(aus)
end subroutine mkorbm











