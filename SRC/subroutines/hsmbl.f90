subroutine hsmbl(p,rsm,e,q,lmax,cy,hsm,hsmp) !,slat
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,lat_plat
  use m_struc_def           !Cgetarg
  use m_lmfinit,only: lat_alat,lat_tol
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_lattic,only: lat_awald
  use m_lattic,only: lat_nkd
  use m_lattic,only: lat_nkq
  !- Bloch-sum of smooth Hankel functions and energy derivatives
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   lmax  :l-cutoff for hsm
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   hsm   :Bloch-summed smoothed Hankels
  !o   hsmp  :Energy derivatives of hsm
  !r Remarks
  !r  Hankel functions evaluated by Ewald summation.
  !r  p in units of alat, qlv in units of 2 pi / alat
  !r  See also hsmq for a vectorized version.
  !u Updates
  !u   26 Jan 07 Works with positive energies e
  !u   1 May 00 Adapted from nfp hsmbl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: lmax
  real(8):: rsm , e , q(3) , p(3) , cy(1)
  !      type(s_lat)::slat

  double complex hsm(1),hsmp(1)
  ! ... Local parameters
  integer:: nkd , nkq , ilm , l , m
  double precision :: alat,p1(3),plat(3,3),qlat(3,3),sp,pi,rwald,awald,asm, &
       tol,vol
  double complex cfac,phase
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  call shorbz(p,p1,plat,qlat)
  pi = 4d0*datan(1d0)
  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  awald=lat_awald
  tol=lat_tol
  vol=lat_vol
  nkd=lat_nkd
  nkq=lat_nkq
  rwald = 1d0/awald
  asm = 1d0/rsm
  if (rsm < rwald) then
     call hsmblq ( p1 , e , q , awald , lmax , alat , rv_a_oqlv , &
          nkq , vol , hsm , hsmp )
     call hsmbld ( p1 , rsm , e , q , awald , lmax , alat , rv_a_odlv &
          , nkd , hsm , hsmp )
  else
     call hsmblq ( p1 , e , q , asm , lmax , alat , rv_a_oqlv , nkq &
          , vol , hsm , hsmp )
  endif
  ! ... Multiply by phase to undo shortening
  cfac = dcmplx(0d0,1d0)*phase
  ilm = 0
  do    l = 0, lmax
     cfac = cfac*dcmplx(0d0,-1d0)
     do    m = 1, 2*l+1
        ilm = ilm+1
        hsm(ilm)  = cfac*cy(ilm)*hsm(ilm)
        hsmp(ilm) = cfac*cy(ilm)*hsmp(ilm)
     enddo
  enddo
end subroutine hsmbl


subroutine hsmblq(p,e,q,a,lmax,alat,qlv,nkq,vol,dl,dlp)
  !- k-space part of smooth hankel bloch sum
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p (units of alat)
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   a     :Ewald smoothing parameter
  !i   lmax  :l-cutoff for dl,dlp
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   qlv   :reciprocal lattice vectors, units of 2pi/alat
  !i   nkq   :number of r.l.v.
  !i   vol   :cell volume
  !o Outputs
  !i   dl    :k-summed smoothed Bloch hankels
  !i   dlp   :energy derivative of dl
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nkq,lmax
  double precision :: a,alat,e,vol
  double precision :: q(3),p(3),qlv(3,nkq)
  double complex dl(1),dlp(1)
  ! ... Local parameters
  integer :: lmxx,nlm,ilm,ir
  parameter (lmxx=11)
  double precision :: r(3),tpi,gamma,fpibv,tpiba,scalp,r2,den0,den1
  double precision :: yl((lmxx+1)**2)
  double complex eiphi

  if (lmax > lmxx) call rxi('hsmblq: increase lmxx to',lmax)
  tpi = 8d0*datan(1d0)
  gamma = .25d0/(a*a)
  fpibv = 2d0*tpi/vol
  tpiba = tpi/alat
  nlm = (lmax+1)**2
  do  10  ilm = 1, nlm
     dl(ilm) = dcmplx(0d0,0d0)
     dlp(ilm) = dcmplx(0d0,0d0)
10 enddo
  do  26  ir = 1, nkq
     r(1) = tpiba*(q(1)+qlv(1,ir))
     r(2) = tpiba*(q(2)+qlv(2,ir))
     r(3) = tpiba*(q(3)+qlv(3,ir))
     scalp = alat*(r(1)*p(1)+r(2)*p(2)+r(3)*p(3))
     eiphi = dcmplx(dcos(scalp),dsin(scalp))
     call sylm(r,yl,lmax,r2)
     den0 = dexp(gamma*(e-r2))/(r2-e)
     den1 = den0/(r2-e)
     do  ilm = 1, nlm
        dl(ilm)  = dl(ilm) + yl(ilm)*eiphi*den0
        dlp(ilm) = dlp(ilm) + yl(ilm)*eiphi*den1
     enddo
26 enddo
  do  35  ilm = 1, nlm
     dl(ilm)  = fpibv*dl(ilm)
     dlp(ilm) = fpibv*dlp(ilm) + gamma*dl(ilm)
35 enddo
end subroutine hsmblq


subroutine hsmbld(p,rsm,e,q,a,lmax,alat,dlv,nkd,dl,dlp)
  !- Adds real space part of reduced structure constants (ewald).
  !u Updates
  !u   10 May 07 New protections to handle error functions underflow
  !     implicit none
  ! ... Passed parameters
  integer :: lmxx,ilm,ir,l,lmax,m,nkd,nm
  double precision :: q(3),p(3),dlv(3,nkd)
  double complex dl(1),dlp(1)
  ! ... Local parameters
  logical :: lpos,lzero
  parameter (lmxx=11)
  double precision :: a,a2,akap,alat,asm,asm2,cc,ccsm,derfc,e,emkr,gl, &
       qdotr,r1,r2,rsm,srpi,ta,ta2,tasm,tasm2,umins,uplus,tpi,kap
  double precision :: yl((lmxx+1)**2),chi1(-1:10),chi2(-1:10),r(3)
  !     double complex zchi1(-1:10),zchi2(-1:10)
  double complex cfac,zikap,expikr,zuplus,zerfc
  double precision :: srfmax,fmax
  parameter (srfmax=16d0, fmax=srfmax*srfmax)

  if (lmax > lmxx) call rxi('hsmbld: increase lmxx to',lmax)
  tpi = 8.d0*datan(1.d0)
  srpi = dsqrt(tpi/2.d0)
  lpos = e .gt. 0d0
  if (lpos) then
     !     Methfessl's 'akap' is sqrt(-e) = i kap by standard kap=sqrt(e)
     kap = dsqrt(e)
     zikap = dcmplx(0d0,1d0)*kap
  else
     akap = dsqrt(-e)
  endif
  ta = 2.d0*a
  a2 = a*a
  ta2 = 2.d0*a2
  cc = 4.d0*a2*a*dexp(e/(ta*ta))/srpi

  asm = 1d0/rsm
  tasm = 2.d0*asm
  asm2 = asm*asm
  tasm2 = 2.d0*asm2
  ccsm = 4d0*asm2*asm*dexp(e/(tasm*tasm))/srpi

  do  20  ir = 1, nkd
     r(1) = alat*(p(1)-dlv(1,ir))
     r(2) = alat*(p(2)-dlv(2,ir))
     r(3) = alat*(p(3)-dlv(3,ir))
     call sylm(r,yl,lmax,r2)
     r1 = dsqrt(r2)
     lzero = .false.

     ! --- Make the xi's from -1 to lmax ---
     if (r1 < 1d-6) then
        do  31  l = 1, lmax
           chi1(l) = 0d0
           chi2(l) = 0d0
31      enddo
        if (lpos) then
           !         do  l = 1, lmax
           !           zchi1(l) = 0d0
           !           zchi2(l) = 0d0
           !         enddo
           zuplus = zerfc(zikap/ta)
           chi1(0) = -zikap*zuplus &
                +2/dsqrt(4*datan(1d0))*a * cdexp(-(zikap/ta)**2)
           zuplus = zerfc(zikap/tasm)
           chi2(0) = -zikap*zuplus &
                +2/dsqrt(4*datan(1d0))*asm * cdexp(-(zikap/tasm)**2)
           chi1(-1) = -zerfc(-zikap/ta)/zikap
           chi2(-1) = -zerfc(-zikap/tasm)/zikap
        else
           chi1(0) = ta*dexp(e/(2d0*ta2))/srpi - akap*derfc(akap/ta)
           chi1(-1) = derfc(akap/ta)/akap
           chi2(0) = tasm*dexp(e/(2d0*tasm2))/srpi - akap*derfc(akap/tasm)
           chi2(-1) = derfc(akap/tasm)/akap
        endif
     else
        if (lpos) then
           expikr = exp(zikap*r1)
           zuplus = zerfc(zikap/ta+r1*a)*expikr
           chi1(0) = expikr/r1 - dble(zuplus)/r1
           chi1(-1) = expikr/zikap + dimag(zuplus)/kap
           !          zchi1(0) = expikr/r1 - dble(zuplus)/r1
           !          zchi1(-1) = expikr/zikap + dimag(zuplus)/kap
        else
           !         If large, these are -log(uplus),-log(umins); then chi->0
           !r        If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
           !r        -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r
           !r                    =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
           !r         u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2] -> 0
           !r        Also if akap*r >> 1,   chi < dexp(-akap*r1) -> 0
           emkr = dexp(-akap*r1)
           if (.5d0*akap/a+r1*a > srfmax .AND. &
                .5d0*akap/a-r1*a > srfmax .OR. akap*r1 > fmax) then
              lzero = .true.
              !            uplus = derfc(.5d0*akap/a+r1*a)/emkr
              !            umins = derfc(.5d0*akap/a-r1*a)*emkr
              !            print *, 'approx',(.5d0*akap/a+r1*a)**2 - akap*r1
              !            print *, '-log up',-dlog(uplus)
              !            print *, 'approx', (.5d0*akap/a-r1*a)**2 + akap*r1
              !            print *, '-log um', -dlog(umins)
              !            print *, r1*a
              !            stop
           else
              uplus = derfc(.5d0*akap/a+r1*a)/emkr
              umins = derfc(.5d0*akap/a-r1*a)*emkr
              chi1(0) = 0.5d0*(umins-uplus)/r1
              chi1(-1) = (umins+uplus)/(2.d0*akap)
           endif
        endif
        if (lzero) then
           do  30  l = -1, lmax
              chi1(l) = 0
30         enddo
           lzero = .false.
        else
           gl = cc*dexp(-a2*r2)/ta2
           do  32  l = 1, lmax
              chi1(l) = ((2*l-1)*chi1(l-1)-e*chi1(l-2)-gl)/r2
              gl = ta2*gl
32         enddo
        endif

        !        chi2 is complex; so is chi1, but the imaginary part
        !        is the same, so the difference is real
        if (lpos) then
           zuplus = zerfc(zikap/tasm+r1*asm)*expikr
           chi2(0) = expikr/r1 - dble(zuplus)/r1
           chi2(-1) = expikr/zikap + dimag(zuplus)/kap
           !          zchi2(0) = expikr/r1 - dble(zuplus)/r1
           !          zchi2(-1) = expikr/zikap + dimag(zuplus)/kap
        else
           if (.5d0*akap/asm+r1*asm > srfmax .AND. &
                .5d0*akap/asm-r1*asm > srfmax .OR. akap*r1 > fmax)then
              lzero = .true.
           else
              uplus = derfc(.5d0*akap/asm+r1*asm)/emkr
              umins = derfc(.5d0*akap/asm-r1*asm)*emkr
              chi2(0) = 0.5d0*(umins-uplus)/r1
              chi2(-1) = (umins+uplus)/(2d0*akap)
           endif
        endif
        if (lzero) then
           do  40  l = -1, lmax
              chi2(l) = 0
40         enddo
           lzero = .false.
        else
           gl = ccsm*dexp(-asm2*r2)/tasm2
           do  33  l = 1, lmax
              chi2(l) = ((2*l-1)*chi2(l-1)-e*chi2(l-2)-gl)/r2
              gl = tasm2*gl
33         enddo
        endif
     endif
     qdotr = tpi*(q(1)*dlv(1,ir)+q(2)*dlv(2,ir)+q(3)*dlv(3,ir))
     cfac = cdexp(dcmplx(0d0,qdotr))
     ilm = 0
     do  38  l = 0, lmax
        nm = 2*l+1
        do  39  m = 1, nm
           ilm = ilm+1
           dl(ilm) = dl(ilm) + yl(ilm)*(chi2(l)-chi1(l))*cfac
           dlp(ilm) = dlp(ilm) + yl(ilm)*0.5d0*(chi2(l-1)-chi1(l-1))*cfac
39      enddo
        cfac = cfac*dcmplx(0d0,1d0)
38   enddo
20 enddo

end subroutine hsmbld

