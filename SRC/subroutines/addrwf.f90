module m_addrwf
contains
  subroutine addrwf(mode,z,l,v,ndg,n1,nr,rofi,rwgt,eadd,ev,fac,gadd, g,s)
    !- Add constant * radial wave function to another radial wave function
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  mode   :0 use both large and small components of radial w.f.
    !i         :1 use both large component of radial w.f. only
    !i   z     :nuclear charge
    !i         :(used only to compute overlap s, mode=0)
    !i   l     :l-quantum number
    !i   v     :spherical potential, without nuclear part
    !i         :(used only to compute overlap s, mode=0)
    !i   ndg   :leading dimension of g and gadd
    !i   n1    :if 0<n1<=nr, rwgt(n1) is scaled by 2
    !i         :(see vxtrap.f)
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights for numerical integration
    !i   ev    :eigenvalue of input wave function g
    !i         :(used only to compute overlap s, mode=0)
    !i   eadd  :eigenvalue of wave function gadd
    !i         :(used only to compute overlap s, mode=0)
    !i   fac   :Add fac*gadd into g
    !i   gadd  :See fac
    !o Inputs/Outputs
    !o   g     :g is overwritten by g + fac*g
    !o   s     :overlap between gadd and new g
    !r Remarks
    !r   Input g and gadd are assumed to be solutions of the Schrodinger
    !r   equation with eigenvalues ev and eadd.  (For the scalar
    !r   relativistic case, the inner product depends slightly
    !r   on z,v, and eigenvalue)
    !u Updates
    !u   04 Sep 04 Adapted to extended local orbitals
    !u   12 Jul 04 ndg,n1 arguments (altered argument list)
    !u   14 Feb 02 New routine
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: l,ndg,n1,nr,mode
    real(8):: fac,rofi(nr),rwgt(nr), &
         gadd(ndg,2),g(ndg,2),s
    real(8),optional:: z,v(nr),ev,eadd
    ! ... Local parameters
    integer :: ir
    double precision :: cc,vi,fllp1,gf11,gf22,gf12,r,tmc
    common /cc/ cc
    fllp1 = l*(l+1)
    s = 0
    if (n1 > 0 .AND. n1 < nr) rwgt(n1) = 2*rwgt(n1)
    if (mode == 0) then
       do  ir = 2, nr
          r = rofi(ir)
          if (fac /= 0) then
             g(ir,1) = g(ir,1) + fac*gadd(ir,1)
             g(ir,2) = g(ir,2) + fac*gadd(ir,2)
          endif
          !       Rest of the loop computes overlap between new g and gadd
          vi = v(ir) - 2d0*z/r
          tmc = cc - (vi-ev)/cc
          gf11 = 1d0 + fllp1/(tmc*r)**2
          tmc = cc - (vi-eadd)/cc
          gf22 = 1d0 + fllp1/(tmc*r)**2
          gf12 = (gf11 + gf22)/2
          s = s + rwgt(ir)*(gf12*g(ir,1)*gadd(ir,1) + g(ir,2)*gadd(ir,2))
       enddo
    else
       do  ir = 2, nr
          r = rofi(ir)
          if (fac /= 0) then
             g(ir,1) = g(ir,1) + fac*gadd(ir,1)
          endif
          s = s + rwgt(ir)*g(ir,1)*gadd(ir,1)
       enddo
    endif
    if (n1 > 0 .AND. n1 < nr) rwgt(n1) = rwgt(n1)/2
  end subroutine addrwf
end module m_addrwf

subroutine wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
     phzp,dphzp,pnz,rsml,ehl,g0,g1,gz)
  !- Add a linear combination of two w.f. to a 3rd to make a local orbital
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   l     :l quantum number
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   rwgt  :radial mesh weights
  !i   phi   :1st wave function at rmax, i.e. r*g0
  !i   dphi  :radial derivative of r*g0
  !i   phip  :2nd wave function at rmax, i.e. r*g1
  !i   dphip :radial derivative of r*g1
  !i   g0    :1st radial w.f. to which to orthogonalize gz
  !i   g1    :2nd radial w.f. to which to orthogonalize gz
  !i   ehl   :energy of smoothed Hankel tail for extended local orbital
  !i   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
  ! o Inputs/Outputs
  ! o  phz   :3rd wave function at rmax, i.e. r*gz
  ! o        :Input for standard local orbital; scaled for extended loc. orb
  ! o  dphz  :radial derivative of r*gz
  ! o        :Input for standard local orbital; scaled for extended loc. orb
  ! o  phzp  :energy derivative of phz
  ! o        :scaled for extended local orbital
  ! o  dphzp :energy derivative of dphz
  ! o        :scaled for extended local orbital
  ! o  gz    :on input, radial w.f.
  ! o        :on output, gz is overwritten by:
  ! o        :(gz - alpha g0 - beta g1) so that value and slope
  ! o        :of the result are zero at rmt, standard local orbital
  ! o        :scaled to match sm. hankel at rmt, extended loc. orb.
  !u Updates
  !u   04 Sep 04 Adapted to extended local orbitals
  !u   06 Mar 02 New routine
  ! ----------------------------------------------------------------------
  !     implicit none
  use m_addrwf
  ! ... Passed parameters
  integer :: l,nr
  double precision :: a,rofi(nr),rwgt(nr),rsml(*),ehl(*)
  double precision :: phi,dphip,dphi,phip,phz,dphz,phzp,dphzp,pnz(*)
  double precision :: g0(nr,2),g1(nr,2),gz(nr,2)
  ! ... Local parameters
  integer :: k,lpzi,nrmx
  parameter (nrmx=1501)
  double precision :: det,au,bu,as,bs,fac,x,xx
  double precision :: gzbig(nrmx*2)
  k = l+1
  lpzi = 0
  if (pnz(k) >  0) lpzi = 1
  if (pnz(k) >= 10) lpzi = 2
  if (pnz(k) >= 20) lpzi = 3
  if (lpzi == 0) then
     return
  elseif (lpzi == 1) then
     det = phi*dphip - dphi*phip
     au = dphip/det
     bu = -dphi/det
     as = -phip/det
     bs = phi/det
     fac = phz*au + dphz*as
     x = 0
     call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g0(1,2),g=gz(1,2),s=xx)
     call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g0,g=gz,s=xx)
     fac = phz*bu + dphz*bs
     call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g1(1,2),g=gz(1,2),s=xx)
     call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g1,g=gz,s=xx)
  elseif (lpzi == 2 .OR. lpzi == 3) then
     call dcopy(nr,gz,1,gzbig,1)
     call dcopy(nr,gz(1,2),1,gzbig(1+nr),1)
     call rwftai(5,rofi(nr),a,nr,nr,rofi,phz,dphz,xx,l, &
          ehl(k),rsml(k),gzbig)
     !       If rwftai scales gzbig, rescale phz,gz
     if (gzbig(nr) /= gz(nr,1)) then
        xx = gzbig(nr)/gz(nr,1)
        phz   = phz*xx
        dphz  = dphz*xx
        phzp  = phzp*xx
        dphzp = dphzp*xx
        call dscal(nr,xx,gz(1,1),1)
        call dscal(nr,xx,gz(1,2),1)
     endif
  endif
end subroutine wf2lo

subroutine ortrwf(mode,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,ez,g0,g1,gz,D)
  !- Orthogonalize a radial wave function gz to a pair of other functions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :0 do not change gz, but return scaling factor D that would
  !i         :  normalize gz after orthogonalization to (g0,g1)
  !i         :1 orthonormalize gz.
  !i         :  NB: this routine assumes g0 and g1 are orthogonal
  !i         :2 orthonormalize g0,g1; do not change gz or compute D
  !i         :3 orthogonalize g1 to g0; do not normalize
  !i         :4 orthogonalize g0 and g1 to gz
  !i         :10s digit
  !i         :0 use both large and small components of radial w.f.
  !i         :1 use large component of radial w.f. only.
  !i            In this case, z,v,e0,e1,ez are not used
  !i   z     :nuclear charge
  !i   l     :l quantum number
  !i   v     :spherical potential (atomsr.f)
  !i   ng    :leading dimension of g and gadd
  !i   n1    :if 0<n1<=nr, rwgt(n1) is scaled by 2
  !i         :(see vxtrap.f)
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   rwgt  :radial mesh weights
  !i   e0    :energy eigenvalue of g0
  !i   e1    :energy eigenvalue of g1
  !i   ez    :energy eigenvalue of gz
  !i   g0    :1st radial w.f. to which to orthogonalize gz
  !i   g1    :2nd radial w.f. to which to orthogonalize gz
  !i   gz    :radial w.f. to orthogonalize
  !o Outputs
  !o   D     :scaling factor that normalizes the orthogonalized gz
  !l Local variables
  !l         :
  !r Remarks
  !r
  !b Bug
  !b   for 1s digit mode=1, this routine assumes g0 and g1 are orthogonal
  !u Updates
  !u   12 Jul 04 Add option 3 to mode.  New argument list
  !u   06 Mar 02 New routine
  ! ----------------------------------------------------------------------
  !     implicit none
  use m_addrwf
  ! ... Passed parameters
  integer :: mode,l,ng,n1,nr
  double precision :: z,v(nr),rofi(nr),rwgt(nr),e0,e1,ez,D
  double precision :: g0(ng,2),g1(ng,2),gz(ng,2)
  ! ... Local parameters
  integer :: mode0,mode1
  double precision :: s00,s01,s11,s0z,s1z,szz,x,s01hat,s11hat,s1zhat
  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)
  ! --- mode 2 : orthonormalize g0 and g1 ---
  if (mode0 == 2) then
     !       <g0 g0>
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
     call dscal(nr,1/sqrt(s00),g0(1,1),1)
     call dscal(nr,1/sqrt(s00),g0(1,2),1)
     !       <g0 g1>
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
     call daxpy(nr,-s01,g0(1,1),1,g1(1,1),1)
     call daxpy(nr,-s01,g0(1,2),1,g1(1,2),1)
     !       <g1 g1>
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
     call dscal(nr,1/sqrt(s11),g1(1,1),1)
     call dscal(nr,1/sqrt(s11),g1(1,2),1)
     !       Check
     !        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
     !        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
     !        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
     return
     ! --- mode 3 : orthogonalize g1 to g0 ---
  elseif (mode0 == 3) then
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
     call daxpy(nr,-s01/s00,g0(1,1),1,g1(1,1),1)
     call daxpy(nr,-s01/s00,g0(1,2),1,g1(1,2),1)
     !       Check
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
     return
  endif
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,ez,ez,0d0,gz,gz,szz)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,0d0,g0,gz,s0z)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g1,gz,s1z)
  ! --- mode 4 : orthogonalize g0 and g1 to gz ---
  if (mode0 == 4) then
     call daxpy(nr,-s0z/szz,gz(1,1),1,g0(1,1),1)
     call daxpy(nr,-s0z/szz,gz(1,2),1,g0(1,2),1)
     call daxpy(nr,-s1z/szz,gz(1,1),1,g1(1,1),1)
     call daxpy(nr,-s1z/szz,gz(1,2),1,g1(1,2),1)
     !       Check
     !        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g0,gz,s0z)
     !        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g1,gz,s1z)
     return
  endif
  s01hat = s01/sqrt(s00)
  s11hat = s11 - s01hat**2
  s1zhat = s1z - s01*s0z/s00
  !     Scaling factor that normalizes the orthogonalized gz
  !     D = sqrt(szz - s0z**2/s00 - s1z**2/s11)
  D = sqrt(szz - s0z**2/s00 - s1zhat**2/s11hat)
  if (mode0 == 0) return
  !     Orthogonalize
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,-s0z/s00,g0,gz,x)
  call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,-s1z/s11,g1,gz,x)
  !     Normalize
  call dscal(nr,1/D,gz(1,1),1)
  call dscal(nr,1/D,gz(1,2),1)
end subroutine ortrwf

