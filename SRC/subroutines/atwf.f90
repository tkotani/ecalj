module m_atwf !- Make properties related to core for one sphere
  public atwf
  private
contains
  subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0, & 
       nphimx,ncore,konfig,ecore,gcore,gval,nmcore)
    use m_ftox
    use m_rhocor,only: getcor
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 return ncore, and konfig, and nphimx only;
    !i         :  see description below for contents of nphimx
    !i         :1s digit
    !i         :1 return valence wave functions
    !i         :2 return core wave functions
    !i         :3 combination of 1+2
    !i       using large component only
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   lmxa  :augmentation l-cutoff
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pnz   :pnu por local orbitals
    !i   rmt   :MT boundary
    !i   z     :nuclear charge      (not used if mode=0)
    !i   v0    :spherical potential (not used if mode=0)
    !i   ehl   :energy of smoothed Hankel tail for extended local orbital
    !i   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
    ! o Inputs/Outputs
    ! o  nphimx:dimensions gval.  Must be at least as large as the
    ! o        :number of valence wave functions
    ! o        :For mode=0, nphimx is output and is assigned to
    !i         :maximum number radial wave functions for any l channel.
    !o Outputs
    !o   ncore :number of core levels
    !o   konfig:1s digit contains core configuration
    !o         :10s digit:
    !o         : 0 -> no local orbitals
    !o         : 1 -> local orbital with p.q.n. < pnu
    !o         : 2 -> local orbital with p.q.n. > pnu
    !o   ... The following are not used if mode=0
    !o   ecore :core eigenvalues
    !o   gcore :core wave functions
    !o   gval  :valence wave functions
    !o          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
    !o            i=0 : phi
    !o            i=1 : phidot
    !o            i=2 : local orbital
    !r Remarks
    !u Updates
    !u    4 Sep 04 Adapted to extended local orbitals
    !u   22 Dec 01 Adjustments to accomodate changes in phidx
    !u   22 Apr 01 Created by MvS
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,nr,nsp,lmxa,ncore,konfig(1+lmxa),n0,nrmx,nphimx
    parameter (n0=10)
    double precision :: rmt,z,a,v0(nr,nsp),pnu(n0,nsp),pnz(n0,nsp), &
         gval(nr*2,0:lmxa,nphimx,nsp),ecore(*),gcore(nr,2,*), rsml(n0),ehl(n0)
    logical :: lpz
    integer :: l,isp,konf,konfz,k,mode0,  nmcore !,mode1
    double precision :: sumtc,sumec,e,ez,xx
    double precision :: rofi(nr),rwgt(nr),rhoc(nr,2),gp(2*nr*4)
    double precision :: phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp
    logical:: isanrg, l_dummy_isanrg
    mode0 = mod(mode,10)
    ! --- Count number of core states ---
    lpz = .false.
    ncore = 0
    do  l = 0, lmxa
       k = l+1
       konfig(k) = pnu(k,1)
       konfz = mod(pnz(k,1),10d0)
       if (konfz == 0) konfz = konfig(k)
       l_dummy_isanrg=isanrg(konfz,konfig(k)-1,konfig(k)+1,'atwf:','pnuz',.true.)
       do  konf = l+1, min(konfz,konfig(k))-1
          ncore = ncore+nsp
       enddo
       if (konfz < konfig(k)) then
          konfig(k) = konfz + 10
          lpz = .true.
       elseif (konfz > konfig(k)) then
          konfig(k) = konfig(k) + 20
          lpz = .true.
       endif
    enddo
    if (mode0 == 0) then
       nphimx = 2
       if (lpz) nphimx = 3
       return
    endif
    call radmsh(rmt,a,nr,rofi)
    call radwgt(rmt,a,nr,rwgt)
    ! --- Valence wave functions ---
    if (mod(mode0,2) == 1) then
       do  l = 0, lmxa
          k = l+1
          do  isp = 1, nsp
             konf = pnu(k,1)
             call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(1,isp),4, gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
             gval(:,l,2,isp)=gp(1:2*nr) 
             if (konf /= konfig(k)) then
                call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(1,isp),2, gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)
                call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
                     phzp,dphzp,pnz(1,isp),rsml,ehl, gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
             endif
          enddo
       enddo
    endif
    if (mode0 >= 2) call getcor(1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,[0d0,0d0],sumec,sumtc, &
         rhoc,ncore,ecore,gcore,nmcore) !nmcore jun2012 !Core eigenfunctions and eigenvalues --- 
  end subroutine atwf
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
    use m_addrwf,only: addrwf
    use m_vxtrap,only: rwftai
    implicit none
    integer :: l,nr
    double precision :: a,rofi(nr),rwgt(nr),rsml(*),ehl(*)
    double precision :: phi,dphip,dphi,phip,phz,dphz,phzp,dphzp,pnz(*)
    double precision :: g0(nr,2),g1(nr,2),gz(nr,2)
    integer :: k,lpzi
    double precision :: det,au,bu,as,bs,fac,x,xx
    double precision :: gzbig(nr*2)
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
       call rwftai(rofi(nr),a,nr,nr,rofi,phz,dphz,xx,l, ehl(k),rsml(k),gzbig)
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
end module m_atwf
