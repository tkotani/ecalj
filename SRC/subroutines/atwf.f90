subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0, &
     nphimx,ncore,konfig,ecore,gcore,gval,nmcore)
  use m_ftox
  !- Make properties related to core for one sphere
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 return ncore, and konfig, and nphimx only;
  !i         :  see description below for contents of nphimx
  !i         :1s digit
  !i         :1 return valence wave functions
  !i         :2 return core wave functions
  !i         :3 combination of 1+2
  !i         :10s digit concerns orthogonalization
  !i         :0 do not orthogonalize
  !i         :1 return orthogonalized to valence orbitals
  !i         :2 return orthogonalized to valence orbitals
  !i         :  using large component only
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
  integer :: l,isp,konf,konfz,k,mode0,mode1,  nmcore
  double precision :: sumtc,sumec,e,ez,xx
  double precision :: rofi(nr),rwgt(nr),rhoc(nr,2),gp(2*nr*4)
  double precision :: phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp
  logical:: isanrg, l_dummy_isanrg
  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)
  ! --- Count number of core states ---
  lpz = .false.
  ncore = 0
  do  l = 0, lmxa
     k = l+1
     konfig(k) = pnu(k,1)
     konfz = mod(pnz(k,1),10d0)
     if (konfz == 0) konfz = konfig(k)
     l_dummy_isanrg=isanrg(konfz,konfig(k)-1,konfig(k)+1,'atwf:','pnuz',.true.)
     !       lpz = konfz .ne. konfig(k)
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
           call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(1,isp),4, &
                gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
           gval(:,l,2,isp)=gp(1:2*nr) !call dcopy(2*nr,gp,1,gval(1,l,2,isp),1)
           !         phi,phidot already orthogonal if mode1=1
           if (mode1 == 2) &
                call ortrwf(10*(mode1-1)+2,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt, &
                e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)
           !     ... Make local orbital
           if (konf /= konfig(k)) then
              l_dummy_isanrg=isanrg(nphimx,3,3,'atwf:','nphimx',.true.)
              call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(1,isp),2, &
                   gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)
              l_dummy_isanrg=isanrg(mode1,0,2,'atwf:','10s digit mode',.true.)
              if (mode1 == 0) then
                 call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
                      phzp,dphzp,pnz(1,isp),rsml,ehl, &
                      gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
              elseif (pnz(l+1,isp) < 10) then
                 call ortrwf(10*(mode1-1)+1,z,l,v0(1,isp),nr,nr,nr,rofi, &
                      rwgt,e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp), &
                      gval(1,l,3,isp),xx)
              endif
           endif
        enddo
     enddo
  endif
  ! --- Core eigenfunctions and eigenvalues ---
  if (mode0 >= 2) then
     call getcor(1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,0d0,sumec,sumtc, &
          rhoc,ncore,ecore,gcore,nmcore) !nmcore jun2012
  endif
end subroutine atwf

