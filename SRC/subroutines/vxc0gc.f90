subroutine vxc0gc(nr,nsp,rofi,rwgt,rho,vxc,exc,rep,rmu,lxcfun)
  !- Gradient-corrected part of vxc and exc in a spherical potential
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   rofi  :radial mesh points
  !i   rwgt  :radial mesh weights
  !i   rho   :(true rho)*(4*pi*r**2)
  !i   lxcfun:specifies exchange-correlation functional
  !i         :1s digit sets local xc functional
  !i         :  1    Ceperly-Alder
  !i         :  2    Barth-Hedin (ASW fit)
  !i         :  3,4  LD part of PW91 and PBE
  !i         :10s digit sets gradient corrections
  !i         :  0    LSDA
  !i         :  1    Langreth-Mehl
  !i         :  2    PW91
  !i         :  3    PBE
  !i         :  4    PBE with Becke exchange
  !i         :100s digit sets flag to treat core in perturbation th.
  !o Outputs
  !o   vxc   :nonlocal XC potential
  !o   exc   :nonlocal XC energy
  !o   rep   :int rho * exc
  !o   rmu   :int rho * vxc
  !l Local variables
  !r Remarks
  !u Updates
  !u   05 Apr 09  Use rwgt for mesh weights
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lxcfun,nr,nsp
  double precision :: rofi(nr),rwgt(nr),vxc(nr,nsp),rho(nr,nsp), &
       rep(nsp),rmu(nsp),exc(nr)
  ! ... Local parameters
  integer :: nrmx
  parameter (nrmx=1501)
  double precision :: pi,rho2,rho3,ub4pi,rp(nrmx*2),rho0(2)
  integer :: i,ir,isp,ogrh,oagrh,oggrh,ogrgag

  pi = 4d0*datan(1d0)
  ub4pi = 1d0/(4d0*pi)
  call dpzero(vxc, nr*nsp)
  call dpzero(exc, nr)
  call dpzero(rep, nsp)
  call dpzero(rmu, nsp)
  if (iabs(lxcfun)/100 == 0) return
  if (nr > nrmx) call rx('vxc0gc: nr > nrmx')

  ! --- Make true rho ---
  do  10  isp = 1, nsp
     rho2 = rho(2,isp)/rofi(2)**2
     rho3 = rho(3,isp)/rofi(3)**2
     rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
     rp(1+nr*(isp-1)) = rho0(isp)
     do  20  ir = 2, nr
        rp(ir+nr*(isp-1)) = rho(ir,isp)*ub4pi/rofi(ir)**2
20   enddo
10 enddo

  !      print *, 'rl,l=0'
  !      call prmr(21,rofi,rp,1)

  ! --- Gradient correction ---
  !      call defrr(ogrh , nrmx*nsp)
  !      call defrr(oggrh, nrmx*nsp)
  !      call defrr(oagrh, nrmx*(3*nsp-2))
  !      call defrr(ogrgag,nrmx*(2*nsp-1))
  !      call vxcgr2(nr,nsp,nr,rofi,rp,w(ogrh),w(oggrh),
  !     .w(oagrh),w(ogrgag),exc,vxc)
  call vxcgr2(nr,nsp,nr,rofi,rp, exc,vxc)
  !      call rlse(ogrh)
  do  24  i  = 1, nsp
     rep(i) = 0d0
     rmu(i) = 0d0
     do  22  ir = 1, nr
        rep(i) = rep(i) + rwgt(ir)*rho(ir,i)*exc(ir)
        rmu(i) = rmu(i) + rwgt(ir)*rho(ir,i)*vxc(ir,i)
22   enddo
24 enddo

end subroutine vxc0gc

