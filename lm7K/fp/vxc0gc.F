      subroutine vxc0gc(nr,nsp,rofi,rwgt,rho,vxc,exc,rep,rmu,lxcfun)
C- Gradient-corrected part of vxc and exc in a spherical potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   rho   :(true rho)*(4*pi*r**2)
Ci   lxcfun:specifies exchange-correlation functional
Ci         :1s digit sets local xc functional
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci         :10s digit sets gradient corrections
Ci         :  0    LSDA
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Ci         :100s digit sets flag to treat core in perturbation th.
Co Outputs
Co   vxc   :nonlocal XC potential
Co   exc   :nonlocal XC energy
Co   rep   :int rho * exc
Co   rmu   :int rho * vxc
Cl Local variables
Cr Remarks
Cu Updates
Cu   05 Apr 09  Use rwgt for mesh weights
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lxcfun,nr,nsp
      double precision rofi(nr),rwgt(nr),vxc(nr,nsp),rho(nr,nsp),
     .rep(nsp),rmu(nsp),exc(nr)
C ... Local parameters
      integer nrmx
      parameter (nrmx=1501)
      double precision pi,rho2,rho3,ub4pi,rp(nrmx*2),rho0(2)
      integer i,ir,isp,ogrh,oagrh,oggrh,ogrgag
      integer w(1)
      common /w/ w

      pi = 4d0*datan(1d0)
      ub4pi = 1d0/(4d0*pi)
      call dpzero(vxc, nr*nsp)
      call dpzero(exc, nr)
      call dpzero(rep, nsp)
      call dpzero(rmu, nsp)
      if (iabs(lxcfun)/100 .eq. 0) return
      if (nr .gt. nrmx) call rx('vxc0gc: nr > nrmx')

C --- Make true rho ---
      do  10  isp = 1, nsp
        rho2 = rho(2,isp)/rofi(2)**2
        rho3 = rho(3,isp)/rofi(3)**2
        rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
        rp(1+nr*(isp-1)) = rho0(isp)
        do  20  ir = 2, nr
   20   rp(ir+nr*(isp-1)) = rho(ir,isp)*ub4pi/rofi(ir)**2
   10 continue

C      print *, 'rl,l=0'
C      call prmr(21,rofi,rp,1)

C --- Gradient correction ---
c      call defrr(ogrh , nrmx*nsp)
c      call defrr(oggrh, nrmx*nsp)
c      call defrr(oagrh, nrmx*(3*nsp-2))
c      call defrr(ogrgag,nrmx*(2*nsp-1))
c      call vxcgr2(nr,nsp,nr,rofi,rp,w(ogrh),w(oggrh),
c     .w(oagrh),w(ogrgag),exc,vxc)
      call vxcgr2(nr,nsp,nr,rofi,rp, exc,vxc)
c      call rlse(ogrh)
      do  24  i  = 1, nsp
        rep(i) = 0d0
        rmu(i) = 0d0
        do  22  ir = 1, nr
          rep(i) = rep(i) + rwgt(ir)*rho(ir,i)*exc(ir)
          rmu(i) = rmu(i) + rwgt(ir)*rho(ir,i)*vxc(ir,i)
   22   continue
   24 continue

      end

