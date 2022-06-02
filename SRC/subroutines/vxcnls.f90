subroutine vxcnls(a,ri,lcut,nr,np,nlm,nsp, yl,gyl,rwgt,wp,rl,lxcg, vl,rep,rmu) 
  use m_xcpbe,  only: xcpbe
  !!= Gradient correction to nspher. density on a radial and angular mesh =
  !!*NOTE: seeing tol=1d-10 is changed on Dec1st 2010. But this is empirically determined for Co atom.
  !! If GGA is unstable than LDA, Check this routine. I am not so definite for tol and so.
  !!----------------------------------------------------------------------
  !!*Inputs
  !!   ri    :mesh of points
  !!   lcut  :1 if cutoff exc for small rho to avoid blowup in exc
  !!   nr    :number of radial mesh points
  !!   np    :number of points for angular integration
  !!   nlm   :maximum (l+1)**2
  !!   nsp   :2 for spin-polarized case, otherwise 1
  !!   yl    :Ylm's tabulated on angular mesh
  !!   gyl   :Gradient of Ylm's tabulated on angular mesh
  !!   rwgt  :radial mesh weights
  !!   wp    :angular mesh weights (not needed unless debugging)
  !!   rl    :density on radial mesh
  !!   agrl  :|rl|
  !!   lxcg : unused in new mode (PBE-GGA lxcg=103 only).
  !!*Outputs
  !!   vl    :vxc in GGA
  !!   rep   :\int rho * exc
  !!   rmu   :\int rho * vxc
  !!
  !!*Local variables
  !!   rp    :spin pol density on the combined radial and angular mesh
  !!   ylwp  :Ylm*wp, where wp are weights for angular mesh
  !!   grp   :grad rp
  !!   ggrp  :Laplacian rp
  !!** these are not used in new mode
  !!   agrp   :|grad rp| Not used in new mode
  !!
  !!*Updates
  !!  Dec1st 10 tol=1d-10.
  !!  Previous setting has a problem for Kino's Co MMOM= 0 0 4 0 case.
  !!      Sep 10 takao
  !!   05 Apr 09 reduced the calling arguments for better portability
  !!   29 Apr 05 (ATP) adaped to new vxcnsp
  !!----------------------------------------------------------------------
  implicit none
  integer :: lcut,nr,np,nlm,nsp,lxcg
  double precision :: ri(nr),gyl(np,nlm,3),ylwp(np,nlm), &
       wp(np),rep(2),rmu(2),rwgt(nr),vxcnl(nr,nlm,nsp),excnl(nr,nlm), &
       rl(nr,nlm,nsp),agrl(nr,nlm,nsp),agrp(nr,np,nsp), &
       vl(nr,nlm,nsp),rp(nr,np,nsp),grp(nr,np,3,nsp),ggrp(nr,np,nsp)

  double precision :: pi,vhold,sumnl(0:20),repnl(4),rmunl(4),weight,weightirr
  double precision :: wk(nr,nsp),wk2(nr,nsp)
  integer :: ilm,ip,ipr,ir,i,l,ll,lmax,nn,oagrp,ogagrp,lopgfl, &
       ovxcp,oexcp,owk,owkl,lmx,jx
  integer :: ir0,ir1,nri,incn
3584 logical lz
  parameter (incn=50)
  real(8):: tol=1d-10

  integer:: irr,j,idiv,isp,icenter,k
  real(8):: a,bb, & ! & takao used to get derivetive of ri. See vxcnls.F
  rdir(3,np),cy1,dydummy,drdp(nr),rmax,ddr,fac,dgr2dn(1:np),polinta,grptot(3)
  integer,parameter:: ndiv=3
  ! ino delete target      real(8),target :: yl(np,nlm)
  !      real(8) :: yl(np,nlm)
  real(8) :: yl(np,4)
  ! ino Dec.12.2011 delete drp_drl      real(8),pointer:: drp_drl(:,:) =>NULL()
  real(8),allocatable:: dgrp_drl(:,:,:,:,:), &
       exci(:,:),vxci(:,:,:),dvxcdgr(:,:,:),grho2_updn(:,:,:),excl(:,:),vlxc(:,:,:)
  real(8),allocatable:: grho2_updn_forcall(:,:,:)

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCC
  logical:: newmode=.true., &
       debug=.false.,plotstop=.false., grpzerotest=.false.
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCC

  real(8),allocatable:: dvxcdgr_dr(:,:,:),nabla_dvxcdgr(:,:,:,:),vxcnlnp(:,:,:)
  real(8)::dvxcdgr_ilm,ddd(np),grpt(3),ggrpt
  logical::  lrat
  integer::lerr

  !      integer,allocatabel:: irr(:,:) ! irr(0:2,ir) specify radial index related to dgrp_drl.
  !     call pshpr(80)
  call tcn('vxcnls')
  if(debug) print *,'vxcnls: goto vxcnls 11111'
  if (lxcg == 0) return
  call getpr(ipr)
  pi = 4d0*datan(1d0)
  nn = 4
  call dpzero(repnl,4)
  call dpzero(rmunl,4)
  call dpzero(excnl,nr*nlm)
  call dpzero(vxcnl,nr*nlm*nsp)
  lz = ri(1) .eq. 0

  ! --- Make ylwp = yl*wp for fast multiplication ---
  do  ilm = 1, nlm
     do  ip = 1, np
        ylwp(ip,ilm) = yl(ip,ilm)*wp(ip)
     enddo
  enddo

  ! --- Generate density point-wise through sphere ---
  do  i = 1, nsp
     call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0, &
          rp(1,1,i),nr)
  enddo


  ! akao developing ccccccccccccccccccccccccccccccccccccccccccccc
  if(newmode) then
     ! --- Gradient of density point-wise through sphere ---
     !       if (lz) lopgfl = 10
     lopgfl=10
     do i = 1, nsp
        call gradfl(ll(nlm),nlm,nr,np,1,nr,1,lopgfl,nn,ri,yl,gyl, &
             rl(1,1,i),grp(1,1,1,i),ggrp(1,1,i))
     enddo

     !! -- If negative density, set to tol ---
     !! In the folloing of call xcpbe, corresponding vxc and so become zero. So they are dummy.
     do  i = 1, nsp
        do  ip = 1, np
           do  ir = 1, nr
              if (rp(ir,ip,i) <= tol) rp(ir,ip,i)   = tol
              if (rp(ir,ip,i) <= tol) grp(ir,ip,:,i)= tol*10
              if (rp(ir,ip,i) <= tol) ggrp(ir,ip,i) = tol*10
              !  takao modified this on Dec1st 2010 for stability. Now tol=1d-10
              ! Previous version is with tol=1d-20 grp=ggrp=tol=10*tol caused problem.
              ! I think we need to examine xcpbe.F90 in future. (because we suppose
              !   grp=ggrp=0 means LDA. and shuold be safer. But it is opposite).

           enddo
        enddo
     enddo


     !$$$C --- drpdrl, dgrpdrl !takao
     !$$$      if(lxcg/=3) call rx('vxcnls: GGA=3(PBE) is only implemented yet.')
     !$$$      drp_drl=>yl  !diagonal for radial mesh
     !$$$      bb = ri(nr)/(exp(a*(nr-1))-1d0) ! drdp = drofi/di, where rofi(i) = b [e^(a(i-1)) -1]
     !$$$      do ir=1,nr
     !$$$        drdp(ir)=a*bb*exp(a*(ir-1))
     !$$$      enddo
     ! f(p) = p(p-1)/2 * f_-1 + (1-p**2) *f_0 + p(p+1)/2 * f_1
     ! dfdp = (p-.5)* f_-1  - 2p * f_0 + (p+.5) * f_1
     ! d_dr = dfdp* dp/dr
     !      do i=-1,1
     !        d_dp(i,-1)= (i-.5d0) !dpdr(ir)
     !        d_dp(i,0) = (-2d0*i)
     !        d_dp(i,1) = i+.5d0
     !      enddo
     ! -- radial part of dgrp_drl


     !$$$      cy1 = dsqrt(3d0/(16*datan(1d0)))
     !$$$      do ip = 1, np
     !$$$         rdir(1,ip)=yl(ip,4)/cy1 !x
     !$$$         rdir(2,ip)=yl(ip,2)/cy1 !y
     !$$$         rdir(3,ip)=yl(ip,3)/cy1 !z
     !$$$      enddo
     !$$$      allocate( dgrp_drl(nr,np,3,-1:1,nlm ))
     !$$$      dgrp_drl=0d0
     !$$$      do irr=2,nr !1,nr
     !$$$      do idiv=-1,1,2
     !$$$        ddr= idiv*0.5d0/drdp(irr) ! -0.5d0 for idiv=-1, 0.5 for idiv=1
     !$$$        ir=irr+idiv
     !$$$        if(idiv>nr) cycle
     !$$$        do ip=1,np
     !$$$        do j=1,3
     !$$$          dgrp_drl(irr,ip,j,idiv,1:nlm)=
     !$$$     &    dgrp_drl(irr,ip,j,idiv,1:nlm) + rdir(j,ip)*ddr*drp_drl(ip,1:nlm)  !contribution from radial derivative. 2 is center
     !$$$          !Note : ir = irr+ idiv, where idiv=-1,0,1.
     !$$$        enddo
     !$$$        enddo
     !$$$      enddo
     !$$$      enddo
     !$$$C -- spherical part of dgdp. So smaller contribution than the radial parts.
     !$$$      idiv=0
     !$$$      do irr=2,nr
     !$$$        do ip=1,np
     !$$$        do j=1,3
     !$$$        dgrp_drl(irr,ip,j,idiv,1:nlm)=
     !$$$     &  dgrp_drl(irr,ip,j,idiv,1:nlm) + gyl(ip,1:nlm,j)/ri(irr) ! see gradfl. 1st digit of lopgfl=lx=0
     !$$$        enddo
     !$$$        enddo
     !$$$      enddo
     !$$$      if(debug) print *,'dgrp_drl sumcheck=',sum(abs(dgrp_drl))

     !      ilm=2
     !      ip=2
     !      do ir=2,nr
     !        write(6,"(i5,3(3d12.4,3x))")ir, dgrp_drl(ir,ip,1:3,-1:1,ilm)
     !      enddo

     ! -- grp is calculated (simple. but consistent with the definition of dgrp_drl, rather than
     !    gradfl. It might be better to use this.
     if(debug) print *,'vxcnls: aaa111'
     allocate(grho2_updn(nr,np,2*nsp-1) )
     do isp=1,nsp
        do ip=1,np

           !$$$        do j=1,3
           !$$$          grp(1,ip,j,isp) = 0d0 !dummy
           !$$$          do ir=4,nr-1
           !$$$            grp(ir,ip,j,isp) = sum( dgrp_drl(ir,ip,j,-1:1,1:nlm)*rl(ir-1:ir+1,1:nlm,isp) )
           !$$$          enddo
           !$$$          do ir=1,3  ! interpolation
           !$$$            call polinx(ri(4), grp(4,ip,j,isp),   4,ri(ir), 0d0,grp(ir,ip,j,isp) ,dydummy)
           !$$$          enddo
           !$$$          call polinx(ri(nr-4),grp(nr-4,ip,j,isp),4,ri(nr), 0d0,grp(nr,ip,j,isp) ,dydummy)
           !$$$        enddo
           if(grpzerotest) then
              grp=0d0
              ggrp=0d0
           endif
           do ir=1,nr
              grho2_updn(ir,ip,isp) = sum(grp(ir,ip,:,isp)**2)
              if(nsp==2) grho2_updn(ir,ip,3) = &
                   sum(grp(ir,ip,1,1:2))**2 + sum(grp(ir,ip,2,1:2))**2 +sum(grp(ir,ip,3,1:2))**2
           enddo
        enddo
     enddo

     ! ccccccccccccccccccccc
     !        isp=1
     !        do ir=1,nr
     !c          write(6,"(i5,8d12.4)") ir, grho2_updn(ir,1:8,isp)
     !        write(6,"('qqq ',i5,8d12.4)") ir, sum(abs(grp(ir,1,1:3,isp))),
     !     &  abs(ggrp(ir,1,isp))
     !        enddo
     ! cccccccccccccccccccc


     ! --- Non local part of exchange correlation, and its derivatives with respect to rp,grp.
     !      call xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci, & !Mandatory Arguments
     !     &           d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments
     !    option =2 ->PBE
     ! NOTE: rp(ir,np,isp=1) contains total electron density when nsp=1.
     if(debug) print *,'vxcnls: aaa222'
     allocate(exci(nr,np),excl(nr,nlm),vxci(nr,np,nsp),dvxcdgr(nr,np,3))
     fac=1d0/2d0  !This fac is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
     if(nsp==2) fac=1d0

     !! == call xcpbe, which is the exchange-correlation kernel for PBE-GGA related functionals ==
     ! Be careful  1:dexdn(up) 2:dexdn(dn) 3:decdn(tot) for
     ! input
     !  rho_updn: density for up and down. For nsp=1, only up is used--->total dnsity is twice of the up density.
     !  grho2_updn:square of gradient of up and down densities.
     !             When nsp=2, grho2_upsn(ir,ip,3) should contans the square for the total density.

     ! allocate for calling xcpbe
     allocate( grho2_updn_forcall(nr,np,2*nsp-1)  )
     do isp=1,2*nsp-1
        do ip=1,np
           do ir=1,nr
              grho2_updn_forcall(ir,ip,isp)=fac**2*grho2_updn(ir,ip,isp)
           enddo
        enddo
     enddo
     call xcpbe(exci=exci,npts=nr*np,nspden=nsp, &
          option=2,&! & Choice of the functional =2:PBE-GGA
     order=1, &!  order=1 means we only calculate first derivative of rho*exc(rho,\nable rho).
     rho_updn=fac*rp,vxci=vxci,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0, & ! & Mandatory Arguments
     !     &           dvxcdgr=dvxcdgr, grho2_updn=fac**2*grho2_updn)   !Optional Arguments
     dvxcdgr=dvxcdgr, grho2_updn=grho2_updn_forcall)   !Optional Arguments
     deallocate(grho2_updn_forcall)
     !! Output:
     exci = 2d0*exci !in Ry.
     vxci = 2d0*vxci !in Ry.
     dvxcdgr= 2d0*dvxcdgr !in Ry.

     !! ==== checkwrite ====
     if(debug) print *,'vxcnls: aaa333'
     if(plotstop) then
        isp=1
        do ir=1,nr
           write(2106,"(i5,18e12.4)") ir,ri(ir),vxci(ir,1:8,isp) !np-4:np,isp)
           write(3106,"(i5,18e12.4)") ir,ri(ir),dvxcdgr(ir,1:8,3) !np-4:np,isp)
           write(4106,"(i5,18e12.4)") ir,ri(ir),dvxcdgr(ir,1:8,1) !np-4:np,isp)
           write(5106,"(i5,18e12.4)") ir,ri(ir),rp(ir,1:8,1) !np-4:np,isp)
           write(6106,"(i5,18e12.4)") ir,ri(ir),sqrt(grho2_updn(ir,1:8,isp))/sqrt(rp(ir,1:8,1))
        enddo
     endif

     !$$$cccccccccccccccccccccccc
     !$$$c      print *,' exci sumcheck=',sum(abs(exci)),sum(abs(vxci)),sum(abs(dvxcdgr))
     !$$$C --- Calculate
     !$$$c      exl= exci
     !$$$c      vlxc = vxci*drp_drl + dvxcdgr*dgrp_drl
     !$$$Co   vlxc    :vxcnl is added to vlxc
     !$$$Co   rep   :int rho * excnl added to rep
     !$$$Co   rmu   :int rho * vxcnl added to rmu

     !      print *,'sumcheck sum(wp)=',sum(wp) !---> this gives =4*pi

     !! == vlxc is given ==
     ! definition of Exc
     ! Exc= \sum_{ir,ip} exci(ir,ip)* ri(ir)**2*rwgt(ip)*wp(ip)
     ! Normalization: Sphere Volume (in Bohr**3) = \sum_{ir,ip} ri(ir)**2*rwgt(ip)*wp(ip)
     ! Potential:   v(ir,theta,phi,isp)= \sum_ilm vlxc(ir,ilm,isp)*Y_ilm(theta,phi)
     !              ip=(theta phi)_ip
     !              1= \sum_ip w(ip) Y_ilm(ip) Y_ilm(ip)
     allocate(vlxc(nr,nlm,nsp))
     do ilm=1,nlm
        ! print *,'norm  check1=',ilm ,sum(yl(1:np,ilm)*ylwp(1:np,ilm))
        ! sum(yl*ylwp)=1 for each ilm with double precision accuracy
        do ir=1,nr
           excl(ir,ilm)   = sum( exci(ir,1:np)*ylwp(1:np,ilm) ) !drp_drl(1:np,ilm) )
           do isp=1,nsp
              vlxc(ir,ilm,isp) = sum( vxci(ir,1:np,isp)*ylwp(1:np,ilm) ) !drp_drl(1:np,ilm) )
           enddo
        enddo
     enddo

     !$$$cccccccccccccccccccccccc
     !$$$c      if(plotstop) then
     !$$$c      isp=1
     !$$$c      do ir=1,nr
     !$$$c        write(6,"(i5,100d12.4)") ir, dvxcdgr(ir,1:3,3),dvxcdgr(ir,np/2:np/2+3,3)
     !$$$cxxxxxxxx        write(6,"(i5,100d12.4)") ir, bsum(grp(irr,ip,j,1:nsp))vxcdgr(ir,1:3,3),dvxcdgr(ir,np/2:np/2+3,3)
     !$$$c      enddo
     !$$$c      stop 'test xxx xxxxxxxxxxxxxxx '
     !$$$c      endif
     !$$$ccccccccccccccccccccccccccccccccccccccc
     !$$$        vxcnl=0d0
     !$$$      do isp=1,nsp
     !$$$      do ilm=1,nlm
     !$$$        do irr=2,nr
     !$$$        do idiv=-1,1
     !$$$          ir = irr + idiv
     !$$$          if(ir==1) cycle
     !$$$          if(ir>nr) cycle
     !$$$          do ip=1,np
     !$$$            do j=1,3
     !$$$            grptot(j) =  sum(grp(irr,ip,j,1:nsp))
     !$$$            enddo
     !$$$            dgr2dn(ip) = sum(dgrp_drl(irr,ip,1:3,idiv,ilm) * grptot(1:3))
     !$$$          enddo
     !$$$          weight = ri(ir)**2*drdp(ir)
     !$$$          weightirr = ri(irr)**2*drdp(irr)
     !$$$          ddd= dvxcdgr(irr,1:np,isp)+dvxcdgr(irr,1:np,3)
     !$$$          vxcnl(ir,ilm,isp) = vxcnl(ir,ilm,isp)
     !$$$     &    + sum(ddd*dgr2dn(1:np)*wp(1:np))* weightirr/weight
     !$$$        enddo
     !$$$        enddo
     !$$$      enddo
     !$$$      enddo
     !$$$!! === interpolation to 1,2,3 and nr-1,nr ===
     !$$$      do isp=1,nsp
     !$$$      do ilm=1,nlm
     !$$$        do ir=1,3
     !$$$          call polinx(ri(4), vxcnl(4,ilm,isp), 4,ri(ir), 0d0,vxcnl(ir,ilm,isp) ,dydummy)
     !$$$        enddo
     !$$$        call polinx(ri(nr-5),vxcnl(nr-5,ilm,isp),4,ri(nr-1), 0d0,vxcnl(nr-1,ilm,isp) ,dydummy)
     !$$$        call polinx(ri(nr-5),vxcnl(nr-5,ilm,isp),4,ri(nr), 0d0,vxcnl(nr,ilm,isp) ,dydummy)
     !$$$      enddo
     !$$$      enddo
     !$$$      if(debug) print *, 'vxcnls: goto vxcnls 99999xxx'


     !! == a new way to calculate vxcnlc ==
     !!  vxcnl= - grad [ frac{\partial rho exc}{\partial |nabla n|}  * frac{\nabla n}{nabla n}  ]
     !!  This expresssion is obtained by a partial derivative of rho exc with respect to (nabla n).
     !!  We neglect surface terms since true and counter components are cancelled out.
     !!   dvxcdgr(ir,ip,3): index=3 means frac{\partial rho exc}{\partial |nabla n|}/|nabla n|. See xcpbe.F90

     !! == dvxcdgr(ir,1:np,3) ---> nabla_dvxcdgr(:,ir,1:np) ==
     allocate(dvxcdgr_dr(nr,np,3),vxcnlnp(nr,np,nsp))
     do isp=1,3 ! 1:dexdn(up) 2:dexdn(dn) 3:decdn(tot)
        if(nsp==1 .AND. isp==2) cycle !dvxcddgr(:,:,isp=2) is dummy if nsp=1
        do ip=1,np
           !        call poldvm(x=ri,y=dvxcdgr(1:nr,ip,isp),np=nr,n=nn,
           !     &  lrat=.false.,tol=1d-12,lerr=lerr,
           !     &  yp=dvxcdgr_dr(1,ip,isp) ) !radial derivarive of dvxcdgr
           call poldvm(ri,dvxcdgr(1,ip,isp),nr,nn, &
                .false.,1d-12,lerr, &
                dvxcdgr_dr(1,ip,isp) ) !radial derivarive of dvxcdgr

           !! tol is twiced to avoid strange interpolation in poldvm
           do ir=1,nr
              if(isp==3) then
                 if(sum(rp(ir,ip,1:nsp)) <= 4*tol) dvxcdgr_dr(ir,ip,isp) = 0d0
              else
                 if(rp(ir,ip,isp) <= 2*tol) dvxcdgr_dr(ir,ip,isp) = 0d0
              endif
           enddo
           if (lerr /= 0) call rxi('vxcnls: poldvm 111',1000*isp+ip)
        enddo
     enddo

     !        print *,'44444444444 vxcnls'
     !! ==== Split grad r- into x,y,z- components ====
     cy1 = dsqrt(3d0/(16*datan(1d0)))
     do ip = 1, np
        rdir(1,ip)=yl(ip,4)/cy1 !x
        rdir(2,ip)=yl(ip,2)/cy1 !y
        rdir(3,ip)=yl(ip,3)/cy1 !z
        !        print *,' rdir**2 sum=',sum(rdir(:,ip)**2)
     enddo
     !        print *,'5555555555555 vxcnls'
     !! convert dvxcdgr ->  dvxcdgr_ilm  np to ilm representation.
     allocate( nabla_dvxcdgr(1:3,nr,np,3) )
     nabla_dvxcdgr = 0d0
     do isp=1,3 ! 1:dexdn(up) 2:dexdn(dn) 3:decdn(tot)
        if(nsp==1 .AND. isp==2) cycle !dvxcddgr(:,:,isp=2) is dummy
        do ir=2,nr
           do ip=1,np
              nabla_dvxcdgr(1:3,ir,ip,isp) = dvxcdgr_dr(ir,ip,isp)*rdir(1:3,ip)
           enddo
           do ilm =1,nlm
              ddd= dvxcdgr(ir,1:np,isp) + dvxcdgr(ir,1:np,3)
              dvxcdgr_ilm = sum( ddd*ylwp(1:np,ilm) ) !(ir,ilm) component
              do ip=1,np
                 nabla_dvxcdgr(1:3,ir,ip,isp) = &
                      nabla_dvxcdgr(1:3,ir,ip,isp) + dvxcdgr_ilm/ri(ir)* gyl(ip,ilm,1:3)
              enddo
           enddo
        enddo
     enddo

     !! == Calculate vxcnl in np represetation, and then it is converted to nlm representation. ==
     ! ccccccccccccccccccc
     !      ggrp=0d0
     ! ccccccccccccccccccc

     fac=1d0/2d0
     if(nsp==2) fac=1d0
     do isp=1,nsp
        do ir=2,nr
           do ip=1,np
              do j=1,3
                 grpt(j)  = sum( grp(ir,ip,j,1:nsp))  !gradient for total density
              enddo
              ggrpt    = sum( ggrp(ir,ip,1:nsp) )  !laplacian for total density
              vxcnlnp(ir,ip,isp) = &
                   - fac*ggrp(ir,ip,isp) * dvxcdgr(ir,ip,isp) &
                   - fac*sum( grp(ir,ip,:,isp) * nabla_dvxcdgr(:,ir,ip,isp) ) &
                   - ggrpt * dvxcdgr(ir,ip,3) &
                   - sum( grpt(:) * nabla_dvxcdgr(:,ir,ip,3) )
           enddo
           do ilm=1,nlm
              vxcnl(ir,ilm,isp) = sum( vxcnlnp(ir,1:np,isp)*ylwp(1:np,ilm) )
           enddo
        enddo
     enddo
     ir=1
     do isp=1,nsp
        do ilm=1,nlm
           call polinx(ri(2),vxcnl(2,ilm,isp), 4,ri(ir), 0d0,vxcnl(ir,ilm,isp) ,dydummy)
        enddo
     enddo
     deallocate(nabla_dvxcdgr,dvxcdgr_dr,vxcnlnp)
     !!
     vl=vlxc + vxcnl

     !C$$$ccccccccccccccccccccccc
     if(plotstop) then
        isp=1
        do ir=1,nr
           !        write(6,"(i5,4d12.4,3x,4d12.4,3x,4d12.4)") ir,
           !     &  vlxc(ir,1:4,isp), vxcnl(ir,1:4,isp),vlxc(ir,1:4,isp)+vxcnl(ir,1:4,isp)
           write(106,"(i5,e12.4,2x,4(3e12.4,2x))") ir,ri(ir), &
                vl(ir,1,isp),  vlxc(ir,1 ,isp) ,vxcnl(ir,1 ,isp)
           !     &  ,vl(ir,11,isp), vlxc(ir,11 ,isp),vxcnl(ir,11,isp),
           !     &  ,vl(ir,21,isp), vlxc(ir,21 ,isp),vxcnl(ir,21,isp)
           write(1106,"(i5,e12.4,2x,4(3e12.4,2x))") ir,ri(ir), &
                sum(grp(ir,1,1:3,isp)**2) ,(rl(ir,1,isp))**2 , sum(grp(ir,1,1:3,isp)**2)/(rl(ir,1,isp))**2
           !     &  ,sum(grp(ir,11,1:3,isp)**2),(rl(ir,11,isp))**2,sum(grp(ir,11,1:3,isp)**2)/(rl(ir,11,isp))**2
           !     &  ,sum(grp(ir,21,1:3,isp)**2),(rl(ir,21,isp))**2,sum(grp(ir,21,1:3,isp)**2)/(rl(ir,21,isp))**2
        enddo
        write(106,*)
        write(1106,*)
        !        isp=1
        !        do ilm=1,nlm
        !        write(6,"(i5,100d12.4)") ilm, sum(vl(:,ilm,isp)), sum(abs(vl(:,ilm,isp)))
        !        enddo
     endif
     !$$$cccccccccccccccccccccc

     !! == rho*exc, rho*vxc ==
     rmu=0d0
     rep=0d0
     do i = 1, nsp
        do ilm = 1, nlm
           do ir = 1, nr
              weight = ri(ir)**2*rwgt(ir)
              rmu(i) = rmu(i) + rl(ir,ilm,i)*vl(ir,ilm,i)*weight
              rep(i) = rep(i) + rl(ir,ilm,i)*excl(ir,ilm)*weight
           enddo
        enddo
        if(ipr>41) then
           if ( i == 1) print 7256, rmu(i),rep(i)
           if ( i == 2) print 7266, rmu(i),rep(i), rmu(1)+rmu(2),rep(1)+rep(2)
        endif
7256    format(' vxc_gga  : rho*vxc=',f11.6,' rho*exc=',f11.6,a)
7266    format('    spin 2:         ',f11.6,'         ',f11.6/ &
             '     total:         ',f11.6,'         ',f11.6)
     enddo
     if(plotstop) stop 'xxxxxxxxxx testing yyy new xxxxxxxx 111 xcpbe'
     call tcx('vxcnls')

     ! cccccccccccccccccccccccccccccccc
     !      print *,' sum wp=',sum(wp)/pi/4d0,1d0/sqrt(4d0*pi)
     !        do ilm =1,3
     !          do ip=1,np
     !            write(6,"(2i5,13d13.5)") ip,ilm,yl(ip,ilm),wp(ip),sum(abs(gyl(ip,ilm,1:3)))
     !          enddo
     !        enddo
     !      stop 'test end of vxcnlsxxx'
     ! ccccccccccccccccccccccccccccccc

     return
  endif
  !! = This is the end of new mode =












  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Mark's old version. Not go into hereafter when newmode=T
  ! --- If negative density, set to tol ---
  do  i = 1, nsp
     do  ip = 1, np
        do  ir = 1, nr
           !            if (rp(ir,ip,i) .le. 0d0) rp(ir,ip,i) = tol
           if (rp(ir,ip,i) <= tol) rp(ir,ip,i) = tol
        enddo
     enddo
  enddo
  ! --- Gradient of density point-wise through sphere ---
  if(debug) print *,'vxcnls: goto vxcnls 44444 lz=',lz
  lopgfl = 0
  if (lz) lopgfl = 10
  do  10  i = 1, nsp
     call gradfl(ll(nlm),nlm,nr,np,1,nr,1,lopgfl,nn,ri,yl,gyl, &
          rl(1,1,i),grp(1,1,1,i),ggrp(1,1,i))
10 enddo

  ! ccccccccccccccccccccccccccccccccc
  !      print *,'ttt: test grp=ggrp=0 tttttttttttttt'
  if(grpzerotest) then
     grp=0d0
     ggrp=0d0
  endif
  ! ccccccccccccccccccccccccccccccccc

  ! --- Potential from spherical part of density (valid for small r) ---
  if(debug) print *,'vxcnls: goto vxcnls 3333 lz=',lz
  if (lz) then
     do   i = 1, nsp
        call dpcopy(rl(1,1,i),wk(1,i),1,nr,dsqrt(4*pi))
        do   ir = 1, nr
           wk(ir,i) = wk(ir,i)*ri(ir)**2
        enddo
     enddo
     call dpzero(wk2,nr*nsp)
     call vxc0gc(nr,nsp,ri,rwgt,wk,wk2,excnl,repnl(3),rmunl(3), &
          100*lxcg)
     call dscal(nr,dsqrt(4*pi),excnl,1)
     do  6  i = 1, nsp
        if (ipr >= 30 .AND. i == 1) &
             print 725, rmunl(i+2),repnl(i+2),'  (l=0 rho)'
        if (ipr >= 30 .AND. i == 2) &
             print 726, rmunl(i+2),repnl(i+2), &
             rmunl(3)+rmunl(4),repnl(3)+repnl(4)
        call dpcopy(wk2(1,i),vxcnl(1,1,i),1,nr,dsqrt(4*pi))
6    enddo
     !       call prrmsh('l=0 vxcnl in vxcnls',ri,vxcnl,nr,nr,1)
     !       call prrmsh('rl in vxcnls',ri,rl,nr,nr,1)
     !       call prrmsh('wl in vxcnls',ri,rwgt,nr,nr,1)
  endif

  ! --- agrp, agrl = abs grad rho and its Yl-projection ---
  if(debug) print *,'vxcnls: goto vxcnls 55555'
  do    i = 1, nsp
     do    ip = 1, np
        do    ir = 1, nr
           agrp(ir,ip,i) = &
                dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
        enddo
     enddo
  enddo
  do   i = 1, nsp
     call dgemm('N','N',nr,nlm,np,1d0,agrp(1,1,i),nr,ylwp,np,0d0, &
          agrl(1,1,i),nr)
  enddo

  ! --- Do gradient in blocks nr-incn..nr, nr-2*incn..nr-incn ... ---
  ir1  = nr
  lmax = ll(nlm)
30 continue
  if(debug) print *,'vxcnls: goto vxcnls 6666'
  ir0 = max(ir1-incn,1)
  if (lz) ir0 = max(ir1-incn,2)
  nri = ir1-ir0+1
  if (nri <= 0) goto 32
  vhold = (vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2

  ! ...   Gradient-corrected Vxc for points between ir0 and ir1
  !      call defrr(oagrp,nri*np*(3*nsp-2))
  !      call defrr(ogagrp,nri*np*nsp*3)
  !      call defrr(ovxcp,nri*np*nsp)
  !      call defrr(oexcp,nri*np*nsp)
  !      call defrr(owk,nr)
  !      call defrr(owkl,nri*nlm*nsp)
  !      call xxcnls(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,ylwp,
  !     .wp,w(owkl),rp,ggrp,grp,agrl,w(oagrp),w(ogagrp),rl,rwgt,lcut,
  !     .w(ovxcp),w(oexcp),vxcnl,excnl,sumnl)
  if(debug) print *,'vxcnls: goto vxcnls 77777a',lxcg,lmax,ir0,ir1,nr
  if(debug) print *,'vxcnls: goto vxcnls 77777b',np,nlm,nn,lcut

  call xxcnls2(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,ylwp &
       ,wp,rp,ggrp,grp,agrl,rl,rwgt,lcut, &
       vxcnl,excnl,sumnl)

  if(debug) print *,'vxcnls: goto vxcnls 7777711111'
  !      stop 'xxxxxxxxxxxxxxxxxx test end'

  !      call rlse(oagrp)
  ir1 = ir0-1

  ! ... Check rmu to determine largest lmax next pass
  do  34  l = lmax, 0, -1
     lmx = l
     if (dabs(sumnl(l)) > 1d-7) goto 35
34 enddo
35 lmax = lmx
  if (dabs(vhold-(vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2) < 1d-6 &
       .AND. dabs(vhold) > 1) goto 32
  goto 30
32 continue
  if(debug) print *,'vxcnls: goto vxcnls 88888'

  ! --- Add nonlocal vxc into vl ----
  call daxpy(nr*nlm*nsp,1d0,vxcnl,1,vl,1)
  if (lz) then
     do  66  i = 1, nsp
        vl(1,1,i) = (vl(2,1,i)*ri(3)-vl(3,1,i)*ri(2))/(ri(3)-ri(2))
        jx = 1
        call polint(ri(2),vl(2,1,i),nr-1,nn,ri,0d0,0,jx,vl(1,1,i),vhold)
        if (ipr >= 50 .AND. dabs(vhold) > dabs(vl(1,1,i)/100)) &
             print 345, vl(1,1,i), vhold/vl(1,1,i)*100
345     format(' vxcnls (warning): expect error in V at origin:', &
             'V=',1pe10.3,' est err=',0pf7.1,'%')
66   enddo
  endif
  !     call prrmsh('nlocal v(l=0)',ri,vxcnl,nr*nlm,nr,1)
  if(debug) print *, 'vxcnls: goto vxcnls 99999'

  ! --- Nonlocal rho*exc, rho*vxc ---
  do  60  i = 1, nsp
     do   ilm = 1, nlm
        do   ir = 1, nr
           weight = ri(ir)**2*rwgt(ir)
           rmunl(i) = rmunl(i) + rl(ir,ilm,i)*vxcnl(ir,ilm,i)*weight
           repnl(i) = repnl(i) + rl(ir,ilm,i)*excnl(ir,ilm)*weight
        enddo
     enddo
     if (ipr >-1 .AND. i == 1) print 725, rmunl(i),repnl(i)
     if (ipr >-1 .AND. i == 2) print 726, rmunl(i),repnl(i), &
          rmunl(1)+rmunl(2),repnl(1)+repnl(2)
725  format(' vxcnls: nlc rmu=',f11.6,'  rep=',f11.6,a)
726  format(' spin 2:         ',f11.6,'      ',f11.6/ &
          '  total:         ',f11.6,'      ',f11.6)
     rep(i) = rep(i) + repnl(i)
     rmu(i) = rmu(i) + rmunl(i)
60 enddo
  ! --- Print out rmu by angular momentum ---
  if (ipr >= 35) then
     lmax = ll(nlm)
     do  42  i = 1, nsp
        do    l = 0, lmax
           sumnl(l) = 0d0
        enddo
        do   ilm = 1, nlm
           l = ll(ilm)
           do   ir = 1, nr
              sumnl(l) =sumnl(l)+ri(ir)**2*rl(ir,ilm,i)*vxcnl(ir,ilm,i)*rwgt(ir)
           enddo
        enddo
        if (i == 1) print 341, (sumnl(l),l=0,lmax)
        if (i == 2) print 342, (sumnl(l),l=0,lmax)
341     format(' rvnlc by L: ',f12.6,4f10.6:/(15x,4f10.6))
342     format('     spin 2: ',f12.6,4f10.6:/(15x,4f10.6))
42   enddo
  endif
  !      if (ipr .ge. 40) print 887,
  !     .  vl(1,1,1), vxcnl(1,1,1), vl(nr,1,1), vxcnl(nr,1,1)
  !  887 format(' V_0(0)=',f15.6,'  nloc VXC_0(0)=',f12.6/
  !     .       ' V_0(R)=',f15.6,'  nloc VXC_0(R)=',f12.6)

  !     call poppr

  !$$$ccccccccccccccccccccccc
  if(plotstop) then
     isp=1
     do ir=1,nr
        !         write(6,"(i5,4d12.4,3x,4d12.4,3x,4d12.4)") ir,
        !     & vl(ir,1:4,isp), vxcnl(ir,1:4,isp),vl(ir,1:4,isp)+vxcnl(ir,1:4,isp)
        !     & vl(ir,21,isp), vxcnl(ir,21,isp),vl(ir,21,isp)+vxcnl(ir,21,isp)
        write(107,"(i5,e12.4,2x,4(3e12.4,2x))") ir,ri(ir), &
             vl(ir,1,isp),  vxcnl(ir,1 ,isp),vl(ir,1,isp) -vxcnl(ir,1 ,isp), &
             vl(ir,11,isp), vxcnl(ir,11,isp),vl(ir,11,isp)-vxcnl(ir,11,isp), &
             vl(ir,21,isp), vxcnl(ir,21,isp),vl(ir,21,isp)-vxcnl(ir,21,isp)
     enddo
  endif
  if(plotstop) then
     isp=1
     do ilm=1,nlm
        write(6,"(i5,100d12.4)") ilm, sum(vl(:,ilm,isp)), sum(abs(vl(:,ilm,isp)))
     enddo
  endif
  !$$$cccccccccccccccccccccc

  ! -- rho*exc, rho*vxc ---
  do i = 1, nsp
     if ( i == 1) print 7251, rmu(i),rep(i)
     if ( i == 2) print 7261, rmu(i),rep(i),rmu(1)+rmu(2),rep(1)+rep(2)
  enddo
7251 format(' origmode vxcnls:     rmu=',f11.6,'  rep=',f11.6,a)
7261 format(' origmode spin 2:         ',f11.6,'      ',f11.6/ &
       '  total:         ',f11.6,'      ',f11.6)
  !      if(plotstop) stop 'xxxxxxxxxx testing old yyy xcfun xxxxxxxxxxx'
  call tcx('vxcnls')
end subroutine vxcnls




!      subroutine xxcnls(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,
!     .ylwp,wp,wkl,rp,ggrp,grp,agrl,agrp,gagrp,rl,rwgt,lcut,vxcp,excp,
!     .vxcnl,excnl,sumnl)
subroutine xxcnls2(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,ylwp &
     ,wp,rp,ggrp,grp,agrl,rl,rwgt,lcut, &! & fixed bug (remove wkl) 11July2010.takao
     vxcnl,excnl,sumnl)
  !     implicit none
  integer :: ir0,ir1,nr,np,nlm,nsp,nn,lmax,lcut,lxcg
  double precision :: rp(nr,np,nsp),ggrp(nr,np,nsp),ri(1),wp(np), &
       grp(nr,np,3,nsp),agrp(ir0:ir1,np,3*nsp-2),vxcnl(nr,nlm,nsp), &
       excnl(nr,nlm),agrl(nr,nlm,nsp),gagrp(ir0:ir1,np,3,nsp), &
       !     .yl(np,nlm),ylwp(np,nlm),gyl(np,nlm,3),wkl(ir0:ir1,nlm,nsp),
       yl(np,4),ylwp(np,nlm),gyl(np,nlm,3),wkl(ir0:ir1,nlm,nsp), &
       vxcp(ir0:ir1,np,nsp),excp(ir0:ir1,np),rl(nr,nlm,nsp),rwgt(nr), &
       sumnl(0:20)
  integer :: ll,ir,ip,i,nri,ilm,l,iprint,mlm,lopgfl
  logical :: debug=.false.
  !      real(8),allocatable:: wkl(:,:,:)
  !      allocate(wkl(ir0:ir1,nlm,nsp))

  if(debug) print *,'xxcnls2 00000'
  !     call pshpr(80)
  nri = ir1-ir0+1

  ! --- gagrp(store in vxcp) = grad rho . grad abs grad rho ---
  lopgfl = 0
  if (ri(ir0) == 0) lopgfl = 10
  if(debug) print *,'xxcnls2 1111'
  do  30  i = 1, nsp
     call gradfl(lmax,nlm,nr,np,ir0,ir1,0,lopgfl,nn,ri,yl,gyl, &
          agrl(1,1,i),gagrp(ir0,1,1,i),0d0)
30 enddo
  do    i  = 1, nsp
     do    ip = 1, np
        do    ir = ir0, ir1
           vxcp(ir,ip,i) = &
                gagrp(ir,ip,1,i)*grp(ir,ip,1,i) + &
                gagrp(ir,ip,2,i)*grp(ir,ip,2,i) + &
                gagrp(ir,ip,3,i)*grp(ir,ip,3,i)
        enddo
     enddo
  enddo

  ! --- store in agrp:  grad total rho . grad abs grad total rho ---
  if(debug) print *,'xxcnls2 2222'
  if (nsp == 2) then
     do    ip = 1, np
        do    ir = ir0, ir1
           agrp(ir,ip,1) = &
                (grp(ir,ip,1,1)+grp(ir,ip,1,2))* &
                (gagrp(ir,ip,1,1)+gagrp(ir,ip,1,2)) + &
                (grp(ir,ip,2,1)+grp(ir,ip,2,2))* &
                (gagrp(ir,ip,2,1)+gagrp(ir,ip,2,2)) + &
                (grp(ir,ip,3,1)+grp(ir,ip,3,2))* &
                (gagrp(ir,ip,3,1)+gagrp(ir,ip,3,2))
        enddo
     enddo
  endif

  ! --- Copy grad rho . grad abs grad rho into gagrp ---
  if(debug) print *,'xxcnls2 333333'
  do   i  = 1, nsp
     do   ip = 1, np
        do   ir = ir0, ir1
           gagrp(ir,ip,i,1) = vxcp(ir,ip,i)
        enddo
     enddo
  enddo
  if (nsp == 2) then
     do   ip = 1, np
        do   ir = ir0, ir1
           gagrp(ir,ip,3,1) = agrp(ir,ip,1)
        enddo
     enddo
  endif
  !     call px('gr.gagr',nri,nlm,1,np,ri(ir0),wp,gagrp(ir0,1,3,1),yl,wkl)


  ! --- Make agrp+,agrp- for ir0 .. ir1 ---
  if(debug) print *,'xxcnls2 44444'
  do    i = 1, nsp
     do    ip = 1, np
        do    ir = ir0, ir1
           agrp(ir,ip,i) = &
                dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
        enddo
     enddo
  enddo

  ! --- Make agrp (total rho) agrp+.agrp-  for ir0 .. ir1 ---
  if(debug) print *,'xxcnls2 55555'

  if (nsp == 2) then
     ! cccccccccccccccccccccccccccc
     !      return
     ! ccccccccccccccccccccccccccc
     do  ip = 1, np
        do  ir = ir0, ir1
           agrp(ir,ip,3) = &
                dsqrt((grp(ir,ip,1,1)+grp(ir,ip,1,2))**2 + &
                (grp(ir,ip,2,1)+grp(ir,ip,2,2))**2 + &
                (grp(ir,ip,3,1)+grp(ir,ip,3,2))**2)
           agrp(ir,ip,4) = &
                grp(ir,ip,1,1)*grp(ir,ip,1,2) + &
                grp(ir,ip,2,1)*grp(ir,ip,2,2) + &
                grp(ir,ip,3,1)*grp(ir,ip,3,2)
        enddo
     enddo
     ! cccccccccccccccccccccccccccc
     !      return
     ! ccccccccccccccccccccccccccc
     !       call px('x',nri,nlm,nsp,np,ri(ir0),wp,agrp(ir0,1,3),yl,wkl)
  endif


  ! --- Make nonlocal potential for points ir0 .. ir1 ---
  call dpzero(vxcp,nri*np*nsp)
  call dpzero(excp,nri*np)
  if(debug) print *,'xxcnls2 66666'
  do  50  ip = 1, np
     if (lxcg > 2) then
        call vxcgga(lxcg,nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp), &
             agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1), &
             ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4), &
             gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1), &
             gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp), &
             excp(ir0,ip))
     else
        if (lcut == 0) then
           call vxnloc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1), &
                agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp), &
                agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1), &
                gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1), &
                vxcp(ir0,ip,nsp),excp(ir0,ip))
        else
           call vxnlcc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1), &
                agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp), &
                agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1), &
                gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1), &
                vxcp(ir0,ip,nsp),excp(ir0,ip))
        endif
     endif
50 enddo
  if(debug) print *,'xxcnls2 77777'

  ! ... (test): yl projection of various quantities'
  !      call px('rho',nr,nlm,nsp,np,ri,wp,rp,yl,wkl)
  !      call px('ggrh',nr,nlm,nsp,np,ri,wp,ggrp,yl,wkl)
  !      call px('agrh',nri,nlm,nsp,np,ri(ir0),wp,agrp,yl,wkl)
  !      call px('gr.gagr',nri,nlm,nsp,np,ri(ir0),wp,gagrp,yl,wkl)
  !      call px('vxc',nri,nlm,nsp,np,ri(ir0),wp,vxcp,yl,wkl)
  !      call px('exc',nri,nlm,nsp,np,ri(ir0),wp,excp,yl,wkl)

  ! --- Yl-projection of vxc,exc into vxcnl,excnl ---
  mlm = (lmax+1)**2
  call dgemm('N','N',nri,mlm,np,1d0,excp(ir0,1),nri,ylwp,np,0d0, &
       wkl,nri)
  do    ilm = 1, mlm
     do    ir = ir0, ir1
        excnl(ir,ilm) = wkl(ir,ilm,1)
     enddo
  enddo
  do  62  i = 1, nsp
     call dgemm('N','N',nri,mlm,np,1d0,vxcp(ir0,1,i),nri,ylwp,np,0d0, &
          wkl,nri)
     do    ilm = 1, mlm
        do    ir = ir0, ir1
           vxcnl(ir,ilm,i) = wkl(ir,ilm,1)
        enddo
     enddo
62 enddo
  !     call prmr(nr,ri,vxcnl,nlm)

  ! --- Estimate rmu in this ir0 ir1 interval by angular momentum ---
  if(debug) print *,'xxcnls2 88888'
  do  82  i = 1, nsp
     do  83  l = 0, lmax
        sumnl(l) = 0d0
83   enddo
     do    ilm = 1, nlm
        l = ll(ilm)
        do    ir = ir0, ir1
           sumnl(l) =sumnl(l)+rl(ir,ilm,i)*vxcnl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
        enddo
     enddo
     if (iprint() >= 80) then
        if (i == 1) print 341, ri(ir0),(sumnl(l),l=0,lmax)
        if (i == 2) print 342,         (sumnl(l),l=0,lmax)
341     format(' R>',f8.6,': ',f12.6,4f10.6:/(15x,4f10.6))
342     format('     spin 2: ',f12.6,4f10.6:/(15x,4f10.6))
     endif
82 enddo
  if(debug) print *,'xxcnls2 99999'
  !     call poppr
end subroutine xxcnls2






!      subroutine px(strn,nr,nlm,nsp,np,ri,wp,fp,yl,fl)
!      implicit none
!      character *(*) strn
!      integer nr,nlm,nsp,np
!      double precision ri(1),wp(np),fp(nr,np,nsp),yl(np,nlm),
!     .  fl(nr,nlm,nsp)

!      call dpzero(fl,nr*nlm*nsp)
!      call fp2yl(nr,nlm,nsp,np,wp,fp,yl,fl)
!      print *, fl(1,1,1), fl(1,1,2)
!      print *, fl(591,1,1), fl(591,1,2)
!      print *, strn
!      call prmr(nr,ri,fl,nlm)

!      end

!      subroutine prmr(nr,ri,f,nl)
!      implicit none
!      integer nr,nl,ir,j,fopna,ifi,ir0
!      double precision ri(nr),f(nr,nl)
!      character*(10) fmt
!      ifi = fopna('out',19,0)
!      ir0 = 1
!C ... first nonzero l=0 point ...
!      do  20  ir = 1, nr
!        ir0 = ir
!        if (dabs(f(ir,1)) .gt. 1d-12) goto 21
!   20 continue
!   21 continue

!      write(ifi,*) nr-ir0+1, nl+1
!      do  10  ir = ir0, nr
!C        write(ifi,333) ri(ir), (f(ir,3*j-2), j=1, nl)
!        write(ifi,333) ri(ir), (f(ir,j), j=1, nl)
!C 333   format(f12.7,(7g18.10:/12x))
!  333   format(f12.9,(9g18.10:/12x))
!   10 continue
!      call fclose(ifi)
!      print *, 'prmr:'
!      pause
!      end


