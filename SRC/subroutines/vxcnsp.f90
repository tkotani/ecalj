subroutine vxcnsp(isw,a,ri,nr,rwgt,nlm,nsp,rl,lxcfun,rc, &
     focexc,focex,focec,focvxc,reps,repsx,repsc,rmu,vl,fl,qs)
  use m_lgunit,only:stdo
  !- Add vxc to potential in sphere and make integrals rho*vxc,rho*exc
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isw   :1s digit
  !i         : 1, make fl, otherwise not
  !i         :10s digit
  !i         : 0 calculated LDA for density as is.
  !i         : 1 reset points of negative density to positive density
  !i         : 2 for any point where rho<0 or rho_isp<0, zero potential
  !i   ri    :radial mesh points
  !i   nr    :number of radial mesh points
  !i   rwgt  :radial mesh weights
  !i   nlm   :L-cutoff for density expansion
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   rl    :full charge density * r**2
  !i   lxcfun:specifies exchange-correlation functional
  !!          1:VWN
  !!          2:barth-hedin
  !!         103:GGA-PBE
  !!         lpert=T  is not supported now.
  !!
  !o Outputs
  !o   focex :(100s digit lxcfun): integral rc * vx
  !o   focec :(100s digit lxcfun): integral rc * vc
  !o   focexc:(100s digit lxcfun): integral rc * vxc
  !o   focvxc:(100s digit lxcfun): integral rc * (dvxc/drho * rho)
  !o         : for these integrals, see Remarks
  !o   reps  :integral rl*exc
  !o   repsx :integral rl*ex
  !o   repsc :integral rl*ec
  !o   rmu   :integral rl*vxc
  !o   vl    :local exchange-correlation potential added to vl
  !o   fl    :fl(:,:,1:nsp) = vxc by L
  !o         :fl(:,:,1+nsp) = exc by L
  !o         :fl is made only if isw is nonzero
  !o   qs    :spin-resolved charge
  !l Local variables
  !l   yl    :YL(rp)
  !l   lxcf  :LD part of LDA+GGA functional
  !r Remarks
  !r   Perturbation treatment:
  !r     vxc[rho + rc] = vxc[rho] + rc * (dvxc/drho)
  !r   Correction to int rho * vxc:
  !r     int rho * vxc = int rho*vxc[rho] + int rho*rc*(dvxc/drho)
  !r     Correction = focvxc = int rho * rc * (dvxc/drho)
  !r   Correction to int rho * exc:
  !r     Corr = int rho * rc * (dexc/drho) = int rho * rc * (vxc-exc)/rho
  !r          = int rc * (vxc-exc) = focexc - int rc exc
  !r     Second term is not retained.
  !u Updates
  !u   06 Apr 09 Routine calls vxcnls to handle some GGA potentials
  !u   14 Mar 04 Optionally Makes fl.  New argument list
  !u   14 Jun 02 repsx and repsc (T. Miyake)
  !u    8 Feb 02 focex and focec (T. Miyake)
  !u   13 Jun 00 spin polarized
  !u    3 Jun 00 Adapted from old FP vxcnsp; pert. corr from nfp vxcdnsp.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: isw,nr,nlm,lxcfun,nsp
  double precision :: ri(nr),reps(2),rmu(2),repsx(2),repsc(2), &
       rwgt(nr),rc(nr),rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,1+nsp)
  double precision :: focexc(2),focex(2),focec(2),focvxc(2)
  double precision :: qs(2)
  ! ... Local parameters
  integer :: nnn,nlmx
  parameter(nnn=122,nlmx=64)
  logical :: lpert
  integer :: ipr,ll,lmax,np,nph,nth,nxl(0:7),lxcf,lxcg
  double precision :: p(3,nnn),wp(nnn),p2(nnn,3),r2(nnn),fpi
  !     double precision yl(nnn*nlmx)
  !      real(8), allocatable :: ylwp(:)
  real(8), allocatable :: yl(:),rp(:,:,:),gyl(:,:),grp(:,:),agrp(:), &
       ggrp(:),agrl(:),vxcnl(:),excnl(:)

  data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
  logical:: debug=.false.
  real(8):: a !takao used to get derivetive of ri. See vxcnls.F
  real(8),allocatable::vlxc(:,:,:)

  !      stdo = lgunit(1)
  call getpr(ipr)
  lxcf = mod(lxcfun,100)
  lxcg = mod(lxcfun/100,100)
  lpert = lxcfun .ge. 10000
  fpi = 16d0*datan(1d0)

  ! ... Create angular integration mesh
  lmax = ll(nlm)
  !     if (lxcg .ne. 0) lmax = lmax+1
  if (lmax > 6) then
     nth = 2*lmax+2
     nph = 0
  else
     nth = nxl(lmax)
     nph = 0
  endif
  call fpiint(nth,nph,np,p,wp)
  if(ipr >= 30) write(stdo,200) nth,nph,np,nr
200 format(' mesh:   nth,nph=',2i4,'   gives',i4,'  angular points,', &
       '   nrad=',i4)
  if (np > nnn) call rxi('vxcnsp: increase nnn, need',np)
  if ((lmax+2)**2*np > nlmx*nnn) call rx('vxcnsp: increase nlm')

  rmu(2)  = 0
  reps(2) = 0
  repsx(2) = 0
  repsc(2) = 0

  !!== Scale rl to true density ==
  call vxcns3(nr,nlm,nsp,ri,rl,1)
  if (lpert) then
     call vxcns3(nr,1,1,ri,rc,1)
     call dscal(nr,1/fpi,rc,1)
  endif

  p2=transpose(p)
  !!== XC part ==
  allocate(vlxc(nr,nlm,nsp))
  !!=== LDA exchange-correlation potential ===
  if(lxcg==0) then
     allocate (yl(np*nlm))
     call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax,np,yl,r2)
     vlxc=0d0
     call vxcns2(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,rc,lxcf, &
          lpert,focexc,focex,focec,focvxc,reps,repsx,repsc,rmu, &
          vlxc,fl,qs)
     deallocate (yl)
  else
     !!=== GGA ===
     !!    GGA case:  Need both YL to lmax+1 to make grad YL
     allocate (yl(np*(lmax+2)**2)) !,ylwp(np*(lmax+2)**2))
     call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax+1,np,yl,r2)
     allocate (gyl(np*nlm,3))
     call ropylg(1,lmax,nlm,np,np,p2(1,1),p2(1,2),p2(1,3),r2,yl,gyl)
     !        allocate(grp(nr*np*nsp,3),agrp(nr*np*nsp),ggrp(nr*np*nsp))
     !        allocate(rp(nr,np,nsp))
     !        allocate(agrl(nr*nlm*nsp))!,vxcnl(nr*nlm*nsp),excnl(nr*nlm*nsp))
     if(debug) print *,'goto vxcnls'
     call vxcnls(a,ri,0,nr,np,nlm,nsp, yl,gyl,rwgt,wp,rl,lxcg, vlxc,reps,rmu)
     deallocate(yl,gyl) !,grp,agrp,ggrp,agrl)!,vxcnl,excnl,rp)
  endif
  vl=vl+vlxc
  deallocate(vlxc)

  ! ... Undo the rl scaling
  call vxcns3(nr,nlm,nsp,ri,rl,0)
  if (lpert) then
     call vxcns3(nr,1,1,ri,rc,0)
     call dscal(nr,fpi,rc,1)
  endif

end subroutine vxcnsp

subroutine vxcns2(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,rc,lxcf, &
     lpert,focexc,focex,focec,focvxc,rep,repx,repc,rmu,vl,fl,qs)
  use m_lgunit,only:stdo
  !- Make vxc, rep and rmu in sphere, for local XC functional.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isw   :1s digit
  !i         : 1, make fl, otherwise not
  !i         :10s digit
  !i         : 0 calculated LDA for density as is.
  !i         : 1 reset points of negative density to positive density
  !i         : 2 for any point where rho<0 or rho_isp<0, zero potential
  !i   ri    :radial mesh points
  !i   nr    :number of radial mesh points
  !i   rwgt  :radial mesh weights
  !i   np    :number of spherical mesh points
  !i   wp    :spherical mesh weights
  !i   nlm   :L-cutoff for density expansion
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   rl    :charge density
  !i   rc    :l=0 charge density*4*pi*ri**2 to be added in pert. theory
  !i   lxcf  :specifies exchange-correlation functional
  !i         :  1    Ceperly-Alder
  !i         :  2    Barth-Hedin (ASW fit)
  !i         :  3,4  LD part of PW91 and PBE
  !i   lpert :Make perturbation integrals focexc and focvxc
  !i   yl    :Spherical harmonics YL(rp) !sep2010 takao. yl is not overwritten by ylwp now.
  !o Outputs
  !o   focexc:(lpert only): integral rc * vxc
  !o   focex :(lpert only): integral rc * vx
  !o   focec :(lpert only): integral rc * vc
  !o   focvxc:(lpert only): integral rc * (dvxc/drho * rho)
  !o   rep   :integral rl*exc
  !o   repx  :integral rl*ex
  !o   repc  :integral rl*ec
  !o   rmu   :integral rl*vxc
  !o   vl    :local exchange-correlation potential added to vl
  !o   fl    :fl(:,:,1:nsp) = vxc by L
  !o         :fl(:,:,1+nsp) = exc by L
  !o         :fl is only made if 1s digit isw is nonzero
  !o   qs    :spin-resolved charge
  !l Local variables
  !l   rp    :list of points on the spherical mesh
  !l   exc   :exchange density on (nr,np) mesh
  !l   vxc   :exchange potential on (nr,np) mesh
  !r Remarks
  !r   For perturbation treatment, take numerical derivatives
  !r   df/dr = d/dr (vxc*r**alfa) instead of d/dr (vxc) because
  !r   f is nearly linear for alpha=2/3.
  !r
  !r   In the spin polarized case, the perturbation density rc is not
  !r   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
  !r   we use dm=0 and thus drho1 = drho2 = drho/2; thus
  !r     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
  !u Updates
  !u   20 Nov 09 New 20s digit for isw
  !u   02 Apr 09 New option (10's digit isw)
  !u   14 Jun 02 repx and repc (T. Miyake)
  !u    8 Feb 02 focex and focec (T. Miyake)
  !u   13 Jun 00 spin polarized
  !u    3 Jun 00 adapted from old FP code
  !u    10.07.96 dln: modified for Perdew GGA
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  logical :: lpert
  integer :: nr,nlm,nsp,lxcf,np,isw
  double precision :: ri(1),yl(np,nlm),wp(np),rep(2),repx(2),repc(2), &
       rmu(2),rwgt(nr),focexc(2),focex(2),focec(2),focvxc(2), &
       rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,nsp+1),rc(nr) &
       ,ylwp(np,nlm)  !to avoid confusion. Sep2010 takao
  ! ... Local parameters
  integer :: ilm,ip,ipr,ir,i,n1,ineg(2),isw1
  ! ino size 41->42      double precision weight,suml(0:40),fac,dmach,alfa,sum(2)
  double precision :: weight,suml(0:41),fac,dmach,alfa,sum(2)
  double precision :: rhot,f1,f2,dfdr,dvdr,f
  double precision :: qs(2),tol
  real(8), allocatable :: rp(:,:),rps(:,:,:),exc(:,:), &
       excx(:,:),excc(:,:)
  real(8), allocatable :: vxc(:,:,:),vx(:,:,:), &
       vc(:,:,:),vxc2(:,:,:),vx2(:,:,:),vc2(:,:,:)
  real(8),allocatable:: vlxc(:,:,:)


  allocate (rp(nr,np),rps(nr,np,nsp),exc(nr,np), &
       excx(nr,np),excc(nr,np))
  allocate (vxc(nr,np,nsp),vx(nr,np,nsp),vc(nr,np,nsp), &
       vxc2(nr,np,nsp),vx2(nr,np,nsp),vc2(nr,np,nsp))

  !     stdo = lgunit(1)
  call getpr(ipr)

  !     'small' density
  tol = 1d-15
  isw1 = mod(isw/10,10)

  !     for numerical differentiation
  fac = dmach(1)**(1d0/3d0)
  alfa = 2d0/3d0
  n1 = nr*np
  !     call prrmsh('rl',rl,rl,nr,nr,nlm*2)

  ! --- Generate density point-wise through sphere ---
  call dpzero(rp,n1)
  ineg(1) = 0
  ineg(2) = 0
  do  i = 1, nsp
     call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0, &
          rps(1,1,i),nr)

     ! cccccccccccccccccccccccccc
     do  ip = 1, np
        do  ir = 1, nr
           if (rps(ir,ip,i) < 0d0) then
              !                print *,'xxx vvv', ir,ip,rps(ir,ip,i)
              ineg(i) = ineg(i) + 1
           endif
        enddo
     enddo
     ! cccccccccccccccccccccccccc

     !   ... Counts number of points with negative density; reset to pos.
     if (isw1 == 1) then
        do  ip = 1, np
           do  ir = 1, nr
              if (rps(ir,ip,i) < 0d0) then
                 ineg(i) = ineg(i) + 1
                 rps(ir,ip,i) = tol
              endif
           enddo
        enddo
     endif

     call daxpy(n1,1d0,rps(1,1,i),1,rp,1)
  enddo


  if (ineg(1)+ineg(2) /= 0) then
     write(6,"(a,i10,a,2i10)") ' vxcns2 (warning): nr*np=',nr*np, &
          "  negative density # of points=", ineg(1:nsp)
     !        call info5(5,0,0,' vxcns2 (warning):  negative density,'
     !     .  //' %i points %?#n==2#(spin 1) %i points (spin 2)##',
     !     .  ineg(1),nsp,ineg(2),0,0)
  endif

  ! --- Perturbation treatment: df/dr in vxc2 ---
  if (lpert) then

     ! --- Stop if GGA ---
     if (lxcf > 2) &
          call rx('vxcnsp: Perturbation treatment'// &
          ' is not implemented for GGA')

     !       Add rp*fac/2 into rp+, rp- and fac*rp into rp
     if (nsp == 2) then
        do  i = 1, nsp
           call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
        enddo
     endif
     call dscal(n1,1+fac,rp,1)
     !       Exchange potential at rp+drho
     do  i = 1, nsp
        call evxcv(rp,rps(1,1,i),n1,1,lxcf, &
             exc,excx,excc,vxc2(1,1,i),vx2(1,1,i),vc2(1,1,i))
     enddo
     !       Restore rp,rps; also add -drho to rp and -drho/2 to rps
     call dscal(n1,1/(1+fac),rp,1)
     if (nsp == 2) then
        do  i = 1, nsp
           call daxpy(n1,-fac,rp,1,rps(1,1,i),1)
        enddo
     endif
     call dscal(n1,(1-fac),rp,1)
     !       Exchange potential at rp-drho
     do  i = 1, nsp
        call evxcv(rp,rps(1,1,i),n1,1,lxcf, &
             exc,excx,excc,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
     enddo
     !       Restore rp,rps
     call dscal(n1,1/(1-fac),rp,1)
     if (nsp == 2) then
        do  i = 1, nsp
           call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
        enddo
     endif

     do  i = 1, nsp
        do  ip = 1, np
           do  ir = 1, nr
              rhot = rp(ir,ip)
              if (rhot > 0) then
                 f1 = vxc (ir,ip,i)*(rhot*(1-fac))**alfa
                 f2 = vxc2(ir,ip,i)*(rhot*(1+fac))**alfa
                 dfdr = (f2-f1)/(2d0*fac*rhot)
                 vxc2(ir,ip,i) = dfdr
              else
                 vxc2(ir,ip,i) = 0
              endif
           enddo
        enddo
     enddo
  else
     call dpzero(excx,n1)
     call dpzero(excc,n1)
  endif

  ! --- vxc, exc for unperturbed density ---
  !     print *, '!!'; lxcf=1
  if (lxcf > 2) then
     call evxcp(rps,rps(1,1,nsp),n1,nsp,lxcf,excx,excc,exc, &
          vx(1,1,1),vx(1,1,2),vc(1,1,1),vc(1,1,2), &
          vxc,vxc(1,1,nsp))
  else
     do  i = 1, nsp
        call evxcv(rp,rps(1,1,i),n1,nsp,lxcf, &
             exc,excx,excc,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
     enddo
  endif
  if (isw1 == 2) then
     do  i = 1, nsp
        do  ip = 1, np
           do  ir = 1, nr
              if (rp(ir,ip) <= 0 .OR. rps(ir,ip,i) <= 0) then
                 vxc(ir,ip,i) = 0
                 vx(ir,ip,i) = 0
                 vc(ir,ip,i) = 0
              endif
           enddo
        enddo
     enddo
  endif

  ! cccccccccccccccccccccc
  !      print *,' np=',np
  !      do ir=1,nr
  !        write(6,"(i5,5d13.5)") ir,exc(ir,np-4:np)
  !      enddo
  !      stop 'xxxxxxxxxxxx test exc vxcnsp'
  ! cccccccccccccccccccccc

  ! --- Integrals rep, rmu ---
  !      rpneg = 0
  !      do  24  i = 1, nsp
  !        qs(i) = 0d0
  !        rep(i) = 0d0
  !        repx(i) = 0d0
  !        repc(i) = 0d0
  !        rmu(i) = 0d0
  !        do  22  ip = 1, np
  !        do  22  ir = 1, nr
  !        rpneg = min(rpneg,rps(ir,ip,i))
  !        weight = (ri(ir)**2*rwgt(ir))*wp(ip)
  !        qs(i)  = qs(i)  + rps(ir,ip,i)*weight
  !        rep(i) = rep(i) + exc(ir,ip)*rps(ir,ip,i)*weight
  !        repx(i) = repx(i) + excx(ir,ip)*rps(ir,ip,i)*weight
  !        repc(i) = repc(i) + excc(ir,ip)*rps(ir,ip,i)*weight
  !   22   rmu(i) = rmu(i) + vxc(ir,ip,i)*rps(ir,ip,i)*weight
  !        if (ipr.ge.30 .and. i.eq.1) write(stdo,725) rmu(i),rep(i),qs(i)
  !        if (ipr.ge.30 .and. i.eq.2) write(stdo,726) rmu(i),rep(i),qs(i),
  !     .    rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2)
  !  725   format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
  !  726   format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/
  !     .         '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)
  !   24 continue
  !      if (rpneg .lt. 0 .and. ipr .ge. 20) write(stdo,727) rpneg
  !  727 format(' vxcnsp (warning): negative rho: min val = ',1pe10.2)
  call vxcns4(0,ri,nr,rwgt,np,wp,nsp,rps,exc,excx,excc,vxc, &
       rep,repx,repc,rmu,qs)

  ! --- Add perturbation to vxc ---
  !     Integrals focexc = int rc vxc, focvxc= rc * dvxc/dr * rhot
  if (lpert) then
     focvxc(1) = 0
     focvxc(2) = 0
     focexc(1) = 0
     focexc(2) = 0
     focex(1)  = 0
     focex(2)  = 0
     focec(1)  = 0
     focec(2)  = 0
     do  i  = 1, nsp
        do  ip = 1, np
           do  ir = 1, nr
              rhot = rp(ir,ip)
              if (rhot <= 0 .OR. rps(ir,ip,i) <= 0) then
                 vxc2(ir,ip,i) = 0
              endif
              !             Debugging
              !              if (rps(ir,ip,i) .lt. 0 .or. rhot .lt. 0) then
              !                if (vxc(ir,ip,i) .ne. 0 .or. vxc2(ir,ip,i) .ne. 0) then
              !                  print *, vxc(ir,ip,i), vxc2(ir,ip,i)
              !                  stop 'oops'
              !                endif
              !              endif
              if (rps(ir,ip,i) > 0 .AND. rhot > 0) then
                 f  = vxc(ir,ip,i)*rhot**alfa
                 dfdr = vxc2(ir,ip,i)
                 dvdr = (dfdr - alfa*f/rhot) / rhot**alfa
                 weight = (ri(ir)**2*rwgt(ir))*wp(ip) * rc(ir)
                 focvxc(i) = focvxc(i) + weight*dvdr*rps(ir,ip,i)
                 focexc(i) = focexc(i) + weight*vxc(ir,ip,i)/nsp
                 focex(i)  = focex(i)  + weight*vx(ir,ip,i)/nsp
                 focec(i)  = focec(i)  + weight*vc(ir,ip,i)/nsp
                 vxc(ir,ip,i) = vxc(ir,ip,i) + dvdr*rc(ir)
              endif
           enddo
        enddo
     enddo
  endif

  ! --- Scale yl by wp for fast multiplication --- ---
  do   ilm = 1, nlm
     do   ip = 1, np
        ylwp(ip,ilm) = yl(ip,ilm)*wp(ip) ! it was yl=yl*wp
     enddo
  enddo
  ! --- Add Yl-projection of vxc into vl ---
  do  30  i = 1, nsp
     call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,ylwp,np,1d0, &
          vl(1,1,i),nr)
30 enddo


  !$$$ccccccccccccccccccccccc
  !      allocate(vlxc(nr,nlm,nsp))
  !      do   i = 1, nsp
  !        call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,0d0,
  !     .  vlxc(1,1,i),nr)
  !      enddo
  !      isp=1
  !      do ir=1,nr
  !       write(6,"(i5,8d12.4)") ir, vxc(ir,1:8,isp) !vlxc(ir,1:8,isp)
  !      enddo
  !$$$      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxx test end vxcnsp 123'
  !$$$cccccccccccccccccccccc

  ! --- Optionally make l-resolved vxc_L,exc_L ---
  if (mod(isw,10) == 1) then
     do  i = 1, nsp
        call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,ylwp,np,0d0, fl(1,1,i),nr)
     enddo
     call dgemm('N','N',nr,nlm,np,1d0,exc,nr,ylwp,np,0d0, & ! & yl-->ylwp
     fl(1,1,1+nsp),nr)
  endif

  ! --- Print out int (rl*vl) resolved by l ---
  if (ipr > 30) then
     call vxcns5(0,ipr,'rho*vtot',nlm,nsp,nr,ri,rwgt,rl,vl,suml,sum)
     !      lmax = ll(nlm)
     !      do  42  i = 1, nsp
     !      do  43  l = 0, lmax
     !   43 suml(l) = 0d0
     !      do  40  ilm = 1, nlm
     !      l = ll(ilm)
     !      do  40  ir = 1, nr
     !   40 suml(l) = suml(l) + rl(ir,ilm,i)*vl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
     !      if (i .eq. 1) write(stdo,341) (suml(l),l = 0,lmax)
     !      if (i .eq. 2) write(stdo,342) (suml(l),l = 0,lmax)
     !  341 format(' rho*vxc by l: ',f13.6,4f10.6:/(18x,4f10.6))
     !  342 format('       spin 2: ',f13.6,4f10.6:/(18x,4f10.6))
     !   42 continue
  endif
  deallocate (rp,rps,exc,excx,excc)
  deallocate (vxc,vx,vc,vxc2,vx2,vc2)
end subroutine vxcns2

subroutine vxcns3(nr,nlm,nsp,ri,rl,isgn)
  !- Scales rho by r**2, or undoes scaling
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nr    :number of radial mesh points
  !i   nlm   :L-cutoff for density expansion
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   ri    :radial mesh points
  !i   isgn  :1, scale rl by 1/r**2
  !i         :else scale rl by r**2
  !o Outputs
  !o   rl   :rl is scaled by r**2 or 1/r**2, depending on isgn
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nr,nlm,nsp,isgn,i,ilm,ir
  double precision :: rl(nr,nlm,nsp),ri(nr),rho2,rho3

  ! --- Scale rho by 1/r**2 ---
  if (isgn == 1) then
     do    i = 1, nsp
        do    ilm = 1, nlm
           rl(1,ilm,i) = 0d0
           do    ir = 2, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)/ri(ir)**2
           enddo
        enddo
     enddo
     !  ...  Extrapolate rho to origin
     do  20  i = 1, nsp
        rho2 = rl(2,1,i)
        rho3 = rl(3,1,i)
        rl(1,1,i) = (rho2*ri(3)-rho3*ri(2))/(ri(3)-ri(2))
20   enddo
  else
     do   i = 1, nsp
        do   ilm = 1, nlm
           do   ir = 1, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)*ri(ir)**2
           enddo
        enddo
     enddo
  endif
end subroutine vxcns3

subroutine vxcns4(isw,ri,nr,rwgt,np,wp,nsp,rp,exc,excx,excc,vxc, &
     rep,repx,repc,rmu,qs)
  use m_lgunit,only:stdo
  !- Integrals of vxc, reps and rmu in sphere.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isw
  !i   ri    :radial mesh points
  !i   nr    :number of radial mesh points
  !i   rwgt  :radial mesh weights
  !i   np    :number of spherical mesh points
  !i   wp    :spherical mesh weights
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !l   rp    :list of points on the spherical mesh
  !l   exc   :exchange-correlation energy density on (nr,np) mesh
  !l   excx  :exchange             energy density on (nr,np) mesh
  !l   excc  :correlation          energy density on (nr,np) mesh
  !l   exc   :exchange-correlation energy density on (nr,np) mesh
  !l   vxc   :exchange correlation potential on (nr,np) mesh
  !o Outputs
  !o   rep   :integral rl*exc
  !o   repx  :integral rl*ex
  !o   repc  :integral rl*ec
  !o   rmu   :integral rl*vxc
  !o   qs    :spin-resolved charge
  !l Local variables
  !l   rpneg :number of points at which density < 0
  !r Remarks
  !u Updates
  !u   12 Mar 04 created from vxcns2
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nr,nsp,np,isw
  double precision :: ri(1),wp(np),rep(2),repx(2),repc(2),rmu(2), &
       rwgt(nr),rp(nr,np,nsp),exc(nr,np),excx(nr,np),excc(nr,np), &
       vxc(nr,np,nsp)
  ! ... Local parameters
  integer :: ip,ipr,ir,i,isw1
  double precision :: weight,rpneg,qs(2)

  !      stdo = lgunit(1)
  call getpr(ipr)
  isw1 = mod(isw/10,10)

  ! --- Integrals reps, rmu ---
  rpneg = 0
  do  24  i = 1, nsp
     qs(i) = 0d0
     rep(i) = 0d0
     repx(i) = 0d0
     repc(i) = 0d0
     rmu(i) = 0d0
     do  ip = 1, np
        do  ir = 1, nr
           rpneg = min(rpneg,rp(ir,ip,i))
           weight = (ri(ir)**2*rwgt(ir))*wp(ip)
           qs(i)  = qs(i)  + rp(ir,ip,i)*weight
           if (isw1 == 2) then
              if (rp(ir,ip,1)+rp(ir,ip,nsp) <= 0 .OR. rp(ir,ip,i) < 0) then
                 exc(ir,ip) = 0
                 excx(ir,ip) = 0
                 excc(ir,ip) = 0
                 vxc(ir,ip,i) = 0
              endif
           endif
           !         Debugging
           !          if (rp(ir,ip,1) .lt. 0 .or. rp(ir,ip,nsp) .lt. 0) then
           !            if (vxc(ir,ip,i) .ne. 0) then
           !              print *, vxc(ir,ip,i)
           !              stop 'oops'
           !            endif
           !          endif
           rep(i) = rep(i) + exc(ir,ip)*rp(ir,ip,i)*weight
           repx(i) = repx(i) + excx(ir,ip)*rp(ir,ip,i)*weight
           repc(i) = repc(i) + excc(ir,ip)*rp(ir,ip,i)*weight
           rmu(i) = rmu(i) + vxc(ir,ip,i)*rp(ir,ip,i)*weight
        enddo
     enddo
     if (ipr >= 30 .AND. i == 1) write(stdo,725) rmu(i),rep(i),qs(i)
     if (ipr >= 30 .AND. i == 2) write(stdo,726) rmu(i),rep(i),qs(i), &
          rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2)
725  format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
726  format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/ &
          '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)
24 enddo
  if (rpneg < 0 .AND. ipr >= 10) write(stdo,727) rpneg
727 format(' vxcnsp (warning): negative rho: min val = ',1pe10.2)

end subroutine vxcns4

subroutine vxcns5(isw,ipr,strn,nlml,nsp,nr,ri,rwgt,rl,vl,suml,sum)
  use m_lgunit,only:stdo
  !- Integrals of rl*vl resolved by l
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isw   :0 rl is true density
  !i         :1 rl is true density * ri**2
  !i   strn  :used in printout
  !i   ri    :radial mesh points
  !i   rwgt  :radial mesh weights
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nlml  :number of Ylm's in rl and vl
  !o Outputs
  !o   suml  :l-and spin-resolved integrals of rl*vl
  !o   sum   :spin-resolved integrals of rl*vl
  !r Remarks
  !u Updates
  !u   12 Mar 04 created from vxcns2
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character strn*(*)
  integer :: isw,ipr,nlml,nsp,nr
  double precision :: ri(nr),rwgt(nr),rl(nr,nlml,nsp),vl(nr,nlml,nsp), &
       suml(0:20,2),sum(2)
  ! ... Local parameters
  integer :: ir,lmax,ilm,isp,l,ll
  double precision :: dsum,dot3

  !      stdo = lgunit(1)
  lmax = ll(nlml)

  do   isp = 1, nsp
     call dpzero(suml(0,isp),lmax+1)
     do  ilm = 1, nlml
        l = ll(ilm)
        if (isw == 1) then
           suml(l,isp) = suml(l,isp) + &
                dot3(nr,rl(1,ilm,isp),vl(1,ilm,isp),rwgt)
        else
           do  ir = 1, nr
              suml(l,isp) = suml(l,isp) + &
                   rl(ir,ilm,isp)*vl(ir,ilm,isp)*ri(ir)**2*rwgt(ir)
           enddo
        endif
     enddo

     sum(isp) = dsum(lmax+1,suml(0,isp),1)

     if (ipr > 30) then
        if (isp == 1) &
             write(stdo,341) strn,sum(isp),(suml(l,isp),l = 0,lmax)
        if (isp == 2) &
             write(stdo,342) sum(isp),(suml(l,isp),l = 0,lmax)
341     format(1x,a8,f13.6,' ... by l: ',f12.6,4f10.6:/(18x,4f10.6))
342     format('  spin 2:',f13.6,' ... by l: ',f12.6, &
             4f10.6:/(18x,4f10.6))
     endif

  enddo

end subroutine vxcns5

!      subroutine fp2yl(nr,nlm,nsp,np,wp,fp,yl,fl)
!C- Add Yl-projection of function tabulated on a mesh
!      implicit none
!      integer nr,nlm,np,ip,ilm,ir,i,nsp
!      double precision fl(nr,nlm,nsp),fp(nr,np,nsp),yl(np,nlm),wp(np),xx

!      do  20  i = 1, nsp
!      do  20  ip = 1, np
!      do  20  ilm = 1, nlm
!      xx = wp(ip)*yl(ip,ilm)
!      do  20  ir = 1, nr
!   20 fl(ir,ilm,i) = fl(ir,ilm,i) + fp(ir,ip,i)*xx
!      end

