module m_vxcatom !Vxc LDA for sites (spherical expansion)
  public vxcnsp,vxc0sp 
  private
  contains
subroutine vxcnsp(isw,a,ri,nr,rwgt,nlm,nsp,rl,lxcfun, reps,repsx,repsc,rmu,vl,fl,qs) 
  use m_lgunit,only: stdo
  use m_ropyln,only: ropyln,ropylg
!  use m_vxcfunc,only: vxcgga
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
  !ox   focex :(100s digit lxcfun): integral rc * vx
  !ox   focec :(100s digit lxcfun): integral rc * vc
  !ox   focexc:(100s digit lxcfun): integral rc * vxc
  !ox   focvxc:(100s digit lxcfun): integral rc * (dvxc/drho * rho)
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
  !rxxx   Perturbation treatment:
  !rxxx     vxc[rho + rc] = vxc[rho] + rc * (dvxc/drho)
  !rxxx   Correction to int rho * vxc:
  !rxxx     int rho * vxc = int rho*vxc[rho] + int rho*rc*(dvxc/drho)
  !rxxx     Correction = focvxc = int rho * rc * (dvxc/drho)
  !rxxx   Correction to int rho * exc:
  !rxxx     Corr = int rho * rc * (dexc/drho) = int rho * rc * (vxc-exc)/rho
  !rxxx          = int rc * (vxc-exc) = focexc - int rc exc
  !r     Second term is not retained.
  !u Updates
  !u   06 Apr 09 Routine calls vxcnls to handle some GGA potentials
  !u   14 Mar 04 Optionally Makes fl.  New argument list
  !u   14 Jun 02 repsx and repsc (T. Miyake)
  !u    8 Feb 02 focex and focec (T. Miyake)
  !u   13 Jun 00 spin polarized
  !u    3 Jun 00 Adapted from old FP vxcnsp; pert. corr from nfp vxcdnsp.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: isw,nr,nlm,lxcfun,nsp
  double precision :: ri(nr),reps(nsp),rmu(nsp),repsx(nsp),repsc(nsp), &
       rwgt(nr),rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,1+nsp) !,rc(nr)
!  double precision :: focexc(nsp),focex(nsp),focec(nsp),focvxc(nsp)
  double precision :: qs(nsp)
  integer :: nnn,nlmx
  parameter(nnn=122,nlmx=64)
!  logical :: lpert
  integer :: ipr,ll,lmax,np,nph,nth,nxl(0:7),lxcf,lxcg
  double precision :: p(3,nnn),wp(nnn),p2(nnn,3),r2(nnn),fpi
  real(8), allocatable :: yl(:),rp(:,:,:),gyl(:,:),grp(:,:),agrp(:), &
       ggrp(:),agrl(:),vxcnl(:),excnl(:)
  data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
  logical:: debug=.false.
  real(8):: a !takao used to get derivetive of ri. See vxcnls.F
  real(8),allocatable::vlxc(:,:,:)
  call getpr(ipr)
  lxcf = mod(lxcfun,100)
  lxcg = mod(lxcfun/100,100) !GGA
!  lpert = .false. !lxcfun .ge. 10000
  fpi = 16d0*datan(1d0)
  ! ... Create angular integration mesh
  lmax = ll(nlm)
  if (lmax > 6) then
     nth = 2*lmax+2
     nph = 0
  else
     nth = nxl(lmax)
     nph = 0
  endif
  call fpiint(nth,nph,np,p,wp)
  if(ipr >= 30) write(stdo,200) nth,nph,np,nr
200 format(' mesh:   nth,nph=',2i4,'   gives',i4,'  angular points,','   nrad=',i4)
  if (np > nnn) call rxi('vxcnsp: increase nnn, need',np)
  if ((lmax+2)**2*np > nlmx*nnn) call rx('vxcnsp: increase nlm')
  rmu  = 0
  reps = 0
  repsx = 0
  repsc = 0
  !!== Scale rl to true density ==
  call vxcns3(nr,nlm,nsp,ri,rl,1)
!  if (lpert) then
!     call vxcns3(nr,1,1,ri,rc,1)
!     call dscal(nr,1/fpi,rc,1)
!  endif
  p2=transpose(p)
  !!== XC part ==
  allocate(vlxc(nr,nlm,nsp))
  !!=== LDA exchange-correlation potential ===
  if(lxcg==0) then
     allocate (yl(np*nlm))
     call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax,np,yl,r2)
     vlxc=0d0
     call vxcns2(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,lxcf, & 
          reps,repsx,repsc,rmu, & !,focexc,focex,focec,focvxc lpert,,rc
          vlxc,fl,qs)
     deallocate (yl)
  else
     !!    GGA case:  Need both YL to lmax+1 to make grad YL
     allocate (yl(np*(lmax+2)**2)) 
     call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax+1,np,yl,r2)
     allocate (gyl(np*nlm,3))
     call ropylg(1,lmax,nlm,np,np,p2(1,1),p2(1,2),p2(1,3),r2,yl,gyl)
     if(debug) print *,'goto vxcnls'
     call vxcnls(a,ri,nr,np,nlm,nsp, yl,gyl,rwgt,wp,rl,lxcg, vlxc,reps,rmu)
     deallocate(yl,gyl) 
  endif
  vl=vl+vlxc
  deallocate(vlxc)
  ! ... Undo the rl scaling
  call vxcns3(nr,nlm,nsp,ri,rl,0)
!  if (lpert) then
!     call vxcns3(nr,1,1,ri,rc,0)
!     call dscal(nr,fpi,rc,1)
!  endif
end subroutine vxcnsp
subroutine vxcns2(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,lxcf, & 
     rep,repx,repc,rmu,vl,fl,qs) !,focexc,focex,focec,focvxc lpert, ,rc
  use m_lgunit,only:stdo
  use m_xclda,only:evxcv,evxcp
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
  !i   lpert=F :Make perturbation integrals focexc and focvxc
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
  implicit none
!  logical :: lpert
  integer :: nr,nlm,nsp,lxcf,np,isw
  double precision :: ri(nr),yl(np,nlm),wp(np),rep(nsp),repx(nsp),repc(nsp), &
       rmu(nsp),rwgt(nr), & !,focexc(nsp),focex(nsp),focec(nsp),focvxc(nsp)
       rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,nsp+1) & !,rc(nr)
       ,ylwp(np,nlm)  !to avoid confusion. Sep2010 takao
  integer :: ilm,ip,ipr,ir,i,n1,ineg(nsp),isw1
  double precision :: weight,suml(0:41),fac,dmach,alfa,ssum(nsp)
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
  call getpr(ipr)
  tol = 1d-15
  isw1 = mod(isw/10,10)
  !     for numerical differentiation
  fac = dmach(1)**(1d0/3d0)
  alfa = 2d0/3d0
  n1 = nr*np
  ! --- Generate density point-wise through sphere ---
  call dpzero(rp,n1)
  ineg = 0
  do  i = 1, nsp
     call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0, rps(1,1,i),nr)
     do  ip = 1, np
        do  ir = 1, nr
           if (rps(ir,ip,i) < 0d0) then
              ineg(i) = ineg(i) + 1
           endif
        enddo
     enddo
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
  if (sum(ineg) /= 0) then
     write(6,"(a,i10,a,2i10)") ' vxcns2 (warning): nr*np=',nr*np, &
          "  negative density # of points=", ineg(1:nsp)
  endif
  ! --- Perturbation treatment: df/dr in vxc2 ---
  ! if (lpert) then
  !    ! --- Stop if GGA ---
  !    if (lxcf > 2) &
  !         call rx('vxcnsp: Perturbation treatment is not implemented for GGA')
  !    !       Add rp*fac/2 into rp+, rp- and fac*rp into rp
  !    if (nsp == 2) then
  !       do  i = 1, nsp
  !          call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
  !       enddo
  !    endif
  !    call dscal(n1,1+fac,rp,1)
  !    !       Exchange potential at rp+drho
  !    do  i = 1, nsp
  !       call evxcv(rp,rps(1,1,i),n1,1,lxcf, &
  !            exc,excx,excc,vxc2(1,1,i),vx2(1,1,i),vc2(1,1,i))
  !    enddo
  !    !       Restore rp,rps; also add -drho to rp and -drho/2 to rps
  !    call dscal(n1,1/(1+fac),rp,1)
  !    if (nsp == 2) then
  !       do  i = 1, nsp
  !          call daxpy(n1,-fac,rp,1,rps(1,1,i),1)
  !       enddo
  !    endif
  !    call dscal(n1,(1-fac),rp,1)
  !    !       Exchange potential at rp-drho
  !    do  i = 1, nsp
  !       call evxcv(rp,rps(1,1,i),n1,1,lxcf, &
  !            exc,excx,excc,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
  !    enddo
  !    !       Restore rp,rps
  !    call dscal(n1,1/(1-fac),rp,1)
  !    if (nsp == 2) then
  !       do  i = 1, nsp
  !          call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
  !       enddo
  !    endif

  !    do  i = 1, nsp
  !       do  ip = 1, np
  !          do  ir = 1, nr
  !             rhot = rp(ir,ip)
  !             if (rhot > 0) then
  !                f1 = vxc (ir,ip,i)*(rhot*(1-fac))**alfa
  !                f2 = vxc2(ir,ip,i)*(rhot*(1+fac))**alfa
  !                dfdr = (f2-f1)/(2d0*fac*rhot)
  !                vxc2(ir,ip,i) = dfdr
  !             else
  !                vxc2(ir,ip,i) = 0
  !             endif
  !          enddo
  !       enddo
  !    enddo
  ! else
     call dpzero(excx,n1)
     call dpzero(excc,n1)
  !endif
  ! --- vxc, exc for unperturbed density ---
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
  call vxcns4(0,ri,nr,rwgt,np,wp,nsp,rps,exc,excx,excc,vxc, &
       rep,repx,repc,rmu,qs)
  ! --- Add perturbation to vxc ---
  !     Integrals focexc = int rc vxc, focvxc= rc * dvxc/dr * rhot
  ! if (lpert) then
  !    focvxc = 0
  !    focexc = 0
  !    focex  = 0
  !    focec  = 0
  !    do  i  = 1, nsp
  !       do  ip = 1, np
  !          do  ir = 1, nr
  !             rhot = rp(ir,ip)
  !             if (rhot <= 0 .OR. rps(ir,ip,i) <= 0) then
  !                vxc2(ir,ip,i) = 0
  !             endif
  !             !             Debugging
  !             !              if (rps(ir,ip,i) .lt. 0 .or. rhot .lt. 0) then
  !             !                if (vxc(ir,ip,i) .ne. 0 .or. vxc2(ir,ip,i) .ne. 0) then
  !             !                  print *, vxc(ir,ip,i), vxc2(ir,ip,i)
  !             !                  stop 'oops'
  !             !                endif
  !             !              endif
  !             if (rps(ir,ip,i) > 0 .AND. rhot > 0) then
  !                f  = vxc(ir,ip,i)*rhot**alfa
  !                dfdr = vxc2(ir,ip,i)
  !                dvdr = (dfdr - alfa*f/rhot) / rhot**alfa
  !                weight = (ri(ir)**2*rwgt(ir))*wp(ip) * rc(ir)
  !                focvxc(i) = focvxc(i) + weight*dvdr*rps(ir,ip,i)
  !                focexc(i) = focexc(i) + weight*vxc(ir,ip,i)/nsp
  !                focex(i)  = focex(i)  + weight*vx(ir,ip,i)/nsp
  !                focec(i)  = focec(i)  + weight*vc(ir,ip,i)/nsp
  !                vxc(ir,ip,i) = vxc(ir,ip,i) + dvdr*rc(ir)
  !             endif
  !          enddo
  !       enddo
  !    enddo
  ! endif
  ! --- Scale yl by wp for fast multiplication --- ---
  do   ilm = 1, nlm
     do   ip = 1, np
        ylwp(ip,ilm) = yl(ip,ilm)*wp(ip) ! it was yl=yl*wp
     enddo
  enddo
  ! --- Add Yl-projection of vxc into vl ---
  do  30  i = 1, nsp
     call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,ylwp,np,1d0, vl(1,1,i),nr)
30 enddo
  ! --- Optionally make l-resolved vxc_L,exc_L ---
  if (mod(isw,10) == 1) then
     do  i = 1, nsp
        call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,ylwp,np,0d0, fl(1,1,i),nr)
     enddo
     call dgemm('N','N',nr,nlm,np,1d0,exc,nr,ylwp,np,0d0, & ! & yl-->ylwp
     fl(1,1,1+nsp),nr)
  endif
  ! --- Print out int (rl*vl) resolved by l ---
  if (ipr > 30) call vxcns5(0,ipr,'rho*vtot',nlm,nsp,nr,ri,rwgt,rl,vl,suml,ssum)
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
  implicit none
  integer :: nr,nsp,np,isw
  double precision :: ri(nr),wp(np),rep(nsp),repx(nsp),repc(nsp),rmu(nsp), &
       rwgt(nr),rp(nr,np,nsp),exc(nr,np),excx(nr,np),excc(nr,np), vxc(nr,np,nsp)
  integer :: ip,ipr,ir,i,isw1
  double precision :: weight,rpneg,qs(nsp)
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
           rep(i) = rep(i) + exc(ir,ip)*rp(ir,ip,i)*weight
           repx(i) = repx(i) + excx(ir,ip)*rp(ir,ip,i)*weight
           repc(i) = repc(i) + excc(ir,ip)*rp(ir,ip,i)*weight
           rmu(i) = rmu(i) + vxc(ir,ip,i)*rp(ir,ip,i)*weight
        enddo
     enddo
     if(ipr>=30.AND.i == 1)write(stdo,725)rmu(i),rep(i),qs(i)
     if(ipr>=30.AND.i == 2)write(stdo,726)rmu(i),rep(i),qs(i),sum(rmu),sum(rep),sum(qs)
725  format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
726  format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/ &
          '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)
24 enddo
  if (rpneg < 0 .AND. ipr >= 10) write(stdo,727) rpneg
727 format(' vxcnsp (warning): negative rho: min val = ',1pe10.2)
end subroutine vxcns4
subroutine vxcns5(isw,ipr,strn,nlml,nsp,nr,ri,rwgt,rl,vl,suml,ssum)
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
  implicit none
  character strn*(*)
  integer :: isw,ipr,nlml,nsp,nr
  double precision :: ri(nr),rwgt(nr),rl(nr,nlml,nsp),vl(nr,nlml,nsp), suml(0:20,nsp),ssum(nsp)
  integer :: ir,lmax,ilm,isp,l,ll
  double precision :: dot3
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
     ssum(isp) = sum(suml(0:lmax,isp))
     if (ipr > 30) then
        if (isp == 1) &
             write(stdo,341) strn,ssum(isp),(suml(l,isp),l = 0,lmax)
        if (isp == 2) &
             write(stdo,342) ssum(isp),(suml(l,isp),l = 0,lmax)
341     format(1x,a8,f13.6,' ... by l: ',f12.6,4f10.6:/(18x,4f10.6))
342     format('  spin 2:',f13.6,' ... by l: ',f12.6, &
             4f10.6:/(18x,4f10.6))
     endif
  enddo
end subroutine vxcns5
subroutine vxcnls(a,ri,nr,np,nlm,nsp, yl,gyl,rwgt,wp,rl,lxcg, vl,rep,rmu) 
  use m_xcpbe,  only: xcpbe
  !!= Gradient correction to nspher. density on a radial and angular mesh =
  !!*NOTE: seeing tol=1d-10 is changed on Dec1st 2010. But this is empirically determined for Co atom.
  !! If GGA is unstable than LDA, Check this routine. I am not so definite for tol and so.
  !!----------------------------------------------------------------------
  !!*Inputs
  !!   ri    :mesh of points
  !!   lcut=0 !!!  :1 if cutoff exc for small rho to avoid blowup in exc
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
  integer :: lcut=0,nr,np,nlm,nsp,lxcg
  double precision :: ri(nr),gyl(np,nlm,3),ylwp(np,nlm), &
       wp(np),rep(nsp),rmu(nsp),rwgt(nr),vxcnl(nr,nlm,nsp),excnl(nr,nlm), &
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
           if ( i == 2) print 7266, rmu(i),rep(i), sum(rmu),sum(rep)
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
  if(debug) print *,'vxcnls: goto vxcnls 77777a',lxcg,lmax,ir0,ir1,nr,np,nlm,nn
  call xxcnls2(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,ylwp &
       ,wp,rp,ggrp,grp,agrl,rl,rwgt, & !lcut, &
       vxcnl,excnl,sumnl)
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
     if (ipr >-1 .AND. i == 2) print 726, rmunl(i),repnl(i), sum(rmunl),sum(repnl)
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
     if ( i == 2) print 7261, rmu(i),rep(i),sum(rmu),sum(rep)
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
     ,wp,rp,ggrp,grp,agrl,rl,rwgt,& !lcut, &! & fixed bug (remove wkl) 11July2010.takao
     vxcnl,excnl,sumnl)
  !     implicit none
  use m_xclda,only: vxnloc
!  use m_vxcfunc,only: vxcgga
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
!     if (lxcg > 2) then
!        call vxcgga(lxcg,nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp), &
!             agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1), &
!             ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4), &
!             gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1), &
!             gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp), &
!             excp(ir0,ip))
!     else
!        if (lcut == 0) then
           call vxnloc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1), &
                agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp), &
                agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1), &
                gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1), &
                vxcp(ir0,ip,nsp),excp(ir0,ip))
!        else
!           call vxnlcc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1), &
!                agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp), &
!                agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1), &
!                gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1), &
!                vxcp(ir0,ip,nsp),excp(ir0,ip))
!        endif
!     endif
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


subroutine vxc0gc(nr,nsp,rofi,rwgt,rho,vxc,exc,rep,rmu,lxcfun)
  use m_xclda,only: vxcgr2
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
  double precision :: pi,rho2,rho3,ub4pi,rp(nrmx*nsp),rho0(nsp)
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

subroutine vxc0sp(a,b,rofi,rho,nr,v,rho0,rep,rmu,nsp,exrmx)
  use m_lmfinit,only: lxcf_g=>lxcf, stdo
  use m_ropyln,only: ropyln
  use m_xclda,only: evxcv,evxcp,vxcgr2
  !!- Adds xc part to spherical potential, makes integrals rmu and rep
  !! ----------------------------------------------------------------------
  ! i Inputs
  ! i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! i   rofi  :radial mesh points
  ! i   rho   :density = (true rho)*(4*pi*r**2)
  ! i   nr    :number of radial mesh points
  ! i   nsp   :2 for spin-polarized case, otherwise 1
  ! i  globalvariables%lxcf  :type of xc potential
  !!    1 : VWN
  !!    2 : Barth-Hedin
  !!   103: GGA PBE
  ! o Outputs
  ! o   v     :vxc is added to v
  ! o   rho0  :density extrapolated to origin
  ! o   rep   :integral rho * exc.
  ! o   rmu   :integral rho * vxc.
  ! o   exrmx :exchange energy density at rmax
  !!   vxc0sp contained in lm7.0 beta compiled by M.van Schilfgaarde.
  !! ----------------------------------------------------------------------
  implicit none
  integer :: nr,nsp
  double precision :: a,b,rofi(nr),v(nr,nsp),rho(nr,nsp), &
       rep(nsp),rmu(nsp),rho0(nsp),qs(nsp),exrmx
  double precision :: pi,rho2,rho3, &
       ub4pi,wgt,exc(nr),vxc(nr,nsp),rp(nr,nsp),repl(nsp),rmul(nsp) !exc(nrmx)-->exc(nr) automatic array.
  integer:: lx ,  i , ir , isp , iprint , nglob, lxcfun,lxcf
  real(8) ,allocatable :: grh_rv(:)
  real(8) ,allocatable :: agrh_rv(:)
  real(8) ,allocatable :: ggrh_rv(:)
  real(8) ,allocatable :: grgag_rv(:)
  real(8) ,allocatable :: excx_rv(:)
  real(8) ,allocatable :: excc_rv(:)
  real(8) ,allocatable :: vxcx_rv(:)
  real(8) ,allocatable :: vxcc_rv(:),vx2(:),vc2(:)
  character (2) :: st

  integer:: nth,nph,lxcg,np,lmax,nlm
  integer,parameter:: nnn=122
  real(8):: p(3,nnn),wp(nnn),rwgt(nr) ,p2(nnn,3),r2(nnn)
  real(8),allocatable:: gyl(:,:),yl(:),rll(:,:),vxc_(:,:)
  logical:: newv=.true.
  pi = 4d0*datan(1d0)
  ub4pi = 1d0/(4d0*pi)
  lxcfun = lxcf_g
  lxcf = mod(lxcfun,100)
  ! --- Add background rho to calculate vxc ---
  !      rhobg = 0d0
  !      call getsyv('rhobg',rhobg,i)
  !      call addzbk(rofi,nr,1,nsp,rho,rhobg,1d0)

  !! === Extrapolate rho to origin ===
  do  10  isp = 1, nsp
     rep(isp) = 0d0
     rmu(isp) = 0d0
     rho2 = rho(2,isp)/rofi(2)**2
     rho3 = rho(3,isp)/rofi(3)**2
     rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
10 enddo
  !! === Make true rho ===
  do   isp = 1, nsp
     rp(1,isp) = rho0(isp)
     do   ir = 2, nr
        rp(ir,isp) = rho(ir,isp)*ub4pi/rofi(ir)**2
     enddo
  enddo

  !      print *,'end of do 20'
  !! === Generate vxc,exc on a mesh ===
  !      print *,'end of do 20 xxxx lxcf2=',lxcf2
  if (lxcfun==103) then
     if( .NOT. newv) then
        allocate(excx_rv(nr))
        allocate(excc_rv(nr))
        allocate(vxcx_rv(nr))
        allocate(vxcc_rv(nr))
        allocate(vx2(nr),vc2(nr))
        call evxcp(rp,rp(1,2),nr,nsp,lxcf,excx_rv,excc_rv,exc, &
             vxcx_rv,vx2, vxcc_rv,vc2,vxc,vxc(1,2))
        deallocate(vx2,vc2)
        do isp = 1, nsp
           vxc(1,isp) = (vxc(2,isp)*rofi(3)-vxc(3,isp)*rofi(2))/(rofi(3)-rofi(2))
        enddo
        deallocate(vxcc_rv,vxcx_rv,excc_rv,excx_rv)
     else
        nph=0
        call fpiint(-4,nph,np,p,wp)
        p2=transpose(p)
        call radwgt(rofi(nr),a,nr,rwgt)
        lmax=0
        nlm=(lmax+1)**2 !yl(l=1) is needed for gradfl called in vxcnls
        allocate (yl(np*(lmax+2)**2),gyl(np*nlm,3)) !nlm=1 Only s channel. Sphelical atom.
        call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax+1,np,yl,r2)
        !         call ropylg(1,lmax,nlm,np,np,p2,p2(1+np),p2(1+2*np),r2,yl,gyl)
        gyl=0d0
        !         print *,'goto vxcnls only s channel',sum(abs(gyl))
        !         nlm=1
        !         allocate(rll(nr,nsp))
        !         rll(:,1:nsp)=
        !         print *,'goto vxcnls',sum(rll)
        call vxcnls(a,rofi,nr,np,nlm,nsp, &
             yl,gyl,rwgt,wp,sqrt(4d0*pi)*rp(:,1:nsp),lxcfun-100,  vxc,rep,rmu)
        deallocate(yl,gyl)
        vxc= vxc/sqrt(4d0*pi)
     endif
  elseif(lxcfun==1 .OR. lxcfun==2) then
     allocate(excx_rv(nr))
     allocate(excc_rv(nr))
     allocate(vxcx_rv(nr))
     allocate(vxcc_rv(nr))
     if (nsp == 1) then
        call evxcv ( rp , rp , nr , nsp , lxcf, &
             exc , excx_rv , excc_rv , vxc , vxcx_rv , vxcc_rv )
     else
        rp(1:nr,2)= rp(1:nr,1)+rp(1:nr,2)
        call evxcv ( rp ( 1 , 2 ) , rp , nr , 2 , lxcf, exc , excx_rv &
             , excc_rv , vxc , vxcx_rv , vxcc_rv )
        rp(1:nr,2) = rp(1:nr,2) - rp(1:nr,1)
        rp(1:nr,1) = rp(1:nr,1) + rp(1:nr,2)
        call evxcv ( rp , rp ( 1 , 2 ) , nr , 2 , lxcf, exc , excx_rv &
             , excc_rv , vxc ( 1 , 2 ) , vxcx_rv , vxcc_rv )
        rp(1:nr,1) = rp(1:nr,1) - rp(1:nr,2)
     endif
     deallocate(vxcc_rv,vxcx_rv,excc_rv,excx_rv)
     !! === Integrals ===
     do  14  i  = 1, nsp
        rep(i) = 0d0
        rmu(i) = 0d0
        do  12  ir = 1, nr
           wgt = 2*(mod(ir+1,2)+1)/3d0
           if (ir == 1 .OR. ir == nr) wgt = 1d0/3d0
           wgt = wgt * a*(rofi(ir)+b)
           !          qs(i)  = qs(i)  + wgt*rho(ir,i)
           rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
           rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
12      enddo
        !        repl(i) = rep(i)
        !        rmul(i) = rmu(i)
14   enddo
  endif
  if (lxcfun==103 .AND. ( .NOT. newv)) then
     call vxcgr2 ( nr , nsp , nr , rofi , rp , exc , vxc )
     ! ...   Redo integrals, with gradient correction
     do  24  i  = 1, nsp
        repl(i) = rep(i)
        rmul(i) = rmu(i)
        rep(i) = 0d0
        rmu(i) = 0d0
        do  22  ir = 1, nr
           wgt = 2*(mod(ir+1,2)+1)/3d0
           if (ir == 1 .OR. ir == nr) wgt = 1d0/3d0
           wgt = wgt * a*(rofi(ir)+b)
           rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
           rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
22      enddo
24   enddo
  endif
  !! --- Add to V ---
  v(:,1:nsp)=v(:,1:nsp)+vxc(:,1:nsp)
  exrmx = exc(nr)
  !! --- Undo background rho for purposes of calculating vxc !call addzbk(rofi,nr,1,nsp,rho,rhobg,-1d0)
  if (iprint() < 80) return
  write(stdo,"(/' vxc0sp: reps        rmu   ')")
  do i = 1, nsp
     st = ' '
     if (i < nsp) st = 'up'
     if (i == 2)   st = 'dn'
     write(stdo,335) st, rep(i),  rmu(i)
  enddo
  if(nsp==2) write(stdo,335) '  ', sum(rep), sum(rmu)
335 format(1x,a2,2x,4f12.6)
end subroutine vxc0sp
end module m_vxcatom
