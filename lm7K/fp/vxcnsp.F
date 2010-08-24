      subroutine vxcnsp(isw,a,ri,nr,rwgt,nlm,nsp,rl,lxcfun,rc,
     .focexc,focex,focec,focvxc,reps,repsx,repsc,rmu,vl,fl,qs)
C- Add vxc to potential in sphere and make integrals rho*vxc,rho*exc
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1s digit
Ci         : 1, make fl, otherwise not
Ci         :10s digit
Ci         : 0 calculated LDA for density as is.
Ci         : 1 reset points of negative density to positive density
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rl    :full charge density * r**2
Ci   lxcfun:specifies exchange-correlation functional
Ci         :1s digit sets local xc functional
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci         :100s digit sets gradient corrections
Ci         :  0    LSDA
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Ci         :10000s digit sets flag to treat core in perturbation th.
Ci   rc    :l=0 charge density*4*pi*ri**2 to be added in pert. theory
Ci         :Not used unless 100s digit of lxcfun is set.
Co Outputs
Co   focex :(100s digit lxcfun): integral rc * vx
Co   focec :(100s digit lxcfun): integral rc * vc
Co   focexc:(100s digit lxcfun): integral rc * vxc
Co   focvxc:(100s digit lxcfun): integral rc * (dvxc/drho * rho)
Co         : for these integrals, see Remarks
Co   reps  :integral rl*exc
Co   repsx :integral rl*ex
Co   repsc :integral rl*ec
Co   rmu   :integral rl*vxc
Co   vl    :local exchange-correlation potential added to vl
Co   fl    :fl(:,:,1:nsp) = vxc by L
Co         :fl(:,:,1+nsp) = exc by L
Co         :fl is made only if isw is nonzero
Co   qs    :spin-resolved charge
Cl Local variables
Cl   yl    :YL(rp)
Cl   lxcf  :LD part of LDA+GGA functional
Cr Remarks
Cr   Perturbation treatment:
Cr     vxc[rho + rc] = vxc[rho] + rc * (dvxc/drho)
Cr   Correction to int rho * vxc:
Cr     int rho * vxc = int rho*vxc[rho] + int rho*rc*(dvxc/drho)
Cr     Correction = focvxc = int rho * rc * (dvxc/drho)
Cr   Correction to int rho * exc:
Cr     Corr = int rho * rc * (dexc/drho) = int rho * rc * (vxc-exc)/rho
Cr          = int rc * (vxc-exc) = focexc - int rc exc
Cr     Second term is not retained.
Cu Updates
Cu   06 Apr 09 Routine calls vxcnls to handle some GGA potentials
Cu   14 Mar 04 Optionally Makes fl.  New argument list
Cu   14 Jun 02 repsx and repsc (T. Miyake)
Cu    8 Feb 02 focex and focec (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    3 Jun 00 Adapted from old FP vxcnsp; pert. corr from nfp vxcdnsp.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isw,nr,nlm,lxcfun,nsp
      double precision ri(nr),reps(2),rmu(2),repsx(2),repsc(2),
     .rwgt(nr),rc(nr),rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,1+nsp)
      double precision focexc(2),focex(2),focec(2),focvxc(2)
      double precision qs(2)
C ... Local parameters
      integer nnn,nlmx
      parameter(nnn=122,nlmx=64)
      logical lpert
      integer stdo,lgunit,ipr,ll,lmax,np,nph,nth,nxl(0:7),lxcf,lxcg
      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn),fpi
C     double precision yl(nnn*nlmx)
      real(8), allocatable :: ylwp(:)
      real(8), allocatable :: yl(:),rp(:,:,:),gyl(:,:),grp(:,:),agrp(:),
     .ggrp(:),agrl(:),vxcnl(:),excnl(:)

      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
      logical:: debug=.false.
      real(8):: a !takao used to get derivetive of ri. See vxcnls.F
      real(8),allocatable::vlxc(:,:,:)

      stdo = lgunit(1)
      call getpr(ipr)
      lxcf = mod(lxcfun,100)
      lxcg = mod(lxcfun/100,100)
      lpert = lxcfun .ge. 10000
      fpi = 16d0*datan(1d0)

C ... Create angular integration mesh
      lmax = ll(nlm)
C     if (lxcg .ne. 0) lmax = lmax+1
      if (lmax .gt. 6) then
        nth = 2*lmax+2
        nph = 0
      else
        nth = nxl(lmax)
        nph = 0
      endif
      call fpiint(nth,nph,np,p,wp)
      if(ipr .ge. 30) write(stdo,200) nth,nph,np,nr
  200 format(' mesh:   nth,nph=',2i4,'   gives',i4,'  angular points,',
     .'   nrad=',i4)
      if (np .gt. nnn) call rxi('vxcnsp: increase nnn, need',np)
      if ((lmax+2)**2*np .gt. nlmx*nnn) call rx('vxcnsp: increase nlm')

      rmu(2)  = 0
      reps(2) = 0
      repsx(2) = 0
      repsc(2) = 0

C ... Scale rl to true density
      call vxcns3(nr,nlm,nsp,ri,rl,1)
      if (lpert) then
        call vxcns3(nr,1,1,ri,rc,1)
        call dscal(nr,1/fpi,rc,1)
      endif

C ... Spherical harmonic polynomials ... to be overwritten by Yl * wp
      call dmcpy(p,1,3,p2,np,1,np,3)
C     GGA case:  Need both yl*wp and YL to lmax+1 to make grad YL
      if (lxcg .ne. 0) then
        allocate (yl(np*(lmax+2)**2),ylwp(np*(lmax+2)**2))
        call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax+1,np,yl,r2)
        call dcopy(np*nlm,yl,1,ylwp,1)
      else
        allocate (ylwp(np*nlm))
        call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,ylwp,r2)
      endif

      allocate(vlxc(nr,nlm,nsp))
C ... Generate the local exchange-correlation potential
C     On output, yl is overwritten by YL * wp  <-------Confusing!!! remove this in future.
      vlxc=0d0
      call vxcns2(isw,ri,nr,rwgt,np,wp,ylwp,nlm,nsp,rl,rc,lxcf,
     .lpert,focexc,focex,focec,focvxc,reps,repsx,repsc,rmu,
     .vlxc,fl,qs)

C --- Gradient correction to LD potential ---
      if(lxcg/=0) then
        allocate (gyl(np*nlm,3))
        call ropylg(1,lmax,nlm,np,np,p2,p2(1+np),p2(1+2*np),r2,yl,gyl)
        allocate(grp(nr*np*nsp,3),agrp(nr*np*nsp),ggrp(nr*np*nsp))
        allocate(rp(nr,np,nsp))
        allocate(agrl(nr*nlm*nsp),vxcnl(nr*nlm*nsp),excnl(nr*nlm*nsp))
        if(debug) print *,'goto vxcnls'
        call vxcnls(a,ri,0,nr,np,nlm,nsp,rp,grp,ggrp,agrp,yl,gyl,ylwp,
     .  rwgt,wp,rl,agrl,lxcg,vxcnl,excnl,vlxc,reps,rmu)
        deallocate(yl,gyl,grp,agrp,ggrp,agrl,vxcnl,excnl,rp)
      endif
      vl=vl+vlxc
      deallocate(vlxc)

C ... Undo the rl scaling
      call vxcns3(nr,nlm,nsp,ri,rl,0)
      if (lpert) then
        call vxcns3(nr,1,1,ri,rc,0)
        call dscal(nr,fpi,rc,1)
      endif

      deallocate(ylwp)

      end

      subroutine vxcns2(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,rc,lxcf,
     .lpert,focexc,focex,focec,focvxc,rep,repx,repc,rmu,vl,fl,qs)
C- Make vxc, rep and rmu in sphere, for local XC functional.
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1s digit
Ci         : 1, make fl, otherwise not
Ci         :10s digit
Ci         : 0 calculated LDA for density as is.
Ci         : 1 reset points of negative density to positive density
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   np    :number of spherical mesh points
Ci   wp    :spherical mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rl    :charge density
Ci   rc    :l=0 charge density*4*pi*ri**2 to be added in pert. theory
Ci   lxcf  :specifies exchange-correlation functional
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci   lpert :Make perturbation integrals focexc and focvxc
Cio Inputs/Outputs
Cio  yl    :Spherical harmonics YL(rp)
Cio        :yl is OVERWRITTEN by YL * wp on output
Co Outputs
Co   focexc:(lpert only): integral rc * vxc
Co   focex :(lpert only): integral rc * vx
Co   focec :(lpert only): integral rc * vc
Co   focvxc:(lpert only): integral rc * (dvxc/drho * rho)
Co   rep   :integral rl*exc
Co   repx  :integral rl*ex
Co   repc  :integral rl*ec
Co   rmu   :integral rl*vxc
Co   vl    :local exchange-correlation potential added to vl
Co   fl    :fl(:,:,1:nsp) = vxc by L
Co         :fl(:,:,1+nsp) = exc by L
Co         :fl is only made if 1s digit isw is nonzero
Co   qs    :spin-resolved charge
Cl Local variables
Cl   rp    :list of points on the spherical mesh
Cl   exc   :exchange density on (nr,np) mesh
Cl   vxc   :exchange potential on (nr,np) mesh
Cr Remarks
Cr   For perturbation treatment, take numerical derivatives
Cr   df/dr = d/dr (vxc*r**alfa) instead of d/dr (vxc) because
Cr   f is nearly linear for alpha=2/3.
Cr
Cr   In the spin polarized case, the perturbation density rc is not
Cr   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
Cr   we use dm=0 and thus drho1 = drho2 = drho/2; thus
Cr     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
Cu Updates
Cu   20 Nov 09 New 20s digit for isw
Cu   02 Apr 09 New option (10's digit isw)
Cu   14 Jun 02 repx and repc (T. Miyake)
Cu    8 Feb 02 focex and focec (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    3 Jun 00 adapted from old FP code
Cu    10.07.96 dln: modified for Perdew GGA
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lpert
      integer nr,nlm,nsp,lxcf,np,isw
      double precision ri(1),yl(np,nlm),wp(np),rep(2),repx(2),repc(2),
     .rmu(2),rwgt(nr),focexc(2),focex(2),focec(2),focvxc(2),
     .rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,nsp+1),rc(nr)
C ... Local parameters
      integer ilm,ip,ipr,ir,i,n1,ineg(2),isw1
      double precision weight,suml(0:40),fac,dmach,alfa,sum(2)
      double precision rhot,f1,f2,dfdr,dvdr,f
      double precision qs(2),tol
      real(8), allocatable :: rp(:,:),rps(:,:,:),exc(:,:),
     .excx(:,:),excc(:,:)
      real(8), allocatable :: vxc(:,:,:),vx(:,:,:),
     .vc(:,:,:),vxc2(:,:,:),vx2(:,:,:),vc2(:,:,:)
      real(8),allocatable:: vlxc(:,:,:)


      allocate (rp(nr,np),rps(nr,np,nsp),exc(nr,np),
     .excx(nr,np),excc(nr,np))
      allocate (vxc(nr,np,nsp),vx(nr,np,nsp),vc(nr,np,nsp),
     .vxc2(nr,np,nsp),vx2(nr,np,nsp),vc2(nr,np,nsp))

C     stdo = lgunit(1)
      call getpr(ipr)

C     'small' density
      tol = 1d-15
      isw1 = mod(isw/10,10)

C     for numerical differentiation
      fac = dmach(1)**(1d0/3d0)
      alfa = 2d0/3d0
      n1 = nr*np
C     call prrmsh('rl',rl,rl,nr,nr,nlm*2)

C --- Generate density point-wise through sphere ---
      call dpzero(rp,n1)
      ineg(1) = 0
      ineg(2) = 0
      do  i = 1, nsp
        call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0,
     .  rps(1,1,i),nr)

cccccccccccccccccccccccccccc
          do  ip = 1, np
            do  ir = 1, nr
              if (rps(ir,ip,i) .lt. 0d0) then
c                print *,'xxx vvv', ir,ip,rps(ir,ip,i)
                ineg(i) = ineg(i) + 1
              endif
            enddo
          enddo
cccccccccccccccccccccccccccc

C   ... Counts number of points with negative density; reset to pos.
        if (isw1 .eq. 1) then
          do  ip = 1, np
            do  ir = 1, nr
              if (rps(ir,ip,i) .lt. 0d0) then
                ineg(i) = ineg(i) + 1
                rps(ir,ip,i) = tol
              endif
            enddo
          enddo
        endif

        call daxpy(n1,1d0,rps(1,1,i),1,rp,1)
      enddo


      if (ineg(1)+ineg(2) .ne. 0) then
        write(6,"(a,i10,a,2i10)") ' vxcns2 (warning): nr*np=',nr*np,
     &  "  negative density # of points=", ineg(1:nsp)
c        call info5(5,0,0,' vxcns2 (warning):  negative density,'
c     .  //' %i points %?#n==2#(spin 1) %i points (spin 2)##',
c     .  ineg(1),nsp,ineg(2),0,0)
      endif

C --- Perturbation treatment: df/dr in vxc2 ---
      if (lpert) then

C --- Stop if GGA ---
        if (lxcf .gt. 2)
     .  call rx('vxcnsp: Perturbation treatment'//
     .  ' is not implemented for GGA')

C       Add rp*fac/2 into rp+, rp- and fac*rp into rp
        if (nsp .eq. 2) then
          do  i = 1, nsp
            call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
          enddo
        endif
        call dscal(n1,1+fac,rp,1)
C       Exchange potential at rp+drho
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,1,lxcf,
     .    exc,excx,excc,vxc2(1,1,i),vx2(1,1,i),vc2(1,1,i))
        enddo
C       Restore rp,rps; also add -drho to rp and -drho/2 to rps
        call dscal(n1,1/(1+fac),rp,1)
        if (nsp .eq. 2) then
          do  i = 1, nsp
            call daxpy(n1,-fac,rp,1,rps(1,1,i),1)
          enddo
        endif
        call dscal(n1,(1-fac),rp,1)
C       Exchange potential at rp-drho
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,1,lxcf,
     .    exc,excx,excc,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
        enddo
C       Restore rp,rps
        call dscal(n1,1/(1-fac),rp,1)
        if (nsp .eq. 2) then
          do  i = 1, nsp
            call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
          enddo
        endif

        do  i = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              rhot = rp(ir,ip)
              if (rhot .gt. 0) then
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

C --- vxc, exc for unperturbed density ---
C     print *, '!!'; lxcf=1
      if (lxcf .gt. 2) then
        call evxcp(rps,rps(1,1,nsp),n1,nsp,lxcf,excx,excc,exc,
     .  vx(1,1,1),vx(1,1,2),vc(1,1,1),vc(1,1,2),
     .  vxc,vxc(1,1,nsp))
      else
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,nsp,lxcf,
     .    exc,excx,excc,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
        enddo
      endif
      if (isw1 .eq. 2) then
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              if (rp(ir,ip) .le. 0 .or. rps(ir,ip,i) .le. 0) then
                vxc(ir,ip,i) = 0
                vx(ir,ip,i) = 0
                vc(ir,ip,i) = 0
              endif
            enddo
          enddo
        enddo
      endif

cccccccccccccccccccccccc
c      print *,' np=',np
c      do ir=1,nr
c        write(6,"(i5,5d13.5)") ir,exc(ir,np-4:np)
c      enddo
c      stop 'xxxxxxxxxxxx test exc vxcnsp'
cccccccccccccccccccccccc

C --- Integrals rep, rmu ---
C      rpneg = 0
C      do  24  i = 1, nsp
C        qs(i) = 0d0
C        rep(i) = 0d0
C        repx(i) = 0d0
C        repc(i) = 0d0
C        rmu(i) = 0d0
C        do  22  ip = 1, np
C        do  22  ir = 1, nr
C        rpneg = min(rpneg,rps(ir,ip,i))
C        weight = (ri(ir)**2*rwgt(ir))*wp(ip)
C        qs(i)  = qs(i)  + rps(ir,ip,i)*weight
C        rep(i) = rep(i) + exc(ir,ip)*rps(ir,ip,i)*weight
C        repx(i) = repx(i) + excx(ir,ip)*rps(ir,ip,i)*weight
C        repc(i) = repc(i) + excc(ir,ip)*rps(ir,ip,i)*weight
C   22   rmu(i) = rmu(i) + vxc(ir,ip,i)*rps(ir,ip,i)*weight
C        if (ipr.ge.30 .and. i.eq.1) write(stdo,725) rmu(i),rep(i),qs(i)
C        if (ipr.ge.30 .and. i.eq.2) write(stdo,726) rmu(i),rep(i),qs(i),
C     .    rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2)
C  725   format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
C  726   format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/
C     .         '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)
C   24 continue
C      if (rpneg .lt. 0 .and. ipr .ge. 20) write(stdo,727) rpneg
C  727 format(' vxcnsp (warning): negative rho: min val = ',1pe10.2)
      call vxcns4(0,ri,nr,rwgt,np,wp,nsp,rps,exc,excx,excc,vxc,
     .rep,repx,repc,rmu,qs)

C --- Add perturbation to vxc ---
C     Integrals focexc = int rc vxc, focvxc= rc * dvxc/dr * rhot
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
              if (rhot .le. 0 .or. rps(ir,ip,i) .le. 0) then
                vxc2(ir,ip,i) = 0
              endif
C             Debugging
C              if (rps(ir,ip,i) .lt. 0 .or. rhot .lt. 0) then
C                if (vxc(ir,ip,i) .ne. 0 .or. vxc2(ir,ip,i) .ne. 0) then
C                  print *, vxc(ir,ip,i), vxc2(ir,ip,i)
C                  stop 'oops'
C                endif
C              endif
              if (rps(ir,ip,i) .gt. 0 .and. rhot .gt. 0) then
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

C --- Scale yl by wp for fast multiplication --- ---
      do  12  ilm = 1, nlm
        do  12  ip = 1, np
   12 yl(ip,ilm) = yl(ip,ilm)*wp(ip)

C --- Add Yl-projection of vxc into vl ---
      do  30  i = 1, nsp
   30 call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,1d0,
     .vl(1,1,i),nr)


C$$$ccccccccccccccccccccccc
c      allocate(vlxc(nr,nlm,nsp))
c      do   i = 1, nsp
c        call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,0d0,
c     .  vlxc(1,1,i),nr)
c      enddo
c      isp=1
c      do ir=1,nr
c       write(6,"(i5,8d12.4)") ir, vxc(ir,1:8,isp) !vlxc(ir,1:8,isp)
c      enddo
C$$$      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxx test end vxcnsp 123'
C$$$cccccccccccccccccccccc

C --- Optionally make l-resolved vxc_L,exc_L ---
      if (mod(isw,10) .eq. 1) then
        do  i = 1, nsp
          call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,0d0,
     .    fl(1,1,i),nr)
        enddo
        call dgemm('N','N',nr,nlm,np,1d0,exc,nr,yl,np,0d0,
     .  fl(1,1,1+nsp),nr)
      endif

C --- Print out int (rl*vl) resolved by l ---
      if (ipr .gt. 30) then
        call vxcns5(0,ipr,'rho*vtot',nlm,nsp,nr,ri,rwgt,rl,vl,suml,sum)
C      lmax = ll(nlm)
C      do  42  i = 1, nsp
C      do  43  l = 0, lmax
C   43 suml(l) = 0d0
C      do  40  ilm = 1, nlm
C      l = ll(ilm)
C      do  40  ir = 1, nr
C   40 suml(l) = suml(l) + rl(ir,ilm,i)*vl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
C      if (i .eq. 1) write(stdo,341) (suml(l),l = 0,lmax)
C      if (i .eq. 2) write(stdo,342) (suml(l),l = 0,lmax)
C  341 format(' rho*vxc by l: ',f13.6,4f10.6:/(18x,4f10.6))
C  342 format('       spin 2: ',f13.6,4f10.6:/(18x,4f10.6))
C   42 continue
      endif
      deallocate (rp,rps,exc,excx,excc)
      deallocate (vxc,vx,vc,vxc2,vx2,vc2)
      end

      subroutine vxcns3(nr,nlm,nsp,ri,rl,isgn)
C- Scales rho by r**2, or undoes scaling
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ri    :radial mesh points
Ci   isgn  :1, scale rl by 1/r**2
Ci         :else scale rl by r**2
Co Outputs
Co   rl   :rl is scaled by r**2 or 1/r**2, depending on isgn
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer nr,nlm,nsp,isgn,i,ilm,ir
      double precision rl(nr,nlm,nsp),ri(nr),rho2,rho3

C --- Scale rho by 1/r**2 ---
      if (isgn .eq. 1) then
        do  10  i = 1, nsp
          do  10  ilm = 1, nlm
            rl(1,ilm,i) = 0d0
            do  10  ir = 2, nr
   10   rl(ir,ilm,i) = rl(ir,ilm,i)/ri(ir)**2
C  ...  Extrapolate rho to origin
        do  20  i = 1, nsp
          rho2 = rl(2,1,i)
          rho3 = rl(3,1,i)
   20   rl(1,1,i) = (rho2*ri(3)-rho3*ri(2))/(ri(3)-ri(2))
      else
        do  110  i = 1, nsp
          do  110  ilm = 1, nlm
            do  110  ir = 1, nr
  110   rl(ir,ilm,i) = rl(ir,ilm,i)*ri(ir)**2
      endif
      end

      subroutine vxcns4(isw,ri,nr,rwgt,np,wp,nsp,rp,exc,excx,excc,vxc,
     .rep,repx,repc,rmu,qs)
C- Integrals of vxc, reps and rmu in sphere.
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   np    :number of spherical mesh points
Ci   wp    :spherical mesh weights
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cl   rp    :list of points on the spherical mesh
Cl   exc   :exchange-correlation energy density on (nr,np) mesh
Cl   excx  :exchange             energy density on (nr,np) mesh
Cl   excc  :correlation          energy density on (nr,np) mesh
Cl   exc   :exchange-correlation energy density on (nr,np) mesh
Cl   vxc   :exchange correlation potential on (nr,np) mesh
Co Outputs
Co   rep   :integral rl*exc
Co   repx  :integral rl*ex
Co   repc  :integral rl*ec
Co   rmu   :integral rl*vxc
Co   qs    :spin-resolved charge
Cl Local variables
Cl   rpneg :number of points at which density < 0
Cr Remarks
Cu Updates
Cu   12 Mar 04 created from vxcns2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nsp,np,isw
      double precision ri(1),wp(np),rep(2),repx(2),repc(2),rmu(2),
     .rwgt(nr),rp(nr,np,nsp),exc(nr,np),excx(nr,np),excc(nr,np),
     .vxc(nr,np,nsp)
C ... Local parameters
      integer ip,ipr,ir,i,stdo,lgunit,isw1
      double precision weight,rpneg,qs(2)

      stdo = lgunit(1)
      call getpr(ipr)
      isw1 = mod(isw/10,10)

C --- Integrals reps, rmu ---
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
            if (isw1 .eq. 2) then
              if (rp(ir,ip,1)+rp(ir,ip,nsp).le.0 .or. rp(ir,ip,i).lt.0) then
                exc(ir,ip) = 0
                excx(ir,ip) = 0
                excc(ir,ip) = 0
                vxc(ir,ip,i) = 0
              endif
            endif
C         Debugging
C          if (rp(ir,ip,1) .lt. 0 .or. rp(ir,ip,nsp) .lt. 0) then
C            if (vxc(ir,ip,i) .ne. 0) then
C              print *, vxc(ir,ip,i)
C              stop 'oops'
C            endif
C          endif
            rep(i) = rep(i) + exc(ir,ip)*rp(ir,ip,i)*weight
            repx(i) = repx(i) + excx(ir,ip)*rp(ir,ip,i)*weight
            repc(i) = repc(i) + excc(ir,ip)*rp(ir,ip,i)*weight
            rmu(i) = rmu(i) + vxc(ir,ip,i)*rp(ir,ip,i)*weight
          enddo
        enddo
        if (ipr.ge.30 .and. i.eq.1) write(stdo,725) rmu(i),rep(i),qs(i)
        if (ipr.ge.30 .and. i.eq.2) write(stdo,726) rmu(i),rep(i),qs(i),
     .  rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2)
  725   format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
  726   format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/
     .  '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)
   24 continue
      if (rpneg .lt. 0 .and. ipr .ge. 10) write(stdo,727) rpneg
  727 format(' vxcnsp (warning): negative rho: min val = ',1pe10.2)

      end

      subroutine vxcns5(isw,ipr,strn,nlml,nsp,nr,ri,rwgt,rl,vl,suml,sum)
C- Integrals of rl*vl resolved by l
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :0 rl is true density
Ci         :1 rl is true density * ri**2
Ci   strn  :used in printout
Ci   ri    :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :number of Ylm's in rl and vl
Co Outputs
Co   suml  :l-and spin-resolved integrals of rl*vl
Co   sum   :spin-resolved integrals of rl*vl
Cr Remarks
Cu Updates
Cu   12 Mar 04 created from vxcns2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character strn*(*)
      integer isw,ipr,nlml,nsp,nr
      double precision ri(nr),rwgt(nr),rl(nr,nlml,nsp),vl(nr,nlml,nsp),
     .suml(0:20,2),sum(2)
C ... Local parameters
      integer ir,stdo,lgunit,lmax,ilm,isp,l,ll
      double precision dsum,dot3

      stdo = lgunit(1)
      lmax = ll(nlml)

      do   isp = 1, nsp
        call dpzero(suml(0,isp),lmax+1)
        do  ilm = 1, nlml
          l = ll(ilm)
          if (isw .eq. 1) then
            suml(l,isp) = suml(l,isp) +
     .      dot3(nr,rl(1,ilm,isp),vl(1,ilm,isp),rwgt)
          else
            do  ir = 1, nr
              suml(l,isp) = suml(l,isp) +
     .        rl(ir,ilm,isp)*vl(ir,ilm,isp)*ri(ir)**2*rwgt(ir)
            enddo
          endif
        enddo

        sum(isp) = dsum(lmax+1,suml(0,isp),1)

        if (ipr .gt. 30) then
          if (isp .eq. 1)
     .    write(stdo,341) strn,sum(isp),(suml(l,isp),l = 0,lmax)
          if (isp .eq. 2)
     .    write(stdo,342) sum(isp),(suml(l,isp),l = 0,lmax)
  341     format(1x,a8,f13.6,' ... by l: ',f12.6,4f10.6:/(18x,4f10.6))
  342     format('  spin 2:',f13.6,' ... by l: ',f12.6,
     .    4f10.6:/(18x,4f10.6))
        endif

      enddo

      end

C      subroutine fp2yl(nr,nlm,nsp,np,wp,fp,yl,fl)
CC- Add Yl-projection of function tabulated on a mesh
C      implicit none
C      integer nr,nlm,np,ip,ilm,ir,i,nsp
C      double precision fl(nr,nlm,nsp),fp(nr,np,nsp),yl(np,nlm),wp(np),xx
C
C      do  20  i = 1, nsp
C      do  20  ip = 1, np
C      do  20  ilm = 1, nlm
C      xx = wp(ip)*yl(ip,ilm)
C      do  20  ir = 1, nr
C   20 fl(ir,ilm,i) = fl(ir,ilm,i) + fp(ir,ip,i)*xx
C      end

