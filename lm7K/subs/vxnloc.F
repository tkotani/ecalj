      subroutine vxnloc(n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom,
     .  grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
C- Langreth-Mehl-Hu gradient correction to exc and vxc
C ---------------------------------------------------------------
Ci Inputs:
Ci   rhop,rhom, grhop,grhom, ggrhop,ggrhom (rho, grad and Laplacian rho)
Ci   grggr,grggrp,grggrm (grad rho . grad abs grad rho)
Ci   (nsp=2) grho: grad total rho, grpgrm: grad rho+ . grad rho-
Co Outputs:
Co   vxc, exc
Cr Remarks:
Cr   References PRB28,1809(1983) and PRB40,1997(1989).
Cr   If nsp=1, rhop, grhop, etc are for total rho.
Cr   factor f is empirically determined cutoff.  f=0 for pure gradient.
Cr   cutoff eliminates the blowup of vxc for small rho,
Cr   for numerical convenience.  Should be no significant change if 0.
Cr   Langreth-Mehl form does not use grggrp,grggrm
Cr   Should be used together with Barth-Hedin functional.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsp,n,i
      double precision rhop(1),rhom(1),grhop(1),grhom(1),grho(1),
     .  ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1),
     .  vxc1(1),vxc2(1),exc(1)
C Local parameters
      logical warned
      double precision pi,fivth,th4,th2,th,sth,sevni,f8,f9,rho,gp2,gm2,
     .                 bigf,aa,d,polar,f,cutoff,cutof1,gro,ggro,g2,
     .                 cutof2,hh,rp,rm,grp,grm,ggrp,ggrm,rmin,xd,expf,xx
      parameter (f=0.15d0,hh=1d-3,fivth=5d0/3d0,th4=4d0/3d0,
     .  th2=2d0/3d0,th=1d0/3d0,sth=7d0/3d0,sevni=7d0/9d0, rmin=1d-15)

      warned = .false.
      pi = 4d0*datan(1d0)
      f8 = (9d0*pi)**(1d0/6d0)*f
      f9 = 5d0/6d0*2d0**th2
      aa = (pi/(8d0*(3d0*pi*pi)**th4))
      if (nsp .eq. 1) then
        do  10  i = 1, n
        rho = rhop(i)
        if (rho .gt. rmin) then
        gro  = grhop(i)
        ggro = ggrhop(i)

        bigf = f8*gro/(rho**(7d0/6d0))
        cutoff = dexp(-hh*(gro*gro/(rho*rho))*rho**(-th2))
        expf = 2d0*dexp(-bigf)
        exc(i) = exc(i) + aa*(gro*gro/(rho**sth))*(expf-sevni)
        vxc1(i) = vxc1(i) + 2d0*aa*rho**(-th)*
     .    (sevni*(ggro/rho-th2*gro*gro/(rho*rho)) - expf*(
     .      ggro*(1d0-bigf/2d0)/rho -
     .      (th2+bigf*(-11d0/6d0+bigf*7d0/12d0))*gro*gro/(rho*rho) +
     .      bigf*(bigf-3d0)*grggr(i)/(2*gro*rho)))*cutoff
        endif
   10   continue
      else
        do  20  i = 1, n
        rho = rhop(i)+rhom(i)
        if (rho .gt. rmin) then
        rp   = rhop(i)
        rm   = rhom(i)
        if (rp .le. 0d0) then
          if (.not. warned)
     .      print *, 'vxnloc (warning) rho+ not positive but rho is'
          rp = rho/1000
          warned = .true.
        endif
        if (rm .le. 0d0) then
          if (.not. warned)
     .    print *, 'vxnloc (warning) rho- not positive but rho is'
          rm = rho/1000
          warned = .true.
        endif
        grp  = grhop(i)
        grm  = grhom(i)
        gro  = grho(i)
        g2   = gro*gro
        gp2  = grp*grp
        gm2  = grm*grm
        ggrp = ggrhop(i)
        ggrm = ggrhom(i)
        ggro = ggrp+ggrm
        polar = (rp-rm)/rho
        bigf = f8*gro/(rho**(7d0/6d0))
        d = dsqrt(((1d0+polar)**fivth+(1d0-polar)**fivth)/2d0)
        expf = 2d0/d*dexp(-bigf)

        exc(i) = exc(i) + aa/rho*(expf*g2/(rho**th4)
     .   - sevni/(2d0**th)*(gp2/(rp**th4)+gm2/(rm**th4)))
        cutof1 = dexp(-hh*(gp2/(rp*rp))*(2*rp)**(-th2))
        xd = f9*rho**(th-4d0)/(d*d)*(rp**th2-rm**th2)
        xx = (2d0-bigf)*ggro/rho -
     .    (th4+bigf*(-11d0/3d0+bigf*7d0/6d0))*g2/(rho*rho) +
     .    bigf*(bigf-3d0)*grggr(i)/(gro*rho)
        vxc1(i) = vxc1(i) + aa/rho**th*(
     .   -sevni/(2d0**th)*(rho/rp)**th*(th4*gp2/(rp*rp)-2d0*ggrp/rp)
     .    -expf*(xx -
     .    xd*((1d0-bigf)*rm*g2-(2d0-bigf)*rho*(grpgrm(i)+gm2))))
     .    *cutof1
        cutof2=dexp(-hh*(gm2/(rm*rm))*(2*rm)**(-th2))
        vxc2(i) = vxc2(i) + aa/rho**th*(
     .   -sevni/(2d0**th)*(rho/rm)**th*(th4*gm2/(rm*rm)-2d0*ggrm/rm)
     .    -expf*(xx +
     .    xd*((1d0-bigf)*rp*g2-(2d0-bigf)*rho*(gp2+grpgrm(i)))))
     .    *cutof2
        endif
   20   continue
      endif
      end
      subroutine pbevxc(nr,nsp,nrx,rofi,rp,exc,vxc,mode,units,
     . grh,ggrh,agrh,grgagr)
C- Interface for Perdew-Burke-Ernzerhof generalized gradient approximation
c  10.07.96  dln
C ---------------------------------------------------------------
Ci Inputs:
Ci   nr   - No. of points in the radial mesh
Ci   nsp  - 1 spin-restricted, 2- spin-pol.
Ci   nrx  - Max. dimension of arrays (> or = nr)
Ci   rofi - radial points
Ci   rp   - rho
Ci   mode - 0 : LSDA
Ci          1 : PW91
Ci          2 : PBE
Ci         -1 : PW91 correction only
Ci         -2 : PBE  correction only
Ci         10 : LSDA, same as 0, but no final extrapolation
Ci   units - Ry : output in Rydbergs
Ci           not Ry: Hartree
Co Outputs:
Co   exc(nrx) -  exchange-corr. energy
Co   vxc(nrx) -  exchange-corr. potential
Cr Remarks:
Cr   Gradient arrays are NOT used if no GGA is required.
Cr   This saves memory
Cr   working arrays arrays:
Cr   grh -
Cr   ggrh -
Cr   agrh -
Cr   grgagrh -
Cr   References:
C ----------------------------------------------------------------
C     implicit none

c - Input variables
      integer nr,nsp,nrx,mode,isp
      double precision rp(nrx,nsp),
     .  exc(nrx),vxc(nrx,nsp),rofi(nr),
     .  Xagrh1,Xagrh2,Xgrgagr1,Xgrgagr2,Xggrh1,Xggrh2,XXagrh,
     .  XXgrgagr

c - Work arrays and internal variables
      double precision grh(nrx,2),ggrh(nrx,2),
     .  agrh(nrx,4),grgagr(nrx,3)
      double precision exlsd, vxuplsd, vxdnlsd, eclsd,
     .  vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91,
     .  ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe,
     .  vxdnpbe, ecpbe, vcuppbe, vcdnpbe
      integer ir,i,iprint
      character*2, units

      if(mode.eq.0.or.mode.eq.10) GOTO 9999 !! Do not calculate gradients

      call pshpr(iprint()-20)

C --- grad(rho), laplacian rho ---
      call radgrx(nr,nrx,nsp,rofi,rp,grh)
      call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
      do  20  i  = 1, nsp
      do  24  ir = 2, nr
   24 ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
      ggrh(1,i) =(rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))

C --- grad rho . grad abs grad rho ---
      do  26  ir = 1, nr
   26 agrh(ir,i) = dabs(grh(ir,i))
      call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
      do  28  ir = 1, nr
   28 grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
   20 continue

C --- Extra terms g(n), g(n+).g(n-), g(n).g(abs(g(n))) if spin pol ---
      if (nsp .eq. 2) then
        do  32  ir = 1, nr
   32   agrh(ir,3) = dabs(grh(ir,1)+grh(ir,2))
        call radgrx(nr,nrx,1,rofi,agrh(1,3),grgagr(1,3))
        do  34  ir = 1, nr
   34   grgagr(ir,3) = (grh(ir,1)+grh(ir,2))*grgagr(ir,3)
        do  36  ir = 1, nr
   36   agrh(ir,4) = grh(ir,1)*grh(ir,2)
      endif

      call poppr

C --- Call PBE subroutine for each radial point ---
 9999 continue
      do  40  ir=1,nr
      if (mode.ne.0.and.mode.ne.10) then
        Xagrh1  =agrh  (ir,1)
        Xgrgagr1=grgagr(ir,1)
        Xggrh1  =ggrh  (ir,1)
        Xagrh2  =agrh  (ir,nsp)
        Xgrgagr2=grgagr(ir,nsp)
        Xggrh2  =ggrh  (ir,nsp)
        XXagrh  =agrh  (ir,2*nsp-1)
        XXgrgagr=grgagr(ir,2*nsp-1)
      else
        Xagrh1  =0d0
        Xgrgagr1=0d0
        Xggrh1  =0d0
        Xagrh2  =0d0
        Xgrgagr2=0d0
        Xggrh2  =0d0
        XXagrh  =0d0
        XXgrgagr=0d0
      endif
      call easypbe(rp    (ir,1),   ! rho up
     .             Xagrh1,         ! grad (rho_up)
     .             Xgrgagr1,       ! (grad up).(grad |grad up|)
     .             Xggrh1,         ! Laplacian of up
     .             rp    (ir,nsp), ! rho dwn
     .             Xagrh2,         ! grad (rho_dwn)
     .             Xgrgagr2,       !(grad dwn).(grad |grad dwn|)
     .             Xggrh2,         ! Laplacian of dwn
     .             XXagrh,         ! |grad rho|
     .             XXgrgagr,       ! (grad rho).(grad |grad rho|)
     .             1,                ! do correlations
     .             1,                ! do potential
     .             exlsd,  ! LSD exchange energy density
     .             vxuplsd,! up LSD exchange potential
     .             vxdnlsd,! down LSD exchange potential
     .             eclsd, ! LSD exchange-correlation energy density
     .             vcuplsd,! up LSD exchange-correlation potential
     .             vcdnlsd,! down LSD exchange-correlation potential
     .             expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, !PW91
     .             expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)       !PBE

c-  LSDA part -- if mode >= 0 !!! Remember output in Hartree!!
      if(mode.eq.0.or.mode.eq.10) then
        exc(ir)=exlsd+eclsd
        vxc(ir,1)=vXuplsd+vCuplsd
        if(nsp.eq.2) vxc(ir,2)=vXdnlsd+vCdnlsd
      endif

c- PW91 part
      if(mode.eq.1) then                ! PW91 GGA
        exc(ir)=expw91+ecpw91
        vxc(ir,1)=vxuppw91+vcuppw91
        if(nsp.eq.2) vxc(ir,2)=vxdnpw91+vcdnpw91
      endif

c- PBE part
      if (mode.eq.2) then          !PBE GGA
        exc(ir)=expbe+ecpbe
        vxc(ir,1)=vxuppbe+vcuppbe
        if(nsp.eq.2) vxc(ir,2)=vxdnpbe+vcdnpbe
      endif

c- PW correction
      if (mode.eq.-1) then
        exc(ir)=expw91+ecpw91-exlsd-eclsd
        vxc(ir,1)=vxuppw91+vcuppw91-vXuplsd-vCuplsd
        if(nsp.eq.2) vxc(ir,2)=vxdnpw91+vcdnpw91-vXdnlsd-vCdnlsd
      endif

c- PBE correction
      if (mode.eq.-2) then
        exc(ir)=expbe+ecpbe-exlsd-eclsd
        vxc(ir,1)=vxuppbe+vcuppbe-vXuplsd-vCuplsd
        if(nsp.eq.2) vxc(ir,2)=vxdnpbe+vcdnpbe-vXdnlsd-vCdnlsd
      endif

   40 continue

      if (units .eq. 'Ry') then
        do  50  ir=1,nr
          exc(ir)=exc(ir)*2
          do  52  isp=1,nsp
   52     vxc(ir,isp)=vxc(ir,isp)*2
   50   continue
      endif

      if (mode.ne.10) then
        do  60  i = 1, nsp
   60   vxc(1,i)=(vxc(2,i)*rofi(3)-vxc(3,i)*rofi(2))/(rofi(3)-rofi(2))
      endif
      end

