      subroutine vxnlcc(n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom,
     .grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
C- Same as vxnloc but exc also cut off to avoid blowup in exc
C     implicit none
C Passed parameters
      integer nsp,n,i
      double precision rhop(1),rhom(1),grhop(1),grhom(1),grho(1),
     .ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1),
     .vxc1(1),vxc2(1),exc(1)
C Local parameters
C     logical warned
      double precision pi,fivth,th4,th2,th,sth,sevni,f8,f9,rho,gp2,gm2,
     .bigf,aa,d,polar,f,cutoff,cutof1,gro,ggro,g2,
     .cutof2,hh,rp,rm,grp,grm,ggrp,ggrm,rmin,xd,expf,xx
      parameter (f=0.15d0,hh=1d-3,fivth=5d0/3d0,th4=4d0/3d0,
     .th2=2d0/3d0,th=1d0/3d0,sth=7d0/3d0,sevni=7d0/9d0, rmin=1d-15)

      call tcn('vxnlcc')
C     warned = .false.
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
            exc(i) = exc(i) + aa*(gro*gro/(rho**sth))*(expf-sevni)*cutoff
            vxc1(i) = vxc1(i) + 2d0*aa*rho**(-th)*
     .      (sevni*(ggro/rho-th2*gro*gro/(rho*rho)) - expf*(
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
C       if (rp .le. 0d0)
C    .    print *, 'vxnloc (warning) rho+ not positive but rho is'
C       if (rm .le. 0d0)
C    .    print *, 'vxnloc (warning) rho- not positive but rho is'
            if (rp .le. 0d0) rp = rho/1000
            if (rm .le. 0d0) rm = rho/1000
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

            cutof1 = dexp(-hh*(gp2/(rp*rp))*(2*rp)**(-th2))
            xd = f9*rho**(th-4d0)/(d*d)*(rp**th2-rm**th2)
            xx = (2d0-bigf)*ggro/rho -
     .      (th4+bigf*(-11d0/3d0+bigf*7d0/6d0))*g2/(rho*rho) +
     .      bigf*(bigf-3d0)*grggr(i)/(gro*rho+1d-16)
            vxc1(i) = vxc1(i) + aa/rho**th*(
     .      -sevni/(2d0**th)*(rho/rp)**th*(th4*gp2/(rp*rp)-2d0*ggrp/rp)
     .      -expf*(xx -
     .      xd*((1d0-bigf)*rm*g2-(2d0-bigf)*rho*(grpgrm(i)+gm2))))
     .      *cutof1
            cutof2=dexp(-hh*(gm2/(rm*rm))*(2*rm)**(-th2))
            vxc2(i) = vxc2(i) + aa/rho**th*(
     .      -sevni/(2d0**th)*(rho/rm)**th*(th4*gm2/(rm*rm)-2d0*ggrm/rm)
     .      -expf*(xx +
     .      xd*((1d0-bigf)*rp*g2-(2d0-bigf)*rho*(gp2+grpgrm(i)))))
     .      *cutof2
            exc(i) = exc(i) + aa/rho*(expf*g2/(rho**th4)
     .      - sevni/(2d0**th)*(gp2/(rp**th4)+gm2/(rm**th4)))
     .      *dsqrt(cutof1)*dsqrt(cutof2)
          endif
   20   continue
      endif
      call tcx('vxnlcc')
      end

