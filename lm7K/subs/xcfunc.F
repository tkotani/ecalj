      logical function testwritexcfun()
      testwritexcfun=.false.
      end

C
C This file includes kinds of the xc functionals in LDA and GGA.
C
C I think no other files contained such functionals in lm7K package.
C
C For given input density (and their gradients), this returns Exc and its gradients.
C The easypbe.F, K. Burke, May 14, 1996. is contained at the bottom of this file.
C
C Two kinds of evxcv is included. One is evxcv_new (lmf7.2), the other is evxcv_old (lmf7.1).
C They caused a little idfferences.
C===================================================================
      subroutine evxcv(rho,rhosp,n,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
c      subroutine evxcv_new(rho,rhosp,n,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
C- XC energy density and potential for a vector of points.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rho   :spin-1 + spin-2 density
Ci   rhosp :spin-1 density (unused for nsp=1)
Ci   n     :number of points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lxcf  :local exchange-correlation functional index
Ci         :1= Ceperly Alder
Ci         :2= Hedin-Lundqvist
Ci         :3,4= PW91
Co Outputs
Co   exc   :local exchange energy density for the n points
Co   vxc   :local exchange potential for the n points
Co         := d (rho exc) /rho)
Co   ex    :exchange-only part of exc
Co   vx    :exchange-only part of vxc
Co   ec    :correlation-only part of exc
Co   vc    :correlation-only part of vxc
Cr Remarks
Cu Updates
Cu   21 Nov 09 criterion to evaluate exc,vxc: rho>rhomin
Cu             No longer checks individual spin parts
Cu   21 Apr 09 Calls evxcp for lxcf=3,4
Cu    8 Feb 02 vx and vc (T. Miyake)
Cu   14 Jan 02 ex and ec (T. Miyake)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,nsp,lxcf
      double precision rho(n),rhosp(n),exc(n),ex(n),ec(n),
     .vxc(n),vx(n),vc(n)
C ... Local parameters
C --- Vosko-Ceperley-Alder ---
      double precision ap,xp0,bp,cp,qp,cp1,cp2,cp3,
     .af,xf0,bf,cf,qf,cf1,cf2,cf3
      double precision wk(n),wk2(n)
      parameter (ap=.0621814d0,bp=3.72744d0,qp=6.1519908d0,
     .cp=12.9352d0,cp1=1.2117833d0,cp2=1.1435257d0,
     .cp3=-.031167608d0,xp0=-.10498d0)
      parameter (af=.0310907d0,bf=7.06042d0,qf=4.7309269d0,
     .cf1=2.9847935d0,cf2=2.7100059d0,cf3=-.1446006d0,
     .cf=18.0578d0,xf0=-.32500d0)
      integer i
      double precision atnf,atnp,beta,dbeta,dfs,duc,vxc1,exc0,excf,excp,
     .fs,rs,s,s4,tf1,th,th4,thx,tp1,ucf,ucp,vxc0,x,xfx,xpx,sth
      parameter (th=1d0/3d0, th4=4d0/3d0, thx=0.519842099789746d0)

C --- Hedin-Lundqvist ---
      double precision mucf,mucp,nuce,unthrd,fothrd,fpi3,
     .cph,cfh,aph,afh,fh,fhofxi,z,xi,cmuxp,efmep,epscf,epscp,
     .epsx0,epsx1,gamma,tauce,rmin
      parameter(unthrd=1d0/3d0,fothrd=4d0/3d0,rmin=1d-21,
     .fpi3=12.566370614d0/3,gamma=5.1297628d0,cmuxp=1.2217741d0,
     .epsx0=.9163305866d0,epsx1=.23817361d0)
C    .  fpi3=4.1887902d0,gamma=5.1297628d0,cmuxp=1.2217741d0,
C ... V. Barth, Hedin parametrization ---
C     parameter CPH=-.0504D0,CFH=-.0254D0,APH=30d0,AFH=75D0)
C --- Taken from ASW ---
      parameter (cph=-.045d0,cfh=-.0225d0,aph=21d0,afh=52.9167d0)
      fh(z)     = (1d0+z*z*z)*dlog(1d0+1d0/z)+.5d0*z-z*z-unthrd

      if (lxcf .gt. 1) goto 200
C --- Vosko-Ceperley-Alder, spin polarized case ---
      if (nsp .eq. 2) then
        do  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = dsqrt(rs)
            xpx = x*x+bp*x + cp
            atnp = datan(qp/(2d0*x+bp))
            excp = ap*(dlog(x*x/xpx) + cp1*atnp
     .      -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
            tp1  = (x*x+bp*x)/xpx
            ucp  = excp - ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))

            xfx = x*x+bf*x + cf
C           s = 2*rhosp/rho-1 
            s = min(max(2*rhosp(i)/rho(i)-1d0,-1d0),1d0)
C           sth -> (2*rhosp/rho)^1/3 as rhosp->0
            sth = (s+1d0)**th
            s4 = s**4 - 1d0
            fs = (sth**4 + (1d0-s)**th4 - 2d0) / thx
            beta = 1d0/(1.74208d0 +x*(3.182d0 +x*(.09873d0+x*.18268d0)))
            dfs = th4*(sth - (1d0-s)**th)/thx
            dbeta = -(.27402d0*x + .09873d0 + 1.591d0/x)*beta**2
            atnf = datan(qf/(2d0*x + bf))
            excf = af*(dlog(x*x/xfx) + cf1*atnf
     .      -cf3*(dlog((x-xf0)**2/xfx) + cf2*atnf))
            exc0 = excp + fs*(excf-excp)*(1d0+s4*beta)
            tf1  = (x*x+bf*x)/xfx
            ucf  = excf-af/3d0*(1d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
            vxc0 = (ucf-ucp)*fs - (excf-excp)*(s-1d0)*dfs
            duc  = (ucf-ucp)*beta*s4*fs +
     .      (excf-excp)*(-rs/3d0)*dbeta*s4*fs
            vxc1 = duc - (excf-excp)*beta*(s-1)*(4d0*s**3*fs+s4*dfs)
            vxc(i) = ucp - 1.2217736d0/rs*sth + vxc0 + vxc1
            vx(i)  =     - 1.2217736d0/rs*sth
            vc(i)  = ucp                      + vxc0 + vxc1
            exc(i) = exc0 - 0.9163306d0/rs - 0.2381735d0/rs*fs
            ex(i)  =      - 0.9163306d0/rs
            ec(i)  = exc0                  - 0.2381735d0/rs*fs
          else
            vxc(i) = 0d0
            vx(i)  = 0d0
            vc(i)  = 0d0
            exc(i) = 0d0
            ex(i)  = 0d0
            ec(i)  = 0d0
          endif
        enddo
      else
        do  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = dsqrt(rs)
            xpx = x*x+bp*x + cp
            atnp = datan(qp/(2d0*x+bp))
            excp = ap*(dlog(x*x/xpx) + cp1*atnp
     .      -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
            tp1  = (x*x+bp*x)/xpx
            ucp  = excp - ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
            vxc(i) = ucp - 1.2217736d0/rs
            vx(i)  =     - 1.2217736d0/rs
            vc(i)  = ucp
            exc(i) = excp - 0.9163306d0/rs
            ex(i)  =      - 0.9163306d0/rs
            ec(i)  = excp
          else
            vxc(i) = 0d0
            vx(i)  = 0d0
            vc(i)  = 0d0
            exc(i) = 0d0
            ex(i)  = 0d0
            ec(i)  = 0d0
          endif
        enddo
      endif
      return

C --- Barth-Hedin ---
  200 continue
      if (lxcf .gt. 2) goto 300
      if (nsp .eq. 2) then
        do  i = 1, n
C         if (rho(i) .gt. rmin .and. rho(i) .ge. rhosp(i)) then
          if (rho(i) .gt. rmin) then
C            if (rho(i) .lt. rhosp(i)) then
C              print *, 'hi'
C            endif
            rs = (rho(i)*fpi3)**(-th)
            x = rs/aph
            mucp = cph*dlog(1d0 + 1d0/x)
            efmep = -cph*fh(x)
            epscp = -efmep
            x = rs/afh
            epscf = cfh*fh(x)
            efmep = efmep + epscf
            nuce = gamma*efmep
            mucf = cfh*dlog(1d0 + 1d0/x)
            tauce = mucf - mucp - fothrd*efmep
            xi = min(max(rhosp(i)/rho(i),0d0),1d0)
            fhofxi = (xi**fothrd + (1d0-xi)**fothrd -.79370052598d0)/
     .      .206299474d0
            vxc(i)  =  mucp + tauce*fhofxi
     .      + (nuce - cmuxp/rs)*((2d0*xi)**unthrd) - nuce
            vx(i) =  - cmuxp/rs*((2d0*xi)**unthrd)
            vc(i) =  mucp + tauce*fhofxi
     .      +  nuce*((2d0*xi)**unthrd)             - nuce
            exc(i)  = epscp - epsx0/rs + fhofxi*(epscf-epscp-epsx1/rs)
            ex(i) =         - epsx0/rs
            ec(i) = epscp + fhofxi*(epscf-epscp-epsx1/rs)
          else
            vxc(i) = 0d0
            vx(i)  = 0d0
            vc(i)  = 0d0
            exc(i) = 0d0
            ex(i)  = 0d0
            ec(i)  = 0d0
          endif
        enddo
      else
        do  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = rs/aph
            mucp = cph*dlog(1d0 + 1d0/x)
            efmep = -cph*fh(x)
            epscp = -efmep
            vxc(i) = mucp - cmuxp/rs
            vx(i)  =      - cmuxp/rs
            vc(i)  = mucp
            exc(i) = epscp - epsx0/rs
            ex(i)  = - epsx0/rs
            ec(i)  = epscp
          else
            vxc(i) = 0d0
            vx(i)  = 0d0
            vc(i)  = 0d0
            exc(i) = 0d0
            ex(i)  = 0d0
            ec(i)  = 0d0
          endif
        enddo
      endif
      return

C --- PW91, PBE ---
  300 continue
      if (lxcf .gt. 4) goto 400
      if (nsp .eq. 1) then
        call evxcp(rho,rho,n,nsp,lxcf,ex,ec,exc,vx,wk,vc,wk,vxc,wk)
      else
        wk2 = rho - rhosp
        call evxcp(rhosp,wk2,n,nsp,lxcf,ex,ec,exc,vx,wk,vc,wk,vxc,wk)
      endif
      return

C --- No local functional ---
  400 continue
      call setpr(30)
      call rxi('evxcv: no functional, lxcf =',lxcf)
      end

C      integer function excsan(lxcfun,ifi)
CC- Sanity checks for a specified xc switch, optional printout
CC     implicit none
C      integer lxcfun,ifi
C      integer lxcf,lxcg,k,i1,i2
C      character *80 sout
C
C      lxcf = mod(mod(lxcfun,100),3)
C      lxcg = lxcfun/100
C
C      k = 0
C      if (lxcf .lt. 1 .or. lxcf .gt. 2) k = -1
C      if (lxcg .ne. 0) k = k - 10
C      excsan = k
C      if (ifi .le. 0) return
C
C      if (k .lt. 0) call rxi('evxcv: unknown XC code:',lxcfun)
C      if (lxcf .eq. 1) write(sout,101)
C      if (lxcf .eq. 2) write(sout,102)
C  101 format(' XC potential is Ceperly-Alder/Vosko')
C  102 format(' XC potential is Barth-Hedin')
C      call strip(sout,i1,i2)
C      if (lxcg .eq. 0) write(sout(i2+1:),201) 'no'
C      if (lxcg .gt. 0) write(sout(i2+1:),201) 'with'
C  201 format(' (',a,' gradient corrections)')
C      call strip(sout,i1,i2)
C      write(ifi,'(a)') sout(1:i2)
C
C      end



c for older version of smvxcm.F ==============================
c      subroutine evxcv(rho,rhosp,n,nsp,lxcf,exc,excx,excc,vxc,vxcx,vxcc)
      subroutine evxcv_old(rho,rhosp,n,nsp,lxcf,exc,excx,excc,vxc,vxcx,vxcc)
C- XC energy density and potential for a vector of points.
C ----------------------------------------------------------------------
C!!! This routine has a bug when rho>0 and rhosp=0  for XCFUN=2
Ci Inputs
Ci   rho   :spin-1 + spin-2 density
Ci   rhosp :spin-1 density (unused for nsp=1)
Ci   n     :number of points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lxcf  :local exchange-correlation functional index
Ci         :1= Ceperly Alder
Ci         :2= Hedin-Lundqvist
Co Outputs
Co   exc   :local exchange energy density for the n points
Co   vxc   :local exchange potential for the n points
Co         := d (rho exc) /rho)
Co   excx  :exchange-only part of exc
Co   vxcx  :exchange-only part of vxc
Co   excc  :correlation-only part of exc
Co   vxcc  :correlation-only part of vxc
Cr Remarks
Cu Updates
Cu   25 Jul 08 Fixed bug in CA functional when rho>0 and rhosp=0
Cu    8 Feb 02 vxcx and vxcc (T. Miyake)
Cu   14 Jan 02 excx and excc (T. Miyake)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,nsp,lxcf
      double precision rho(n),rhosp(n),exc(n),excx(n),excc(n),
     .vxc(n),vxcx(n),vxcc(n)
C ... Local parameters
C --- Vosko-Ceperley-Alder ---
      double precision ap,xp0,bp,cp,qp,cp1,cp2,cp3,
     .af,xf0,bf,cf,qf,cf1,cf2,cf3
      parameter (ap=.0621814d0,bp=3.72744d0,qp=6.1519908d0,
     .cp=12.9352d0,cp1=1.2117833d0,cp2=1.1435257d0,
     .cp3=-.031167608d0,xp0=-.10498d0)
      parameter (af=.0310907d0,bf=7.06042d0,qf=4.7309269d0,
     .cf1=2.9847935d0,cf2=2.7100059d0,cf3=-.1446006d0,
     .cf=18.0578d0,xf0=-.32500d0)
      integer i
      double precision atnf,atnp,beta,dbeta,dfs,duc,vxc1,exc0,excf,excp,
     .fs,rs,s,s4,tf1,th,th4,thx,tp1,ucf,ucp,vxc0,x,xfx,xpx,sth
      parameter (th=1d0/3d0, th4=4d0/3d0, thx=0.519842099789746d0)

C --- Hedin-Lundqvist ---
      double precision mucf,mucp,nuce,unthrd,fothrd,fpi3,
     .cph,cfh,aph,afh,fh,fhofxi,z,xi,cmuxp,efmep,epscf,epscp,
     .epsx0,epsx1,gamma,tauce,rmin
      parameter(unthrd=1d0/3d0,fothrd=4d0/3d0,rmin=1d-21,
     .fpi3=12.566370614d0/3,gamma=5.1297628d0,cmuxp=1.2217741d0,
     .epsx0=.9163305866d0,epsx1=.23817361d0)
C    .  fpi3=4.1887902d0,gamma=5.1297628d0,cmuxp=1.2217741d0,
C ... V. Barth, Hedin parametrization ---
C     parameter CPH=-.0504D0,CFH=-.0254D0,APH=30d0,AFH=75D0)
C --- Taken from ASW ---
      parameter (cph=-.045d0,cfh=-.0225d0,aph=21d0,afh=52.9167d0)
      fh(z)     = (1d0+z*z*z)*dlog(1d0+1d0/z)+.5d0*z-z*z-unthrd

      if (lxcf .gt. 1) goto 200
C --- Vosko-Ceperley-Alder, spin polarized case ---
      if (nsp .eq. 2) then
        do  110  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = dsqrt(rs)
            xpx = x*x+bp*x + cp
            atnp = datan(qp/(2d0*x+bp))
            excp = ap*(dlog(x*x/xpx) + cp1*atnp
     .      -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
            tp1  = (x*x+bp*x)/xpx
            ucp  = excp - ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))

            xfx = x*x+bf*x + cf
C           s = 2*rhosp/rho-1
            s = min(max(2*rhosp(i)/rho(i)-1d0,-1d0),1d0)
C           sth -> (2*rhosp/rho)^1/3 as rhosp->0
            sth = (s+1d0)**th
            s4 = s**4 - 1d0
            fs = (sth**4 + (1d0-s)**th4 - 2d0) / thx
            beta = 1d0/(1.74208d0 +x*(3.182d0 +x*(.09873d0+x*.18268d0)))
            dfs = th4*(sth - (1d0-s)**th)/thx
            dbeta = -(.27402d0*x + .09873d0 + 1.591d0/x)*beta**2
            atnf = datan(qf/(2d0*x + bf))
            excf = af*(dlog(x*x/xfx) + cf1*atnf
     .      -cf3*(dlog((x-xf0)**2/xfx) + cf2*atnf))
            exc0 = excp + fs*(excf-excp)*(1d0+s4*beta)
            tf1  = (x*x+bf*x)/xfx
            ucf  = excf-af/3d0*(1d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
            vxc0 = (ucf-ucp)*fs - (excf-excp)*(s-1d0)*dfs
            duc  = (ucf-ucp)*beta*s4*fs +
     .      (excf-excp)*(-rs/3d0)*dbeta*s4*fs
            vxc1 = duc - (excf-excp)*beta*(s-1)*(4d0*s**3*fs+s4*dfs)
            vxc(i)  = ucp - 1.2217736d0/rs*sth + vxc0 + vxc1
            vxcx(i) =     - 1.2217736d0/rs*sth
            vxcc(i) = ucp                      + vxc0 + vxc1
            exc(i)  = exc0 - 0.9163306d0/rs - 0.2381735d0/rs*fs
            excx(i) =      - 0.9163306d0/rs
            excc(i) = exc0                  - 0.2381735d0/rs*fs
          else
            vxc(i) = 0d0
            vxcx(i) = 0d0
            vxcc(i) = 0d0
            exc(i) = 0d0
            excx(i) = 0d0
            excc(i) = 0d0
          endif
  110   continue
      else
        do  120  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = dsqrt(rs)
            xpx = x*x+bp*x + cp
            atnp = datan(qp/(2d0*x+bp))
            excp = ap*(dlog(x*x/xpx) + cp1*atnp
     .      -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
            tp1  = (x*x+bp*x)/xpx
            ucp  = excp - ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
            vxc(i)  = ucp - 1.2217736d0/rs
            vxcx(i) =     - 1.2217736d0/rs
            vxcc(i) = ucp
            exc(i)  = excp - 0.9163306d0/rs
            excx(i) =      - 0.9163306d0/rs
            excc(i) = excp
          else
            vxc(i) = 0d0
            vxcx(i) = 0d0
            vxcc(i) = 0d0
            exc(i) = 0d0
            excx(i) = 0d0
            excc(i) = 0d0
          endif
  120   continue
      endif
      return

C --- Barth-Hedin ---
  200 continue
      if (nsp .eq. 2) then
        do  210  i = 1, n
          if (rho(i) .gt. rmin .and. rho(i) .ge. rhosp(i)) then
            rs = (rho(i)*fpi3)**(-th)
            x = rs/aph
            mucp = cph*dlog(1d0 + 1d0/x)
            efmep = -cph*fh(x)
            epscp = -efmep
            x = rs/afh
            epscf = cfh*fh(x)
            efmep = efmep + epscf
            nuce = gamma*efmep
            mucf = cfh*dlog(1d0 + 1d0/x)
            tauce = mucf - mucp - fothrd*efmep
            xi = max(rhosp(i)/rho(i),0d0)
            fhofxi = (xi**fothrd + (1d0-xi)**fothrd -.79370052598d0)/
     .      .206299474d0
            vxc(i)  =  mucp + tauce*fhofxi
     .      + (nuce - cmuxp/rs)*((2d0*xi)**unthrd) - nuce
            vxcx(i) =  - cmuxp/rs*((2d0*xi)**unthrd)
            vxcc(i) =  mucp + tauce*fhofxi
     .      +  nuce*((2d0*xi)**unthrd)             - nuce
            exc(i)  = epscp - epsx0/rs + fhofxi*(epscf-epscp-epsx1/rs)
            excx(i) =       - epsx0/rs
            excc(i) = epscp + fhofxi*(epscf-epscp-epsx1/rs)
          else
            vxc(i) = 0d0
            vxcx(i) = 0d0
            vxcc(i) = 0d0
            exc(i) = 0d0
            excx(i) = 0d0
            excc(i) = 0d0
          endif
  210   continue
      else
        do  220  i = 1, n
          if (rho(i) .gt. rmin) then
            rs = (rho(i)*fpi3)**(-th)
            x = rs/aph
            mucp = cph*dlog(1d0 + 1d0/x)
            efmep = -cph*fh(x)
            epscp = -efmep
            vxc(i)  = mucp - cmuxp/rs
            vxcx(i) =      - cmuxp/rs
            vxcc(i) = mucp
            exc(i)  = epscp - epsx0/rs
            excx(i) = - epsx0/rs
            excc(i) = epscp
          else
            vxc(i) = 0d0
            vxcx(i) = 0d0
            vxcc(i) = 0d0
            exc(i) = 0d0
            excx(i) = 0d0
            excc(i) = 0d0
          endif
  220   continue
      endif

      end

C$$$      integer function excsan(lxcfun,ifi)
C$$$C- Sanity checks for a specified xc switch, optional printout
C$$$C     implicit none
C$$$      integer lxcfun,ifi
C$$$      integer ltype,lgrad,k,i1,i2
C$$$      character *80 sout

C$$$      ltype = mod(mod(lxcfun,10),3)
C$$$      lgrad = lxcfun/10

C$$$      k = 0
C$$$      if (ltype .lt. 1 .or. ltype .gt. 2) k = -1
C$$$      if (lgrad .ne. 0) k = k - 10
C$$$      excsan = k
C$$$      if (ifi .le. 0) return

C$$$      if (k .lt. 0) call rxi('evxcv: unknown XC code:',lxcfun)
C$$$      if (ltype .eq. 1) write(sout,101)
C$$$      if (ltype .eq. 2) write(sout,102)
C$$$  101 format(' XC potential is Ceperly-Alder/Vosko')
C$$$  102 format(' XC potential is Barth-Hedin')
C$$$      call strip(sout,i1,i2)
C$$$      if (lgrad .eq. 0) write(sout(i2+1:),201) 'no'
C$$$      if (lgrad .gt. 0) write(sout(i2+1:),201) 'with'
C$$$  201 format(' (',a,' gradient corrections)')
C$$$      call strip(sout,i1,i2)
C$$$      write(ifi,'(a)') sout(1:i2)

C$$$      end



C===================================================================
      subroutine evxcp(r1,r2,n,nsp,lxcf,ex,ec,exc,vxup,vxdn,vcup,vcdn,
     .vxcup,vxcdn)
C-Local part xc energy and potential for PW91 or PBE (in LDA, PBE=PW91)
C ----------------------------------------------------------------------
Ci Inputs: r1, r2 charge density in a mesh of n points. If nsp=1, r1
Ci         is total charge; if nsp=2 r1 and r2 are up and 
Ci         down charge.
Ci         lxcf=3,4 for PBE (AKA PW91)
Co Outputs:
Co   ex    :exchange-only part of exc
Co   ec    :correlation-only part of exc
Co   exc   :local exchange energy density for the n points
Co   vxup  :exchange-only part of vxc, spin1
Co   vcup  :correlation-only part of vxc, spin1
Co   vxcup :local exchange potential for the n points, spin1
Co   vxdn  :exchange-only part of vxc, spin2
Co   vcdn  :correlation-only part of vxc, spin2
Co   vxcdn :local exchange potential for the n points, spin2
Co         If nsp=1, only spin numbers are returned
Cr Remarks: 
Cr         If nsp=1 r2 points to an arbitrary address. NB evxcp calls
Cr         easypbe written by K. Burke, which uses Hartree atomic units.
Cr         This routine adapted from pbevxc in FP-LMTO and uses names
Cr         longer than six characters.
Cu Updates
Cu   21 Apr 09 Returns exchange and correlation parts separately
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsp,n,lxcf
      double precision r1(n),r2(n),ex(n),ec(n),exc(n)
      double precision vxup(n),vxdn(n),vcup(n),vcdn(n),vxcup(n),vxcdn(n)
C Local Variables
      integer i,iprint
      double precision up, dn, exlsd, vxuplsd, vxdnlsd, eclsd,
     .vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91,
     .ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe,
     .vxdnpbe, ecpbe, vcuppbe, vcdnpbe

      if (lxcf .ne. 3 .and. lxcf .ne. 4) then
        if (iprint() .lt. 10) call pshpr(10)
        call rxi('evxcp cannot handle lxcf =',lxcf)
      endif

      do  i = 1, n
        if (nsp .eq. 1) then
          up = r1(i) / 2
          dn = up
        else
          up = r1(i)
          dn = r2(i)
        endif
        call easypbe(up,0d0,0d0,0d0,dn,0d0,0d0,0d0,0d0,0d0,1,1,0,
     .  exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe) 
C   ... Times two to convert Ha to Ry
        ex(i) = exlsd * 2d0
        ec(i) = eclsd * 2d0
        exc(i) = (exlsd + eclsd) * 2d0
        vxup(i)  = vxuplsd * 2d0
        vcup(i)  = vcuplsd * 2d0
        vxcup(i) = (vxuplsd + vcuplsd) * 2d0
        if (nsp .eq. 2) then
          vxdn(i)  = vxdnlsd * 2d0
          vcdn(i)  = vcdnlsd * 2d0
          vxcdn(i) = (vxdnlsd + vcdnlsd) * 2d0
        endif
      enddo
      end


      subroutine vxcgga(lxcg,n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom,
     .grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
C- PW91 and PBE gradient corrections to exc and vxc for a set of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rhop  :spin up density if nsp=2; otherwise total density
Ci   rhom  :spin down density if nsp=2; otherwise not used
Ci   grhop :|grad rhop| or |grad rho| if nsp=1
Ci   grhom :|grad rhom| (nsp=2)
Ci   ggrhop:Laplacian of rhop
Ci   ggrhom:Laplacian of rhom
Ci   grho  :|grad total rho| if nsp=2; otherwise not used
Ci   grpgrm:grad rho+ . grad rho- (not used here)
Ci   grggr :(grad rho).(grad |grad rho|)
Ci   grggrp:(grad up).(grad |grad up|) if nsp=2; otherwise ggrgr
Ci   grggrm:(grad dn).(grad |grad dn|) if nsp=2; otherwise not used
Co Outputs
Co   vxc1  :GGA contr. to vxc+ if nsp=2; otherwise GGA contr. to vxc
Co   vxc2  :GGA contr. to vxc- if nsp=2; otherwise GGA contr. to vxc
Co    exc  :GGA contr. to exc
Cl Local variables
Cl   symbols match easypbe
Cb Bugs
Cr Remarks
Cu Updates
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsp,n,lxcg
      double precision rhop(1),rhom(1),grhop(1),grhom(1),grho(1),
     .ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1),
     .vxc1(1),vxc2(1),exc(1)
C Local Variables
      integer i,isw
      double precision up,dn,agrup,delgrup,uplap,agrdn,delgrdn,dnlap,
     .agr,delgr
      double precision exlsd,vxuplsd,vxdnlsd,eclsd,
     .vcuplsd,vcdnlsd,expw91,vxuppw91,vxdnpw91,
     .ecpw91,vcuppw91,vcdnpw91,expbe,vxuppbe,
     .vxdnpbe,ecpbe,vcuppbe,vcdnpbe

      logical:: testwritexcfun
      if (lxcg .ne. 3 .and. lxcg .ne. 4) call rx('evxcp: bad lxcg')

      do  i = 1, n
        if (nsp .eq. 1) then
          up = rhop(i) / 2
          dn = up
          agrup = grhop(i) / 2
          agrdn = agrup
          uplap = ggrhop(i) / 2
          dnlap = uplap
          delgrup = grggrp(i) / 4
          delgrdn = delgrup
          agr = grhop(i)
          delgr = grggrp(i)
        else
          up = rhop(i)
          dn = rhom(i)
          agrup = grhop(i)
          agrdn = grhom(i)
          uplap = ggrhop(i)
          dnlap = ggrhom(i)
          delgrup = grggrp(i)
          delgrdn = grggrm(i)
          agr = grho(i)
          delgr = grggr(i)
        endif
cccccccccccccccccccccc
c        write(1219,*)' goto easypbe i lxcg=',i,lxcg
cccccccccccccccccccccc
        call easypbe(up,agrup,delgrup,uplap,
     .  dn,agrdn,delgrdn,dnlap,agr,delgr,1,1,
     .  isw(lxcg.eq.4),
     .  exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

        if (lxcg .eq. 2) then
          exc(i) = exc(i) + (expw91 + ecpw91 - exlsd - eclsd) * 2d0
          vxc1(i) = vxc1(i) +
     .    (vxuppw91 + vcuppw91 - vxuplsd - vcuplsd) * 2d0
          if (nsp .eq. 2) vxc2(i) = vxc2(i) +
     .    (vxdnpw91 + vcdnpw91 - vxdnlsd - vcdnlsd) * 2d0
        else
          exc(i) = exc(i) + (expbe + ecpbe - exlsd - eclsd) * 2d0
          vxc1(i) = vxc1(i) +
     .    (vxuppbe + vcuppbe - vxuplsd - vcuplsd) * 2d0

ccccccccccccccccccccc
          if(testwritexcfun()) then
          aaa=(vxuppbe + vcuppbe - vxuplsd - vcuplsd)
          bbb=(vxdnpbe + vcdnpbe - vxdnlsd - vcdnlsd)
c          if(abs(aaa)>1.or.abs(bbb)>1) then
             write(1219,"('ggaxxx ',i8,255d15.6)")i,up,agrup,delgrup,uplap
             write(1219,"('ggaxxx ',i8,255d15.6)")i,dn,agrdn,delgrdn,dnlap
             write(1219,"('ggaxxx ',i8,255d15.6)")i,agr,delgr
             write(1219,"('ggaxxx ',i8,255f12.3)")i,aaa, vxuppbe, vcuppbe, vxuplsd, vcuplsd
             write(1219,"('ggaxxx ',i8,255f12.3)")i,bbb, vxdnpbe, vcdnpbe, vxdnlsd, vcdnlsd
             write(1219,*)
c          endif  
          endif
ccccccccccccccccccccccccc

          if (nsp .eq. 2) vxc2(i) = vxc2(i) +
     .    (vxdnpbe + vcdnpbe - vxdnlsd - vcdnlsd) * 2d0
        endif
      enddo

ccccccccccccccccccccccccc
      if(testwritexcfun()) then
      stop 'xxxxxxxxx end of vxcgga xxxxxxxxxxxxxxxxxxxxxxxx'
      endif
ccccccccccccccccccccccccc

      end


      subroutine vxnloc(n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom,
     .grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
C- Langreth-Mehl-Hu gradient correction to exc and vxc
C ---------------------------------------------------------------
Ci Inputs:
Ci   rhop  :spin up density if nsp=2; otherwise total density
Ci   rhom  :spin down density if nsp=2; otherwise not used
Ci   grhop :|grad rhop| or |grad rho| if nsp=1
Ci   grhom :|grad rhom| (nsp=2)
Ci   grho  :|grad total rho| if nsp=2; otherwise not used
Ci   ggrhop:Laplacian of rhop
Ci   ggrhom:Laplacian of rhom
Ci   grggr :(grad rho).(grad |grad rho|)
Ci   grggrp:(grad up).(grad |grad up|) (not used here)
Ci   grggrm:(grad dn).(grad |grad dn|) (not used here)
Ci   grpgrm: grad rho+ . grad rho- (nsp=2)
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
     .ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1),
     .vxc1(1),vxc2(1),exc(1)
C Local parameters
      logical warned
      double precision pi,fivth,th4,th2,th,sth,sevni,f8,f9,rho,gp2,gm2,
     .bigf,aa,d,polar,f,cutoff,cutof1,gro,ggro,g2,
     .cutof2,hh,rp,rm,grp,grm,ggrp,ggrm,rmin,xd,expf,xx
      parameter (f=0.15d0,hh=1d-3,fivth=5d0/3d0,th4=4d0/3d0,
     .th2=2d0/3d0,th=1d0/3d0,sth=7d0/3d0,sevni=7d0/9d0, rmin=1d-15)

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
            if (rp .le. 0d0) then
              if (.not. warned)
     .        print *, 'vxnloc (warning) rho+ not positive but rho is'
              rp = rho/1000
              warned = .true.
            endif
            if (rm .le. 0d0) then
              if (.not. warned)
     .        print *, 'vxnloc (warning) rho- not positive but rho is'
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
     .      - sevni/(2d0**th)*(gp2/(rp**th4)+gm2/(rm**th4)))
            cutof1 = dexp(-hh*(gp2/(rp*rp))*(2*rp)**(-th2))
            xd = f9*rho**(th-4d0)/(d*d)*(rp**th2-rm**th2)
            xx = (2d0-bigf)*ggro/rho -
     .      (th4+bigf*(-11d0/3d0+bigf*7d0/6d0))*g2/(rho*rho) +
     .      bigf*(bigf-3d0)*grggr(i)/(gro*rho)
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
          endif
   20   continue
      endif
      end
      subroutine pbevxc(nr,nsp,nrx,rofi,rp,exc,vxc,mode,units,
     .grh,ggrh,agrh,grgagr)
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
     .exc(nrx),vxc(nrx,nsp),rofi(nr),
     .Xagrh1,Xagrh2,Xgrgagr1,Xgrgagr2,Xggrh1,Xggrh2,XXagrh,
     .XXgrgagr

c - Work arrays and internal variables
      double precision grh(nrx,2),ggrh(nrx,2),
     .agrh(nrx,4),grgagr(nrx,3)
      double precision exlsd, vxuplsd, vxdnlsd, eclsd,
     .vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91,
     .ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe,
     .vxdnpbe, ecpbe, vcuppbe, vcdnpbe
      integer ir,i,iprint
      character*2, units

      if(mode.eq.0.or.mode.eq.10) GOTO 9999 !! Do not calculate gradients

      call pshpr(iprint()-20)

C --- grad(rho), laplacian rho ---
      call radgrx(nr,nrx,nsp,rofi,rp,grh)
      call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
      do  20  i  = 1, nsp
        do  24  ir = 2, nr
   24   ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
        ggrh(1,i) =(rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))

C --- grad rho . grad abs grad rho ---
        do  26  ir = 1, nr
   26   agrh(ir,i) = dabs(grh(ir,i))
        call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
        do  28  ir = 1, nr
   28   grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
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
     .  Xagrh1,         ! grad (rho_up)
     .  Xgrgagr1,       ! (grad up).(grad |grad up|)
     .  Xggrh1,         ! Laplacian of up
     .  rp    (ir,nsp), ! rho dwn
     .  Xagrh2,         ! grad (rho_dwn)
     .  Xgrgagr2,       !(grad dwn).(grad |grad dwn|)
     .  Xggrh2,         ! Laplacian of dwn
     .  XXagrh,         ! |grad rho|
     .  XXgrgagr,       ! (grad rho).(grad |grad rho|)
     .  1,                ! do correlations
     .  1,                ! do potential
     .  0,                ! Never do Becke exchange
     .  exlsd,  ! LSD exchange energy density
     .  vxuplsd,! up LSD exchange potential
     .  vxdnlsd,! down LSD exchange potential
     .  eclsd, ! LSD exchange-correlation energy density
     .  vcuplsd,! up LSD exchange-correlation potential
     .  vcdnlsd,! down LSD exchange-correlation potential
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, !PW91
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)       !PBE

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




C===================================================================
C===================================================================
      subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     1           agr,delgr,lcor,lpot,lbecke,
     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c EASYPBE is a driver for the PBE subroutines, using simple inputs
c K. Burke, May 14, 1996.
c inputs: up=up density
c       : agrup=|grad up|
c       : delgrup=(grad up).(grad |grad up|) 
c       : uplap=grad^2 up=Laplacian of up
c       : dn,agrdn,delgrdn,dnlap=corresponding down quantities
c       : agr=|grad rho|
c       : delgr=(grad rho).(grad |grad rho|) 
c       : lcor=flag to do correlation(=0=>don't)
c       : lpot=flag to do potential(=0=>don't)
c       : lbecke=flag to return Becke exchange instead of PW91 (MvS)
c outputs: exlsd=LSD exchange energy density, so that
c               ExLSD=int d^3r rho(r) exlsd(r)
c        : vxuplsd=up LSD exchange potential
c        : vxdnlsd=down LSD exchange potential
c        : exclsd=LSD exchange-correlation energy density
c        : vxcuplsd=up LSD exchange-correlation potential
c        : vxcdnlsd=down LSD exchange-correlation potential
c        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
c        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
Cr DLN,ATP,MvS have made a few trivial changes (see below)
Cr NB: in the absence of gradient corrections, LSD=PW91=PBE
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c needed constants:
c pi32=3 pi**2
c alpha=(9pi/4)**thrd
      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE exchange
c use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
c do up exchange
c fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3) 
c s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
c u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
c v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
      eps = dsqrt(d1mach(3))
      rho2=2.d0*up
      if(rho2.gt.eps)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        call exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
        call exchpw91(lbecke,rho2,s,u,v,exuppw91,vxuppw91)
        call exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
      else
        exuplsd=0.d0
        vxuplsd=0.d0
        exuppw91=0.d0
        vxuppw91=0.d0
        exuppbe=0.d0
        vxuppbe=0.d0
      endif
c repeat for down
      rho2=2.d0*dn
      if(rho2.gt.eps)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        call exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
        call exchpw91(lbecke,rho2,s,u,v,exdnpw91,vxdnpw91)
        call exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
      else
        exdnlsd=0.d0
        vxdnlsd=0.d0
        exdnpw91=0.d0
        vxdnpw91=0.d0
        exdnpbe=0.d0
        vxdnpbe=0.d0
      endif
   10 continue 
c construct total density and contribution to ex
      if(up.gt.eps.and.dn.gt.eps)then
        rho=up+dn
        exlsd=(exuplsd*up+exdnlsd*dn)/rho
        expw91=(exuppw91*up+exdnpw91*dn)/rho
        expbe=(exuppbe*up+exdnpbe*dn)/rho
      else
        exlsd=0.d0
        expw91=0.d0
        expbe=0.d0
      endif
      if(lcor.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Now do correlation
c zet=(up-dn)/rho
c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c uu=delgrad/(rho^2*twoksg^3)
c rholap=Laplacian
c vv=Laplacian/(rho*twoksg^2)
c ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
c ec=lsd correlation energy
c vcup=lsd up correlation potential
c vcdn=lsd down correlation potential
c h=gradient correction to correlation energy
c dvcup=gradient correction to up correlation potential
c dvcdn=gradient correction to down correlation potential
C atp added this line ..
      rho = up + dn
C mvs added check that both up and dn must be > eps
      if(rho.lt.eps .or. up.lt.eps .or. dn.lt.eps) then
        eclsd = 0
        vcuplsd = 0
        vcdnlsd = 0
        ecpw91 = 0
        vcuppw91 = 0
        vcdnpw91 = 0
        ecpbe = 0
        vcuppbe = 0
        vcdnpbe = 0
        return
      endif
      zet=(up-dn)/rho
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rho)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rho)
      uu=delgr/(rho*rho*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rho*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn,
     1                  H,DVCUP,DVCDN)
      eclsd=ec
      ecpbe=ec+h
      vcuplsd=vcup
      vcdnlsd=vcdn
      vcuppbe=vcup+dvcup
      vcdnpbe=vcdn+dvcdn
      call CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      call CORPW91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,DVCUP,DVCDN)
      ecpw91=ec+h
      vcuppw91=vcup+dvcup
      vcdnpw91=vcdn+dvcdn
      return
      end
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c       e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c       e_x[PBE]=e_x[unif]*FxPBE(s)
c       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
cdln  parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
cdln  parameter(um=0.2195149727645171d0,uk=0.4040d0,ul=um/uk)
      parameter(um=0.2195149727645171d0)
      save uk
      data uk/0d0/
      logical testwritexcfun

C.....dln: now can change uk
      if (uk .eq. 0d0) call setuk(0.8040d0)
      call getuk(uk)
      ul=um/uk
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
        ex=exunif
        vx=ex*thrd4
        return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24) 
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
cccccccccccccccccccc
      if(testwritexcfun()) then
      write(1219,"('ggazzz ',255d15.6)")exunif,FxPBE,U*FSS,s2*s*FSS,v*FS,s,u,v,vx
      endif
cccccccccccccccccccc
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,
     1                  H,DVCUP,DVCDN)
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c EU=unpolarized LSD correlation energy
c EURS=dEU/drs
c EP=fully polarized LSD correlation energy
c EPRS=dEP/drs
c ALFM=-spin stiffness, [c](3).
c ALFRSM=-dalpha/drs
c F=spin-scaling factor from [c](9).
c construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm-
     1((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPW91(lbecke,D,S,U,V,EX,VX)
C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  INPUT lbecke : use Becke exchange instead of PW91
C  INPUT D : DENSITY
C  INPUT S:  ABS(GRAD D)/(2*KF*D)
C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
C Burke had
C      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
C      parameter(a=7.7956D0,b1=0.004d0)
C ATP put a3 and b1 into data statements
      parameter(a1=0.19645D0,a2=0.27430D0,a4=100.d0)
      parameter(a=7.7956D0)
C Burke had
C      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
C which seems inconsistent with ax=.. in EXCHPBE
C ATP changed to
      parameter(ax=-0.738558766382022405884230032680836d0)
C ---
C      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      data a3/0.15084D0/ b1/0.004d0/
C for Becke exchange, set a3=b1=0
      if (lbecke .eq. 1) then
        a3 = 0d0
        b1 = 0d0
      endif
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
C  LOCAL EXCHANGE OPTION
C     EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
C  LOCAL EXCHANGE OPTION:
C     VX = FAC*THRD4
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(gam=0.5198421D0,fzz=1.709921D0)
      parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,1.00D0,RS,EU,EURS)
      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,1.00D0,RS,EP,EPRS)
      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
C  CALLED BY SUBROUTINE CORLSD
      IMPLICIT REAL*8 (A-H,O-Z)
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
     1                   DVCUP,DVCDN)
C  pw91 CORRELATION, modified by K. Burke to put all arguments 
c  as variables in calling statement, rather than in common block
c  May, 1996.
C  INPUT RS: SEITZ RADIUS
C  INPUT ZET: RELATIVE SPIN POLARIZATION
C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
      parameter(alf=0.09D0)
      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
C  LOCAL CORRELATION OPTION:
C     H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      if (zet .eq. 1) then
        GZ = ((1.D0+ZET)**THRDM)/3.D0
      else
        GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      endif
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  LOCAL CORRELATION OPTION:
C     DVCUP = 0.0D0
C     DVCDN = 0.0D0
      RETURN
      END

      subroutine setuk(uk)
C     implicit none
      double precision uk,uk_value
      common/ukstor/uk_value
      uk_value=uk
      return
      entry getuk(uk)
      uk=uk_value
      return
      end


