      subroutine evxcv(rho,rhosp,n,nsp,lxcf,exc,excx,excc,vxc,vxcx,vxcc)
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
     .                 vxc(n),vxcx(n),vxcc(n)
C ... Local parameters
C --- Vosko-Ceperley-Alder ---
      double precision ap,xp0,bp,cp,qp,cp1,cp2,cp3,
     .                 af,xf0,bf,cf,qf,cf1,cf2,cf3
      parameter (ap=.0621814d0,bp=3.72744d0,qp=6.1519908d0,
     .           cp=12.9352d0,cp1=1.2117833d0,cp2=1.1435257d0,
     .           cp3=-.031167608d0,xp0=-.10498d0)
      parameter (af=.0310907d0,bf=7.06042d0,qf=4.7309269d0,
     .           cf1=2.9847935d0,cf2=2.7100059d0,cf3=-.1446006d0,
     .           cf=18.0578d0,xf0=-.32500d0)
      integer i
      double precision atnf,atnp,beta,dbeta,dfs,duc,vxc1,exc0,excf,excp,
     .  fs,rs,s,s4,tf1,th,th4,thx,tp1,ucf,ucp,vxc0,x,xfx,xpx,sth
      parameter (th=1d0/3d0, th4=4d0/3d0, thx=0.519842099789746d0)

C --- Hedin-Lundqvist ---
      double precision mucf,mucp,nuce,unthrd,fothrd,fpi3,
     .  cph,cfh,aph,afh,fh,fhofxi,z,xi,cmuxp,efmep,epscf,epscp,
     .  epsx0,epsx1,gamma,tauce,rmin
      parameter(unthrd=1d0/3d0,fothrd=4d0/3d0,rmin=1d-21,
     .  fpi3=12.566370614d0/3,gamma=5.1297628d0,cmuxp=1.2217741d0,
     .  epsx0=.9163305866d0,epsx1=.23817361d0)
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
     .           -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
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
     .           -cf3*(dlog((x-xf0)**2/xfx) + cf2*atnf))
            exc0 = excp + fs*(excf-excp)*(1d0+s4*beta)
            tf1  = (x*x+bf*x)/xfx
            ucf  = excf-af/3d0*(1d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
            vxc0 = (ucf-ucp)*fs - (excf-excp)*(s-1d0)*dfs
            duc  = (ucf-ucp)*beta*s4*fs +
     .             (excf-excp)*(-rs/3d0)*dbeta*s4*fs
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
     .           -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
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
     .        .206299474d0
            vxc(i)  =  mucp + tauce*fhofxi
     .              + (nuce - cmuxp/rs)*((2d0*xi)**unthrd) - nuce
            vxcx(i) =  - cmuxp/rs*((2d0*xi)**unthrd)
            vxcc(i) =  mucp + tauce*fhofxi
     .              +  nuce*((2d0*xi)**unthrd)             - nuce
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

      integer function excsan(lxcfun,ifi)
C- Sanity checks for a specified xc switch, optional printout
C     implicit none
      integer lxcfun,ifi
      integer ltype,lgrad,k,i1,i2
      character *80 sout

      ltype = mod(mod(lxcfun,10),3)
      lgrad = lxcfun/10

      k = 0
      if (ltype .lt. 1 .or. ltype .gt. 2) k = -1
      if (lgrad .ne. 0) k = k - 10
      excsan = k
      if (ifi .le. 0) return

      if (k .lt. 0) call rxi('evxcv: unknown XC code:',lxcfun)
      if (ltype .eq. 1) write(sout,101)
      if (ltype .eq. 2) write(sout,102)
  101 format(' XC potential is Ceperly-Alder/Vosko')
  102 format(' XC potential is Barth-Hedin')
      call strip(sout,i1,i2)
      if (lgrad .eq. 0) write(sout(i2+1:),201) 'no'
      if (lgrad .gt. 0) write(sout(i2+1:),201) 'with'
  201 format(' (',a,' gradient corrections)')
      call strip(sout,i1,i2)
      write(ifi,'(a)') sout(1:i2)

      end

