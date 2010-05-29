      subroutine evxcv(rho,rhosp,n,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
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

