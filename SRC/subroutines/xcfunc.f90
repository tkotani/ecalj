module m_xclda
  public evxcp,evxcv ,vxcgr2,vxnloc !,vxcgga
  private
contains
subroutine vxcgga(lxcg,n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom, &
     grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
  !- PW91 and PBE gradient corrections to exc and vxc for a set of points
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n     :number of points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   rhop  :spin up density if nsp=2; otherwise total density
  !i   rhom  :spin down density if nsp=2; otherwise not used
  !i   grhop :|grad rhop| or |grad rho| if nsp=1
  !i   grhom :|grad rhom| (nsp=2)
  !i   ggrhop:Laplacian of rhop
  !i   ggrhom:Laplacian of rhom
  !i   grho  :|grad total rho| if nsp=2; otherwise not used
  !i   grpgrm:grad rho+ . grad rho- (not used here)
  !i   grggr :(grad rho).(grad |grad rho|)
  !i   grggrp:(grad up).(grad |grad up|) if nsp=2; otherwise ggrgr
  !i   grggrm:(grad dn).(grad |grad dn|) if nsp=2; otherwise not used
  !o Outputs
  !o   vxc1  :GGA contr. to vxc+ if nsp=2; otherwise GGA contr. to vxc
  !o   vxc2  :GGA contr. to vxc- if nsp=2; otherwise GGA contr. to vxc
  !o    exc  :GGA contr. to exc
  !l Local variables
  !l   symbols match easypbe
  !b Bugs
  !r Remarks
  !u Updates
  !-----------------------------------------------------------------------
  implicit none
  ! Passed parameters
  integer :: nsp,n,lxcg
  double precision :: rhop(1),rhom(1),grhop(1),grhom(1),grho(1), &
       ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1), &
       vxc1(1),vxc2(1),exc(1)
  ! Local Variables
  integer :: i,isw
  double precision :: up,dn,agrup,delgrup,uplap,agrdn,delgrdn,dnlap, &
       agr,delgr,aaa,bbb
  double precision :: exlsd,vxuplsd,vxdnlsd,eclsd, &
       vcuplsd,vcdnlsd,expw91,vxuppw91,vxdnpw91, &
       ecpw91,vcuppw91,vcdnpw91,expbe,vxuppbe, &
       vxdnpbe,ecpbe,vcuppbe,vcdnpbe

  !logical:: testwritexcfun
  if (lxcg /= 3 .AND. lxcg /= 4) call rx('evxcp: bad lxcg')

  do  i = 1, n
     if (nsp == 1) then
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
     ! cccccccccccccccccccc
     !        write(1219,*)' goto easypbe i lxcg=',i,lxcg
     ! cccccccccccccccccccc
     call easypbe(up,agrup,delgrup,uplap, &
          dn,agrdn,delgrdn,dnlap,agr,delgr,1,1, &
          isw(lxcg.eq.4), &
          exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
          expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
          expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

     if (lxcg == 2) then
        exc(i) = exc(i) + (expw91 + ecpw91 - exlsd - eclsd) * 2d0
        vxc1(i) = vxc1(i) + &
             (vxuppw91 + vcuppw91 - vxuplsd - vcuplsd) * 2d0
        if (nsp == 2) vxc2(i) = vxc2(i) + &
             (vxdnpw91 + vcdnpw91 - vxdnlsd - vcdnlsd) * 2d0
     else
        ! ccccccccccccccc
        exc(i) = exc(i) + (expbe + ecpbe - exlsd - eclsd) * 2d0
        !          exc(i) = (expbe + ecpbe - exlsd - eclsd) * 2d0
        !          exc(i) = (exlsd + eclsd) * 2d0
        !          exc(i) = (expbe + ecpbe ) * 2d0
        ! ccccccccccccccc
        vxc1(i) = vxc1(i) + &
             (vxuppbe + vcuppbe - vxuplsd - vcuplsd) * 2d0

        ! ccccccccccccccccccc
        if(testwritexcfun()) then
           aaa=(vxuppbe + vcuppbe - vxuplsd - vcuplsd)
           bbb=(vxdnpbe + vcdnpbe - vxdnlsd - vcdnlsd)
           !          if(abs(aaa)>1.or.abs(bbb)>1) then
           write(1219,"('ggaxxx ',i8,255d15.6)")i,up,agrup,delgrup,uplap
           write(1219,"('ggaxxx ',i8,255d15.6)")i,dn,agrdn,delgrdn,dnlap
           write(1219,"('ggaxxx ',i8,255d15.6)")i,agr,delgr
           write(1219,"('ggaxxx ',i8,255f12.3)")i,aaa, vxuppbe, vcuppbe, vxuplsd, vcuplsd
           write(1219,"('ggaxxx ',i8,255f12.3)")i,bbb, vxdnpbe, vcdnpbe, vxdnlsd, vcdnlsd
           write(1219,*)
           !          endif
        endif
        ! ccccccccccccccccccccccc

        if (nsp == 2) vxc2(i) = vxc2(i) + &
             (vxdnpbe + vcdnpbe - vxdnlsd - vcdnlsd) * 2d0
     endif
  enddo

  ! ccccccccccccccccccccccc
  if(testwritexcfun()) then
     stop 'xxxxxxxxx end of vxcgga xxxxxxxxxxxxxxxxxxxxxxxx'
  endif
  ! ccccccccccccccccccccccc

end subroutine vxcgga


subroutine vxcgr2(nr,nsp,nrx,rofi,rp, exc,vxc)
  use m_lmfinit,only: lxcf_g=>lxcf
  !      subroutine vxcgr2(nr,nsp,nrx,rofi,rp,
  !     .grh,ggrh,agrh,grgagr,exc,vxc)
  ! akao automatic array version
  !- Gradient correction to vxc, exc for a mesh of points.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nrx   :leading dimension of the radial function arrays
  !i   rofi  :radial mesh points
  !i   rp    :density rho on a radial mesh
  !i         :the following work arrays are dimensioned (nrx,2)
  !i   grh   :work array : radial grad rho
  !i   ggrh  :work array : laplacian rho
  !i   agrh  :work array : abs(grh)
  !i   grgagr:work array : grad rho . grad abs grad rho
  !o Outputs
  !o   exc   :gradient contribution to energy added to exc
  !o   vxc   :gradient contribution to potential added to vxc
  !l Local variables
  !r Remarks
  !r
  !u Updates
  !u   18 Jun 04 Bug fix
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nr,nsp,nrx
  double precision :: rp(nrx,nsp),grh(nrx,2),ggrh(nrx,2), &
       agrh(nrx,4),grgagr(nrx,3),exc(nrx),vxc(nrx,2),rofi(nr)
  ! ... Local parameters
  integer :: ir,i,lxcf,lxcg,nglob

  ! angenglob      lxcg = mod(nglob('lxcf')/100,100)
  !      lxcg = mod(globalvariables%lxcf/100,100)


  ! ccccccccccccccccccccccc
  lxcg=3
  ! ccccccccccccccccccccccc


  !      integer iprint
  !      if (iprint() .ge. 80) then
  !        call prmr(20,rofi,rp(2,1),1)
  !        call prmr(20,rofi,rp(2,2),1)
  !      endif

  ! --- grad(rho), laplacian rho ---
  call radgrx(nr,nrx,nsp,rofi,rp,grh)
  call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
  do  20  i  = 1, nsp
     do  ir = 2, nr
        ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
     enddo
     ggrh(1,i) =(rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))

     ! --- grad rho . grad abs grad rho ---
     do   ir = 1, nr
        agrh(ir,i) = dabs(grh(ir,i))
     enddo
     call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
     do  ir = 1, nr
        grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
     enddo
20 enddo

  ! --- Extra terms g(n), g(n+).g(n-), g(n).g(abs(g(n))) if spin pol ---
  if (nsp == 2) then
     do   ir = 1, nr
        agrh(ir,3) = dabs(grh(ir,1)+grh(ir,2))
     enddo
     call radgrx(nr,nrx,1,rofi,agrh(1,3),grgagr(1,3))
     do  ir = 1, nr
        grgagr(ir,3) = (grh(ir,1)+grh(ir,2))*grgagr(ir,3)
     enddo
     do   ir = 1, nr
        agrh(ir,4) = grh(ir,1)*grh(ir,2)
     enddo
  endif

  ! --- Gradient term for all points ---
  if (lxcg >= 3) then
     ! angenglob        lxcf = mod(nglob('lxcf'),100)
     lxcf = mod(lxcf_g,100) !globalvariables%lxcf,100)
     if (lxcf /= 3 .AND. lxcf /= 4) call &
          rx('vxcgf2: inconsistent use of local and GGA functionals')
     call vxcgga(lxcg,nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp), &
          ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4), &
          grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp), &
          vxc(1,1),vxc(1,nsp),exc)
  elseif (lxcg == 2) then
     call rx('PW91 no longer implemented')
  else
     call vxnloc(nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp), &
          ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4), &
          grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp), &
          vxc(1,1),vxc(1,nsp),exc)
  endif
  do  i = 1, nsp
     vxc(1,i) = (vxc(2,i)*rofi(3)-vxc(3,i)*rofi(2))/(rofi(3)-rofi(2))
  enddo

end subroutine vxcgr2

subroutine vxnloc(n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom, &
     grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
  !- Langreth-Mehl-Hu gradient correction to exc and vxc
  ! ---------------------------------------------------------------
  !i Inputs:
  !i   rhop  :spin up density if nsp=2; otherwise total density
  !i   rhom  :spin down density if nsp=2; otherwise not used
  !i   grhop :|grad rhop| or |grad rho| if nsp=1
  !i   grhom :|grad rhom| (nsp=2)
  !i   grho  :|grad total rho| if nsp=2; otherwise not used
  !i   ggrhop:Laplacian of rhop
  !i   ggrhom:Laplacian of rhom
  !i   grggr :(grad rho).(grad |grad rho|)
  !i   grggrp:(grad up).(grad |grad up|) (not used here)
  !i   grggrm:(grad dn).(grad |grad dn|) (not used here)
  !i   grpgrm: grad rho+ . grad rho- (nsp=2)
  !o Outputs:
  !o   vxc, exc
  !r Remarks:
  !r   References PRB28,1809(1983) and PRB40,1997(1989).
  !r   If nsp=1, rhop, grhop, etc are for total rho.
  !r   factor f is empirically determined cutoff.  f=0 for pure gradient.
  !r   cutoff eliminates the blowup of vxc for small rho,
  !r   for numerical convenience.  Should be no significant change if 0.
  !r   Langreth-Mehl form does not use grggrp,grggrm
  !r   Should be used together with Barth-Hedin functional.
  ! ----------------------------------------------------------------
  implicit none
  ! Passed parameters
  integer :: nsp,n,i
  double precision :: rhop(1),rhom(1),grhop(1),grhom(1),grho(1), &
       ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1), &
       vxc1(1),vxc2(1),exc(1)
  ! Local parameters
  logical :: warned
  double precision :: pi,fivth,th4,th2,th,sth,sevni,f8,f9,rho,gp2,gm2, &
       bigf,aa,d,polar,f,cutoff,cutof1,gro,ggro,g2, &
       cutof2,hh,rp,rm,grp,grm,ggrp,ggrm,rmin,xd,expf,xx
  parameter (f=0.15d0,hh=1d-3,fivth=5d0/3d0,th4=4d0/3d0, &
       th2=2d0/3d0,th=1d0/3d0,sth=7d0/3d0,sevni=7d0/9d0, rmin=1d-15)

  warned = .false.
  pi = 4d0*datan(1d0)
  f8 = (9d0*pi)**(1d0/6d0)*f
  f9 = 5d0/6d0*2d0**th2
  aa = (pi/(8d0*(3d0*pi*pi)**th4))
  if (nsp == 1) then
     do  10  i = 1, n
        rho = rhop(i)
        if (rho > rmin) then
           gro  = grhop(i)
           ggro = ggrhop(i)

           bigf = f8*gro/(rho**(7d0/6d0))
           cutoff = dexp(-hh*(gro*gro/(rho*rho))*rho**(-th2))
           expf = 2d0*dexp(-bigf)
           exc(i) = exc(i) + aa*(gro*gro/(rho**sth))*(expf-sevni)
           vxc1(i) = vxc1(i) + 2d0*aa*rho**(-th)* &
                (sevni*(ggro/rho-th2*gro*gro/(rho*rho)) - expf*( &
                ggro*(1d0-bigf/2d0)/rho - &
                (th2+bigf*(-11d0/6d0+bigf*7d0/12d0))*gro*gro/(rho*rho) + &
                bigf*(bigf-3d0)*grggr(i)/(2*gro*rho)))*cutoff
        endif
10   enddo
  else
     do  20  i = 1, n
        rho = rhop(i)+rhom(i)
        if (rho > rmin) then
           rp   = rhop(i)
           rm   = rhom(i)
           if (rp <= 0d0) then
              if ( .NOT. warned) &
                   print *, 'vxnloc (warning) rho+ not positive but rho is'
              rp = rho/1000
              warned = .true.
           endif
           if (rm <= 0d0) then
              if ( .NOT. warned) &
                   print *, 'vxnloc (warning) rho- not positive but rho is'
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

           exc(i) = exc(i) + aa/rho*(expf*g2/(rho**th4) &
                - sevni/(2d0**th)*(gp2/(rp**th4)+gm2/(rm**th4)))
           cutof1 = dexp(-hh*(gp2/(rp*rp))*(2*rp)**(-th2))
           xd = f9*rho**(th-4d0)/(d*d)*(rp**th2-rm**th2)
           xx = (2d0-bigf)*ggro/rho - &
                (th4+bigf*(-11d0/3d0+bigf*7d0/6d0))*g2/(rho*rho) + &
                bigf*(bigf-3d0)*grggr(i)/(gro*rho)
           vxc1(i) = vxc1(i) + aa/rho**th*( &
                -sevni/(2d0**th)*(rho/rp)**th*(th4*gp2/(rp*rp)-2d0*ggrp/rp) &
                -expf*(xx - &
                xd*((1d0-bigf)*rm*g2-(2d0-bigf)*rho*(grpgrm(i)+gm2)))) &
                *cutof1
           cutof2=dexp(-hh*(gm2/(rm*rm))*(2*rm)**(-th2))
           vxc2(i) = vxc2(i) + aa/rho**th*( &
                -sevni/(2d0**th)*(rho/rm)**th*(th4*gm2/(rm*rm)-2d0*ggrm/rm) &
                -expf*(xx + &
                xd*((1d0-bigf)*rp*g2-(2d0-bigf)*rho*(gp2+grpgrm(i))))) &
                *cutof2
        endif
20   enddo
  endif
end subroutine vxnloc
  ! This file includes kinds of the xc functionals in LDA
  !===================================================================
  subroutine evxcv(rho,rhosp,n,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
    !      subroutine evxcv_new(rho,rhosp,n,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
    !- XC energy density and potential for a vector of points.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   rho   :spin-1 + spin-2 density
    !i   rhosp :spin-1 density (unused for nsp=1)
    !i   n     :number of points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lxcf  :local exchange-correlation functional index
    !i         :1= Ceperly Alder
    !i         :2= Hedin-Lundqvist
    !i         :3,4= PW91
    !o Outputs
    !o   exc   :local exchange energy density for the n points
    !o   vxc   :local exchange potential for the n points
    !o         := d (rho exc) /rho)
    !o   ex    :exchange-only part of exc
    !o   vx    :exchange-only part of vxc
    !o   ec    :correlation-only part of exc
    !o   vc    :correlation-only part of vxc
    !r Remarks
    !u Updates
    !u   21 Nov 09 criterion to evaluate exc,vxc: rho>rhomin
    !u             No longer checks individual spin parts
    !u   21 Apr 09 Calls evxcp for lxcf=3,4
    !u    8 Feb 02 vx and vc (T. Miyake)
    !u   14 Jan 02 ex and ec (T. Miyake)
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: n,nsp,lxcf
    double precision :: rho(n),rhosp(n),exc(n),ex(n),ec(n), &
         vxc(n),vx(n),vc(n)
    ! ... Local parameters
    ! --- Vosko-Ceperley-Alder ---
    double precision :: ap,xp0,bp,cp,qp,cp1,cp2,cp3, &
         af,xf0,bf,cf,qf,cf1,cf2,cf3
    double precision :: wk(n),wk2(n)
    parameter (ap=.0621814d0,bp=3.72744d0,qp=6.1519908d0, &
         cp=12.9352d0,cp1=1.2117833d0,cp2=1.1435257d0, &
         cp3=-.031167608d0,xp0=-.10498d0)
    parameter (af=.0310907d0,bf=7.06042d0,qf=4.7309269d0, &
         cf1=2.9847935d0,cf2=2.7100059d0,cf3=-.1446006d0, &
         cf=18.0578d0,xf0=-.32500d0)
    integer :: i
    double precision :: atnf,atnp,beta,dbeta,dfs,duc,vxc1,exc0,excf,excp, &
         fs,rs,s,s4,tf1,th,th4,thx,tp1,ucf,ucp,vxc0,x,xfx,xpx,sth
    parameter (th=1d0/3d0, th4=4d0/3d0, thx=0.519842099789746d0)

    ! --- Hedin-Lundqvist ---
    double precision :: mucf,mucp,nuce,unthrd,fothrd,fpi3, &
         cph,cfh,aph,afh,fh,fhofxi,z,xi,cmuxp,efmep,epscf,epscp, &
         epsx0,epsx1,gamma,tauce,rmin
    parameter(unthrd=1d0/3d0,fothrd=4d0/3d0,rmin=1d-21, &
         fpi3=12.566370614d0/3,gamma=5.1297628d0,cmuxp=1.2217741d0, &
         epsx0=.9163305866d0,epsx1=.23817361d0)
    !    .  fpi3=4.1887902d0,gamma=5.1297628d0,cmuxp=1.2217741d0,
    ! ... V. Barth, Hedin parametrization ---
    !     parameter CPH=-.0504D0,CFH=-.0254D0,APH=30d0,AFH=75D0)
    ! --- Taken from ASW ---
    parameter (cph=-.045d0,cfh=-.0225d0,aph=21d0,afh=52.9167d0)
    fh(z)     = (1d0+z*z*z)*dlog(1d0+1d0/z)+.5d0*z-z*z-unthrd

    if (lxcf > 1) goto 200
    ! --- Vosko-Ceperley-Alder, spin polarized case ---
    !      if(i==1) write(*,*)'vwn isp=2' !this is to avoid problem for ifort (probably optimization bug)
    if (nsp == 2) then
       do  i = 1, n
          if (rho(i) > rmin) then
             rs = (rho(i)*fpi3)**(-th)
             x = dsqrt(rs)
             xpx = x*x+bp*x + cp
             atnp = datan(qp/(2d0*x+bp))
             excp = ap*(dlog(x*x/xpx) + cp1*atnp &
                  -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
             tp1  = (x*x+bp*x)/xpx
             ucp  = excp - ap/3d0*(1d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))

             xfx = x*x+bf*x + cf
             !           s = 2*rhosp/rho-1
             s = min(max(2*rhosp(i)/rho(i)-1d0,-1d0),1d0)
             !           sth -> (2*rhosp/rho)^1/3 as rhosp->0
             sth = (s+1d0)**th
             s4 = s**4 - 1d0
             fs = (sth**4 + (1d0-s)**th4 - 2d0) / thx
             beta = 1d0/(1.74208d0 +x*(3.182d0 +x*(.09873d0+x*.18268d0)))
             dfs = th4*(sth - (1d0-s)**th)/thx
             dbeta = -(.27402d0*x + .09873d0 + 1.591d0/x)*beta**2
             atnf = datan(qf/(2d0*x + bf))
             excf = af*(dlog(x*x/xfx) + cf1*atnf &
                  -cf3*(dlog((x-xf0)**2/xfx) + cf2*atnf))
             exc0 = excp + fs*(excf-excp)*(1d0+s4*beta)
             tf1  = (x*x+bf*x)/xfx
             ucf  = excf-af/3d0*(1d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
             vxc0 = (ucf-ucp)*fs - (excf-excp)*(s-1d0)*dfs
             duc  = (ucf-ucp)*beta*s4*fs + &
                  (excf-excp)*(-rs/3d0)*dbeta*s4*fs
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
          if (rho(i) > rmin) then
             rs = (rho(i)*fpi3)**(-th)
             x = dsqrt(rs)
             xpx = x*x+bp*x + cp
             atnp = datan(qp/(2d0*x+bp))
             excp = ap*(dlog(x*x/xpx) + cp1*atnp &
                  -cp3*(dlog((x-xp0)**2/xpx) + cp2*atnp))
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

    ! --- Barth-Hedin ---
200 continue
    if (lxcf > 2) goto 300
    if (nsp == 2) then
       do  i = 1, n
          !         if (rho(i) .gt. rmin .and. rho(i) .ge. rhosp(i)) then
          if (rho(i) > rmin) then
             !            if (rho(i) .lt. rhosp(i)) then
             !              print *, 'hi'
             !            endif
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
             fhofxi = (xi**fothrd + (1d0-xi)**fothrd -.79370052598d0)/ &
                  .206299474d0
             vxc(i)  =  mucp + tauce*fhofxi &
                  + (nuce - cmuxp/rs)*((2d0*xi)**unthrd) - nuce
             vx(i) =  - cmuxp/rs*((2d0*xi)**unthrd)
             vc(i) =  mucp + tauce*fhofxi &
                  +  nuce*((2d0*xi)**unthrd)             - nuce
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
          if (rho(i) > rmin) then
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

    ! --- PW91, PBE ---
300 continue
    if (lxcf > 4) goto 400
    if (nsp == 1) then
       call evxcp(rho,rho,n,nsp,lxcf,ex,ec,exc,vx,wk,vc,wk,vxc,wk)
    else
       wk2 = rho - rhosp
       call evxcp(rhosp,wk2,n,nsp,lxcf,ex,ec,exc,vx,wk,vc,wk,vxc,wk)
    endif
    return

    ! --- No local functional ---
400 continue
    !      call setpr(30)
    call rxi('evxcv: no functional, lxcf =',lxcf)
  end subroutine evxcv

  !===================================================================
  subroutine evxcp(r1,r2,n,nsp,lxcf,ex,ec,exc,vxup,vxdn,vcup,vcdn, &
       vxcup,vxcdn)
    ! Local part xc energy and potential for PW91 or PBE (in LDA, PBE=PW91)
    ! ----------------------------------------------------------------------
    !i Inputs: r1, r2 charge density in a mesh of n points. If nsp=1, r1
    !i         is total charge; if nsp=2 r1 and r2 are up and
    !i         down charge.
    !i         lxcf=3,4 for PBE (AKA PW91)
    !o Outputs:
    !o   ex    :exchange-only part of exc
    !o   ec    :correlation-only part of exc
    !o   exc   :local exchange energy density for the n points
    !o   vxup  :exchange-only part of vxc, spin1
    !o   vcup  :correlation-only part of vxc, spin1
    !o   vxcup :local exchange potential for the n points, spin1
    !o   vxdn  :exchange-only part of vxc, spin2
    !o   vcdn  :correlation-only part of vxc, spin2
    !o   vxcdn :local exchange potential for the n points, spin2
    !o         If nsp=1, only spin numbers are returned
    !r Remarks:
    !r         If nsp=1 r2 points to an arbitrary address. NB evxcp calls
    !r         easypbe written by K. Burke, which uses Hartree atomic units.
    !r         This routine adapted from pbevxc in FP-LMTO and uses names
    !r         longer than six characters.
    !u Updates
    !u   21 Apr 09 Returns exchange and correlation parts separately
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed Parameters
    integer :: nsp,n,lxcf
    double precision :: r1(n),r2(n),ex(n),ec(n),exc(n)
    double precision :: vxup(n),vxdn(n),vcup(n),vcdn(n),vxcup(n),vxcdn(n)
    ! Local Variables
    integer :: i,iprint
    double precision :: up, dn, exlsd, vxuplsd, vxdnlsd, eclsd, &
         vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91, &
         ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe, &
         vxdnpbe, ecpbe, vcuppbe, vcdnpbe

    if (lxcf /= 3 .AND. lxcf /= 4) then
       if (iprint() < 10) call pshpr(10)
       call rxi('evxcp cannot handle lxcf =',lxcf)
    endif

    do  i = 1, n
       if (nsp == 1) then
          up = r1(i) / 2
          dn = up
       else
          up = r1(i)
          dn = r2(i)
       endif
       call easypbe(up,0d0,0d0,0d0,dn,0d0,0d0,0d0,0d0,0d0,1,1,0, &
            exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
            expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
            expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
       !   ... Times two to convert Ha to Ry
       ex(i) = exlsd * 2d0
       ec(i) = eclsd * 2d0
       exc(i) = (exlsd + eclsd) * 2d0
       vxup(i)  = vxuplsd * 2d0
       vcup(i)  = vcuplsd * 2d0
       vxcup(i) = (vxuplsd + vcuplsd) * 2d0
       if (nsp == 2) then
          vxdn(i)  = vxdnlsd * 2d0
          vcdn(i)  = vcdnlsd * 2d0
          vxcdn(i) = (vxdnlsd + vcdnlsd) * 2d0
       endif
    enddo
  end subroutine evxcp
  !===================================================================
  subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
       agr,delgr,lcor,lpot,lbecke, &
       exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
       expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
       expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! EASYPBE is a driver for the PBE subroutines, using simple inputs
    ! K. Burke, May 14, 1996.
    ! inputs: up=up density
    !       : agrup=|grad up|
    !       : delgrup=(grad up).(grad |grad up|)
    !       : uplap=grad^2 up=Laplacian of up
    !       : dn,agrdn,delgrdn,dnlap=corresponding down quantities
    !       : agr=|grad rho|
    !       : delgr=(grad rho).(grad |grad rho|)
    !       : lcor=flag to do correlation(=0=>don't)
    !       : lpot=flag to do potential(=0=>don't)
    !       : lbecke=flag to return Becke exchange instead of PW91 (MvS)
    ! outputs: exlsd=LSD exchange energy density, so that
    !               ExLSD=int d^3r rho(r) exlsd(r)
    !        : vxuplsd=up LSD exchange potential
    !        : vxdnlsd=down LSD exchange potential
    !        : exclsd=LSD exchange-correlation energy density
    !        : vxcuplsd=up LSD exchange-correlation potential
    !        : vxcdnlsd=down LSD exchange-correlation potential
    !        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
    !        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
    !r DLN,ATP,MvS have made a few trivial changes (see below)
    !r NB: in the absence of gradient corrections, LSD=PW91=PBE
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! needed constants:
    ! pi32=3 pi**2
    ! alpha=(9pi/4)**thrd
    implicit real*8(a-h,o-z)
    implicit integer(i-n)
    parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
    parameter(pi32=29.608813203268075856503472999628d0)
    parameter(pi=3.1415926535897932384626433832795d0)
    parameter(alpha=1.91915829267751300662482032624669d0)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! PBE exchange
    ! use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
    ! do up exchange
    ! fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3)
    ! s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
    ! u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
    ! v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
    eps = dsqrt(d1mach(3))
    rho2=2.d0*up
    if(rho2 > eps)then
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
    ! repeat for down
    rho2=2.d0*dn
    if(rho2 > eps)then
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
10  continue
    ! construct total density and contribution to ex
    if(up > eps .AND. dn > eps)then
       rho=up+dn
       exlsd=(exuplsd*up+exdnlsd*dn)/rho
       expw91=(exuppw91*up+exdnpw91*dn)/rho
       expbe=(exuppbe*up+exdnpbe*dn)/rho
    else
       exlsd=0.d0
       expw91=0.d0
       expbe=0.d0
    endif
    if(lcor == 0)return
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! Now do correlation
    ! zet=(up-dn)/rho
    ! g=phi(zeta)
    ! rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
    ! sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
    ! twoksg=2*Ks*phi
    ! t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
    ! uu=delgrad/(rho^2*twoksg^3)
    ! rholap=Laplacian
    ! vv=Laplacian/(rho*twoksg^2)
    ! ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
    ! ec=lsd correlation energy
    ! vcup=lsd up correlation potential
    ! vcdn=lsd down correlation potential
    ! h=gradient correction to correlation energy
    ! dvcup=gradient correction to up correlation potential
    ! dvcdn=gradient correction to down correlation potential
    ! atp added this line ..
    rho = up + dn
    ! mvs added check that both up and dn must be > eps
    if(rho < eps .OR. up < eps .OR. dn < eps) then
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
    call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn, &
         H,DVCUP,DVCDN)
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
  end subroutine easypbe
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
    !----------------------------------------------------------------------
    !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
    !  K Burke's modification of PW91 codes, May 14, 1996
    !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !  INPUT rho : DENSITY
    !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
    !  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
    !  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
    !   (for U,V, see PW86(24))
    !  input lgga:  (=0=>don't put in gradient corrections, just LDA)
    !  input lpot:  (=0=>don't get potential and don't need U and V)
    !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! References:
    ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
    ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
    !     {\bf 40},  3399  (1989) (E).
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! Formulas:
    !       e_x[unif]=ax*rho^(4/3)  [LDA]
    ! ax = -0.75*(3/pi)^(1/3)
    !       e_x[PBE]=e_x[unif]*FxPBE(s)
    !       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
    ! uk, ul defined after [a](13)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    implicit integer(i-n)
    parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
    parameter(pi=3.14159265358979323846264338327950d0)
    parameter(ax=-0.738558766382022405884230032680836d0)
    ! ln  parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
    ! ln  parameter(um=0.2195149727645171d0,uk=0.4040d0,ul=um/uk)
    parameter(um=0.2195149727645171d0)
!    save uk
!    data uk/0d0/
!    logical :: testwritexcfun
!
!.....dln: now can change uk
!    if (uk == 0d0) call setuk(0.8040d0)
    !    call getuk(uk)
    uk=0.8040d0
    ul=um/uk
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! construct LDA exchange energy density
    exunif = AX*rho**THRD
    if(lgga == 0)then
       ex=exunif
       vx=ex*thrd4
       return
    endif
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! construct PBE enhancement factor
    S2 = S*S
    P0=1.d0+ul*S2
    FxPBE = 1d0+uk-uk/P0
    EX = exunif*FxPBE
    if(lpot == 0)return
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !  ENERGY DONE. NOW THE POTENTIAL:
    !  find first and second derivatives of Fx w.r.t s.
    !  Fs=(1/s)*d FxPBE/ ds
    !  Fss=d Fs/ds
    Fs=2.d0*uk*ul/(P0*P0)
    Fss=-4.d0*ul*S*Fs/P0
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! calculate potential from [b](24)
    VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
    ! cccccccccccccccccc
    if(testwritexcfun()) then
       write(1219,"('ggazzz ',255d15.6)")exunif,FxPBE,U*FSS,s2*s*FSS,v*FS,s,u,v,vx
    endif
    ! cccccccccccccccccc
    RETURN
  END SUBROUTINE EXCHPBE
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn, &
       H,DVCUP,DVCDN)
    !----------------------------------------------------------------------
    !  Official PBE correlation code. K. Burke, May 14, 1996.
    !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
    !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
    !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
    !       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
    !       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
    !       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
    !       :  UU,VV,WW, only needed for PBE potential
    !       : lgga=flag to do gga (0=>LSD only)
    !       : lpot=flag to do potential (0=>energy only)
    !  output: ec=lsd correlation energy from [a]
    !        : vcup=lsd up correlation potential
    !        : vcdn=lsd dn correlation potential
    !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
    !        : dvcup=nonlocal correction to vcup
    !        : dvcdn=nonlocal correction to vcdn
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! References:
    ! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
    !     {\sl Generalized gradient approximation made simple}, sub.
    !     to Phys. Rev.Lett. May 1996.
    ! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
    !     construction of a generalized gradient approximation:  The PW91
    !     density functional}, submitted to Phys. Rev. B, Feb. 1996.
    ! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    implicit integer(i-n)
    ! thrd*=various multiples of 1/3
    ! numbers for use in LSD energy spin-interpolation formula, [c](9).
    !      GAM= 2^(4/3)-2
    !      FZZ=f''(0)= 8/(9*GAM)
    ! numbers for construction of PBE
    !      gamma=(1-log(2))/pi^2
    !      bet=coefficient in gradient expansion for correlation, [a](4).
    !      eta=small number to stop d phi/ dzeta from blowing up at
    !          |zeta|=1.
    parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
    parameter(sixthm=thrdm/2.d0)
    parameter(thrd4=4.d0*thrd)
    parameter(GAM=0.5198420997897463295344212145565d0)
    parameter(fzz=8.d0/(9.d0*GAM))
    parameter(gamma=0.03109069086965489503494086371273d0)
    parameter(bet=0.06672455060314922d0,delt=bet/gamma)
    parameter(eta=1.d-12)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! find LSD energy contributions, using [c](10) and Table I[c].
    ! EU=unpolarized LSD correlation energy
    ! EURS=dEU/drs
    ! EP=fully polarized LSD correlation energy
    ! EPRS=dEP/drs
    ! ALFM=-spin stiffness, [c](3).
    ! ALFRSM=-dalpha/drs
    ! F=spin-scaling factor from [c](9).
    ! construct ec, using [c](8)
    rtrs=dsqrt(rs)
    CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
         0.49294D0,rtrs,EU,EURS)
    CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
         0.62517D0,rtRS,EP,EPRS)
    CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
         0.49671D0,rtRS,ALFM,ALFRSM)
    ALFC = -ALFM
    Z4 = ZET**4
    F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
    EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! LSD potential from [c](A1)
    ! ECRS = dEc/drs [c](A2)
    ! ECZET=dEc/dzeta [c](A3)
    ! FZ = dF/dzeta [c](A4)
    ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
    FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
    ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
         -(1.D0-Z4)*ALFM/FZZ)
    COMM = EC -RS*ECRS/3.D0-ZET*ECZET
    VCUP = COMM + ECZET
    VCDN = COMM - ECZET
    if(lgga == 0)return
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! PBE correlation energy
    ! G=phi(zeta), given after [a](3)
    ! DELT=bet/gamma
    ! B=A of [a](8)
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
    if(lpot == 0)return
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
    G4 = G3*G
    T6 = T4*T2
    RSTHRD = RS/3.D0
    GZ=(((1.d0+zet)**2+eta)**sixthm- &
         ((1.d0-zet)**2+eta)**sixthm)/3.d0
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
  END SUBROUTINE CORPBE
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
    ! slimmed down version of GCOR used in PW91 routines, to interpolate
    ! LSD correlation energy, as given by (10) of
    ! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
    ! K. Burke, May 11, 1996.
    IMPLICIT REAL*8 (A-H,O-Z)
    Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
    Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
    Q2 = DLOG(1.D0+1.D0/Q1)
    GG = Q0*Q2
    Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
    GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
    RETURN
  END SUBROUTINE GCOR2
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE EXCHPW91(lbecke,D,S,U,V,EX,VX)
    !  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
    !  INPUT lbecke : use Becke exchange instead of PW91
    !  INPUT D : DENSITY
    !  INPUT S:  ABS(GRAD D)/(2*KF*D)
    !  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
    !  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
    !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
    IMPLICIT REAL*8 (A-H,O-Z)
    implicit integer(i-n)
    ! Burke had
    !      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
    !      parameter(a=7.7956D0,b1=0.004d0)
    ! ATP put a3 and b1 into data statements
    parameter(a1=0.19645D0,a2=0.27430D0,a4=100.d0)
    parameter(a=7.7956D0)
    ! Burke had
    !      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
    ! which seems inconsistent with ax=.. in EXCHPBE
    ! ATP changed to
    parameter(ax=-0.738558766382022405884230032680836d0)
    ! ---
    !      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
    parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
    data a3/0.15084D0/ b1/0.004d0/
    ! for Becke exchange, set a3=b1=0
    if (lbecke == 1) then
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
    !  LOCAL EXCHANGE OPTION
    !     EX = FAC
    !  ENERGY DONE. NOW THE POTENTIAL:
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
    !  LOCAL EXCHANGE OPTION:
    !     VX = FAC*THRD4
    RETURN
  END SUBROUTINE EXCHPW91
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
    !  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
    !  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
    !  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
    !     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
    !  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
    IMPLICIT REAL*8 (A-H,O-Z)
    parameter(gam=0.5198421D0,fzz=1.709921D0)
    parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
    F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
    CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
         0.49294D0,1.00D0,RS,EU,EURS)
    CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
         0.62517D0,1.00D0,RS,EP,EPRS)
    CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
         0.49671D0,1.00D0,RS,ALFM,ALFRSM)
    !  ALFM IS MINUS THE SPIN STIFFNESS ALFC
    ALFC = -ALFM
    Z4 = ZET**4
    EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
    !  ENERGY DONE. NOW THE POTENTIAL:
    ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
    FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
    ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
         -(1.D0-Z4)*ALFM/FZZ)
    COMM = EC -RS*ECRS/3.D0-ZET*ECZET
    VCUP = COMM + ECZET
    VCDN = COMM - ECZET
    RETURN
  END SUBROUTINE CORLSD
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
    !  CALLED BY SUBROUTINE CORLSD
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
  END SUBROUTINE GCOR
  !----------------------------------------------------------------------
  !######################################################################
  !----------------------------------------------------------------------
  SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H, &
       DVCUP,DVCDN)
    !  pw91 CORRELATION, modified by K. Burke to put all arguments
    !  as variables in calling statement, rather than in common block
    !  May, 1996.
    !  INPUT RS: SEITZ RADIUS
    !  INPUT ZET: RELATIVE SPIN POLARIZATION
    !  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
    !  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
    !  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
    !  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
    !  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
    !  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
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
    !  LOCAL CORRELATION OPTION:
    !     H = 0.0D0
    !  ENERGY DONE. NOW THE POTENTIAL:
    CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
    RSTHRD = RS/3.D0
    R4 = RSTHRD*CCRS/COEFF
    if (zet == 1) then
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
    !  LOCAL CORRELATION OPTION:
    !     DVCUP = 0.0D0
    !     DVCDN = 0.0D0
    RETURN
  end SUBROUTINE CORpw91

  logical function testwritexcfun()
    testwritexcfun=.false.
  end function testwritexcfun
end module m_xclda



subroutine radgrx(nr,nrx,nsp,ri,f,gf)
  implicit none
  integer :: nr,nrx,nsp,nn,i,iprint,jx
  double precision :: ri(nr),f(nrx,nsp),gf(nrx,nsp),tol,egf0
  logical :: lerr
  parameter (tol=1d-12,nn=6)
  do  10  i = 1, nsp
     !        call prmr(nr,ri,gf,1)
     call poldvm(ri(2),f(2,i),nr-1,nn,.false.,tol,lerr,gf(2,i))
     jx = 1
     call polint(ri(2),gf(2,i),nr-1,nn,ri,0d0,0,jx,gf(1,i),egf0)
     if (iprint() >= 50 .AND. dabs(egf0/gf(1,i)) > 1d-2) &
          print 345, gf(1,i), egf0/gf(1,i)*100
345  format(' radgrx (warning): expect error in gradient at origin:', &
          'f=',1pe10.3,' est err=',0pf7.1,'%')
10 enddo
end subroutine radgrx

subroutine dsred(nm,n,hr,ar)
  !- Reduction of nonorthogonal symmetric matrix to orthogonal form
  ! ----------------------------------------------------------------
  !i Inputs
  !i   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
  !i   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
  !i   n:  order of h and a
  !o Outputs
  !o   H replaced by H'' = L^-1 H (L+)^-1
  !r Remarks
  !r   Makes h'ij  = (hij  - sum_k<i lik h'kj)/lii
  !r         h''ij = (h'ij - sum_k<j h''ik (l*)jk)/ljj
  !r   This version uses vectorizable BLAS-style daxpy loops.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,nm
  double precision :: hr(nm,n),ar(nm,n)
  ! Local parameters
  integer :: i,j,k

  ! --- Make h' ---
  do  10  i = 1, n
     do  k = 1, i-1
        call daxpy(n,-ar(i,k),hr(k,1),nm,hr(i,1),nm)
     enddo
     call dscal(n,1/ar(i,i),hr(i,1),nm)
10 enddo

  ! --- Make h'' (lower triangle only) ---
  do  30  j = 1, n
     do    k = 1, j-1
        call daxpy(n-j+1,-ar(j,k),hr(j,k),1,hr(j,j),1)
     enddo
     call dscal(n-j+1,1/ar(j,j),hr(j,j),1)

     ! --- Copy lower triangle into upper ---
     do  i = j+1, n
        hr(j,i) =  hr(i,j)
     enddo
30 enddo

  !      print 337, hr,hi
  !      pause
  !  337 format(9f10.6)
end subroutine dsred

