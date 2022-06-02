subroutine sitepot(z,rmt,rg,a,nr,cofg,cofh,ceh,rfoc, &
     lfoc,nlml,qmom, rofi,rwgt,rho1,rho2,rhoc,rhol1,rhol2,v1,v2, &
     wk,val_ves,cpn_ves,rho_exc,rho_vxc,val_vef,xcore,qloc,qlocc, &
     gpotb,focexc,focvxc,   ifigwa,v1_,v2_)
  !- Makes the potential at one site, and associated energy terms.
  !  Arrays:
  !  rho1-rho2  local atomic valence density (input)
  !  rhol1      true electron density, rho1 plus rhoc
  !  rhol2      compensated smooth density, rho2 plus gaussians
  !  v1         true potential at this site
  !  v2         smooth potential at this site
  !  gpotb      integrals of smooth estat pot times gaussians

  use m_lldata,only: ll
  implicit real*8 (a-h,p-z), integer (o)
  parameter ( kmx=20, lmx=6, nlmx=36, nrx=1501, lmc=14 )
  dimension rho1(nr,1),rho2(nr,1),rofi(nr),rwgt(nr),qmom(1), &
       rhol2(nr,1),vval(nlmx),wk(nr),v1(nr,1),v2(nr,1),df(0:20), &
       pkl(0:kmx,0:lmx),xi(0:20),cof(nlmx),cy((lmc+1)*(lmc+1)), &
       gpotb(nlml),rhoc(1),rhol1(nr,1),rhocsm(nrx),rhonsm(nrx), &
       rhochs(nrx)
  ! akao
  real(8) :: v1_(nr,nlml),v2_(nr,nlml)
  real(8),allocatable:: rhol1x(:,:)
  !$$$#ifdef COMMONLL
  !$$$      integer(4) ll(51**2)
  !$$$      common/llblock/ll
  !$$$#else
  !$$$      integer(4):: ll
  !$$$#endif

  allocate(rhol1x(nr,nlml))

  call tcn('sitepot')
  ipr=iprint()
  if (nr > nrx) call rxi('sitepot: need nrx at least',nr)
  pi=4d0*datan(1d0)
  srfpi=dsqrt(4*pi)
  y0=1d0/srfpi
  lmax=ll(nlml)
  if (nlml > nlmx) call rxi('sitepot: increase nlmx, need',nlml)
  call setdfac(20,df)
  call sylmnc(cy,lmc)
  b=rmt/(dexp(a*nr-a)-1.d0)
  ag=1d0/rg
  afoc=1d0/rfoc

  ! --- make core and nucleus pseudodensities
  fac=4*pi*(ag*ag/pi)**1.5d0
  ! ... renormalize gaussian
  sumg=0d0
  do i=2,nr
     r=rofi(i)
     gnu=fac* r*r * dexp(-ag*ag*r*r)
     sumg=sumg+rwgt(i)*gnu
  enddo
  if (dabs(sumg-1d0) > 1d-4) &
       call wgf('sitepot: large gaussian, integral=',sumg)
  sfac=1d0/sumg
  fac=fac*sfac

  sumh=0d0
  qcor2=0d0
  qnuc2=0d0
  qcor1=0d0
  do i=1,nr
     r=rofi(i)
     gnu=fac* r*r * dexp(-ag*ag*r*r)
     call hansmr(r,ceh,afoc,xi,1)
     rhonsm(i)=-z*gnu
     rhocsm(i)=srfpi*cofg*gnu+srfpi*cofh*xi(0)*r*r
     rhochs(i)=srfpi*cofh*xi(0)*r*r
     sumh=sumh+rwgt(i)*rhochs(i)
     qcor2=qcor2+rwgt(i)*rhocsm(i)
     qnuc2=qnuc2+rwgt(i)*rhonsm(i)
     qcor1=qcor1+rwgt(i)*rhoc(i)
  enddo
  samh=-y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
  if (ipr >= 20 .AND. dabs(samh) > 1d-6) &
       write(6,942) samh,sumh,samh-sumh
942 format(' charge in smH core:  total',f12.6,'   in sphere',f12.6/ &
       ' core spill-out charge is  ',f12.6)

  ! --- make full electron density in rhol1
  call dpcopy(rho1,rhol1,1,nr*nlml,1d0)
  call dpadd (rhol1,rhoc,1,nr,y0)

  ! --- assemble the full smooth compensated density rho2+gval+gcor+gnuc
  do ilm=1,nlml
     l=ll(ilm)
     cof(ilm)=qmom(ilm)*4*pi/df(2*l+1)
     fac=sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
     do i=1,nr
        r=rofi(i)
        gnu=cof(ilm)*fac* r**(l+2) * dexp(-ag*ag*r*r)
        rhol2(i,ilm)=rho2(i,ilm)+gnu
     enddo
  enddo

  do i=1,nr
     rhol2(i,1)=rhol2(i,1)+y0*(rhocsm(i)+rhonsm(i))
  enddo

  ! ... check sphere neutrality for safety
  sum1=0d0
  sum2=0d0
  do i=1,nr
     sum2=sum2+rwgt(i)*rhol2(i,1)*srfpi
     sum1=sum1+rwgt(i)*rhol1(i,1)*srfpi
  enddo
  sum1=sum1-z
  if (dabs(sum1-sum2) > 1d-6) &
       call rxf('sitepot: sphere not neutral',sum1-sum2)

  ! --- solve poisson equation for the true and smooth densities
  call dpzero(vval, nlml)
  call poinsp(z,vval,nlml,a,b,v1,rofi,rhol1,wk,nr, &
       rhoves,rhves1,vnucl,vsum)
  call poinsp(0d0,vval,nlml,a,b,v2,rofi,rhol2,wk,nr, &
       rhoves,rhves1,vnucl,vsum)

  ! ccccccccccccccccccccccccccc
  !      print *," v1 v2=",v1(nr,1),v2(nr,1)
  !      stop
  ! ccccccccccccccccccccccccccc

  ! --- gpotb: integrals of compensating gaussians times smooth estat pot
  sgpotb=0d0
  do ilm=1,nlml
     l=ll(ilm)
     cof0=4*pi/df(2*l+1)
     fac=sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
     sum1=0d0
     do i=1,nr
        r=rofi(i)
        gnu=cof0*fac* r**(l+2) * dexp(-ag*ag*r*r)
        sum1=sum1+rwgt(i)*v2(i,ilm)*gnu
     enddo
     gpotb(ilm)=sum1
     sgpotb=sgpotb+qmom(ilm)*gpotb(ilm)
  enddo

  ! --- electrostatic energy integrals involving nonspher terms
  rvs1=0d0
  rvs2=0d0
  do ilm=1,nlml
     do i=2,nr
        v1tot=v1(i,ilm)
        if (ilm == 1) v1tot=v1tot-srfpi*2*z/rofi(i)
        rvs1=rvs1+rwgt(i)*v1tot*rhol1(i,ilm)
        rvs2=rvs2+rwgt(i)*v2(i,ilm)*rhol2(i,ilm)
     enddo
  enddo

  ! --- electrostatic integrals involving spherical terms only
  vesn2=0d0
  vesc2=0d0
  vesc1=0d0
  vnucl=0d0
  do i=2,nr
     v1tot=y0*v1(i,1)-2*z/rofi(i)
     vesc1=vesc1+rwgt(i)*rhoc(i)*v1tot
     vesn2=vesn2+rwgt(i)*rhonsm(i)*y0*v2(i,1)
     vesc2=vesc2+rwgt(i)*rhocsm(i)*y0*v2(i,1)
     vnucl=vnucl+rwgt(i)*rhol1(i,1)*(1d0/rofi(i)-1d0/rmt)
  enddo

  vnucl=2*srfpi*vnucl+2*z/rmt+y0*vval(1)
  vesn1=-z*vnucl

  ! ... valence density times electrostatic potential
  vales1=rvs1-vesc1
  vales2=rvs2-vesn2-vesc2
  val_ves= vales1-vales2

  ! ... core plus nucleus times estatic potential
  vcpn1=vesc1+vesn1
  vcpn2=vesn2+vesc2
  cpn_ves=vcpn1-vcpn2

  ! --- add xc potentials to v1 and v2
  if(ifigwa==332) then !Omit Vxc for ifigwa==332 case
     v1_(1:nr,1:nlml)=v1(1:nr,1:nlml)
     v2_(1:nr,1:nlml)=v2(1:nr,1:nlml)
  endif

  if(ifigwa==333) then !Omit Vxc for ifigwa==332 case
     lxcfun=nglob('lxcfun')
     lxcf=mod(lxcfun,10)
     v1_(1:nr,1:nlml)=v1(1:nr,1:nlml)
     rhol1x(1:nr,1:nlml) = rhol1(1:nr,1:nlml)
     call dpadd (rhol1x,-rhoc(1:nr),1,nr,y0)
     call vxcnsp(rhol1x,v1_,rofi,rwgt,a,nr,nlml,rep1,rmu1,cy,lxcf)
  endif


  lxcfun=nglob('lxcfun')
  lxcf=mod(lxcfun,10)
  call vxcnsp(rhol1,v1,rofi,rwgt,a,nr,nlml,rep1,rmu1,cy,lxcf)

  if (lfoc == 0) then
     call vxcnsp(rho2, v2,rofi,rwgt,a,nr,nlml,rep2,rmu2,cy,lxcf)
     focexc=0d0
     focvxc=0d0
  else if (lfoc == 1) then
     call dpadd (rho2,rhochs,1,nr,y0)
     call vxcnsp(rho2, v2,rofi,rwgt,a,nr,nlml,rep2,rmu2,cy,lxcf)
     focexc=0d0
     focvxc=0d0
     call dpadd (rho2,rhochs,1,nr,-y0)
  else if (lfoc == 2) then
     call vxcdnsp(rho2,v2,rofi,rwgt,a,nr,nlml,rep2,rmu2,cy,lxcf, &
          rhochs,focexc,focvxc)
     if (ipr >= 40) write(6,941) sumh,focexc,focvxc
941  format(' smH core xc integrals',2f12.6)
  else
     call rxi('sitepot: cannot handle lfoc=',lfoc)
  endif

  if(ifigwa==333) then !Omit Vxc for ifigwa==332 case
     v2_(1:nr,1:nlml)=v2(1:nr,1:nlml)
  endif


  ! takao
  !  500 continue


  ! --- integrals over core times effective potential
  vefc1=0d0
  do i=2,nr
     v1tot=y0*v1(i,1)-2*z/rofi(i)
     vefc1=vefc1+rwgt(i)*rhoc(i)*v1tot
  enddo

  ! --- integrals involving the full nonspherical potential
  vefv2=0d0
  vefv1=0d0
  if (ipr >= 40) write (6,351)
  do ilm=1,nlml
     rvsm=0d0
     rvtr=0d0
     do i=2,nr
        vtr=v1(i,ilm)
        if (ilm == 1) vtr=vtr-srfpi*2*z/rofi(i)
        rvtr=rvtr+rwgt(i)*rho1(i,ilm)*vtr
        rvsm=rvsm+rwgt(i)*rho2(i,ilm)*v2(i,ilm)
     enddo
     vefv1=vefv1+rvtr
     vefv2=vefv2+rvsm
     top=dmax1(dabs(rvsm),dabs(rvtr))
     if(ipr >= 40 .AND. top >= 1d-8) write(6,350) ilm,rvtr,rvsm
350  format(i4,3x,2f15.6)
351  format(/' ilm',10x,'rho*vtr',9x,'rho*vsm')
  enddo
  ! ... smooth xc potential includes foca head; undo in integral
  vefv2=vefv2-focvxc

  rho_exc=rep1-rep2
  rho_vxc=rmu1-rmu2
  xcore=vefc1
  vefv2=vefv2+sgpotb
  val_vef=vefv1-vefv2

  ! --- charges, printout
  q1 =0d0
  q2 =0d0
  qcor1=0d0
  do i=1,nr
     q1=q1+rwgt(i)*rho1(i,1)*srfpi
     q2=q2+rwgt(i)*rho2(i,1)*srfpi
     qcor1=qcor1+rwgt(i)*rhoc(i)
  enddo

  qlocc=qcor1-qcor2
  if (ipr >= 40) then
     write (6,251)
     write (6,250) rep1,rep2,rho_exc,rmu1,rmu2,rho_vxc, &
          vefv1,vefv2,val_vef
250  format(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6/ &
          ' val*vef  ',3f15.6)
251  format(/' local terms:     true',11x,'smooth',9x,'local')
  endif
  qloc=q1-q2
  if (ipr >= 40) write (6,252) q1,q2,qloc,qcor1,qcor2,qlocc
252 format(' val chg: ',3f15.6/' core chg:',3f15.6)
  call tcx('sitepot')

  deallocate(rhol1x)
end subroutine sitepot


! --- getqval
subroutine getqval(lmxa,pl,z,kcor,lcor,qcor,qc,qv)
  implicit real*8 (a-h,p-z), integer (o)
  dimension pl(0:8)
  qc=0d0
  do l=0,lmxa
     konfig=pl(l)
     do konf=l+1,konfig-1
        deg=(2*(2*l+1))
        if (konf == kcor .AND. l == lcor) deg=deg+qcor
        qc=qc+deg
     enddo
  enddo
  qv=z-qc
end subroutine getqval





