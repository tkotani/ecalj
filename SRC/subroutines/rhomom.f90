subroutine rhomom (sv_p_orhoat, qmom,vsum) ! Multipole moments Q_L = Q_aL^Zc+ Q_aL^v (See.Eq.(28) JPSJ034702)
  use m_struc_def
  use m_lmfinit,only: nsp,nbas,sspec=>v_sspec,jnlml,ispec
  use m_lgunit,only:stdo
  use m_hansr,only:corprm
  !i   orhoat:vector of offsets containing site density
  !i jnlml(ib): address of (ilm,ib) index. qmom(j:j+nlml-1) is for ib. Search jnlm below.
  !o   qmom  :Multipole moments, stored as a single long vector.
  !o         := integral r^l Y_L (rho1-rho2) + smooth core 
  !o   vsum  :sum over all sites ib of difference vs1_ib - vs2_ib
  !o         :vs1_ib = integral in sphere of estat potential[true rho]
  !o         :vs2_ib = integral in sphere of estat potential[sm rho]
  implicit none
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  integer:: ipr,iprint,j1,ib,is,igetss,lmxl,nr,nlml,ilm,j,lfoc
  real(8) ,allocatable :: rofi(:),rwgt(:), h_rv(:), v_rv(:)
  real(8):: qmom(*) , vsum,vs1(nbas),vs2(nbas), z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,ceh,rfoc
  real(8),parameter:: fpi  = 16d0*datan(1d0), y0 = 1d0/dsqrt(fpi),pi = 4d0*datan(1d0), srfpi = dsqrt(4d0*pi)
  ipr  = iprint()
  if (ipr >= 45) write(stdo,"(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x, 'Qc',8x,'Z')")
  do  ib = 1, nbas
     is = ispec(ib) 
     lmxl=sspec(is)%lmxl
     z  = sspec(is)%z
     a  = sspec(is)%a
     nr = sspec(is)%nr
     rmt= sspec(is)%rmt
     rg = sspec(is)%rg
     j1 = jnlml(ib)! (ilm,ib) index
     if (lmxl == -1) cycle
     call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     nlml = (lmxl+1)**2
     associate(rho1=>sv_p_orhoat(1,ib)%v,rho2=>sv_p_orhoat(2,ib)%v,rhoc=>sv_p_orhoat(3,ib)%v,vsum1=>vs1(ib),vsum2=>vs2(ib))
       call pvrhom(rmt,a,nlml,nr,nsp,rho1,rho2,rhoc, cofg,cofh,rg,ceh,rfoc,z, qmom(j1),vsum1,vsum2)
     endassociate
     if(ipr >= 45) then
        write(stdo,220) ib,1,qmom(j1),qmom(j1)/y0,qc,z
        do ilm = 2, nlml
           if (dabs(qmom(j1+ilm-1)) > 1d-6) write(stdo,220) ib,ilm,qmom(j1+ilm-1)
        enddo
     endif
  enddo
  vsum=sum(vs1-vs2)
  220 format(i13,i6,f12.6,f12.6,2f9.2)
end subroutine rhomom
subroutine pvrhom(rmt,a,nlml,nr,nsp,rho1,rho2,rhoc,cofg,cofh,rg,ceh,rfoc,z,qmomj,vsum1,vsum2)
  !Multipole moments for one site. Q_L = Q_aL^Zc+ Q_aL^v (See.Eq.(28) JPSJ034702)
  use m_hansr,only: hansmr
  use m_ll,only:ll
  !i   nlml  :L-cutoff for charge
  !i   nr    :number of radial mesh points
  !i   nsp   :number of spins
  !i   rofi  :radial mesh points
  !i   rwgt  :radial integration weights
  !i   rho1  :local true density*r**2, tabulated on a radial mesh
  !i   rho2  :local smoothed density*r**2, tabulated on a radial mesh
  !i   rhoc  :core density
  !i   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
  !i   cofh  :coefficient to Hankel part of pseudocore density (corprm)
  !i   rg    :smoothing radius for compensating gaussians
  !i   ceh   :energy of hankel function to fit core tail
  !i   rfoc  :smoothing radius for hankel head fitted to core tail
  !i   z     :nuclear charge
  !o   qmomj  :multipole moments for one site
  !r       qmomj= Q_aL^Zc+ Q_aL^v (See.Eq.(28) JPSJ034702)
  !             = Q_aL^Zc + \integral r^l (rho1-rho2)
  !Warn qmomj(ilm=1) may diffrer from this definition. Add + cofg - y0*z for true qmom. See code below.
  implicit none
  integer :: nlml,nr,nsp,i,ilm,l,m,lmxl,isp
  real(8) :: ceh,cofg,cofh,rfoc,rg,z,rmt,rofi(nr),rwgt(nr),h(nr),qmomj(nlml),rhoc(nr,nsp), &
       rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),xi0(nr),a,vsum1,vsum2,&
       ag,delq,fac,qcor1,qcor2,qnuc2,b,q1,facs,vhrho,vsum,cg,af,q2,v(nr),facc,&
       rhochs,rhocsm,rhonsm,ssum,sumg,xi(0:0),qcor1s!,gnu(nr)
  real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi,fpi = 4d0*pi
  call radmsh( rmt , a , nr , rofi )
  call radwgt( rmt , a , nr , rwgt )
  do ilm=1,nlml 
     l=ll(ilm)
     qmomj(ilm)= sum([ (sum(rwgt*rofi**l*(rho1(:,ilm,isp)-rho2(:,ilm,isp))), isp=1,nsp) ]) !Q_aL^v Eq.(28)
  enddo
  ag  = 1d0/rg
  do i = 1, nr
     call hansmr(rofi(i),ceh,1d0/rfoc,xi,0)
     xi0(i)=xi(0)
  enddo
  qcor1 = sum([(sum(rwgt*rhoc(:,isp)),isp=1,nsp)])      ! component1 True core charge inside rmax
  qcor2 = srfpi*cofh*sum(rwgt*xi0*rofi**2) + srfpi*cofg ! component2 Smooth core charge inside rmax
  qmomj(1) = qmomj(1) + (qcor1 - qcor2)/srfpi  !Add (true core q) - (sm core q) !qmomj is defined with coff is comp2
  b = rofi(nr)/(dexp(a*nr-a)-1d0)
  facs = 1d0/(3-nsp)
  q1=sum( rwgt* (facs*(srfpi*(rho1(:,1,1)+rho1(:,1,nsp)) + rhoc(:,1) + rhoc(:,nsp))))
  h=0d0 
  call poiss0(z,a,b,rofi,h,nr,0d0,   v, vhrho,vsum,1)
  vsum1 = sum(rwgt(2:nr)*rofi(2:nr)**2*(v(2:nr)-2d0*z/rofi(2:nr)))
  ! ... Smooth density, including compensating gaussians
  cg = qmomj(1) + cofg - y0*z ! this is true Q^Zc_aL in Eq.(25) for ilm=1.
  facc = fpi*(ag*ag/pi)**1.5d0
!  q2= srfpi*sum(rwgt*( facs*(rho2(:,1,1)+rho2(:,1,nsp)) + cg*facc*gnu + cofh*xi0*rofi**2))
  call poiss0(0d0,a,b,rofi,h,nr,0d0, v, vhrho,vsum,1)
  vsum2 = fpi*sum(rwgt*rofi**2*v)
end subroutine pvrhom
