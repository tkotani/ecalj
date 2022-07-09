subroutine rhomom (sv_p_orhoat, qmom,vsum)
  use m_struc_def
  use m_lmfinit,only: nsp,nbas,sspec=>v_sspec,jnlml,ispec
  use m_lgunit,only:stdo
  !- Multipole moments of valence sphere densities
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   orhoat:vector of offsets containing site density
  !i jnlml(ib): address of (ilm,ib) index. qmom(j:j+nlml-1) is for ib. Search jnlm below.
  !o Outputs
  !o   qmom  :Multipole moments, stored as a single long vector.
  !o         := integral r^l Y_L (rho1-rho2) + core spillout
  !o   vsum  :sum over all sites ib of difference vs1_ib - vs2_ib
  !o         :vs1_ib = integral in sphere of estat potential[true rho]
  !o         :vs2_ib = integral in sphere of estat potential[sm rho]
  ! ----------------------------------------------------------------------
  implicit none
  type(s_rv1) :: sv_p_orhoat(3,*)
  real(8):: qmom(1) , vsum,vs1(nbas),vs2(nbas)
  integer:: ipr,iprint,j1,ib,is,igetss,lmxl,nr,nlml,ilm,j,lfoc
  real(8) ,allocatable :: rofi(:),rwgt(:), h_rv(:), v_rv(:)
  real(8):: z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,ceh,rfoc
  real(8),parameter:: fpi  = 16d0*datan(1d0), y0 = 1d0/dsqrt(fpi)
  real(8),parameter:: pi = 4d0*datan(1d0), srfpi = dsqrt(4d0*pi)
  ipr  = iprint()
  if (ipr >= 40) write(stdo,"(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x, 'Qc',8x,'Z')")
  do  ib = 1, nbas
     is = ispec(ib) !ssite(ib)%spec
     lmxl=sspec(is)%lmxl
     z  = sspec(is)%z
     a  = sspec(is)%a
     nr = sspec(is)%nr
     rmt= sspec(is)%rmt
     rg = sspec(is)%rg
     j1 = jnlml(ib)! (ilm,ib) index
     if (lmxl == -1) cycle
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     nlml = (lmxl+1)**2
     call pvrhom(rmt,a,nlml, nr, nsp, &
          sv_p_orhoat(1,ib)%v, sv_p_orhoat(2,ib)%v, sv_p_orhoat(3,ib)%v,&
          cofg, cofh, rg, ceh, rfoc, z, qmom(j1), vs1(ib),vs2(ib))
     if (ipr >= 40) then
        write(stdo,220) ib,1,qmom(j1),qmom(j1)/y0,qc,z
        do  ilm = 2, nlml
           if (dabs(qmom(j1+ilm-1)) > 1d-6) write(stdo,220) ib,ilm,qmom(j1+ilm-1)
        enddo
     endif
  enddo
  vsum=sum(vs1-vs2)
  220 format(i13,i6,f12.6,f12.6,2f9.2)
end subroutine rhomom
subroutine pvrhom(rmt,a,nlml,nr,nsp,rho1,rho2,rhoc,cofg,cofh,rg,ceh,rfoc,z,qmom,vsum1,vsum2)
  use m_hansr,only: hansmr
  !- Multipole moments for one site
  ! ----------------------------------------------------------------------
  !i Inputs
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
  !o Outputs
  !o   qmom  :multipole moments for one site
  !r Remarks
  !r   Q_L = integral r^l (rho1-rho2) + l=0 contr. from core spillout
  !r   The core spillout term is:   qcore(rhoc)-z  - sm_qcore-sm_qnuc
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlml,nr,nsp,i,ilm,l,m,lmxl,ll,isp
  real(8) :: ceh,cofg,cofh,rfoc,rg,z,rmt
  real(8) :: rofi(nr),rwgt(nr),h(nr),qmom(nlml),rhoc(nr,nsp), &
       rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),xi0(nr),a,vsum1,vsum2
  real(8):: ag,delq,fac,qcor1,qcor2,qnuc2,b,q1,facs,vhrho,vsum,cg,af,q2,v(nr),facc,&
       rhochs,rhocsm,rhonsm,ssum,sumg,xi(0:0),qcor1s,gnu(nr)
  real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi,fpi = 4d0*pi
  call radmsh( rmt , a , nr , rofi )
  call radwgt( rmt , a , nr , rwgt )
  do ilm=1,nlml 
     l=ll(ilm)
     qmom(ilm)= sum([ (sum(rwgt*rofi**l*(rho1(:,ilm,isp)-rho2(:,ilm,isp))), isp=1,nsp) ])
  enddo
  ag  = 1d0/rg
  fac = 1d0/sum(rwgt*rofi**2*exp(-ag**2*rofi**2)) 
  gnu = rofi**2 * exp(-ag**2*rofi**2) ! ... Renormalize gaussian
  do  i = 1, nr
     call hansmr(rofi(i),ceh,1d0/rfoc,xi,0)
     xi0(i)=xi(0)
  enddo
  qcor2 = srfpi*cofh*sum(rwgt*xi0*rofi**2) + srfpi*cofg*fac*sum(rwgt*gnu) !smooth core charge inside rmax*fac
  qnuc2 = -z*fac *sum(rwgt*gnu) !smooth nuclear charge inside rmax * fac
  qcor1 = sum(rwgt*rhoc(:,1))  ! true core charge inside rmax
  qcor1s= sum(rwgt*rhoc(:,nsp))
  if (nsp == 2) qcor1 = qcor1+qcor1s
  delq = qcor1-z - qcor2-qnuc2 !(true core q-z) - (sm core q - sm z) * fac
  qmom(1) = qmom(1) + delq/srfpi !! ... l=0 includes core; some might be spilled out
  !------------------------------------------
  b = rofi(nr)/(dexp(a*nr-a)-1d0)
  facs = 1d0/(3-nsp)
  q1=sum( rwgt* (facs*(srfpi*(rho1(:,1,1)+rho1(:,1,nsp)) + rhoc(:,1) + rhoc(:,nsp))))
  call poiss0(z,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum1 = sum(rwgt*rofi**2*(v-2d0*z/rofi))
  ! ... Smooth density, including compensating gaussians
  cg = qmom(1) + cofg - y0*z
  facc = fpi*(ag*ag/pi)**1.5d0
  q2= srfpi*sum(rwgt*( facs*(rho2(:,1,1)+rho2(:,1,nsp)) + cg*facc*gnu + cofh*xi0*rofi**2))
  call poiss0(0d0,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum2 = fpi*sum(rwgt*rofi**2*v)
end subroutine pvrhom
