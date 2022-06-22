subroutine rhomom (sv_p_orhoat, qmom,vsum)
  use m_struc_def
  use m_lmfinit,only: nsp,nbas,ssite=>v_ssite,sspec=>v_sspec
  use m_lgunit,only:stdo
  !- Multipole moments of valence sphere densities
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i     Elts read: spec
  !i   sspec :struct containing species-specific information
  !i     Elts read: lmxl z qc a nr rmt rg
  !i     Passed to: corprm
  !i   orhoat:vector of offsets containing site density
  !o Outputs
  !o   qmom  :Multipole moments, stored as a single long vector.
  !o         := integral r^l Y_L (rho1-rho2) + core spillout
  !o   vsum  :sum over all sites ib of difference vs1_ib - vs2_ib
  !o         :vs1_ib = integral in sphere of estat potential[true rho]
  !o         :vs2_ib = integral in sphere of estat potential[sm rho]
  !r Remarks
  !u Updates
  !u   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
  !u   11 Jun 00 spin polarized
  !u   22 Apr 00 Adapted from nfp rhomom.f
  ! ----------------------------------------------------------------------
  implicit none
  type(s_rv1) :: sv_p_orhoat(3,*)
  real(8):: qmom(1) , vsum
  integer:: ipr,iprint,j1,ib,is,igetss,lmxl,nr,nlml,ilm,j,lfoc
  real(8) ,allocatable :: rofi(:),rwgt(:), h_rv(:), v_rv(:)
  real(8):: z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,ceh,rfoc,vs1,vs2
  real(8),parameter:: fpi  = 16d0*datan(1d0), y0 = 1d0/dsqrt(fpi)
  real(8),parameter:: pi = 4d0*datan(1d0), srfpi = dsqrt(4d0*pi)
  ipr  = iprint()
  if (ipr >= 40) write(stdo,221)
  j1 = 1
  vsum = 0d0
  do  ib = 1, nbas
     is = ssite(ib)%spec
     lmxl=sspec(is)%lmxl
     z=sspec(is)%z
     qc=sspec(is)%qc
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     rg=sspec(is)%rg
     if (lmxl == -1) cycle
     allocate(rofi(nr),rwgt(nr))
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     nlml = (lmxl+1)**2
     call radmsh( rmt , a , nr , rofi )
     call radwgt( rmt , a , nr , rwgt )
     call pvrhom(a,nlml, nr, nsp, rofi, rwgt,&
          sv_p_orhoat(1,ib)%v, sv_p_orhoat(2,ib)%v, sv_p_orhoat(3,ib)%v,&
          cofg, cofh, rg, ceh, rfoc, z, qmom(j1), vs1,vs2)
     if (ipr >= 40) then
        write(stdo,220) ib,1,qmom(j1),qmom(j1)/y0,qc,z
        do  ilm = 2, nlml
           j = j1+ilm-1
           if (dabs(qmom(j)) > 1d-6) write(stdo,220) ib,ilm,qmom(j)
        enddo
     endif
!     call pvrhm2(z,a,nr, nlml , nsp , rofi , rwgt , sv_p_orhoat( 1 , ib )%v &
!          , sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v , qmom ( &
!          j1 ) , cofg , cofh , rg , ceh , rfoc, vs1 , vs2   )
     vsum = vsum+vs1-vs2
     j1 = j1+nlml
     deallocate(rofi,rwgt)
  enddo
220  format(i13,i6,f12.6,f12.6,2f9.2)
221  format(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x, 'Qc',8x,'Z')
end subroutine rhomom

subroutine pvrhom(a,nlml,nr,nsp,rofi,rwgt,rho1,rho2,rhoc,cofg,cofh,rg,ceh,rfoc,z,qmom,vsum1,vsum2)
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
  !i   h     :work array of dimension nr
  !o Outputs
  !o   qmom  :multipole moments for one site
  !r Remarks
  !r   Q_L = integral r^l (rho1-rho2) + l=0 contr. from core spillout
  !r   The core spillout term is:
  !r      qcore(rhoc)-z  - sm_qcore-sm_qnuc
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlml,nr,nsp
  real(8) :: ceh,cofg,cofh,rfoc,rg,z
  real(8) :: rofi(nr),rwgt(nr),h(nr),qmom(nlml),rhoc(nr,nsp), &
       rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),xi0(nr),a,vsum1,vsum2
  integer :: n0,i,ilm,l,m,lmxl,ll,isp
  parameter (n0=10)
  real(8):: ag,delq,fac,qcor1,qcor2,qnuc2,r, &
       rhochs,rhocsm,rhonsm,ssum,sumg,xi(0:n0),qcor1s,gnu(nr)
  real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi,fpi = 4d0*pi
  ! ... l=0 includes core; some might be spilled out
  !     qcor1 = true core charge inside rmax
  !     qcor2 = smooth core charge inside rmax * fac
  !     qnuc2 = smooth nuclear charge inside rmax * fac
  !     delq = (true core q-z) - (sm core q - sm z) * fac
  do ilm=1,nlml 
     l=ll(ilm)
     qmom(ilm)= sum([ (sum(rwgt*rofi**l*(rho1(:,ilm,isp)-rho2(:,ilm,isp))), isp=1,nsp) ])
  enddo
  ag  = 1d0/rg
  fac = 1d0/sum(rwgt*rofi**2*exp(-ag**2*rofi**2)) 
  gnu = fac * rofi**2 * dexp(-ag**2*rofi**2) ! ... Renormalize gaussian
  xi0(1)=0d0
  do  i = 2, nr
     call hansmr(rofi(i),ceh,1d0/rfoc,xi,1)
     xi0(i)=xi(0)
  enddo
  qcor2 = srfpi*cofh*sum(rwgt*xi0*rofi**2) + srfpi*cofg*sum(rwgt*gnu)
  qnuc2 = -z*sum(rwgt*gnu)
  qcor1 = sum(rwgt*rhoc(:,1))
  qcor1s= sum(rwgt*rhoc(:,nsp))
  if (nsp == 2) qcor1 = qcor1+qcor1s
  delq = qcor1-z - qcor2-qnuc2
  qmom(1) = qmom(1) + delq/srfpi
  !------------------------------------------
  block
    real(8):: b,q1,facs,vhrho,vsum,r,cg,ag,af,fac,q2,v(nr)
  b = rofi(nr)/(dexp(a*nr-a)-1d0)
  q1 = 0d0
  facs = 1d0/(3-nsp)
  do  i = 1, nr
     h(i) = facs*(srfpi*(rho1(i,1,1)+rho1(i,1,nsp)) + rhoc(i,1) + rhoc(i,nsp))
     q1 = q1 + rwgt(i)*h(i)
  enddo
  call poiss0(z,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum1 = 0d0
  do  i = 2, nr
     r = rofi(i)
     vsum1 = vsum1 + 4*pi*rwgt(i)*r*r*(v(i)-2*z/r)
  enddo
  ! ... Smooth density, including compensating gaussians
  cg = qmom(1) + cofg - y0*z
  ag = 1d0/rg
  af = 1d0/rfoc
  fac = fpi*(ag*ag/pi)**1.5d0
  q2 = 0d0
  facs = 1d0/(3-nsp)
  do  i = 1, nr
     r = rofi(i)
     !gnu =  r*r * dexp(-ag*ag*r*r)
     call hansmr(r,ceh,af,xi,1)
     h(i) = srfpi*(facs*(rho2(i,1,1)+rho2(i,1,nsp)) &
          + cg*fac*gnu(i) + cofh*xi(0)*r*r)
     q2 = q2 + rwgt(i)*h(i)
  enddo
  call poiss0(0d0,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum2 = 0d0
  do  i = 2, nr
     r = rofi(i)
     vsum2 = vsum2 + 4*pi*rwgt(i)*r*r*v(i)
  enddo
  endblock
end subroutine pvrhom

subroutine pvrhm2(z,a,nr,nlml,nsp,rofi,rwgt,rho1,rho2,rhoc,qmom, &
     cofg,cofh,rg,ceh,rfoc,vsum1,vsum2)
  !- Integral over electrostatic potential, to fix energy origin
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   z     :nuclear charge
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   nlml  :L-cutoff for charge
  !i   nsp   :number of spins
  !i   rofi  :radial mesh points
  !i   rwgt  :radial integration weights
  !i   rho1  :local true density, tabulated on a radial mesh
  !i   rho2  :local smoothed density, tabulated on a radial mesh
  !i   rhoc  :core density
  !i   qmom  :multipole moments for one site (only l=0 term used)
  !i   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
  !i         := y0 * pseudocore charge
  !i   cofh  :coefficient to Hankel part of pseudocore density (corprm)
  !i   rg    :smoothing radius for compensating gaussians
  !i   ceh   :energy of hankel function to fit core tail
  !i   rfoc  :smoothing radius for hankel head fitted to core tail
  !i   h     :work array of dimension nr
  !i   v     :work array (holds spherical estat potential - poiss0.f)
  !o Outputs
  !o   vsum1 :integral in sphere of electrostatic potential for true rho
  !o   vsum2 :integral in sphere of electrostatic potential for smooth rho
  !o         :including compensating gaussians
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nr,nlml,nsp
  double precision :: a,ceh,cofg,cofh,rfoc,rg,vsum1,vsum2,z,rofi(nr), &
       rwgt(nr),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),h(nr),v(nr),qmom(*), &
       rhoc(nr,nsp)
  integer :: n0,i
  parameter (n0=10)
  double precision :: af,ag,b,cg,fac,facs,fpi,gnu,pi,q1,q2,r,srfpi, vhrho,vsum,y0,xi(0:n0)
  pi = 4d0*datan(1d0)
  fpi = 4d0*pi
  srfpi = dsqrt(fpi)
  y0 = 1d0/srfpi
  b = rofi(nr)/(dexp(a*nr-a)-1d0)

  ! ... True density
  q1 = 0d0
  facs = 1d0/(3-nsp)
  do  i = 1, nr
     h(i) = facs*(srfpi*(rho1(i,1,1)+rho1(i,1,nsp)) + rhoc(i,1) + rhoc(i,nsp))
     q1 = q1 + rwgt(i)*h(i)
  enddo
  call poiss0(z,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum1 = 0d0
  do  i = 2, nr
     r = rofi(i)
     vsum1 = vsum1 + 4*pi*rwgt(i)*r*r*(v(i)-2*z/r)
  enddo

  ! ... Smooth density, including compensating gaussians
  cg = qmom(1) + cofg - y0*z
  ag = 1d0/rg
  af = 1d0/rfoc
  fac = fpi*(ag*ag/pi)**1.5d0
  q2 = 0d0
  facs = 1d0/(3-nsp)
  do  i = 1, nr
     r = rofi(i)
     gnu = fac* r*r * dexp(-ag*ag*r*r)
     call hansmr(r,ceh,af,xi,1)
     h(i) = srfpi*(facs*(rho2(i,1,1)+rho2(i,1,nsp)) &
          + cg*gnu + cofh*xi(0)*r*r)
     q2 = q2 + rwgt(i)*h(i)
  enddo
  call poiss0(0d0,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
  vsum2 = 0d0
  do  i = 2, nr
     r = rofi(i)
     vsum2 = vsum2 + 4*pi*rwgt(i)*r*r*v(i)
  enddo

  !      write (stdo,400) q1,q2,vsum1,vsum2
  !  400 format(' qsph=',2f12.6,'   vsum=',2f12.6)
end subroutine pvrhm2


