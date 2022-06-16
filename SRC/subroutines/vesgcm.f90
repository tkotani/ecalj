subroutine vesgcm(nbas,ssite,sspec,cy,qmom,ng,gv,&
  kv,cv,cg1,cgsum,k1,k2,k3,smpot,f,gpot0,hpot0,qsmc,zsum,vrmt)
  use m_struc_def 
  use m_lmfinit,only:lat_alat
  use m_lattic,only: lat_vol,rv_a_opos
  use m_lgunit,only:stdo
  !- Adds contribution from compensating gaussians to smooth estat pot.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl rg rmt
  !i     Stored:    *
  !i     Passed to: corprm
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat vol
  !i     Stored:    *
  !i     Passed to: *
  !i   cy    :Normalization constants for spherical harmonics
  !i   qmom  :multipole moments of on-site the local density:
  !i         :integral r^l (rho1-rho2) + l=0 contr. from core spillout
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   cv    :work array
  !i   cg1   :work array holding compensating gaussian density
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  ! o Inputs/Outputs
  ! o  smpot :On input sm estat potential, without compensating
  ! o        :gaussians + hankels
  ! o        :On output, estat potential from local charges is added
  !o Outputs
  !o   cgsum :density from compensating gaussians
  !o   f     :contribution from compensating gaussians is added to force
  !o   gpot0 :vector of integrals of compensating gaussians * (phi[n0])
  !o         :where phi[n0] is the electrostatic potential without
  !o         :without the contribution from the local charge.
  !o         :To improve accuracy, the latter part is computed
  !o         :analytically; see ugcomp.f
  !o      = Gaussian v n0(r), where n0 is  the last term in Eq.28 in JPSJ.84.034702
  !o
  !o   hpot0 :integrals of phi0 * smooth Hankels (part of contr. to core)
  !o      = Hankel v n0(r), where n0 is the last term in Eq.28 in JPSJ.84.034702
  !o
  !o   qsmc  :smoothed core charge
  !o   zsum  :total nuclear charge
  !l Local variables
  !l   qcorg :gaussian part of the pseudocore charge; see corprm.f
  !l   qcorh :hankel part of the pseudocore charge; see corprm.f
  !r Remarks
  !r   Local charge consists of a sum of gaussians that compensate for
  !r   the difference in multipole moments of true and smooth local charge
  !r   and a contribution from the smooth core charge.
  !r     g(qmpol) + g(qcore-z) + h(ncore)
  !u Updates
  !u   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
  !u   20 Apr 01 Generates vrmt
  !u   23 Apr 00 Adapted from nfp ves_gcomp.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: k1,k2,k3,nbas,ng,kv(ng,3)
  double precision :: qsmc,zsum
  real(8):: qmom(*) , gv(ng,3) , cy(*) , f(3,nbas) , gpot0(*) , hpot0(nbas) , vrmt(nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex smpot(k1,k2,k3),cv(ng),cg1(ng),cgsum(ng)
  integer :: k0,nlmx,i,ib,ilm,iprint,is,iv0,kb,kmax,l,ll, lmxl,m,nlm,lfoc
  parameter (k0=3, nlmx=64)
  double precision :: alat,ceh,cofg,cofh,g2,pi,qc,qcorg,qcorh,qsc,rfoc, &
       rg,sum1,sum2,sum3,tpiba,vol,xx,y0,z
  double precision :: cof(nlmx),df(0:20),tau(3),v(3),ddot,gvr,rmt,fac, gvb
  double complex gkl(0:k0,nlmx),img
  call tcn('vesgcm')
  call stdfac(20,df)
  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4d0*pi)
  alat=lat_alat
  vol=lat_vol
  tpiba = 2*pi/alat
  call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)
  ! --- Accumulate FT of Gaussian + Hankel density for listed vectors ---
  !     and make integrals g*phi0, h*phi0
  call dpzero(cgsum,2*ng)
  iv0 = 0
  kmax = 0
  qsmc = 0d0
  zsum = 0d0
  do  ib = 1, nbas
     is=ssite(ib)%spec
     tau=rv_a_opos(:,ib) !ssite(ib)%pos
     lmxl=sspec(is)%lmxl
     rg=sspec(is)%rg
     if (lmxl == -1) goto 10
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     qsmc = qsmc+qc
     zsum = zsum+z
     nlm = (lmxl+1)**2
     if (nlm > nlmx) call rxi('vesgcm: increase nlmx to',nlm)
     ilm = 0
     do  l = 0, lmxl
        do m = -l,l
           ilm = ilm+1
           cof(ilm) = qmom(ilm+iv0)*4*pi/df(2*l+1)
           gpot0(ilm+iv0) = 0d0
        enddo
     enddo
     hpot0(ib) = 0d0
     cof(1) = cof(1) + 4*pi*y0*(qcorg-z)
     call dpzero(cg1,2*ng)
     do  i = 1, ng
        v(1) = gv(i,1)
        v(2) = gv(i,2)
        v(3) = gv(i,3)
        call gklft(v,rg,0d0,tau,alat,kmax,nlm,k0,cy,gkl)

        do  ilm = 1, nlm
           cg1(i) = cg1(i) + cof(ilm)*gkl(0,ilm)/vol
           gpot0(ilm+iv0) = gpot0(ilm+iv0) + dconjg(cv(i))*gkl(0,ilm)
        enddo
        call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy,gkl)
        cg1(i) = cg1(i) + cofh*gkl(0,1)/vol
        hpot0(ib) = hpot0(ib) + dconjg(cv(i))*gkl(0,1)
     enddo
     !        call dpadd(cgsum,cg1,1,2*ng,1d0)
     cgsum= cgsum+cg1
     !   ... Multiply factors into gpot0
     do  ilm = 1, nlm
        l = ll(ilm)
        gpot0(ilm+iv0) = gpot0(ilm+iv0)*4*pi/df(2*l+1)
     enddo
     !   ... Force of smooth density on the compensating gaussians
     sum1 = 0d0
     sum2 = 0d0
     sum3 = 0d0
     do  i = 1, ng
        xx = -dimag(dconjg(cv(i))*cg1(i))
        sum1 = sum1 + xx*gv(i,1)
        sum2 = sum2 + xx*gv(i,2)
        sum3 = sum3 + xx*gv(i,3)
     enddo
     sum1 = sum1*vol*tpiba
     sum2 = sum2*vol*tpiba
     sum3 = sum3*vol*tpiba
     f(1,ib) = f(1,ib)+sum1
     f(2,ib) = f(2,ib)+sum2
     f(3,ib) = f(3,ib)+sum3
     do  kb = 1, nbas
        f(1,kb) = f(1,kb) - sum1/nbas
        f(2,kb) = f(2,kb) - sum2/nbas
        f(3,kb) = f(3,kb) - sum3/nbas
     enddo
     iv0 = iv0+nlm
10   continue
  enddo
  ! --- Add 8pi/G**2 * (FT gaussian+Hankel density) into smpot ---
  if (iprint() > 40) write(stdo,300) cgsum(1),dble(cgsum(1)*vol)
300 format(/' vesgcm: smooth density G=0 term =', 2f11.6,'   Q = ',f12.6)
  !! Commented by obata (but not packed in git by obata ---fixed by t.kotani)
  do  i = 2, ng
     g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
     cv(i) = cv(i) + (8*pi)*cgsum(i)/g2
  enddo
  ! --- Electrostatic potential at rmt ---
  call dpzero(vrmt,nbas)
  img = (0d0,1d0)
  do  ib = 1, nbas
     is=ssite(ib)%spec
     tau=rv_a_opos(:,ib) !ssite(ib)%pos
     call dscal(3,alat,tau,1)
     rmt=sspec(is)%rmt
     !       Add a negligibly small amount to rmt to handle case rmt=0
     rmt = rmt+1d-32
     do  i = 2, ng
        v(1) = gv(i,1)*tpiba
        v(2) = gv(i,2)*tpiba
        v(3) = gv(i,3)*tpiba
        g2 = dsqrt(ddot(3,v,1,v,1))
        gvr = g2*rmt
        fac = sin(gvr)/gvr
        gvb = v(1)*tau(1) + v(2)*tau(2) + v(3)*tau(3)
        vrmt(ib) = vrmt(ib) + dble(cv(i)*fac*exp(img*gvb))
     enddo
  enddo
  ! --- Put cv back into smpot array ---
  call gvputf(ng,1,kv,k1,k2,k3,cv,smpot)
  call tcx('vesgcm')
end subroutine vesgcm
