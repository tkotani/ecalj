module m_vesgcm
  use m_lmfinit,only: lmxl_i=>lmxl,rmt_i=>rmt,rg_i=>rg,nlmxlx
  use m_ll,only:ll
  private
  public vesgcm
  contains
    subroutine vesgcm(qmom,ng,gv,kv,cv,cg1,cgsum,smpot,f,gpot0,hpot0,qsmc,zsum,vrmt)! Adds contribution from gaussians+smHanmekels to 0th estat pot.
    use m_lmfinit,only:alat=>lat_alat,ispec,nbas,cy=>rv_a_ocy,nlmx,k0
    use m_lattic,only: vol=>lat_vol,rv_a_opos
    use m_lgunit,only:stdo
    use m_supot,only: k1=>n1,k2=>n2,k3=>n3
    use m_hansr,only:corprm
    !i   cy    :Normalization constants for spherical harmonics
    !i   qmom  :multipole moments of on-site the local density:
    !i         :integral r^l (rho1-rho2) + l=0 contr. from core spillout
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !i   kv    :indices for gather/scatter operations (gvlist.f)
    !i   cv    :work array
    !i   cg1   :work array holding compensating gaussian density
    !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
    !io  smpot :On input sm estat potential, without compensating
    !io        :gaussians + hankels
    !io        :On output, estat potential from local charges is added
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
    integer :: ng,kv(ng,3), i,ib,ilm,iprint,is,iv0,kb,kmax,l, lmxl,m,nlm,lfoc
    real(8)::qsmc,zsum,qmom(nlmxlx,nbas),gv(ng,3),f(3,nbas),gpot0(*),hpot0(nbas),vrmt(nbas),ceh,cofg,cofh,g2,qc,&
         qcorg,qcorh,qsc,rfoc, rg,sum1,sum2,sum3,tpiba,xx,z,cof(nlmx),df(0:20),tau(3),v(3),gvr,rmt,fac, gvb,fadd(3)
    complex(8):: smpot(k1,k2,k3),cv(ng),cg1(ng),cgsum(ng),gkl(0:k0,nlmx),img=(0d0,1d0)
    real(8),parameter::pi = 4d0*datan(1d0), y0 = 1d0/dsqrt(4d0*pi)
    call tcn('vesgcm')
    call stdfac(20,df)
  tpiba = 2*pi/alat
  call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)
  ! --- Accumulate FT of Gaussian + Hankel density for listed vectors and make integrals g*phi0, h*phi0
  cgsum=0d0
  iv0 = 0
  kmax = 0
  qsmc = 0d0
  zsum = 0d0
  do  ib = 1, nbas
     is=ispec(ib)
     tau=rv_a_opos(:,ib) 
     lmxl=lmxl_i(is)
     rg=rg_i(is)
     if (lmxl == -1) goto 10
     call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     qsmc = qsmc+qc
     zsum = zsum+z
     nlm = (lmxl+1)**2
     if (nlm > nlmx) call rxi('vesgcm: increase nlmx to',nlm)
     ilm = 0
     do  l = 0, lmxl
        do m = -l,l
           ilm = ilm+1
           cof(ilm) = qmom(ilm,ib)*4d0*pi/df(2*l+1)
           gpot0(ilm+iv0) = 0d0
        enddo
     enddo
     hpot0(ib) = 0d0
!     cof(1) = cof(1) + 4*pi*y0*(qcorg-z)
     cg1=0d0
     do  i = 1, ng
        v(:) = gv(i,:)
        call gklft(v,rg,0d0,tau,alat,kmax,nlm,k0,cy,gkl)

        do  ilm = 1, nlm
           cg1(i) = cg1(i) + cof(ilm)*gkl(0,ilm)/vol
           gpot0(ilm+iv0) = gpot0(ilm+iv0) + dconjg(cv(i))*gkl(0,ilm)
        enddo
        call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy,gkl)
        cg1(i) = cg1(i) + cofh*gkl(0,1)/vol
        hpot0(ib) = hpot0(ib) + dconjg(cv(i))*gkl(0,1)
     enddo
     cgsum= cgsum+cg1
     !   ... Multiply factors into gpot0
     do  ilm = 1, nlm
        l = ll(ilm)
        gpot0(ilm+iv0) = gpot0(ilm+iv0)*4d0*pi/df(2*l+1)
     enddo
     ! ... Force of smooth density on the compensating gaussians
     fadd = [ (-sum(dimag(dconjg(cv(:))*cg1(:))*gv(:,i)), i=1,3)]*vol*tpiba
     f(:,ib) = f(:,ib) + fadd
     do  kb = 1, nbas
        f(:,kb) = f(:,kb) - fadd/nbas
     enddo
     iv0 = iv0+nlm
10   continue
  enddo
  ! --- Add 8pi/G**2 * (FT gaussian+Hankel density) into smpot ---
  if (iprint() > 40) write(stdo,300) cgsum(1),dble(cgsum(1)*vol)
300 format(/' vesgcm: smooth density G=0 term =', 2f11.6,'   Q = ',f12.6)
    do  i = 2, ng
       g2 = tpiba*tpiba*sum(gv(i,:)**2) 
       cv(i) = cv(i) + (8*pi)*cgsum(i)/g2
    enddo
    ! --- Electrostatic potential at rmt ---
    vrmt=0d0
    do ib = 1, nbas
       is=ispec(ib)
       tau= alat * rv_a_opos(:,ib) 
       rmt= rmt_i(is) +1d-32 !  Add a negligibly small amount to rmt to handle case rmt=0
       do i = 2, ng
          v(:) = gv(i,:)*tpiba
          g2 = (sum(v**2))**.5 
          gvr = g2*rmt
          fac = sin(gvr)/gvr
          gvb = sum(v*tau) 
          vrmt(ib) = vrmt(ib) + dble(cv(i)*fac*exp(img*gvb))
       enddo
    enddo
    call gvputf(ng,1,kv,k1,k2,k3,cv,smpot)!Put cv back into smpot array 
    call tcx('vesgcm')
  end subroutine vesgcm
end module m_vesgcm
