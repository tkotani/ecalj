Cgetarg...info...           structure ['vesgcm', 'lat', 'slat']
Cgetarg...info...           structure ['vesgcm', 'site', 'ssite']
Cgetarg...info...           structure ['vesgcm', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['vesgcm', 'm_struc_def']
      subroutine vesgcm(nbas,ssite,sspec,slat,cy,qmom,ng,gv,
     .  kv,cv,cg1,cgsum,k1,k2,k3,smpot,f,gpot0,hpot0,qsmc,zsum,vrmt)
       
       use m_struc_def  !Cgetarg

C- Adds contribution from compensating gaussians to smooth estat pot.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl rg rmt
Ci     Stored:    *
Ci     Passed to: corprm
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat vol
Ci     Stored:    *
Ci     Passed to: *
Ci   cy    :Normalization constants for spherical harmonics
Ci   qmom  :multipole moments of on-site the local density:
Ci         :integral r^l (rho1-rho2) + l=0 contr. from core spillout
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   cv    :work array
Ci   cg1   :work array holding compensating gaussian density
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Cio Inputs/Outputs
Cio  smpot :On input sm estat potential, without compensating
Cio        :gaussians + hankels
Cio        :On output, estat potential from local charges is added
Co Outputs
Co   cgsum :density from compensating gaussians
Co   f     :contribution from compensating gaussians is added to force
Co   gpot0 :vector of integrals of compensating gaussians * (phi[n0])
Co         :where phi[n0] is the electrostatic potential without
Co         :without the contribution from the local charge.
Co         :To improve accuracy, the latter part is computed
Co         :analytically; see ugcomp.f
Co   hpot0 :integrals of phi0 * smooth Hankels (part of contr. to core)
Co   qsmc  :smoothed core charge
Co   zsum  :total nuclear charge
Cl Local variables
Cl   qcorg :gaussian part of the pseudocore charge; see corprm.f
Cl   qcorh :hankel part of the pseudocore charge; see corprm.f
Cr Remarks
Cr   Local charge consists of a sum of gaussians that compensate for
Cr   the difference in multipole moments of true and smooth local charge
Cr   and a contribution from the smooth core charge.
Cr     g(qmpol) + g(qcore-z) + h(ncore)
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   20 Apr 01 Generates vrmt
Cu   23 Apr 00 Adapted from nfp ves_gcomp.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,ng,kv(ng,3)
      double precision qsmc,zsum
Cgetarg       double precision ssite(1),sspec(1),qmom(1),gv(ng,3),cy(1),slat(1),
Cgetarg      .  f(3,nbas),gpot0(1),hpot0(nbas),vrmt(nbas)
       real(8):: qmom(1) , gv(ng,3) , cy(1) , f(3,nbas) , gpot0(1) , 
     . hpot0(nbas) , vrmt(nbas) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex smpot(k1,k2,k3),cv(ng),cg1(ng),cgsum(ng)
C ... Local parameters
      integer stdo,k0,nlmx,i,ib,ilm,iprint,is,iv0,kb,kmax,l,lgunit,ll,
     .  lmxl,m,nlm,lfoc
      parameter (k0=3, nlmx=64)
      double precision alat,ceh,cofg,cofh,g2,pi,qc,qcorg,qcorh,qsc,rfoc,
     .  rg,sum1,sum2,sum3,tpiba,vol,xx,y0,z
      double precision cof(nlmx),df(0:20),tau(3),v(3),ddot,gvr,rmt,fac,
     .  gvb
      double complex gkl(0:k0,nlmx),img

      call tcn('vesgcm')

      call stdfac(20,df)
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack2('lat alat vol',slat,alat,vol)
       
       alat=slat%alat
       vol=slat%vol

      tpiba = 2*pi/alat

      call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)

C --- Accumulate FT of Gaussian + Hankel density for listed vectors ---
C     and make integrals g*phi0, h*phi0
      call dpzero(cgsum,2*ng)
      iv0 = 0
      kmax = 0
      qsmc = 0d0
      zsum = 0d0
      do  ib = 1, nbas
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,ib,is,tau,0,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack('spec lmxl rg',sspec,is,lmxl,rg,0,0)
         
         lmxl=sspec(is)%lmxl
         rg=sspec(is)%rg

        if (lmxl .eq. -1) goto 10
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        qc = qcorg+qcorh
        qsmc = qsmc+qc
        zsum = zsum+z
        nlm = (lmxl+1)**2
        if (nlm .gt. nlmx) call rxi('vesgcm: increase nlmx to',nlm)
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
        call dpadd(cgsum,cg1,1,2*ng,1d0)

C   ... Multiply factors into gpot0
        do  ilm = 1, nlm
          l = ll(ilm)
          gpot0(ilm+iv0) = gpot0(ilm+iv0)*4*pi/df(2*l+1)
        enddo

C   ... Force of smooth density on the compensating gaussians
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

C --- Add 8pi/G**2 * (FT gaussian+Hankel density) into smpot ---
      if (iprint() .gt. 40) write(stdo,300) cgsum(1),dble(cgsum(1)*vol)
  300 format(/' vesgcm: smooth density G=0 term =',
     .   2f11.6,'   Q = ',f12.6)
      cgsum(1) = (0d0,0d0)
      do  i = 2, ng
        g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
        cv(i) = cv(i) + (8*pi)*cgsum(i)/g2
      enddo

C --- Electrostatic potential at rmt ---
      call dpzero(vrmt,nbas)
      img = (0d0,1d0)
      do  ib = 1, nbas
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,ib,is,tau,0,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

        call dscal(3,alat,tau,1)
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack2('spec rmt',sspec,is,rmt)
         
         rmt=sspec(is)%rmt

C       Add a negligibly small amount to rmt to handle case rmt=0
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

C --- Put cv back into smpot array ---
      call gvputf(ng,1,kv,k1,k2,k3,cv,smpot)

      call tcx('vesgcm')
Cgetarg       end
       end subroutine vesgcm 


