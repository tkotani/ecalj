#define F90 1
      subroutine msh21c(ib,ssite,sspec,slat,ng,gv,kv,k1,k2,k3,fn)
C- One-center expansion at one site of a function tabulated on a mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rmt lmxl a nr
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat nabc
Ci     Stored:    *
Ci     Passed to: *
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   cv    :work array
Ci   k1,k2,k3 dimensions of smrho,fn for smooth mesh density
Ci   fn    :F.T. of function for which to make one-center expansion
Co Outputs
Co   The one-center expansion of function fn is written to
Co   disk for one site. This is a debugging routine, and not optimized.
Cr Remarks
Cr   A PW exp(i.q.r) has a one-center expansion at radius r
Cr      sum_L C_L Y_L(r) where C_L = 4 pi i^l j_l(|rq|) Y_L(q)
Cr   Adapted from mshvmt.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ib,k1,k2,k3,ng,kv(ng,3)
      double precision ssite(1),sspec(1),slat(1),gv(ng,3)
      double complex fn(k1,k2,k3)
C ... Local parameters
      integer i,is,lmxx,nlmx,lmxl,nlm,ngabc(3),
     .  n1,n2,n3,m,ilm,l,ipr,nrmx,nr,ir
      parameter (lmxx=5,nrmx=1501)
      double precision alat,pi,tpiba,tau(3),rmt,r,fac,plat(3,3),a
C     double precision q(3),qlat(3,3)
      double precision gv2(ng,3),agv(ng),cgp(ng),sgp(ng)
      double precision phil(ng,0:lmxx)
      double complex vvali,fprli,cv(ng)
      double precision rofi(nrmx),rwgt(nrmx)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
#if F90
      real(8),allocatable:: yl(:,:),f1c(:,:)
#else
      parameter (nlmx=(lmxx+1)**2)
      double precision yl(ng,nlmx),f1c(nrmx,nlmx)
#endif

      call getpr(ipr)
      pi = 4d0*datan(1d0)
      call upack('lat alat plat nabc',slat,alat,plat,ngabc,0,0)
      call upack('site spec pos',ssite,ib,is,tau,0,0)
      call upack('spec rmt lmxl',sspec,is,rmt,lmxl,0,0)
      call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
      nlm = (lmxl+1)**2
      tpiba = 2*pi/alat
      call gvgetf(ng,1,kv,k1,k2,k3,fn,cv)
      if (lmxl .eq. -1) return

C --- YL(G)*G**l, agv=|G| for each g ---
      call dpcopy(gv,gv2,1,3*ng,tpiba)
      call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv)
      do  i = 1, ng
        agv(i) = sqrt(agv(i))
      enddo

      call radmsh(rmt,a,nr,rofi)
      call radwgt(rmt,a,nr,rwgt)

#if F90
      nlmx = nlm
      allocate(yl(ng,nlmx),f1c(nrmx,nlmx))
#endif
      if (nlm .gt. nlmx) call rxi('msh21c: increase nlmx to',nlm)

C --- For each radius, do ---
      call dpzero(f1c,nrmx*nlm)
      do  ir = 1, nr
        r = rofi(ir)

C   ... j_l(|r*q|)/r**l for each G and l=0..lmax
        call ropbes(agv,r**2,lmxl,cgp,sgp,phil,ng)
C   ... Phases calculated straightforwardly.
        call dscal(3,alat,tau,1)
        do  i = 1, ng
          fac = -(tau(1)*gv2(i,1)+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
          cgp(i) = dcos(fac)
          sgp(i) = dsin(fac)
        enddo

C   ... Sum_G 4*pi*(i*r)**l j_l(|r*G|)/(r*G)**l YL(G) G**l
        ilm = 0
        fprli = 4*pi
        do  l  = 0, lmxl
          do  m = -l, l
            ilm = ilm+1
            vvali = 0
            do  i = 1, ng
              vvali = vvali + (phil(i,l)*yl(i,ilm))*
     .                        (cv(i)*dcmplx(cgp(i),-sgp(i)))
            enddo
            f1c(ir,ilm) = fprli*vvali
          enddo
          fprli = fprli*(0d0,1d0)*r
        enddo
      enddo

      call prrmsh('one-center expansion',rofi,f1c,nrmx,nr,nlm)

#if F90
      deallocate(yl,f1c)
#endif

      end

