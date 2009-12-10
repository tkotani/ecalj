      subroutine ftlxp(nbas,ssite,sspec,alat,ng,gv,cv,k0,nlm0,fkl)
C- Pkl expansion around all sites of a function given as Fourier series.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl rsmv kmxv
Ci     Stored:    *
Ci     Passed to: *
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cv    :FT of function, as list of G-vector coefficients
Ci   k0    :leadingd dimension of fkl
Ci   nlm0  :second dimension of fkl
Co Outputs
Co   fkl   :coefficient to P_kL expansion of function
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   31 May 00 Adapted from nfp ftlxp.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,nbas,ng,nlm0
      double precision ssite(1),sspec(1),gv(ng,3),alat
      double complex cv(ng),fkl(0:k0,nlm0,nbas)
C ... Local parameters
      integer ib,ilm,ipr,is,k,kmxv,l,lgunit,lmxl,ltop,m,nlm,ntop,og2,oh,
     .  oyl,stdo,igetss
      double precision a,fpi,pfac,pi,rsmv,tip,top,tpiba,df(0:20),
     .  fact(0:20),tau(3)
      double complex cfac
C ... Heap
      integer w(1)
      common /w/ w

      stdo = lgunit(1)
      call getpr(ipr)
      call setfac(20,fact)
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      fpi = 4*pi
      tpiba = 2d0*pi/alat

C ... Set up spherical harmonics for max l needed
      ltop = -1
      do  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        lmxl = igetss('spec lmxl',is,sspec)
        ltop = max0(ltop,lmxl)
      enddo
      ntop = (ltop+1)**2
      call defrr(oyl, ntop*ng)
      call defrr(og2, ng)
      call defcc(oh,  ng)
      call ropyln(ng,gv(1,1),gv(1,2),gv(1,3),ltop,ng,w(oyl),w(og2))

C ... Scale g2 by (2pi/a)**2
      call dpcopy(w(og2),w(og2),1,ng,tpiba*tpiba)

C --- Start loop over atoms ---
      do  ib = 1, nbas
        call upack('site spec pos',ssite,ib,is,tau,0,0)
        call upack('spec lmxl rsmv kmxv',sspec,is,lmxl,rsmv,kmxv,0)
        if (lmxl .eq. -1) goto 10

        nlm = (lmxl+1)**2
        if (nlm.gt.nlm0) call rxi('ftlcxp: nlm.gt.nlm0, need',nlm)
        if (kmxv.gt.k0)  call rxi('ftlcxp: kmxv.gt.k0, need',kmxv)
        if (ipr .ge. 60) write(stdo,200) ib,is,nlm,rsmv,lmxl,kmxv,tau
  200   format(' ib=',2i3,'  nlm=',i3,'  rsmv=',f6.3,
     .     '  lv,kv=',2i2,'  pos=',3f7.4/
     .     ' Expansion coefficients:')

        call ftlxp2(rsmv,tau,kmxv,k0,nlm,ng,gv,w(og2),w(oyl),w(oh),cv,
     .    fkl(0,1,ib))

C   ... Put in factors independent of G
        a = 1d0/rsmv
        cfac = (1d0,0d0)
        ilm = 0
        do  l = 0, lmxl
          do m = -l,l
            ilm = ilm+1
            do  k = 0, kmxv
              pfac = (4*a*a)**k * a**l * fact(k) * df(2*l+1)
              fkl(k,ilm,ib) = fkl(k,ilm,ib)*cfac*fpi/pfac
            enddo
          enddo
          cfac = cfac*dcmplx(0d0,tpiba)
        enddo

C   ... Printout
        if (ipr .ge. 60) then
          do  ilm = 1, nlm
            top = 0d0
            tip = 0d0
            do  k = 0, kmxv
              top = dmax1(top,dabs(dble(fkl(k,ilm,ib))))
              tip = dmax1(tip,dabs(dimag(fkl(k,ilm,ib))))
            enddo
            if (top .gt. 1d-8)
     .         write(stdo,400) ilm,(dble(fkl(k,ilm,ib)),k = 0,kmxv)
            if (tip .gt. 1d-8)
     .         write(stdo,401)     (dimag(fkl(k,ilm,ib)),k = 0,kmxv)
  400       format(i4,' RE',5f12.6:/(7x,5f12.6))
  401       format(4x,' IM',5f12.6:/(7x,5f12.6))
          enddo
        endif

   10   continue
      enddo

      call rlse(oyl)

      end

      subroutine ftlxp2(rsmv,tau,kmxv,k0,nlm,ng,gv,g2,yl,h,cv,fkl)
C     implicit none
C ... Passed parameters
      integer k0,kmxv,ng,nlm
      double precision rsmv,yl(ng,1),tau(3),gv(ng,3),g2(ng)
      double complex h(ng),cv(ng),fkl(0:k0,nlm)
C ... Local parameters
      integer i,ilm,k
      double precision scalp,sum1,sum2,pi,gam

      pi = 4d0*datan(1d0)
      gam = 0.25d0*rsmv**2

C ... Init h with conjg(phase) * exponential factor * coeff
      do  i = 1, ng
        scalp = -2*pi*(tau(1)*gv(i,1)+tau(2)*gv(i,2)+tau(3)*gv(i,3))
        h(i) = dcmplx(dcos(scalp),-dsin(scalp))*dexp(-gam*g2(i))*cv(i)
      enddo

C ... Loop over k and ilm, do the scalar products with cv
      do  k = 0, kmxv
        do  ilm = 1, nlm
          sum1 = 0d0
          sum2 = 0d0
          do  i = 1, ng
            sum1 = sum1 + dble(h(i))*yl(i,ilm)
            sum2 = sum2 + dimag(h(i))*yl(i,ilm)
          enddo
          fkl(k,ilm) = dcmplx(sum1,sum2)
        enddo
C   ... Multiply h to get next k-th power of (-g2)
        do  i = 1, ng
          h(i) = -h(i)*g2(i)
        enddo
      enddo

      end

