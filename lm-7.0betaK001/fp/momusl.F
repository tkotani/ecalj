      subroutine momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp,
     .  rofi,rwgt,v0,v1,qum,vum)
C- Moments of ul*ul, ul*sl, sl*sl and their integrals with true pot.
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   rmt   :augmentation radius, in a.u.
Ci   lmxa  :augmentation l-cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   lmxl  :l-cutoff for density, potential on the radial mesh
Ci   nlml  :(lmxl+1)*(lmxl+1)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nsp   :number of spin channels
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights: int f(r) dr = sum_ir f(ir) * rwgt(ir)
Ci   v0    :spherical potential defining wave functions, without 2Z/r
Ci   v1    :true potential.  Spherical part need not be v0.
Co Outputs
Co   qum   :moments (u,s,gz)_l1 * (u,s,gz)_l2 * r**l
Co         :qum(l1,l2,l,1) = (u_l1 u_l2) * r**l
Co         :qum(l1,l2,l,2) = (u_l1 s_l2) * r**l
Co         :qum(l1,l2,l,3) = (s_l1 s_l2) * r**l
Co         :qum(l1,l2,l,4) = (u_l1 g_l2) * r**l
Co         :qum(l1,l2,l,5) = (s_l1 g_l2) * r**l
Co         :qum(l1,l2,l,6) = (g_l1 g_l2) * r**l
Co   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
Co         :decomposed by M where v1 = sum_M v1_M Y_M
Co         :vum(l1,l2,M,1) = (u_l1 v1_M u_l2)
Co         :vum(l1,l2,M,2) = (u_l1 v1_M s_l2)
Co         :vum(l1,l2,M,3) = (s_l1 v1_M s_l2)
Co         :vum(l1,l2,M,4) = (u_l1 v1_M g_l2)
Co         :vum(l1,l2,M,5) = (s_l1 v1_M g_l2)
Co         :vum(l1,l2,M,6) = (g_l1 v1_M g_l2)
Co         :Note that
Co         :vum(l2,l1,M,2) = (s_l1 v1_M u_l2)
Co         :vum(l2,l1,M,4) = (g_l1 v1_M u_l2)
Co         :vum(l2,l1,M,5) = (g_l1 v1_M s_l2)
Cl Local variables
Cl   ul    :radial wave functions; see makusp for definition
Cl   sl    :radial wave functions; see makusp for definition
Cl   ruu   :diagonal w.f. products; see makusp for definition
Cl   rus   :diagonal w.f. products; see makusp for definition
Cl   rss   :diagonal w.f. products; see makusp for definition
Cr Remarks
Cr   Diagonal integrals (l1=l2, l=0 potential) are treated specially:
Cr   small component of w.f. is explicitly taken into account.
Cr   This routine does not rely on any special properties of (u,s,gz)
Cu Updates
Cu   31 Jul 04 normalization of local orbitals now depends on type.
Cu             Altered argument list.
Cu   21 Aug 01 Added local orbitals.  Altered argument list.
Cu   13 Jun 00 spin polarized
Cu   16 May 00 Adapted from nfp momusl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmxa,lmxl,nlml,nr,nsp,n0
      parameter (n0=10)
      double precision rmt,a,z,rofi(nr),rwgt(nr),rsml(n0),ehl(n0),
     .  v0(nr,nsp),v1(nr,nlml,nsp),pnu(n0,nsp),pnz(n0,nsp),
     .  qum(0:lmxa,0:lmxa,0:lmxl,6,nsp),
     .  vum(0:lmxa,0:lmxa,nlml,6,nsp)
C ... Local parameters
      logical lpz1,lpz2
      integer l1,l2,lm,i,mlm,isp
      double precision pi,srfpi,xx,
     .  quu,qus,qss,quz,qsz,qzz,vuu,vus,vss,vv,vuz,vsz,vzz,
     .  ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp),
     .  sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp),
     .  gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)

      call tcn('momusl')
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      call dpzero(gz,nr*(lmxa+1)*nsp)
      call makusp(n0,z,nsp,rmt,lmxa,v0,a,nr,xx,xx,pnu,pnz,rsml,ehl,
     .  ul,sl,gz,ruu,rus,rss)

C --- Moments ---
      do  isp = 1, nsp
      do  l1 = 0, lmxa
        do  l2 = 0, lmxa
          lpz1 = pnz(l1+1,1) .ne. 0
          lpz2 = pnz(l2+1,1) .ne. 0
          do  lm = 0, lmxl
            quu = 0d0
            qus = 0d0
            qss = 0d0
            quz = 0d0
            qsz = 0d0
            qzz = 0d0
            if (l1.eq.l2 .and. lm.eq.0) then
              do  i = 2, nr
                xx = rwgt(i)
                quu = quu + xx * ruu(i,l1,1,isp)
                qus = qus + xx * rus(i,l1,1,isp)
                qss = qss + xx * rss(i,l1,1,isp)
              enddo
              if (lpz1) then
              do  i = 2, nr
                xx = rwgt(i)
                quz = quz + xx * ruu(i,l1,2,isp)
                qsz = qsz + xx * rus(i,l1,2,isp)
                qzz = qzz + xx * rss(i,l1,2,isp)
              enddo
              endif
            else
              do  i = 2, nr
                xx = rwgt(i) * rofi(i)**lm
                quu = quu + xx * ul(i,l1,isp) * ul(i,l2,isp)
                qus = qus + xx * ul(i,l1,isp) * sl(i,l2,isp)
                qss = qss + xx * sl(i,l1,isp) * sl(i,l2,isp)
              enddo
              if (lpz2) then
              do  i = 2, nr
                xx = rwgt(i) * rofi(i)**lm
                quz = quz + xx * ul(i,l1,isp) * gz(i,l2,isp)
                qsz = qsz + xx * sl(i,l1,isp) * gz(i,l2,isp)
                qzz = qzz + xx * gz(i,l1,isp) * gz(i,l2,isp)
              enddo
              endif
            endif
            qum(l1,l2,lm,1,isp) = quu
            qum(l1,l2,lm,2,isp) = qus
            qum(l1,l2,lm,3,isp) = qss
            qum(l1,l2,lm,4,isp) = quz
            qum(l1,l2,lm,5,isp) = qsz
            qum(l1,l2,lm,6,isp) = qzz
          enddo
        enddo
      enddo
      enddo

C     These should correspond to sabs, made by potpus
C      do  l1 = 0, lmxa
C        print 333, (qum(l1,l1,0,i,1), i=1,6)
C  333   format(6f12.6)
C      enddo
C      pause

C --- Integrals ul*ul*V_M, ul*sl*V_M, sl*sl*V_M ---
      do  isp = 1, nsp
      do  mlm = 1, nlml
        do  l1 = 0, lmxa
          do  l2 = 0, lmxa
            lpz1 = pnz(l1+1,1) .ne. 0
            lpz2 = pnz(l2+1,1) .ne. 0
            vuu = 0d0
            vus = 0d0
            vss = 0d0
            vuz = 0d0
            vsz = 0d0
            vzz = 0d0
            if (l1.eq.l2 .and. mlm.eq.1) then
              do  i = 2, nr
                vv = v1(i,mlm,isp) - srfpi*2d0*z/rofi(i)
                vuu = vuu + (rwgt(i)*vv) * ruu(i,l1,1,isp)
                vus = vus + (rwgt(i)*vv) * rus(i,l1,1,isp)
                vss = vss + (rwgt(i)*vv) * rss(i,l1,1,isp)
              enddo
              if (lpz1) then
                do  i = 2, nr
                  vv = v1(i,mlm,isp) - srfpi*2d0*z/rofi(i)
                  vuz = vuz + (rwgt(i)*vv) * ruu(i,l1,2,isp)
                  vsz = vsz + (rwgt(i)*vv) * rus(i,l1,2,isp)
                  vzz = vzz + (rwgt(i)*vv) * rss(i,l1,2,isp)
                enddo
              endif
            else
              do  i = 2, nr
                vv = v1(i,mlm,isp)
                if (mlm .eq. 1) vv = vv - srfpi*2d0*z/rofi(i)
                vuu = vuu + (rwgt(i)*vv) * ul(i,l1,isp)*ul(i,l2,isp)
                vus = vus + (rwgt(i)*vv) * ul(i,l1,isp)*sl(i,l2,isp)
                vss = vss + (rwgt(i)*vv) * sl(i,l1,isp)*sl(i,l2,isp)
              enddo
              if (lpz2) then
                do  i = 2, nr
                  vv = v1(i,mlm,isp)
                  if (mlm .eq. 1) vv = vv - srfpi*2d0*z/rofi(i)
                  vuz = vuz + (rwgt(i)*vv) * ul(i,l1,isp)*gz(i,l2,isp)
                  vsz = vsz + (rwgt(i)*vv) * sl(i,l1,isp)*gz(i,l2,isp)
                  vzz = vzz + (rwgt(i)*vv) * gz(i,l1,isp)*gz(i,l2,isp)
                enddo
              endif
            endif
            vum(l1,l2,mlm,1,isp) = vuu
            vum(l1,l2,mlm,2,isp) = vus
            vum(l1,l2,mlm,3,isp) = vss
            vum(l1,l2,mlm,4,isp) = vuz
            vum(l1,l2,mlm,5,isp) = vsz
            vum(l1,l2,mlm,6,isp) = vzz
          enddo
        enddo
      enddo
      enddo

C     These should correspond to vabs, made by potpus
C      do  l1 = 0, lmxa
C        print 333, (vum(l1,l1,1,i,1)/srfpi, i=1,6)
C  333   format(6f12.6)
C      enddo
C      pause

C      if (pnz(2,1) .ne. 0) then
C        call wrupw(1,2,nr,rofi,sl,gz,nlml,v1,rwgt)
C      endif


      call tcx('momusl')
      end
C      subroutine wrupw(lm1,lm2,nr,rofi,f1,f2,nlmv,vlm,rwgt)
CC- Brute force integration of wf(1)*V*wf(2) in sphere
C      implicit none
C      integer lm1,lm2,nr,nlmv
C      double precision f1(nr,0:1),f2(nr,0:1),rofi(nr),rwgt(nr)
C      double precision vlm(nr,nlmv)
C      integer nth,nnn,nlmx
C      parameter(nth=-62,nnn=62,nlmx=64)
C      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn)
C      double precision yl(nnn,nlmx),flm(nr,nlmx),sum
C      double precision f1p(nr,nnn),f2p(nr,nnn),vp(nr,nnn)
C      integer np,l1,l2,ll,nlm1,nlm2,ir,ip,lmax
C
C      l1 = ll(lm1)
C      l2 = ll(lm2)
C      nlm1 = (l1+1)**2
C      nlm2 = (l2+1)**2
C
C      print *, '!!'
C      vlm(:,1) = 0
C
C      call fpiint(nth,0,np,p,wp)
C      if (np .ne. nnn) stop 'oops'
C      call dmcpy(p,1,3,p2,np,1,np,3)
C      lmax = max(l1,l2,ll(nlmv))
C      call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax+1,np,yl,r2)
C
CC ... f1(l) -> f1(lm) -> f1(sphere)
C      flm = 0
C      do  ir = 1, nr
C        flm(ir,lm1) = f1(ir,l1)
C      enddo
C      call dgemm('N','T',nr,np,nlm1,1d0,flm,nr,yl,np,0d0,f1p,nr)
Cc     call wrhomt('f1p','f1(pointwise)',0,f1p,rofi,nr,nlm1,1)
C
C      sum = 0
C      do  ip = 1, np
C      do  ir = 1, nr
C        sum = sum + f1p(ir,ip)**2 *rwgt(ir)*wp(ip)
C      enddo
C      enddo
C      print 333, 'f1*f1 = ', sum
C  333 format(1x,a,f12.6)
C
CC ... f2(l) -> f2(lm) -> f2(sphere)
C      flm = 0
C      do  ir = 1, nr
C        flm(ir,lm2) = f2(ir,l2)
C      enddo
C      call dgemm('N','T',nr,np,nlm2,1d0,flm,nr,yl,np,0d0,f2p,nr)
CC     call wrhomt('f2p','f2(pointwise)',0,f2p,rofi,nr,nlm2,1)
C
CC ... v(lm) > v(sphere)
C      call dgemm('N','T',nr,np,nlmv,1d0,vlm,nr,yl,np,0d0,vp,nr)
C      call wrhomt('vp','v(pointwise)',0,vp,rofi,nr,nlm2,1)
C
C      sum = 0
C      do  ip = 1, np
C      do  ir = 1, nr
C        sum = sum + f2p(ir,ip)**2 *rwgt(ir)*wp(ip)
C      enddo
C      enddo
C      print 333, 'f2*f2 = ', sum
C
C      sum = 0
C      do  ip = 1, np
C      do  ir = 1, nr
C        sum = sum + f1p(ir,ip)*f2p(ir,ip) *rwgt(ir)*wp(ip)
C      enddo
C      enddo
C      print 333, 'f1*f2 = ', sum
C
C      sum = 0
C      do  ip = 1, np
C      do  ir = 1, nr
C        sum = sum + f1p(ir,ip)*vp(ir,ip)*f2p(ir,ip) *rwgt(ir)*wp(ip)
C      enddo
C      enddo
C      print *, 'f1*v*f2 = ', sum
C
C
C
C      end

