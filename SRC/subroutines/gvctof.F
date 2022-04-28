      subroutine gvctof(iopt,alat,plat,pos,n1,n2,n3,gmax,ng)
      use m_ftox
C- Makes k-space cutoff associated with mesh orders n1,n2,n3
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt       0 Use default (smaller of iopt=1,2)
Ci              1 use Nyquist cutoff
Ci              2 use cutoff for largest sphere in BZ
Ci   alat,plat  Real-space lattice vectors
Ci   n1,n2,n3   no. divisions along the three lattice vectors
Co Outputs
Co   gmax       Energy cutoff for these n1,n2,n3
Co   ng         Number of lattice vectors
Cr Remarks
Cr   Adapted from nfp gvcutoff.f
Cu Updates
Cu   07 Feb 01 changed gmax tolerance to be consistent with gvlst2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n1,n2,n3,ng,iopt
      double precision gmax,alat,plat(3,3),pos(3)
C ... Local parameters
      integer ipr,iprint,m,k1,k2,k3,ipv,nv1,nv2,ipr0,ipr1,stdo,lgunit
      parameter (ipr0=30,ipr1=40)
      double precision pi,vol,tpiba,tol,volg,volgs,gvol,h1,
     .h2,h3,gg,h0,voln,gbot,g1,g2,g3,gmax1,vol1,qlat(3,3),
     .g(3),gs(3),plat1(3,3),qlat1(3,3),gmax2

      stdo = lgunit(1)
      ipr = iprint()
      pi = 4*datan(1d0)
      call dinv33(plat,1,qlat,vol)
      vol = dabs(alat**3*vol)

      tpiba = 2*pi/alat
      tol = 1d-8

C --- Basis vectors for real-space mesh and Q-space supercell ---
      do  10  m = 1, 3
        plat1(m,1) = plat(m,1)/n1
        plat1(m,2) = plat(m,2)/n2
        plat1(m,3) = plat(m,3)/n3
        qlat1(m,1) = qlat(m,1)*n1
        qlat1(m,2) = qlat(m,2)*n2
        qlat1(m,3) = qlat(m,3)*n3
   10 continue
C     pl1=alat*dsqrt(plat1(1,1)**2+plat1(2,1)**2+plat1(3,1)**2)
C     pl2=alat*dsqrt(plat1(1,2)**2+plat1(2,2)**2+plat1(3,2)**2)
C     pl3=alat*dsqrt(plat1(1,3)**2+plat1(2,3)**2+plat1(3,3)**2)

C --- Cutoff corresponding to recip supercell volume ---
      volg = (2*pi)**3 / vol
      volgs = n1*n2*n3* volg
      gvol = (volgs*3/(4*pi))**(1d0/3d0)

C --- Get shortest mesh pts where ki .ne. 0, gmax is pi/(longest) ---
      h1 = 1d10
      h2 = 1d10
      h3 = 1d10
      do   k1 = -5, 5
      do   k2 = -5, 5
      do   k3 = -5, 5
        do  22  m = 1, 3
          g(m) = k1*plat1(m,1)+k2*plat1(m,2)+k3*plat1(m,3)
   22   continue
        gg = alat*dsqrt(g(1)**2+g(2)**2+g(3)**2)
        if (k1 .ne. 0) h1 = dmin1(h1,gg)
        if (k2 .ne. 0) h2 = dmin1(h2,gg)
        if (k3 .ne. 0) h3 = dmin1(h3,gg)
      enddo
      enddo
      enddo
      h0 = dmax1(h1,h2,h3)
      gmax = pi/h0

      if(ipr>=ipr0) write(stdo,ftox)' GVCTOF: mesh division',n1,n2,n3,'length',h1,h2,h3
      voln = (4*pi/3)*gmax**3
      ipv = int(100*voln/volgs)
      if (ipr .ge. ipr1) write(stdo,311) gmax,voln,ipv
  311 format('   Nyquist cutoff pi/h    ',f10.3,'    (volume',
     .f10.2,',',i4,'%)')

C --- Alternative: non-overlapping spheres on recip superlattice ---
      gbot = 1d10
      do   k1 = -5, 5
      do   k2 = -5, 5
      do   k3 = -5, 5
        g1 = k1*qlat1(1,1)+k2*qlat1(1,2)+k3*qlat1(1,3) - pos(1)
        g2 = k1*qlat1(2,1)+k2*qlat1(2,2)+k3*qlat1(2,3) - pos(2)
        g3 = k1*qlat1(3,1)+k2*qlat1(3,2)+k3*qlat1(3,3) - pos(3)
        if (k1.ne.0 .or. k2.ne.0 .or. k3.ne.0) then
          gbot = dmin1(gbot,g1*g1+g2*g2+g3*g3)
        endif
      enddo
      enddo
      enddo
      
      gmax1 = 0.5d0 * tpiba*dsqrt(gbot)
      vol1 = (4*pi/3)*gmax1**3
      ipv = int(100*vol1/volgs)
      if (ipr .ge. ipr1) write(stdo,320) gmax1,vol1,ipv
  320 format('   largest sphere within BZ ',f8.3,'    (volume',
     .f10.2,',',i4,'%)')

      if (ipr .ge. ipr1) write(stdo,330) gvol,volgs,100
  330 format('   cutoff for largest vector',f8.3,'    (volume',
     .f10.2,',',i4,'%)')
      if (gmax1 .lt. gmax .and. iopt .eq. 0) gmax = gmax1
      if (iopt .eq. 2) gmax = gmax1

C --- Count g-vectors within sphere of radius gmax ---
      gmax2 = (gmax-tol)**2
      nv1 = 0
      nv2 = 0
      do   k1 = 0, n1-1
      do   k2 = 0, n2-1
      do   k3 = 0, n3-1
        g(1) = k1*qlat(1,1)+k2*qlat(1,2)+k3*qlat(1,3) - pos(1)
        g(2) = k1*qlat(2,1)+k2*qlat(2,2)+k3*qlat(2,3) - pos(2)
        g(3) = k1*qlat(3,1)+k2*qlat(3,2)+k3*qlat(3,3) - pos(3)
        call shorbz(g,gs,qlat1,plat1)
C        if (abs(gs(1)+7) .lt. tol .and. abs(gs(2)+1) .lt. tol) then
C          print *, gs
C        endif
        gg = (tpiba*tpiba)*(gs(1)**2+gs(2)**2+gs(3)**2)
        if (gg .le. gmax2) nv2 = nv2+1
        nv1 = nv1+1
      enddo
      enddo
      enddo
      if(ipr>=ipr0) write(stdo,ftox)' Reciprocal lattice: use sphere of radius',
     &     ftof(gmax),'keeping',nv2,'vectors of',nv1
      ng = nv2
      end

