      subroutine lscoul(ri,nr,nsp,nlml,z,e,rhoL,vL)
C- Thomas-Fermi screened coulomb potential inside a sphere
C ----------------------------------------------------------------
Ci Inputs
Ci   ri    table of points on radial mesh
Ci   z     nuclear charge, inducing -2*z/r*exp(-sqrt(-e)r)
Ci   e     TF screening energy.  Input e should be <0
Ci   rhoL  density, expressed as sum rhoL(r,L) Y_L
Co Outputs
Co   vL    screened potential V = sum vL(r,L) Y_L
Cr Remarks
Cr   lscoul solves the equation for vL
Cr     (-nabla - e) vL = 8*pi*rhoL
Cr   vL is obtained by direct integration, viz
Cr   vL(ir) = 8*pi*(f1(ir) + f2(ir))
Cr       f1 = J(r) int_0^r H(r') n(r')
Cr       f2 = H(r) int_r^R J(r') n(r')
Cr   Solution corresponds to the boundary condition rhoL=0 for r>rmax.
Cr   The total charge corresponding to the screened potential is
Cr     rho-scr  =  1/8pi (-nabla) vL  =  e/8pi vL + rhoL
C ----------------------------------------------------------------
C     implicit none
      integer nr,nsp,nlml
      double precision z,e,ri(nr),rhoL(nr,nlml,nsp),vL(nr,nlml,nsp)
C Local variables
      integer nrmx,npol,ir,isp,lmxx,lmax,ll,l,ilm,iprint
      parameter (nrmx=1501, lmxx=8)
      double precision q,qeff,qscr,atepi,y0,xx,phii(0:lmxx+2),
     .  psii(0:lmxx+2),errmx,bessel(nrmx,0:lmxx),hankel(nrmx,0:lmxx),
     .  bessw(nrmx),hankw(nrmx),bessi(nrmx),hanki(nrmx)
C ... for debugging
C      integer i,nd,np
C      parameter (nd=36,np=62,nth=np)
C      double precision yl(nd,np),gyl(nrmx,nd,3),wk(nr),wk2(nr)
C      double precision grp(nr,np,3),ggrp(nr,np),ggvl(nr,nd)
C      double precision p(3,np),wp(np),xp(np),yp(np),zp(np),r2(np)

      if (nsp .ne. 1) call rx('lscoul: not spin pol')
      if (nr .gt. nrmx) call rxi('lscoul: need nrmx at least',nr)
      isp = 1
      atepi = 32*datan(1d0)
      y0 = 1/dsqrt(atepi/2)
      lmax = ll(nlml)

C --- Tabulate Bessel, Hankel functions on the radial mesh ---
C     phii, psi are Bessel / r**l,  Hankel * r**l+1
      do  21  l = 0, lmax
      bessel(1,l) = 0
      hankel(1,l) = 0
   21 continue
      bessel(1,0) = 1
#if OKA
      bessel(1,0) = 2
#endif
      do  20  ir = 2, nr
        call bessl(e*ri(ir)**2,lmax,phii,psii)
        do  20  l = 0, lmax
#if OKA
          call dscal(lmax+1,2d0,phii,1)
#endif
          xx = ri(ir)**l
          bessel(ir,l) = phii(l)*xx
          hankel(ir,l) = psii(l)/(xx*ri(ir))
   20 continue
      
C --- For each L make integral ---
      do  22  ilm = 1, nlml
        l = ll(ilm)
        bessw(1) = 0
        hankw(1) = 0
        do  24  ir = 2, nr
          bessw(ir) = bessel(ir,l)*rhoL(ir,ilm,isp)
          hankw(ir) = hankel(ir,l)*rhoL(ir,ilm,isp)
   24   continue

C   ... Integrate H_L(r) n_L(r), J_L(R) n_L(r)
        npol = 6
        call politg(ri,bessw,nr,npol,0,errmx,bessi)
        call politg(ri,hankw,nr,npol,0,errmx,hanki)
        if (ilm .eq. 1) qeff = bessi(nr)
        do  26  ir = 2, nr
          vL(ir,ilm,isp) = atepi*(hankel(ir,l)*bessi(ir) +
     .                            bessel(ir,l)*(hanki(nr)-hanki(ir)))
   26   continue
   22 continue
C ... vL(1,1)
      call polinx(ri(2),vL(2,1,isp),npol,0d0,0d0,vL(1,1,isp),errmx)

C ... Screened potential from nucleus
      do  30  ir = 2, nr
   30 vL(ir,1,isp) = vL(ir,1,isp) - 2*z/y0*hankel(ir,0)

C ... Printout 
      if (iprint() .ge. 50) then
        npol = 6
        call politg(ri,rhoL,nr,npol,0,errmx,bessi)
        q = bessi(nr)/y0
        do  40  ir = 1, nr
   40   bessw(ir) = vL(ir,1,isp)*ri(ir)**2*e/atepi
        call politg(ri,bessw,nr,npol,0,errmx,bessi)
        qscr = e*bessi(nr)/y0
        print 333, z,e,q,qeff/y0,qscr
  333   format(' lscoul:  z=',f10.6,'  e-TF=',f8.4,'  quscr=',f10.6,
     .    '  qeff=',f10.6,'  qscr=',f11.6)
C        call prrmsh('vL-scr',ri,vL,nr,nr,nlml)
C        call prrmsh('rho-scr',ri,bessw,nr,nr,1)
      endif

C... Debugging: (-nabla - e) phi = 8 pi rhoL
C      call fpiint(-np,0,i,p,wp)
C      if (i .ne. np) stop 'bug'
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C      call gradfl(lmax,nd,nr,np,1,nr,1,10,npol,ri,yl,gyl,wk,
C     .    wk2,vL,grp,ggrp)
C      call fp2yl(nr,nlml,np,wp,ggrp,yl,ggvl)
C      do  46  ilm = 1, nlml
C      do  46  ir = 1, nr
C   46 ggvl(ir,ilm) =  (-ggvl(ir,ilm) - e*vL(ir,ilm,isp))*ri(ir)**2/atepi
C      call prrmsh('rhoL',ri,rhoL,nr,nr,nlml)
C      call prrmsh('vL-GF',ri,vL,nr,nr,nlml)
C      call prrmsh('(-nabla - e)vL-GF',ri,ggvl,nr,nr,nlml)

C   ... Use a Jacobian for a shifted mesh.  Doesn't improve anything.
C       ri2(ir) = ir-1
C       bessw(ir) = bessw(ir) * a*(ri(ir)+b)
C       hankw(ir) = hankw(ir) * a*(ri(ir)+b)
C     ri2(ir) = 0
C     bessw(1) = bessw(1) * a*(ri(1)+b)
C     hankw(1) = hankw(1) * a*(ri(1)+b)
C     call politg(ri2,bessw,nr,npol,0,errmx,bessi)
      end
      subroutine fp2yl(nr,nlm,np,wp,fp,yl,fl)
C- Yl-projection of function tabulated on a mesh
C     implicit none
      integer nr,nlm,np,ip,ilm,ir
      double precision fl(nr,nlm),fp(nr,np),yl(np,nlm),wp(np),xx

      call dpzero(fl,nr*nlm)
      do  20  ip = 1, np
      do  20  ilm = 1, nlm
      xx = wp(ip)*yl(ip,ilm)
      do  20  ir = 1, nr
   20 fl(ir,ilm) = fl(ir,ilm) + fp(ir,ip)*xx
      end

