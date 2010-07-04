#define VANILLA 1
#define NOQUAD 1
CSFCPP#define F90 1
#define FFTW 1
#define IN_PLACE 1
      subroutine gradfl(lmax,nd,nr,np,ir0,ir1,lgg,lx,nn,ri,yl,gyl,fl,
     .  gp,ggp)
C- Gradient, Laplacian of function point-wise through sphere from YL expansion
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  :density is expanded to l-cutoff lmax
Ci   nd    :dimensions yl,gyl
Ci   nr    :number of radial mesh points
Ci   np    :number of angular mesh points
Ci   ir0   :gp and gpp are made for radial points between ir0,ir1
Ci   ir1   :gp and gpp are made for radial points between ir0,ir1
Ci   lgg   :if zero, make gradient gp only; gpp not addressed.
Ci   lx    :(ones digit) if 1, fl scaled by r**2 
Ci         :(tens digit): extrapolate 1st point (ir0=1) from others
Ci         :(100  digit): rational function interpolation for radial deriv
Ci   nn    :nn: number of points used to differentite radial f
Ci   ri    :vector of radial mesh points
Ci   yl    :Spherical harmonics for L=0:(lmax+1)^2 at each angular mesh point
Ci         :Generate with a call to ropyln
Ci   gyl   :Gradient of YL.
Ci         :Generate with a call to ropylg
Ci   fl    :function to be differentiated, on the combined radial 
Ci         :and angular mesh
Co Outputs
Ci   gp    :gradient of fl, on the combined radial and angular mesh,
Ci         :x,y,z components
Ci   ggp   :Laplacian of fl, on the combined radial and angular mesh
Cl Local variables
Cl   gf    :Work array (used for radial derivative of fl)
Cl   ggf   :Work array (used for 2nd radial derivative of fl)
Cr Remarks
Cr   
Cu Updates
Cu   02 Apr 09 Made gf,ggf local; fixed bug for 10s digit lx case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer np,nr,nd,lx,nn,lmax,ir0,ir1
      double precision fl(nr,1),gp(ir0:ir1,np,3),ggp(ir0:ir1,np),
     .  ri(nr),yl(np,nd),gyl(np,nd,3)
C ... Local parameters
      integer i0,ilm,ip,ir,j,l,m,l2,lerr,lgg,iprint,jx,nx
      double precision xx,cy1,tol,egf0,gf(nr),ggf(nr)
      logical lrat
      parameter (tol=1d-12)

      if (ir0 .lt. 1) call rx('gradfl: illegal value of ir0')
      cy1 = dsqrt(3/(16*datan(1d0)))
      l2 = mod(lx,100)
      lrat = lx/100 .ne. 0
      i0 = 1
      nx = nr
      if (l2/10 .ne. 0) then
        i0 = 2
        nx = nr-1
      endif

C --- Contribution (grad fl) r^-l yl ---
      call dpzero(gp,    (ir1-ir0+1)*np*3)
      if (lgg .ne. 0) call dpzero(ggp,   (ir1-ir0+1)*np)
      ilm = 0
      do  20  l = 0, lmax
      do  20  m = -l, l
      ilm = ilm+1
      if (mod(l2,10) .eq. 0) then
        call poldvm(ri(i0),fl(i0,ilm),nx,nn,lrat,tol,lerr,gf(i0))
        if (lerr .ne. 0) goto 99
      else
        do  22  ir = i0, nr
   22   ggf(ir) = fl(ir,ilm)/ri(ir)**2
        call poldvm(ri(i0),ggf(i0),nx,nn,lrat,tol,lerr,gf(i0))
        if (lerr .ne. 0) goto 99
      endif
C     Extrapolate gf to first point
      if (l2/10 .ne. 0) then
        jx = 1
        call polint(ri(2),gf(2),nx,nn,ri,0d0,0,jx,gf,egf0)
        lerr = 1
        if (iprint() .ge. 40 .and.
     .      dabs(egf0) .gt. 1d-3*max(dabs(gf(1)),dabs(gf(2)))) then
          call info5(40,0,0,' gradfl (warning): uncertainty in grad'//
     .    ' f(r=0,L=%i):  f=%;3g  est err= %;3g',ilm,gf(1),egf0,0,0)
        endif
      endif
      do  24  ip = 1, np
      do  24  ir = ir0, ir1
   24 gp(ir,ip,1) = gp(ir,ip,1) + gf(ir)*yl(ip,ilm)

C --- Laplacian: (nabla fl) Yl + fl (nabla Yl) ---
      if (lgg .ne. 0) then
        call poldvm(ri(i0),gf(i0),nx,nn,lrat,tol,lerr,ggf(i0))
        if (lerr .ne. 0) goto 99
        if (mod(l2,10) .eq. 0) then
          do  25  ir = i0, nr
          xx = 1/ri(ir)
   25     ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx*xx
        else
          do  27  ir = i0, nr
          xx = 1/ri(ir)
   27     ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx**4
        endif
        if (i0 .eq. 2) ggf(1)= (ri(3)*ggf(2)-ri(2)*ggf(3))/(ri(3)-ri(2))
        do  28  ip = 1, np
        do  28  ir = ir0, ir1
   28   ggp(ir,ip) = ggp(ir,ip) + ggf(ir)*yl(ip,ilm)
      endif
   20 continue

C ... Split grad r- into x,y,z- components
      do  26  j = 3, 1, -1
      do  26  ip = 1, np
      xx = yl(ip,j)/cy1
      if (j .eq. 1) xx = yl(ip,4)/cy1
      do  26  ir = ir0, ir1
   26 gp(ir,ip,j) = xx*gp(ir,ip,1)

C --- Contribution fl(r) grad r^-l yl (use gf as work array) ---
      ilm = 0
      do  10  l = 0, lmax
      do  10  m = -l, l
        ilm = ilm+1
C ...   Factor 1/r from grad Yl
        if (mod(l2,10) .eq. 0) then
          do  13  ir = max(i0,ir0), nr
   13     gf(ir) = fl(ir,ilm)/ri(ir)
        else
          do  15  ir = max(i0,ir0), nr
   15     gf(ir) = fl(ir,ilm)/ri(ir)**3
        endif
        if (i0 .gt. ir0) gf(1) = (ri(3)*gf(2)-ri(2)*gf(3))/(ri(3)-ri(2))
        do  18  j = 1, 3
        do  18  ip = 1, np
        xx = gyl(ip,ilm,j)
        do  18  ir = ir0, ir1
   18   gp(ir,ip,j) = gp(ir,ip,j) + gf(ir)*xx
   10 continue
      return

C --- Error handling ---
   99 print *, 'gradfl: stopping at ilm=',ilm,'  point', lerr
      call rx('gradfl: can''t diff radial function')

      end
C$$$#if TEST
C$$$C Test program to check
C$$$      subroutine fmain
C$$$      implicit none
C$$$      integer nr,nlmx,lmax,nlm1,ir,ii,i,l,ilm,i1,i2,nsize,ll
C$$$      integer nlmf,nnn,np,nph,nth,oxp,oyp,ozp,oyl,ogyl,orp,or2,
C$$$     .  oagfl,oggfl,ogfl,ogp,oggp,oggpb,orh,orhol,ofrl,lp,ilm2
C$$$      parameter (nr=250,nlmx=49,nsize=500000,nnn=144)
C$$$      double precision p(3,nnn),wp(nnn),rofi(nr),a,b,rmax,scl,scl2
C$$$      real w(nsize)
C$$$      common /w/ w
C$$$      call wkinit(nsize)

C$$$   99 print *, 'ilm1, scl1, ilm2, scl2:'
C$$$      ilm = 2
C$$$      scl = 1
C$$$      ilm2 = 5
C$$$      scl2 = 0
C$$$      read(*,*) ilm, scl, ilm2, scl2
C$$$      lmax = ll(max(ilm,ilm2))+2

C$$$C ... Angular mesh
C$$$      nth=lmax+1
C$$$      nph=2*nth
C$$$      nlmf = (lmax+1)**2
C$$$      call fpiint(nth,nph,np,p,wp)
C$$$      print *, nr, 'radial points;', np, ' angular points'
C$$$      call defrr(oxp,     np)
C$$$      call defrr(oyp,     np)
C$$$      call defrr(ozp,     np)
C$$$      call defrr(or2,     np)
C$$$C ... 3 necessary if two derivatives taken ?
C$$$      call defrr(oyl,     (lmax+3)**2*np)
C$$$      call dcopy(np,p(1,1),3,w(oxp),1)
C$$$      call dcopy(np,p(2,1),3,w(oyp),1)
C$$$      call dcopy(np,p(3,1),3,w(ozp),1)
C$$$      call defrr(ogyl,     nlmf*np*3)

C$$$C ... Radial mesh
C$$$      rmax = .1d0
C$$$      a = 1d-6
C$$$      a = .001
C$$$      b = rmax/(dexp(a*(nr-1))-1d0)
C$$$      call radmsh(rmax,a,nr,rofi)
C$$$C ... Other setup
C$$$      call defrr(orp,     nr*np)
C$$$      call defrr(ogp,     nr*np*3)
C$$$      call defrr(oggp,    nr*np)
C$$$      call defrr(oggpb,   nr*np*3*3)
C$$$      call defrr(ogfl,    nr)
C$$$      call defrr(oggfl,   nr)
C$$$      call defrr(orhol,   nr*nlmf)
C$$$      call defrr(ofrl,    nr*nlmf)

C$$$      call testg(rofi,nr,np,ilm,scl,ilm2,scl2,nlmf,w(orp),
C$$$     .  w(ogp),w(oggp),w(oggpb),w(ofrl),w(oyl),w(ogyl),wp,w(orhol),
C$$$     .  w(ogfl),w(oggfl),w(oxp),w(oyp),w(ozp),w(or2))

C$$$      end
C$$$      subroutine radmsh(r,a,nr,rofi)
C$$$      implicit real*8 (a-h,p-z), integer (o)
C$$$      dimension rofi(nr)
C$$$      b=r/(dexp(a*nr-a)-1.d0)
C$$$      do 1 ir=1,nr
C$$$    1 rofi(ir)=b*(dexp(a*ir-a)-1d0)
C$$$      end
C$$$      subroutine fp2yl(nr,nlm,np,wp,fp,yl,fl)
C$$$C- Yl-projection of function tabulated on a mesh
C$$$      implicit none
C$$$      integer nr,nlm,np,ip,ilm,ir
C$$$      double precision fl(nr,nlm),fp(nr,np),yl(np,nlm),wp(np),xx

C$$$      call dpzero(fl,nr*nlm)
C$$$      do  20  ip = 1, np
C$$$      do  20  ilm = 1, nlm
C$$$      xx = wp(ip)*yl(ip,ilm)
C$$$      do  20  ir = 1, nr
C$$$   20 fl(ir,ilm) = fl(ir,ilm) + fp(ir,ip)*xx
C$$$      end
C$$$      subroutine prmr(strn,nr,nrx,rofi,pow,f,nl)
C$$$      implicit none
C$$$      integer nr,nrx,nl,ir,j,fopna,ifi
C$$$      double precision rofi(nrx),f(nrx,nl),pow
C$$$      character*(10) fmt, strn*(*)
C$$$      ifi = fopna('out',19,0)
C$$$      write(ifi,*) nr, nl+1,   ' scaled by', pow
C$$$      do  10  ir = 1, nr
C$$$        write(ifi,333) rofi(ir),
C$$$     .    (f(ir,j)*(rofi(ir)+1d-12)**pow, j=1, nl)
C$$$C  333   format(f12.5,(7g18.10:/12x))
C$$$  333   format(f12.5,(17f12.6:/12x))
C$$$   10 continue
C$$$      call fclose(ifi)
C$$$      print *, strn
C$$$      pause
C$$$      end
C$$$      subroutine testg(rofi,nr,np,ilm1,scl1,ilm2,scl2,nlm,rp,gp,ggp,
C$$$     .  ggpb,frl,yl,gyl,wp,rl,gfl,ggfl,xp,yp,zp,r2)

C$$$      implicit none
C$$$      integer nr,np,ilm1,ilm2,nlm
C$$$      double precision rofi(1),yl(np,nlm),gyl(np,nlm,3),
C$$$     .  wp(np),rl(nr,nlm),frl(nr,nlm),gfl(nr),ggfl(nr),
C$$$     .  rp(nr,np),gp(nr,np,3),ggp(nr,np),ggpb(nr,np,3,3),x2,
C$$$     .  xp(1),yp(1),zp(1),r2(1),xx,scl1,scl2,phi(0:20),psi(0:20),e,pow
C$$$      integer ip,ir,ll,lmax,nlmx,lx1,lx2,lr2,ir0,ir1
C$$$      logical lrat

C$$$      ir0 = 247
C$$$      ir1 = 249
C$$$      lrat = .false.
C$$$      lr2 = 10
C$$$      lmax = ll(nlm)
C$$$      lx1 = ll(ilm1)
C$$$      lx2 = ll(ilm2)
C$$$      nlmx = min(nlm,16)
C$$$      if (ll(nlm) .lt. lx1+2) call rx('testg: need bigger nlm')
C$$$      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C$$$      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)

C$$$CC ... Show by brute force for rl = 1, laplacian is -l(l+1)/r**2
C$$$C      print *, 
C$$$C     .'Show by brute force for rl = 1, laplacian is -l(l+1)/r**2'
C$$$C      call dpzero(rl, nr*nlm)
C$$$C      call dcopy(nr,scl1,0,rl(1,ilm1),1)
C$$$C      do  10  ip = 1, np
C$$$C      do  10  ilm = 1, nlm
C$$$C        do  15  ir = 1, nr
C$$$C   15   rp(ir,ip) = rp(ir,ip) + rl(ir,ilm)*yl(ip,ilm)
C$$$C   10 continue
C$$$C      call blap(rofi,nr,np,nlm,rp,gp,ggp,ggpb,
C$$$C     .  lrat,lr2,frl,yl,gyl,wp,gfl,ggfl,0d0,lmax,nlmx)


C$$$C...  Gradient, laplacian of Hankel or Bessel tabulated on a mesh
C$$$      print *, 'nabla of hankel or bessel, brute force'
C$$$      e = -.7d0
C$$$      do  110  ir = 2, nr
C$$$      call bessl(e*rofi(ir)**2,max(lx1,lx2),phi,psi)
C$$$C ... here for Bessel
C$$$      lrat = .false.
C$$$      pow = -lx1
C$$$      xx = rofi(ir)**lx1
C$$$      x2 = rofi(ir)**lx2
C$$$      do  110  ip = 1, np
C$$$  110 rp(ir,ip) = scl1*phi(lx1)*xx*yl(ip,ilm1) +
C$$$     .            scl2*phi(lx2)*x2*yl(ip,ilm2)
C$$$C ... here for Hankel (note for near origin, need rational f interp.
C$$$C      lrat = .true.
C$$$C      pow = lx1+1
C$$$C      xx = rofi(ir)**(-lx1-1)
C$$$C      do  110  ip = 1, np
C$$$C  110 rp(ir,ip) = psi(lx1)*xx*yl(ip,ilm1)
C$$$      call makghl(rofi,nr,np,ilm1,nlm,gp,frl,yl,xp,yp,zp,wp,e,pow)

C$$$      if (mod(lr2,10) .ne. 0) then
C$$$        do  60  ip = 1, np
C$$$        do  60  ir = 1, nr
C$$$   60   rp(ir,ip) = rp(ir,ip)*rofi(ir)**2
C$$$      endif
C$$$      call blap(rofi,nr,np,nlm,rp,gp,ggp,ggpb,
C$$$     .  ir0,ir1,lrat,lr2,frl,yl,gyl,wp,gfl,ggfl,pow,lmax,nlmx)


C$$$      end
C$$$      subroutine blap(ri,nr,np,nlm,rp,gp,ggp,ggpb,
C$$$     .  ir0,ir1,lrat,lr2,fl,yl,gyl,wp,gf,ggf,pow,lmax,nlmx)
C$$$C- Laplacian by brute force of function tabulated on a mesh
C$$$      implicit none
C$$$      integer nr,np,nlm,lr2,ir0,ir1
C$$$      logical lrat
C$$$      double precision ri(1),yl(np,nlm),gyl(np,nlm,3),
C$$$     .  wp(np),fl(nr,nlm),gf(nr),ggf(nr),
C$$$     .  rp(nr,np),gp(nr,np,3),ggp(nr,np),ggpb(ir0:ir1,np,3,3),xx,pow
C$$$      integer ip,ir,lmax,j,nlmx,lx,itwo(2)

C$$$      lx = 0
C$$$      if (lrat) lx = 100
C$$$      call fp2yl(nr,nlm,np,wp,rp,yl,fl)
C$$$      call prmr('fl made from points...',nr,nr,ri,pow,fl,nlmx)
C$$$C ... Gradient and Laplacian
C$$$      call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C$$$      call gradfl(lmax,nlm,nr,np,1,nr,1,lx+lr2,8,ri,yl,gyl,gf,ggf,fl,
C$$$     .  gp,ggp)
C$$$      call fp2yl(nr,nlm,np,wp,ggp,yl,fl)
C$$$      call prmr('Laplacian from gradfl ...',nr,nr,ri,pow,fl,nlmx)
C$$$C ... grad (gradient) ... points ir0:ir1 only
C$$$      do  12  j = 1, 3
C$$$        call fp2yl(nr,nlm,np,wp,gp(1,1,j),yl,fl)
C$$$        call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C$$$          call prmr('gp',nr,nr,ri,pow+1,fl,nlmx)
C$$$        call gradfl(lmax,nlm,nr,np,ir0,ir1,0,lx+10,8,ri,yl,gyl,
C$$$     .    gf,ggf,fl,ggpb(ir0,1,1,j),ggp)
C$$$        call fp2yl(ir1-ir0+1,nlm,np,wp,ggpb(ir0,1,j,j),yl,fl)
C$$$        call prmr('ggpb',ir1-ir0+1,ir1-ir0+1,ri(ir0),pow+2,fl,nlmx)
C$$$   12 continue
C$$$      do  16  ip = 1, np
C$$$      do  16  ir = ir0, ir1
C$$$      xx = ggpb(ir,ip,1,1) + ggpb(ir,ip,2,2) + ggpb(ir,ip,3,3)
C$$$   16 ggpb(ir,ip,1,1) = xx
C$$$      call fp2yl(ir1-ir0+1,nlm,np,wp,ggpb,yl,fl)
C$$$      call prmr('Laplacian by grad(grad)',ir1-ir0+1,ir1-ir0+1,ri(ir0),
C$$$     .  pow,fl,min(nlm,16))

C$$$      end
C$$$      subroutine makghl(ri,nr,np,ilm1,nlm,gp,fl,yl,xp,yp,zp,wp,e,pow)
C$$$      implicit none
C$$$      integer nr,np,nlm,ilm1
C$$$      double precision ri(nr),gp(nr,np,3),
C$$$     .  fl(nr,nlm),yl(np,nlm),xp(np),yp(np),zp(np),wp(np),e,pow
C$$$      integer ndim,ir,ip,j,lmax,ll,nlmx
C$$$      parameter (ndim=200)
C$$$      double precision dr(3),cy(16**2),hl(ndim),ghl(ndim,3),
C$$$     .   hd(ndim),ghd(ndim,3)
C$$$      common /static/ cy


C$$$      nlmx = min(nlm,16)
C$$$      call sylmnc(cy,15)
C$$$      lmax = ll(ilm1)
C$$$      do  10  ir = 2, nr
C$$$      do  10  ip = 1, np
C$$$        dr(1) = xp(ip)*ri(ir)
C$$$        dr(2) = yp(ip)*ri(ir)
C$$$        dr(3) = zp(ip)*ri(ir)
C$$$        call solhpg(e,dr,lmax,ndim,hl,ghl,hd,ghd,cy)
C$$$        do  10  j = 1, 3
C$$$        gp(ir,ip,j) = ghl(ilm1,j)
C$$$   10 continue

C$$$      do  20  j = 1, 3
C$$$        call fp2yl(nr,nlm,np,wp,gp(1,1,j),yl,fl)
C$$$C       call prmr('exact grad hl',nr,nr,ri,pow+1,fl,nlmx)
C$$$   20 continue

C$$$      end
C$$$      subroutine solhpg(e,dr,lmax,ndim,hl,ghl,hd,ghd,cy)
C$$$C- Solid Hankel functions with energy derivatives and gradients
C$$$      implicit real*8 (a-h,p-z), integer (o)
C$$$      dimension cy(1),dr(3),hl(ndim),ghl(ndim,3),phi(0:30),psi(0:30),
C$$$     .   hd(ndim),ghd(ndim,3)
C$$$      nlm=(lmax+1)**2
C$$$      if((lmax+2)**2.gt.ndim) call rx('solhgp: ndim too small')

C$$$C --- Make solid Hankel functions HL ---
C$$$      call sylm(dr,hl,lmax+1,r2)
C$$$      call bessl(e*r2,lmax+2,phi,psi)
C$$$      ilm=0
C$$$      fac=dsqrt(r2)
C$$$      do 10 l=0,lmax+1
C$$$        fac=fac/r2
C$$$        psidot=((l+l+1)*psi(l)-psi(l+1))/(e+e)
C$$$        do 10 m=-l, l
C$$$        ilm=ilm+1
C$$$        hd(ilm)=fac*psidot*cy(ilm)*hl(ilm)
C$$$  10    hl(ilm)=fac*psi(l)*cy(ilm)*hl(ilm)

C$$$C ------ make gradients ----------
C$$$      do 20 m=1,3
C$$$      do 20 ilm=1,nlm
C$$$      ghd(ilm,m)=0d0
C$$$  20  ghl(ilm,m)=0d0

C$$$      nlm1=lmax*lmax
C$$$      do 22 ilm=1,nlm
C$$$      call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
C$$$      ghl(ilm,1)=ghl(ilm,1)-cx1*hl(kx1)-cx2*hl(kx2)
C$$$      ghl(ilm,2)=ghl(ilm,2)-cy1*hl(ky1)-cy2*hl(ky2)
C$$$      ghl(ilm,3)=ghl(ilm,3)-cz*hl(kz)
C$$$      ghd(ilm,1)=ghd(ilm,1)-cx1*hd(kx1)-cx2*hd(kx2)
C$$$      ghd(ilm,2)=ghd(ilm,2)-cy1*hd(ky1)-cy2*hd(ky2)
C$$$      ghd(ilm,3)=ghd(ilm,3)-cz*hd(kz)
C$$$      if(ilm.le.nlm1) then
C$$$        xx=e*hl(ilm)
C$$$        ghl(kx1,1)=ghl(kx1,1)+cx1*xx
C$$$        ghl(kx2,1)=ghl(kx2,1)+cx2*xx
C$$$        ghl(ky1,2)=ghl(ky1,2)+cy1*xx
C$$$        ghl(ky2,2)=ghl(ky2,2)+cy2*xx
C$$$        ghl(kz,3)=ghl(kz,3)+cz*xx
C$$$        xx=hl(ilm)+e*hd(ilm)
C$$$        ghd(kx1,1)=ghd(kx1,1)+cx1*xx
C$$$        ghd(kx2,1)=ghd(kx2,1)+cx2*xx
C$$$        ghd(ky1,2)=ghd(ky1,2)+cy1*xx
C$$$        ghd(ky2,2)=ghd(ky2,2)+cy2*xx
C$$$        ghd(kz,3)=ghd(kz,3)+cz*xx
C$$$        endif
C$$$  22  continue
C$$$      end

C$$$#endif

