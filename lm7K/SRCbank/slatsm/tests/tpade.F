      subroutine emesh(semsh,z,w)
C- Energy mesh for energy contour integration
C ----------------------------------------------------------------
Cr  Structure of semsh:
Cr    1:   nz  number of points
Cr    2:   modec for contour; see zmesh
Cr    3-4: emin,emax
Cr    5-6: ecc,eps: (for modec=10) passed to zmesh
Cr    5:   delta    (for modec=0) passed to zmesh
C ----------------------------------------------------------------
      implicit none
C Passed variables
      double precision semsh(6)
C Local variables
      integer nz,modec
      double precision z(2,1),w(2,1)
      double precision emin,emax,ecc,eps,delta

      nz    = semsh(1)
      modec = semsh(2)
      emin  = semsh(3)
      emax  = semsh(4)
      delta = semsh(5)
      ecc   = semsh(5)
      eps   = semsh(6)
      call zmesh(modec,emin,emax,ecc,eps,delta,nz,z,w)

      end
      subroutine zmesh(modec,eb,ef,ecc,eps,delta,nz,z,w)
C- Makes a complex energy mesh for contour integration
C ----------------------------------------------------------------
Ci modec:  0 uniform mesh shifted off the real axis by i*delta
Ci         1 ditto, but use -i*delta
Ci        10 elliptical contour, top half of complex plane
Ci        11 elliptical contour, bottom half of complex plane
Ci delta: imaginary part for uniform mesh
Ci eb,ef: integration from eb to ef
Ci ecc:   eccentricity for ellipse (0 => circle, 1 => line)
Ci eps:   parameter controlling logarithmic bunching of points near Ef
Ci        eps=0 => no bunching
Ci nz:    number of points
Co Outputs
Co   z,w: complex points and weights for contour integration
Cr Remarks
Cr   Contour integration by Gaussian quadrature of an angle theta on
Cr   an elliptical path, for theta in [-pi,0].  When eps>0, there is
Cr   an additional mapping to a logarithmic mesh to bunch points near
Cr   the real axis.  Choose eps in [0,1], the closer to 1 the more
Cr   bunching:  theta(x)=(1-exp(x))/eps, x(theta)=log(1-theta*eps)
Cr   Points are ReZ = z0 + zr*cos(theta), Im Z = zr*(1-ecc)*sin(theta),
Cr   where z0=(ef+eb)/2 and zr=(ef-eb)/2.
Cr
Cr   For uniform mesh, weights are defined for trapezoidal integration
C ----------------------------------------------------------------
      implicit none
C Passed variables
      integer nz,modec
      double precision z(2,nz),w(2,nz)
C     double complex z(nz),w(nz)
      double precision eb,ef,ecc,eps,delta
C Local variables
      logical lellip
      integer i,lgunit,iprint,signz
      double precision wgt,expx,theta,ct,st,pi,z0,zr,rm,de,del,eta,datan
      character*80 outs

      signz = mod(modec,10)
      lellip = mod(modec/10,10) .eq. 1
      if (modec-signz .ne. 0 .and. modec-signz .ne. 10)
     .  call rx('zmesh: illegal modec')

C --- Setup and default settings ---
      pi = 4*datan(1d0)
      eta = 1-ecc
      if (ecc .le. 0d0) eta = 1d0
      del = delta
*     if (delta .eq. 0d0) del = 1d-2
      if (iprint() .ge. 20) print *, ' '
      call awrit5('%x ZMESH: nz=%i  ecc,eps=%1;4,1d,%1;4d'//
     .  '  e1,e2=%1;4,2d,%1;4,2d',outs,80,0,nz,1-eta,eps,eb,ef)
      if (.not. lellip) call awrit5('%x ZMESH:  nz=%i  delta=%1;4d'//
     .  '%?;n==1;(*-1);;'//
     .  '  e1,e2=%1;4,1d,%1;4,1d',outs,80,0,nz,del,signz,eb,ef)
      if (nz .le. 1 .or. lellip .and. nz .le. 0 .or. iprint() .ge. 20)
     .  call awrit0('%a',outs,-80,-lgunit(1))
      if (nz .le. 1 .or. lellip .and. nz .le. 0)
     .  call rx('... illegal mesh')

C --- Contour case ---
      if (lellip) then
C   ... Gaussian quadrature on interval [0,log(1-eps)] or [0,1] if eps=0
C   ... Use z, and w as a work space
        call mklegw(nz,z,z(1+nz,1),0)
C   ... Flip to  put z of highest last
        call dscal(nz,-1d0,z,1)
C   ... Copy (pts,wgt) to (real, imag) of array w for temp storage
        do  5  i = 1, nz
        w(1,i) = z(i,1)
    5   w(2,i) = z(i+nz,1)
        rm = 1
        if (eps .ne. 0) rm = dlog(1d0-eps)
        do  10  i = 1, nz
        w(1,i) = rm/2*w(1,i) + rm/2
   10   w(2,i) = w(2,i)*rm/2
C   ... Additional mapping to to logarithmic mesh in [0,1]
        if (eps .ne. 0) then
          do  20  i = 1, nz
          expx = dexp(w(1,i))
          w(1,i) = (1d0-expx)/eps
   20     w(2,i) = -expx*w(2,i)/eps
        endif
C   ... Map to points on an ellipse
        zr = (ef-eb)/2
        z0 = (ef+eb)/2
        do  30  i = 1, nz
          theta = -pi + pi*w(1,i)
          wgt = w(2,i)
          ct = dcos(theta)
          st = dsin(theta)
          z(1,i) = z0 + zr*ct
          z(2,i) = -zr*eta*st
          w(1,i) = -pi*zr*st*wgt
          w(2,i) = -pi*zr*eta*ct*wgt
          if (signz .eq. 1) then
            z(2,i) = -z(2,i)
            w(2,i) = -w(2,i)
          endif
   30   continue
C --- Uniformly spaced points shifted by delta off the real axis---
      else
        if (nz .le. 1) nz = 30
        de = (ef-eb)/(nz-1)
        do  40  i = 1, nz
          z(1,i) = eb + (i-1)*de
          z(2,i) = del
          if (signz .eq. 1) z(2,i) = -z(2,i)
          w(1,i) = de
          w(2,i) = 0d0
   40   continue
        w(1,1) = w(1,1)/2
        w(1,nz) = w(1,nz)/2
      endif

C --- Printout ---
      if (iprint() .ge. 60 .and. lellip) then
        print 332
  332   format(10x,'z',21x,'w')
        do  50  i = 1, nz
   50   print 333, z(1,i),z(2,i),w(1,i),w(2,i)
  333   format(2f10.6,2x,2f10.6)
      endif
      end
      subroutine fmain
C Tests zmesh and Pade approximations to GF
C There are three test cases: simple pole, simple s-band,
C and semi-infinite s-band.  Tests for integrated properties
C are written to stdout; also tests for DOS, which are calculated
C directly and estimated by a Pade approximant.  These
C are written to file dos.dat.
C
C To check stdout and dos with reference
C echo / | a.out | tee out
C diff out out.pade; mc -f'(7f12.6/12x,6f12.6)' dos.dat dos.pade -- -px| head
C
C To compare the exact DOS with the Pade approximant for these cases:
C fplot -frme 0,1,1,1.5 -lt 2 -colsy 2 dos.dat -colsy 8 -lt 3,bold=6 dos.dat\
C       -frme 0,1,.5,1 -lt 2 -colsy 4 dos.dat -colsy 10 -lt 3,bold=6 dos.dat\
C       -frme 0,1,0,.5 -lt 2 -colsy 6 dos.dat -colsy 12 -lt 3,bold=6 dos.dat
      implicit none
      integer nzmx,nz,npad
      parameter(nzmx=200)
      double complex zm(nzmx),wz(nzmx),g(nzmx),sum,z3,pdot,dcsqrt
      double complex zpad(nzmx),wpad(nzmx),gp(nzmx),gb(nzmx),gs(nzmx)
      double precision eb,ef,ecc,eps,delta,result,pi,fac,xx
      double precision dosb(nzmx),doss(nzmx),gfrb(nzmx),gfrs(nzmx),
     .  dosp(nzmx),gfrp(nzmx),cpad(nzmx,nzmx+2)
      double precision pole,psmall
      integer i,ifi,jfi,fopng,sgn
      pi = 4*datan(1d0)
      eb = -2.1d0
      ef = 0d0
      nz = 20
      ecc = .5d0
      eps = 0d0
      pole = -.1d0
      delta =  .001d0
      psmall = .02d0
      sgn = -1
      print *, 'input nz,ecc,eps,pole,psmall,delta'
      read(*,*) nz,ecc,eps,pole,psmall,delta
      npad = 200
      i = (1+sgn)/2
      call zmesh(i,1.5d0*eb,-1.5d0*eb,0d0,0d0,delta,npad,zpad,wpad)
      call pshpr(80)
      call zmesh(10+i,eb,ef,ecc,eps,delta,nz,zm,wz)
      ifi = -1
      ifi = fopng('zpoints.dat',ifi,0)
      do  5  i = 1, nz
      write(ifi,100) dble(zm(i)), dimag(zm(i))
    5 continue
  100 format(2e16.7)
      call awrit6(' nz=%i  ecc=%d  eps=%d  pole=%d  psmall=%d'//
     .  '  delta=%d',' ',80,6,nz,ecc,eps,pole,psmall,delta)

c --- test 1: simple pole 1/(z-a) with pole a inside integration interval
c     (Test is strictest when pole close to fermi energy)
c     pole = -0.01d0
      sum = dcmplx(0d0,0d0)
      xx = 0
      do  10  i = 1, nz
        g(i) = 1d0/(zm(i)-pole)
        xx = max(xx,cdabs(g(i)))
        sum = sum + g(i)*wz(i)
   10 continue
      result = sgn*dimag(sum)/pi
      call awrit3(' Simple pole at e=%1;4d: gmax='//
     .  '%1;3g  sum = %1;10d',' ',80,6,pole,xx,result)
      call padcof(nz,zm,g,nz,cpad)
      call pade(npad,zpad,nz,zm,nz,cpad,gp)

c --- test 2: contour integration of full-chain s-band gf x 2 for spin
c     integrated up to half band should give 1
      sum = dcmplx(0d0,0d0)
      xx = 0d0
      do  20  i = 1, nz
        fac = -dsign(1d0,dble(zm(i)))
        z3 = zm(i) + zm(i)**3*psmall
        pdot = 1 + 3*zm(i)*zm(i)*psmall
        g(i) = -2d0*fac*pdot/dcsqrt(z3**2-4d0)
        xx = max(xx,cdabs(g(i)))
        sum = sum + g(i)*wz(i)
   20 continue
      result = sgn*dimag(sum)/pi
      call awrit2(' s band chain: gmax='//
     .  '%1;3g  sum = %1;10d',' ',80,6,xx,result)
      call padcof(nz,zm,g,nz,cpad)
      call pade(npad,zpad,nz,zm,nz,cpad,gb)

c --- test 3: integration of semi-infinite s-band gf x 2 for spin
      sum = dcmplx(0d0,0d0)
      xx = 0d0
      do  30  i = 1, nz
        fac = -dsign(1d0,dble(zm(i)))
        z3 = zm(i) + zm(i)**3*psmall
        pdot = 1 + 3*zm(i)*zm(i)*psmall
        g(i) = pdot*(z3 + fac*dcsqrt(z3**2 - 4))
        xx = max(xx,cdabs(g(i)))
        sum = sum + g(i)*wz(i)
   30 continue
      result = sgn*dimag(sum)/pi
      call awrit2(' semi-infinite s-band: gmax='//
     .  '%1;3g  sum = %1;10d',' ',80,6,xx,result)
      call padcof(nz,zm,g,nz,cpad)
      call pade(npad,zpad,nz,zm,nz,cpad,gs)

C --- DOS on a fine mesh ---
      nz = npad
      do  50  i = 1, nz
        g(i) = 1d0/(zpad(i)-pole)
        dosp(i) = sgn*dimag(g(i))/pi
        gfrp(i) = dble(g(i))/pi
        z3 = zpad(i) + zpad(i)**3*psmall
        pdot = 1 + 3*zpad(i)*zpad(i)*psmall
        fac = -dsign(1d0,dble(zpad(i)))
        g(i) = -2*fac*pdot/dcsqrt(z3**2 - 4)
        dosb(i) = sgn*dimag(g(i))/pi
        gfrb(i) = dble(g(i))/pi
        g(i) = pdot*(z3 + fac*dcsqrt(z3**2 - 4))
        doss(i) = sgn*dimag(g(i))/pi
        gfrs(i) = dble(g(i))/pi
   50 continue
      print *, 'dumping dos in file dos.dat ...'
      jfi = -1
      jfi = fopng('dos.dat',jfi,0)
      rewind jfi
      call awrit2('%% rows %i cols %i',' ',80,jfi,nz,7+6)
      write(jfi,301)
  301 format('#',5x,'z',9x,'dos(pole)  Re(G)',8x,
     .  'dos(sband) Re(G)',8x,'dos(semi-inft) Re(G)/pade')
      do  60  i = 1, nz
   60 write(jfi,300) dble(zpad(i)),
     .    dosp(i),gfrp(i),dosb(i),gfrb(i),doss(i),gfrs(i),
     .    sgn*dimag(gp(i))/pi,dble(gp(i))/pi,
     .    sgn*dimag(gb(i))/pi,dble(gb(i))/pi,
     .    sgn*dimag(gs(i))/pi,dble(gs(i))/pi
      call fclose(ifi)
      call fclose(jfi)
  300 format(7f12.6/12x,6f12.6)
      call cexit(1,1)
      end
