      subroutine gklq(lmax,rsm,q,p,e,kmax,k0,alat,dlv,nkd,nrx,yl,wk,job,
     .gkl)
C- Bloch sum of k,L-dependent gaussians (vectorizes)
C ---------------------------------------------------------------
Ci Inputs:
Ci  lmax   :l-cutoff for gkl
Ci   rsm   :smoothing radius
Ci   q     :wave number for Bloch sum (units of 2*pi/alat)
Ci   p     :connecting vector (units of alat)
Ci   e     :G_kL scaled by exp(e*rsm**2/4)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci  dlv,nkd:direct lattice vectors, and number
Ci   nrx   :leading dimension of wk,yl
Ci   yl    :work array, dimensioned nrx*(lmax+1)**2
Ci   wk    :work array, dimensioned at least nrx*(2*lmax+10)
Ci   nrx   :dimensions work arrays yl, and must be >= max(nkq,nkd)
Ci   job   :1s digit
Ci         :0, generate wk(1,2,3,4)
Ci          1  assume wk(1,3,4) and yl have been generated already
Ci          2  assume wk(1,2,3,4) and yl have been generated already
Ci         :10s digit
Ci         :0  use standard phase convention, do not shorten p
Ci         :1  use standard phase convention, but shorten p
Ci         :2  scale standard phase convention by exp(-i q . p)
Ci         :3  like 2, but also shorten p
Ci   k0    :leading dimension of gkl
Co Outputs:
Co   yl:  ylm(1..nkd,1..(lmax+1)**2) for points alat*(p-dlv(1..nkd))
Co   wk:  (*,1) holds r**2
Co        (*,2) holds Y0 exp(-(r/rsm)**2)
Co        (*,3) holds cos(q.dlv)
Co        (*,4) holds sin(q.dlv)
Co   gkl: G_kL * exp(e*rsm**2/4) generated for (0:kmax,0:lmax)
Ce External routines required: ropyln
Cu Updates
Cu  15 Aug 00 extended to e>0; added 1000 digit job
C ---------------------------------------------------------------
C     implicit none
      integer k0,kmax,lmax,nkd,nrx,job
      double precision alat,rsm,p(3),q(3),dlv(3,nkd),wk(nrx,4),yl(nrx,1)
      double precision gkl(2,0:k0,1),e
C Local variables
      integer ilm,ir,k,l,m,nlm,ik1,ik2,lc1,ls1,lc2,ls2,job0,job1
      double precision qdotr,pi,tpi,y0,ta2,x,y,a2,g0fac,xx1,xx2,x1,x2,
     .y2,p1(3),sp,cosp,sinp

C --- Setup ---
      if (kmax .lt. 0 .or. lmax .lt. 0 .or. rsm .eq. 0d0) return
      job0 = mod(job,10)
      job1 = mod(job/10,10)
      nlm = (lmax+1)**2
      pi  = 4*datan(1d0)
      tpi = 2*pi
      y0  = 1/dsqrt(4*pi)
      a2  = 1/rsm**2
      ta2 = 2*a2
      do  10  ilm = 1, nlm
        do  10  k = 0, kmax
          gkl(1,k,ilm) = 0d0
          gkl(2,k,ilm) = 0d0
   10 continue
C ... Shorten connecting vector; need to adjust phase later
      if (job1 .eq. 1 .or. job1 .eq. 3) then
        call shortn(p,p1,dlv,nkd)
      else
        call dcopy(3,p,1,p1,1)
      endif

C --- Put ylm in yl and alat**2*(p-dlv)**2 in wk(1) ---
      if (job0 .eq. 0) then
        do  20  ir = 1, nkd
          wk(ir,2) = alat*(p1(1)-dlv(1,ir))
          wk(ir,3) = alat*(p1(2)-dlv(2,ir))
          wk(ir,4) = alat*(p1(3)-dlv(3,ir))
   20   continue
        call ropyln(nkd,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
C ...   cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
        do  22  ir = 1, nkd
          qdotr = 2*pi*(q(1)*dlv(1,ir)+ q(2)*dlv(2,ir)+ q(3)*dlv(3,ir))
          wk(ir,3) = dcos(qdotr)
          wk(ir,4) = dsin(qdotr)
          wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
   22   continue
      elseif (job0 .eq. 1) then
        do  24  ir = 1, nkd
          wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
   24   continue
      endif

      lc1 = 5
      ls1 = 6
      lc2 = 7
      ls2 = 8
      ik1 = 9
      ik2 = 10+lmax

C --- Outer loop over k (in blocks of 2), and over l ---
      do  30  k = 0, kmax, 2
        do  30  l = 0, lmax
          g0fac = 1/rsm*ta2**(l+1)/pi * dexp(e*rsm*rsm/4)

C   ... Make radial part of the G_kl(1..nkd) for k= 0, 1
          if (k .eq. 0) then
            do  32  ir = 1, nkd
              xx1 = g0fac*wk(ir,2)
              xx2 = (ta2*wk(ir,1)-3-2*l)* ta2 * xx1
              wk(ir,ik1+l) = xx1
              wk(ir,ik2+l) = xx2
              wk(ir,lc1) = wk(ir,3)*xx1
              wk(ir,ls1) = wk(ir,4)*xx1
              wk(ir,lc2) = wk(ir,3)*xx2
              wk(ir,ls2) = wk(ir,4)*xx2
   32       continue
C   ... Make radial part of the G_kl(1..nkd) for k, k+1 from k-1, k-2
C       and cos(q.dlv) * G_kl and sin(q.dlv) * G_kl
          else
            x = 2*(k-1)*(2*k + 2*l-1)
            y = 4*k + 2*l-1
            x2 = 2*k*(2*(k+1) + 2*l-1)
            y2 = 4*(k+1) + 2*l-1
            do  34  ir = 1, nkd
              xx1 = ta2*((ta2*wk(ir,1)-y)*wk(ir,ik2+l) - x*ta2*wk(ir,ik1+l))
              xx2 = ta2*((ta2*wk(ir,1)-y2)*xx1         -x2*ta2*wk(ir,ik2+l))
              wk(ir,ik1+l) = xx1
              wk(ir,ik2+l) = xx2
              wk(ir,lc1) = wk(ir,3)*xx1
              wk(ir,ls1) = wk(ir,4)*xx1
              wk(ir,lc2) = wk(ir,3)*xx2
              wk(ir,ls2) = wk(ir,4)*xx2
   34       continue
          endif

C   ... For each point, add G_kl Y_L exp(i q.dlv) into Bloch G_kL
          ilm = l*l
          if (k .lt. kmax) then
            do  36  m = -l, l
              ilm = ilm+1
              do  38  ir = nkd, 1, -1
                gkl(1,k,ilm) = gkl(1,k,ilm) + wk(ir,lc1)*yl(ir,ilm)
                gkl(2,k,ilm) = gkl(2,k,ilm) + wk(ir,ls1)*yl(ir,ilm)
                gkl(1,k+1,ilm) = gkl(1,k+1,ilm) + wk(ir,lc2)*yl(ir,ilm)
                gkl(2,k+1,ilm) = gkl(2,k+1,ilm) + wk(ir,ls2)*yl(ir,ilm)
   38         continue
   36       continue
          else
            do  46  m = -l, l
              ilm = ilm+1
              do  48  ir = nkd, 1, -1
                gkl(1,k,ilm) = gkl(1,k,ilm) + wk(ir,lc1)*yl(ir,ilm)
                gkl(2,k,ilm) = gkl(2,k,ilm) + wk(ir,ls1)*yl(ir,ilm)
   48         continue
   46       continue
          endif
   30 continue

C ... Put in phase to undo shortening, or different phase convention
      sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      if (job1 .ge. 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
      if (sp .ne. 0d0) then
        cosp = dcos(sp)
        sinp = dsin(sp)
        do  40  ilm = 1, nlm
          do  40  k   = 0, kmax
            x1 = gkl(1,k,ilm)
            x2 = gkl(2,k,ilm)
            gkl(1,k,ilm) = x1*cosp - x2*sinp
            gkl(2,k,ilm) = x2*cosp + x1*sinp
   40   continue
      endif

      end

