subroutine gklq(lmax,rsm,q,p,e,kmax,k0,alat,dlv,nkd,nrx,yl,wk,job, &
     gkl)
  !- Bloch sum of k,L-dependent gaussians (vectorizes)
  ! ---------------------------------------------------------------
  !i Inputs:
  !i  lmax   :l-cutoff for gkl
  !i   rsm   :smoothing radius
  !i   q     :wave number for Bloch sum (units of 2*pi/alat)
  !i   p     :connecting vector (units of alat)
  !i   e     :G_kL scaled by exp(e*rsm**2/4)
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i  dlv,nkd:direct lattice vectors, and number
  !i   nrx   :leading dimension of wk,yl
  !i   yl    :work array, dimensioned nrx*(lmax+1)**2
  !i   wk    :work array, dimensioned at least nrx*(2*lmax+10)
  !i   nrx   :dimensions work arrays yl, and must be >= max(nkq,nkd)
  !i   job   :1s digit
  !i         :0, generate wk(1,2,3,4)
  !i          1  assume wk(1,3,4) and yl have been generated already
  !i          2  assume wk(1,2,3,4) and yl have been generated already
  !i         :10s digit
  !i         :0  use standard phase convention, do not shorten p
  !i         :1  use standard phase convention, but shorten p
  !i         :2  scale standard phase convention by exp(-i q . p)
  !i         :3  like 2, but also shorten p
  !i   k0    :leading dimension of gkl
  !o Outputs:
  !o   yl:  ylm(1..nkd,1..(lmax+1)**2) for points alat*(p-dlv(1..nkd))
  !o   wk:  (*,1) holds r**2
  !o        (*,2) holds Y0 exp(-(r/rsm)**2)
  !o        (*,3) holds cos(q.dlv)
  !o        (*,4) holds sin(q.dlv)
  !o   gkl: G_kL * exp(e*rsm**2/4) generated for (0:kmax,0:lmax)
  !e External routines required: ropyln
  !u Updates
  !u  15 Aug 00 extended to e>0; added 1000 digit job
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: k0,kmax,lmax,nkd,nrx,job
  double precision :: alat,rsm,p(3),q(3),dlv(3,nkd),wk(nrx,*),yl(nrx,1)
  double precision :: gkl(2,0:k0,1),e
  ! Local variables
  integer :: ilm,ir,k,l,m,nlm,ik1,ik2,lc1,ls1,lc2,ls2,job0,job1
  double precision :: qdotr,pi,tpi,y0,ta2,x,y,a2,g0fac,xx1,xx2,x1,x2, &
       y2,p1(3),sp,cosp,sinp

  ! --- Setup ---
  if (kmax < 0 .OR. lmax < 0 .OR. rsm == 0d0) return
  job0 = mod(job,10)
  job1 = mod(job/10,10)
  nlm = (lmax+1)**2
  pi  = 4*datan(1d0)
  tpi = 2*pi
  y0  = 1/dsqrt(4*pi)
  a2  = 1/rsm**2
  ta2 = 2*a2
  do    ilm = 1, nlm
     do    k = 0, kmax
        gkl(1,k,ilm) = 0d0
        gkl(2,k,ilm) = 0d0
     enddo
  enddo
  ! ... Shorten connecting vector; need to adjust phase later
  if (job1 == 1 .OR. job1 == 3) then
     call shortn(p,p1,dlv,nkd)
  else
     call dcopy(3,p,1,p1,1)
  endif

  ! --- Put ylm in yl and alat**2*(p-dlv)**2 in wk(1) ---
  if (job0 == 0) then
     do  20  ir = 1, nkd
        wk(ir,2) = alat*(p1(1)-dlv(1,ir))
        wk(ir,3) = alat*(p1(2)-dlv(2,ir))
        wk(ir,4) = alat*(p1(3)-dlv(3,ir))
20   enddo
     call ropyln(nkd,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     ! ...   cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
     do  22  ir = 1, nkd
        qdotr = 2*pi*(q(1)*dlv(1,ir)+ q(2)*dlv(2,ir)+ q(3)*dlv(3,ir))
        wk(ir,3) = dcos(qdotr)
        wk(ir,4) = dsin(qdotr)
        wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
22   enddo
  elseif (job0 == 1) then
     do  24  ir = 1, nkd
        wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
24   enddo
  endif

  lc1 = 5
  ls1 = 6
  lc2 = 7
  ls2 = 8
  ik1 = 9
  ik2 = 10+lmax

  ! --- Outer loop over k (in blocks of 2), and over l ---
  do  301  k = 0, kmax, 2
     do  30  l = 0, lmax
        g0fac = 1/rsm*ta2**(l+1)/pi * dexp(e*rsm*rsm/4)

        !   ... Make radial part of the G_kl(1..nkd) for k= 0, 1
        if (k == 0) then
           do  32  ir = 1, nkd
              xx1 = g0fac*wk(ir,2)
              xx2 = (ta2*wk(ir,1)-3-2*l)* ta2 * xx1
              wk(ir,ik1+l) = xx1
              wk(ir,ik2+l) = xx2
              wk(ir,lc1) = wk(ir,3)*xx1
              wk(ir,ls1) = wk(ir,4)*xx1
              wk(ir,lc2) = wk(ir,3)*xx2
              wk(ir,ls2) = wk(ir,4)*xx2
32         enddo
           !   ... Make radial part of the G_kl(1..nkd) for k, k+1 from k-1, k-2
           !       and cos(q.dlv) * G_kl and sin(q.dlv) * G_kl
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
34         enddo
        endif

        !   ... For each point, add G_kl Y_L exp(i q.dlv) into Bloch G_kL
        ilm = l*l
        if (k < kmax) then
           do  36  m = -l, l
              ilm = ilm+1
              do  38  ir = nkd, 1, -1
                 gkl(1,k,ilm) = gkl(1,k,ilm) + wk(ir,lc1)*yl(ir,ilm)
                 gkl(2,k,ilm) = gkl(2,k,ilm) + wk(ir,ls1)*yl(ir,ilm)
                 gkl(1,k+1,ilm) = gkl(1,k+1,ilm) + wk(ir,lc2)*yl(ir,ilm)
                 gkl(2,k+1,ilm) = gkl(2,k+1,ilm) + wk(ir,ls2)*yl(ir,ilm)
38            enddo
36         enddo
        else
           do  46  m = -l, l
              ilm = ilm+1
              do  48  ir = nkd, 1, -1
                 gkl(1,k,ilm) = gkl(1,k,ilm) + wk(ir,lc1)*yl(ir,ilm)
                 gkl(2,k,ilm) = gkl(2,k,ilm) + wk(ir,ls1)*yl(ir,ilm)
48            enddo
46         enddo
        endif
30   enddo
301 enddo

  ! ... Put in phase to undo shortening, or different phase convention
  sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  if (job1 >= 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
  if (sp /= 0d0) then
     cosp = dcos(sp)
     sinp = dsin(sp)
     do    ilm = 1, nlm
        do    k   = 0, kmax
           x1 = gkl(1,k,ilm)
           x2 = gkl(2,k,ilm)
           gkl(1,k,ilm) = x1*cosp - x2*sinp
           gkl(2,k,ilm) = x2*cosp + x1*sinp
        enddo
     enddo
  endif

end subroutine gklq

