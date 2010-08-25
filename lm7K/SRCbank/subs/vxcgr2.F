      subroutine vxcgr2(nr,nsp,nrx,rofi,rp,
     . exc,vxc)
c      subroutine vxcgr2(nr,nsp,nrx,rofi,rp,
c     .grh,ggrh,agrh,grgagr,exc,vxc)
ctakao automatic array version
C- Gradient correction to vxc, exc for a mesh of points.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nrx   :leading dimension of the radial function arrays
Ci   rofi  :radial mesh points
Ci   rp    :density rho on a radial mesh
Ci         :the following work arrays are dimensioned (nrx,2)
Ci   grh   :work array : radial grad rho
Ci   ggrh  :work array : laplacian rho
Ci   agrh  :work array : abs(grh)
Ci   grgagr:work array : grad rho . grad abs grad rho
Co Outputs
Co   exc   :gradient contribution to energy added to exc
Co   vxc   :gradient contribution to potential added to vxc
Cl Local variables
Cr Remarks
Cr   
Cu Updates
Cu   18 Jun 04 Bug fix
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nsp,nrx
      double precision rp(nrx,nsp),grh(nrx,2),ggrh(nrx,2),
     .agrh(nrx,4),grgagr(nrx,3),exc(nrx),vxc(nrx,2),rofi(nr)
C ... Local parameters
      integer ir,i,lxcf,lxcg,nglob

      lxcg = mod(nglob('lxcf')/100,100)
C      integer iprint
C      if (iprint() .ge. 80) then
C        call prmr(20,rofi,rp(2,1),1)
C        call prmr(20,rofi,rp(2,2),1)
C      endif

C --- grad(rho), laplacian rho ---
      call radgrx(nr,nrx,nsp,rofi,rp,grh)
      call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
      do  20  i  = 1, nsp
        do  24  ir = 2, nr
   24   ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
        ggrh(1,i) =(rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))

C --- grad rho . grad abs grad rho ---
        do  26  ir = 1, nr
   26   agrh(ir,i) = dabs(grh(ir,i))
        call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
        do  28  ir = 1, nr
   28   grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
   20 continue

C --- Extra terms g(n), g(n+).g(n-), g(n).g(abs(g(n))) if spin pol ---
      if (nsp .eq. 2) then
        do  32  ir = 1, nr
   32   agrh(ir,3) = dabs(grh(ir,1)+grh(ir,2))
        call radgrx(nr,nrx,1,rofi,agrh(1,3),grgagr(1,3))
        do  34  ir = 1, nr
   34   grgagr(ir,3) = (grh(ir,1)+grh(ir,2))*grgagr(ir,3)
        do  36  ir = 1, nr
   36   agrh(ir,4) = grh(ir,1)*grh(ir,2)
      endif

C --- Gradient term for all points ---
      if (lxcg .ge. 3) then
        lxcf = mod(nglob('lxcf'),100)
        if (lxcf .ne. 3 .and. lxcf .ne. 4) call
     .  rx('vxcgf2: inconsistent use of local and GGA functionals')
        call vxcgga(lxcg,nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp),
     .  ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4),
     .  grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp),
     .  vxc(1,1),vxc(1,nsp),exc)
      elseif (lxcg .eq. 2) then
        call rx('PW91 no longer implemented')
      else
        call vxnloc(nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp),
     .  ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4),
     .  grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp),
     .  vxc(1,1),vxc(1,nsp),exc)
      endif
      do  66  i = 1, nsp
   66 vxc(1,i) = (vxc(2,i)*rofi(3)-vxc(3,i)*rofi(2))/(rofi(3)-rofi(2))

      end
      subroutine radgrx(nr,nrx,nsp,ri,f,gf)
C     implicit none
      integer nr,nrx,nsp,nn,i,iprint,jx
      double precision ri(nr),f(nrx,nsp),gf(nrx,nsp),tol,egf0
      logical lerr
      parameter (tol=1d-12,nn=6)

      do  10  i = 1, nsp
C        call prmr(nr,ri,gf,1)
        call poldvm(ri(2),f(2,i),nr-1,nn,.false.,tol,lerr,gf(2,i))
        jx = 1
        call polint(ri(2),gf(2,i),nr-1,nn,ri,0d0,0,jx,gf(1,i),egf0)
        if (iprint() .ge. 50 .and. dabs(egf0/gf(1,i)) .gt. 1d-2)
     .  print 345, gf(1,i), egf0/gf(1,i)*100
  345   format(' radgrx (warning): expect error in gradient at origin:',
     .  'f=',1pe10.3,' est err=',0pf7.1,'%')
   10 continue
      end

