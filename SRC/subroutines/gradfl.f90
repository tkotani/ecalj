subroutine gradfl(lmax,nd,nr,np,ir0,ir1,lgg,lx,nn,ri,yl,gyl,fl,  gp,ggp)
  !- Gradient, Laplacian of function point-wise through sphere from YL expansion
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lmax  :density is expanded to l-cutoff lmax
  !i   nd    :dimensions yl,gyl
  !i   nr    :number of radial mesh points
  !i   np    :number of angular mesh points
  !i   ir0   :gp and gpp are made for radial points between ir0,ir1
  !i   ir1   :gp and gpp are made for radial points between ir0,ir1
  !i   lgg   :if zero, make gradient gp only; gpp not addressed.
  !i   lx    :(ones digit) if 1, fl scaled by r**2
  !i         :(tens digit): extrapolate 1st point (ir0=1) from others
  !i         :(100  digit): rational function interpolation for radial deriv
  !i   nn    :nn: number of points used to differentite radial f
  !i   ri    :vector of radial mesh points
  !i   yl    :Spherical harmonics for L=0:(lmax+1)^2 at each angular mesh point
  !i         :Generate with a call to ropyln
  !i   gyl   :Gradient of YL.
  !i         :Generate with a call to ropylg
  !i   fl    :function to be differentiated, on the combined radial
  !i         :and angular mesh
  !o Outputs
  !i   gp    :gradient of fl, on the combined radial and angular mesh,
  !i         :x,y,z components
  !i   ggp   :Laplacian of fl, on the combined radial and angular mesh
  !l Local variables
  !l   gf    :Work array (used for radial derivative of fl)
  !l   ggf   :Work array (used for 2nd radial derivative of fl)
  !r Remarks
  !r
  !u Updates
  !u   02 Apr 09 Made gf,ggf local; fixed bug for 10s digit lx case
  ! ----------------------------------------------------------------------
  implicit none
  integer :: np,nr,nd,lx,nn,lmax,ir0,ir1
  double precision :: fl(nr,1),gp(ir0:ir1,np,3),ggp(ir0:ir1,np), ri(nr),yl(np,nd),gyl(np,nd,3)
  integer :: i0,ilm,ip,ir,j,l,m,l2,lerr,lgg,iprint,jx,nx
  double precision :: xx,cy1,tol,egf0,gf(nr),ggf(nr)
  logical :: lrat
  parameter (tol=1d-12)
  if (ir0 < 1) call rx('gradfl: illegal value of ir0')
  cy1 = dsqrt(3/(16*datan(1d0)))
  l2 = mod(lx,100)
  lrat = lx/100 .ne. 0
  i0 = 1
  nx = nr
  if (l2/10 /= 0) then
     i0 = 2
     nx = nr-1
  endif
  ! --- Contribution (grad fl) r^-l yl ---
  call dpzero(gp,    (ir1-ir0+1)*np*3)
  if (lgg /= 0) call dpzero(ggp,   (ir1-ir0+1)*np)
  ilm = 0
  do  201  l = 0, lmax
     do  20  m = -l, l
        ilm = ilm+1
        if (mod(l2,10) == 0) then
           call poldvm(ri(i0),fl(i0,ilm),nx,nn,lrat,tol,lerr,gf(i0))
           if (lerr /= 0) goto 99
        else
           do   ir = i0, nr
              ggf(ir) = fl(ir,ilm)/ri(ir)**2
           enddo
           call poldvm(ri(i0),ggf(i0),nx,nn,lrat,tol,lerr,gf(i0))
           if (lerr /= 0) goto 99
        endif
        !     Extrapolate gf to first point
        if (l2/10 /= 0) then
           jx = 1
           call polint(ri(2),gf(2),nx,nn,ri,0d0,0,jx,gf,egf0)
           lerr = 1
           if (iprint() >= 50 .AND. &! & takao. too noizy.
              dabs(egf0) .gt. 1d-3*max(dabs(gf(1)),dabs(gf(2)))) then
!         call info5(40,0,0,' gradfl (warning): uncertainty in grad'// &
!                 ' f(r=0,L=%i):  f=%;3g  est err= %;3g',ilm,gf(1),egf0,0,0)
!   print *,'TAKAO: this warning is probably not a problem. If you like, plot ri.vs.gf as in gradfl.'
           endif
        endif
        do    ip = 1, np
           do    ir = ir0, ir1
              gp(ir,ip,1) = gp(ir,ip,1) + gf(ir)*yl(ip,ilm)
           enddo
        enddo
        ! --- Laplacian: (nabla fl) Yl + fl (nabla Yl) ---
        if (lgg /= 0) then
           call poldvm(ri(i0),gf(i0),nx,nn,lrat,tol,lerr,ggf(i0))
           if (lerr /= 0) goto 99
           if (mod(l2,10) == 0) then
              do    ir = i0, nr
                 xx = 1/ri(ir)
                 ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx*xx
              enddo
           else
              do    ir = i0, nr
                 xx = 1/ri(ir)
                 ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx**4
              enddo
           endif
           if (i0 == 2) ggf(1)= (ri(3)*ggf(2)-ri(2)*ggf(3))/(ri(3)-ri(2))
           do  ip = 1, np
              do  ir = ir0, ir1
                 ggp(ir,ip) = ggp(ir,ip) + ggf(ir)*yl(ip,ilm)
              enddo
           enddo
        endif
20   enddo
201 enddo
  ! ... Split grad r- into x,y,z- components
  do    j = 3, 1, -1
     do    ip = 1, np
        xx = yl(ip,j)/cy1
        if (j == 1) xx = yl(ip,4)/cy1
        do    ir = ir0, ir1
           gp(ir,ip,j) = xx*gp(ir,ip,1)
        enddo
     enddo
  enddo
  ! --- Contribution fl(r) grad r^-l yl (use gf as work array) ---
  ilm = 0
  do  101  l = 0, lmax
     do  10  m = -l, l
        ilm = ilm+1
        ! ...   Factor 1/r from grad Yl
        if (mod(l2,10) == 0) then
           do   ir = max(i0,ir0), nr
              gf(ir) = fl(ir,ilm)/ri(ir)
           enddo
        else
           do    ir = max(i0,ir0), nr
              gf(ir) = fl(ir,ilm)/ri(ir)**3
           enddo
        endif
        if (i0 > ir0) gf(1) = (ri(3)*gf(2)-ri(2)*gf(3))/(ri(3)-ri(2))
        do    j = 1, 3
           do    ip = 1, np
              xx = gyl(ip,ilm,j)
              do    ir = ir0, ir1
                 gp(ir,ip,j) = gp(ir,ip,j) + gf(ir)*xx
              enddo
           enddo
        enddo
10   enddo
101 enddo
  return
  ! --- Error handling ---
99 print *, 'gradfl: stopping at ilm=',ilm,'  point', lerr
  call rx('gradfl: can''t diff radial function')
end subroutine gradfl
