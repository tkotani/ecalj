module m_lattic
  public m_lattic_init,Setopos
  real(8), allocatable,protected,public ::  rv_a_odlv (:)
  real(8), allocatable,protected,public ::  rv_a_oqlv (:)
  real(8), protected,public:: lat_plat(3,3),lat_qlat(3,3),lat_awald,lat_vol
  real(8), allocatable,protected,public :: rv_a_opos(:,:)
  integer,protected,public:: lat_nkd,lat_nkq
  logical,protected,public:: lattic_init=.false.
  ! real(8),protected:: lat_dist(3,3) !unused because no deformation of cell allowed currently.
  private
contains
  subroutine Setopos(posin) !called from lmfp to revise atomic position by reading rst file.
    use m_lmfinit,only: ctrl_nbas
    integer:: i
    real(8):: posin(:,:)
    do i=1,ctrl_nbas
       rv_a_opos(:,i)= posin(:,i) !v_ssite(i)%pos
    enddo
  end subroutine Setopos
  subroutine m_lattic_init() !no shear now  ldist=0
    use m_lmfinit,only:ctrl_nbas,lat_alat,lat_as,lat_tol,lat_rpad,lat_nkdmx,lat_nkqmx,lat_gam,&
        lat_platin, poss=>pos
    ! Sets up the real and rmeciprocal space lattice vectors
    ! ----------------------------------------------------------------------
    !u Updates
    !u   2 Mar 04 Pass rpad to lattc
    !u   5 Jun 01 (ATP) Now calls lattc after lattice transformation
    !u  19 Apr 00 Fixed rotations; new argument list
    ! ----------------------------------------------------------------------
    implicit none
    integer::  lmxst , nkd , nkdmx , nkq , nkqmx , nbas,i_data_size,ib
    real(8),allocatable:: rv_a_tmp(:)
    real(8):: alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol, &
         xx1,xx2,dotprd,pi,rpad, plat0(3,3),plat(3,3),qlat(3,3) 
    equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
    call tcn('m_lattic_init')
    lattic_init=.true.
    alat=lat_alat
    awald0=lat_as
    tol=lat_tol
    rpad=lat_rpad
    nkdmx=lat_nkdmx
    nkqmx=lat_nkqmx
    gam = lat_gam
    alat = lat_alat
    plat0=lat_platin
    nbas=ctrl_nbas
    allocate(rv_a_opos(3,nbas))
    do ib=1,nbas
       rv_a_opos(:,ib)= poss(:,ib) !v_ssite(ib)%pos
    enddo
    if (abs(gt-1d0) > 1d-10) then
       call rdistn ( rv_a_opos , rv_a_opos , nbas , gx , gy , gz , gt )
       !        call rdistn ( rv_a_opos , rv_a_opos , nbaspp , gx , gy , gz , gt )
       !      elseif (ldist .ne. 0) then
       !        call lattdf ( ldist , dist , plat0 , nbaspp , rv_a_opos , 0 ,  0d0 )
       !      else
       !        dist=0d0
       !        dist(1,1) = 1
       !        dist(2,2) = 1
       !        dist(3,3) = 1
    endif
    allocate(rv_a_odlv(abs(3*nkdmx)))
    allocate(rv_a_oqlv(abs(3*nkqmx)))
    lmxst = 6
    call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx , &
         gy , gz , gt , plat , qlat , lmxst , vol , awald , rv_a_odlv &
         , nkd , rv_a_oqlv , nkq , nkdmx , nkqmx )
    lat_vol  =vol
    !      lat_plat0=plat0
    lat_plat =plat
    lat_qlat =qlat
    !! reduce size. necessary?
    i_data_size=size(rv_a_oqlv)
    allocate(rv_a_tmp(i_data_size))
    rv_a_tmp=rv_a_oqlv
    deallocate(rv_a_oqlv)
    i_data_size=min(i_data_size,3*nkq)
    allocate(rv_a_oqlv(3*nkq))
    rv_a_oqlv(:i_data_size)=rv_a_tmp(:i_data_size)
    deallocate(rv_a_tmp)
    !! reduce size. necessary?
    i_data_size=size(rv_a_odlv)
    allocate(rv_a_tmp(i_data_size))
    rv_a_tmp=rv_a_odlv
    deallocate(rv_a_odlv)
    i_data_size=min(i_data_size,3*nkd)
    allocate(rv_a_odlv(3*nkd))
    rv_a_odlv(:i_data_size)=rv_a_tmp(:i_data_size)
    deallocate(rv_a_tmp)
    lat_awald=awald
    lat_nkd=nkd
    lat_nkq=nkq
    !      lat_dist=dist
    call tcx('m_lattic_init')
  end subroutine m_lattic_init
  subroutine lattdf(ldist,defgrd,plat,nbas,bas,ngen,gen) !- Rotates or deforms lattice
    ! ----------------------------------------------------------------
    !i Inputs
    !i   ldist: 0, no deformations.  For abs(ldist):
    !i          1: defgrd holds rot about spec'd angle
    !i          2, lattice deformed with a general linear transformation
    !i          3, lattice deformed by a shear.
    !i          SIGN ldist <0 => suppress rotation of plat
    !i   defgrd:transformation matrix, whose form is specified by ldist
    !i   nbas  :size of basis
    !i   ngen  :number of generators of the symmetry group
    !o Outputs
    ! o  plat  :primitive lattice vectors, in units of alat
    ! o        :On output, plat is transformed by defgrd
    ! o  bas   :basis vectors, in units of alat
    ! o        :On output, bas is transformed by defgrd
    ! o  gen   :On input, generators of the symmetry group
    ! o        :On output, generators are transformed by defgrd
    ! o        :Note: gen can only be transformed by a rotation
    ! o  defgrd:Output defgrd is actual transformation matrix.
    !u Updates
    !u   19 Mar 06 Blend Voigt strains into other types
    ! ----------------------------------------------------------------
    implicit none
    integer :: ldist,nbas,ngen
    double precision :: defgrd(3,3), plat(3,3), bas(3,1), gen(3,3,1)
    double precision :: work(3,3),rinv(3,3),det,gold(3,3)
    integer :: ipr,i,j,ib
    real(8) ,allocatable :: bwk_rv(:)
    if (ldist == 0) return
    call getpr(ipr)
    allocate(bwk_rv(3*nbas))
    call dcopy(9,defgrd,1,work,1)
    if (iabs(ldist) == 1) call makrot(work,defgrd)
    if (iabs(ldist) == 3) then
       det = defgrd(1,3)
       if ( det == 0 ) then
          if (allocated(bwk_rv)) deallocate(bwk_rv)
          return
       endif
       call shear(0,plat,bas,det,defgrd,defgrd)
    endif
    if (ipr >= 30) then
       print 333, ldist, ((defgrd(i,j), j=1,3), i=1,3)
333    format(/' LATTDF:  deformation matrix for mode',i3,':'/ &
            (3f12.7))
    endif
    ! ... Rotate or shear plat
    if (ldist >= 1 .AND. ldist <= 3) then
       !call dcopy(9,plat,1,work,1)
       call dinv33(defgrd,0,rinv,det)
       plat=matmul(defgrd,plat) !call dmpy(defgrd,3,1,work,3,1,plat,3,1,3,3,3)
       if (ipr >= 30) then
          print 334, ((work(i,j), i=1,3), (plat(i,j), i=1,3),j=1,3)
334       format(10x,' Lattice vectors:',25x,'Transformed to:'/ &
               (3f12.7,2x,3f12.7))
       endif
    endif
    if (iabs(ldist) >= 1 .AND. iabs(ldist) <= 3) then
       if (ipr >= 30 .AND. nbas > 0) &
            print '(10x,''  Basis vectors:'',25x,''Transformed to:'')'
       do  10  ib = 1, nbas
          !call dcopy(3,bas(1,ib),1,work,1)
          bas(:,ib)=matmul(defgrd,bas(:,ib)) !call dmpy(defgrd,3,1,work,3,1,bas(1,ib),3,1,3,1,3)
          if (ipr >= 30) print '(3f12.7,2x,3f12.7)', &
               (work(i,1), i=1,3), (bas(i,ib), i=1,3)
10     enddo
       if (ipr >= 30 .AND. ngen > 0) &
            print '(15x,''Group ops:'',26x,''Rotated to:'')'
       call dinv33(defgrd,0,rinv,det)
       do  20  ib = 1, ngen
          gold=gen(:,:,ib) !call dcopy(9,gen(1,1,ib),1,gold,1)
          work=matmul(defgrd,gen(:,:,ib)) !call dmpy(defgrd,3,1,gen(1,1,ib),3,1,work,3,1,3,3,3)
          gen(:,:,ib)=matmul(work,rinv) !call dmpy(work,3,1,rinv,3,1,gen(1,1,ib),3,1,3,3,3)
          if (ipr >= 30) print '(/(3f12.7,2x,3f12.7))', &
               ((gold(j,i), i=1,3), (gen(j,i,ib), i=1,3),j=1,3)
20     enddo
    endif
    if (allocated(bwk_rv)) deallocate(bwk_rv)
  end subroutine lattdf
  subroutine shear(nbas,plat,tau,alpha,eps,s)!- Apply a pure strain to lattice and basis vectors
    use m_lmfinit,only: stdo
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis.  nbas=0 => plat, tau not sheared
    !i   plat  :basis vectors, in units of alat
    !i   tau   :position vectors, in units of alat
    !i   alpha :tight-binding screening parameters
    !o   eps   :alpha*eps(1..6):  Voigt tensor strains
    !o Outputs
    !o   s     :deformation matrix
    !o         :s and eps can occupy the same address space
    !l Local variables
    !l         :
    !r Remarks
    !u Updates
    !u   19 Mar 06 Added output s
    ! ----------------------------------------------------------------------
    integer :: nbas
    double precision :: plat(3,3),tau(3,nbas),eps(6),alpha
    integer :: i,j,ind(2,21),n,iprint !,nbmx
    !parameter (nbmx = 256)
    double precision :: s(3,3),b(3,3),x(21),t(3,nbas),det,e(6)
    if (dabs(alpha) < 1d-8) return
    !if (nbas > nbmx) call fexit(-1,9,'Increase nbmx in SHEAR',0)
    e=alpha*eps !call dcopy(6,eps,1,e,1)  !call dscal(6,alpha,e,1)
    call xxxes(e,s,det)
    ! ... Shear plat, tau
    if (nbas > 0) then
       if (iprint() >= 40) then
          write (stdo,10) plat
          write (stdo,20)
          write (stdo,30) ((tau(i,j),i=1,3),j=1,nbas)
          write (stdo,40) s
       endif
       plat=matmul(s,plat) !call dmpy(s,3,1,plat,3,1,b,3,1,3,3,3)  !call dcopy(9,b,1,plat,1)
       tau =matmul(s,tau) !call dmpy(s,3,1,tau,3,1,t,3,1,3,nbas,3)   !call dcopy(3*nbas,t,1,tau,1)
       if (iprint() >= 40) then
          write (stdo,10) plat
          write (stdo,20)
          write (stdo,30) ((tau(i,j),i=1,3),j=1,nbas)
       endif
    endif
    call xxxse(s,e)
    call dscal(6,1d0/alpha,e,1)
    n = 0
    do  1  i = 1, 6
       call xxxadd(i,i,n,e,ind,x)
1   enddo
    do   i = 1, 5
       do    j = i+1, 6
          call xxxadd(i,j,n,e,ind,x)
       enddo
    enddo
    if (n == 0 .OR. iprint() < 30) return
    write (stdo,50) alpha,det-1
    write (stdo,60) ((ind(i,j),i=1,2),j=1,n)
    write (stdo,70) (x(i),i=1,n)
10  format (' SHEAR: Lattice vectors:'/3(8x,3f10.6/))
20  format ('        Basis atoms:')
30  format (8x,3f10.6)
40  format ('       Lattice and basis sheared by'/3(8x,3f10.6/))
50  format (/ &
         1x,'SHEAR: distortion amplitude =',f9.6, &
         '  Vol. dilatation =',f9.6/ &
         8x,'The second derivative E/vol w.r.t alpha = ', &
         'W'''' = sum_ij x_ij c_ij'/ &
         8x,'has coefficients x_ij to elastic constants c_ij as shown:')
60  format (6(5x,2i1,5x))
70  format (6f12.6)
  end subroutine shear
  subroutine xxxes(e,s,det)  ! Make deformation tensor
    double precision :: e(6),s(3,3),det
    s(1,1) = 1 + e(1)
    s(2,2) = 1 + e(2)
    s(3,3) = 1 + e(3)
    s(1,2) = e(6)
    s(2,1) = e(6)
    s(1,3) = e(5)
    s(3,1) = e(5)
    s(2,3) = e(4)
    s(3,2) = e(4)
    det=s(1,1)*s(2,2)*s(3,3)+s(1,2)*s(2,3)*s(3,1) &
         +s(1,3)*s(2,1)*s(3,2)-s(1,3)*s(2,2)*s(3,1) &
         -s(1,2)*s(2,1)*s(3,3)-s(1,1)*s(2,3)*s(3,2)
  end subroutine xxxes
  subroutine xxxse(s,e) ! Make engineering strains
    double precision :: s(3,3),e(6)
    integer :: i
    do  1  i = 1, 3
       e(i) = s(i,i) - 1d0
1   enddo
    e(4) = 2*s(2,3)
    e(5) = 2*s(1,3)
    e(6) = 2*s(1,2)
  end subroutine xxxse
  subroutine xxxadd(i,j,n,e,ind,x)! Add to list of non-zero elastic constant coefficients
    implicit none
    integer :: i,j,n,ind(2,21)
    double precision :: e(6),x(21),xx
    xx = e(i)*e(j)
    if (i /= j) xx = 2*xx
    if (dabs(xx) > 1d-8) then
       n = n+1
       ind(1,n) = i
       ind(2,n) = j
       x(n) = xx
    endif
  end subroutine xxxadd
  subroutine makrot(rot,r)
    !- Make rotation vector and angle into deformation gradient tensor
    !     implicit none
    double precision :: rot(4),r(3,3)
    double precision :: n(3),cost,sint,rnorm,dcos,dsin,dsqrt
    integer :: i,j,k,e,mod
    cost = dcos(rot(4))
    sint = dsin(rot(4))
    rnorm = dsqrt(rot(1)**2 + rot(2)**2 + rot(3)**2)
    do  1  k = 1, 3
       n(k) = rot(k) / rnorm
1   enddo
    do  21  i = 1, 3
       do  2  j = 1, 3
          r(i,j) = (1 - cost)*n(i)*n(j)
          if (i == j) then
             r(i,j) = r(i,j) + cost
          else
             k = 3 - mod(i+j,3)
             if (i < k) then
                e = 1
                if (k-i == 2) e = -1
             else
                e = -1
                if (i-k == 2) e = 1
             endif
             r(i,j) = r(i,j) + sint*n(k)*e
          endif
2      enddo
21  enddo
  end subroutine makrot
end module m_lattic

subroutine lattc(as,tol,rpad,alat,alat0,platin,g1,g2,g3,gt,plat,qlat,lmax,vol,awald,dlat,nkd,glat,nkg,nkdmx,nkgmx)  ! Sets up the real and reciprocal space lattice vectors for Ewald
  use m_ftox
  use m_lgunit,only:stdo
  ! ----------------------------------------------------------------
  !i Inputs
  !i   as    :dimensionless Ewald parameter (2 is suggested).
  !i         :Ewald parameter awald scales with the lattice as
  !i         :as/(vol)**(1/3)
  !i   tol   :tolerance for ewald sums
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   platin:primitive lattice translation vectors, in units of alat
  !i   g1,g2,g3:x,y,z distortions of platin
  !i   gt:    multiplier of g1,g2,g3.  gt=1 => no distortion
  !i   lmax:  Ewald sums will be taken to maximum L.
  !i   nkdmx: maximum number of direct lattice vectors
  !i   nkgmx: maximum number of reciprocal lattice vectors
  !o Outputs
  !o    plat:   distorted lattice vectors
  !o    qlat:   distorted reciprocal vectors
  !o     vol:   cell volume
  !o   awald:   ewald parameter
  !o   dlat,nkd: direct lattice vectors and number
  !o   glat,nkg: reciprocal lattice vectors and number
  !r Remarks
  !r   awald is in (atomic units)^-1.
  !r   Tolerance tol is the estimated error for a lattice of unit volume.
  !r   plat(*,k) holds lattice vector k
  !r   plat can point to same address as platin; then platin => plat
  !r
  !r   Local variables r0 and q0 are the ranges for
  !r   real-space and reciprocal lattice vectors for a unit lattice.
  !u Updates
  !u   2 Mar 04 New rpad: truncate radius of lattice vectors to rpad*rmax
  !u            when list has to be padded in order to include at
  !r            least one lattice vector.
  ! ----------------------------------------------------------------
  implicit none
  integer :: lmax,nkd,nkg,nkdmx,nkgmx
  double precision :: as,tol,alat,g1,g2,g3,gt,vol,awald,alat0,rpad, &
       glat(3,nkgmx),dlat(3,nkdmx),platin(3,3),plat(3,3),qlat(3,3)
  integer :: k,iprint,m,i1mach,modeg(3),isw
  double precision :: qlat0(3,3),vol0,plat0(3,3),radd,qadd
  double precision :: qdist0,a0,rdist0,tol1,r0,q0,one(3,3),oned(3,3)
  integer:: ifile_handle,ifp,i
  real(8):: rxx,qxx
  plat0=platin !call dcopy(9,platin,1,plat0,1)
  call dinv33(plat0,1,qlat0,vol0)
  vol0 = dabs(vol0)
  call rdistn(plat0,plat,3,g1,g2,g3,gt)
  call dinv33(plat,1,qlat,vol)
  vol = dabs(vol)*(alat**3)
  if (iprint() < 20) goto 20
  write(stdo,351)
351 format(/t17,'Plat',t55,'Qlat')
  write (stdo,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
  open(newunit=ifp,file='PlatQlat.chk')
  write (ifp,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
  write (ifp,"('             PLAT              and         QLAT    ')")
  close(ifp)
350 format(3f11.6,5x,3f11.6)
  if (dabs(gt-1d0) > 1.d-5) then
     write(stdo,ftox)' Distorted with gx,y,z=',ftof(g1),ftof(g2),ftof(g3),'gt=:',ftof(gt)
     write(stdo,350) ((plat(m,k),m=1,3),(qlat(m,k),m=1,3),k=1,3)
     one=0
     one(1,1) = 1
     one(2,2) = 1
     one(3,3) = 1
     call rdistn(one,oned,3,g1,g2,g3,gt)
     if (iprint() > 40) write(stdo,352) ((oned(m,k),m=1,3), k=1,3)
352  format(t14,'shear matrix'/(3f11.7))
  endif
  write(stdo,ftox)'  Cell vol= ',ftof(vol)
  if((dabs(vol-vol0*(alat**3)) > 1d-9))write(stdo,ftox)'(undistorted vol=',ftof(vol0*(alat**3))
20 continue
  ! --- Set up real and reciprocal vectors ---
  ! The errors are estimated making a continuum approximation to a
  ! discrete set of lattice sums.  Adding .7, slightly more than
  ! half the average spacing makes the continuum approximation err
  ! on the safe side.
  rdist0 = vol0**(1d0/3d0)
  qdist0 = 1d0/rdist0
  radd = 0.7d0*rdist0 
  qadd = 0.7d0*qdist0 
  a0 = as/rdist0
  awald = a0/alat
  tol1 = tol*alat0**(lmax+1)
  call lctoff(a0,vol0,lmax,tol1,r0,q0)
  modeg(1) = 2
  modeg(2) = 2
  modeg(3) = 2
!  rxx=r0+radd
!  qxx=q0+qadd
  rxx= maxval([r0+radd,(sum(plat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12 for very anisotropic cases safer.
  qxx= maxval([q0+qadd,(sum(qlat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12
  call xlgen(plat0,rxx,rpad*(r0+radd),nkdmx,11,modeg,nkd,dlat)
  call xlgen(qlat0,qxx,rpad*(q0+qadd),nkgmx,11,modeg,nkg,glat)
  call rdistn(dlat,dlat,nkd,g1,g2,g3,gt)
  call qdistn(glat,glat,nkg,g1,g2,g3,gt)
  ! --- Printout ---
  if (iprint() < 30) goto 60
  write (stdo,340) as,tol,alat,awald
  write (stdo,342) r0+radd,nkd,q0+qadd,nkg
340 format(/' LATTC:  as=',f6.3,'   tol=',1p,e9.2, &
       '   alat=',0p,f8.5,'   awald=',f6.3)
342 format(9x,'r1=',f7.3,'   nkd=',i4,'      q1=',f7.3,'   nkg=',i4)
60 if(dabs(alat0/alat-1d0) > 0.04d0) call rx('lattc: alat and alat0 deviate by more than 4 %')
end subroutine lattc
subroutine lctoff(a0,v0,lmax,tol,r0,q0) !- makes limits r0,q0 for sums in real and recip space for a lattice
  !  with lattice constant 1.
  !u Updates
  !u   25 Jun 03 (Kino) bug fix in dimension of f and g
  !     implicit none
  integer :: lmax,i
  double precision :: a0,q0,r0,tol,v0
  double precision :: gq0,gq1,pi,q1,q2,r1,r2
  double precision :: f(0:lmax),g(0:lmax)
  pi = 4d0*datan(1d0)
  q1 = 0.001d0
  if (lmax > 2) q1 = dsqrt(.5d0*(lmax-2))*a0/pi
  gq1 = (2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
  if (tol > gq1) call info0(10,0,0,' lctoff (warning): tol gt gq1')
  q2 = 50d0
  q0 = 5d0
  do  33  i = 1, 25
     gq0 = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
     if(gq0 > tol) q1 = q0
     if(gq0 < tol) q2 = q0
     q0 = .5d0*(q1+q2)
33 enddo
  r1 = 0.1d0
  r2 = 50d0
  r0 = 5d0
  do  15  i = 1, 25
     call dlmtor(r0,a0,lmax,f,g)
     if(f(lmax) > tol) r1 = r0
     if(f(lmax) <= tol) r2 = r0
     r0 = .5d0*(r1+r2)
15 enddo
  !|    try = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
  !|    write(6,957) q0,try,r0,f(lmax)
  ! 957 format(' lcut: q0=',f12.6,'   try=',f12.6,'   r0=',f12.6,
  !|   .  '   f=',f12.6)
end subroutine lctoff
subroutine dlmtor(r,a,lmax,f,fbar)  !- Radial part of damped lmtos f and fbar, l=0 to lmax
  !     implicit none
  integer :: l,lmax
  double precision :: a,f(0:lmax),fbar(0:lmax),r
  double precision :: derfc,emz2,erfc0,erfc1,erfc2,fbsrpi, flm2,g,ta2r,z
  fbsrpi = 0.564189835d0
  z = a*r
  emz2 = dexp(-z*z)
  erfc0 = derfc(z)
  erfc1 = -z*erfc0 + fbsrpi*emz2
  erfc2 = -0.5d0*z*erfc1 + 0.25d0*erfc0
  f(0) = erfc0/r
  fbar(0) = -erfc2/(a*a*r)
  ta2r = 2d0*a*a*r
  g = 2d0*a*emz2*fbsrpi/r
  flm2 = fbsrpi*emz2/z - erfc0
  do  10  l = 1, lmax
     f(l) = ((l+l-1)/r)*f(l-1) + g
     fbar(l) = ((l+l-1)/r)*fbar(l-1) - flm2
     flm2 = f(l-1)
     g = g*ta2r
10 enddo
end subroutine dlmtor
subroutine qdist(q,q1,ux,uy,uz,gam) !  (from msm) new version of 03.10.89. u=(ux,uy,uz) gives direction,
  !  gam is multiplier in real space along u, volume is conserved.
  !     implicit none
  double precision :: q(3),q1(3),ux,uy,uz,gam
  double precision :: g,u2,a,xxx,yyy
  g = 1d0/gam
  u2 = ux*ux+uy*uy+uz*uz
  a = (ux*q(1)+uy*q(2)+uz*q(3))/u2
  xxx = 1d0/dsqrt(dabs(g))
  if (g < 0d0) xxx = 1
  yyy = a*(g-xxx)
  q1(1) = xxx*q(1)+yyy*ux
  q1(2) = xxx*q(2)+yyy*uy
  q1(3) = xxx*q(3)+yyy*uz
end subroutine qdist
subroutine rdist(v,v1,ux,uy,uz,g)
  !     implicit none
  double precision :: v(3),v1(3),ux,uy,uz,g,u2,a,xxx,yyy
  u2 = ux*ux+uy*uy+uz*uz
  a = (ux*v(1)+uy*v(2)+uz*v(3))/u2
  xxx = 1d0/dsqrt(dabs(g))
  if (g < 0d0) xxx=1d0
  yyy = a*(g-xxx)
  v1(1) = xxx*v(1)+yyy*ux
  v1(2) = xxx*v(2)+yyy*uy
  v1(3) = xxx*v(3)+yyy*uz
end subroutine rdist
subroutine rdistn(a1,a2,na,gx,gy,gz,gt)  !- distort na real-space vectors in array a1 into array a2
  !     implicit none
  integer :: na,ia
  double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
  do 10 ia = 1, na
     call rdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10 enddo
end subroutine rdistn
subroutine qdistn(a1,a2,na,gx,gy,gz,gt)  !- distort na q-space vectors in array a1 into array a2
  !     implicit none
  integer :: na,ia
  double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
  do 10 ia = 1,na
     call qdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10 enddo
end subroutine qdistn

