subroutine grfmsh(isw,alat,ng,gv,kv,k1,k2,k3,n1,n2,n3,fn, &
     fgrd,flap)
  !- Gradient and Laplacian of a function tabulated on a uniform mesh
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isw   :1s digit
  !i         : 0 input function fn is real
  !i         : 1 input function fn is complex
  !i         :10s digit
  !i         : 0 Input fn is in real space
  !i         : 1 Input fn is in reciprocal space
  !i         :100s digit handles output
  !i         : 0 No output
  !i         : 1 Output fgrd, in reciprocal space, no flap
  !i         : 2 Output fgrd, in real space, no flap
  !i         : 3 Output flap, in reciprocal space, no fgrd
  !i         : 4 Output flap, in real space, no fgrd
  !i         : 5 Output fgrd and flap, in reciprocal space
  !i         : 6 Output fgrd and flap, in real space
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   vol   :cell volume
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   k1..k3:dimensions of fn,fgrd,flap
  !i   n1..n3:size of uniform mesh for which fn,fgrd,flap are tabulated
  !i   fn    :function on uniform mesh, either in real or recip space,
  !i         :depending on isw
  !o Outputs
  !o   fgrd  :gradient of fn either in real or recip space, depending on isw
  !o   flap  :laplacian of fn either in real or recip space, depending on isw
  !l Local variables
  !l   gi    :i * G (for gradient)
  !l   g2    :(i * G)^2 (for Laplacian)
  !l   fg    :gradient of function, G space
  !l   fl    :Laplacian of function, G space
  !r Remarks
  !r Sample test:
  !r   cmfft -f9f15.10 -248 rho | mc -f9f15.10  . igx -xe | cmfft -f9f15.10 -i -248 .
  !u Updates
  !u   08 Apr 09 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ng,isw,k1,k2,k3,n1,n2,n3
  integer :: kv(ng,3)
  double precision :: alat,gv(ng,3)
  double complex fn(k1,k2,k3),fgrd(k1,k2,k3,3),flap(k1,k2,k3)
  ! ... Local parameters
  integer :: i,isw0,isw1,isw2
  double precision :: pi,tpiba,g2
  double complex tpibai,gi(3)
  complex(8),allocatable:: fg(:),fgg(:,:),fg2(:)

  call tcn('grfmsh')
  pi   = 4d0*datan(1d0)
  tpiba=2*pi/alat
  isw0 = mod(isw,10)
  isw1 = mod(isw/10,10)
  isw2 = mod(isw/100,10)
  if (isw0 /= 1) &
       call rx('grfmsh: gradient of real function not implemented')

  ! ... FT of smooth function to reciprocal space
  if (isw1 == 0) then
     !       call zprm3('fn(r)',0,fn,n1,n2,n3)
     call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)
     !       call zprm3('fn(G)',0,fn,n1,n2,n3)
  endif

  ! ... Gather function G coefficients
  allocate(fg(ng))
  call gvgetf(ng,1,kv,k1,k2,k3,fn,fg)

  ! ... Restore given function to real space
  if (isw1 == 0) then
     call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,1)
  endif

  ! ... Make iG * f(G)  and (iG)^2 * f(G)
  fg(1) = 0
  tpibai = dcmplx(0d0,1d0)*tpiba
  allocate(fgg(ng,3),fg2(ng))
  do  i = 1, ng

     gi(1) = tpibai*(gv(i,1))
     gi(2) = tpibai*(gv(i,2))
     gi(3) = tpibai*(gv(i,3))
     g2 = -tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)

     fgg(i,1) = gi(1) * fg(i)
     fgg(i,2) = gi(2) * fg(i)
     fgg(i,3) = gi(3) * fg(i)
     !       fgg(i,1) = gi(1)
     !       fgg(i,2) = gi(2)
     !       fgg(i,3) = gi(3)
     fg2(i)   = g2    * fg(i)

  enddo

  ! ... Scatter gradient, laplacian into respective arrays
  if (isw2 == 1 .OR. isw2 == 2 .OR. isw2 == 5 .OR. isw2 == 6) then
     call gvputf(ng,1,kv,k1,k2,k3,fgg(1,1),fgrd(1,1,1,1))
     call gvputf(ng,1,kv,k1,k2,k3,fgg(1,2),fgrd(1,1,1,2))
     call gvputf(ng,1,kv,k1,k2,k3,fgg(1,3),fgrd(1,1,1,3))
     !       call zprm3('i Gx fn(G)',0,fgrd(1,1,1,1),n1,n2,n3)
     !       call zprm3('i Gy fn(G)',0,fgrd(1,1,1,2),n1,n2,n3)
     !       call zprm3('i Gz fn(G)',0,fgrd(1,1,1,3),n1,n2,n3)
  endif
  if (isw2 == 3 .OR. isw2 == 4 .OR. isw2 == 5 .OR. isw2 == 6) then
     call gvputf(ng,1,kv,k1,k2,k3,fg2,     flap)
  endif

  ! ... Gradient, laplacian in real space
  if (isw2 == 2 .OR. isw2 == 6) then
     call fftz3(fgrd(1,1,1,1),n1,n2,n3,k1,k2,k3,1,0,1)
     call fftz3(fgrd(1,1,1,2),n1,n2,n3,k1,k2,k3,1,0,1)
     call fftz3(fgrd(1,1,1,3),n1,n2,n3,k1,k2,k3,1,0,1)
     !       call zprm3('gradx fn(r)',0,fgrd(1,1,1,1),n1,n2,n3)
     !       call zprm3('grady fn(r)',0,fgrd(1,1,1,2),n1,n2,n3)
     !       call zprm3('gradz fn(r)',0,fgrd(1,1,1,3),n1,n2,n3)
  endif
  if (isw2 == 4 .OR. isw2 == 6) then
     call fftz3(flap         ,n1,n2,n3,k1,k2,k3,1,0,1)
     !       call zprm3('lap fn(r)',0,flap,n1,n2,n3)
  endif

  deallocate(fg,fgg,fg2)

  call tcx('grfmsh')
end subroutine grfmsh

