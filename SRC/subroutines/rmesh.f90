subroutine rmesh(z,rmax,lrel,lgrad,nrmx,a,nr) !- Generate parameters for shifted logarithmic radial mesh
  use m_lgunit,only:stdo
  !i   z     :nuclear charge
  !i   rmax  :augmentation radius, in a.u.
  !i   lrel  :0 for non-relativistic
  !i         :1 for scalar relativistic
  !i         :2 for Dirac equation
  !i   lgrad :0 for LDA, nonzero for gradient corrections
  !i   nrmx  :maximum allowed number of points
  ! o Inputs/Outputs
  ! o  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! o        :a is not altered if input a>0; otherwise a is set here.
  ! o        :When a is set, it is independent of rmax and nr
  ! o  nr    :number of radial mesh points
  ! o        :nr is not altered if input nr>0; otherwise nr is set here.
  ! o        :The calculated value of nr depends on both a and z
  !l Local variables
  !l         :
  !r Remarks
  !r   Uses input values for a,nr if >0; otherwise rmesh sets them
  !u Updates
  !u   18 Mar 03 Default parameters for fully relativistic case
  !u   11 Oct 02 No longer uses a smaller number of points for
  !u             the nonrelativistic case.
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nrmx,nr,lgrad,lrel,i1mach,iprint
  double precision :: z,rmax,a,b
  integer :: nrmax
  parameter (nrmax=1501)

  b = 1d0/(2*z+1)
  if (lrel == 2) then
     nr = nrmax
     a = .01d0
  elseif (lgrad /= 0) then
     if (a < 1d-6) a = 0.015d0
     if (nr <= 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
     !      No longer treat nonrelativistic case separately
     !      elseif (lrel .ne. 0) then
     !        if (a .lt. 1d-6) a = 0.03d0
     !        if (nr .le. 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
     !      else
     !        if (a .lt. 1d-6) a = 0.02d0
     !        if (nr .le. 0) nr = .5d0+dlog(1+rmax/b)/a
  else
     if (a < 1d-6) a = 0.03d0
     if (nr <= 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
  endif
  nr = max0(51,((nr-1)/2)*2+1)
  if (nrmx > 0) nr = min0(nr,nrmx)
  !     b = rmax/(dexp(a*(nr-1)) - 1d0)
  !      if (iprint() .ge. 50) print 333, z,a,nr,rmax
  !  333 format(' rmesh:  Z=',f5.1,'  a=',f6.3,'  nr=',i4,
  !     .  '  rmax=',f8.5)
  if (iprint() >= 50)write(stdo,"( &
       ' RMESH:  Z=',f9.3,'  a=',f9.3,'  nr=',i6,' rmax=',f9.4)")z,a,nr,rmax
end subroutine rmesh

