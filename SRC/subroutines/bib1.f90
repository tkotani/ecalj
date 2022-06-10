
subroutine freq01 (nx,ua,   freqx,freqw,wx,expa)
  ! generates a gaussian point x between (0,1) and w = (1-x)/x
  ! and the weights in x
  ! also generates expa = exp(-ua^2 w^2)
  ! nx    = no. gaussian points
  !  ua   = s.o.
  ! freqx = gaussian points
  ! freqw = (1-x)/x
  ! wx    = weights of gaussian points x
  ! expa  = s.o.
  ! originally 92.02.27 Ferdi.Aryasetiawan
  implicit real*8 (a-h,o-z)
  integer:: nx,ix
  real(8):: freqx(nx),freqw(nx),wx(nx),expa(nx)
  call gauss   (nx,0.d0,1.d0,freqx,wx)! generate gaussian points
  ! calculate w = 1/(1+x)
  ua2        = ua*ua
  do      ix = 1,nx
     freqw(ix)  = (1.d0 - freqx(ix)) / freqx(ix)
     expa(ix)   = dexp(-ua2*freqw(ix)*freqw(ix))
  end do
  return
end subroutine freq01
!--------------------------------------------------------------------
subroutine freq01x (nx,    freqx,freqw,wx) !,expa)
  ! 92.02.27
  ! generates a gaussian point x between (0,1) and w = (1-x)/x
  ! and the weights in x
  ! also generates expa = exp(-ua^2 w^2)
  ! nx    = no. gaussian points
  !  ua   = s.o.
  ! freqx = gaussian points
  ! freqw = (1-x)/x
  ! wx    = weights of gaussian points x
  ! expa  = s.o.
  implicit real*8 (a-h,o-z)
  real(8):: freqx(nx),freqw(nx),wx(nx),expa(nx)
  integer:: nx,ix
  ! generate gaussian points
  call gauss   (nx,0.d0,1.d0,freqx,wx)
  ! calculate w = 1/(1+x)
  !      ua2        = ua*ua
  write(6,"(' --- freq01x:  ix    x    freqw(a.u.)---')")
  do      ix = 1,nx
     freqw(ix)  = (1.d0 - freqx(ix)) / freqx(ix)
     write(6,"('            ',i4,2f9.4)") ix,freqx(ix),freqw(ix)
  end do
  return
end subroutine freq01x
