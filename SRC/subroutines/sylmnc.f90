subroutine sylmnc(c,lmx)
  !- Normalization constants for the spherical harmonics
  ! ----------------------------------------------------------------
  !i Inputs
  !i   lmx
  !o Outputs
  !o   C
  !r Remarks
  !r   use together with sylm (from ASW package)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: lmx
  double precision :: c(1)
  ! Local parameters
  integer :: i,i1,i2,l,lav,lp1,m,n1,n2,n3
  double precision :: fn2,fpi,tlp1,tpi,y0
  tpi = 8.d0*datan(1.d0)
  fpi = 2.d0*tpi
  y0 = 1.d0/dsqrt(fpi)
  c(1) = y0
  do  21  l = 1, lmx
     lp1 = l+1
     tlp1 = l+lp1
     lav = l*lp1 + 1
     c(lav) = dsqrt(tlp1/fpi)
     do  2  m = 1, l
        n2 = lp1-m
        n1 = n2+1
        n3 = l+m
        fn2 = n2
        do  1  i = n1, n3
           fn2 = fn2*i
1       enddo
        i1 = lav+m
        i2 = lav-m
        c(i1) = dsqrt(tlp1/(fn2*tpi))
        c(i2) = c(i1)
2    enddo
21 enddo
end subroutine sylmnc


