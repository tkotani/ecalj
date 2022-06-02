subroutine dfphi(de1,del1,de2,del2,nr,g,gp,fivep)
  !- Numerically differentiates g using five point finite difference
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   de1   :see Remarks
  !i   del1  :see Remarks
  !i   de2   :see Remarks
  !i   del2  :see Remarks
  !i   nr    :number of mesh points
  !i   g     :normalized wave function times r
  !i   fivep :if true, five-point formula; else three-point formula
  ! o Inputs/Outputs
  ! o  gp    :On input,
  ! o        :contains g at de1+del1, -de1+del1, de2+del2, -de2+del2
  ! o        :On output, contains energy derivatives of g:
  ! o        :two derivatives or four derivatives, depending on fivep
  !r Remarks
  !r   A simple five-point estimate for the numerical differentiation for
  !r   the first four derivatives of phi is based on integration of phi
  !r   at four energies de1+del1, -de1+del1, de2+del2, -de2+del2.  The
  !r   offset del1 and del2 are present because (for consistency's sake)
  !r   RSEQ is called which takes evenly spaced increments in slopes
  !r   about the average.  But the deviations del1 and del2 are only
  !r   nonzero to second order in the increment de, and this can be
  !r   exploited to obtain accurate five point estimates without having
  !r   to solve for the four simultaneous equations.
  !r
  !r   A three point for energy differences about Enu of 0, ep = de+del,
  !r   em = -de+del gives
  !r
  !r   (1)  gp =  2 del/(de^2-del^2) g(0) +
  !r              (de-del)/(de+del)/(2 de) g(ep) -
  !r              (de+del)/(de-del)/(2 de) g(em) +
  !r
  !r              1/6 gppp (de^2-del^2) + 1/12 gpppp del (de^2-del^2) + ...
  !r
  !r
  !r   (2) gpp = -2/(de^2-del^2) g(0) +
  !r              1/(de+del)/de g(ep) +
  !r              1/(de-del)/de g(em) -
  !r
  !r              2/3 gppp del + 1/12 gpppp del (de^2 + 3 del^2) + ...
  !r
  !r
  !r   The gppp term in (1) can be knocked out by taking two three point
  !r   formulas in linear combination (de1^2-del1^2) and (de2^2-del2^2)
  !r   leaving only fourth and higher order terms.  Because del is of
  !r   order de^2, the fourth order term is of the same order as the
  !r   fifth order term and there is no advantage in eliminating it.
  !r   Also from the difference between this more accurate (5-point)
  !r   estimate for gp, an estimate for gppp can be made as
  !r
  !r   (3) gppp = -6 ( gp(five point) - gp(three point) ) /(de1^2-del1^2)
  !r
  !r             + 1/2 gpppp del1
  !r
  !r   which is again accurate to order de^2.  Once gppp is known to this
  !r   order the term proportional to gppp in three point estimate for
  !r   gpp can be subtracted out directly and the gpppp term can be
  !r   eliminated by taking three point formulas (with the gppp term
  !r   removed) in linear combinations (de1^2+3*del1^2) and
  !r   (de2^2+3*del2^2), and finally the fourth order derivative can be
  !r   estimated from the difference in the five-point estimate for gpp
  !r   and the three point estimate.
  !  ----------------------------------------------------------------
  !     implicit none
  !     Passed parameters
  integer :: nr
  logical :: fivep
  double precision :: de1,del1,de2,del2,g(nr),gp(nr,4)
  !     Local parameters
  integer :: i
  double precision :: gp5p,gpp5p,gppp,gpp3p,gpp32
  double precision :: xx1,xx2,xx3,xx4,gp3p,w01,w11,w21,w02,w12,w22, &
       w01d,w11d,w21d,w02d,w12d,w22d,wp1d,wp2d,gpppp

  ! --- Constants common to 3-point and 5-point formulas ---
  w01 = 2*del1/(de1**2-del1**2)
  w11 = (de1-del1)/(de1+del1)/(2*de1)
  w21 = (de1+del1)/(de1-del1)/(2*de1)
  w01d = -2/(de1**2-del1**2)
  w11d = 1/(de1+del1)/de1
  w21d = 1/(de1-del1)/de1

  if ( .NOT. fivep) goto 20
  if (dabs(del1/de1) > .1 .OR.  dabs(del1/de1) > .1) then
     !       if (iprint() .ge. 20) print *, 'dfphi:  large del; use 3 point'
     goto 20
  endif

  ! --- Extra constants for 5-point formula ---
  xx1 = de1**2 - del1**2
  xx2 = de2**2 - del2**2
  xx3 = de1**2 + 3*del1**2
  xx4 = de2**2 + 3*del2**2
  w02 = 2*del2/(de2**2-del2**2)
  w12 = (de2-del2)/(de2+del2)/(2*de2)
  w22 = (de2+del2)/(de2-del2)/(2*de2)
  wp1d = 2d0*del1/3
  w02d = -2/(de2**2-del2**2)
  w12d = 1/(de2+del2)/de2
  w22d = 1/(de2-del2)/de2
  wp2d = 2d0*del2/3

  do  10  i = 1, nr

     ! Three point formula for gp; store in temporary gp3p
     gp3p = w01*g(i) + w11*gp(i,1) - w21*gp(i,2)

     ! Five point estimate for gp
     gp5p = (xx2*gp3p - xx1*(w02*g(i) + w12*gp(i,3) - w22*gp(i,4))) &
          /(xx2-xx1)

     ! Difference between five point and three point gives estimate for gppp
     gppp = -6/xx1*(gp5p - gp3p)

     ! Three point estimates for gpp with correction for gppp
     gpp3p = w01d * g(i) + w11d * gp(i,1) + w21d*gp(i,2) &
          -wp1d * gppp
     gpp32 = w02d * g(i) + w12d * gp(i,3) + w22d*gp(i,4) &
          -wp2d * gppp

     ! Five point estimate for gpp with correction for gppp
     gpp5p = (gpp3p*xx4 - gpp32*xx3) / (xx4 - xx3)

     ! Difference between five point and three point gives est for gpppp
     gpppp = -12/xx3*(gpp5p - gpp3p)

     gp(i,1) = gp5p
     gp(i,2) = gpp5p
     gp(i,3) = gppp
     gp(i,4) = gpppp
10 enddo
  return

  ! Three point formulae:  only gp, gpp calculated
20 continue
  do  30  i = 1, nr
     gp3p    = w01*g(i)  + w11*gp(i,1)  - w21*gp(i,2)
     gp(i,2) = w01d*g(i) + w11d*gp(i,1) + w21d*gp(i,2)
     gp(i,1) = gp3p
30 enddo

end subroutine dfphi

