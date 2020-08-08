      subroutine dfphi(de1,del1,de2,del2,nr,g,gp,fivep)
C- Numerically differentiates g using five point finite difference
C ----------------------------------------------------------------------
Ci Inputs
Ci   de1   :see Remarks
Ci   del1  :see Remarks
Ci   de2   :see Remarks
Ci   del2  :see Remarks
Ci   nr    :number of mesh points
Ci   g     :normalized wave function times r
Ci   fivep :if true, five-point formula; else three-point formula
Cio Inputs/Outputs
Cio  gp    :On input,
Cio        :contains g at de1+del1, -de1+del1, de2+del2, -de2+del2
Cio        :On output, contains energy derivatives of g:
Cio        :two derivatives or four derivatives, depending on fivep
Cr Remarks
Cr   A simple five-point estimate for the numerical differentiation for
Cr   the first four derivatives of phi is based on integration of phi
Cr   at four energies de1+del1, -de1+del1, de2+del2, -de2+del2.  The
Cr   offset del1 and del2 are present because (for consistency's sake)
Cr   RSEQ is called which takes evenly spaced increments in slopes
Cr   about the average.  But the deviations del1 and del2 are only
Cr   nonzero to second order in the increment de, and this can be
Cr   exploited to obtain accurate five point estimates without having
Cr   to solve for the four simultaneous equations.
Cr
Cr   A three point for energy differences about Enu of 0, ep = de+del,
Cr   em = -de+del gives
Cr
Cr   (1)  gp =  2 del/(de^2-del^2) g(0) +
Cr              (de-del)/(de+del)/(2 de) g(ep) -
Cr              (de+del)/(de-del)/(2 de) g(em) +
Cr
Cr              1/6 gppp (de^2-del^2) + 1/12 gpppp del (de^2-del^2) + ...
Cr
Cr
Cr   (2) gpp = -2/(de^2-del^2) g(0) +
Cr              1/(de+del)/de g(ep) +
Cr              1/(de-del)/de g(em) -
Cr
Cr              2/3 gppp del + 1/12 gpppp del (de^2 + 3 del^2) + ...
Cr
Cr
Cr   The gppp term in (1) can be knocked out by taking two three point
Cr   formulas in linear combination (de1^2-del1^2) and (de2^2-del2^2)
Cr   leaving only fourth and higher order terms.  Because del is of
Cr   order de^2, the fourth order term is of the same order as the
Cr   fifth order term and there is no advantage in eliminating it.
Cr   Also from the difference between this more accurate (5-point)
Cr   estimate for gp, an estimate for gppp can be made as
Cr
Cr   (3) gppp = -6 ( gp(five point) - gp(three point) ) /(de1^2-del1^2)
Cr
Cr             + 1/2 gpppp del1
Cr
Cr   which is again accurate to order de^2.  Once gppp is known to this
Cr   order the term proportional to gppp in three point estimate for
Cr   gpp can be subtracted out directly and the gpppp term can be
Cr   eliminated by taking three point formulas (with the gppp term
Cr   removed) in linear combinations (de1^2+3*del1^2) and
Cr   (de2^2+3*del2^2), and finally the fourth order derivative can be
Cr   estimated from the difference in the five-point estimate for gpp
Cr   and the three point estimate.
C  ----------------------------------------------------------------
C     implicit none
C     Passed parameters
      integer nr
      logical fivep
      double precision de1,del1,de2,del2,g(nr),gp(nr,4)
C     Local parameters
      integer i
      double precision gp5p,gpp5p,gppp,gpp3p,gpp32
      double precision xx1,xx2,xx3,xx4,gp3p,w01,w11,w21,w02,w12,w22,
     .w01d,w11d,w21d,w02d,w12d,w22d,wp1d,wp2d,gpppp

C --- Constants common to 3-point and 5-point formulas ---
      w01 = 2*del1/(de1**2-del1**2)
      w11 = (de1-del1)/(de1+del1)/(2*de1)
      w21 = (de1+del1)/(de1-del1)/(2*de1)
      w01d = -2/(de1**2-del1**2)
      w11d = 1/(de1+del1)/de1
      w21d = 1/(de1-del1)/de1

      if (.not. fivep) goto 20
      if (dabs(del1/de1) .gt. .1 .or.  dabs(del1/de1) .gt. .1) then
C       if (iprint() .ge. 20) print *, 'dfphi:  large del; use 3 point'
        goto 20
      endif

C --- Extra constants for 5-point formula ---
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

C Three point formula for gp; store in temporary gp3p
        gp3p = w01*g(i) + w11*gp(i,1) - w21*gp(i,2)

C Five point estimate for gp
        gp5p = (xx2*gp3p - xx1*(w02*g(i) + w12*gp(i,3) - w22*gp(i,4)))
     .  /(xx2-xx1)

C Difference between five point and three point gives estimate for gppp
        gppp = -6/xx1*(gp5p - gp3p)

C Three point estimates for gpp with correction for gppp
        gpp3p = w01d * g(i) + w11d * gp(i,1) + w21d*gp(i,2)
     .  -wp1d * gppp
        gpp32 = w02d * g(i) + w12d * gp(i,3) + w22d*gp(i,4)
     .  -wp2d * gppp

C Five point estimate for gpp with correction for gppp
        gpp5p = (gpp3p*xx4 - gpp32*xx3) / (xx4 - xx3)

C Difference between five point and three point gives est for gpppp
        gpppp = -12/xx3*(gpp5p - gpp3p)

        gp(i,1) = gp5p
        gp(i,2) = gpp5p
        gp(i,3) = gppp
        gp(i,4) = gpppp
   10 continue
      return

C Three point formulae:  only gp, gpp calculated
   20 continue
      do  30  i = 1, nr
        gp3p    = w01*g(i)  + w11*gp(i,1)  - w21*gp(i,2)
        gp(i,2) = w01d*g(i) + w11d*gp(i,1) + w21d*gp(i,2)
        gp(i,1) = gp3p
   30 continue

      end

