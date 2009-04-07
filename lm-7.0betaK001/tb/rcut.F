      subroutine rcut(val,slo,curv,r1,rc,c)
C - Make coefficients to order-5 polynomial to cut off smoothly to zero
C-----------------------------------------------------------------------
Ci Input: val,slo,curve,r1,rc
Co Output: c, coefficients to fifth order polynomial
Cr Remarks
Cr For a scaling law or pair potential (or any function), given
Cr its value, slope and curvature at some r1 and a cut-off distance
Cr rc, this returns a polynomial, valid in the range r1 < r < rc which
Cr is to augment (replace) the function in this range so it goes
Cr smoothly to zero at rc
C-----------------------------------------------------------------------
C     implicit none
C Passed paramters
      double precision r1,rc,val,slo,curv,c(0:5)

C Local paramters
      integer i,j,iprint,i1mach,ierr
      double precision r,a(6,6),w(36),w2(7)
      double precision p,pp,ppp
C --- polynomial and its derivatives ---
      p(r) = c(0)+c(1)*r+c(2)*r**2+c(3)*r**3+c(4)*r**4+c(5)*r**5 
      pp(r) = c(1)+2d0*c(2)*r+3d0*c(3)*r**2+4d0*c(4)*r**3+5d0*c(5)*r**4
      ppp(r) = 2d0*c(2)+6d0*c(3)*r+12d0*c(4)*r**2+2d1*c(5)*r**3

C --- make vector ---
      c(0) = val
      c(1) = slo
      c(2) = curv
      call dcopy(3,0d0,0,c(3),1)
C --- make matrix ---
      a(1,1) =  1d0
      a(1,2) =  r1
      a(1,3) =  r1**2
      a(1,4) =  r1**3
      a(1,5) =  r1**4
      a(1,6) =  r1**5
      a(2,1) =  0d0
      a(2,2) =  1d0
      a(2,3) =  2d0 * r1
      a(2,4) =  3d0 * r1**2
      a(2,5) =  4d0 * r1**3
      a(2,6) =  5d0 * r1**4
      a(3,1) =  0d0
      a(3,2) =  0d0
      a(3,3) =  2d0
      a(3,4) =  6d0 * r1
      a(3,5) = 12d0 * r1**2
      a(3,6) = 20d0 * r1**3
      a(4,1) =  1d0
      a(4,2) =  rc
      a(4,3) =  rc**2
      a(4,4) =  rc**3
      a(4,5) =  rc**4
      a(4,6) =  rc**5
      a(5,1) =  0d0
      a(5,2) =  1d0
      a(5,3) =  2d0 * rc
      a(5,4) =  3d0 * rc**2
      a(5,5) =  4d0 * rc**3
      a(5,6) =  5d0 * rc**4
      a(6,1) =  0d0
      a(6,2) =  0d0
      a(6,3) =  2d0
      a(6,4) =  6d0 * rc
      a(6,5) = 12d0 * rc**2
      a(6,6) = 20d0 * rc**3
      
      if (iprint() .gt. 120) then
        call awrit0(' RCUT: matrix:',' ',128,i1mach(2))
        do  i = 1, 6
          write (*, 1) (a(i,j), j = 1, 6)
        enddo
        call awrit0(' RCUT: vector:',' ',128,i1mach(2))
        write (*,1) (c(i), i = 0, 5)
      endif

C --- solve the linear problem ---
      call dqinvb(' ',a,6,2,6,1,w,6,w2,c,6,ierr)

      if (ierr .ne. 0) then
        if (iprint() .gt. 10) then
          call awrit1(' RCUT failed to invert matrix, ierr=%i',' ',128
     .      ,i1mach(2),ierr)
        endif
        call rx(' RCUT exiting ..')
      else
        if (iprint() .gt. 120) then
          call awrit0(' RCUT: coefficients:',' ',128,i1mach(2))
          write (*,1) (c(i), i = 0, 5)
          write (*,2)
          write (*,3) val,p(r1),slo,pp(r1),curv,ppp(r1)
          write (*,4)
          write (*,5) p(rc),pp(rc),ppp(rc)
        endif
      endif
    1 format (6g15.4)
    2 format (8x,'h(r1)',8x,'p(r1)',9x,'h''(r1)',9x,'p''(r1)',9x,
     .        'h''''(r1)',9x,'p''''(r1)')
    3 format (6g15.4)
    4 format (8x,'p(rc)',8x,'p''(rc)',8x,'p''''(rc)')
    5 format (3g15.4)
      end
