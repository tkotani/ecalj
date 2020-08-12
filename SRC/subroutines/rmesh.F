      subroutine rmesh(z,rmax,lrel,lgrad,nrmx,a,nr)
C- Generate parameters for shifted logarithmic radial mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   lrel  :0 for non-relativistic
Ci         :1 for scalar relativistic
Ci         :2 for Dirac equation
Ci   lgrad :0 for LDA, nonzero for gradient corrections
Ci   nrmx  :maximum allowed number of points
Cio Inputs/Outputs
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio        :a is not altered if input a>0; otherwise a is set here.
Cio        :When a is set, it is independent of rmax and nr
Cio  nr    :number of radial mesh points
Cio        :nr is not altered if input nr>0; otherwise nr is set here.
Cio        :The calculated value of nr depends on both a and z
Cl Local variables
Cl         :
Cr Remarks
Cr   Uses input values for a,nr if >0; otherwise rmesh sets them
Cu Updates
Cu   18 Mar 03 Default parameters for fully relativistic case
Cu   11 Oct 02 No longer uses a smaller number of points for
Cu             the nonrelativistic case.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrmx,nr,lgrad,lrel,i1mach,iprint
      double precision z,rmax,a,b
      integer nrmax
      parameter (nrmax=1501)

      b = 1d0/(2*z+1)
      if (lrel .eq. 2) then
        nr = nrmax
        a = .01d0
      elseif (lgrad .ne. 0) then
        if (a .lt. 1d-6) a = 0.015d0
        if (nr .le. 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
C      No longer treat nonrelativistic case separately
C      elseif (lrel .ne. 0) then
C        if (a .lt. 1d-6) a = 0.03d0
C        if (nr .le. 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
C      else
C        if (a .lt. 1d-6) a = 0.02d0
C        if (nr .le. 0) nr = .5d0+dlog(1+rmax/b)/a
      else
        if (a .lt. 1d-6) a = 0.03d0
        if (nr .le. 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
      endif
      nr = max0(51,((nr-1)/2)*2+1)
      if (nrmx .gt. 0) nr = min0(nr,nrmx)
C     b = rmax/(dexp(a*(nr-1)) - 1d0)
C      if (iprint() .ge. 50) print 333, z,a,nr,rmax
C  333 format(' rmesh:  Z=',f5.1,'  a=',f6.3,'  nr=',i4,
C     .  '  rmax=',f8.5)
      if (iprint() .ge. 50)
     .call awrit4(' RMESH:  Z=%d  a=%1;4d  nr=%i  rmax=%1;6d',
     .' ',80,i1mach(2),z,a,nr,rmax)
      end

