      subroutine pade(nz,z,nc,zc,ndc,c,f)
C- Generate Pade approximation to a function at a list of points
C ----------------------------------------------------------------
Ci Inputs
Ci   z,nz:  points at which to evaluate Pade approximant, and no.
Ci   zc,nc: points used to generate Pade coefficients, and number
Ci   c,ndc: Pade coefficients; see padcof for their generation,
Ci          and leading dimension of c
Ci          c must be dimensioned c(ndc,nc+2)
Co Outputs
Co   f:     function values at the points zc.
C ----------------------------------------------------------------
C     implicit none
      integer i,nc,nz,im,ndc
      double complex f(nz),z(nz),zc(*),c(ndc,nc+2)
      c(1,nc+1) = c(1,1)
      c(1,nc+2) = 1
      do  10  im = 1, nz
        c(2,nc+1) = c(1,nc+1)
        c(2,nc+2) = c(1,nc+2) + (z(im)-zc(1))*c(2,2)
        do  20  i = 2, nc-1
          c(i+1,nc+1) = c(i,nc+1) + (z(im)-zc(i))*c(i+1,i+1)*c(i-1,nc+1)
          c(i+1,nc+2) = c(i,nc+2) + (z(im)-zc(i))*c(i+1,i+1)*c(i-1,nc+2)
   20   continue
        f(im) = c(nc,nc+1)/c(nc,nc+2)
   10 continue
      end
      subroutine padcof(n,z,f,ndc,c)
C- Generate by recursion coefficients for Pade approximation
C ----------------------------------------------------------------
Ci Inputs
Ci   n:    number of points at which function is known
Ci   z,f:  points and function values
Ci   ndc:  leading dimension of c array
Co Outputs
Co   c:    Pade coefficients.  
Co         NB: padcof requires that c be dimensioned (ndc,n), but
Co         pade requires the dimensions to be (ndc,n+2)
Cu Updates
Cu   9 Mar 00 Handle case all f are zero
C ----------------------------------------------------------------
C     implicit none
      integer i,j,n,ndc
      double precision f(2,*),z(2,*),c(2,ndc,n)
      logical lnz
#if NOQUAD
      double precision zn(2),xd(2),yd(2),zd(2)
#else
      real(16) zn(2),xd(2),yd(2),zd(2)
#endif
      lnz = .false.
      do  10  j = 1, n
        c(1,1,j) = f(1,j)
        c(2,1,j) = f(2,j)
        lnz = lnz .or. f(1,j) .ne. 0 .or. f(2,j) .ne. 0
   10 continue
      do  20  j = 2, n
C ... Do the arithmetic with quad precision
        do  20  i = 2, j
C       c(i,j) = (c(i-1,i-1)-c(i-1,j)) / ( (z(j)-z(i-1)) * c(i-1,j) )
          zn(1) = c(1,i-1,i-1)-c(1,i-1,j)
          zn(2) = c(2,i-1,i-1)-c(2,i-1,j)
          xd(1) = z(1,j)-z(1,i-1)
          xd(2) = z(2,j)-z(2,i-1)
          yd(1) = c(1,i-1,j)
          yd(2) = c(2,i-1,j)
          zd(1) = xd(1)*yd(1) - xd(2)*yd(2) 
          zd(2) = xd(1)*yd(2) + xd(2)*yd(1) 
          if (lnz) call qdiv(zn(1),zn(2),zd(1),zd(2),zn(1),zn(2))
          c(1,i,j) = zn(1)
          c(2,i,j) = zn(2)
   20 continue

      end

