      subroutine shear(nbas,plat,tau,alpha,eps,s)
C- Apply a pure strain to lattice and basis vectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis.  nbas=0 => plat, tau not sheared
Ci   plat  :basis vectors, in units of alat
Ci   tau   :position vectors, in units of alat
Ci   alpha :tight-binding screening parameters
Co   eps   :alpha*eps(1..6):  Voigt tensor strains
Co Outputs
Co   s     :deformation matrix
Co         :s and eps can occupy the same address space
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   19 Mar 06 Added output s
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas
      double precision plat(3,3),tau(3,nbas),eps(6),alpha
C ... Local parameters
      integer i,j,ind(2,21),n,iprint,nbmx,nglob,stdo
      parameter (nbmx = 256)
      double precision s(3,3),b(3,3),x(21),t(3,nbmx),det,e(6)

      stdo = nglob('stdo')
      if (dabs(alpha) .lt. 1d-8) return
      if (nbas .gt. nbmx)
     .  call fexit(-1,9,'Increase nbmx in SHEAR',0)
      call dcopy(6,eps,1,e,1)
      call dscal(6,alpha,e,1)
      call xxxes(e,s,det)

C ... Shear plat, tau
      if (nbas .gt. 0) then
        if (iprint() .ge. 40) then
          write (stdo,10) plat
          write (stdo,20)
          write (stdo,30) ((tau(i,j),i=1,3),j=1,nbas)
          write (stdo,40) s
        endif
        call dmpy(s,3,1,plat,3,1,b,3,1,3,3,3)
        call dcopy(9,b,1,plat,1)
        call dmpy(s,3,1,tau,3,1,t,3,1,3,nbas,3)
        call dcopy(3*nbas,t,1,tau,1)
        if (iprint() .ge. 40) then
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
    1 continue
      do  2  i = 1, 5
      do  2  j = i+1, 6
        call xxxadd(i,j,n,e,ind,x)
    2 continue
      if (n .eq. 0 .or. iprint() .lt. 30) return
      write (stdo,50) alpha,det-1
      write (stdo,60) ((ind(i,j),i=1,2),j=1,n)
      write (stdo,70) (x(i),i=1,n)
   10 format (' SHEAR: Lattice vectors:'/3(8x,3f10.6/))
   20 format ('        Basis atoms:')
   30 format (8x,3f10.6)
   40 format ('       Lattice and basis sheared by'/3(8x,3f10.6/))
   50 format (/
     .   1x,'SHEAR: distortion amplitude =',f9.6,
     .      '  Vol. dilatation =',f9.6/
     .   8x,'The second derivative E/vol w.r.t alpha = ',
     .      'W'''' = sum_ij x_ij c_ij'/
     .   8x,'has coefficients x_ij to elastic constants c_ij as shown:')
   60 format (6(5x,2i1,5x))
   70 format (6f12.6)
      end
      subroutine xxxes(e,s,det)
C Make deformation tensor
      double precision e(6),s(3,3),det
      s(1,1) = 1 + e(1)
      s(2,2) = 1 + e(2)
      s(3,3) = 1 + e(3)
      s(1,2) = e(6)
      s(2,1) = e(6)
      s(1,3) = e(5)
      s(3,1) = e(5)
      s(2,3) = e(4)
      s(3,2) = e(4)
      det=s(1,1)*s(2,2)*s(3,3)+s(1,2)*s(2,3)*s(3,1)
     .   +s(1,3)*s(2,1)*s(3,2)-s(1,3)*s(2,2)*s(3,1)
     .   -s(1,2)*s(2,1)*s(3,3)-s(1,1)*s(2,3)*s(3,2)
      end
      subroutine xxxse(s,e)
C Make engineering strains
      double precision s(3,3),e(6)
      integer i
      do  1  i = 1, 3
        e(i) = s(i,i) - 1d0
    1 continue
      e(4) = 2*s(2,3)
      e(5) = 2*s(1,3)
      e(6) = 2*s(1,2)
      end
      subroutine xxxadd(i,j,n,e,ind,x)
C Add to list of non-zero elastic constant coefficients
C     implicit none
      integer i,j,n,ind(2,21)
      double precision e(6),x(21),xx
      xx = e(i)*e(j)
      if (i .ne. j) xx = 2*xx
      if (dabs(xx) .gt. 1d-8) then
        n = n+1
        ind(1,n) = i
        ind(2,n) = j
        x(n) = xx
      endif
      end


