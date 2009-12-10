      integer function chkhss(hess,n,w5,evtol,iopt,z,e)
C- Returns the number of eigenvalues greater than a specified value
Co z: eigenvectors if iopt > 1; otherwise, z unchanged
Co e: eigenvalues
Co hess: DESTROYED on output
C     implicit none
      integer n,iopt
      double precision hess(n,n),w5(n,5),evtol,z(n,n),e(n)
      double precision dum
      integer i,j

      j = 0
      chkhss = j
      if (iopt .eq. 0) return
      if (iopt .gt. 1) j = n
      call pshpr(0)
C     call prm(.false.,1,6,' ',w5,n,5)

      call dsev1(n,hess,dum,w5,0,.false.,.false.,.false.,
     .  j,9d9,i,z,e)
      call poppr
      j = 0
      do  60  i = 1, n
   60 if (e(i) .lt. evtol) j = j+1
      if (j .gt. 0 .and. iopt .eq. 1) then
        print 333, j,e(1)
  333   format(/' chkhss: Hessian has',i3,' negative evals; e(1)=',
     .    f8.4,'.  Not updated.')
      elseif (j .gt. 0) then
        print 332, j,e(1)
  332   format(/' chkhss: Hessian has',i3,' negative evals; e(1)=',
     .    f8.4,'.  Project unwanted part.')
C        call prmx('e',e(1),n,n,1)
C        call prmx('z',z,n,n,n)
        call rx('not implemented')
      endif
      chkhss = j

      end

