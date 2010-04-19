      subroutine efgrad(z,efg)
C- Computation of electric field gradient
C     implicit none
      double precision z,efg(5)

      integer i,j,n,ier
      double precision v(3,3),vi(3,3),d(3),e(3),e2(3),tau(2,3)
      double precision tr(3,3),ti(3,3)
      double precision conv1,conv2,pi,f0,s3,
     .  ax,bx,cx,dx,ex
      data conv1 /162.1/
      pi=4.0*atan(1.0)
      f0=sqrt(15.0/16.0/pi)
      s3=sqrt(3.0)
      conv2=conv1*3.0
      
      write(*,'(11x,"Tensor axes",5x,"esu/cm^2",2x,
     .          "  V/m^2")')
      write(*,'("      ",5x,"           ",5x," x10^13 ",2x,
     .          "  x10^19 ")')
      
      ax = efg(1)
      bx = efg(2)
      cx = efg(3)
      dx = efg(4)
      ex = efg(5)
      v(1,1) = 2*ex - 2/s3*cx
      v(2,2) = -2*ex - 2/s3*cx
      v(3,3) = 4/s3*cx
      v(1,2) = 2*ax
      v(1,3) = 2*dx
      v(2,3) = 2*bx
      v(2,1) = v(1,2)
      v(3,1) = v(1,3)
      v(3,2) = v(2,3)
      do  2  i = 1,3
        do  1  j = 1,3
          v(i,j) = f0*v(i,j)
          vi(i,j) = 0.d0
          tr(i,j) = 0.d0
          ti(i,j) = 0.d0
    1   continue 
        tr(i,i) = 1.d0
    2 continue 
      n = 3
      call htridi(n,n,v,vi,d,e,e2,tau)
      call imtql2(n,n,d,e,tr,ier)
      call rxx (ier.gt.0,' efgrad : IER ne 0')
      call htribk(n,n,v,vi,tau,n,tr,ti)
      do  5  i = 1,3
        write(*,'(4x,3x,3f6.2,2x,f8.2,2x,f8.2,5x)')
     .           (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)          
    5 continue 
      end

