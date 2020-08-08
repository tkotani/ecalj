      implicit none
      integer(4):: iq,nq,iline
      real (8):: j1,j2,jpp0,jppf,jpmf,jpp,jpm,h0,q(3)
     & ,qinit(3),qend(3),jpm0,omega,mmom,pi,x,dx
      pi =4d0*atan(1d0)
c      print *,pi
      j1= -13.2   !Solovyev
      j2= -23.5   !97.86 at .25 .25 .25
      mmom = 4.50d0

c      j1=0d0
c      j2=0d0

c      j1= -30
c      j2= -30  ! 160.0

c      j1= -9.8  /11.605*2.5**2 ! Pask
c      j2= -24.5 /11.605*2.5**2  

c      j1=-0.77* 2.5**2 !Kohgi
c      j2=-0.89* 2.5**2
      jpp0 = jppf((/0d0,0d0,0d0/),j1,j2)  !NN   intra
      jpm0 = jpmf((/0d0,0d0,0d0/),j1,j2)  !NNN  inter
      h0 =  jpp0-jpm0

      x=0
      do iline = 1,5
       nq=30
       if(iline==1) then
          qinit= (/0d0,0d0,0d0/)
          qend = (/0.25d0,0.25d0,0.25d0/) 
       elseif(iline==2) then
         qinit = qend
         qend  = (/.5,.5,-.5/) 
       elseif(iline==3) then
         qinit  = qend
         qend  =  (/0d0,0d0,0d0/)
       elseif(iline==4) then
         qinit  = qend
         qend =  (/-0.25,-0.25,.75/)
       elseif(iline==5) then
         qinit = (/0d0,0d0,0d0/)
         qend  = (/.5,.5,-.5/)
       endif
       dx = sqrt(sum((qend-qinit)**2))/dble(nq)
       do iq= 0,nq
        if(iq/=0) x=x+dx
        q= qinit + (qend-qinit) *(iq)/dble(nq)
        jpp = jppf( 2*pi*q, j1,j2)
        jpm = jpmf( 2*pi*q, j1,j2)
c        write(6,"('xxx: ',5d13.5)") q, jpp-h0,jpm
        omega = 2d0/mmom *( (jpp-h0)**2 - jpm**2 )**0.5d0
c        omega = 2d0/mmom *( jpp**2 - jpm**2 )**0.5d0
        write(6,"(f9.5, ' q=',3f9.5,' omg=',f13.5)") x,q, omega
       enddo
      enddo
      end
C--- NNN
      subroutine setNNN(a,nnn)
      real(8) :: a(3,6),pi
      integer(4):: nnn
      a(:,1) = (/ 1, 0,0/)
      a(:,2) = (/-1, 0,0/)
      a(:,3) = (/ 0, 1,0/)
      a(:,4) = (/ 0,-1,0/)
      a(:,5) = (/ 0, 0, 1/)
      a(:,6) = (/ 0, 0,-1/)
      nnn=6
      end
C--- NN
      subroutine setNN(a,nn)
      real(8) :: a(3,12)
      integer(4):: nn
      a(:,1) = (/ 1, 1,0/)
      a(:,2) = (/-1, 1,0/)
      a(:,3) = (/ 1, -1,0/)
      a(:,4) = (/-1, -1,0/)
      a(:,5) = (/0, 1, 1/)
      a(:,6) = (/0,-1, 1/)
      a(:,7) = (/0, 1, -1/)
      a(:,8) = (/0,-1, -1/)
      a(:,9)  = (/ 1, 0,  1/)
      a(:,10) = (/-1,0,  1/)
      a(:,11) = (/ 1,0, -1/)
      a(:,12) = (/-1,0, -1/)
      a = a*0.5d0
      nn=12
      end

C---------------------
      real(8) function jpmf(q,j1,j2)
! Inter sublattice
      implicit none
      integer(4):: i,nn,nnn
      complex(8):: img=(0d0,1d0)
      real(8) :: a(3,12),pi,q(3),j1,j2,aaf(3)
c      aaf(1:3)= 0d0 !(/1d0,1d0,1d0/)
      call setNNN(a,nnn)
      jpmf=0d0
      do i=1,nnn
         jpmf    = jpmf + j2* exp ( img* sum(a(1:3,i)*q(1:3)) )
      enddo
c      write(6,"(' q=',3f9.5,' jpmf=',f13.5)") q,jpmf
      call setNN(a,nn)
      do i=1,nn
         if( abs( sum(a(:,i)) )>1d-8) then !not in 111 plane
           jpmf = jpmf + j1* exp ( img* sum (a(1:3,i)*q(1:3)) )
         endif
      enddo
      end
c---------------------
      real(8) function jppf(q,j1,j2)
!  intra subrattice
      implicit none
      integer(4):: i,nn
      complex(8):: img=(0d0,1d0)
      real(8) :: a(3,12),pi,q(3),j1,j2
      call setNN(a,nn)
      jppf=0d0
      do i=1,nn
         if( abs( sum(a(:,i)) )<1d-8) then  !in plane
           jppf = jppf + j1* exp( img* sum(a(1:3,i)*q(1:3)) )
         endif
      enddo
      end
