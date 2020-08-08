! This is an example, to determine J1+,J1-,J2
! From the date of SW at M/2,M,X.
!
!         if(abs(d1-121.93)<.05) then
!         if(abs(d2-119.96)<.05) then
!         if(abs(d3-24.15) <.05) then
!      121.93,... are SWenergy at M/2,M,X.
! 
      implicit none
      integer(4):: iq,nq,iline,ifmat,ix,iy,nmbas,iffi,ixx,
     & i1p,i1m,i2
      real (8):: j1,j2,jpp0,jpp,jpm,h0,q(3),j1p,j1m
     & ,qinit(3),qend(3),jpm0,omega,mmom,pi,x,dx,meff,omega2
     & ,d1,d2,d3
! ,jpmf
      complex(8):: jpmf,jppf,phase,img=(0d0,1d0)
      complex(8),allocatable:: x0matinv0(:,:)
      pi =4d0*atan(1d0)

!--- We obtained 
cQSGW
c j1p j1m j2 mmom=     -2.0310     -4.5850     -5.7470      5.3600
c d1 d2 d3=     23.1313     18.8992     14.2174
c
cLDA
c j1p j1m j2 mmom=    -12.1200    -12.1200    -18.7600      4.2500
c d1 d2 d3=     87.1906     60.1984      0.0000

cQSGW
c j1p j1m j2 mmom=     -0.4500     -0.7550    -19.0525      1.9500
c d1 d2 d3=    121.8923    120.0062     24.1939

cLDA
c j1p j1m j2 mmom=      0.2650      0.2650    -28.2700      1.0900
c d1 d2 d3=    308.3119    310.2508      0.0000


c xxxxxxxxxxxxxxxx
c      mmom= 4.25

      do i1m = 1,50
      do i1p = 1,50
      do i2  = 1,100 

c         mmom=5.36
c         j1p = -3
c         j1m = -3
c         j2  = -7

c         j1p= -2*117*13.605/1000
c         j1m= -2*235*13.605/1000
c         j2=  -2*242*13.605/1000
c         j1p=   -0 !.265
c         j2 =   -20
c         j1m=j1p

         mmom=1.95
         j1p =  -0.4  - i1p*0.0025  !trial values
         j1m =  -0.7  - i1m*0.0025
         j2  = -18.9  - i2 *0.0025

c         write(6,"(' j1p j1m j2 mmom=',4f12.4)") j1p,j1m,j2,mmom
         call calsw(j1p,j1m,j2,mmom,d1,d2,d3)
c         write(6,"(' d1 d2 d3=',4f12.4)") d1,d2,d3

c         if(abs(d1-23.18)<0.05) then
c         if(abs(d2-18.85)<0.05) then
c         if(abs(d3-14.18)<0.05) then

c         if(abs(d1-87.19) <.02) then
c         if(abs(d2-60.20) <.02) then
c           write(6,"(' j1p j1m j2 mmom=',4f12.4)") j1p,j1m,j2,mmom
c           write(6,"(' d1 d2 d3=',4f12.4)") d1,d2,d3
c            stop "xxxxxxxxxxxxxxxx"
c         endif
c         endif

c         if(abs(d1-308.40) <.1) then
c         if(abs(d2-310.18) <.1) then
c           write(6,"(' j1p j1m j2 mmom=',4f12.4)") j1p,j1m,j2,mmom
c           write(6,"(' d1 d2 d3=',4f12.4)") d1,d2,d3
c            stop "xxxxxxxxxxxxxxxx"
c         endif
c         endif
c         endif
c      enddo

         if(abs(d1-121.93)<.05) then
         if(abs(d2-119.96)<.05) then
         if(abs(d3-24.15) <.05) then
           write(6,"(' j1p j1m j2 mmom=',4f12.4)") j1p,j1m,j2,mmom
           write(6,"(' d1 d2 d3=',4f12.4)") d1,d2,d3
            stop "xxxxxxxxxxxxxxxx"
         endif
         endif
         endif

c           write(6,"(' j1p j1m j2 mmom=',4f12.4)") j1p,j1m,j2,mmom
c           write(6,"(' d1 d2 d3=',4f12.4)") d1,d2,d3
c            stop "xxxxxxxxxxxxxxxx"
      enddo      
      enddo
      enddo
      end

      subroutine calsw(j1p,j1m,j2,mmom,d1,d2,d3)
      implicit none
      integer(4):: iq,nq,iline,ifmat,ix,iy,nmbas,iffi,ixx,iqx
      real (8):: j1,j2,jpp0,jpp,jpm,h0,q(3),j1p,j1m
     & ,qinit(3),qend(3),jpm0,omega,mmom,pi,x,dx,meff,omega2,
     &  d1,d2,d3
! ,jpmf
      complex(8):: jpmf,jppf,phase,img=(0d0,1d0)
      complex(8),allocatable:: x0matinv0(:,:)
      pi =4d0*atan(1d0)
c      print *,pi
c      j1= -13.2   !Solovyev
c      j2= -23.5   !97.86 at .25 .25 .25

c      j1= -13.2   !Solovyev
c      j2= -23.5   !97.86 at .25 .25 .25

c      j1p=j1
c      j1m=j1
c      j1=0d0
c      j2=0d0

c      j1= -30
c      j2= -30  ! 160.0

c      j1= -9.8  /11.605*2.5**2 ! Pask
c      j2= -24.5 /11.605*2.5**2  

c      j1=-0.77* 2.5**2 !Kohgi
c      j2=-0.89* 2.5**2

c      mmom = 4.50d0

      jpp0 = jppf((/0d0,0d0,0d0/),j1p,j2)  !NN   intra
      jpm0 = jpmf((/0d0,0d0,0d0/),j1m,j2)  !NNN  inter
      h0 =  jpp0-jpm0

      x=0

      iqx=0
      do iline = 1,2
c       print *
       nq=1
       if(iline==1) then
         qinit= (/0d0,0d0,0d0/)
         qend = (/0.25d0,0.25d0,0.25d0/) 
       elseif(iline==2) then
         qinit = (/.25d0,.25d0,-.25d0/)
         qend  = (/.5,.5,-.5/) 
c       elseif(iline==3) then
cc         qinit  = qend
c         qend  =  (/0d0,0d0,0d0/)
c       elseif(iline==4) then
c         qinit  = qend
c         qend =  (/-0.25,-0.25,.75/)
c       elseif(iline==5) then
c         qinit = (/0d0,0d0,0d0/)
c         qend  = (/.5,.5,-.5/)
       endif
       dx = sqrt(sum((qend-qinit)**2))/dble(nq)
       do iq= 0,nq
        iqx=iqx+1
        if(iq/=0) x=x+dx
        q= qinit + (qend-qinit) *(iq)/dble(nq)
c      iffi=30
c      open(iffi,file='qlist')
c      do  ixx=1,1 !dummy
c      do 
c        read(iffi,*,end=1010) q

        jpp = jppf( 2*pi*q, j1p,j2)
        jpm = jpmf( 2*pi*q, j1m,j2)

c matrix element -----------
        nmbas=2
        allocate(x0matinv0(nmbas,nmbas))
        phase = exp( 2*pi*img*sum(q*(/0d0,-0.5d0,-0.5d0/)) )
! 0d0,-0.5d0,-0.5d0 is R(site2) - R(site1). See LMTO file.
!
        x0matinv0(1,1) = jpp-h0
        x0matinv0(1,2) = jpm * dconjg( phase )
        x0matinv0(2,1) = jpm * phase 
        x0matinv0(2,2) = jpp-h0
        meff = mmom
        x0matinv0 = 2*x0matinv0/meff**2
        omega = meff*sqrt(
     &   x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1))

c----------------------
c        ifmat= 1017
c        open(ifmat, file="JmatTest")
c        write(ifmat,"(3d18.10, 3x, 255d18.10)") 
c     &        q, ((x0matinv0(ix,iy),ix=1,nmbas),iy=1,nmbas) 
c        write(6,"(' q=',3f9.5,' omg=',2f13.5)") q, omega !,omega2
        
        if(iqx==2) d1=omega
        if(iqx==3) d2=omega
        if(iqx==4) d3=omega
       enddo
      enddo
 1010 end


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
c      real(8) function jpmf(q,j1,j2)
      complex(8) function jpmf(q,j1,j2)
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
      call setNN(a,nn)
      do i=1,nn
         if( abs( sum(a(:,i)) )>1d-8) then !not in 111 plane
           jpmf = jpmf + j1* exp ( img* sum (a(1:3,i)*q(1:3)) )
         endif
      enddo
c      write(6,"(' q=',3f9.5,' jpmf=',5f13.5)") q,jpmf
      end
c---------------------
!      real(8) function jppf(q,j1,j2)
      complex(8) function jppf(q,j1,j2)
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
c      write(6,"(' q=',3f9.5,' jppf=',5f13.5)") q,jppf
      end
