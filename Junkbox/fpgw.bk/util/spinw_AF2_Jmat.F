! This generate the exchange interaction J for given J1 J2 J3(=
! Need qlist file
      implicit none
      integer(4):: iq,nq,iline,ifmat,ix,iy,nmbas,iffi,ixx
      real (8):: j1,j2,j3,jpp0,jpp,jpm,h0,q(3)
     & ,qinit(3),qend(3),jpm0,omega,mmom,pi,x,dx,meff,omega2
! ,jpmf
      complex(8):: jpmf,jppf,phase,img=(0d0,1d0)
      complex(8),allocatable:: x0matinv0(:,:)
      pi =4d0*atan(1d0)
c      print *,pi

C     These J are defined with H = 1/2 sum J'_ij ... = sum J_ij
      j1= -13.2   !Solovyev  PRB 58, 15496
      j2= -23.5   !97.86 meV = 7.195 mRy at .25 .25 .25  Z point
      j3= j1
      mmom = 4.50d0

C     j1 = J connecting NN,  d=1/sqrt(2) type A coupled to B  ! J1-
C     j2 = J connecting 2NN, d=1, type A coupled to B
C     j3 = J connecting NN,  d=1/sqrt(2) type A coupled to A  ! J1+ FM alligend 
C      j1=1d0/8 * 1
C      j2=1d0/8 * 2
C      j3=1d0/8 * 3
C      mmom = 1d0


C      j1=0d0
C      j2=0d0
C      j3=1d0
C      mmom = 0.5d0


c      j1= -30
c      j2= -30  ! 160.0

c      j1= -9.8  /11.605*2.5**2 ! Pask
c      j2= -24.5 /11.605*2.5**2  

c      j1=-0.77* 2.5**2 !Kohgi
c      j2=-0.89* 2.5**2

C     These J are defined with H = 1/2 sum J_ij ...
      print *, 'scale J by 1/2 / 1000 /13.6 ...'
      j1 = j1/2 / 1000 / 13.6d0
      j2 = j2/2 / 1000 / 13.6d0
      j3 = j3/2 / 1000 / 13.6d0

      print "('j1,j2,j3 (mRy)=',3f12.6,' mom=',f12.6)",
     .  j1*1000,j2*1000,j3*1000,mmom

      jpp0 = jppf((/0d0,0d0,0d0/),j3,j2)  !NN   intra  j2 not used
      jpm0 = jpmf((/0d0,0d0,0d0/),j1,j2)  !NNN  inter
      h0 =  jpp0-jpm0

      x=0

c      do iline = 1,2
c       nq=30
c       if(iline==1) then
c         qinit= (/0d0,0d0,0d0/)
c         qend = (/0.25d0,0.25d0,0.25d0/) 
c       elseif(iline==2) then
c         qinit = qend
c         qend  = (/.5,.5,-.5/) 
c       elseif(iline==3) then
cc         qinit  = qend
c         qend  =  (/0d0,0d0,0d0/)
c       elseif(iline==4) then
c         qinit  = qend
c         qend =  (/-0.25,-0.25,.75/)
c       elseif(iline==5) then
c         qinit = (/0d0,0d0,0d0/)
c         qend  = (/.5,.5,-.5/)
c       endif
c       dx = sqrt(sum((qend-qinit)**2))/dble(nq)
c       do iq= 0,nq
c        if(iq/=0) x=x+dx
c        q= qinit + (qend-qinit) *(iq)/dble(nq)
      iffi=30
      open(iffi,file='qlist')

C     print *, '!! removing phase ...'

      do  ixx=1,1 !dummy
      do 
        read(iffi,*,end=1010) q

        jpp = jppf( 2*pi*q, j3,j2)
        jpm = jpmf( 2*pi*q, j1,j2)

c matrix element -----------
        nmbas=2
        allocate(x0matinv0(nmbas,nmbas))
        phase = exp( 2*pi*img*sum(q*(/0d0,-0.5d0,-0.5d0/)) )
! 0d0,-0.5d0,-0.5d0 is R(site2) - R(site1). See LMTO file.
!
C!!     phase = 1

        x0matinv0(1,1) = jpp-h0
        x0matinv0(1,2) = jpm * dconjg( phase )
        x0matinv0(2,1) = jpm * phase 
        x0matinv0(2,2) = jpp-h0
        meff = mmom
C!! PRB 58, 15496
C!!     Scaling between Jij(takao) S_i S_j -> Jij(mark) e_i e_j
        x0matinv0 = x0matinv0/(meff/2)**2
C!!     To match Takao's susceptibility calculation, need scale A-B by (-1)
        x0matinv0(1,2) = -x0matinv0(1,2)
        x0matinv0(2,1) = -x0matinv0(2,1)
C!!     Spin wave spectra in when J is in units of S
C       When J is in units of e, need to scale by 4/m_i/m_j
        omega = 1*meff*sqrt(
     &   x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1))
c----------------------
        ifmat= 1017
        open(ifmat, file="JmatTest")
        write(ifmat,"(3d18.10, 3x, 255d18.10)") 
     &        q, ((x0matinv0(ix,iy),ix=1,nmbas),iy=1,nmbas) 
c        omega2 = 2d0/mmom *( (jpp-h0)**2 - jpm**2 )**0.5d0
c        write(6,"(f9.5, ' q=',3f9.5,' omg=',2f13.5)") x,q, omega !,omega2
        write(6,"(' q=',3f9.5,' omg=',2f13.5)") q, omega !,omega2

        deallocate(x0matinv0)

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

C      do  j = 1, nnn
C        print *, j, sngl (a(:,j) - (/0d0,-0.5d0,-0.5d0/))
C      enddo
C      stop

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
! A sublattice coupled to B sublattice
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
!      real(8) function jppf(q,j3,j2)
      complex(8) function jppf(q,j3,j2)
! A sublattice coupled to A sublattice
      implicit none
      integer(4):: i,nn
      complex(8):: img=(0d0,1d0)
      real(8) :: a(3,12),pi,q(3),j3,j2
      call setNN(a,nn)
      jppf=0d0
      do i=1,nn
         if( abs( sum(a(:,i)) )<1d-8) then  !in plane
           jppf = jppf + j3* exp( img* sum(a(1:3,i)*q(1:3)) )
         endif
      enddo
c      write(6,"(' q=',3f9.5,' jppf=',5f13.5)") q,jppf
      end
