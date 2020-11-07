      program gaunt
C- Calculate core-core exchange energy.
      implicit none
c takao obtain cg upto 2*l
c     double precision cg( (2*nl-1)**2,(2*nl-1)**2, 4*(nl-1) )
      integer:: jmax,ngaunt,ngautx,jm1,jm2,jm3,jm4,m1,m2,m3,m4
     & ,j1,j2,j3,j4,ll,ngau,ngau2,ngaut
      integer,parameter:: lmax=3, nl=lmax+1
      double precision cg( (2*nl-1)**2,(2*nl-1)**2, 0:4*(nl-1) )
     &     , cgau(0:2*(nl-1),nl**2, nl**2 )
      
c- CG coeff. generation
      jmax=  2*(nl-1)
      call clebsh_t( cg,jmax )
      print *,' end of clebsh'
c---------------------------------------------------------------
c     make Gaunt coefficient. Def. by eq(A2), and Rose eq.4.34
c     cgau(k,L',L)  ; k=j2 L'=jm3 L=jm1
c     cgau(k,j3m3,j1m1)=(-1)**(m1-m3)*cgau(k,j1m1,j3m3)
c
c  Gaunt(or CG) coefficients are not zero when the conditions 1,2,
c     and 3 are hold.
c  1.   (-1)**(j1+j2+j3).eq.1.
c  2.   abs(m1-m3).le.j2
c  3.   j2.ge.abs(j3-j1).and.j2.le.(j3+j1)
c
      ngaut=0
      ngautx=0
      do 10 j3=0,nl-1
      do 10 j1=0,nl-1
c takao  max of j2 is 2*lmx, because max of j3 and j1 is lmx. cgau(j2,jm3,jm1)
      do 10 j2=0,2*(nl-1)
      do 10 m1=-j1,j1
      do 10 m3=-j3,j3
        jm1= j1**2+m1+j1+1
        jm3= j3**2+m3+j3+1
        m2 = m3 - m1
        jm2= j2**2+m2+j2+1

        cgau(j2,jm3,jm1)=0.0d0
        if(abs(m2).le.j2. and.
     &    mod(j1+j2+j3,2).eq.0. and.
     &    j2.ge.abs(j3-j1).and. j2.le.(j3+j1) ) then
          ngaut=ngaut+1
cccc
          if( j2.le.nl-1) ngautx=ngautx+1
ccc
          cgau(j2,jm3,jm1)=cg(jm1,jm2,j3)
     &      *cg(j1**2+j1+1,j2**2+j2+1,j3)
     &      *sqrt( (2.0d0*j1+1.0d0)/(2.0d0*j3+1.0d0) )
        endif
   10 continue
      PRINT *,' * Gaunt coef. end;  num of Gaunt; nl  ='
     & , ngaut ,nl
      print *,' ngautx=',ngautx

c$$$c----check write--------------------------------------c
c$$$      if(iprint().ge.120) then
        ngau=0
        ngau2=0
        do 31 j3=0,nl-1
        do 31 j1=0,nl-1
        do 31 j2=0,2*(nl-1)
        do 31 m1=-j1,j1
        do 31 m3=-j3,j3
          jm1= j1**2+m1+j1+1
          jm3= j3**2+m3+j3+1
          m2 = m3 - m1
          jm2= j2**2+m2+j2+1
          if(jm3.ge.jm1) then
          if( (-1)**(j1+j2+j3).eq.1.and.abs(m1-m3).le.j2.and.
     &            j2.ge.abs(j3-j1).and.j2.le.(j3+j1) ) then
ccccccccccccccccccccccccccc             
            if(j3==3.and.j1==3) then
               write(6,119)j2,j3,m3,j1,m1,
     &        cgau(j2,jm3,jm1)
c     &        cgau(j2,jm1,jm3)*(-1)**(m1-m3),(-1)**(j1+j2+j3)
              ngau=ngau+1
            endif
cccccccccccccccccccccccccc           
         else
            if(abs(cgau(j2,jm1,jm3))>1d-8) then
               stop 'unexpected: something strange?'
            endif
            ngau2=ngau2+1
          endif
          endif
  119     format('   gaunt j2 j3m3 j1m1 parity= '
     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
  129     format('                   gaunt j2 j3m3 j1m1 parity= '
     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
   31   continue
c$$$  endif
        print *,' ngau=',ngau
        end
      
        subroutine clebsh_t(cg,j1mx)
c- generate crebsh gordon coefficient
c  takao 1993 8/17
c     eq.3.18 Rose's book , 'elementary theory of the angular momentum
c
c    C(j1,j2,j; m1,m2,m1+m2) = cg(jm1,jm2,j)
c      jm1=j1**2+(j1+1+m1)
c
c    C(j1,j2,j;0,0,0) at j1+j2+j=even is exactly 0  (eq.3.22),
c      however this routine gives the value oder 1.0d-16, according to
c      numerical cancellation is not complete.
c----------------------------------------------------------------c
c      implicit double precision (a-h,o-z)
      implicit none
      integer j1mx,ll
      integer
     &    j1,j2,jm1,jm2,m1,m2,m3,j3,nu,
     &    id1,id2,id3
      double precision
     &    k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
      double precision
     &    fac,fac2, igan,
     &    cg( (j1mx+1)**2, (j1mx+1)**2, 0:2*j1mx)
      external igan,ll
      print *, ' go into clebsh j1mx=',j1mx
      do 403   jm1=1, (j1mx+1)**2
        do 405 jm2=1, (j1mx+1)**2
          j1 = ll(jm1)
          m1 = jm1-(j1**2+j1+1)
          j2 = ll(jm2)
          m2 = jm2-(j2**2+j2+1)
          m3  =  m1+m2

          do 303 j3=0, 2*j1mx

c            write(6,309) j1,m1,jm1,j2,m2,jm2,j3,m3
c  309       format(' j1 m1 j1m1=',3i4,'  j2 m2 jm2',3i4,'  j3 m3=',2i4)
            cg(jm1,jm2,j3)=0.0d0
cccc
            if( j3.gt.j1+j2.or.j3.lt.abs(j1-j2) )  go to 303
            if(abs(m3).gt.j3) goto 303
c              write(6,*) '  goto calculation'
cccc
            k1= igan(j3+j1-j2)
            k2= igan(j3-j1+j2)
            k3= igan(j1+j2-j3)
            k4= igan(j3+m3)
            k5= igan(j3-m3)

            k6= igan(j1+j2+j3+1)
            k7= igan(j1-m1)
            k8= igan(j1+m1)
            k9= igan(j2-m2)
            k10=igan(j2+m2)

            fac2 =  k6*k7*k8*k9*k10
c            if(fac2.eq.0.0d0) call rx( ' k6k7k8k9k10=fac2=0')
            fac  = sqrt( (2*j3+1) *k1*k2*k3*k4*k5 /fac2 )

            do 36 nu =0, j2+j3+m1
              id1=j3-j1+j2-nu
              id2=j3+m3-nu
              id3=nu+j1-j2-m3
              if(id1.ge.0.and.id2.ge.0.and.id3.ge.0.and.
     &          j2+j3+m1-nu.ge.0.and. j1-m1+nu.ge.0  )  then
c
                k1=igan(j2+j3+m1-nu)
                k2=igan(j1-m1+nu)
                k3=igan(nu)
                k4=igan(id1)
                k5=igan(id2)
                k6=igan(id3)
                fac2= k3*k4*k5*k6
                cg(jm1,jm2,j3)=cg(jm1,jm2,j3)  
     &            + fac*(-1)**(nu+j2+m2)
     &              *k1*k2 / fac2
              endif
   36       continue
cccccccccccccccccccccccccccccccccccccccccccccccc
c--test write
c      write(105,105) j1,m1,j2,m2,j3,m3,cg(jm1,jm2,j3)
c  105 format( ' ###  j1m1 j2m2 jm = ',2i3,' ',2i3,' ',2i3, d23.16)
ccccccccccccccccccccccccccccccccccccccccccccccccc
cx            endif
  303     continue
  405   continue
  403 continue
      end
      
      double precision function igan(i)
      integer:: i,ix
      igan  = 1d0
      do ix =1,i
        igan=igan*ix
      enddo
      end

      integer function ll(ilm)
C- Returns l, given lm index
C     implicit none
C Passed parameters
      integer ilm
      integer lmaxx
      parameter (lmaxx=17)
      integer lla(lmaxx**2)
      data lla/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,
     .21*10,23*11,25*12,27*13,29*14,31*15,33*16 /
      ll = lla(ilm)
      end
