c--------------------------------------------------------------------
      subroutine rotcg(lmxax,symops,ng,cgr)
c --- rotated CG coefficients.
      implicit none
      integer(4) :: lmxax, ng, nlmxa,
     &              lnjcg, lnxcg, 
     &              ilma,la,ilmb,lh,ii,indx,icg1,icg2,icg,
     &              ig,lm1,lm2,lm,l,m,md,lmd,lmxcg,ilm ,lx
      real(8) :: 
     &   cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng), 
     &   symops(9,ng)        ,sumr
      real(8),allocatable:: cg(:),dlmm(:,:,:,:),cgn(:,:,:)
      integer(4),allocatable :: jcg(:),indxcg(:)
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4):: ll
#endif
c --- CG coefficienets. <LM3|lm1 lm2>
c inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
c Injcg = indxcg(inxcg) to indxcg(inxcg)-1
c cg(inxcg)  : = <lm3|lm1 lm2>
c jcg(lnjcg) : = lm3
      print *,' rotcg:'
c      do ig=1,ng
c      write(6,*)' transposed symope ig  =',ig
c     write(6,'(3f12.6)') symops(1:3,ig)
c     write(6,'(3f12.6)') symops(4:6,ig)
c     write(6,'(3f12.6)') symops(7:9,ig)
c      enddo
      allocate(cgn((2*lmxax+1)**2,(lmxax+1)**2,(lmxax+1)**2) )
      cgn = 0d0
      nlmxa = (lmxax+1)**2
      lmxcg = lmxax
      if (lmxcg .le. 6) then
        lnjcg = 6500
        lnxcg = 1300
      else if (lmxcg .le. 8) then
        lnjcg = 22700
        lnxcg = 3400
      else if (lmxcg .le. 10) then
        lnjcg = 62200
        lnxcg = 7400
      else
        write(6,*) 'rotcg: cannot handle lmxcg=',lmxcg
        stop       'rotcg: cannot handle lmxcg'
      endif
      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
      call scg(lmxcg,cg,indxcg,jcg)
c----
      do ilma = 1, nlmxa
        la = ll(ilma)
        do ilmb = 1, nlmxa
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do icg = icg1, icg2
            ilm  = jcg(icg)
            cgn(ilm, ilma,ilmb)  = cg(icg) ! ilm is move to 1st argument.!
          enddo
        enddo
      enddo
ccccccccccccccccccccccccccccccccc
c      print *,' cgn=',cgn(1,1,1),1/sqrt(4*3.1415926) !=1d0
c      stop
ccccccccccccccccccccccccccccccccc

c --- Rotation matrix
      allocate(
     & dlmm( -2*lmxax:2*lmxax, -2*lmxax:2*lmxax, 0:2*lmxax,ng))
      call rotdlmm(symops, ng, 2*lmxax+1,dlmm)
c --- Rotated CG
      do lm =  1, (2*lmxax+1)**2
        l = ll(lm)
        m = lm - l**2 -l -1
        lx = l**2 +l +1
        do ig  = 1, ng
          do lm2 = 1, nlmxa
            do lm1 = 1, nlmxa
! ilm is move to 1st argument.!
              cgr(lm1,lm2,lm,ig) = sum(cgn(lx-l:lx+l,lm1,lm2)*dlmm(-l:l,m,l,ig))
c        sumr = 0d0
c        do md = -l,l
c          lmd = l**2 +l +1 + md
c          sumr = sumr + cgn(lm1,lm2,lmd)*dlmm(md,m,l,ig)
c        enddo
c        cgr(lm1,lm2,lm,ig) = sumr
            enddo
          enddo
        enddo
      enddo
      deallocate( cg,dlmm,cgn,jcg,indxcg)
      print *,' rotcg end:'
      end

      subroutine rotdlmm(symops,ng,nl ,dlmm)
c- Generate rotation matrix D^l_{m,m'} for L-representaiton, corresponding
c  to points group operations.
Ci symops(9,ng),ng; point ops.
Ci nl; num.of l =lmax+1
Co dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2); D^l_{m,m'}. Indexes are for Real harmonics.
Cr   dlmmc is used as work area about 200kbyte used for  s,p,d,f -> nl=4
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ng,nl
      double precision SYMOPS(9,ng),
     .        am(3,3) ,fac1,fac2
      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)

      double complex   dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
      double precision det,igann,osq2
      double complex   msc(0:1,2,2), mcs(0:1,2,2),Img
     &            ,dum(2)
      parameter(Img=(0d0,1d0))

      print *; print *,' ROTDLMM:'
      do 10 ig =1,ng
        do 20 i=1,3
        do 20 j=1,3
          am(i,j) = symops(i+3*(j-1),ig)
   20   continue
c calculate determinant(signature)
        det= am(1,1)*am(2,2)*am(3,3)
     &        -am(1,1)*am(3,2)*am(2,3)
     &        -am(2,1)*am(1,2)*am(3,3)
     &        +am(2,1)*am(3,2)*am(1,3)
     &        +am(3,1)*am(1,2)*am(2,3)
     &        -am(3,1)*am(2,2)*am(1,3)
        if(abs(abs(det)-1d0).ge.1d-10) then
          print *,' rotdlmm: det/=1 ig and det=',ig,det
          stop
        endif
c seek Euler angle   print *,' goto cbeta',ig,det
        cbeta = am(3,3)/det
c added region correction so as to go beyond domain error for functions, dsqrt and acos.
        if(abs(cbeta-1d0).le.1d-6) cbeta= 1d0
        if(abs(cbeta+1d0).le.1d-6) cbeta=-1d0
        beta = dacos(cbeta)
        sbeta= sin(beta)
c beta= 0~pi
        if(sbeta.le.1.0d-6) then
          calpha= 1d0
          salpha= 0d0
          alpha = 0d0
          cgamma= am(2,2)/det
          sgamma= am(2,1)/det
        else
          salpha =  am(2,3)/sbeta/det
          calpha =  am(1,3)/sbeta/det
          sgamma =  am(3,2)/sbeta/det
          cgamma = -am(3,1)/sbeta/det
        endif
        co2 = dcos(beta/2)
        so2 = dsin(beta/2)
c         print *,' calpha=',calpha
        if(abs(calpha-1.0d0).le.1.0d-6) calpha= 1.0d0
        if(abs(calpha+1.0d0).le.1.0d-6) calpha=-1.0d0
        if(abs(cgamma-1.0d0).le.1.0d-6) cgamma= 1.0d0
        if(abs(cgamma+1.0d0).le.1.0d-6) cgamma=-1.0d0
        alpha=dacos(calpha)
        if(salpha.lt.0d0) alpha=-alpha
        gamma=dacos(cgamma)
        if(sgamma.lt.0d0) gamma=-gamma
c         print *,'alpha beta gamma det=',alpha,beta,gamma,det
        do 30 l =  0, nl-1
        do 30 md= -l, l
        do 30 m = -l, l
c  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wisley and Sons.  p.52 (4.13)
          fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
          fac2 = 0d0
          do 40 ikap=0,2*l
            if(l-md-ikap.ge.0 .and. l+m-ikap.ge.0
     &          .and.ikap+md-m.ge.0) then
              add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap)
     &           *igann(ikap+md-m)*igann(ikap) )
              if(2*l+m-md-2*ikap.ne.0) add=add*co2**(2*l+m-md-2*ikap)
              if(md-m+2*ikap.ne.0)     add=add*(-so2)**(md-m+2*ikap)
              fac2 = fac2+add
            endif
   40     continue
c l-th rep. is odd or even according to (det)**l
          dlmmc(md,m,l,ig) = fac1*fac2*det**l*
     &        cdexp( -Img*(alpha*md+gamma*m) )
   30   continue

        am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
        am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
        am(1,3)= sin(beta)*cos(alpha)
        am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
        am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
        am(2,3)= sin(beta)*sin(alpha)
        am(3,1)=-sin(beta)*cos(gamma)
        am(3,2)= sin(beta)*sin(gamma)
        am(3,3)= cos(beta)

        if(abs(am(1,1)*det-symops(1,ig)).gt.1.0d-8.or.
     &    abs(am(2,1)*det-symops(2,ig)).gt.1.0d-8.or.
     &    abs(am(3,1)*det-symops(3,ig)).gt.1.0d-8.or.
     &    abs(am(1,2)*det-symops(4,ig)).gt.1.0d-8.or.
     &    abs(am(2,2)*det-symops(5,ig)).gt.1.0d-8.or.
     &    abs(am(3,2)*det-symops(6,ig)).gt.1.0d-8.or.
     &    abs(am(1,3)*det-symops(7,ig)).gt.1.0d-8.or.
     &    abs(am(2,3)*det-symops(8,ig)).gt.1.0d-8.or.
     &    abs(am(3,3)*det-symops(9,ig)).gt.1.0d-8) then
          print *,' rotdlmm: not agree. symgrp and one by eular angle'
          stop
        endif
cccccccccccccccccccccccc
c        if(iprint().ge.140) then
        if(.false.) then
          print *;print *;print *,' **** group ops no. ig=', ig
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
          print *,' by Eular angle '
          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
        endif
 1731   format (' ',3f9.4)
cccccccccccccccccccccccc
   10 continue
c conversion to cubic rep. Belows are from csconvs
c  msc mcs conversion matrix generation 2->m 1->-m for m>0
      osq2 = 1d0/sqrt(2d0)
      do m = 0,1
        Msc(m,1,1)= osq2*(-1)**m
        Msc(m,1,2)=-osq2*Img*(-1)**m
        Msc(m,2,1)= osq2
        Msc(m,2,2)= osq2*Img

        Mcs(m,1,1)= osq2*(-1)**m
        Mcs(m,1,2)= osq2
        Mcs(m,2,1)= osq2*Img*(-1)**m
        Mcs(m,2,2)=-osq2*Img
      enddo
c
      print * ,' goto do 23'
      do 23 is=1,ng
        if(.false.) then
c        if(iprint().ge.150) then
          print *; print *,' **** group ops no. ig=', is
          write(6,1731) symops(1,is),symops(4,is),symops(7,is)
          write(6,1731) symops(2,is),symops(5,is),symops(8,is)
          write(6,1731) symops(3,is),symops(6,is),symops(9,is)
        endif
c convert to cubic rep.
      do 23   l =0,nl-1
        do 33 m2=-l,l
        do 33 m1= 1,l
          dum(1)= dlmmc(m2, m1,l,is)
          dum(2)= dlmmc(m2,-m1,l,is)
          mx    = mod(m1,2)
          dlmmc(m2,  m1,l,is)=
     &                       dum(1)*msc(mx,1,1)
     &                      +dum(2)*msc(mx,2,1)
          dlmmc(m2, -m1,l,is)=
     &                       dum(1)*msc(mx,1,2)
     &                      +dum(2)*msc(mx,2,2)
   33   continue
        do 43 m2=  1,l
        do 43 m1= -l,l
          dum(1)=dlmmc( m2, m1,l,is)
          dum(2)=dlmmc(-m2, m1,l,is)
          mx=mod(m2,2)
          dlmmc( m2, m1,l,is)=
     &                       mcs(mx,1,1)*dum(1)
     &                      +mcs(mx,1,2)*dum(2)
          dlmmc(-m2, m1,l,is)=
     &                       mcs(mx,2,1)*dum(1)
     &                      +mcs(mx,2,2)*dum(2)
   43   continue
        do 53 m2=-l,l
        do 53 m1=-l,l
          dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
          if( abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12 ) stop
     &      ' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12'
   53   continue
ccccccccccccccccccccc
        if(.false.) then
c        if(.true.) then
c        if(iprint().ge.41) then
          print *; print *,'  points ops  ig, l=', is,l,' cubic   '
          do m2=-l,l
            write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
c    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
          enddo
        endif
cccccccccccccccccccccc
   23 continue
      print *,' end of rotdlmm'
      end
c--------------------------------------------
      double precision function igann(i)
      igann  = 1d0
      do ix =1,i
        igann=igann*ix
      enddo
      end
