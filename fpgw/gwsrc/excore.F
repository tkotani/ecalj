      subroutine excore(nrmx,nl,nnc,nclass,nspn,nbas, 
     &  phic,nindxc,iclass,
     &  a,b,nr,rofi)
C- Calculate core-core exchange energy.
      implicit none
      integer(4):: nrmx,nnc,nclass,nspn,ncmx,nl,nbas
      integer(4):: nindxc(0:nl-1,nclass),nr(nclass),iclass(nclass)
      real(8) :: a(nclass),b(nclass),rofi(nrmx,nclass)
      real(8) :: phic (nrmx,0:nl-1,nnc,nclass,nspn),wgtx

      real(8):: exacct,exacc(nbas,nspn),rydberg
c-local
      integer(4),parameter :: nsxamx=10 000 000
      real(8),allocatable :: 
     & gc(:,:,:,:), sxadata(:),sigkcc(:,:,:,:,:)
      integer(4),allocatable:: indxsxa(:,:),lcr(:,:,:),ncr(:,:)
      integer(4):: ir,l,n,ic,ncrr,ncmxx(nclass),kmax,k,it,isp,
     & isx,lm1,lm2,lm3,lm4,ibas,icr,icrd,nsxatot,ifexcore
      real(8):: rkp(nrmx,0:2*(nl-1)),rkm(nrmx,0:2*(nl-1))
#ifdef COMMONLL
      integer(4) ll(51**2)
      common/llblock/ll
#else
      integer(4):: ll
      external ll
#endif

c-----------------------------------------------------------
c ncmx
      do ic = 1,nclass
        ncmxx(ic) = sum(nindxc(0:nl-1,ic))
      enddo
      ncmx = maxval(ncmxx(1:nclass))

      print *,' ncmx nl nspn=',ncmx,nl,nspn

c --- convert format of core function ---
c nindxc(l,ic) phic ---> ncr lcr gc
c
      allocate( gc(nrmx,ncmx,nspn,nclass)
     & ,lcr(ncmx,nspn,nclass),ncr(nspn,nclass) )
      do ic = 1,nclass
        ncrr = 0
        do  l = 0,nl-1
          do  n = 1,nindxc(l,ic)
            print *, ' l n nindx=',l,n,nindxc(l,ic)
            ncrr = ncrr + 1
            do isp=1,nspn
              gc(1:nrmx,ncrr,isp,ic) = phic(1:nrmx,l,n,ic,isp)
            enddo
            lcr(ncrr,1:nspn,ic) = l
          end do
        end do
        ncr(1:nspn,ic) = ncrr
      enddo
c---------------------------------------------------------------
      print *,' goto alloc'
      allocate( sxadata(nsxamx),indxsxa(6,nsxamx),
     &  sigkcc(0:2*(nl-1), ncmx, ncmx, nclass,nspn) )
      print *,' end of alloc'
c
      kmax  = 2*(nl-1)
      do ic = 1,nclass
        print *,' make rkp rkm ic=',ic
c- rkp,rkm ----------- This is from subroutine bess
        do k=0, kmax
          rkp(1,k)=0d0
          rkm(1,k)=0d0
          do ir=2,nr(ic)
            rkp(ir,k) = rofi(ir,ic)**k
            rkm(ir,k) = rofi(ir,ic)**(-k-1)
          enddo
        enddo
        print *,' end of make rkp rkm ic=',ic
c- radial integrals
        do isp = 1,nspn
          call intsigkcc(sigkcc(0,1,1,ic,isp), ncmx,
     &      gc(1,1,isp,ic), 
     &      a(ic),b(ic),nr(ic),rofi(1,ic),nl,
     &      ncr(isp,ic),lcr(1,isp,ic),
     &      kmax, rkp, rkm )
        enddo
      enddo

c- make structure spherical part --- sxadata
      call mksxa(nl,nsxamx,
     o   sxadata,indxsxa,nsxatot)

c- core-core exchange =sxadata * radial integral
      print * !; print *,' go into EXEX.CORE part'
      exacc = 0d0
      do 300 isx=1,nsxatot
c	print *,' do 300 isx isp=',isx, isp
      do 300 isp=1,nspn
        k   = indxsxa (1,isx)
        lm1 = indxsxa (2,isx)
        lm3 = indxsxa (3,isx)
        lm2 = indxsxa (4,isx)
        lm4 = indxsxa (5,isx)
c	print *,' do 300 isx isp= xxx'
        if(lm3.ne.lm2 .or. lm1.ne.lm4) cycle
c- icr->l3l2 icrd->l1l4
        do 350 ibas = 1,nbas
          ic = iclass(ibas)
c	  print *,' ibas ic=',nbas,ibas,ic
          do icr = 1, ncr(isp,ic)
            do icrd= 1, ncr(isp,ic)
ccccccccccccccccccccccc
c	print *,' isx=',isx,ic,isp,icr,icrd,lm1,lm3
c	print *,' lcr=',lcr(icr, isp,ic),lcr(icrd, isp,ic)
c	print *,' ll=',ll(lm1),ll(lm3)
ccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccc
c     test
c            if(ll(lm3)==1.or.ll(lm1)==1) cycle
ccccccccccccccccccccccccccccccccccccc


              if( lcr(icr, isp,ic) /= ll(lm3)) cycle
              if( lcr(icrd,isp,ic) /= ll(lm1)) cycle
              exacc(ibas,isp) = exacc(ibas,isp) 
     &       - sxadata(isx) * sigkcc(k,icr,icrd,ic,isp)
            enddo
          enddo
  350   continue
c---------------------
  300 continue
      wgtx=1d0
      if(nspn==1) wgtx=2d0
      exacct = sum(exacc(1:nbas,1:nspn))*wgtx

      ifexcore=9051
      open (ifexcore,file='TEEXXcc')
      write(6,*) '==== EXCORE ==> TEEXXcc ============'
      write(ifexcore,*) '======================================='
      write(ifexcore,*) '  Exchange energy core-core   Exx (eV) '
      write(ifexcore,*) '======================================='
      write(ifexcore,*) ' *** '
      write(ifexcore,
     &   "( f20.10,2i4,' ! Ex core-core total (eV) --  nbas nsp ')")
     &   exacct*rydberg(), nbas, nspn
      write(6,"(' excore total =',f13.6)") exacct*rydberg()
      do isp =1,nspn
        do ibas=1,nbas
          write(6,"(' ibas isp=',2i3,' Ex core-core (eV) =',f13.6)")
     &  ibas,isp,exacc(ibas,isp)*rydberg()
c	write(ifexcore,"(' ibas isp=',2i3,' Ex core-core (Ry) =',f13.6)")
c     &  ibas,isp,exacc(ibas,isp)
          write(ifexcore,"( f20.10,2i4,' ! Ex core-core (eV)  ibas isp ')")
     &  exacc(ibas,isp)*rydberg(),ibas,isp
        enddo
      enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine intsigkcc(sigkcc,ncmx,
     &         gc, a,b,nr,rofi, nl,
     &         ncr,lcr ,
     &         kmax, rkp, rkm )
C- integral \sigma^k(l1,l3,ip1,ip3,l2,l4,ip2,ip4) in the Phys.Rev.B34, 5512 ,Simpson
C---------------------------------------------------------------------
Ci ncmx,nrec, gl,gpl,gc , a,b,nr,rofi,nl,ncr,lcr
Co  sigkcc
C---------------------------------------------------------------------
      implicit none
      integer nr, ncmx, nl,l1,l2,l3,l4,l1l3,l2l4,k
     &   ,ip1,ip3,ip2,ip4,ip2ip4,icr,l1l4,icrd,ip1ip3,ir,lr0
     &   ,ncr,lcr(ncmx)
      double precision a,b,rofi(nr)
     &   , gc(nr,ncmx)
     &   , sum1,sum2   
     &   , sigkcc(0:2*(nl-1), ncmx, ncmx  ) 
     &   , int1(nr),int2(nr),a1(nr),a2(nr),b1(nr),f13(nr)
      integer kmax
      double precision rkp(nr,0:kmax),rkm(nr,0:kmax)
      print *,' intsigkcc:'
c--------------
      print *,' normcheck:'
      do icr  = 1,ncr; l2 =lcr(icr);  l3=l2
        do icrd = 1,ncr; l1 =lcr(icrd); l4=l1
          call gintxx(gc(1,icr),gc(1,icrd),a,b,nr,sum1)
          write(6,"(' norm check ='2i3,2x,f13.6)") icr,icrd,sum1
        enddo
      enddo
c      do ir=1,nr
c        write(6,"(' rofi ='i3,2x,f13.6)") ir,rofi(ir)
c      enddo
c
c ---------------between core and core->  sigkc
      do 2215 icr  = 1,ncr; l2 =lcr(icr);  l3=l2
      do 2215 icrd = 1,ncr; l1 =lcr(icrd); l4=l1
        do 2210 k  = abs(l1-l3),l1+l3
          if( k<abs(l2-l4) .or. k>l2+l4 ) cycle
          if(mod(k+l1+l3,2)==1 .or. mod(k+l2+l4,2)==1) cycle
          a1(1)    = 0d0
          a1(2:nr) = rkp(2:nr,k) *gc(2:nr,icrd)
          a2(1)    = 0d0
          a2(2:nr) = rkm(2:nr,k) *gc(2:nr,icrd)
          b1(1:nr) = gc(1:nr,icr)
          call intn_smp_g(a1,b1,int1,a,b,rofi,nr,lr0)
          call intn_smp_g(a2,b1,int2,a,b,rofi,nr,lr0)

c check write ---------------------------------------------
          call gintxx(a1,b1,a,b,nr,sum1)
          call gintxx(a2,b1,a,b,nr,sum2)
          write(6,"(' integral ='2d13.6)") sum1,sum2
          write(6,"(' integral ='2d13.6)") int1(1),int2(1)
c        endif
c---
          f13(1)    = 0d0
          f13(2:nr) = rkm(2:nr,k) *( int1(1)-int1(2:nr) )
     &            + rkp(2:nr,k) * int2(2:nr)
          a1(1:nr)  = gc(1:nr,icr) *f13(1:nr)
          b1(1:nr)  = gc(1:nr, icrd)
          call gintxx(a1,b1,A,B,NR, sigkcc(k, icr, icrd ) )
c chekc write
c        if(iprint().ge.130)
          write(6,"( ' k icr icrd =',3i5,'  sigkcc=',d15.8)")  
     &    k, icr,icrd ,sigkcc(k,icr,icrd)
 2210   continue
 2215 continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mksxa(nl,nsxamx,
     o   sxadata,indxsxa,nsxatot)
C- make sx (spherical-dependent part for exx energy.) ------------c
Ci INPUT
Ci   nl
Ci   nsxamx
Co   indxsxa, sxdata, nsxatot
Ci
Cr       ngaut is number for non-zero cgau.
Cr  Generation of
Cr     cgau(k,L',L) for R=R'
C-------------------------------------------------------------------------c
      implicit none
      integer nl,nlx,
     &    m,istcl,lrd,lr,j1,j2,j3,m1,m2,m3,jm1,jm2,jm3,
     &    ngau,ngaut,ngau2,k1,k2,j4,m4,jm1m,
     &    jm2m,jm4,km1,km2,
     &    i,j,k ,jmax
      integer nsxamx, nsxatot
      integer indxsxa(6,nsxamx)
      double precision sxadata(nsxamx)
c
      double complex  
     &  sxx(-(nl-1):(nl-1),-(nl-1):(nl-1),-(nl-1):(nl-1),-(nl-1):(nl-1))
c     & ,syy(-3:3,-3:3,-3:3,-3:3)
      double precision osq2

c takao obtain cg upto 2*l
c      double precision cg( (2*nl-1)**2,(2*nl-1)**2, 4*(nl-1) )
      double precision cg( (2*nl-1)**2,(2*nl-1)**2, 0:4*(nl-1) )
     &  , cgau(0:2*(nl-1),nl**2, nl**2 ) ,dum
c      integer iprint

      external clebsh ,mscmul
      double complex msc(0:1,2,2),mcs(0:1,2,2),Img
     &          ,ap,am,mscmul
c
      integer ngautx
      data Img/(0.0d0,1.0d0)/

c----------------------------------------------------------------------
      print *,' goto mksxa'

c  msc, conversion matrix generation 1->m 2->-m for m>0
      osq2 = 1d0/sqrt(2d0) !sq2=1/sqrt(2)
      do m=0,1
        Msc(m,1,1)= osq2*(-1)**m
        Msc(m,1,2)=-osq2*Img*(-1)**m
        Msc(m,2,1)= osq2
        Msc(m,2,2)= osq2*Img

        Mcs(m,1,1)= osq2*(-1)**m
        Mcs(m,1,2)= osq2
        Mcs(m,2,1)= osq2*Img*(-1)**m
        Mcs(m,2,2)=-osq2*Img
      enddo

c- CG coeff. generation
      jmax=  2*(nl-1)
      call clebsh( cg,jmax )
      print *,' end of clebsh'
c- check write
c      if(iprint().ge.120) then
c        do 106 j1=0,nl-1
c        do 106 j2=0,nl-1
c        do 106 j3=abs(j1-j2),j1+j2
c        do 106 m1=-j1,j1
c        do 106 m2=-j2,j2
c          jm1=j1**2+j1+1+m1
c          jm2=j2**2+j2+1+m2
c          m3=m1+m2
c          if(abs(m3).gt.j3) go to 106
c          write(6,105) j1,m1,j2,m2,j3,m3,cg(jm1,jm2,j3)
c  105     format( ' ###  j1m1 j2m2 jm = ',2i3,' ',2i3,' ',2i3, d12.5)
c  106   continue
c      endif

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
c$$$        ngau=0
c$$$        ngau2=0
c$$$        do 31 j3=0,nl-1
c$$$        do 31 j1=0,nl-1
c$$$        do 31 j2=0,2*(nl-1)
c$$$        do 31 m1=-j1,j1
c$$$        do 31 m3=-j3,j3
c$$$          jm1= j1**2+m1+j1+1
c$$$          jm3= j3**2+m3+j3+1
c$$$          m2 = m3 - m1
c$$$          jm2= j2**2+m2+j2+1
c$$$          if(jm3.ge.jm1) then
c$$$          if( (-1)**(j1+j2+j3).eq.1.and.abs(m1-m3).le.j2.and.
c$$$     &      j2.ge.abs(j3-j1).and.j2.le.(j3+j1) ) then
c$$$
c$$$            write(6,119)j2,j3,m3,j1,m1, cgau(j2,jm3,jm1),
c$$$     &        cgau(j2,jm1,jm3)*(-1)**(m1-m3),(-1)**(j1+j2+j3)
c$$$            ngau=ngau+1
c$$$          else
c$$$            write(6,129)j2,j3,m3,j1,m1, cgau(j2,jm3,jm1),
c$$$     &        cgau(j2,jm1,jm3)*(-1)**(m1-m3),(-1)**(j1+j2+j3)
c$$$            ngau2=ngau2+1
c$$$          endif
c$$$          endif
c$$$  119     format('   gaunt j2 j3m3 j1m1 parity= '
c$$$     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
c$$$  129     format('                   gaunt j2 j3m3 j1m1 parity= '
c$$$     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
c$$$   31   continue
c$$$      endif

c-- sxa cal. for R=R', see eq.(A11) --------------------------------c
      print *
      print *
      print *, '  *** Go into SXA cal. for R=Rdash pair '
      nsxatot=0
c      nsxatot2=0
      do 107 k=0,2*(nl-1)
      do 107 j1=0,nl-1
      do 107 j2=0,nl-1
      do 107 j3=0,nl-1
      do 107 j4=0,nl-1
        do 207 m1=-j1,j1
        do 207 m2=-j2,j2
        do 207 m3=-j3,j3
        do 207 m4=-j4,j4
          jm1 =j1**2+j1+1+m1
          jm1m=j1**2+j1+1-m1

          jm2 =j2**2+j2+1+m2
          jm2m=j2**2+j2+1-m2

          jm3 =j3**2+j3+1+m3
          jm4 =j4**2+j4+1+m4

          sxx(m1,m2,m3,m4)=0.0d0
          if(m1+m2+m3+m4.eq.0) then
            sxx(m1,m2,m3,m4)=2*(-1)**(m1+m2)*cgau(k,jm3,jm1m)
     &           *cgau(k,jm2m,jm4)
          endif
cccccccccccccccccccccccccc
c            syy(m1,m2,m3,m4)= sxx(m1,m2,m3,m4)
cccccccccccccccccccccccccc
  207   continue

c convert to real harmonics rep.-------------
        call convsx(sxx,j1,j2,j3,j4,nl,msc,mcs)
c
        do 407 m1=-j1,j1
        do 407 m2=-j2,j2
        do 407 m3=-j3,j3
        do 407 m4=-j4,j4
          jm1 =j1**2+j1+1+m1
          jm2 =j2**2+j2+1+m2
          jm3 =j3**2+j3+1+m3
          jm4 =j4**2+j4+1+m4
c
c          if(iprint().ge.130.and.
c     &       DREAL( sxx(m1,m2,m3,m4) ).gt.1.0d-6 ) write(6,1007)
c     &      k,j1,m1,j2,m2,j3,m3,j4,m4,syy(m1,m2,m3,m4),
c     &      sxx(m1,m2,m3,m4)

          dum= DREAL( sxx(m1,m2,m3,m4) )

          if(dabs(dimag(sxx(m1,m2,m3,m4))) >1d-8)
     &      stop ' MAKSX; im part of sxx .ne. 0'
c
          if(abs(dum)> 1d-8) then
            if(nsxatot.ge.nsxamx) 
     &        stop ' MAKSX: enlarge the size of nsxamx '
            nsxatot=nsxatot+1
c note its oeder! jm1.ge.jm3
            indxsxa(1,nsxatot)=k
            indxsxa(2,nsxatot)=jm1
            indxsxa(3,nsxatot)=jm3
            indxsxa(4,nsxatot)=jm2
            indxsxa(5,nsxatot)=jm4
            sxadata(nsxatot  )=dum
          endif
c
  407   continue
  107 continue
c
      sxadata(1:nsxatot)=sxadata(1:nsxatot)/2d0
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine clebsh()
      subroutine clebsh(cg,j1mx)
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
      integer j1mx
      integer
     &    j1,j2,jm1,jm2,m1,m2,m3,j3,nu,
     &    id1,id2,id3
      double precision
     &    k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
      double precision
     &    fac,fac2, igan,
     &    cg( (j1mx+1)**2, (j1mx+1)**2, 0:2*j1mx)
#ifdef COMMONLL
      integer(4) ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
      external ll
#endif
      external igan
c
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
            if(fac2.eq.0.0d0) stop  ' k6k7k8k9k10=fac2=0'
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine convsx(sxx,j1,j2,j3,j4,nl,msc,mcs)
C- sx based on spherical har. rep. is converted to cub rep.
Ci input
Ci     sxx; in the harmonic rep.
Co output
Co     sxx; in the sperical rep.
c
Cr inversion test was included; uncomment the lines commented by ct.
Cr
C------------------------------------------------
      implicit none
c      implicit double precision(a-h,o-z)
c      parameter(nlx=3)
      integer
     &    nl,j1,j2,j3,j4,m1,m2,m3,m4
      double complex   
     &   sxx (-(nl-1):(nl-1),-(nl-1):(nl-1),
     &       -(nl-1):(nl-1),-(nl-1):(nl-1))
c     &   sxxx(-(nlx-1):(nlx-1),-(nlx-1):(nlx-1),
c     &        -(nlx-1):(nlx-1),-(nlx-1):(nlx-1)),
c     &   syyy(-(nlx-1):(nlx-1),-(nlx-1):(nlx-1),
c     &        -(nlx-1):(nlx-1),-(nlx-1):(nlx-1))

      double complex 
     &    msc(0:1,2,2),Img,ap,am,mscmul
c  inversion test
      double complex 
     &    mcs(0:1,2,2)

c      if(nlx.ne.nl) stop 'CONVSX: nlx.ne.nl'
      do 310 m1=   1, j1
      do 310 m2= -j2, j2
      do 310 m3= -j3, j3
      do 310 m4= -j4, j4
        ap=sxx( m1,m2,m3,m4)
        am=sxx(-m1,m2,m3,m4)
        sxx(  m1,m2,m3,m4) = mscmul(1,m1,msc,ap,am)
        sxx( -m1,m2,m3,m4) = mscmul(2,m1,msc,ap,am)
  310 continue
      do 320 m1= -j1, j1
      do 320 m2=   1, j2
      do 320 m3= -j3, j3
      do 320 m4= -j4, j4
        ap=sxx( m1, m2,m3,m4)
        am=sxx( m1,-m2,m3,m4)
        sxx(  m1, m2,m3,m4) = mscmul(1,m2,msc,ap,am)
        sxx(  m1,-m2,m3,m4) = mscmul(2,m2,msc,ap,am)
  320 continue
      do 330 m1= -j1, j1
      do 330 m2= -j2, j2
      do 330 m3=   1, j3
      do 330 m4= -j4, j4
        ap=sxx( m1, m2, m3,m4)
        am=sxx( m1, m2,-m3,m4)
        sxx(  m1, m2, m3,m4) = mscmul(1,m3,msc,ap,am)
        sxx(  m1, m2,-m3,m4) = mscmul(2,m3,msc,ap,am)
  330 continue
      do 340 m1= -j1, j1
      do 340 m2= -j2, j2
      do 340 m3= -j3, j3
      do 340 m4=   1, j4
        ap=sxx( m1, m2, m3, m4)
        am=sxx( m1, m2, m3,-m4)
        sxx(  m1, m2, m3, m4) = mscmul(1,m4,msc,ap,am)
        sxx(  m1, m2, m3,-m4) = mscmul(2,m4,msc,ap,am)
  340 continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double complex function mscmul(i,m,msc,ap,am)
c----multiple mat matrix
      integer mx,m,i
      double complex msc(0:1,2,2),ap,am
c  msc, conversion matrix generation 1->m 2->-m for m>0
      mx=mod(m,2)
      mscmul=ap*msc(mx,1,i)+am*msc(mx,2,i)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine intn_smp_g(g1,g2,int,a,b,rofi,nr,lr0)
c-- intergral of two wave function.
c   This is for true g1(ir) = r phi(ir), where phi is the true wave function.
c
c int(r) = \int_(r)^(rmax) u1(r') u2(r') dr'
c
c lr0 dummy index, now not used.
c simpson rule ,and with higher rule for odd devision.
c --------------------------------------------------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT none
      integer nr,ir,lr0
      double precision g1(nr),g2(nr),int(nr),a,b,rofi(nr),w1,w2,w3
     &    ,ooth,foth
      integer(4):: mx=1
C Local parameters
c one over three-> oot,  and so on.
      data ooth,foth/0.33333333333333333,1.3333333333333333333/
      data w1,w2,w3/0.41666666666666666,0.6666666666666666666,
     &              -0.083333333333333333/
C                        xxxxx
      if(mod(nr,2).eq.0) 
     &  stop ' intn_smp_g: nr should be odd. '
c
      int(1)=0.0d0
c l00 means u1(r)u2(r)~r**lr0 near r~0
cc simplest formula
cc      do 10 ir=3,nr
cc        int(ir)=int(ir-1)
cc     &         +0.5d0*G1(ir-1)*G2(ir-1)*( a*(b+rofi(ir-1)) )**mx
cc     &         +0.5d0*G1(ir)*G2(ir)*( a*(b+rofi(ir)) )**mx
cc   10 continue
c simpson rule
      DO  10  IR = 3,NR,2
        int(ir)=int(ir-2)
     &             + ooth*G1(IR-2)*G2(IR-2)*( a*(b+rofi(ir-2)) )**mx
     &             + foth*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )**mx
     &             + ooth*G1(IR)*G2(IR)*( a*(b+rofi(ir)) )**mx
   10 CONTINUE

c At the value for odd points, use the same interpolation above
      do 20 ir = 2,nr-1,2
        int(ir)=int(ir-1)
     &             + w1*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )**mx 
     &             + w2*G1(IR)  *G2(IR)*  ( a*(b+rofi(ir)  ) )**mx
     &             + w3*G1(IR+1)*G2(IR+1)*( a*(b+rofi(ir+1)) )**mx
   20 continue
c
      do 30 ir=1,nr
        int(ir)=int(nr)-int(ir)
   30 continue
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function igan(i)
      igan  = 1d0
      do ix =1,i
        igan=igan*ix
      enddo
      end











