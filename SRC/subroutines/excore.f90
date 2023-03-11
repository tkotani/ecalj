module m_excore
  public excore
  contains
subroutine excore(nrmx,nl,nnc,nclass,nspn,nbas, &
     phic,nindxc,iclass, &
     a,b,nr,rofi)
  !- Calculate core-core exchange energy.
  use m_lldata,only: ll
  implicit none
  integer(4):: nrmx,nnc,nclass,nspn,ncmx,nl,nbas
  integer(4):: nindxc(0:nl-1,nclass),nr(nclass),iclass(nclass)
  real(8) :: a(nclass),b(nclass),rofi(nrmx,nclass)
  real(8) :: phic (nrmx,0:nl-1,nnc,nclass,nspn),wgtx

  real(8):: exacct,exacc(nbas,nspn),rydberg
  integer,parameter :: nsxamx=10000000
  real(8),allocatable :: &
       gc(:,:,:,:), sxadata(:),sigkcc(:,:,:,:,:)
  integer(4),allocatable:: indxsxa(:,:),lcr(:,:,:),ncr(:,:)
  integer(4):: ir,l,n,ic,ncrr,ncmxx(nclass),kmax,k,it,isp, &
       isx,lm1,lm2,lm3,lm4,ibas,icr,icrd,nsxatot,ifexcore
  real(8):: rkp(nrmx,0:2*(nl-1)),rkm(nrmx,0:2*(nl-1))
  do ic = 1,nclass
     ncmxx(ic) = sum(nindxc(0:nl-1,ic))
  enddo
  ncmx = maxval(ncmxx(1:nclass))

  print *,' ncmx nl nspn=',ncmx,nl,nspn

  ! --- convert format of core function ---
  ! nindxc(l,ic) phic ---> ncr lcr gc

  allocate( gc(nrmx,ncmx,nspn,nclass) &
       ,lcr(ncmx,nspn,nclass),ncr(nspn,nclass) )
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
  !---------------------------------------------------------------
  print *,' goto alloc'
  allocate( sxadata(nsxamx),indxsxa(6,nsxamx), &
       sigkcc(0:2*(nl-1), ncmx, ncmx, nclass,nspn) )
  print *,' end of alloc'

  kmax  = 2*(nl-1)
  do ic = 1,nclass
     print *,' make rkp rkm ic=',ic
     !- rkp,rkm ----------- This is from subroutine bess
     do k=0, kmax
        rkp(1,k)=0d0
        rkm(1,k)=0d0
        do ir=2,nr(ic)
           rkp(ir,k) = rofi(ir,ic)**k
           rkm(ir,k) = rofi(ir,ic)**(-k-1)
        enddo
     enddo
     print *,' end of make rkp rkm ic=',ic
     !- radial integrals
     do isp = 1,nspn
        call intsigkcc(sigkcc(0,1,1,ic,isp), ncmx, &
             gc(1,1,isp,ic), &
             a(ic),b(ic),nr(ic),rofi(1,ic),nl, &
             ncr(isp,ic),lcr(1,isp,ic), &
             kmax, rkp, rkm )
     enddo
  enddo

  !- make structure spherical part --- sxadata
  call mksxa(nl,nsxamx, &
       sxadata,indxsxa,nsxatot)

  !- core-core exchange =sxadata * radial integral
  print * !; print *,' go into EXEX.CORE part'
  exacc = 0d0
  do 301 isx=1,nsxatot
     !	print *,' do 300 isx isp=',isx, isp
     do 300 isp=1,nspn
        k   = indxsxa (1,isx)
        lm1 = indxsxa (2,isx)
        lm3 = indxsxa (3,isx)
        lm2 = indxsxa (4,isx)
        lm4 = indxsxa (5,isx)
        !	print *,' do 300 isx isp= xxx'
        if(lm3 /= lm2 .OR. lm1 /= lm4) cycle
        !- icr->l3l2 icrd->l1l4
        do 350 ibas = 1,nbas
           ic = iclass(ibas)
           !	  print *,' ibas ic=',nbas,ibas,ic
           do icr = 1, ncr(isp,ic)
              do icrd= 1, ncr(isp,ic)
                 ! ccccccccccccccccccccc
                 !	print *,' isx=',isx,ic,isp,icr,icrd,lm1,lm3
                 !	print *,' lcr=',lcr(icr, isp,ic),lcr(icrd, isp,ic)
                 !	print *,' ll=',ll(lm1),ll(lm3)
                 ! ccccccccccccccccccccc

                 ! ccccccccccccccccccccccccccccccccccc
                 !     test
                 !            if(ll(lm3)==1.or.ll(lm1)==1) cycle
                 ! ccccccccccccccccccccccccccccccccccc


                 if( lcr(icr, isp,ic) /= ll(lm3)) cycle
                 if( lcr(icrd,isp,ic) /= ll(lm1)) cycle
                 exacc(ibas,isp) = exacc(ibas,isp) &
                      - sxadata(isx) * sigkcc(k,icr,icrd,ic,isp)
              enddo
           enddo
350     enddo
        !---------------------
300  enddo
301 enddo
  wgtx=1d0
  if(nspn==1) wgtx=2d0
  exacct = sum(exacc(1:nbas,1:nspn))*wgtx

  open (newunit=ifexcore,file='TEEXXcc')
  write(6,*) '==== EXCORE ==> TEEXXcc ============'
  write(ifexcore,*) '======================================='
  write(ifexcore,*) '  Exchange energy core-core   Exx (eV) '
  write(ifexcore,*) '======================================='
  write(ifexcore,*) ' *** '
  write(ifexcore, &
       "( f20.10,2i4,' ! Ex core-core total (eV) --  nbas nsp ')") &
       exacct*rydberg(), nbas, nspn
  write(6,"(' excore total =',f13.6)") exacct*rydberg()
  do isp =1,nspn
     do ibas=1,nbas
        write(6,"(' ibas isp=',2i3,' Ex core-core (eV) =',f13.6)") &
             ibas,isp,exacc(ibas,isp)*rydberg()
        !	write(ifexcore,"(' ibas isp=',2i3,' Ex core-core (Ry) =',f13.6)")
        !     &  ibas,isp,exacc(ibas,isp)
        write(ifexcore,"( f20.10,2i4,' ! Ex core-core (eV)  ibas isp ')") &
             exacc(ibas,isp)*rydberg(),ibas,isp
     enddo
  enddo
end subroutine excore

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine intsigkcc(sigkcc,ncmx, &
     gc, a,b,nr,rofi, nl, &
     ncr,lcr , &
     kmax, rkp, rkm )
  !- integral \sigma^k(l1,l3,ip1,ip3,l2,l4,ip2,ip4) in the Phys.Rev.B34, 5512 ,Simpson
  !---------------------------------------------------------------------
  !i ncmx,nrec, gl,gpl,gc , a,b,nr,rofi,nl,ncr,lcr
  !o  sigkcc
  !---------------------------------------------------------------------
  implicit none
  integer :: nr, ncmx, nl,l1,l2,l3,l4,l1l3,l2l4,k &
       ,ip1,ip3,ip2,ip4,ip2ip4,icr,l1l4,icrd,ip1ip3,ir,lr0 &
       ,ncr,lcr(ncmx)
  double precision :: a,b,rofi(nr) &
       , gc(nr,ncmx) &
       , sum1,sum2 &
       , sigkcc(0:2*(nl-1), ncmx, ncmx  ) &
       , int1(nr),int2(nr),a1(nr),a2(nr),b1(nr),f13(nr)
  integer :: kmax
  double precision :: rkp(nr,0:kmax),rkm(nr,0:kmax)
  print *,' intsigkcc:'
  !--------------
  print *,' normcheck:'
  do icr  = 1,ncr; l2 =lcr(icr);  l3=l2
     do icrd = 1,ncr; l1 =lcr(icrd); l4=l1
        call gintxx(gc(1,icr),gc(1,icrd),a,b,nr,sum1)
        write(6,"(' norm check ='2i3,2x,f13.6)") icr,icrd,sum1
     enddo
  enddo
  !      do ir=1,nr
  !        write(6,"(' rofi ='i3,2x,f13.6)") ir,rofi(ir)
  !      enddo

  ! ---------------between core and core->  sigkc
  do 2216 icr  = 1,ncr; l2 =lcr(icr);  l3=l2
     do 2215 icrd = 1,ncr; l1 =lcr(icrd); l4=l1
        do 2210 k  = abs(l1-l3),l1+l3
           if( k<abs(l2-l4) .OR. k>l2+l4 ) cycle
           if(mod(k+l1+l3,2)==1 .OR. mod(k+l2+l4,2)==1) cycle
           a1(1)    = 0d0
           a1(2:nr) = rkp(2:nr,k) *gc(2:nr,icrd)
           a2(1)    = 0d0
           a2(2:nr) = rkm(2:nr,k) *gc(2:nr,icrd)
           b1(1:nr) = gc(1:nr,icr)
           call intn_smp_g(a1,b1,int1,a,b,rofi,nr,lr0)
           call intn_smp_g(a2,b1,int2,a,b,rofi,nr,lr0)

           ! check write ---------------------------------------------
           call gintxx(a1,b1,a,b,nr,sum1)
           call gintxx(a2,b1,a,b,nr,sum2)
           write(6,"(' integral ='2d13.6)") sum1,sum2
           write(6,"(' integral ='2d13.6)") int1(1),int2(1)
           !        endif
           !---
           f13(1)    = 0d0
           f13(2:nr) = rkm(2:nr,k) *( int1(1)-int1(2:nr) ) &
                + rkp(2:nr,k) * int2(2:nr)
           a1(1:nr)  = gc(1:nr,icr) *f13(1:nr)
           b1(1:nr)  = gc(1:nr, icrd)
           call gintxx(a1,b1,A,B,NR, sigkcc(k, icr, icrd ) )
           ! chekc write
           !        if(iprint().ge.130)
           write(6,"( ' k icr icrd =',3i5,'  sigkcc=',d15.8)") &
                k, icr,icrd ,sigkcc(k,icr,icrd)
2210    enddo
2215 enddo
2216 enddo
end subroutine intsigkcc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mksxa(nl,nsxamx, &
     sxadata,indxsxa,nsxatot)
  !- make sx (spherical-dependent part for exx energy.) ------------c
  !i INPUT
  !i   nl
  !i   nsxamx
  !o   indxsxa, sxdata, nsxatot
  !i
  !r       ngaut is number for non-zero cgau.
  !r  Generation of
  !r     cgau(k,L',L) for R=R'
  !-------------------------------------------------------------------------c
  implicit none
  integer :: nl,nlx, &
       m,istcl,lrd,lr,j1,j2,j3,m1,m2,m3,jm1,jm2,jm3, &
       ngau,ngaut,ngau2,k1,k2,j4,m4,jm1m, &
       jm2m,jm4,km1,km2, &
       i,j,k ,jmax
  integer :: nsxamx, nsxatot
  integer :: indxsxa(6,nsxamx)
  double precision :: sxadata(nsxamx)

  double complex &
       sxx(-(nl-1):(nl-1),-(nl-1):(nl-1),-(nl-1):(nl-1),-(nl-1):(nl-1))
  !     & ,syy(-3:3,-3:3,-3:3,-3:3)
  double precision :: osq2

  ! takao obtain cg upto 2*l
  !      double precision cg( (2*nl-1)**2,(2*nl-1)**2, 4*(nl-1) )
  double precision :: cg( (2*nl-1)**2,(2*nl-1)**2, 0:4*(nl-1) ) &
       , cgau(0:2*(nl-1),nl**2, nl**2 ) ,dum
  !      integer iprint

  !      external clebsh ,mscmul
  double complex msc(0:1,2,2),mcs(0:1,2,2),Img &
       ,ap,am!,mscmul

  integer :: ngautx
  data Img/(0.0d0,1.0d0)/

  !----------------------------------------------------------------------
  print *,' goto mksxa'

  !  msc, conversion matrix generation 1->m 2->-m for m>0
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

  !- CG coeff. generation
  jmax=  2*(nl-1)
  call clebsh_t( cg,jmax )
  print *,' end of clebsh'
  !- check write
  !      if(iprint().ge.120) then
  !        do 106 j1=0,nl-1
  !        do 106 j2=0,nl-1
  !        do 106 j3=abs(j1-j2),j1+j2
  !        do 106 m1=-j1,j1
  !        do 106 m2=-j2,j2
  !          jm1=j1**2+j1+1+m1
  !          jm2=j2**2+j2+1+m2
  !          m3=m1+m2
  !          if(abs(m3).gt.j3) go to 106
  !          write(6,105) j1,m1,j2,m2,j3,m3,cg(jm1,jm2,j3)
  !  105     format( ' ###  j1m1 j2m2 jm = ',2i3,' ',2i3,' ',2i3, d12.5)
  !  106   continue
  !      endif

  !---------------------------------------------------------------
  !     make Gaunt coefficient. Def. by eq(A2), and Rose eq.4.34
  !     cgau(k,L',L)  ; k=j2 L'=jm3 L=jm1
  !     cgau(k,j3m3,j1m1)=(-1)**(m1-m3)*cgau(k,j1m1,j3m3)

  !  Gaunt(or CG) coefficients are not zero when the conditions 1,2,
  !     and 3 are hold.
  !  1.   (-1)**(j1+j2+j3).eq.1.
  !  2.   abs(m1-m3).le.j2
  !  3.   j2.ge.abs(j3-j1).and.j2.le.(j3+j1)

  ngaut=0
  ngautx=0
  do  j3=0,nl-1
     do  j1=0,nl-1
        ! takao  max of j2 is 2*lmx, because max of j3 and j1 is lmx. cgau(j2,jm3,jm1)
        do  j2=0,2*(nl-1)
           do  m1=-j1,j1
              do  m3=-j3,j3
                 jm1= j1**2+m1+j1+1
                 jm3= j3**2+m3+j3+1
                 m2 = m3 - m1
                 jm2= j2**2+m2+j2+1

                 cgau(j2,jm3,jm1)=0.0d0
                 if(abs(m2) <= j2 .AND. &
                      mod(j1+j2+j3,2) == 0 .AND. &
                      j2 >= abs(j3-j1) .AND. j2 <= (j3+j1) ) then
                    ngaut=ngaut+1
                    ! cc
                    if( j2 <= nl-1) ngautx=ngautx+1
                    ! c
                    cgau(j2,jm3,jm1)=cg(jm1,jm2,j3) &
                         *cg(j1**2+j1+1,j2**2+j2+1,j3) &
                         *sqrt( (2.0d0*j1+1.0d0)/(2.0d0*j3+1.0d0) )
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  PRINT *,' * Gaunt coef. end;  num of Gaunt; nl  =' &
       , ngaut ,nl
  print *,' ngautx=',ngautx

  !$$$c----check write--------------------------------------c
  !$$$      if(iprint().ge.120) then
  !$$$        ngau=0
  !$$$        ngau2=0
  !$$$        do 31 j3=0,nl-1
  !$$$        do 31 j1=0,nl-1
  !$$$        do 31 j2=0,2*(nl-1)
  !$$$        do 31 m1=-j1,j1
  !$$$        do 31 m3=-j3,j3
  !$$$          jm1= j1**2+m1+j1+1
  !$$$          jm3= j3**2+m3+j3+1
  !$$$          m2 = m3 - m1
  !$$$          jm2= j2**2+m2+j2+1
  !$$$          if(jm3.ge.jm1) then
  !$$$          if( (-1)**(j1+j2+j3).eq.1.and.abs(m1-m3).le.j2.and.
  !$$$     &      j2.ge.abs(j3-j1).and.j2.le.(j3+j1) ) then
  !$$$
  !$$$            write(6,119)j2,j3,m3,j1,m1, cgau(j2,jm3,jm1),
  !$$$     &        cgau(j2,jm1,jm3)*(-1)**(m1-m3),(-1)**(j1+j2+j3)
  !$$$            ngau=ngau+1
  !$$$          else
  !$$$            write(6,129)j2,j3,m3,j1,m1, cgau(j2,jm3,jm1),
  !$$$     &        cgau(j2,jm1,jm3)*(-1)**(m1-m3),(-1)**(j1+j2+j3)
  !$$$            ngau2=ngau2+1
  !$$$          endif
  !$$$          endif
  !$$$  119     format('   gaunt j2 j3m3 j1m1 parity= '
  !$$$     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
  !$$$  129     format('                   gaunt j2 j3m3 j1m1 parity= '
  !$$$     &              ,i3,' ',2i3,' ',2i3,' ',2d23.16,'  ',i3)
  !$$$   31   continue
  !$$$      endif

  !-- sxa cal. for R=R', see eq.(A11) --------------------------------c
  print *
  print *
  print *, '  *** Go into SXA cal. for R=Rdash pair '
  nsxatot=0
  !      nsxatot2=0
  do 147 k=0,2*(nl-1)
     do 137 j1=0,nl-1
        do 127 j2=0,nl-1
           do 117 j3=0,nl-1
              do 107 j4=0,nl-1
                 do 237 m1=-j1,j1
                    do 227 m2=-j2,j2
                       do 217 m3=-j3,j3
                          do 207 m4=-j4,j4
                             jm1 =j1**2+j1+1+m1
                             jm1m=j1**2+j1+1-m1

                             jm2 =j2**2+j2+1+m2
                             jm2m=j2**2+j2+1-m2

                             jm3 =j3**2+j3+1+m3
                             jm4 =j4**2+j4+1+m4

                             sxx(m1,m2,m3,m4)=0.0d0
                             if(m1+m2+m3+m4 == 0) then
                                sxx(m1,m2,m3,m4)=2*(-1)**(m1+m2)*cgau(k,jm3,jm1m) &
                                     *cgau(k,jm2m,jm4)
                             endif
                             ! cccccccccccccccccccccccc
                             !            syy(m1,m2,m3,m4)= sxx(m1,m2,m3,m4)
                             ! cccccccccccccccccccccccc
207                       enddo
217                    enddo
227                 enddo
237              enddo

                 ! convert to real harmonics rep.-------------
                 call convsx(sxx,j1,j2,j3,j4,nl,msc,mcs)

                 do 437 m1=-j1,j1
                    do 427 m2=-j2,j2
                       do 417 m3=-j3,j3
                          do 407 m4=-j4,j4
                             jm1 =j1**2+j1+1+m1
                             jm2 =j2**2+j2+1+m2
                             jm3 =j3**2+j3+1+m3
                             jm4 =j4**2+j4+1+m4

                             !          if(iprint().ge.130.and.
                             !     &       DREAL( sxx(m1,m2,m3,m4) ).gt.1.0d-6 ) write(6,1007)
                             !     &      k,j1,m1,j2,m2,j3,m3,j4,m4,syy(m1,m2,m3,m4),
                             !     &      sxx(m1,m2,m3,m4)

                             dum= DREAL( sxx(m1,m2,m3,m4) )

                             if(dabs(dimag(sxx(m1,m2,m3,m4))) >1d-8) &
                                ! top2rx 2013.08.09 kino     &      stop ' MAKSX; im part of sxx .ne. 0'
                                  call rx( ' MAKSX; im part of sxx /= 0')

                             if(abs(dum)> 1d-8) then
                                if(nsxatot >= nsxamx) &
                                ! top2rx 2013.08.09 kino     &        stop ' MAKSX: enlarge the size of nsxamx '
                                     call rx( ' MAKSX: enlarge the size of nsxamx ')
                                nsxatot=nsxatot+1
                                ! note its oeder! jm1.ge.jm3
                                indxsxa(1,nsxatot)=k
                                indxsxa(2,nsxatot)=jm1
                                indxsxa(3,nsxatot)=jm3
                                indxsxa(4,nsxatot)=jm2
                                indxsxa(5,nsxatot)=jm4
                                sxadata(nsxatot  )=dum
                             endif

407                       enddo
417                    enddo
427                 enddo
437              enddo
107           enddo
117        enddo
127     enddo
137  enddo
147 enddo

  sxadata(1:nsxatot)=sxadata(1:nsxatot)/2d0
end subroutine mksxa

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      subroutine clebsh()
subroutine clebsh_t(cg,j1mx)
  !- generate crebsh gordon coefficient
  !  takao 1993 8/17
  !     eq.3.18 Rose's book , 'elementary theory of the angular momentum

  !    C(j1,j2,j; m1,m2,m1+m2) = cg(jm1,jm2,j)
  !      jm1=j1**2+(j1+1+m1)

  !    C(j1,j2,j;0,0,0) at j1+j2+j=even is exactly 0  (eq.3.22),
  !      however this routine gives the value oder 1.0d-16, according to
  !      numerical cancellation is not complete.
  !----------------------------------------------------------------c
  !      implicit double precision (a-h,o-z)
  use m_lldata,only: ll
  implicit none
  integer :: j1mx
  integer :: &
       j1,j2,jm1,jm2,m1,m2,m3,j3,nu, &
       id1,id2,id3
  double precision :: &
       k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
  double precision :: &
       fac,fac2, &!@igan, &
       cg( (j1mx+1)**2, (j1mx+1)**2, 0:2*j1mx)
  print *, ' go into clebsh j1mx=',j1mx
  do 403   jm1=1, (j1mx+1)**2
     do 405 jm2=1, (j1mx+1)**2
        j1 = ll(jm1)
        m1 = jm1-(j1**2+j1+1)
        j2 = ll(jm2)
        m2 = jm2-(j2**2+j2+1)
        m3  =  m1+m2

        do 303 j3=0, 2*j1mx

           !            write(6,309) j1,m1,jm1,j2,m2,jm2,j3,m3
           !  309       format(' j1 m1 j1m1=',3i4,'  j2 m2 jm2',3i4,'  j3 m3=',2i4)
           cg(jm1,jm2,j3)=0.0d0
           ! cc
           if( j3 > j1+j2 .OR. j3 < abs(j1-j2) )  go to 303
           if(abs(m3) > j3) goto 303
           !              write(6,*) '  goto calculation'
           ! cc
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
           ! top2rx 2013.08.09 kino            if(fac2.eq.0.0d0) stop  ' k6k7k8k9k10=fac2=0'
           if(fac2 == 0.0d0) call rx( ' k6k7k8k9k10=fac2=0')
           fac  = sqrt( (2*j3+1) *k1*k2*k3*k4*k5 /fac2 )

           do 36 nu =0, j2+j3+m1
              id1=j3-j1+j2-nu
              id2=j3+m3-nu
              id3=nu+j1-j2-m3
              if(id1 >= 0 .AND. id2 >= 0 .AND. id3 >= 0 .AND. &
                   j2+j3+m1-nu >= 0 .AND. j1-m1+nu >= 0  )  then

                 k1=igan(j2+j3+m1-nu)
                 k2=igan(j1-m1+nu)
                 k3=igan(nu)
                 k4=igan(id1)
                 k5=igan(id2)
                 k6=igan(id3)
                 fac2= k3*k4*k5*k6
                 cg(jm1,jm2,j3)=cg(jm1,jm2,j3) &
                      + fac*(-1)**(nu+j2+m2) &
                      *k1*k2 / fac2
              endif
36         enddo
           ! cccccccccccccccccccccccccccccccccccccccccccccc
           !--test write
           !      write(105,105) j1,m1,j2,m2,j3,m3,cg(jm1,jm2,j3)
           !  105 format( ' ###  j1m1 j2m2 jm = ',2i3,' ',2i3,' ',2i3, d23.16)
           ! ccccccccccccccccccccccccccccccccccccccccccccccc
           !x            endif
303     enddo
405  enddo
403 enddo
end subroutine clebsh_t

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine convsx(sxx,j1,j2,j3,j4,nl,msc,mcs)
  !- sx based on spherical har. rep. is converted to cub rep.
  !i input
  !i     sxx; in the harmonic rep.
  !o output
  !o     sxx; in the sperical rep.

  !r inversion test was included; uncomment the lines commented by ct.
  !r
  !------------------------------------------------
  implicit none
  !      implicit double precision(a-h,o-z)
  !      parameter(nlx=3)
  integer :: &
       nl,j1,j2,j3,j4,m1,m2,m3,m4
  double complex &
       sxx (-(nl-1):(nl-1),-(nl-1):(nl-1), &
       -(nl-1):(nl-1),-(nl-1):(nl-1))
  !     &   sxxx(-(nlx-1):(nlx-1),-(nlx-1):(nlx-1),
  !     &        -(nlx-1):(nlx-1),-(nlx-1):(nlx-1)),
  !     &   syyy(-(nlx-1):(nlx-1),-(nlx-1):(nlx-1),
  !     &        -(nlx-1):(nlx-1),-(nlx-1):(nlx-1))

  double complex &
       msc(0:1,2,2),Img,ap,am!,mscmul
  !  inversion test
  double complex &
       mcs(0:1,2,2)

  !      if(nlx.ne.nl) stop 'CONVSX: nlx.ne.nl'
  do m1=   1, j1
     do m2= -j2, j2
        do m3= -j3, j3
           do m4= -j4, j4
              ap=sxx( m1,m2,m3,m4)
              am=sxx(-m1,m2,m3,m4)
              sxx(  m1,m2,m3,m4) = mscmul(1,m1,msc,ap,am)
              sxx( -m1,m2,m3,m4) = mscmul(2,m1,msc,ap,am)
           enddo
        enddo
     enddo
  enddo
  do m1= -j1, j1
     do m2=   1, j2
        do m3= -j3, j3
           do m4= -j4, j4
              ap=sxx( m1, m2,m3,m4)
              am=sxx( m1,-m2,m3,m4)
              sxx(  m1, m2,m3,m4) = mscmul(1,m2,msc,ap,am)
              sxx(  m1,-m2,m3,m4) = mscmul(2,m2,msc,ap,am)
           enddo
        enddo
     enddo
  enddo
  do  m1= -j1, j1
     do  m2= -j2, j2
        do  m3=   1, j3
           do  m4= -j4, j4
              ap=sxx( m1, m2, m3,m4)
              am=sxx( m1, m2,-m3,m4)
              sxx(  m1, m2, m3,m4) = mscmul(1,m3,msc,ap,am)
              sxx(  m1, m2,-m3,m4) = mscmul(2,m3,msc,ap,am)
           enddo
        enddo
     enddo
  enddo
  do  m1= -j1, j1
     do  m2= -j2, j2
        do  m3= -j3, j3
           do  m4=   1, j4
              ap=sxx( m1, m2, m3, m4)
              am=sxx( m1, m2, m3,-m4)
              sxx(  m1, m2, m3, m4) = mscmul(1,m4,msc,ap,am)
              sxx(  m1, m2, m3,-m4) = mscmul(2,m4,msc,ap,am)
           enddo
        enddo
     enddo
  enddo
end subroutine convsx
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function mscmul(i,m,msc,ap,am)
  !----multiple mat matrix
  integer :: mx,m,i
  double complex msc(0:1,2,2),ap,am
  !  msc, conversion matrix generation 1->m 2->-m for m>0
  mx=mod(m,2)
  mscmul=ap*msc(mx,1,i)+am*msc(mx,2,i)
END function mscmul

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine intn_smp_g(g1,g2,int,a,b,rofi,nr,lr0)
  !-- intergral of two wave function.
  !   This is for true g1(ir) = r phi(ir), where phi is the true wave function.

  ! int(r) = \int_(r)^(rmax) u1(r') u2(r') dr'

  ! lr0 dummy index, now not used.
  ! simpson rule ,and with higher rule for odd devision.
  ! --------------------------------------------------------------
  !      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  IMPLICIT none
  integer :: nr,ir,lr0
  double precision :: g1(nr),g2(nr),int(nr),a,b,rofi(nr),w1,w2,w3 &
       ,ooth,foth
  integer(4):: mx=1
  ! Local parameters
  ! one over three-> oot,  and so on.
  data ooth,foth/0.33333333333333333,1.3333333333333333333/
  data w1,w2,w3/0.41666666666666666,0.6666666666666666666, &
       -0.083333333333333333/
  !                        xxxxx
  if(mod(nr,2) == 0) &
       ! top2rx 2013.08.09 kino     &  stop ' intn_smp_g: nr should be odd. '
       call rx( ' intn_smp_g: nr should be odd. ')

  int(1)=0.0d0
  ! l00 means u1(r)u2(r)~r**lr0 near r~0
  !c simplest formula
  !c      do 10 ir=3,nr
  !c        int(ir)=int(ir-1)
  !c     &         +0.5d0*G1(ir-1)*G2(ir-1)*( a*(b+rofi(ir-1)) )**mx
  !c     &         +0.5d0*G1(ir)*G2(ir)*( a*(b+rofi(ir)) )**mx
  !c   10 continue
  ! simpson rule
  DO  10  IR = 3,NR,2
     int(ir)=int(ir-2) &
          + ooth*G1(IR-2)*G2(IR-2)*( a*(b+rofi(ir-2)) )**mx &
          + foth*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )**mx &
          + ooth*G1(IR)*G2(IR)*( a*(b+rofi(ir)) )**mx
10 enddo

  ! At the value for odd points, use the same interpolation above
  do 20 ir = 2,nr-1,2
     int(ir)=int(ir-1) &
          + w1*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )**mx &
          + w2*G1(IR)  *G2(IR)*  ( a*(b+rofi(ir)  ) )**mx &
          + w3*G1(IR+1)*G2(IR+1)*( a*(b+rofi(ir+1)) )**mx
20 enddo

  do 30 ir=1,nr
     int(ir)=int(nr)-int(ir)
30 enddo
end subroutine intn_smp_g

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real(8) function igan(i)
  integer:: i,ix
  igan= product([(ix,ix=1,i)])
END function igan

end module m_excore
