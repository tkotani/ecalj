c-------------------------------------------------------------
      subroutine sigint(rkp,rkm,kmx,a,b,nr,phi1,phi2,k,rofi, sig)
      implicit none
      integer(4) :: nr,kmx,k,ir
      real(8):: a,b, a1(nr),a2(nr),b1(nr),rkp(nr,0:kmx),rkm(nr,0:kmx),
     &    int1x(nr),int2x(nr), phi1(nr), phi2(nr),rofi(nr),sig
      real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
c
      a1(1) = 0d0;  a1(2:nr) = rkp(2:nr,k)
      a2(1) = 0d0;  a2(2:nr) = rkm(2:nr,k)
      b1(1:nr) = phi1(1:nr)
      call intn_smpxxx(a1,b1,int1x,a,b,rofi,nr,0)
      call intn_smpxxx(a2,b1,int2x,a,b,rofi,nr,0)
c
      a1(1) = 0d0; a1(2:nr) = 
     &  rkm(2:nr,k) *( int1x(1)-int1x(2:nr) )+ rkp(2:nr,k) * int2x(2:nr)
      b1(1:nr) = phi2(1:nr)
      call gintxx(a1,b1,A,B,NR, sig )
      sig = sig/(2*k+1)*fpi
      end



      subroutine mkjp3( q,ngc,ngvecc, alat, qlat, lxx,lx,nxx,nx,  
     &                   bas, a,b,rmax,nr,nrx,rprodx,
     o        rojp,sgpb,fouvb)
c     i    ibas,
c     o        sgpp,fouvp)
C- Make integrals in each MT. and the Fourier matrix.
Cr the integrals rojp, fouvb,fouvp
Cr are for  J_L(r)= j_l(sqrt(e) r)/sqrt(e)**l Y_L,
Cr which behaves as r^l/(2l+1)!! near r=0.
Cr
Cr oniste integral is based on
Cr 1/|r-r'| = \sum 4 pi /(2k+1) \frac{r_<^k }{ r_>^{k+1} } Y_L(r) Y_L(r')
Cr See PRB34 5512(1986) for sigma type integral
Cr
      implicit none
      integer(4) :: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx
      real(8)    :: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat,
     &              qlat(3,3)
ci rho-type onsite integral
      complex(8) :: rojp(ngc, (lxx+1)**2)
ci sigma-type onsite integral
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2)
c     &              sgpp(ngc,  ngc,  (lxx+1)**2)
      real(8),allocatable::cy(:),yl(:)
ci Fourier
      complex(8) :: 
     &              fouvb(ngc,  nxx, (lxx+1)**2)
c     &              fouvp(ngc,  ngc, (lxx+1)**2)
c internal
      integer(4) :: nlx,ig1,ig2,l,n,ir,n1,n2,lm !, ibas
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      real(8)    :: pi,fpi,tpiba, qg1(3),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1,absqg2,
     &  fac,radint,radsigo(0:lx),radsig(0:lx),phi(0:lx),psi(0:lx)
     &   ,rofi(nr),r2s,sig,sig1,sig2,sigx(0:lx),sig0(0:lx) ,qg2(3)
      complex(8) :: img =(0d0,1d0),phase
      complex(8),allocatable :: pjyl(:,:)
      real(8),allocatable ::ajr(:,:,:),a1(:,:,:),rkpr(:,:),rkmr(:,:),
     &  qg(:,:),absqg(:)
c-------------------------------------------------
      print *,' mkjp2:'
      nlx = (lx+1)**2
      allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),
     &  rkpr(nr,0:lx),rkmr(nr,0:lx), qg(3,ngc),absqg(ngc),
     &  pjyl((lx+1)**2,ngc) )
c
      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
      tpiba = 2*pi/alat
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
c      print *,' mkjp: end of sylmnc'
c q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
      do ig1 = 1,ngc
        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        qg1(1:3) = qg(1:3,ig1)
        absqg(ig1)  = sqrt(sum(qg1(1:3)**2))
        absqg1   = absqg(ig1)
        phase = exp( img*sum(qg1(1:3)*bas(1:3))*alat  )
        call sylm(qg1/absqg1,yl,lx,r2s) !spherical factor Y( q+G )
        do lm =1,nlx
          l = ll(lm)
          pjyl(lm,ig1) = fpi*img**l *cy(lm)*yl(lm) *phase  *absqg1**l
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
      enddo

c rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
      rofi(1) = 0d0
      do ir   = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      do l = 0,lx
        rkpr(1:nr,l) = rofi(1:nr)**(l      +1 )
        rkmr(2:nr,l) = rofi(2:nr)**(-l-1   +1 )
        rkmr(1,l)    = rkmr(2,l)
      enddo

c rojp
      print *,' mkjp: rojp'
      do ig1 = 1,ngc
        call wronkj( absqg(ig1)**2, 0d0, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
        do lm = 1,nlx
          l = ll(lm)
          rojp(ig1,lm) = (-fjj(l))* pjyl(lm,ig1)
        enddo
      enddo

c ajr
      do ig1 = 1,ngc
        do ir  = 1,nr
          call bessl(absqg(ig1)**2*rofi(ir)**2,lx,phi,psi)
          do l   = 0,lx
            ajr(ir,l,ig1) = phi(l)* rofi(ir) **(l +1 )
            ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
          enddo
ccccccccccccccccccccccccccc
c        write(116,'(i3,10d13.6)') ir, rofi(ir), ajr(ir,0:lx,ig1)
ccccccccccccccccccccccccccc
        enddo
ccccccccccccccccccccccccccc
c        write(6,*) ig1,sum(ajr(1:nr,0:lx,ig1))
ccccccccccccccccccccccccccc
      enddo

c-------------------------
      do ig1 = 1,ngc
        call sigintAn1( absqg(ig1), lx, rofi, nr
     o                  ,a1(1:nr, 0:lx,ig1) )
      enddo
c-------------------------
      print *,' mkjp: sgpb'
c sgpb
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)                      ! r jl        , r B(r)
c        call sigint(rkpr,rkmr, lx,a,b,nr,  ajr(1:nr,l,ig1),rprodx(1,n,l)
c     &            ,l, rofi, radint)
            call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
            sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
ccccccccccccccccccc
          enddo
        enddo
      enddo

c---------------------------------------
c sgpp block------->removed
c---------------------------------------

c Fourier
c fouvb
      print *,' mkjp: Four'
      fouvb=0d0
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)
cccccccccccccccccccccccccccccccccccccccccccccc
c         print *,' ig1 lm l n=',ig1,lm,l,n
cccccccccccccccccccccccccccccccccccccccccccccc
            call gintxx(ajr(1,l,ig1), rprodx(1,n,l), a,b,nr,
     o                 radint )
cccccccccccccccccccccccccccccccccccccccccccccc
c         print *,' radint=',radint
cccccccccccccccccccccccccccccccccccccccccccccc
            fouvb(ig1, n, lm) =
     &    fpi/(absqg(ig1)**2) *dconjg(pjyl(lm,ig1))*radint

          enddo
        enddo
      enddo
ccccccccccccccccccccc
c        write(6,*)' fourvb sum=',sum (fouvb)
ccccccccccccccccccccc

c-----------------------------
c fouvp block --->removed
c-----------------------------

      deallocate(ajr,a1,   rkpr,rkmr, qg,absqg,   pjyl)
      if (allocated( cy )) deallocate(cy)
      if (allocated( yl )) deallocate(yl)
      end





c--------------------------------------------------------------------------
      subroutine mkjp2( q,ngc,ngvecc, alat, qlat, lxx,lx,nxx,nx,  
     &                   bas, a,b,rmax,nr,nrx,rprodx,
     o        rojp,sgpb,fouvb,
c     i    ibas,  
     o        sgpp,fouvp)
C- Make integrals in each MT. and the Fourier matrix.
Cr the integrals rojp, fouvb,fouvp
Cr are for  J_L(r)= j_l(sqrt(e) r)/sqrt(e)**l Y_L,
Cr which behaves as r^l/(2l+1)!! near r=0.
Cr
Cr oniste integral is based on
Cr 1/|r-r'| = \sum 4 pi /(2k+1) \frac{r_<^k }{ r_>^{k+1} } Y_L(r) Y_L(r')
Cr See PRB34 5512(1986) for sigma type integral
Cr
      implicit none
      integer(4) :: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx
      real(8)    :: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat,
     &              qlat(3,3)
ci rho-type onsite integral
      complex(8) :: rojp(ngc, (lxx+1)**2)
ci sigma-type onsite integral
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2),  
     &              sgpp(ngc,  ngc,  (lxx+1)**2)
      real(8),allocatable::cy(:),yl(:)
ci Fourier
      complex(8) :: 
     &              fouvb(ngc,  nxx, (lxx+1)**2),
     &              fouvp(ngc,  ngc, (lxx+1)**2)
c internal
      integer(4) :: nlx,ig1,ig2,l,n,ir,n1,n2,lm !, ibas
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      real(8)    :: pi,fpi,tpiba, qg1(3),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1,absqg2,
     &  fac,radint,radsigo(0:lx),radsig(0:lx),phi(0:lx),psi(0:lx)
     &   ,rofi(nr),r2s,sig,sig1,sig2,sigx(0:lx),sig0(0:lx) ,qg2(3)
      complex(8) :: img =(0d0,1d0),phase
      complex(8),allocatable :: pjyl(:,:)
      real(8),allocatable ::ajr(:,:,:),a1(:,:,:),rkpr(:,:),rkmr(:,:),
     &  qg(:,:),absqg(:)
c      logical  :: newsg=.true. !new version (Mar2001) for sgpp sgpb or not
c-------------------------------------------------
      print *,' mkjp2:'
      nlx = (lx+1)**2
      allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),
     &  rkpr(nr,0:lx),rkmr(nr,0:lx), qg(3,ngc),absqg(ngc),
     &  pjyl((lx+1)**2,ngc) )
c
      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
      tpiba = 2*pi/alat
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
c      print *,' mkjp: end of sylmnc'
c q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
      do ig1 = 1,ngc
        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        qg1(1:3) = qg(1:3,ig1)
        absqg(ig1)  = sqrt(sum(qg1(1:3)**2))
        absqg1   = absqg(ig1)
cccccccccccccccccccccccccccccccccccccccccccccc
c        write(6,'(i4,3f12.6,3x,f12.6)')ig1, qg1,absqg1
cccccccccccccccccccccccccccccccccccccccccccccc
        phase = exp( img*sum(qg1(1:3)*bas(1:3))*alat  )
        call sylm(qg1/absqg1,yl,lx,r2s) !spherical factor Y( q+G )
        do lm =1,nlx
          l = ll(lm)
          pjyl(lm,ig1) = fpi*img**l *cy(lm)*yl(lm) *phase  *absqg1**l
ccccccccccccccccccc
c         write(6,*) ' pjyl=',lm,ig1,pjyl(lm,ig1)
ccccccccccccccccccc
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
      enddo

c rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
      rofi(1) = 0d0
      do ir   = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      do l = 0,lx
        rkpr(1:nr,l) = rofi(1:nr)**(l      +1 )
        rkmr(2:nr,l) = rofi(2:nr)**(-l-1   +1 )
        rkmr(1,l)    = rkmr(2,l)
      enddo

      print *,' mkjp: rojp'
c rojp
      do ig1 = 1,ngc
        call wronkj( absqg(ig1)**2, 0d0, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
        do lm = 1,nlx
          l = ll(lm)
          rojp(ig1,lm) = (-fjj(l))* pjyl(lm,ig1)
        enddo
      enddo

c ajr
      do ig1 = 1,ngc
        do ir  = 1,nr
          call bessl(absqg(ig1)**2*rofi(ir)**2,lx,phi,psi)
          do l   = 0,lx
            ajr(ir,l,ig1) = phi(l)* rofi(ir) **(l +1 )
            ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
          enddo
ccccccccccccccccccccccccccc
c        write(116,'(i3,10d13.6)') ir, rofi(ir), ajr(ir,0:lx,ig1)
ccccccccccccccccccccccccccc
        enddo
ccccccccccccccccccccccccccc
c        write(6,*) ig1,sum(ajr(1:nr,0:lx,ig1))
ccccccccccccccccccccccccccc
      enddo

c-------------------------
c      if(newsg) then  ! New version Feb 2001.
      do ig1 = 1,ngc
        call sigintAn1( absqg(ig1), lx, rofi, nr
     o                  ,a1(1:nr, 0:lx,ig1) )
      enddo
c      else            ! Old version
c        do l   = 0,lx
c        do ig1 = 1,ngc
c        call sigintA(rkpr(1:nr,0:lx),rkmr(1:nr,0:lx),lx,a,b,nr,ajr(1:nr,l,ig1),a1(1:nr,l,ig1)
c     &            ,l, rofi)
c        enddo
c        enddo
c      endif

c-------------------------
      print *,' mkjp: sgpb'
c sgpb
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)                      ! r jl        , r B(r)
c        call sigint(rkpr,rkmr, lx,a,b,nr,  ajr(1:nr,l,ig1),rprodx(1,n,l)
c     &            ,l, rofi, radint)
            call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
            sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
ccccccccccccccccccc
          enddo
        enddo
      enddo


c-----------------------------------------------------------------------
c sgpp
      print *,' mkjp: sgpp lx=',lx
      do ig1 = 1,ngc
        do ig2 = 1,ig1
          if(newsg) then ! NewVersion
            call sigintpp( absqg(ig1), absqg(ig2), lx, rmax,
     o             sigx)
            do l = 0,lx
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo
          else           ! OldVersion
            do l = 0,lx
              call gintxx(a1(1,l,ig1), ajr(1,l,ig2),A,B,NR, sig1 )
              call gintxx(a1(1,l,ig2), ajr(1,l,ig1),A,B,NR, sig2 )
              sig = (sig1+sig2)/2d0
              radsig(l) = fpi/(2*l+1) *sig
c         radsigo(l) = fpi/(2*l+1) *sig
            enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do l = 0,lx
c         write(6,"( ' l ig1 ig2 sig=',3i5,2d14.6,' diff= ',d14.6)")
c     &     l,ig1,ig2, radsigo(l),radsig(l), radsigo(l)-radsig(l)
c        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do lm = 1,nlx
            l = ll(lm)
            sgpp(ig1,ig2,lm)= dconjg(pjyl(lm,ig1))*radsig(l)* pjyl(lm,ig2)
            sgpp(ig2,ig1,lm)= dconjg(sgpp(ig1,ig2,lm))
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write(6,"('ig1 ig2 lm sgpp=', 3i5,12d13.6)") ig1,ig2,lm,
c     &        dconjg(pjyl(lm,ig1)), radsig(l), pjyl(lm,ig2),sgpp(ig1,ig2,lm)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
          enddo
        enddo
      enddo

c Fourier
c fouvb
      print *,' mkjp: Four'
      fouvb=0d0
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)
            call gintxx(ajr(1,l,ig1), rprodx(1,n,l), a,b,nr,
     o                 radint )
            fouvb(ig1, n, lm) =
     &    fpi/(absqg(ig1)**2) *dconjg(pjyl(lm,ig1))*radint

cccccccccccccccccccccccccccccccccccccc
c      if( abs(fouvb(ig1, n, lm))>1d4 ) then
c        write(6,*)' large fouvb=',ig1,n,lm,fouvb(ig1, n, lm)
c      endif
cccccccccccccccccccccccccccccccccccccc

          enddo
        enddo
      enddo
ccccccccccccccccccccc
c        write(6,*)' fourvb sum=',sum (fouvb)
ccccccccccccccccccccc

c fouvp
      fouvp=0d0
      do ig1 = 1,ngc
        do ig2 = 1,ngc
          absqg1=absqg(ig1)
          absqg2=absqg(ig2)
          call wronkj( absqg1**2, absqg2**2, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
          do lm = 1,nlx
            l = ll(lm)
cccccccccccccccccccccccccccccccccccc
c        call gintxx(ajr(1:nr,l,ig1), ajr(1:nr,l,ig2), a,b,nr,
c     o                 radint )
c        if(abs(ig1-ig2)<5) then
c        write(116,'(3i3,2d13.6)') ig1,ig2,l, radint, radint- (-fjj(l))
c        endif
cccccccccccccccccccccccccccccccccccccc
            fouvp(ig1, ig2, lm) = fpi/(absqg1**2) *dconjg(pjyl(lm,ig1))
     &    * (-fjj(l)) * pjyl(lm,ig2)

cccccccccccccccccccccccccccccccccccccc
c      if( abs(fouvp(ig1, ig2, lm))>1d4 ) then
c        write(6,*)' large fourvp=',ig1,ig2,fouvp(ig1, ig2, lm)
c      endif
cccccccccccccccccccccccccccccccccccccc
          enddo
        enddo
      enddo

cccccccccccccccccccccc1
c        write(6,*)' fouvp sum=',sum (fouvp)
cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
c      do ig1 = 1,ngc
c      do lm = 1,nlx
c        write(6,*)' fourvp=',ig1, lm, fouvp(ig1, ig1, lm)
c      enddo
c      enddo
cccccccccccccccccccccccccccccccccccccccc
      deallocate(ajr,a1,   rkpr,rkmr, qg,absqg,   pjyl)
      if (allocated( cy )) deallocate(cy)
      if (allocated( yl )) deallocate(yl)
      end




c-----------------------------------------------------------------------------------
      subroutine mkjp( q,ngc,ngvecc, alat, qlat, lxx,lx,nxx,nx,  
     &                   bas, a,b,rmax,nr,nrx,rprodx,
     o        rojp,sgpb,sgpp,fouvb,fouvp,ibas)
C- Make integrals in each MT. and the Fourier matrix.
Cr the integrals rojp, fouvb,fouvp
Cr are for  J_L(r)= j_l(sqrt(e) r)/sqrt(e)**l Y_L,
Cr which behaves as r^l/(2l+1)!! near r=0.
Cr
Cr oniste integral is based on
Cr 1/|r-r'| = \sum 4 pi /(2k+1) \frac{r_<^k }{ r_>^{k+1} } Y_L(r) Y_L(r')
Cr See PRB34 5512(1986) for sigma type integral
Cr
      implicit none
      integer(4) :: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx
      real(8)    :: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat,
     &              qlat(3,3)
ci rho-type onsite integral
      complex(8) :: rojp(ngc, (lxx+1)**2)
ci sigma-type onsite integral
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2),  
     &              sgpp(ngc,  ngc,  (lxx+1)**2)
      real(8),allocatable::cy(:),yl(:)
ci Fourier
      complex(8) :: 
     &              fouvb(ngc,  nxx, (lxx+1)**2),
     &              fouvp(ngc,  ngc, (lxx+1)**2)
c internal
      integer(4) :: nlx,ig1,ig2,l,n,ir,n1,n2,lm, ibas
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      real(8)    :: pi,fpi,tpiba, qg1(3),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1,absqg2,
     &  fac,radint,radsigo(0:lx),radsig(0:lx),phi(0:lx),psi(0:lx)
     &   ,rofi(nr),r2s,sig,sig1,sig2,sigx(0:lx),sig0(0:lx) ,qg2(3)
      complex(8) :: img =(0d0,1d0),phase
      complex(8),allocatable :: pjyl(:,:)
      real(8),allocatable ::ajr(:,:,:),a1(:,:,:),rkpr(:,:),rkmr(:,:),
     &  qg(:,:),absqg(:)
      logical  :: newsg=.true. !new version (Mar2001) for sgpp sgpb or not
c-------------------------------------------------
      print *,' mkjp:'
      nlx = (lx+1)**2
      allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),
     &  rkpr(nr,0:lx),rkmr(nr,0:lx), qg(3,ngc),absqg(ngc),
     &  pjyl((lx+1)**2,ngc) )
c
      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
      tpiba = 2*pi/alat
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
c      print *,' mkjp: end of sylmnc'
c q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
      do ig1 = 1,ngc
        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        qg1(1:3) = qg(1:3,ig1)
        absqg(ig1)  = sqrt(sum(qg1(1:3)**2))
        absqg1   = absqg(ig1)
cccccccccccccccccccccccccccccccccccccccccccccc
c        write(6,'(i4,3f12.6,3x,f12.6)')ig1, qg1,absqg1
cccccccccccccccccccccccccccccccccccccccccccccc
        phase = exp( img*sum(qg1(1:3)*bas(1:3))*alat  )
        call sylm(qg1/absqg1,yl,lx,r2s) !spherical factor Y( q+G )
        do lm =1,nlx
          l = ll(lm)
          pjyl(lm,ig1) = fpi*img**l *cy(lm)*yl(lm) *phase  *absqg1**l
ccccccccccccccccccc
c         write(6,*) ' pjyl=',lm,ig1,pjyl(lm,ig1)
ccccccccccccccccccc
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
      enddo

c rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
      rofi(1) = 0d0
      do ir   = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      do l = 0,lx
        rkpr(1:nr,l) = rofi(1:nr)**(l      +1 )
        rkmr(2:nr,l) = rofi(2:nr)**(-l-1   +1 )
        rkmr(1,l)    = rkmr(2,l)
      enddo

      print *,' mkjp: rojp'
c rojp
      do ig1 = 1,ngc
        call wronkj( absqg(ig1)**2, 0d0, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
        do lm = 1,nlx
          l = ll(lm)
          rojp(ig1,lm) = (-fjj(l))* pjyl(lm,ig1)
        enddo
      enddo

c ajr
      do ig1 = 1,ngc
        do ir  = 1,nr
          call bessl(absqg(ig1)**2*rofi(ir)**2,lx,phi,psi)
          do l   = 0,lx
            ajr(ir,l,ig1) = phi(l)* rofi(ir) **(l +1 )
            ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
          enddo
ccccccccccccccccccccccccccc
c        write(116,'(i3,10d13.6)') ir, rofi(ir), ajr(ir,0:lx,ig1)
ccccccccccccccccccccccccccc
        enddo
ccccccccccccccccccccccccccc
c        write(6,*) ig1,sum(ajr(1:nr,0:lx,ig1))
ccccccccccccccccccccccccccc
      enddo

c-------------------------
      if(newsg) then  ! New version Feb 2001.
        do ig1 = 1,ngc
          call sigintAn1( absqg(ig1), lx, rofi, nr
     o                  ,a1(1:nr, 0:lx,ig1) )
        enddo
      else            ! Old version
        do l   = 0,lx
          do ig1 = 1,ngc
            call sigintA(rkpr,rkmr,lx,a,b,nr,ajr(1:nr,l,ig1),a1(1:nr,l,ig1)
     &            ,l, rofi)
          enddo
        enddo
      endif

c-------------------------
      print *,' mkjp: sgpb'
c sgpb
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)                      ! r jl        , r B(r)
c        call sigint(rkpr,rkmr, lx,a,b,nr,  ajr(1:nr,l,ig1),rprodx(1,n,l)
c     &            ,l, rofi, radint)
            call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
            sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
ccccccccccccccccccc
          enddo
        enddo
      enddo


c-----------------------------------------------------------------------
c sgpp
      print *,' mkjp: sgpp lx=',lx
      do ig1 = 1,ngc
        do ig2 = 1,ig1
          if(newsg) then ! NewVersion
            call sigintpp( absqg(ig1), absqg(ig2), lx, rmax,
     o             sigx)
            do l = 0,lx
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo
          else           ! OldVersion
            do l = 0,lx
              call gintxx(a1(1,l,ig1), ajr(1,l,ig2),A,B,NR, sig1 )
              call gintxx(a1(1,l,ig2), ajr(1,l,ig1),A,B,NR, sig2 )
              sig = (sig1+sig2)/2d0
              radsig(l) = fpi/(2*l+1) *sig
c         radsigo(l) = fpi/(2*l+1) *sig
            enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do l = 0,lx
c         write(6,"( ' l ig1 ig2 sig=',3i5,2d14.6,' diff= ',d14.6)")
c     &     l,ig1,ig2, radsigo(l),radsig(l), radsigo(l)-radsig(l)
c        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do lm = 1,nlx
            l = ll(lm)
            sgpp(ig1,ig2,lm)= dconjg(pjyl(lm,ig1))*radsig(l)* pjyl(lm,ig2)
            sgpp(ig2,ig1,lm)= dconjg(sgpp(ig1,ig2,lm))
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write(6,"('ig1 ig2 lm sgpp=', 3i5,12d13.6)") ig1,ig2,lm,
c     &        dconjg(pjyl(lm,ig1)), radsig(l), pjyl(lm,ig2),sgpp(ig1,ig2,lm)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
          enddo
        enddo
      enddo

c Fourier
c fouvb
      print *,' mkjp: Four'
      fouvb=0d0
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)
            call gintxx(ajr(1,l,ig1), rprodx(1,n,l), a,b,nr,
     o                 radint )
            fouvb(ig1, n, lm) =
     &    fpi/(absqg(ig1)**2) *dconjg(pjyl(lm,ig1))*radint

cccccccccccccccccccccccccccccccccccccc
c      if( abs(fouvb(ig1, n, lm))>1d4 ) then
c        write(6,*)' large fouvb=',ig1,n,lm,fouvb(ig1, n, lm)
c      endif
cccccccccccccccccccccccccccccccccccccc

          enddo
        enddo
      enddo
ccccccccccccccccccccc
c        write(6,*)' fourvb sum=',sum (fouvb)
ccccccccccccccccccccc

c fouvp
      fouvp=0d0
      do ig1 = 1,ngc
        do ig2 = 1,ngc
          absqg1=absqg(ig1)
          absqg2=absqg(ig2)
          call wronkj( absqg1**2, absqg2**2, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
          do lm = 1,nlx
            l = ll(lm)
cccccccccccccccccccccccccccccccccccc
c        call gintxx(ajr(1:nr,l,ig1), ajr(1:nr,l,ig2), a,b,nr,
c     o                 radint )
c        if(abs(ig1-ig2)<5) then
c        write(116,'(3i3,2d13.6)') ig1,ig2,l, radint, radint- (-fjj(l))
c        endif
cccccccccccccccccccccccccccccccccccccc
            fouvp(ig1, ig2, lm) = fpi/(absqg1**2) *dconjg(pjyl(lm,ig1))
     &    * (-fjj(l)) * pjyl(lm,ig2)

cccccccccccccccccccccccccccccccccccccc
c      if( abs(fouvp(ig1, ig2, lm))>1d4 ) then
c        write(6,*)' large fourvp=',ig1,ig2,fouvp(ig1, ig2, lm)
c      endif
cccccccccccccccccccccccccccccccccccccc
          enddo
        enddo
      enddo

cccccccccccccccccccccc1
c        write(6,*)' fouvp sum=',sum (fouvp)
cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
c      do ig1 = 1,ngc
c      do lm = 1,nlx
c        write(6,*)' fourvp=',ig1, lm, fouvp(ig1, ig1, lm)
c      enddo
c      enddo
cccccccccccccccccccccccccccccccccccccccc
      deallocate(ajr,a1,   rkpr,rkmr, qg,absqg,   pjyl)
      if (allocated( cy )) deallocate(cy)
      if (allocated( yl )) deallocate(yl)
      end



c------------------------------------------
      subroutine mkjb( lxx,lx,nxx,nx,  
     &                 a,b,nr,nrx,rprodx,
     o        rojb,sgbb)
C--make integrals in each MT. and the Fourier matrix.
      implicit none
      integer(4) :: lxx, lx, nxx, nx(0:lxx),nr,nrx
      real(8)    :: q(3), rprodx(nrx,nxx,0:lxx),a,b
ci rho-type onsite integral
      real(8)    :: rojb(nxx, 0:lxx)
ci sigma-type onsite integral
      real(8)    :: sgbb(nxx, nxx, 0:lxx)
c internal
      integer(4) :: l,n,ir,n1,n2,l1
      real(8)    :: 
     &  fac,
     &  rofi(nr),xxx
      real(8),allocatable :: rkpr(:,:),rkmr(:,:)
c
      allocate(rkpr(nr,0:lx),rkmr(nr,0:lx))
c-------------------------------------------------
c rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
ccccccccccccccccccccccccccccccccc
c      do l = 0,lx
c      do n = 1,nx(l)
c      do n1 = 1,nx(l)
c        call gintxx(rprodx(1:nr,n,l), rprodx(1:nr,n1,l), a,b,nr,
c     o                 xxx )
c      write(6,*)' check rprodx =',l,n,n-n1,xxx
c      enddo
c      enddo
c      enddo
c      stop 'xxx'
cccccccccccccccccccccccccccccccc

ccccccccccccc
      rojb=0d0
      sgbb=0d0
ccccccccccccc
      rofi(1)    = 0d0
      do ir      = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      do l = 0,lx
        rkpr(1:nr,l) = rofi(1:nr)**(l +1)
        rkmr(2:nr,l) = rofi(2:nr)**(-l-1) *rofi(2:nr)
        rkmr(1,l)    = rkmr(2,l)
      enddo
c rojb
      fac = 1d0
      do l = 0,lx
        fac = fac/(2*l+1)
        do n = 1,nx(l)
          call gintxx(rkpr(1,l), rprodx(1,n,l), a,b,nr,
     o                 rojb(n,l) )
        enddo
        rojb(1:nx(l),l) = fac*rojb(1:nx(l),l)
      enddo
c sgbb
      do l  = 0,lx
        do n1 = 1,nx(l)
          do n2 = 1,nx(l)
            call sigint(rkpr,rkmr, lx,a,b,nr,rprodx(1,n1,l),rprodx(1,n2,l)
     &            ,l, rofi, sgbb(n1, n2, l))
          enddo
        enddo
      enddo
c      write(6,*) ' rojbsum=', sum(rojb(:,:)),   sum(abs(rojb(:,:)))
c      write(6,*) ' sgbbsum=', sum(sgbb(:,:,:)), sum(abs(sgbb(:,:,:)))

ccccccccccccccccccccccccccccccccccccc
c      write(6,*)' sigint 1 1 0=',sgbb(1, 1, 0) !/(16d0*datan(1d0))
c      sgbb(1, 1, 0) =0d0
cccccccccccccccccccccccccccccccccccccccc
      deallocate(rkpr,rkmr)


      end

c-----------------------------------------------------------
      subroutine vcoulq2(q,nbloch, ngc,  
     &                  nbas, lx,lxx, nx,nxx,
     &                  alat, qlat, vol, ngvecc,
     &          strx,rojp,rojb, sgbb,sgpb, fouvb,   !sgpp,fouvp,
     i          nblochpmx,bas,rmax,
     o          vcoul)
Co Coulmb matrix for each q. -----------------------------------------------
Ci strx:  Structure factors
Ci nlx corresponds to (lx+1)**2 . lx corresponds to 2*lmxax.
Ci rho-type integral
Ci  ngvecc     : q+G vector
Ci  rojp rojb  : rho-type integral
ci  sigma-type onsite integral
ci  Fourier
Ci  nx(l,ibas) : max number of radial function index for each l and ibas.
Ci               Note that the definition is a bit different from nx in basnfp.
ci  nxx        : max number of nx among all l and ibas.
ci  lx(nbas)   : max number of l for each ibas.
ci  lxx        :
ci
ci  vol : cell vol
c
Co Vcoul
cr vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
cr    vcoul = 2d0*vcoul !  in Ry unit.
c---------------------------------------------------------------------------
c  rojp = <j_aL(r) | P(q+G)_aL > where
c         |P(q+G)_aL> : the projection of exp(i (q+G) r) to aL channnel.
c         |j_aL>      : \def r^l/(2l+1)!! Y_L.  The spherical bessel functions near r=0.  Energy-dependence is omitted.
c
      implicit none
      integer(4) :: nbloch, nblochpmx, nbas, 
     &              lxx,lx(nbas), nxx, nx(0:lxx,nbas)
      real(8)    :: egtpi,vol,q(3),fpi

ci structure con
      complex(8) :: strx((lxx+1)**2, nbas, (lxx+1)**2,nbas)
ci |q+G|**2
      integer(4) :: ngc, ngvecc(3,ngc)
      real(8)    :: qlat(3,3),alat,absqg2(ngc),qg(3)

ci rho-type onsite integral
      complex(8) ::   rojp(ngc, (lxx+1)**2, nbas)
      real(8)    ::   rojb(nxx, 0:lxx, nbas)
ci sigma-type onsite integral
      real(8)    :: sgbb(nxx,  nxx,  0:lxx,      nbas)
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2, nbas)
c     &              ,sgpp(ngc,  ngc,  (lxx+1)**2, nbas)
ci Fourier
      complex(8) :: 
     &             fouvb(ngc,  nxx, (lxx+1)**2, nbas)
Co
     &             ,vcoul(nblochpmx, nblochpmx)
c     &             ,fouvp(ngc,  ngc, (lxx+1)**2, nbas)

cinternals
      integer(4) :: ibl1, ibl2,ig1,ig2,ibas,ibas1,ibas2, 
     &              l,m,n, n1,l1,m1,lm1,n2,l2,m2,lm2,ipl1,ipl2
      integer(4) :: ibasbl(nbloch), nbl(nbloch), lbl(nbloch),
     &              mbl(nbloch), lmbl(nbloch)
      real(8) :: pi, fpivol,tpiba
      complex(8) :: rojpstrx((lxx+1)**2,nbas)

c check
      complex(8),allocatable :: hh(:,:),oo(:,:),zz(:,:)
      real(8),allocatable    :: eb(:)

      complex(8),allocatable :: matp(:),matp2(:)
      complex(8) :: xxx
      integer(4) :: nblochngc,nev,nmx,ix
      logical :: ptest=.false. ! See ptest in basnfp.f

c-------------------
      real(8),   allocatable :: cy(:),yl(:)
      complex(8),allocatable :: pjyl_(:,:),phase(:,:)
      complex(8) :: img=(0d0,1d0)
      real(8):: bas(3,nbas),r2s,rmax(nbas)
      integer(4):: lm
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      real(8)::  fkk(0:lxx),fkj(0:lxx),fjk(0:lxx),fjj(0:lxx),sigx(0:lxx),radsig(0:lxx)
      complex(8):: fouvp_ig1_ig2, fouvp_ig2_ig1, sgpp_ig1_ig2
c---------------------------------------------------------------
      write(6,'(" vcoul: nblochpmx  nbloch ngc=",3i6)') nblochpmx,nbloch,ngc
c     print *, ' sum fouvp=',sum(fouvp(:,:,:,1))
c     print *, ' sum fouvb=',sum(fouvb(:,:,:,1))
      pi    = 4d0*datan(1d0)
      fpi    = 4*pi
      fpivol = 4*pi*vol

c---for sgpp fouvp
      allocate(  !ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),rkpr(nr,0:lx),rkmr(nr,0:lx), 
     &  pjyl_((lxx+1)**2,ngc),phase(ngc,nbas) )
      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
      call sylmnc(cy,lxx)

c=======================================================
      vcoul = 0d0
c-gvec
      tpiba = 2*pi/alat
      do ig1 = 1,ngc
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        absqg2(ig1)  = sum(qg(1:3)**2)
c---for spgg fourvp ----------
        do ibas=1,nbas
          phase(ig1,ibas) = exp( img*sum(qg(1:3)*bas(1:3,ibas))*alat  )
        enddo
        call sylm(qg/sqrt(absqg2(ig1)),yl,lxx,r2s) !spherical factor Y( q+G )
        do lm =1,(lxx+1)**2
          l = ll(lm)
          pjyl_(lm,ig1) = fpi*img**l *cy(lm)*yl(lm)  * sqrt(absqg2(ig1))**l  !*phase
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
c--------------
      enddo
c


c-- index (mx,nx,lx,ibas) order.
      ibl1 = 0
      do ibas= 1, nbas
        do l   = 0, lx(ibas)
c        write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
          do n   = 1, nx(l,ibas)
            do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              mbl   (ibl1) = m
              lmbl  (ibl1) = l**2 + l+1 +m
c        write(6,*)ibl1,n,l,m,lmbl(ibl1)
            enddo
          enddo
        enddo
      enddo
      if(ibl1/= nbloch) then
        write(6,*)' ibl1 nbloch',ibl1, nbloch
Cstop2rx 2013.08.09 kino        stop ' vcoulq: error ibl1/= nbloch'
        call rx( ' vcoulq: error ibl1/= nbloch')
      endif
c-- <B|v|B> block
c      write(6,*)' vcoulq: bvb block xxx rojbsum='
c      write(6,*) sum(rojb(:,:,1))
c      write(6,*) sum(rojb(:,:,2))
c      write(6,*) sum(rojb(:,:,3))
c      write(6,*) sum(rojb(:,:,4))
c      write(6,*)' vcoulq: bvb block xxx sgbbbsum='
c      write(6,*) sum(sgbb(:,:,:,1))
c      write(6,*) sum(sgbb(:,:,:,2))
c      write(6,*) sum(sgbb(:,:,:,3))
c      write(6,*) sum(sgbb(:,:,:,4))
      do ibl1= 1, nbloch
        ibas1= ibasbl(ibl1)
        n1   = nbl (ibl1)
        l1   = lbl (ibl1)
        m1   = mbl (ibl1)
        lm1  = lmbl(ibl1)
        do ibl2= 1, ibl1
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          vcoul(ibl1,ibl2) =
     &     rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2)
     &    *rojb(n2, l2, ibas2)
          if(ibas1==ibas2 .and. lm1==lm2) then
            vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1)
            ! sigma-type contribution. onsite coulomb
          endif
        enddo
      enddo

c <P_G|v|B>
      write(6,*)' vcoulq: pgvb block'
      do ibl2= 1, nbloch
        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        m2   = mbl (ibl2)
        lm2  = lmbl(ibl2)
        do ig1 = 1,ngc
          ipl1 = nbloch + ig1
          vcoul(ipl1,ibl2) = fouvb(ig1,  n2, lm2, ibas2)

          do ibas1= 1, nbas
            do lm1  = 1, (lx(ibas1)+1)**2
              vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
     &      *rojb(n2, l2, ibas2)
              if(ibas1==ibas2 .and.lm1==lm2) then
                vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &        sgpb(ig1, n2, lm2, ibas2)
              endif
            enddo
          enddo
        enddo
      enddo

c <P_G|v|P_G>
      write(6,*)' vcoulq: pgvpg block'
      do ig1 = 1,ngc
        ipl1 = nbloch + ig1
        rojpstrx = 0d0
        do ibas1= 1, nbas
          do lm1  = 1, (lx(ibas1)+1)**2
            do ibas2= 1, nbas
              do lm2  = 1, (lx(ibas2)+1)**2
                rojpstrx(lm2, ibas2) = rojpstrx(lm2, ibas2)+
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
              enddo
            enddo
          enddo
        enddo
        do ig2 = 1,ig1
          ipl2 = nbloch + ig2

          if(ig1==ig2) then
            vcoul(ipl1,ipl2) = fpivol/absqg2(ig1)
          endif

          do ibas2= 1, nbas
c... for fouvp and sgpp -------
            call wronkj( absqg2(ig1), absqg2(ig2), rmax(ibas2),lx(ibas2),
     o                 fkk,fkj,fjk,fjj)
            call sigintpp( absqg2(ig1)**.5, absqg2(ig2)**.5, lx(ibas2), rmax(ibas2),
     o             sigx)
            do l = 0,lx(ibas2)
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo
c------------------------------
            do lm2  = 1, (lx(ibas2)+1)**2
              l= ll(lm2)
c...fouvp sgpp-----------
              fouvp_ig1_ig2 = fpi/absqg2(ig1) *dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2))
     &      * (-fjj(l)) * pjyl_(lm2,ig2)*phase(ig2,ibas2)
              fouvp_ig2_ig1 = fpi/absqg2(ig2) *dconjg(pjyl_(lm2,ig2)*phase(ig2,ibas2))
     &      * (-fjj(l)) * pjyl_(lm2,ig1)*phase(ig1,ibas2)
              sgpp_ig1_ig2  = dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2))*radsig(l)
     &                      * pjyl_(lm2,ig2)*phase(ig2,ibas2)
c------------------------
              vcoul(ipl1,ipl2) = vcoul(ipl1,ipl2) 
     &       +  rojpstrx(lm2,ibas2)*rojp(ig2, lm2, ibas2)
c     &      -  dconjg( fouvp(ig2,  ig1, lm2, ibas2)) !BugFix Mar5-01 It was dcmplx.
c     &      -          fouvp(ig1,  ig2, lm2, ibas2)
c     &      +  sgpp(ig1, ig2, lm2, ibas2)
     &       -  dconjg( fouvp_ig2_ig1 )
     &       -          fouvp_ig1_ig2
     &       +  sgpp_ig1_ig2
            enddo
          enddo
        enddo
      enddo

c-- Right-Upper part of vcoul.
      write(6,*)' vcoulq: right-upper'
      do ipl1=1, nbloch+ngc
        do ipl2=1, ipl1-1
          vcoul(ipl2,ipl1) = dconjg(vcoul(ipl1,ipl2))
        enddo
      enddo

c vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
c     vcoul = 2d0*vcoul !  in Ry unit.
c

c check write
      do ix = 1,nbloch+ngc,20
        write(6,"(' Diagonal Vcoul =',i5,2d18.10)") ix,vcoul(ix,ix)
      enddo

C... deallocate
      if( allocated(yl)   ) deallocate(yl)
      if( allocated(cy)   ) deallocate(cy)
      if( allocated(phase)) deallocate(phase)
      if( allocated(pjyl_)) deallocate(pjyl_)


c-------------------------------------------------------
      if(.not.ptest) return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! Below ia a plane-wave test.
c--- check! Coulomb by plane wave expansion.
      write(6,*) ' --- plane wave Coulomb matrix check 1---- '
      write(197,*) ' --- off diagonal ---- '
      nblochngc = nbloch+ngc
      allocate(matp(nblochngc),matp2(nblochngc))
      do ig1 = 1,ngc
        matp = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp(ibl2) = fouvb(ig1, n2, lm2, ibas2)*absqg2(ig1)/fpi
        enddo
        matp(nbloch+ig1) = 1d0
        ig2=ig1
c      do ig2 = 1,ngc !off diagnal
        matp2 = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp2(ibl2) = fouvb(ig2, n2, lm2, ibas2)*absqg2(ig2)/fpi
        enddo
        matp2(nbloch+ig2) = 1d0
        xxx= sum(
     &    matmul(matp(1:nblochngc),vcoul(1:nblochngc,1:nblochngc))
     &             *dconjg(matp2(1:nblochngc))  )
        if(ig1/=ig2) then  !off diagnal
          if(abs(xxx)>1d-1 ) then
            write(197,'(2i5, 2d13.6)') ig1,ig2, xxx
            write(197,'("    matpp ", 2d13.6)') 
     &      vcoul(nbloch+ig1,nbloch+ig2)
            write(197,*)
          endif
        else
          write(196,'(2i5," exact=",3d13.6,"q ngsum=",3f8.4,i5)')
     &       ig1,ig2,fpi*vol/absqg2(ig1)
     &    , fpi*vol/absqg2(ig2),absqg2(ig1), q(1:3)
     &    , sum(ngvecc(1:3,ig1)**2)
          write(196,'("           cal  =", 2d13.6)') xxx
          write(196,'("           vcoud=", 2d13.6)') 
     &       vcoul(nbloch+ig1,nbloch+ig2)
          write(196,*)
        endif
c      enddo !off diagnal
      enddo
c
      deallocate(matp,matp2)
c      stop ' *** ptest end *** See fort.196 and 197'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end

















c-----------------------------------------------------------
      subroutine vcoulq_not_used(q,nbloch, ngc,  
     &                  nbas, lx,lxx, nx,nxx,
     &                  alat, qlat, vol, ngvecc,
     &          strx,rojp,rojb, sgbb,sgpb,sgpp, fouvb,fouvp,
     i          nblochpmx,
     o          vcoul)
Co Coulmb matrix for each q. -----------------------------------------------
Ci strx:  Structure factors
Ci nlx corresponds to (lx+1)**2 . lx corresponds to 2*lmxax.
Ci rho-type integral
Ci  ngvecc     : q+G vector
Ci  rojp rojb  : rho-type integral
ci  sigma-type onsite integral
ci  Fourier
Ci  nx(l,ibas) : max number of radial function index for each l and ibas.
Ci               Note that the definition is a bit different from nx in basnfp.
ci  nxx        : max number of nx among all l and ibas.
ci  lx(nbas)   : max number of l for each ibas.
ci  lxx        :
ci
ci  vol : cell vol
c
Co Vcoul
cr vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
cr    vcoul = 2d0*vcoul !  in Ry unit.
c---------------------------------------------------------------------------
c  rojp = <j_aL(r) | P(q+G)_aL > where
c         |P(q+G)_aL> : the projection of exp(i (q+G) r) to aL channnel.
c         |j_aL>      : \def r^l/(2l+1)!! Y_L.  The spherical bessel functions near r=0.  Energy-dependence is omitted.
c
      implicit none
      integer(4) :: nbloch, nblochpmx, nbas, 
     &              lxx,lx(nbas), nxx, nx(0:lxx,nbas)
      real(8)    :: egtpi,vol,q(3),fpi

ci structure con
      complex(8) :: strx((lxx+1)**2, nbas, (lxx+1)**2,nbas)
ci |q+G|**2
      integer(4) :: ngc, ngvecc(3,ngc)
      real(8)    :: qlat(3,3),alat,absqg2(ngc),qg(3)

ci rho-type onsite integral
      complex(8) ::   rojp(ngc, (lxx+1)**2, nbas)
      real(8)    ::   rojb(nxx, 0:lxx, nbas)
ci sigma-type onsite integral
      real(8)    :: sgbb(nxx,  nxx,  0:lxx,      nbas)
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2, nbas),  
     &              sgpp(ngc,  ngc,  (lxx+1)**2, nbas)
ci Fourier
      complex(8) :: 
     &             fouvb(ngc,  nxx, (lxx+1)**2, nbas),
     &             fouvp(ngc,  ngc, (lxx+1)**2, nbas),
Co
     &           vcoul(nblochpmx, nblochpmx)

cinternals
      integer(4) :: ibl1, ibl2,ig1,ig2,ibas,ibas1,ibas2, 
     &              l,m,n, n1,l1,m1,lm1,n2,l2,m2,lm2,ipl1,ipl2
      integer(4) :: ibasbl(nbloch), nbl(nbloch), lbl(nbloch),
     &              mbl(nbloch), lmbl(nbloch)
      real(8) :: pi, fpivol,tpiba
      complex(8) :: rojpstrx((lxx+1)**2,nbas)

c check
      complex(8),allocatable :: hh(:,:),oo(:,:),zz(:,:)
      real(8),allocatable    :: eb(:)

      complex(8),allocatable :: matp(:),matp2(:)
      complex(8) :: xxx
      integer(4) :: nblochngc,nev,nmx,ix
      logical :: ptest=.false. ! See ptest in basnfp.f
c---------------------------------------------------------------
      write(6,'(" vcoul: nblochpmx  nbloch ngc=",3i6)') 
     & nblochpmx,nbloch,ngc
c     print *, ' sum fouvp=',sum(fouvp(:,:,:,1))
c     print *, ' sum fouvb=',sum(fouvb(:,:,:,1))
      pi    = 4d0*datan(1d0)
      fpi    = 4*pi
      fpivol = 4*pi*vol
      vcoul = 0d0
c-gvec
      tpiba = 2*pi/alat
      do ig1 = 1,ngc
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        absqg2(ig1)  = sum(qg(1:3)**2)
      enddo
c
c-- index (mx,nx,lx,ibas) order.
      ibl1 = 0
      do ibas= 1, nbas
        do l   = 0, lx(ibas)
c        write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
          do n   = 1, nx(l,ibas)
            do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              mbl   (ibl1) = m
              lmbl  (ibl1) = l**2 + l+1 +m
c        write(6,*)ibl1,n,l,m,lmbl(ibl1)
            enddo
          enddo
        enddo
      enddo
      if(ibl1/= nbloch) then
        write(6,*)' ibl1 nbloch',ibl1, nbloch
Cstop2rx 2013.08.09 kino        stop ' vcoulq: error ibl1/= nbloch'
        call rx( ' vcoulq: error ibl1/= nbloch')
      endif
c-- <B|v|B> block
c      write(6,*)' vcoulq: bvb block xxx rojbsum='
c      write(6,*) sum(rojb(:,:,1))
c      write(6,*) sum(rojb(:,:,2))
c      write(6,*) sum(rojb(:,:,3))
c      write(6,*) sum(rojb(:,:,4))
c      write(6,*)' vcoulq: bvb block xxx sgbbbsum='
c      write(6,*) sum(sgbb(:,:,:,1))
c      write(6,*) sum(sgbb(:,:,:,2))
c      write(6,*) sum(sgbb(:,:,:,3))
c      write(6,*) sum(sgbb(:,:,:,4))

CCC!$OMP parallel do private( ibl1, ibas1,n1,l1,m1,lm1,ibl2)
CCC!$OMP& private( ibas2, n2,l2, m2,lm2)
CCC!$OMP& reduction(+:vcoul)
      do ibl1= 1, nbloch
        ibas1= ibasbl(ibl1)
        n1   = nbl (ibl1)
        l1   = lbl (ibl1)
        m1   = mbl (ibl1)
        lm1  = lmbl(ibl1)
        do ibl2= 1, ibl1
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          vcoul(ibl1,ibl2) =
     &     rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2)
     &    *rojb(n2, l2, ibas2)
          if(ibas1==ibas2 .and. lm1==lm2) then
            vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1)
            ! sigma-type contribution. onsite coulomb
          endif
        enddo
      enddo
CCC!$OMP end parallel do

c <P_G|v|B>
      write(6,*)' vcoulq: pgvb block'

CCC!$OMP parallel do private( ibl2, ibas2, n2,l2,m2, lm2, ig1,ipl1, ibas1, lm1 )
CCC!$OMP& reduction(+:vcoul)
      do ibl2= 1, nbloch
        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        m2   = mbl (ibl2)
        lm2  = lmbl(ibl2)
        do ig1 = 1,ngc
          ipl1 = nbloch + ig1
          vcoul(ipl1,ibl2) = fouvb(ig1,  n2, lm2, ibas2)

          do ibas1= 1, nbas
            do lm1  = 1, (lx(ibas1)+1)**2
              vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
     &      *rojb(n2, l2, ibas2)
              if(ibas1==ibas2 .and.lm1==lm2) then
                vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &        sgpb(ig1, n2, lm2, ibas2)
              endif
            enddo
          enddo
        enddo
      enddo
CCC!$OMP end parallel do

c <P_G|v|P_G>
      write(6,*)' vcoulq: pgvpg block'


      do ig1 = 1,ngc
        ipl1 = nbloch + ig1
        rojpstrx = 0d0
        do ibas1= 1, nbas
          do lm1  = 1, (lx(ibas1)+1)**2
CCC!$OMP parallel do  private(ibas2,lm2)  shared(rojpstrx)
            do ibas2= 1, nbas
              do lm2  = 1, (lx(ibas2)+1)**2
                rojpstrx(lm2, ibas2) = rojpstrx(lm2, ibas2)+
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
              enddo
            enddo
          enddo
        enddo

CCC!$OMP parallel do private(ig2,ipl2,ibas2,lm2)  shared(vcoul)
        do ig2 = 1,ig1
          ipl2 = nbloch + ig2

          if(ig1==ig2) then
            vcoul(ipl1,ipl2) = fpivol/absqg2(ig1)
          endif
c
          do ibas2= 1, nbas
            do lm2  = 1, (lx(ibas2)+1)**2
              vcoul(ipl1,ipl2) = vcoul(ipl1,ipl2) 
     &      -  dconjg( fouvp(ig2,  ig1, lm2, ibas2)) !BugFix Mar5-01 It was dcmplx.
     &      -          fouvp(ig1,  ig2, lm2, ibas2)
     &      +  rojpstrx(lm2,ibas2)*rojp(ig2, lm2, ibas2)
     &      +  sgpp(ig1, ig2, lm2, ibas2)
            enddo
          enddo
        enddo

      enddo ! ig1


c-- Right-Upper part of vcoul.
      write(6,*)' vcoulq: right-upper'
      do ipl1=1, nbloch+ngc
        do ipl2=1, ipl1-1
          vcoul(ipl2,ipl1) = dconjg(vcoul(ipl1,ipl2))
        enddo
      enddo

c vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
c     vcoul = 2d0*vcoul !  in Ry unit.
c

c check write
      do ix = 1,nbloch+ngc,20
        write(6,"(' Diagonal Vcoul =',i5,2d18.10)") ix,vcoul(ix,ix)
      enddo

c-------------------------------------------------------
      if(.not.ptest) return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! Below ia a plane-wave test.
c--- check! Coulomb by plane wave expansion.
      write(6,*) ' --- plane wave Coulomb matrix check 1---- '
      write(197,*) ' --- off diagonal ---- '
      nblochngc = nbloch+ngc
      allocate(matp(nblochngc),matp2(nblochngc))
      do ig1 = 1,ngc
        matp = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp(ibl2) = fouvb(ig1, n2, lm2, ibas2)*absqg2(ig1)/fpi
        enddo
        matp(nbloch+ig1) = 1d0
        ig2=ig1
c      do ig2 = 1,ngc !off diagnal
        matp2 = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp2(ibl2) = fouvb(ig2, n2, lm2, ibas2)*absqg2(ig2)/fpi
        enddo
        matp2(nbloch+ig2) = 1d0
        xxx= sum(
     &    matmul(matp(1:nblochngc),vcoul(1:nblochngc,1:nblochngc))
     &             *dconjg(matp2(1:nblochngc))  )
        if(ig1/=ig2) then  !off diagnal
          if(abs(xxx)>1d-1 ) then
            write(197,'(2i5, 2d13.6)') ig1,ig2, xxx
            write(197,'("    matpp ", 2d13.6)') 
     &      vcoul(nbloch+ig1,nbloch+ig2)
            write(197,*)
          endif
        else
          write(196,'(2i5," exact=",3d13.6,"q ngsum=",3f8.4,i5)')
     &       ig1,ig2,fpi*vol/absqg2(ig1)
     &    , fpi*vol/absqg2(ig2),absqg2(ig1), q(1:3)
     &    , sum(ngvecc(1:3,ig1)**2)
          write(196,'("           cal  =", 2d13.6)') xxx
          write(196,'("           vcoud=", 2d13.6)') 
     &       vcoul(nbloch+ig1,nbloch+ig2)
          write(196,*)
        endif
c      enddo !off diagnal
      enddo
c
      deallocate(matp,matp2)
c      stop ' *** ptest end *** See fort.196 and 197'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end

c-------------------------------------------------------------
      subroutine sigintA(rkp,rkm,kmx,a,b,nr,phi1,a1,k,rofi)
      implicit none
      integer(4) :: nr,kmx,k,ir
      real(8):: a,b, a1(nr),a2(nr),b1(nr),rkp(nr,0:kmx),rkm(nr,0:kmx),
     &    int1x(nr),int2x(nr), phi1(nr), phi2(nr),rofi(nr),sig
      real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
c
      a1(1) = 0d0;  a1(2:nr) = rkp(2:nr,k)
      a2(1) = 0d0;  a2(2:nr) = rkm(2:nr,k)
      b1(1:nr) = phi1(1:nr)
      call intn_smpxxx(a1,b1,int1x,a,b,rofi,nr,0)
      call intn_smpxxx(a2,b1,int2x,a,b,rofi,nr,0)
c
      a1(1) = 0d0; a1(2:nr) = 
     &  rkm(2:nr,k) *( int1x(1)-int1x(2:nr) )+ rkp(2:nr,k) * int2x(2:nr)
      end


