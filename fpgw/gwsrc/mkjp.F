      subroutine vcoulq_4(q,nbloch, ngc,  
     &                  nbas, lx,lxx, nx,nxx,
     &                  alat, qlat, vol, ngvecc,
     &          strx,rojp,rojb, sgbb,sgpb, fouvb,   !sgpp,fouvp,
     i          nblochpmx,bas,rmax,
     i    eee, aa,bb,nr,nrx,rkpr,rkmr,rofi,
!         These inputs are to generate sgpp on the fly.
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

      integer(4):: nrx,nr(nbas),ir,ig
      real(8):: eee , int1x(nrx),int2x(nrx),phi(0:lxx),psi(0:lxx)
     &  ,aa(nbas),bb(nbas),rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas)
     &  ,rofi(nrx,nbas)
      real(8), allocatable:: ajr(:,:,:,:), a1(:,:,:)
      logical :: debug=.false.
c---------------------------------------------------------------
      write(6,'(" vcoulq_4: nblochpmx  nbloch ngc=",3i6)') nblochpmx,nbloch,ngc
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

ccccccccccccccccccccccccccc
c      goto 1112
ccccccccccccccccccccccccccc

c <P_G|v|B>
      if(debug) write(6,*)' vcoulq_4: pgvb block 1111'
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

      if(debug) write(6,*)' vcoulq_4: ajr allocate'
C... prepare funciton ajr and a1.
C... ajr:spherical bessel, a1: integral of (sperical bseel)*(rkp rkm)
c------------------
      allocate( ajr(nrx,0:lxx,nbas,ngc), a1(nrx,0:lxx,nbas) )
      if(debug) write(6,*)' vcoulq_4: end ajr allocate'
      do ig1 = 1,ngc
        do ibas= 1,nbas
          if(debug) write(6,"('ccc: ',10i15)")ig1,ibas
          do ir = 1,nr(ibas)
            call bessl(absqg2(ig1)*rofi(ir,ibas)**2,lxx,phi,psi)
            do l  = 0,lx(ibas)

              if(debug.and.ig==162.and.ibas==8) then
                write(6,"('ccc: ',10i15)")ig1,ibas,ir,l
                write(6,*)"ccc:", phi(l)
                write(6,*)"ccc:", rofi(ir,ibas)
              endif

              ajr(ir,l,ibas,ig1) = phi(l)* rofi(ir,ibas) **(l +1 )
              ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
            enddo
          enddo
        enddo
      enddo
c------------------

c <P_G|v|P_G>
      if(debug) write(6,*)' vcoulq_4: pgvpg block'
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

c----------------------
        do ibas=1,nbas
          do l = 0,lx(ibas)
            call intn_smpxxx( rkpr(1,l,ibas), ajr(1,l,ibas,ig1),int1x
     &         ,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas),0)
            call intn_smpxxx( rkmr(1,l,ibas), ajr(1,l,ibas,ig1),int2x
     &         ,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas),0)
            a1(1,         l,ibas) = 0d0
            a1(2:nr(ibas),l,ibas) = 
     &             rkmr(2:nr(ibas),l,ibas) *( int1x(1)-int1x(2:nr(ibas)) )
     &           + rkpr(2:nr(ibas),l,ibas) *  int2x(2:nr(ibas))
          enddo
        enddo
c---------------------

        do ig2 = 1,ig1
          ipl2 = nbloch + ig2
          if(ig1==ig2) vcoul(ipl1,ipl2) = fpivol/(absqg2(ig1) -eee) !eee is negative
          do ibas2= 1, nbas
c... for fouvp and sgpp -------
            call wronkj( absqg2(ig1), absqg2(ig2), rmax(ibas2),lx(ibas2),
     o                 fkk,fkj,fjk,fjj)

            if(eee==0d0) then
              call sigintpp( absqg2(ig1)**.5, absqg2(ig2)**.5, lx(ibas2), rmax(ibas2),
     o             sigx)
            else
              do l = 0,lx(ibas2)
                call gintxx(a1(1,l,ibas2), ajr(1,l,ibas2,ig2)
     &                   ,aa(ibas2),bb(ibas2),NR(ibas2), sigx(l))
              enddo
            endif
            do l = 0,lx(ibas2)
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo

c------------------------------
            do lm2  = 1, (lx(ibas2)+1)**2
              l= ll(lm2)
c...fouvp sgpp-----------
              fouvp_ig1_ig2 = fpi/(absqg2(ig1)-eee) *dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2))
     &      * (-fjj(l)) * pjyl_(lm2,ig2)*phase(ig2,ibas2)
              fouvp_ig2_ig1 = fpi/(absqg2(ig2)-eee) *dconjg(pjyl_(lm2,ig2)*phase(ig2,ibas2))
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
ccccccccccccccccccccccccccccccc
c 1112 continue
ccccccccccccccccccccccccccccccc


c-- Right-Upper part of vcoul.
      if(debug) write(6,*)' vcoulq_4: right-upper'
      do ipl1=1, nbloch+ngc
        do ipl2=1, ipl1-1
          vcoul(ipl2,ipl1) = dconjg(vcoul(ipl1,ipl2))
        enddo
      enddo

cccccccccccccccccccccccccccccc
c test.xxxxxxxxxx
c$$$      do ibl2= 1, nbloch
c$$$        ibas2= ibasbl(ibl2)
c$$$        n2   = nbl (ibl2)
c$$$        l2   = lbl (ibl2)
c$$$        m2   = mbl (ibl2)
c$$$        lm2  = lmbl(ibl2)
c$$$        if(l2==1.and.ibas2>2) then
c$$$          vcoul(nbloch+1:nbloch+ngc, ibl2) = 0d0
c$$$          vcoul(ibl2, nbloch+1:nbloch+ngc) = 0d0
c$$$        endif
c$$$      enddo
ccccccccccccccccccccccccccc

c vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
c     vcoul = 2d0*vcoul !  in Ry unit.
c

c check write
      do ix = 1,nbloch+ngc,20
        write(6,"(' Diagonal Vcoul =',i5,2d18.10)") ix,vcoul(ix,ix)
      enddo
      if( allocated(yl)   ) deallocate(yl)
      if( allocated(cy)   ) deallocate(cy)
      if( allocated(phase)) deallocate(phase)
      if( allocated(pjyl_)) deallocate(pjyl_)
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




c==========================================================================
      subroutine mkjp_4( q,ngc,ngvecc, alat, qlat, lxx,lx,nxx,nx,  
     i                   bas, a,b,rmax,nr,nrx,rprodx,
     i    eee,rofi,rkpr,rkmr,
     o        rojp,sgpb,fouvb)
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
     &   ,r2s,sig,sig1,sig2,sigx(0:lx),sig0(0:lx) ,qg2(3)
      complex(8) :: img =(0d0,1d0),phase
      complex(8),allocatable :: pjyl(:,:)
      real(8),allocatable ::ajr(:,:,:),a1(:,:,:), !rkpr(:,:),rkmr(:,:),
     &  qg(:,:),absqg(:)


      real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),eee
      logical :: debug=.false.
c rkpr(nr,0:lx),rkmr(nr,0:lx),
c-------------------------------------------------
      if(debug) print *,' mkjp_4:'
      nlx = (lx+1)**2
      allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),
     &  qg(3,ngc),absqg(ngc),
     &  pjyl((lx+1)**2,ngc) )

      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
      tpiba = 2*pi/alat
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
!      print *,' mkjp_4: end of sylmnc'
C... q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
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

cc rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
c      rofi(1) = 0d0
c      do ir   = 1, nr
c        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
c      enddo
c      do l = 0,lx
c        rkpr(1:nr,l) = rofi(1:nr)**(l      +1 )
c        rkmr(2:nr,l) = rofi(2:nr)**(-l-1   +1 )
c        rkmr(1,l)    = rkmr(2,l)
c      enddo

c rojp
      if(debug) print *,' mkjp_4: rojp'
      do ig1 = 1,ngc
        call wronkj( absqg(ig1)**2, eee, rmax,lx,
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
      if(eee==0d0) then
c        print *,' mkjp_4: use sigintAn1 eee=0(r0c=infty) mode'
        do ig1 = 1,ngc
          call sigintAn1( absqg(ig1), lx, rofi, nr
     o                  ,a1(1:nr, 0:lx,ig1) )
        enddo
c      else
c We need to impliment a version of sigintAn1 to treat eee/=0 case...
      endif

c-------------------------
c sgpb
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)                      ! r jl        , r B(r)
            if(eee==0d0) then
              call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
ccccccccccccccccccc
            else !for a while, we use this version of sgpb
              call sigint_4(rkpr(1,l),rkmr(1,l), lx,a,b,nr, ajr(1,l,ig1),rprodx(1,n,l)
     &            , rofi, sig)
            endif
            sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
c        write(6,*)
ccccccccccccccccccc
          enddo
        enddo
      enddo
ccccccccccccccccccccccccc
c      stop 'test end======================'
cccccccccccccccccccccccc

c---------------------------------------
c sgpp block------->removed
c---------------------------------------

c Fourier
c fouvb
      if(debug) print *,' mkjp_4: Four'
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
     &    fpi/(absqg(ig1)**2-eee) *dconjg(pjyl(lm,ig1))*radint !eee is supposed to be negative

          enddo
        enddo
      enddo
ccccccccccccccccccccc
c        write(6,*)' fourvb sum=',sum (fouvb)
ccccccccccccccccccccc

c-----------------------------
c fouvp block --->removed
c-----------------------------

      deallocate(ajr,a1,   qg,absqg,   pjyl)
      if (allocated( cy )) deallocate(cy)
      if (allocated( yl )) deallocate(yl)
      end






      real(8) function fac2m(i)
cC A table of (2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      logical,save::  init=.true.
      real(8),save:: fac2mm(0:100)
      if(init) then
        fac2mm(0)=1d0
        do l=1,100
          fac2mm(l)=fac2mm(l-1)*(2*l-1)
        enddo
      endif
      fac2m=fac2mm(i)
      end
c=====================================================================
      subroutine genjh(eee,nr,a,b,lx, nrx,lxx,
     o      rofi,rkpr,rkmr)
C-- Generate radial mesh rofi, spherical bessel, and hankel functions
Cr  rkpr, rkmr are real fucntions --
ci eee=E= -kappa**2 <0 
cr      rkpr = (2l+1)!! * j_l(i sqrt(abs(E)) r) * r / (i sqrt(abs(E)))**l
cr      rkmr = (2l-1)!! * h_l(i sqrt(abs(E)) r) * r * i*(i sqrt(abs(E)))**(l+1)
cr rkpr reduced to be r**l*r      at E \to 0
cr rkmr reduced to be r**(-l-1)*r at E \to 0
C-----------------------------------------------------------
      implicit none
      integer(4):: nr,lx, nrx,lxx,ir,l
      real(8):: a,b,eee,psi(0:lx),phi(0:lx)
      real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),fac2m
      rofi(1)    = 0d0
      do ir      = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      if(eee==0d0) then
        do l = 0,lx
          rkpr(1:nr,l) = rofi(1:nr)**(l +1)
          rkmr(2:nr,l) = rofi(2:nr)**(-l-1 +1)
          rkmr(1,l)    = rkmr(2,l)
        enddo
      else
        do ir  = 1, nr
          call bessl(eee*rofi(ir)**2,lx,phi(0:lx),psi(0:lx))
          do l = 0,lx    !fac2m(l)= (2l-1)!!
c            print *,' phi=',l,phi(l),phi(l)*fac2m(l+1)
c            print *,' psi=',l,psi(l),psi(l)/fac2m(l)
            rkpr(ir,l) = phi(l)* rofi(ir)**(l +1) *fac2m(l+1)
            if(ir/=1) rkmr(ir,l) = psi(l)* rofi(ir) **(-l ) /fac2m(l)
          enddo
        enddo
        rkmr(1,0:lx) = rkmr(2,0:lx)
      endif
      end
c=============================================================
      subroutine mkjb_4( lxx,lx,nxx,nx,  
     i                 a,b,nr,nrx,rprodx,
     i      rofi,rkpr,rkmr,
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
     &  fac, xxx,fpi,pi,sig
      real(8) :: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx)
      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
c      real(8),allocatable :: rkpr(:,:),rkmr(:,:)
c
c      allocate(rkpr(nr,0:lx),rkmr(nr,0:lx))
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

c      rofi(1)    = 0d0
c      do ir      = 1, nr
c        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
c      enddo
c      do l = 0,lx
c        rkpr(1:nr,l) = rofi(1:nr)**(l +1)
c        rkmr(2:nr,l) = rofi(2:nr)**(-l-1) *rofi(2:nr)
c        rkmr(1,l)    = rkmr(2,l)
c      enddo

C... initialize
      rojb=0d0
      sgbb=0d0
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
            call sigint_4(rkpr(1,l),rkmr(1,l),lx,a,b,nr,rprodx(1,n1,l),rprodx(1,n2,l)
     &            , rofi,sig )
            sgbb(n1, n2, l)=sig/(2*l+1)*fpi
          enddo
        enddo
      enddo
c      write(6,*) ' rojbsum=', sum(rojb(:,:)),   sum(abs(rojb(:,:)))
c      write(6,*) ' sgbbsum=', sum(sgbb(:,:,:)), sum(abs(sgbb(:,:,:)))
ccccccccccccccccccccccccccccccccccccc
c      write(6,*)' sigint 1 1 0=',sgbb(1, 1, 0) !/(16d0*datan(1d0))
c      sgbb(1, 1, 0) =0d0
cccccccccccccccccccccccccccccccccccccccc
c      deallocate(rkpr,rkmr)
      end


c-------------------------------------------------------------
      subroutine sigint_4(rkp,rkm,kmx,a,b,nr,phi1,phi2,rofi, sig)
      implicit none
      integer(4) :: nr,kmx,k,ir
      real(8):: a,b, a1(nr),a2(nr),b1(nr),rkp(nr),rkm(nr),
     &    int1x(nr),int2x(nr), phi1(nr), phi2(nr),rofi(nr),sig
      real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
c
      a1(1) = 0d0;  a1(2:nr) = rkp(2:nr)
      a2(1) = 0d0;  a2(2:nr) = rkm(2:nr)
      b1(1:nr) = phi1(1:nr)
      call intn_smpxxx(a1,b1,int1x,a,b,rofi,nr,0)
      call intn_smpxxx(a2,b1,int2x,a,b,rofi,nr,0)
c
      a1(1) = 0d0; a1(2:nr) = 
     &  rkm(2:nr) *( int1x(1)-int1x(2:nr) )+ rkp(2:nr) * int2x(2:nr)
      b1(1:nr) = phi2(1:nr)
      call gintxx(a1,b1,A,B,NR, sig )
      end

c---------------------------------------------------------------
      subroutine intn_smpxxx(g1,g2,int,a,b,rofi,nr,lr0)
c-- intergral of two wave function. used in ppdf
c
c int(r) = \int_(r)^(rmax) u1(r') u2(r') dr'
c
c lr0 dummy index, now not used.
c simpson rule ,and with higher rule for odd devision.
c --------------------------------------------------------------
      IMPLICIT none
      integer nr,ir,lr0
      double precision g1(nr),g2(nr),int(nr),a,b,rofi(nr),w1,w2,w3
     &    ,ooth,foth
      data ooth,foth/0.33333333333333333,1.3333333333333333333/
      data w1,w2,w3/0.41666666666666666,0.6666666666666666666,
     &              -0.083333333333333333/
      if(mod(nr,2).eq.0) 
Cstop2rx 2013.08.09 kino     &  stop ' INTN: nr should be odd for simpson integration rule'
     &  call rx( ' INTN: nr should be odd for simpson integration rule')
c
      int(1)=0.0d0
      DO  10  IR = 3,NR,2
        int(ir)=int(ir-2)
     &             + ooth*G1(IR-2)*G2(IR-2)*( a*(b+rofi(ir-2)) )
     &             + foth*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )
     &             + ooth*G1(IR)*G2(IR)*( a*(b+rofi(ir)) )
   10 CONTINUE

c At the value for odd points, use the same interpolation above
      do 20 ir = 2,nr-1,2
        int(ir)=int(ir-1)
     &             + w1*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )
     &             + w2*G1(IR)  *G2(IR)*  ( a*(b+rofi(ir)  ) )
     &             + w3*G1(IR+1)*G2(IR+1)*( a*(b+rofi(ir+1)) )
   20 continue
      do ir=1,nr
        int(ir)=int(nr)-int(ir)
      enddo
      END

c-------------------------------------------------------------------
      subroutine sigintAn1( absqg, lx, rofi, nr,
     o             a1int)
c a1int(r') = r' * \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} *
c                j_l(absqg r)/absqg**l
      implicit none
      integer(4) :: nr,l,ir,lx
      real(8):: a1int(nr,0:lx), rofi(nr),absqg
      real(8)::  
     &   ak(0:lx) ,aj(0:lx), dk(0:lx), dj(0:lx),
     &   aknr(0:lx),ajnr(0:lx),dknr(0:lx),djnr(0:lx),
     &   phi(0:lx),psi(0:lx)
c---
c      print *,' sigintAn1: absqg=',absqg
      if(absqg<1d-10) then
c      if(absqg<1d-6) then !23jan2004 1d-10 ok?
Cstop2rx 2013.08.09 kino        stop "sigintAn1: absqg=0 is not supported yet. Improve here."
        call rx( "sigintAn1: absqg=0 is not supported yet. Improve here.")
c This part for absqg=0 has not been checked yet!
c       call bessl(0d0,lx,phi,psi)
c        do ir = 1,nr
c        do l  = 0,lx
c          a1int(ir,l) = .5d0* rofi(nr)**2     * rofi(ir)**l     * phi(l)
c     &                +(1d0/(2d0*l+3d0)-.5d0) * rofi(ir)**(l+2) * phi(l)
c       enddo
c        enddo
      else
        call  radkj(absqg**2, rofi(nr),lx,aknr,ajnr,dknr,djnr,0)
        a1int(1,:) = 0d0
        do ir = 2,nr
          call radkj(absqg**2, rofi(ir),lx,ak,aj,dk,dj,0)
          do l  = 0,lx
            a1int(ir,l) = ( (2*l+1)* aj(l) 
     &     -((l+1)* ajnr(l)+ rofi(nr)*djnr(l) )*(rofi(ir)/rofi(nr))**l)
     &    /absqg**2 
     &   *rofi(ir)
          enddo
        enddo
      endif
c      print *,' sigintAn1: end'
      end

c-------------------------------------------------
      subroutine sigintpp( absqg1, absqg2, lx, rmax,
     o             sig)
c sig(l)   =  \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} *
c               j_l(absqg1 r)/absqg1**l
c               j_l(absqg2 r)/absqg2**l
c e1\ne0 e2\ne0
      implicit none
      integer(4) :: l,lx
      real(8)::  rmax,sig(0:lx), absqg1,absqg2, e1,e2,
     &   ak1(0:lx) ,aj1(0:lx), dk1(0:lx), dj1(0:lx),
     &   ak2(0:lx) ,aj2(0:lx), dk2(0:lx), dj2(0:lx),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx)
c---
      e1 = absqg1**2
      e2 = absqg2**2
c
c      print *," sigintpp",e1,e2
c
      call wronkj( e1,e2, rmax,lx,   fkk,fkj,fjk,fjj )
      call  radkj( e1,    rmax,lx,   ak1,aj1,dk1,dj1,0)
      call  radkj( e2,    rmax,lx,   ak2,aj2,dk2,dj2,0)
c
      do l = 0,lx
        sig(l)= ( -l*(l+1)*rmax*aj1(l)*aj2(l) 
     &            + rmax**3 * dj1(l)*dj2(l)
     &            + 0.5d0*rmax**2* (aj1(l)*dj2(l)+aj2(l)*dj1(l))
     &            - fjj(l)*(2*l+1)*(e1+e2)/2d0
     &           ) /(e1*e2)
      enddo
      end

