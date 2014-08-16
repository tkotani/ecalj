      subroutine ppbafp_v2 (ig,ng,isp,nspin,
     i                 il,in,im, mnl,
     d                 nl,nn,nclass,mnlx,
     i                 mdimx,lx,nx,nxx, !Bloch wave    
     i                 cgr,lmxax,   !rotated CG
     i                 ppbrd,           !radial integrals
     o  ppb)
c ---
c calculates for all R (class)
c <Phi(RLn) Phi(RL'n') B(R,i)>
c   n differenciate core phi phidot localOrbital.
c
c B(R,i)   = orthonormal basis, i=1,mdim
c
c---------------------
c iln            = s. indxln.f
c nn,nl          = number of n and l
c
c mdim           = dimension of orthonormal basis
c mnl   =  (m,n,l) composite idex --- numbering.
c in,il,im      = index for n,l,m s. indxlnm.f
c
c ppb            = <Phi(RLn) Phi(RL'n') B(R,i)>
c
      implicit none
      integer(4):: ic,nclass,nl,nn,mnlx,mdimx,nlmx,nlxx,
     &             isp, nxx,lmxax,
     &             i,lb,nb,mb,lmb,i1,ibas,i2,ig,
     &             np,lp,mp,lmp,n,l,m,lm,  nlmxa,nspin,ng
      integer(4):: 
     i          il(mnlx,nclass),in(mnlx,nclass),im(mnlx,nclass),
     i          mnl(nclass),
c
     i    lx(nclass),nx(0: 2*(nl-1),nclass),    j1,j2
c
      real(8):: ppb(mnlx,mnlx,mdimx,nclass)
c     o          pdb(mnlx,mnlx,mdimx,nclass),
c     o          dpb(mnlx,mnlx,mdimx,nclass),
c     o          ddb(mnlx,mnlx,mdimx,nclass)

      real(8) ::      
     &   cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng), 
     &           ppbrd( 0:nl-1,nn, 0:nl-1,nn, 0:2*(nl-1), nxx, 
     &                 nclass*nspin)
     & ,sum1,xxx
      integer(4):: zvztest
ccccccccccccccccccc
c      write(6,*)'    ppbafp: sum(ppbrd)',sum(ppbrd)
c      ppb=0d0;pdb=0d0;dpb=0d0;ddb=0d0
c      ppb=1d99;pdb=1d99;dpb=1d99;ddb=1d99
cccccccccccccccccccc

c loop over class  ! nclass=nbas
CCC!$OMP parallel do  private ( ic,ibas, i, ib,nb,mb, lmb, i2,np,kp,mp, lmp)
CCC!$OMP& private( i1,n,l,m,lm, lp )
      do ic  = 1, nclass
        ibas = ic
c        do ip1=1,2
c        do ip2=1,2
c          ioff(ip1,ip2) =  ip1+ 2*(ip2-1)+ 4*(ic-1) + 4*nclass*(isp-1)
c        enddo
c        enddo
c     print *, ' ic lx=',ic, lx(ibas)

!! loop over product basis
        i = 0 !i = product basis index.
        do lb  = 0, lx (ibas)
          do nb  = 1, nx (lb,ibas)
            do mb  = -lb, lb
              i    = i+1  !The number of product basis is  =(i at the end of loop).
              lmb  = lb*lb + lb + mb + 1

cccccccccccccccccc
c              write(6,"(' === product basis index: i ibas l n m',10i4)")i,ic,lb,nb,mb
cccccccccccccccccc

c loop over LMTO basis
              do  i2 = 1,mnl(ic)
                np   = in(i2,ic)
                lp   = il(i2,ic)
                mp   = im(i2,ic)
                lmp  = lp*lp + lp + mp + 1

                do  i1 = 1,mnl(ic)
                  n    = in(i1,ic)
                  l    = il(i1,ic)
                  m    = im(i1,ic)
                  lm   = l*l + l + m + 1
cccccccccccccccccccccccccccccccc
c        do j1=1,2;doj2=1,2
c        if(abs(ppbrd(l,n, lp,np, lb,nb, ioff(j1,j2))) >1d10) then
c          print *, abs(ppbrd( l,n, lp,np, lb,nb, ioff(j1,j2)))
c          print *, l,n, lp,np, lb,nb, ioff(j1,j2)
c          stop
c        endif
c       enddo; enddo
ccccccccccccccccccccccccccccccccccccc
                  ppb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
     &                  * ppbrd( l,n, lp,np, lb,nb, isp+nspin*(ic-1))

cccccccccccccccccccccccccccccccccccccccccccccccccc
                  if(zvztest()==2) then
                    if(abs(ppb(i1,i2,i,ic))/=0d0  )then
                      write(6,"('i1 i2 i ic=',4i4,' ppb=',f10.6)") i1,i2,i,ic, ppb(i1,i2,i,ic)
                    endif
                  endif
cccccccccccccccccccccccccccccccccccccccccccccccccc

c        pdb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(1,2))
c        dpb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,1))
c        ddb(i1,i2,i,ic) = cgr(lm,lmp, lmb,ig)
c     &                  * ppbrd( l,n, lp,np, lb,nb, ioff(2,2))
ccccccccccccccccccc
c         xxx = ppb(5+1,1,i,1)
c         ppb(i1,i2,i,ic) = 0d0
c         ppb(5+1,1,i,1) = xxx
c         if(i1==6.and.i2==1.and.ic==1.and.ig==1) then
c          write(6,'(4i4,2d15.6,7i2)')lm,lmp,i,lmb, cgr(lm,lmp, lmb,ig),
c     &       ppbrd( l,n, lp,np, lb,nb, ioff(1,1))
c     &       ,l,n, lp,np, lb,nb, ioff(1,1)
c         endif
ccccccccccccccccccc
c      if(ppb(i1,i2,i,ic)/=0d0)
c     & write(6,'(4i4,d15.6)'),i1,i2,i,ic,ppb(i1,i2,i,ic)
ccccccccccccccccccc
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
CCC!$OMP end parallel do
cccccccccccccccccccccccccccccccccccccccccccc
c      pdb =0d0
c     dpb =0d0
c     ddb =0d0
ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
c     print *,' ppbafp x sum ppb',sum(ppb)
c     stop
cccccccccccccccccccccccccccccccccccccccccccc
      end

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
      integer:: ll1,ll2
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
      ll2 = (2*lmxax+1)**2
      ll1 = (lmxax+1)**2
      allocate( cgn(ll2,ll1,ll1) )
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
Cstop2rx 2013.08.09 kino        stop       'rotcg: cannot handle lmxcg'
        call rx( 'rotcg: cannot handle lmxcg')
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

c$$$c-------------------------------------------------------
c$$$      subroutine drvmelp( iq, ntp0, kp,nt0, k,
c$$$     i        ngc,ngcmx,ngpn, ngpmx,ngvecpB,ngveccBr,geigB,nband,itp0,
c$$$     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c$$$     i       rmel, cmel, nbloch,noccx,nctot,
c$$$     o       rmelt,cmelt)
c$$$c ! q1=q    ---> iq 1:ntp0 q-point
c$$$c ! q2=q-rk ---> kp 1:nt0  occupied
c$$$c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c$$$c  zmelpl
c$$$c-----------------------------------------------------
c$$$      implicit none
c$$$      integer(4):: iq,kp,k,ngp1, ngp2, ngpmx,nqbz,nqibz,
c$$$     &      ngcmx ,nctot,nband,
c$$$     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
c$$$     &             ngvecpB  (3,ngpmx,nqbz),
c$$$     &             ngveccBr (3,ngcmx,nqibz), itp0(ntp0),
c$$$     &     nadd(3),ngpn(nqbz)
c$$$      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
c$$$      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
c$$$     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
c$$$     &         ,rmel(nbloch,noccx,ntp0)
c$$$     &         ,cmel(nbloch,noccx,ntp0)
c$$$     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c$$$     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c$$$     & ,pi=3.1415926535897932D0
c$$$      complex(8) :: geigB  (ngpmx,nband,nqbz)
c$$$c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
c$$$      ngp1 = ngpn(iq)  ! q        in FBZ
c$$$      ngp2 = ngpn(kp)  ! q-rk     in FBZ
c$$$c     ngc  = ngcni(k)  ! k-points in IBZ
c$$$      allocate( geigq(ngp1,ntp0), zmelpl(ngc,nt0,ntp0) )
c$$$      do itx = 1,ntp0
c$$$        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c$$$      enddo
c$$$c     print *,' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c$$$      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp)
c$$$      ! rk          -q          +(q-rk) is not zero.
c$$$c     print *,' drvmelp qdiff',qdiff
c$$$      add  = matmul(qbasinv,qdiff)
c$$$      nadd = idint( add + dsign(.5d0,add))
c$$$c      write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)')
c$$$c     &  add,nadd,abs(add-nadd)
c$$$      if(sum(abs(add-nadd))>1d-7) stop "drvmlp: abs(add-nadd))>1d-10"
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      if(.false.) then !old version
c$$$c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c$$$c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c$$$c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c$$$c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c$$$c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c$$$c     i       shtv, matmul(symope,qibz(1:3,k)), qbas,
c$$$c     o       zmelpl)
c$$$c      else             !new version from Sep-10 2000
c$$$      call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c$$$     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c$$$     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c$$$     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c$$$     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c$$$     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas,
c$$$     o       zmelpl)
c$$$c       endif
c$$$      deallocate(geigq)
c$$$      rmelt=0d0
c$$$      cmelt=0d0
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      zmelpl=0d0
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c$$$     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c$$$      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c$$$     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      print *,' !!!!!! NO PLANE WAVE !!!'
c$$$c      return
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c extensiton to nbloch+ngc
c$$$      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c$$$     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
c$$$      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c$$$     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
c$$$      deallocate(zmelpl)
c$$$
c$$$      end

c------------------------------------------------------
      subroutine drvmelp2( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx_in,nband,itp0, 
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       rmelt,cmelt)
      use m_readqg, only : readqg
      use m_readeigen,only : readgeig
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
c-----------------------------------------------------
      implicit none
      real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4):: ngp1, ngp2, ngpmx_in,nqbz,nqibz,  !iq,kp,k,
     &      ngcmx ,nctot,nband,
     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
     &             ngvecpB1(3,ngpmx_in),  
     &             ngvecpB2(3,ngpmx_in),  
     &             ngveccBr(3,ngcmx), itp0(ntp0),
     &     nadd(3),isp  !,ngpn(nqbz)
      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0) 
     &         ,cmel(nbloch,noccx,ntp0) 
     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     & ,pi=3.1415926535897932D0
!      complex(8) :: geigB  (ngpmx_in,nband,nqbz)
      complex(8):: geig1(ngpmx_in,nband),geig2(ngpmx_in,nband)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
      logical:: debug=.false.
      if(verbose()>=90) debug=.true.
!TIME0
c
c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) print *,' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,   ginv, qt,   ngp1, ngvecpB1) !qt is used q in practice.
      call readqg('QGpsi', q_rk,ginv, q_rkt,ngp2, ngvecpB2)
c      write(*,"('ddd1: ',i4,3f10.6,2x,3f10.6,3x,i15)")ngp1,q,qt      ,sum(abs(ngvecpB1(:,1:ngp1)))
c      write(*,"('ddd1: ',i4,3f10.6,2x,3f10.6,3x,i15)")ngp2,q_rk,q_rkt,sum(abs(ngvecpB2(:,1:ngp2)))
      if(debug) print *,'drvmelp2: goto readgeig'
      call readgeig(q,    ngpmx_in, isp, qu1, geig1)
      call readgeig(q_rk, ngpmx_in, isp, qu2, geig2)
!TIME1 "after readgeig q q_rk"
      if(debug) print *,'drvmelp2: endof readgeig'
c      qu1=qt
c      qu2=q_rkt
c      geig1=0d0
c      geig2=0d0
c      write(*,"('ddd2: ',i4,3f10.6,2x,3f10.6,3x,d23.16,i15)")ngp1,q,qt      ,sum(abs(
c     &       geig1(1:ngp1,itp0(1:ntp0)))), ntp0  ! q1=q    ---> iq ngp1 1:ntp0 q-point
c      write(*,"('ddd2: ',i4,3f10.6,2x,3f10.6,3x,d23.16,i15)")ngp2,q_rk,q_rkt,sum(abs(
c     &       geig2(1:ngp2,1:nt0))), nt0     ! q2=q-rk ---> kp ngp2 1:nt0  occupied
cc      allocate( zmelpl(ngc,nt0,ntp0) )
c      print *,'test 888888888888888888888',ngc,nt0,ntp0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!TIME0
Cstop2rx 2013.08.09 kino      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(qt-qu1))>1d-10) call rx( 'drvmelp2;qu1/=qu1x')
Cstop2rx 2013.08.09 kino      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
      if(sum(abs(q_rkt-qu2))>1d-10) call rx( 'drvmelp2;qu2/=qu2x')
      if(debug) print *,' drvmelp2:end of read geig '

c     ngc  = ngcni(k)  ! k-points in IBZ
c      allocate( geigq(ngp1,ntp0) )
c      do itx = 1,ntp0
c        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c      enddo

c     print *,' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp)
      qdiff = matmul(symope,qik)  - qt + q_rkt
      ! rk    -q  +(q-rk) is not zero.
!! <rk q-rk |q>
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) then
        print *,' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
     &  add,nadd,abs(add-nadd)
Cstop2rx 2013.08.09 kino        stop "drvmlp: abs(add-nadd))>1d-10"
        call rx( "drvmlp: abs(add-nadd))>1d-10")
      endif

      if(debug) print *,' goto alloc zmelpl'
      allocate( zmelpl(ngc,nt0,ntp0) )
cccccccccccccccccccccccccccccccccccccccccccccccc
c      if(.false.) then !old version
c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)), qbas,
c     o       zmelpl)
c      else             !new version from Sep-10 2000
c       call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas,
c     o       zmelpl)
c       endif
      if(debug) print *, 'drvmelp2:goto melpln2'
ccccccccccccccccccc
c      print *, 'ssss: sum check geig1=',sum(abs(geig1(1:ngp1,itp0(1:ntp0))))
c      print *, 'ssss: sum check geig2=',ngp2,nt0,sum(abs(geig2(1:ngp2,1:nt0)))
ccccccccccccccccccccccccccccccccccc
!TIME1 "before drv melpln2t"
!TIME0 
      call melpln2t(ngp1, ngvecpB1
     &           ,  ngp2, ngvecpB2
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,! q1=(shifted q) ->iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,1:nt0), nt0, ! q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied 
     i       shtv, matmul(symope,qik),qik, symope, qbas,
     i       qt, !qt oct2013 means q1
     o       zmelpl)
!TIME1 "after drv melpln2t"
!TIME0 
      if(debug) print *,'end of melpln2t'

 8888 continue
c     deallocate(geig1,geig2)
      rmelt=0d0
      cmelt=0d0
      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =  
     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
cc extensiton to nbloch+ngc
      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
      deallocate(zmelpl)
      if(debug) print *,'end of drvmelp2'
!TIME1 "end drvmelpln2t"
      end



c------------------------------------------------------
      subroutine rotgvec(symops, nqibz, 
     i  ngcmx,ngcn, qbas,ngveccB, 
     o  ngveccBr)
c Rotated ngveccB
c G' = R(G), where R denotes rotation and
cr We determine ngveccBr so that
cr      matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))
cr      =  matmul(qbas, ngveccBr(1:3, igc, iq)).
cr See the variable sumchck.
      implicit none
      integer(4) :: irot, ngpmx, nqbz,nqibz, ngcmx, iqbz,
     i            ngcn(nqibz), idxk,
     &    ngveccB(3,ngcmx,nqibz), ngveccBr(3,ngcmx,nqibz),iq,igc,kp,verbose
      real(8)    :: symops(3,3), qq(3), qrot(3),qbas(3,3)
     &             ,qbasinv(3,3),
     &     qfbz(3),rotnvg(3,3),vec(3),xsign(1:3),det,sumchk !,sumchk2
      logical:: debug=.false.
c
      if(verbose()>=100) debug=.true.
      if(debug) print *,' rotgvec: '
      call dinv33(qbas,0,qbasinv,det)
!      symops = grp(:,:,irot)
      if(debug) write(6,'(3f12.6)') symops(1,1:3)
      if(debug) write(6,'(3f12.6)') symops(2,1:3)
      if(debug) write(6,'(3f12.6)') symops(3,1:3)

      rotnvg = matmul(qbasinv,matmul(symops,qbas))
      if(debug) print *, 'sum rotnvg=',sum(rotnvg),ngcmx,ngcn(1:nqibz)
      sumchk  =0d0
c      sumchk2 =0d0
      do iq  = 1, nqibz
        if(debug) print *,' iq=',iq, ' sum ngveccB=', sum(abs(ngveccB(1:3,1:ngcn(iq), iq)))
        do igc = 1, ngcn(iq)
          vec  = matmul( rotnvg, ngveccB(1:3,igc, iq))
! vec should be the almost integer and ngveccBr = vec.
! But we need this procedure in order to get correct integer value.
          ngveccBr(1:3, igc, iq) = idint( vec + dsign(.5d0,vec))
        enddo
        do igc= 1, ngcn(iq)
          sumchk = sumchk + 
     &    sum(abs( 
     &       matmul(qbas, ngveccBr(1:3, igc, iq))
     &     - matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))))
        enddo
cccccccccccccccc
c        do igc= 1, ngcn(iq)
c          write(6,*) " igc=",igc
c          write(6,'(3f13.5)')
c     &       matmul(qbas, ngveccBr(1:3, igc, iq))
c          write(6,'(3f13.5)')
c     &      matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))
c          sumchk2 = sumchk2 +
c     &    sum(abs(
c     &     qrot + matmul(qbas, ngveccBr(1:3, igc, iq)) ))
c        enddo
cccccccccccccccc
ccccccccccccccc
        if(debug) print *," rotgvec: nmin nmax="
     & ,minval(ngveccBr(1:3, 1:ngcn(iq), iq))
     & ,maxval(ngveccBr(1:3, 1:ngcn(iq), iq))
ccccccccccccccc
      enddo
      if(abs(sumchk)/nqibz/minval(ngcn)>1d-4) then
        print *," rotgvec: sum chk error sumchk=",sumchk
Cstop2rx 2013.08.09 kino        stop "rotgvec: sum chk error >1d-4"
        call rx( "rotgvec: sum chk error >1d-4")
      endif
c      print *,' rotgvec: sumcheck',sumchk !,sumchk2
c      stop "rotgvec: test end"
      end

c
c-----------------------------------------------------
      subroutine q0iwgt(symops,ngrp,q0i,nq0i, 
     o                  wgt0)
c Get weight for each k-points near 0.
c wgt0(irreducible-k, irotation)
      implicit none
      integer(4) :: ixx,ix,i, ngrp,ig,nq0i
      real(8)     :: q0(3,6),q0i(3,6),symops(3,3,ngrp),sym(3,3),
     &   qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp)
c      print *,' q0iwgt:'
      wgt0 = 0d0
      ixx = 0
      do i = 1,nq0i
        qt = q0i(:,i)
c equivalence check
        do ig = 1,ngrp
          sym = symops(:,:,ig)
c          write(66,'(2i3,3f12.5)')i,ig,matmul(sym,qt)
ccccccccccccccccccccccccccccccccccccccc
c          do ix = 1,ixx
c            if(sum(abs(q0in(:,ix)-matmul(sym,qt)))<1d-10) then
c              goto 1111
c            endif
c          enddo
ccccccccccccccccccccccccccccccccccccccc
          ixx = ixx+1
          q0in (:,ixx) = matmul(sym,qt)
          wgt0 (i,ig)  = 1d0
c          write(6,'(2i3,3f12.5)')i,ig,q0in (:,ixx)
 1111     continue
        enddo
      enddo
      wgt0  = wgt0 /ixx
ccccccccccccccccccccc
c     do ig=1,ngrp
c     write(6,*)' sym ig  =',ig
c     write(6,'(3f12.6)') symops(1,1:3,ig)
c     write(6,'(3f12.6)') symops(2,1:3,ig)
c     write(6,'(3f12.6)') symops(3,1:3,ig)
c      enddo
ccccccccccccccccccccc
      end

c--------------------------------------------------
      subroutine melpln
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc,ngvecc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qbas, 
     o       zmelp)
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
co rmelp, cmelp
c-----------------------------------------------------
      implicit none
      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
     &     itp,it, ngc,ngvecc(3,ngc),nnum,inum,ig1,ig2,igc,ntp0,nt0,
     &     nadd(3)
      complex(8):: geigq1(ngp1,ntp0),
     &             geigq2(ngp2,nt0), zmelp(ngc,nt0,ntp0),phase(ngc),
     &             img=(0d0,1d0)
      real(8)::  q(3),qg(3),qbas(3,3),shtv(3),qdummy(3)
     & ,pi=3.1415926535897932D0
c      inum =0
c
c      print *,' melpln: sht nadd=',shtv,nadd
      do igc = 1,ngc
        qg(1:3) = 2*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
        phase(igc) = exp(img*sum(qg*shtv))
      enddo

      zmelp=0d0
      do ig1=1,ngp1
        do ig2=1,ngp2
          do igc=1,ngc
            if(sum(abs( ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2)
     &             -ngvecc(1:3,igc) -nadd))  ==0) then
c          inum = inum+1
! igdelta(3,inum) means the delta function delta_G1,G2+G.
c          igdelta(1:3,inum) = (/ig1,ig2,igc/)
              do itp=1,ntp0
                do it =1,nt0
                  zmelp(igc,it,itp) = zmelp(igc,it,itp)
     &        + dconjg(geigq1(ig1,itp))*geigq2(ig2,it) * phase(igc)

cccccccccccccccccccccccccccccccccccccccccccc
c      if(sum(abs( ngvecp1(1:3,ig1)))<3.and.
c     &   igc==4.and.it==32.and.itp==1) then
c         write(6,'("  g1 g2 zmeladd=",3i3,1x,3i3,1x,3i3,12d12.3)')
c     &   ngvecp1(1:3,ig1),ngvecp2(1:3,ig2),ngvecc(1:3,igc),
c     &   dconjg(geigq1(ig1,itp))*geigq2(ig2,it) * phase(igc),
c     &    geigq1(ig1,itp),geigq2(ig2,it) ,phase(igc)
c      endif
cccccccccccccccccccccccccccccccccccccccccc
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
c      nnum=inum
c      do inum=1,nnum
c        (/ig1,ig2,igc/) = igdelta(1:3,inum)
c      enddo
ccccccccccccccccccccccccccccccccccccc
c      igc=4
c      it = 32
c      itp=1
c      write(6,'("  igc it itp zmel=",3i4,2d15.6)')
c     &  igc,it,itp,zmelp(igc,it,itp)
c      stop 'TEST end melpln:'
ccccccccccccccccccccccccccccccccccccc
      end

c--------------------------------------------------------------------
      subroutine melpln2t
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qi, symope, qbas, !
     i       qt, !qt oct2013 for G1
     o       zmelp)
!! == Mattrix elements <Plane psi |psi> from interstitial plane wave. ==
!! zmelp(igc(qi),it(q2),itp(q1)) = <itp(for q1+G1)| it(for q2+G2) igc> matrix element
!! NOTE: shtv = g(delta_{g^-1})
!! ----------------------------------------------------------
c      use m_read_ppovl,only: getppx2,getppx,
c     &  ppx,ngc2,ngvecc2, nvggg,nvgcgp2,ngvecc, nggg,ngcgp2,ngcread,ggg,ppovlinv
      use m_mpi, only:mpi__rank
      use m_read_ppovl,only: getppx2,
     &  nvggg,nvgcgp2,ngvecc, nggg,ngcgp,ngcread, ggg,ppovlinv,
     &  ngc2,ngvecc2
      implicit none
      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
     &     itp,it, ngc,nnum,inum,ig1,ig2,igc,ntp0,nt0,
     &     nadd(3)
      complex(8):: geigq1(ngp1,ntp0),phas1,phas2,phasqadd,
     &          geigq2(ngp2,nt0),ge2,zmelp(ngc,nt0,ntp0),  !it was zmelp(ngc,nt0*ntp0)
     &             img=(0d0,1d0) ,phas         ,phase2
      real(8)::  q(3),qi(3),qg(3),qbas(3,3),shtv(3)
     & ,pi=3.1415926535897932D0, symope(3,3)   ,qx(3)
      integer(4),allocatable :: ipp(:,:,:)
      complex(8),allocatable :: gg(:,:),gg_(:,:,:),xxx(:,:)
      integer(4) :: iqi,igc2,ngcx,nn(1:3),nn1(3),
     &       nxmin,nymin,nzmin,nxmax,nymax,nzmax,ix,iy
      integer(4)::verbose,nomatm
      logical:: debug=.false.
      real(8):: qu1(3),qu2(3)
      complex(8),allocatable:: geigq1_(:,:)
      integer::igcgp2,iggg,igp1,igp2
      integer:: nxminggg,nxmaxggg,nyminggg,nymaxggg,nzminggg,nzmaxggg
      integer:: nxmingcgp2,nxmaxgcgp2,nymingcgp2,nymaxgcgp2,nzmingcgp2,nzmaxgcgp2
      complex(8):: phas1a
      integer,allocatable:: igggi(:,:,:),igcgp2i(:,:,:),ngveccR(:,:)
      complex(8),allocatable::ggitp(:,:),gp2phas2(:),phase(:)
      integer:: ngcgp2,ngcs(1)
      real(8):: qt(3)
      complex(8)::ppp,zdotc
      complex(8),allocatable:: mmm(:,:),zmelp0(:,:,:),ggitp_(:,:)
#ifdef USE_GEMM_FOR_SUM
      complex(8),allocatable:: z2(:,:)
#endif
c------------------------------
cc getppx allocate and return ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2).
cc      call getppx(qbas,  ngcx,qi,symope) ! rppovl.f

!TIME0
!! allocate and return foobar decleared in use m_read_ppovl,only: foobar
c      if(verbose()>=100) 
cccccccccccccccccccccc
c      nadd=0
cccccccccccccccccccccc
      if(debug) print *,' melpln2t: nadd=',nadd
      call getppx2(qbas,qi) ! rppovl.f
      if(verbose()>=100) print *,' melpln2t: end of getppx2'
      if(ngc/=ngcread) then
         print *,qi,ngcread,ngc
         call rx( 'melpln2: ngc/= ngcx by getppx:PPOVLG')
      endif   
      ngcgp2=ngcgp
!! --- new version --- nov2013takao
!! zmelp(igc,it,itp) = <itp(for G1)|it(for G2) igc> matrix element.
!!   zmelp0(igc,it,itp) = <G1|G2 Gc'> geig^*(G1,itp) geig(G2,it) 
!!   zmelp(igc,it,itp) =   = zmelp0(Gc',it,itp) <Gc'|Gc>^-1
!!   (<Gc'|Gc>^-1 is dconjg(ppovlinv)
!! New ggg matrix <Gc |G1 G2> is introduced.
!!
!!    <Gc G2|G1> is equivalent to <-Gc+G1-G2>; described by ggg
!! Readin input
!!    ggg(1:nggg) = <Gc+G2-G1>
!!    nvggg(3,1:nggg)   for Gc+G2-G1
!!    nvgcgp2(3,ngcgp2) for Gc+G2
!!    ppovlinv(ngc,ngc) <Gc|Gc> matrix
!!
!! Old version used ppx= <Gc|Gc>^-1 <Gc|Gc2> (Gc2 span G1-G2 ) in rdata4gw_v2
!! -------------------------------------
!!  Make igggi inversion table
      nxminggg =minval(nvggg(1,1:nggg)) ;nxmaxggg=maxval(nvggg(1,1:nggg))
      nyminggg =minval(nvggg(2,1:nggg)) ;nymaxggg=maxval(nvggg(2,1:nggg))
      nzminggg =minval(nvggg(3,1:nggg)) ;nzmaxggg=maxval(nvggg(3,1:nggg))
      if(verbose()>100) write(6,"(' nggg nvggg range: ',i6, 6i4)") 
     & nggg, nxminggg,nxmaxggg,nyminggg,nymaxggg,nzminggg,nzmaxggg
      allocate( igggi(nxminggg:nxmaxggg,nyminggg:nymaxggg,nzminggg:nzmaxggg) )
      igggi = -100000
      if(verbose()>=100) print *,' melpln2t: goto do igggi'
      do iggg =1,nggg
        igggi(nvggg(1,iggg),nvggg(2,iggg),nvggg(3,iggg)) = iggg
      enddo
!! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
!! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
      allocate(ngveccR(1:3,1:ngc))
      ngcs(1) = ngc
      if(verbose()>=100) print *,' melpln2t: goto rotation of ngvecc xxxx'
      call rotgvec(symope, 1, ngc, ngcs, qbas, ngvecc, 
     o    ngveccR) 
!TIME1 "end of rotgvec"
!TIME0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! simple case commented out. This worked (but slow) nov2013takao.
c$$$!!      print *,'--- goto simple case for debug---'
c$$$      print *,'fffffff: ngc ngp1 ngp2=',ngc,ngp1,ngp2
c$$$      zmelp0 = 0d0
c$$$      do igc=1,ngc
c$$$        phas =exp( img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shtv) ) !conjg
c$$$      do igp1=1,ngp1     !for ngp1
c$$$      do igp2=1,ngp2
c$$$c        nn = ngvecp1(:,igp1)-ngvecp2(:,igp2)-ngveccR(:,igc) - nadd 
c$$$        nn = -ngvecp1(:,igp1) + ngvecp2(:,igp2) + ngveccR(:,igc) + nadd 
c$$$        iggg = igggi(nn(1),nn(2),nn(3))
c$$$        do it=1,nt0
c$$$        do itp=1,ntp0
c$$$          zmelp0(igc,it,itp) = zmelp0(igc,it,itp)
c$$$     &      +  phas*ggg(iggg)* dconjg(geigq1(igp1,itp))* geigq2(igp2,it)
c$$$        enddo
c$$$        enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
c$$$      return
ccccccccccccccccccccccccccccccccccccccccccccccccc
      if(verbose()>=100) print *,' melpln2t: godo igp1 loop'
      allocate(ggitp(ntp0,ngcgp2))
      ggitp = 0d0
      do igp1  = 1,ngp1   !for ngp1
      do igcgp2= 1,ngcgp2 !for ngc+ngp2 
         nn = ngvecp1(:,igp1)- nvgcgp2(:,igcgp2) - nadd 
         ! G1 -(Gc+G2) - Gadd !  -Gadd= -rk + qt -q_rk
c         print *,' nn=',nn
         if(nn(1)<nxminggg .or. nxmaxggg<nn(1) ) cycle
         if(nn(2)<nyminggg .or. nymaxggg<nn(2) ) cycle
         if(nn(3)<nzminggg .or. nzmaxggg<nn(3) ) cycle
         iggg = igggi(nn(1),nn(2),nn(3))
         if(iggg<0) cycle
c         print *,'iggg=',iggg,nggg
!! ggg(iggg) = <qt+G1 -(rk+Gc) -(q_rk+G2) >, where 
!!       gggx(1:ntp0,igcg2) = gggx(1:ntp0,igcg2)+ ggg(iggg)*geigq1_(1:ntp0,igp1)
         call zaxpy(ntp0,ggg(iggg), geigq1(igp1,1),ngp1, ggitp(1,igcgp2),1)
      enddo
      enddo
!TIME1 "end of ggitp loop"
!TIME0
      if(verbose()>=100) print *,' melpln2t: 11111111'
      deallocate(igggi)
      ggitp = dconjg(ggitp)
!! igcp2 inverse table
      nxmingcgp2 =minval(nvgcgp2(1,1:ngcgp2)) ;nxmaxgcgp2=maxval(nvgcgp2(1,1:ngcgp2))
      nymingcgp2 =minval(nvgcgp2(2,1:ngcgp2)) ;nymaxgcgp2=maxval(nvgcgp2(2,1:ngcgp2))
      nzmingcgp2 =minval(nvgcgp2(3,1:ngcgp2)) ;nzmaxgcgp2=maxval(nvgcgp2(3,1:ngcgp2))
      allocate(igcgp2i(nxmingcgp2:nxmaxgcgp2,nymingcgp2:nymaxgcgp2,nzmingcgp2:nzmaxgcgp2))
      igcgp2i = -100000
      do igcgp2 =1,ngcgp2
        igcgp2i(nvgcgp2(1,igcgp2),nvgcgp2(2,igcgp2),nvgcgp2(3,igcgp2))=igcgp2
      enddo
      if(verbose()>=100) write(6,"(' nggg nvggg range: ',i6, 6i4)") 
     &  ngcgp2, nxmingcgp2,nxmaxgcgp2,nymingcgp2,nymaxgcgp2,nzmingcgp2,nzmaxgcgp2
!! zmelp <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2))
!! note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
!! mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
      allocate(gp2phas2(nt0),phase(ngc))
      do igc=1,ngc
        phase(igc)=exp( img*2d0*pi*sum((q+matmul(qbas,ngveccR(:,igc)))*shtv) )
      enddo  
!! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
      allocate(ggitp_(ngp2,ngcgp2))
      allocate(zmelp0(ngc,nt0,ntp0))
      zmelp0=0d0
      if(verbose()>=100) print *,' melpln2t: 33333333'
!TIME1 "before zmelp0 loop"
!TIME0
!! If range of nvgcgp2 is not large enough, igcgp2 can be -100000.
#ifdef USE_GEMM_FOR_SUM
      write(*,*)'info: USE GEMM FOR SUM (zmelp0=geigq2*ggitp_) in ppbafp.fal.F'
      allocate(z2(ngc,nt0))
      do itp= 1,ntp0
      do igc=1,ngc
         do igp2=1,ngp2
            nn = ngveccR(:,igc) + ngvecp2(:,igp2)
            igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
            ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
         enddo
      enddo
         call zgemm('T','N',ngc,nt0,ngp2,(1.0d0,0.0d0),ggitp_,ngp2,geigq2,ngp2,(0.0d0,0.0d0), 
     &        z2,ngc)
         zmelp0(:,:,itp)= z2
      enddo
      deallocate(z2)
#else
      do igc=1,ngc
      do itp= 1,ntp0
         do igp2=1,ngp2
            nn = ngveccR(:,igc) + ngvecp2(:,igp2) 
            igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
            ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
         enddo
         if(verbose()>=100) print *,' melpln2t: 33333333 22222222222'
!     ! This do block is time-consuming when we use large QpGcut_psi _cou.
!     ! But this version (nov2013) is not so bad, maybe.
         do it=1,nt0
            zmelp0(igc,it,itp) = zmelp0(igc,it,itp) + sum( geigq2(:,it)*ggitp_(:,igc) )
         enddo
      enddo
      enddo
#endif
!TIME1 "end of zmelp0 loop"
!TIME0
c$$$      do itp=1,ntp0
c$$$         do igc =1,ngc
c$$$            do igp2=1,ngp2
c$$$               nn = ngveccR(:,igc) + ngvecp2(:,igp2) 
c$$$               igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
c$$$               ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
c$$$            enddo
c$$$         enddo
c$$$         if(verbose()>=100) print *,' melpln2t: 33333333 22222222222'
c$$$!     ! This do block is time-consuming when we use large QpGcut_psi _cou.
c$$$!     ! But this version (nov2013) is not so bad, maybe.
c$$$         do igc=1,ngc
c$$$            do it=1,nt0
c$$$               zmelp0(igc,it,itp) = zmelp0(igc,it,itp) + sum( geigq2(:,it)*ggitp_(:,igc) )
c$$$            enddo
c$$$         enddo
c$$$      enddo
      if(verbose()>=100) print *,' melpln2t: 4444444444'
      deallocate(ngveccR,phase,ggitp_)
c      print *,'pppppppp',shape(ppovlinv),sum(abs(ppovlinv))
c      print *,'pppppppp',shape(zmelp0),sum(abs(zmelp0))
c      print *,'pppppppp',ngc,ntp0*nt0
c      print *,'pppppppp',shape(zmelp)
c      do ix=1,nt0
c         print *,'ix=',ix, sum(abs(matmul(dconjg(ppovlinv),zmelp0(:,:,ix))))
c         zmelp(:,:,ix)=matmul(dconjg(ppovlinv),zmelp0(:,:,ix))
c      enddo
      call matm(dconjg(ppovlinv),zmelp0,zmelp,ngc,ngc,ntp0*nt0)
c      print *,'pppppppp end of matm' 
      deallocate(igcgp2i,ggitp,gp2phas2,zmelp0)
      if(verbose()>=100) print *,' melpln2t: end'
!TIME1 "end of melpln2t"
      end subroutine melpln2t 


c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$!! old version --------------------------------------------------------
c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c$$$      if(sum(abs(qu1-(/9.6225044864937631E-002,   9.6225044864937631E-002, 0.10000000000000001/)))<1d-8) then
c$$$c$$$         debug=.true.
c$$$c$$$      endif
c$$$c$$$      if(debug) print *,'melpln2: qu1',qu1
c$$$c$$$      if(debug) print *,'melpln2: qu2',qu2
c$$$c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      do ig1=1,ngp1
c$$$        do itp=1,ntp0
c$$$          geigq1_(itp,ig1)=geigq1(ig1,itp)
c$$$        enddo
c$$$      enddo
c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c$$$      if(debug) then
c$$$c$$$       print *,' melpln2 old sumcheck ngc2 ngvec ppx ='
c$$$c$$$     &  ,ngc2,sum(abs(ngvecc2(1:3,1:ngc2))),sum(abs(ppx(1:ngc,1:ngc2)))
c$$$c$$$       print *,'melpln2: 1111111111'
c$$$c$$$       endif
c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      nxmin =minval(ngvecc2(1,1:ngc2)) ;nxmax=maxval(ngvecc2(1,1:ngc2))
c$$$      nymin =minval(ngvecc2(2,1:ngc2)) ;nymax=maxval(ngvecc2(2,1:ngc2))
c$$$      nzmin =minval(ngvecc2(3,1:ngc2)) ;nzmax=maxval(ngvecc2(3,1:ngc2))
c$$$c$$$ccccccccccccccccccccccccccccccccccccccccccc
c$$$c$$$      if(debug) then
c$$$c$$$      write(6,"(' iqi ngc ngc2=',3i8)") iqi,ngc,ngc2
c$$$c$$$c      print *,' scheck2 ppx=',sum(ppovlinv(1:ngcx,1:ngcx))
c$$$c$$$c      print *,' scheck2 ppovls   =',sum(ppovl   (1:ngcx,1:ngc2))
c$$$c$$$      print *,' scheck2 ngvecc2s =',sum(ngvecc2(1:3,1:ngc2))
c$$$c$$$      print *," nxmin:nxmax,nymin:nymax,nzmin:nzmax=",
c$$$c$$$     & nxmin,nxmax,nymin,nymax,nzmin,nzmax
c$$$c$$$      endif
c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccc
c$$$      allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
c$$$      if(debug) print *,'melpln2: 1111111____'
c$$$      ipp = -10000000
c$$$      do igc2 =1,ngc2
c$$$        ipp(ngvecc2(1,igc2),ngvecc2(2,igc2),ngvecc2(3,igc2))=igc2
c$$$c        nn(1:3) = ngvecc2(1:3,igc2)
c$$$c        alat =10.26
c$$$c        print *,'  nn =',nn(1:3)
c$$$c        print *,'  abs=',2*pi/alat
c$$$c     &     *sqrt(sum( (qi(1:3)+ matmul(qbas, nn(1:3)))**2 ) )
c$$$      enddo
c$$$      if(debug) print *,'melpln2: 2222222222222'
c$$$c      inum =0
c$$$c
c$$$c      print *,' melpln2: sht nadd=',shtv,nadd
c$$$c      do igc = 1,ngc
c$$$c        qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
c$$$c        phase(igc) = exp(-img*sum(qg*shtv))
c$$$c
c$$$c      print *,' igc phase=', phase(igc)
c$$$c      enddo
c$$$c      print *," melpln2:goto ig1 ig2 loop"
c$$$c #ifdef EXPAND_MELPLN2
c$$$      allocate(gg(ngc2,nt0*ntp0),gg_(ntp0,nt0,ngc2)) !sf
c$$$      gg_=0d0 !sf
c$$$      do ig1 = 1,ngp1
c$$$        do ig2 = 1,ngp2
c$$$          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
c$$$          if(nn(1)<nxmin.or.nn(1)>nxmax) then
c$$$c$$$ccccccccccccccccccccccccc
c$$$c$$$             do igc2 =1,ngc2
c$$$c$$$                print *,'ttt ngvecc2:',ngvecc2(1:3,igc2)
c$$$c$$$             enddo
c$$$c$$$             print *,'ttt ngc2=',ngc2
c$$$c$$$             print *,'ttt mepln2: nn=',nn(1),nxmin,nxmax,ngp1,ngp2
c$$$c$$$c             print *,'ttt qu1=',qu1,ngvecp1(1:3,ig1)
c$$$c$$$c             print *,'ttt qu2=',qu2,ngvecp2(1:3,ig2)
c$$$c$$$             print *,'ttt qi =',qi
c$$$c$$$             print *,'ttt nadd=',nadd
c$$$c$$$cccccccccccccccccccccccccc
c$$$Cstop2rx 2013.08.09 kino             stop "melpln2:nn(1) out of range"
c$$$             call rx( "melpln2:nn(1) out of range")
c$$$          endif   
c$$$Cstop2rx 2013.08.09 kino          if(nn(2)<nymin.or.nn(2)>nymax) stop "melpln2:nn(2) out of range"
c$$$          if(nn(2)<nymin.or.nn(2)>nymax) call rx( "melpln2:nn(2) out of range")
c$$$Cstop2rx 2013.08.09 kino          if(nn(3)<nzmin.or.nn(3)>nzmax) stop "melpln2:nn(3) out of range"
c$$$          if(nn(3)<nzmin.or.nn(3)>nzmax) call rx( "melpln2:nn(3) out of range")
c$$$        enddo
c$$$      enddo
c$$$      if(debug) print *,'melpln2: 3333333333333  ngp nadd',ngp1,ngp2,ntp0,nadd
c$$$!! This loop is time consuming for SiO2c even ngp1=ngp2=283 ntp0 is alrge?
c$$$!$OMP parallel private(nn,igc2,qg,phas,ge2) 
c$$$      qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, -nadd))
c$$$      phasqadd = exp(-img*sum(qg*shtv))
c$$$      do ig1 = 1,ngp1
c$$$          qg(1:3) = 2d0*pi * ( matmul(qbas, ngvecp1(1:3,ig1)) )
c$$$          phas1 = exp(-img*sum(qg*shtv))
c$$$        do ig2 = 1,ngp2
c$$$          qg(1:3) = 2d0*pi * ( matmul(qbas, -ngvecp2(1:3,ig2)) ) 
c$$$          phas2 = exp(-img*sum(qg*shtv))
c$$$c$$$cccccccccccccccccccc
c$$$c$$$          if(ig1==2) then
c$$$c$$$             print *,' start nn=',ig1,ig2,nadd
c$$$c$$$             print *,' ngvecp1(1:3,ig1)',ngvecp1(1:3,ig1)
c$$$c$$$             print *,' ngvecp2(1:3,ig2)',ngvecp2(1:3,ig2)
c$$$c$$$          endif
c$$$c$$$cccccccccccccccccccc
c$$$          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
c$$$          igc2    = ipp(nn(1),nn(2),nn(3))
c$$$c          qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, nn(1:3)) )
c$$$c          phas = exp(-img*sum(qg*shtv))
c$$$          phas= phasqadd*phas1*phas2
c$$$!$OMP do 
c$$$          do it =1,nt0
c$$$            ge2=dconjg(geigq2(ig2,it)) *phas
c$$$            call zaxpy(ntp0,ge2,geigq1_(1,ig1),1,gg_(1,it,igc2),1)
c$$$c            do itp=1,ntp0
c$$$c              gg_(itp,it,igc2) =gg_(itp,it,igc2) 
c$$$c     &        +  geigq1_(itp,ig1)*ge2
c$$$c            enddo
c$$$          enddo
c$$$!$OMP end do 
c$$$        enddo
c$$$      enddo
c$$$!$OMP end parallel 
c$$$      if(debug) print *,'melpln2: 444444444444444 '
c$$$      do itp=1,ntp0
c$$$        do it =1,nt0
c$$$          do igc2=1,ngc2
c$$$            gg(igc2,it+(itp-1)*nt0)=gg_(itp,it,igc2)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      deallocate(gg_)
c$$$c      if(nomatm()==1) then !sep2009 for case to avoid bug in zgemm (called in matm)
c$$$c        zmelp = matmul(ppx(1:ngc,1:ngc2),gg)
c$$$c      else
c$$$      call matm(ppx,gg,zmelp,ngc,ngc2,ntp0*nt0) !ppx= matmul(ppovlinv(ngc,ngc),ppovl(ngc,ngc2)
c$$$c      endif
c$$$      deallocate(ppx,gg,ngvecc2,ipp,geigq1_)
c$$$      zmelp = dconjg(zmelp) ! <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
c      end subroutine melpln2t 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(6, '(" *** q=",(3f10.4,2x))'),q(1:3)
c      do ig1 = 1,ngp1
c        write(6, '(" ig geig =",i2,2x,3i3,4(2f10.4,2x))')
c     &   ig1,ngvecp1(1:3,ig1), ( geigq1(ig1,itp),itp=1,4)
c      enddo
ccccccccccccccccccccccccccccccccc
cccccccccccccccccccc
c      write(6,"( ' --- sum check gg=',10d13.6)")
c     &  sum(gg(1:ngc2,1,1)),sum(abs(gg(1:ngc2,1,1)))
cccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
c      allocate(xxx(ngc,ngc2))
c      call matm(ppovlinv,ppovl,xxx,ngc,ngc,ngc2)
c     print *,'invtest'
c      ix=1
c     do iy=1,ngc2
c       write(6,'("  ix iy xxx=",2i6,2f12.6)')  ix,iy,xxx(ix,iy)
c      enddo
c     stop ' invtest end'
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc
c      do itp=1,ntp0
c      do it =1,nt0
c      do igc = 1,ngc
c      do igc2= 1,ngc2
c        if( sum(abs(ngvecc(:,igc)-ngvecc2(:,igc2)))==0 ) then
c          zmelp(igc, it,itp) = gg(igc2,it,itp)
c        endif
c      enddo
c      enddo
c      enddo
c      enddo
c      zmelp = dconjg(zmelp) ! <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
c      deallocate(ppx,gg,ngvecc2)
c      return
ccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccc
c      write(6, '(" *** q=",(3f10.4,2x))'),q(1:3)
c      do it = 1,nt0
c      do igc = 1,ngc
c        write(6, '(" it igc nvec gg =",2i4,2x,3i3,4(2f10.4,2x))')
c     &   it ,igc, ngvecc(1:3,igc), ( zmelp(igc,it,itp),itp=1,4)
c      enddo
c      enddo
ccccccccccccccccccccccccccccccccc

