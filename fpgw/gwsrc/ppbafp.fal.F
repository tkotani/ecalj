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
      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
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
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) then
        print *,' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
     &  add,nadd,abs(add-nadd)
        stop "drvmlp: abs(add-nadd))>1d-10"
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
!TIME1 "before drv melpln2"
!TIME0 
      call melpln2(ngp1, ngvecpB1
     &           ,  ngp2, ngvecpB2
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,1:nt0), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
     i       shtv, matmul(symope,qik),qik, symope, qbas, 
     o       zmelpl)
!TIME1 "after drv melpln2"
!TIME0 
      if(debug) print *,'end of melpln2'

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
!TIME1 "end drvmelpln2"
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
     &    ngveccB(3,ngcmx,nqibz), ngveccBr(3,ngcmx,nqibz),iq,igc,kp
      real(8)    :: symops(3,3), qq(3), qrot(3),qbas(3,3)
     &             ,qbasinv(3,3),
     &     qfbz(3),rotnvg(3,3),vec(3),xsign(1:3),det,sumchk !,sumchk2
      logical:: debug=.false.
c
c      print *,' rotgvec: '
      call dinv33(qbas,0,qbasinv,det)
!      symops = grp(:,:,irot)
      if(debug) write(6,'(3f12.6)') symops(1,1:3)
      if(debug) write(6,'(3f12.6)') symops(2,1:3)
      if(debug) write(6,'(3f12.6)') symops(3,1:3)

      rotnvg = matmul(qbasinv,matmul(symops,qbas))
c      print *, 'sum rotnvg=',sum(rotnvg)
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
        stop "rotgvec: sum chk error >1d-4"
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
      subroutine melpln2
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
c     &       (qu1,qu2, ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied                        
     i       shtv,q, qi, symope, qbas,
     o       zmelp)
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
co zmelp
cr use PPOVL through getppx in m_read_ppovl
c-------------------------------------------------
      use m_read_ppovl,only:getppx, ppx,ngc2,ngvecc2
      implicit none
      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
     &     itp,it, ngc,nnum,inum,ig1,ig2,igc,ntp0,nt0,
     &     nadd(3)   !ngvecc(3,ngc)
      complex(8):: geigq1(ngp1,ntp0),
     &          geigq2(ngp2,nt0),ge2,zmelp(ngc,nt0*ntp0),  !it was zmelp(ngc,nt0,ntp0); but no effect sep2009.
     &             img=(0d0,1d0) ,phas         ,phase2
      real(8)::  q(3),qi(3),qg(3),qbas(3,3),shtv(3)
     & ,pi=3.1415926535897932D0, symope(3,3)   ,qx(3)
      integer(4),allocatable :: ipp(:,:,:)
      complex(8),allocatable :: gg(:,:),gg_(:,:,:),xxx(:,:)
      integer(4) :: iqi,igc2,ngcx,nn(1:3),
     &       nxmin,nymin,nzmin,nxmax,nymax,nzmax,ix,iy
      integer(4)::verbose,nomatm
      logical:: debug=.false.
      real(8):: qu1(3),qu2(3)
      complex(8),allocatable:: geigq1_(:,:)
      allocate(geigq1_(ntp0,ngp1))
c------------------------------
      if(verbose()>=100) print *,' melpln2: '

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      if(sum(abs(qu1-(/9.6225044864937631E-002,   9.6225044864937631E-002, 0.10000000000000001/)))<1d-8) then
c$$$         debug=.true.
c$$$      endif
c$$$      if(debug) print *,'melpln2: qu1',qu1
c$$$      if(debug) print *,'melpln2: qu2',qu2
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccc

c- Get ppx and ngvecc2
c   ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
ccc!$OMP parallel
ccc!$OMP master
      call getppx(qbas,  ngcx,qi,symope) ! rppovl.f
ccc!$OMP end master
ccc!$OMP end parallel 
      if(ngc/=ngcx) then
        print *, qi, ngcx, ngc
        stop 'melpln2: ngc/= ngcx by getppx:PPOVL'
      endif
      do ig1=1,ngp1
        do itp=1,ntp0
          geigq1_(itp,ig1)=geigq1(ig1,itp)
        enddo
      enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      if(debug) then
c$$$       print *,' melpln2 old sumcheck ngc2 ngvec ppx ='
c$$$     &  ,ngc2,sum(abs(ngvecc2(1:3,1:ngc2))),sum(abs(ppx(1:ngc,1:ngc2)))
c$$$       print *,'melpln2: 1111111111'
c$$$       endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nxmin =minval(ngvecc2(1,1:ngc2)) ;nxmax=maxval(ngvecc2(1,1:ngc2))
      nymin =minval(ngvecc2(2,1:ngc2)) ;nymax=maxval(ngvecc2(2,1:ngc2))
      nzmin =minval(ngvecc2(3,1:ngc2)) ;nzmax=maxval(ngvecc2(3,1:ngc2))

c$$$ccccccccccccccccccccccccccccccccccccccccccc
c$$$      if(debug) then
c$$$      write(6,"(' iqi ngc ngc2=',3i8)") iqi,ngc,ngc2
c$$$c      print *,' scheck2 ppx=',sum(ppovlinv(1:ngcx,1:ngcx))
c$$$c      print *,' scheck2 ppovls   =',sum(ppovl   (1:ngcx,1:ngc2))
c$$$      print *,' scheck2 ngvecc2s =',sum(ngvecc2(1:3,1:ngc2))
c$$$      print *," nxmin:nxmax,nymin:nymax,nzmin:nzmax=",
c$$$     & nxmin,nxmax,nymin,nymax,nzmin,nzmax
c$$$      endif
c$$$cccccccccccccccccccccccccccccccccccccccccccc

      allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
      if(debug) print *,'melpln2: 1111111____'
      ipp = -10000000
      do igc2 =1,ngc2
        ipp(ngvecc2(1,igc2),ngvecc2(2,igc2),ngvecc2(3,igc2))=igc2
c        nn(1:3) = ngvecc2(1:3,igc2)
c        alat =10.26
c        print *,'  nn =',nn(1:3)
c        print *,'  abs=',2*pi/alat
c     &     *sqrt(sum( (qi(1:3)+ matmul(qbas, nn(1:3)))**2 ) )
      enddo
      if(debug) print *,'melpln2: 2222222222222'

c      inum =0
c
c      print *,' melpln2: sht nadd=',shtv,nadd
c      do igc = 1,ngc
c        qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
c        phase(igc) = exp(-img*sum(qg*shtv))
c
c      print *,' igc phase=', phase(igc)
c      enddo
c      print *," melpln2:goto ig1 ig2 loop"
c #ifdef EXPAND_MELPLN2

      allocate(gg(ngc2,nt0*ntp0),gg_(ntp0,nt0,ngc2)) !sf
      gg_=0d0 !sf
      do ig1 = 1,ngp1
        do ig2 = 1,ngp2
          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
          if(nn(1)<nxmin.or.nn(1)>nxmax) then
c$$$ccccccccccccccccccccccccc
c$$$             do igc2 =1,ngc2
c$$$                print *,'ttt ngvecc2:',ngvecc2(1:3,igc2)
c$$$             enddo
c$$$             print *,'ttt ngc2=',ngc2
c$$$             print *,'ttt mepln2: nn=',nn(1),nxmin,nxmax,ngp1,ngp2
c$$$c             print *,'ttt qu1=',qu1,ngvecp1(1:3,ig1)
c$$$c             print *,'ttt qu2=',qu2,ngvecp2(1:3,ig2)
c$$$             print *,'ttt qi =',qi
c$$$             print *,'ttt nadd=',nadd
c$$$cccccccccccccccccccccccccc
             stop "melpln2:nn(1) out of range"
          endif   
          if(nn(2)<nymin.or.nn(2)>nymax) stop "melpln2:nn(2) out of range"
          if(nn(3)<nzmin.or.nn(3)>nzmax) stop "melpln2:nn(3) out of range"
        enddo
      enddo
      if(debug) print *,'melpln2: 3333333333333  ngp nadd',ngp1,ngp2,nadd

!$OMP parallel private(nn,igc2,qg,phas,ge2) 
      do ig1 = 1,ngp1
        do ig2 = 1,ngp2
c$$$cccccccccccccccccccc
c$$$          if(ig1==2) then
c$$$             print *,' start nn=',ig1,ig2,nadd
c$$$             print *,' ngvecp1(1:3,ig1)',ngvecp1(1:3,ig1)
c$$$             print *,' ngvecp2(1:3,ig2)',ngvecp2(1:3,ig2)
c$$$          endif
c$$$cccccccccccccccccccc
          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
          igc2    = ipp(nn(1),nn(2),nn(3))
          qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, nn(1:3)) )
          phas = exp(-img*sum(qg*shtv))
!$OMP do 
          do it =1,nt0
            ge2=dconjg(geigq2(ig2,it)) *phas
            do itp=1,ntp0
              gg_(itp,it,igc2) =gg_(itp,it,igc2) 
     &        +  geigq1_(itp,ig1)*ge2
            enddo
          enddo
!$OMP end do 
        enddo
      enddo
!$OMP end parallel 

      if(debug) print *,'melpln2: 444444444444444 '
      do itp=1,ntp0
        do it =1,nt0
          do igc2=1,ngc2
            gg(igc2,it+(itp-1)*nt0)=gg_(itp,it,igc2)
          enddo
        enddo
      enddo
      deallocate(gg_)
c      if(nomatm()==1) then !sep2009 for case to avoid bug in zgemm (called in matm)
c        zmelp = matmul(ppx(1:ngc,1:ngc2),gg)
c      else
      call matm(ppx,gg,zmelp,ngc,ngc2,ntp0*nt0) !ppx= matmul(ppovlinv(ngc,ngc),ppovl(ngc,ngc2)
c      endif
      deallocate(ppx,gg,ngvecc2,ipp,geigq1_)
      zmelp = dconjg(zmelp) ! <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
      end subroutine melpln2 


ccccccccccccccccccccccccccccccccc
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

