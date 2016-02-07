      subroutine mkppmt(alat,plat,qlat, q, ng1,ngvec1, 
     i     rmax,  nbas,  bas, lmxa, lmxax, 
     o         ppmt)
C- ppmt value and slope at MT bounadaries for q+G plane waves.
Co  ppmt(1,lm,nbas): value for (lm,ibas)
Co  ppmt(2,lm,nbas):
Cr Oct2005
C------------------
      use m_lldata,only: ll
      implicit none
      integer(4) ::  nbas,  ng1,lmxax
     &       ,ig1,ibas,la,lma,lmh, 
     &        ngvec1(3,ng1),lmxa(nbas)
      real(8) :: absqg, rmax(nbas),pi4,r2s,
     &  q(3),plat(3,3),qlat(3,3),qgg(3),
     &  pi,alat,tpiba, bas(3,nbas) ,facl
      complex(8) :: ppmt(2,(lmxax+1)**2,nbas,ng1)
      integer(4) :: verbose,it,ip
      real(8),allocatable:: cy(:),yl(:),ylr(:)
      real(8) :: ak1(200),aj1(200), dk1(200),dj1(200),tpi,absqg2,qqq(3)
      complex(8) :: img =(0d0,1d0),fac,phase

      real(8):: theta,phi,rdir(3)
      complex(8)::valx1,valx2

c$$$#ifdef COMMONLL
c$$$      integer(4):: ll(51**2)
c$$$      common/llblock/ll
c$$$#else
c$$$      integer(4)::ll
c$$$#endif
      logical ::debug=.false.
c-----------------------------------------------------
      if(debug) allocate(ylr((lmxax+1)**2) )
      if(verbose()>50) print *,' mkppmt:'
      pi  = 4d0*datan(1d0)
      tpi = 2d0*pi
      pi4 = 4d0*pi
      tpiba = tpi/alat
c      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
      call sylmnc(cy,lmxax)
      ppmt = 0d0
      do ig1 =1,ng1
        qgg(1:3)=  tpiba*( q(1:3)+ matmul(qlat, ngvec1(1:3,ig1)) )
!   exp(i qgg *\bfr)

        absqg2  = sum(qgg(1:3)**2)
        absqg = sqrt(absqg2)
        if(absqg==0d0) then
          qqq = (/0d0,0d0,1d0/)
        else
          qqq = qgg/absqg
        endif
        call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )

        do ibas = 1,nbas
          phase = exp( img*sum( qgg*bas(:,ibas) )*alat  )
          call radkj(absqg2, rmax(ibas),lmxa(ibas),ak1,aj1,dk1,dj1,0)
          do lma = 1,(lmxa(ibas)+1)**2
            la = ll(lma)
            fac = pi4* img**la * cy(lma)*yl(lma) * facl(absqg,la)
            ppmt(1,lma,ibas,ig1) = fac*phase* aj1(la+1)
            ppmt(2,lma,ibas,ig1) = fac*phase* dj1(la+1)
          enddo
        enddo

c--- debug check ! val at mt in two kinds of calculations.
        if(debug) then
          do it= 0,11
            do ip= 0,13
              do ibas = 1,nbas
                phase = exp( img*sum( qgg*bas(:,ibas) )*alat  )
                theta = pi/2d0* it/19d0
                phi   = 2*pi  * ip/17d0
                rdir  = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
                call sylm(rdir,ylr,lmxax,r2s) !spherical factor Y( q+G )
                rdir  = rdir*rmax(ibas)
                valx1 = exp(img* sum( qgg*( bas(:,ibas)*alat+rdir) ) )
                valx2 = 0d0
                do lma = 1,(lmxa(ibas)+1)**2
                  valx2= valx2 + ppmt(1,lma,ibas,ig1)*cy(lma) *ylr(lma)
                enddo
                if(abs(valx2 - valx1)/abs(valx1)>0) !.05) 
     &      write(6,"(' ig itip ibas valx1 err=',4i3,2d11.3,3x,d11.2)")
     &       ig1,it,ip,ibas, valx1, abs(valx2 - valx1)/abs(valx1)
              enddo
            enddo
          enddo
          write(6,*)
        endif
c--- debug check end
      enddo
Cstop2rx 2013.08.09 kino      if(debug) stop 'test end ----------'
      if(debug) call rx( 'test end ----------')
      deallocate(yl,cy)
c          call wronkj( absqg1**2, absqg2**2, rmax(ibas),
c     &                 lmxa(ibas),fkk,fkj,fjk,fjj)
c
cyy        do lma  = 1,(lmxa(ibas)+1)**2
cyy          la = ll(lma)
cyy          ppovl(ig1,ig2) = ppovl(ig1,ig2)
cyy     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
cyy     &       * pi4* pi4 * cy(lma)**2*yl1(lma)*yl2(lma)
cyy     &       * (-fjj(la))*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
c        do la  = 0,lmxa(ibas)
c          ppmt(ig1,la,ibas) = aj1(la+1)
c          ppmt(ig1,la,ibas) = dj1(la+1)
c          ppovl(ig1,ig2) = ppovl(ig1,ig2)
c     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
c     &        * pi4 *(2*la+1d0) * plegn(la,cost)
c     &        * (-fjj(la))  *facl(absqg1*absqg2,la)
c
c                           !*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
cxxx     &      * ajaj(ig1,ig2,la,ibas)        ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
cccccccccccc
c        if( abs(absqg1)==0d0)
c     &    write(1126,'(4i3,2d15.6,3x,d15.6)')ig1,ig2,la,ibas
c     &    , -fjj(la),ajaj(ig1,ig2,la,ibas)
c     &    , fjj(la) +ajaj(ig1,ig2,la,ibas)
cccccccccccc
      end
C-------------------------------------------------------------------------

c$$$      subroutine pplmat2(s_lat, q, ng, ngvec, rmax, 
c$$$     &       bmat,kmaxx, nlmax, nbas, ndimh, bas, lmxa, lmxax,
c$$$     &       lh1,lh2,rsmh, eh,ntorb,ntorbx, kmax,
c$$$     &       plhd, plpkl, QpGcutHankel, 
c$$$     o       ppovl,ppsh)
c$$$c <P^q1_G1 | P^q2_G2 > matrix where P^q1_G1 denotes projected plane wave, zero within sphere.
c$$$c <P^q1_G1 | smooth Hankel > matrix
c$$$c
c$$$c Mar-2001 from pplmat
c$$$      implicit none
c$$$      integer(4) ::  nbas, ndimh, ng,kmaxx,nlmax,lmxax,ntorbx
c$$$     &       ,ig1,ig2,ibas,ih,la,ihoff,itorb,k,lma,lmh,lh,ntorb(nbas)
c$$$     &      , ngvec(3,ng),lmxa(nbas),kmax(nbas),
c$$$     &      lh1(ntorbx,nbas),lh2(ntorbx,nbas)
c$$$      real(8)    :: absqg1,absqg2,tripl,rmax(nbas),pi4,
c$$$     &              rsmh(ntorbx,nbas), eh(ntorbx,nbas),r2s,denom,gan
c$$$      complex(8) :: bmat(0:kmaxx, nlmax, nbas, ndimh)
c$$$      real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qg(3,ng)    ,facl
c$$$     & ,pi,alat,tpiba,cost,voltot,plegn,qg1(3),qg2(3),absqg(ng),qqq(3),
c$$$     &  bas(3,nbas), fkk(0:lmxax),fkj(0:lmxax),fjk(0:lmxax),fjj(0:lmxax)
c$$$      complex(8) :: img =(0d0,1d0),phase
c$$$c
c$$$      complex(8) ::  plhd (ng, ndimh)        ! integral plane \times head ???
c$$$     &             , plpkl(ng, nlmax, 0:kmaxx,nbas)  ! integral plane \times poly ???
c$$$      complex(8) :: ppovl(ng,ng), ppsh(ng,ndimh)
c$$$      real(8),allocatable::cy(:),yl(:),yl1(:),yl2(:)
c$$$      integer(4),allocatable :: ngvecx(:,:)
c$$$      complex(8),allocatable :: ppshx(:,:),ppox(:,:,:)
c$$$      complex(8),allocatable :: ppovlx(:,:)
c$$$      real(8) :: gv1(3),gv1x(3),QpGcutHankel,dummy,tpibaqlat(3,3)
c$$$      integer(4) :: ngmx,ig1x,nx(3),n1x,n2x,n3x,n1m,n2m,n3m
c$$$#ifdef COMMONLL
c$$$      integer(4):: ll(51**2)
c$$$      common/llblock/ll
c$$$#else
c$$$      integer(4)::ll
c$$$#endif
c$$$c-----------------------------------------------------------------------
c$$$      print *,' pplmat2:'
c$$$      pi  = 4d0*datan(1d0)
c$$$      pi4 = 16d0*datan(1d0)
c$$$      call u_lat_vecs(s_lat,alat,plat,qlat)
c$$$      tpiba  = 2*pi/alat
c$$$      tpibaqlat = tpiba *qlat
c$$$      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
c$$$      allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
c$$$      call sylmnc(cy,lmxax)
c$$$
c$$$ccc      goto 1201 ! if you want to use old version ---
c$$$c
c$$$c <P^q_G1 | P^q_G2 >
c$$$ccctest      call mkppovl2test(alat,plat,qlat,
c$$$      call mkppovl2(alat,plat,qlat, 
c$$$     &    ng, ngvec, !G1
c$$$     &    ng, ngvec, !G2
c$$$     &    nbas, rmax, bas,  
c$$$     o    ppovl)
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$c      write(3001,"('q=',3d16.8)") q
c$$$c      do ig1 = 1, ng
c$$$c      do ig2 = 1, ng
c$$$c        write(3001,"(2i3,'  ',2d16.8)") ig1,ig2,ppovl(ig1,ig2)
c$$$c      enddo
c$$$c      enddo
c$$$c      stop 'test1 end'
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$c Get ngvecx for the interstitial part of the smooth hankel function.
c$$$      call getgv2(alat,plat,qlat, q, QpGcutHankel, 1, ngmx,dummy)
c$$$      allocate( ngvecx(3,ngmx),ppshx(ngmx,ndimh) )
c$$$      call getgv2(alat,plat,qlat, q, QpGcutHankel, 2, ngmx, ngvecx)  ! for eigenfunction
c$$$      print *,  ' pplmat2: ngmx=',ngmx
c$$$
c$$$c Expansion coefficients of |P^q_G1x> in the expansion .. |Smooth Hankel>
c$$$c ---> ppshx(igx1,ndimh) matrix
c$$$      ppshx = 0d0
c$$$      do ig1x = 1, ngmx
c$$$        qg1(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecx(1:3,ig1x)))
c$$$        absqg1   = sqrt(sum(qg1(1:3)**2))
c$$$        if(absqg1==0d0) then
c$$$          qqq = (/0d0,0d0,1d0/)
c$$$        else
c$$$          qqq = qg1/absqg1
c$$$        endif
c$$$        call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )
c$$$        ihoff = 0
c$$$        do ibas = 1,nbas
c$$$          if(ibas >=2) then
c$$$            ihoff  = ihoff  
c$$$     &      +  sum ( (  lh2(1:ntorb(ibas-1),ibas-1)+1)**2 
c$$$     &              - lh1(1:ntorb(ibas-1),ibas-1)**2   )
c$$$          endif
c$$$          ih = ihoff
c$$$          phase = exp( -img*sum( qg1*bas(:,ibas)*alat )  )
c$$$          do itorb = 1,ntorb(ibas)
c$$$            do lmh   = lh1(itorb,ibas)**2+1, (lh2(itorb,ibas)+1)**2
c$$$              lh = ll(lmh)
c$$$              ih = ih + 1
c$$$              denom = eh(itorb,ibas)-absqg1**2
c$$$              gan   = 1d0/4d0*rsmh(itorb,ibas)**2
c$$$              ppshx(ig1x,ih) = - pi4/denom
c$$$     &      * (-img)**lh * facl(absqg1,lh)
c$$$     &      * cy(lmh)* yl(lmh)* exp( gan*denom)
c$$$     &      * phase  /voltot
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$c
c$$$      allocate(ppovlx(ng,ngmx))
c$$$ccctest      call mkppovl2test(alat,plat,qlat,
c$$$      call mkppovl2(alat,plat,qlat, 
c$$$     &    ng,   ngvec,  !G1
c$$$     &    ngmx, ngvecx, !G1x
c$$$     &    nbas, rmax, bas,  
c$$$     o    ppovlx)
c$$$      call matm(ppovlx,ppshx,ppsh,ng,ngmx,ndimh)
c$$$      deallocate( ppovlx )
c$$$      if (allocated(cy)) deallocate(cy)
c$$$      if (allocated(yl)) deallocate(yl)
c$$$      if (allocated(ngvecx)) deallocate(ngvecx)
c$$$      if (allocated(ppshx)) deallocate(ppshx)
c$$$      return
c$$$
c$$$c <P^q_G1 | Smooth Hankel > = ppsh matrix
c$$$c     n1x = maxval( ngvecx(1,:)) - minval( ngvec(1,:))
c$$$c     n1m = minval( ngvecx(1,:)) - maxval( ngvec(1,:))
c$$$c     n2x = maxval( ngvecx(2,:)) - minval( ngvec(2,:))
c$$$c     n2m = minval( ngvecx(2,:)) - maxval( ngvec(2,:))
c$$$c     n3x = maxval( ngvecx(3,:)) - minval( ngvec(3,:))
c$$$c     n3m = minval( ngvecx(3,:)) - maxval( ngvec(3,:))
c$$$c      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
c$$$c      ppox = 1d99
c$$$c   <G1|G1X>
c$$$c      do ig1  = 1, ng
c$$$c      do ig1x = 1, ngmx
c$$$c       nx(1:3) = ngvecx(1:3,ig1x) -ngvec(1:3,ig1)
c$$$c        if( ppox(nx(1),nx(2),nx(3))==1d99 ) then
c$$$c          call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
c$$$c     i    nx(1:3), ! G1x -G1
c$$$c     o    ppox( nx(1),nx(2),nx(3)))
c$$$c       endif
c$$$c      enddo
c$$$c     enddo
c$$$c      ppsh = 0d0
c$$$c      do ih   = 1, ndimh
c$$$c      do ig1  = 1, ng
c$$$c      do ig1x = 1, ngmx
c$$$c       nx(1:3) = ngvecx(1:3,ig1x) - ngvec(1:3,ig1)
c$$$c       ppovlx  = ppox( nx(1),nx(2),nx(3))
c$$$c       ppsh(ig1,ih) = ppsh(ig1,ih) + ppovlx * ppshx(ig1x,ih)
c$$$c      enddo
c$$$c     enddo
c$$$c     enddo
c$$$c      deallocate( ppox )
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$
c$$$
c$$$c=======================================================================
c$$$ 1201 continue
c$$$c      print *,  ' pplmat: xxx voltot=',voltot
c$$$c <P^q_G1 | P^q_G2 > matrix
c$$$      do ig1 = 1,ng
c$$$        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvec(1:3,ig1)))
c$$$        absqg(ig1)  = sqrt(sum(qg(1:3,ig1)**2))
c$$$      enddo
c$$$
c$$$c --- OLD codes from here ----------------------------------------------
c$$$      ppovl = 0d0
c$$$      do ig1 = 1,ng
c$$$        ppovl(ig1,ig1)= voltot
c$$$      enddo
c$$$
c$$$cyy      allocate(yl1((lmxax+1)**2),yl2((lmxax+1)**2) )
c$$$
c$$$      do ig1 =1,ng
c$$$        qg1(1:3) = qg(1:3,ig1)
c$$$        absqg1   = absqg(ig1)
c$$$cyy        call sylm(qg1/absqg1,yl1,lmxax,r2s) !spherical factor Y( q+G )
c$$$        do ig2 =1,ng
c$$$          qg2(1:3) = qg(1:3,ig2)
c$$$          absqg2   = absqg(ig2)
c$$$cyy        call sylm(qg2/absqg2,yl2,lmxax,r2s) !spherical factor Y( q+G )
c$$$
c$$$          if( absqg1*absqg2 == 0d0) then
c$$$            cost = 1d0
c$$$          else
c$$$            cost = sum(qg1*qg2)/absqg1/absqg2
c$$$          endif
c$$$
c$$$          do ibas = 1,nbas
c$$$            call wronkj( absqg1**2, absqg2**2, rmax(ibas),
c$$$     &                 lmxa(ibas),fkk,fkj,fjk,fjj)
c$$$
c$$$cyy        do lma  = 1,(lmxa(ibas)+1)**2
c$$$cyy          la = ll(lma)
c$$$cyy          ppovl(ig1,ig2) = ppovl(ig1,ig2)
c$$$cyy     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
c$$$cyy     &       * pi4* pi4 * cy(lma)**2*yl1(lma)*yl2(lma)
c$$$cyy     &       * (-fjj(la))*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
c$$$
c$$$            do la  = 0,lmxa(ibas)
c$$$              ppovl(ig1,ig2) = ppovl(ig1,ig2)   
c$$$     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
c$$$     &        * pi4 *(2*la+1d0) * plegn(la,cost) 
c$$$     &        * (-fjj(la))  *facl(absqg1*absqg2,la)
c$$$              !*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
c$$$cxxx     &      * ajaj(ig1,ig2,la,ibas)        ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
c$$$
c$$$cccccccccccc
c$$$c        if( abs(absqg1)==0d0)
c$$$c     &    write(1126,'(4i3,2d15.6,3x,d15.6)')ig1,ig2,la,ibas
c$$$c     &    , -fjj(la),ajaj(ig1,ig2,la,ibas)
c$$$c     &    , fjj(la) +ajaj(ig1,ig2,la,ibas)
c$$$cccccccccccc
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$c      write(3002,"('q=',3d16.8)") q
c$$$c      do ig1 = 1, ng
c$$$c      do ig2 = 1, ng
c$$$c        write(3002,"(2i3,'  ',2d16.8)") ig1,ig2,ppovl(ig1,ig2)
c$$$c      enddo
c$$$c      enddo
c$$$c      stop 'test2 end'
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$cccccccccccccc
c$$$c      do ig1=1,ng-1
c$$$c       ig2=ig1+1
c$$$c        write(6,"(' old ppovl(ig1,ig2)=',i4,2d13.6)")ig1,ppovl(ig1,ig2)
c$$$c      enddo
c$$$c     stop "------testen--------d"
c$$$cccccccccccccc
c$$$
c$$$ 1111 continue
c$$$c <P^q_G1 | Smooth Hankel > = ppsh matrix
c$$$      ppsh=0d0
c$$$      do ig1 = 1, ng
c$$$        qg1(1:3) = qg(1:3,ig1)
c$$$        absqg1   = absqg(ig1)
c$$$        if(absqg1==0d0) then
c$$$          qqq = (/0d0,0d0,1d0/)
c$$$        else
c$$$          qqq = qg1/absqg1
c$$$        endif
c$$$        call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )
c$$$
c$$$        ihoff  = 0
c$$$        do ibas = 1,nbas
c$$$          if(ibas >=2) then
c$$$            ihoff  = ihoff  
c$$$     &      +  sum ( (  lh2(1:ntorb(ibas-1),ibas-1)+1)**2 
c$$$     &                - lh1(1:ntorb(ibas-1),ibas-1)**2   )
c$$$          endif
c$$$c from head
c$$$          phase = exp( -img*sum( qg1*bas(:,ibas)*alat )  )
c$$$          ih = ihoff
c$$$          do itorb = 1,ntorb(ibas)
c$$$            do lmh   = lh1(itorb,ibas)**2+1, (lh2(itorb,ibas)+1)**2
c$$$              lh = ll(lmh)
c$$$              ih = ih + 1
c$$$c           print *,ibas,ip,itorb,lmh,ias,ih,lh
c$$$c           print *,lmh,lh
c$$$              denom = eh(itorb,ibas)-absqg1**2
c$$$              gan   = 1d0/4d0*rsmh(itorb,ibas)**2
c$$$              ppsh(ig1,ih) = ppsh(ig1,ih) - pi4/denom
c$$$c     &         * (-img*absqg1)**lh * cy(lmh)* yl(lmh)* exp( gan*denom)
c$$$     &         * (-img)**lh * facl(absqg1,lh)
c$$$     &         * cy(lmh)* yl(lmh)* exp( gan*denom)
c$$$     &         * phase
c$$$     &         - phase * plhd(ig1, ih)  ! integral of plane \times hd
c$$$            enddo
c$$$          enddo
c$$$c from tail
c$$$          do lma= 1, (lmxa(ibas)+1)**2
c$$$            do k  = 0, kmax(ibas)
c$$$              la  = ll(lma)
c$$$c           print *, ibas,ip,lma,k
c$$$              ppsh(ig1,1:ndimh) = ppsh(ig1,1:ndimh)  
c$$$     &        - phase * plpkl(ig1, lma, k,ibas)
c$$$     &                * bmat(k, lma, ibas, 1:ndimh) ! integral of plane \times hd
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      if (allocated(cy)) deallocate(cy)
c$$$      if (allocated(yl)) deallocate(yl)
c$$$      if (allocated(ngvecx)) deallocate(ngvecx)
c$$$      if (allocated(ppshx)) deallocate(ppshx)
c$$$      end


c--------------------------------------------------------------
      subroutine mkppovl(alat,plat,qlat, q, ng1,ngvec1,ng2,ngvec2, rmax, 
     &         nbas,  bas, lmxa, lmxax, 
     o         ppovl)
c <P^q1_G1 | P^q2_G2 > matrix where P^q1_G1 denotes IPW, zero within sphere.
      implicit none
      integer(4) ::  nbas,  ng1,ng2,lmxax
     &       ,ig1,ig2,ibas,la,lma,lmh, 
     &        ngvec1(3,ng1),ngvec2(3,ng2),lmxa(nbas),ll
      real(8) :: absqg1,absqg2,tripl,rmax(nbas),pi4,r2s
      real(8) :: q(3),plat(3,3),qlat(3,3),qgg1(3,ng1),qgg2(3,ng2), facl
     & ,pi,alat,tpiba,cost,voltot,plegn,qg1(3),qg2(3),qqq(3),
     &  bas(3,nbas), fkk(0:lmxax),fkj(0:lmxax),fjk(0:lmxax),fjj(0:lmxax)
     &  ,absqg1x(ng1),absqg2x(ng2)
      real(8),allocatable::cy(:),yl(:),yl1(:),yl2(:)
      complex(8) :: img =(0d0,1d0),phase
c
      complex(8) :: ppovl(ng1,ng2)
      integer(4):: verbose
c-----------------------------------------------------
      if(verbose()>50) print *,' mkppovl:'
      pi = 4d0*datan(1d0)
      pi4=16d0*datan(1d0)
c      call u_lat_vecs(s_lat,alat,plat,qlat)
      tpiba=2*pi/alat
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
      call sylmnc(cy,lmxax)

c <P^q_G1 | P^q_G2 > matrix
      do ig1 = 1,ng1
ccccccccccccccccccccccccccccc
c       print *, 'mkppovl:xx1 ', tpiba ,q(1:3),ig1
c       print *, 'mkppovl:xx2 ', qlat, ngvec1(1:3,ig1)
cccccccccccccccccccccccccccc
        qgg1(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvec1(1:3,ig1)))
        absqg1x(ig1)  = sqrt(sum(qgg1(1:3,ig1)**2))
c        ppovl(ig1,ig1)= voltot
      enddo

c      print *,  ' mkppovl: 3'

      do ig2 = 1,ng2
        qgg2(1:3,ig2) = tpiba * (q(1:3)+ matmul(qlat, ngvec2(1:3,ig2)))
        absqg2x(ig2)  = sqrt(sum(qgg2(1:3,ig2)**2))
c        ppovl(ig1,ig1)= voltot
      enddo

cyy      allocate(yl1((lmxax+1)**2),yl2((lmxax+1)**2) )

c      print *,  ' pplmat: xxx voltot=',voltot
c
      ppovl = 0d0
      do ig1 =1,ng1
        qg1(1:3) = qgg1(1:3,ig1)
        absqg1   = absqg1x(ig1)
cyy        call sylm(qg1/absqg1,yl1,lmxax,r2s) !spherical factor Y( q+G )
        do ig2 =1,ng2
          qg2(1:3) = qgg2(1:3,ig2)
          absqg2   = absqg2x(ig2)
cyy        call sylm(qg2/absqg2,yl2,lmxax,r2s) !spherical factor Y( q+G )
c
c        print *,  ' ig1 ig2=',ig1,ig2
c
          if(sum(abs(ngvec1(:,ig1)-ngvec2(:,ig2)))==0) 
     &    ppovl(ig1,ig2)= voltot
c
          if( absqg1*absqg2 == 0d0) then
            cost = 1d0
          else
            cost = sum(qg1*qg2)/absqg1/absqg2
          endif
c
          do ibas = 1,nbas
            call wronkj( absqg1**2, absqg2**2, rmax(ibas),
     &                 lmxa(ibas),fkk,fkj,fjk,fjj)
cyy        do lma  = 1,(lmxa(ibas)+1)**2
cyy          la = ll(lma)
cyy          ppovl(ig1,ig2) = ppovl(ig1,ig2)
cyy     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
cyy     &       * pi4* pi4 * cy(lma)**2*yl1(lma)*yl2(lma)
cyy     &       * (-fjj(la))*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)

            do la  = 0,lmxa(ibas)
              ppovl(ig1,ig2) = ppovl(ig1,ig2)   
     &      - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )  !???
     &        * pi4 *(2*la+1d0) * plegn(la,cost) 
     &        * (-fjj(la))  *facl(absqg1*absqg2,la)
              !*(absqg1*absqg2)**la  ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
cxxx     &      * ajaj(ig1,ig2,la,ibas)        ! radjjint =\int_0^a r^2 dr j_l(absqg1 r) j_l(absqg2 r)
cccccccccccc
c        if( abs(absqg1)==0d0)
c     &    write(1126,'(4i3,2d15.6,3x,d15.6)')ig1,ig2,la,ibas
c     &    , -fjj(la),ajaj(ig1,ig2,la,ibas)
c     &    , fjj(la) +ajaj(ig1,ig2,la,ibas)
cccccccccccc
            enddo
          enddo
        enddo
      enddo
      if (allocated(cy)) deallocate(cy)
      if (allocated(yl)) deallocate(yl)
      end

c--------------------------------------------------------------
      subroutine mkppovl2(alat,plat,qlat, ng1,ngvec1,ng2,ngvec2,
     i         nbas, rmax, bas, 
     o         ppovl)
c < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
c
      implicit none
      integer(4) ::  nbas, ng1,ng2,nx(3),
     &        ig1,ig2,ibas, ngvec1(3,ng1),ngvec2(3,ng2),
     &         n1x,n2x,n3x,n1m,n2m,n3m,verbose
      real(8) :: absqg1,absqg2,tripl,rmax(nbas),pi
      real(8) :: plat(3,3),qlat(3,3),
     &  alat,tpiba,tpibaqlat(3,3),voltot, bas(3,nbas)
      complex(8) :: img =(0d0,1d0)
      complex(8) :: ppovl(ng1,ng2)
      complex(8),allocatable :: ppox(:,:,:)
      logical :: debug=.false.
c-----------------------------------------------------
      if(verbose()>50) debug=.true.
      if(debug) print *,' mkppovl2:'
      pi        = 4d0*datan(1d0)
      voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      tpibaqlat =  2*pi/alat *qlat
c < G1 | G2 >
      n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
      n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
      n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
      n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
      n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
      n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
c
      if(debug) print *,' mkppovl2: 1 ',n1x,n1m,n2x,n2m,n3x,n3m
      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
      ppox = 1d99
      do ig1  = 1, ng1
        do ig2  = 1, ng2
          nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
          if( ppox(nx(1),nx(2),nx(3))==1d99 ) then
            call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
     i    nx(1:3), ! G2 -G1
     o    ppox( nx(1),nx(2),nx(3)))
          endif
        enddo
      enddo
      if(debug) print *,' mkppovl2: 2 ',n1x,n1m,n2x,n2m,n3x,n3m
      do ig1 = 1,ng1
        do ig2 = 1,ng2
          nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
          ppovl(ig1,ig2) = ppox( nx(1),nx(2),nx(3) )
        enddo
      enddo
      deallocate(ppox)
c      print *,' mkppovl2: 3 ',n1x,n1m,n2x,n2m,n3x,n3m
      end

c---------------------
      subroutine matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
     i   ngvec,
     o   ppovl)
c ppovl = \int_{cell} IPW(G) d^3r
c G = ggvec = g2 -g1 <g1|g2>
c
      implicit none
      integer(4) :: nbas,ibas,ngvec(3)
      real(8) :: ggvec(3), alat,bas(3,nbas),rmax(nbas),voltot,absg,grmx
     &           ,tpibaqlat(3,3)
      complex(8) :: ppovl
      real(8)    :: pi4 =4d0*3.1415926535897932D0
      complex(8) :: img =(0d0,1d0)
      ggvec(1:3) =  matmul(tpibaqlat, ngvec(1:3))
      absg       =  sqrt(sum(ggvec(1:3)**2))
      ppovl = 0d0
      if(absg==0d0) ppovl =  voltot
      do ibas = 1,nbas
        if(absg==0d0) then
          ppovl = ppovl - pi4*rmax(ibas)**3/3d0
        else

ccccccccccccccccccccccccccc
ctest
c          ppovl = 0d0
c         cycle
ccccccccccccccccccccccccccc

          grmx  = absg* rmax(ibas)
          ppovl = ppovl
     &      - exp( img* sum(ggvec*bas(1:3,ibas))*alat )  
     &        * pi4/absg**3 * ( -grmx * cos(grmx) + sin(grmx) )
        endif
      enddo
      end

c---------------------
      real(8) function facl(a,l)
      integer l
      real(8) a
      if(l==0) then
        facl=1d0
      else
        facl=a**l
      endif
      end















c--------------------------------------------------------------
c--------------------------------------------------------------
c--------------------------------------------------------------
c--------------------------------------------------------------
      subroutine mkppovl2test(alat,plat,qlat, ng1,ngvec1,ng2,ngvec2,
     i         nbas, rmax, bas, 
     o         ppovl)
c < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
c
      use m_lldata,only: ll
      implicit none
      integer(4) ::  nbas, ng1,ng2,nx(3),
     &        ig1,ig2,ibas, ngvec1(3,ng1),ngvec2(3,ng2),
     &         n1x,n2x,n3x,n1m,n2m,n3m
      real(8) :: absqg1,absqg2,tripl,rmax(nbas),pi
      real(8) :: plat(3,3),qlat(3,3),
     &  alat,tpiba,tpibaqlat(3,3),voltot, bas(3,nbas)
      complex(8) :: img =(0d0,1d0)
      complex(8) :: ppovl(ng1,ng2)
      complex(8),allocatable :: ppox(:,:,:)
c$$$#ifdef COMMONLL
c$$$      integer(4) ll(51**2)
c$$$      common/llblock/ll
c$$$#else
c$$$      integer(4) ll
c$$$#endif
c-----------------------------------------------------
      print *,' mkppovl2test:'
c      pi        = 4d0*datan(1d0)
c      voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
c      tpibaqlat =  2*pi/alat *qlat
c < G1 | G2 >
c     n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
c     n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
c     n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
c     n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
c     n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
c     n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
c
c      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
c      ppox = 1d99
c      do ig1  = 1, ng1
c      do ig2  = 1, ng2
c       nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
c        if( ppox(nx(1),nx(2),nx(3))==1d99 ) then
c          call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
c     i    nx(1:3), ! G2 -G1
c     o    ppox( nx(1),nx(2),nx(3)))
c       endif
c      enddo
c     enddo
      ppovl = 0d0
      do ig1 = 1,ng1
        do ig2 = 1,ng2
          nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
          if(sum(nx**2)==0)  ppovl(ig1,ig2) = 1d0
        enddo
      enddo
c      deallocate(ppox)
      end
