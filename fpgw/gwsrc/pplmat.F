      subroutine mkppmt(alat,plat,qlat, q, ng1,ngvec1, 
     i     rmax,  nbas,  bas, lmxa, lmxax, 
     o         ppmt)
C- ppmt value and slope at MT bounadaries for q+G plane waves.
Co  ppmt(1,lm,nbas): value for (lm,ibas)
Co  ppmt(2,lm,nbas):
Cr Oct2005
C------------------
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

#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4)::ll
#endif
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
      if(debug) stop 'test end ----------'
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

      subroutine pplmat2(s_lat, q, ng, ngvec, rmax, 
     &       bmat,kmaxx, nlmax, nbas, ndimh, bas, lmxa, lmxax,
     &       lh1,lh2,rsmh, eh,ntorb,ntorbx, kmax,
     &       plhd, plpkl, QpGcutHankel, 
     o       ppovl,ppsh)
c <P^q1_G1 | P^q2_G2 > matrix where P^q1_G1 denotes projected plane wave, zero within sphere.
c <P^q1_G1 | smooth Hankel > matrix
c
c Mar-2001 from pplmat
      implicit none
      integer(4) ::  nbas, ndimh, ng,kmaxx,nlmax,lmxax,ntorbx
     &       ,ig1,ig2,ibas,ih,la,ihoff,itorb,k,lma,lmh,lh,ntorb(nbas)
     &      , ngvec(3,ng),lmxa(nbas),kmax(nbas),
     &      lh1(ntorbx,nbas),lh2(ntorbx,nbas)
      real(8)    :: absqg1,absqg2,tripl,rmax(nbas),pi4,
     &              rsmh(ntorbx,nbas), eh(ntorbx,nbas),r2s,denom,gan
      complex(8) :: bmat(0:kmaxx, nlmax, nbas, ndimh)
      real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qg(3,ng)    ,facl
     & ,pi,alat,tpiba,cost,voltot,plegn,qg1(3),qg2(3),absqg(ng),qqq(3),
     &  bas(3,nbas), fkk(0:lmxax),fkj(0:lmxax),fjk(0:lmxax),fjj(0:lmxax)
      complex(8) :: img =(0d0,1d0),phase
c
      complex(8) ::  plhd (ng, ndimh)        ! integral plane \times head ???
     &             , plpkl(ng, nlmax, 0:kmaxx,nbas)  ! integral plane \times poly ???
      complex(8) :: ppovl(ng,ng), ppsh(ng,ndimh)
      real(8),allocatable::cy(:),yl(:),yl1(:),yl2(:)
      integer(4),allocatable :: ngvecx(:,:)
      complex(8),allocatable :: ppshx(:,:),ppox(:,:,:)
      complex(8),allocatable :: ppovlx(:,:)
      real(8) :: gv1(3),gv1x(3),QpGcutHankel,dummy,tpibaqlat(3,3)
      integer(4) :: ngmx,ig1x,nx(3),n1x,n2x,n3x,n1m,n2m,n3m
#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4)::ll
#endif
c-----------------------------------------------------------------------
      print *,' pplmat2:'
      pi  = 4d0*datan(1d0)
      pi4 = 16d0*datan(1d0)
      call u_lat_vecs(s_lat,alat,plat,qlat)
      tpiba  = 2*pi/alat
      tpibaqlat = tpiba *qlat
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
      call sylmnc(cy,lmxax)

ccc      goto 1201 ! if you want to use old version ---
c
c <P^q_G1 | P^q_G2 >
ccctest      call mkppovl2test(alat,plat,qlat,
      call mkppovl2(alat,plat,qlat, 
     &    ng, ngvec, !G1
     &    ng, ngvec, !G2
     &    nbas, rmax, bas,  
     o    ppovl)

ccccccccccccccccccccccccccccccccccccccc
c      write(3001,"('q=',3d16.8)") q
c      do ig1 = 1, ng
c      do ig2 = 1, ng
c        write(3001,"(2i3,'  ',2d16.8)") ig1,ig2,ppovl(ig1,ig2)
c      enddo
c      enddo
c      stop 'test1 end'
ccccccccccccccccccccccccccccccccccccccc

c Get ngvecx for the interstitial part of the smooth hankel function.
      call getgv2(alat,plat,qlat, q, QpGcutHankel, 1, ngmx,dummy)
      allocate( ngvecx(3,ngmx),ppshx(ngmx,ndimh) )
      call getgv2(alat,plat,qlat, q, QpGcutHankel, 2, ngmx, ngvecx)  ! for eigenfunction
      print *,  ' pplmat2: ngmx=',ngmx

c Expansion coefficients of |P^q_G1x> in the expansion .. |Smooth Hankel>
c ---> ppshx(igx1,ndimh) matrix
      ppshx = 0d0
      do ig1x = 1, ngmx
        qg1(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecx(1:3,ig1x)))
        absqg1   = sqrt(sum(qg1(1:3)**2))
        if(absqg1==0d0) then
          qqq = (/0d0,0d0,1d0/)
        else
          qqq = qg1/absqg1
        endif
        call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )
        ihoff = 0
        do ibas = 1,nbas
          if(ibas >=2) then
            ihoff  = ihoff  
     &      +  sum ( (  lh2(1:ntorb(ibas-1),ibas-1)+1)**2 
     &              - lh1(1:ntorb(ibas-1),ibas-1)**2   )
          endif
          ih = ihoff
          phase = exp( -img*sum( qg1*bas(:,ibas)*alat )  )
          do itorb = 1,ntorb(ibas)
            do lmh   = lh1(itorb,ibas)**2+1, (lh2(itorb,ibas)+1)**2
              lh = ll(lmh)
              ih = ih + 1
              denom = eh(itorb,ibas)-absqg1**2
              gan   = 1d0/4d0*rsmh(itorb,ibas)**2
              ppshx(ig1x,ih) = - pi4/denom
     &      * (-img)**lh * facl(absqg1,lh)
     &      * cy(lmh)* yl(lmh)* exp( gan*denom)
     &      * phase  /voltot
            enddo
          enddo
        enddo
      enddo
c
      allocate(ppovlx(ng,ngmx))
ccctest      call mkppovl2test(alat,plat,qlat,
      call mkppovl2(alat,plat,qlat, 
     &    ng,   ngvec,  !G1
     &    ngmx, ngvecx, !G1x
     &    nbas, rmax, bas,  
     o    ppovlx)
      call matm(ppovlx,ppshx,ppsh,ng,ngmx,ndimh)
      deallocate( ppovlx )
      if (allocated(cy)) deallocate(cy)
      if (allocated(yl)) deallocate(yl)
      if (allocated(ngvecx)) deallocate(ngvecx)
      if (allocated(ppshx)) deallocate(ppshx)
      return

c <P^q_G1 | Smooth Hankel > = ppsh matrix
c     n1x = maxval( ngvecx(1,:)) - minval( ngvec(1,:))
c     n1m = minval( ngvecx(1,:)) - maxval( ngvec(1,:))
c     n2x = maxval( ngvecx(2,:)) - minval( ngvec(2,:))
c     n2m = minval( ngvecx(2,:)) - maxval( ngvec(2,:))
c     n3x = maxval( ngvecx(3,:)) - minval( ngvec(3,:))
c     n3m = minval( ngvecx(3,:)) - maxval( ngvec(3,:))
c      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
c      ppox = 1d99
c   <G1|G1X>
c      do ig1  = 1, ng
c      do ig1x = 1, ngmx
c       nx(1:3) = ngvecx(1:3,ig1x) -ngvec(1:3,ig1)
c        if( ppox(nx(1),nx(2),nx(3))==1d99 ) then
c          call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
c     i    nx(1:3), ! G1x -G1
c     o    ppox( nx(1),nx(2),nx(3)))
c       endif
c      enddo
c     enddo
c      ppsh = 0d0
c      do ih   = 1, ndimh
c      do ig1  = 1, ng
c      do ig1x = 1, ngmx
c       nx(1:3) = ngvecx(1:3,ig1x) - ngvec(1:3,ig1)
c       ppovlx  = ppox( nx(1),nx(2),nx(3))
c       ppsh(ig1,ih) = ppsh(ig1,ih) + ppovlx * ppshx(ig1x,ih)
c      enddo
c     enddo
c     enddo
c      deallocate( ppox )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc




c=======================================================================
 1201 continue
c      print *,  ' pplmat: xxx voltot=',voltot
c <P^q_G1 | P^q_G2 > matrix
      do ig1 = 1,ng
        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvec(1:3,ig1)))
        absqg(ig1)  = sqrt(sum(qg(1:3,ig1)**2))
      enddo

c --- OLD codes from here ----------------------------------------------
      ppovl = 0d0
      do ig1 = 1,ng
        ppovl(ig1,ig1)= voltot
      enddo

cyy      allocate(yl1((lmxax+1)**2),yl2((lmxax+1)**2) )

      do ig1 =1,ng
        qg1(1:3) = qg(1:3,ig1)
        absqg1   = absqg(ig1)
cyy        call sylm(qg1/absqg1,yl1,lmxax,r2s) !spherical factor Y( q+G )
        do ig2 =1,ng
          qg2(1:3) = qg(1:3,ig2)
          absqg2   = absqg(ig2)
cyy        call sylm(qg2/absqg2,yl2,lmxax,r2s) !spherical factor Y( q+G )

          if( absqg1*absqg2 == 0d0) then
            cost = 1d0
          else
            cost = sum(qg1*qg2)/absqg1/absqg2
          endif

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

ccccccccccccccccccccccccccccccccccccccc
c      write(3002,"('q=',3d16.8)") q
c      do ig1 = 1, ng
c      do ig2 = 1, ng
c        write(3002,"(2i3,'  ',2d16.8)") ig1,ig2,ppovl(ig1,ig2)
c      enddo
c      enddo
c      stop 'test2 end'
ccccccccccccccccccccccccccccccccccccccc

cccccccccccccc
c      do ig1=1,ng-1
c       ig2=ig1+1
c        write(6,"(' old ppovl(ig1,ig2)=',i4,2d13.6)")ig1,ppovl(ig1,ig2)
c      enddo
c     stop "------testen--------d"
cccccccccccccc

 1111 continue
c <P^q_G1 | Smooth Hankel > = ppsh matrix
      ppsh=0d0
      do ig1 = 1, ng
        qg1(1:3) = qg(1:3,ig1)
        absqg1   = absqg(ig1)
        if(absqg1==0d0) then
          qqq = (/0d0,0d0,1d0/)
        else
          qqq = qg1/absqg1
        endif
        call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )

        ihoff  = 0
        do ibas = 1,nbas
          if(ibas >=2) then
            ihoff  = ihoff  
     &      +  sum ( (  lh2(1:ntorb(ibas-1),ibas-1)+1)**2 
     &                - lh1(1:ntorb(ibas-1),ibas-1)**2   )
          endif
c from head
          phase = exp( -img*sum( qg1*bas(:,ibas)*alat )  )
          ih = ihoff
          do itorb = 1,ntorb(ibas)
            do lmh   = lh1(itorb,ibas)**2+1, (lh2(itorb,ibas)+1)**2
              lh = ll(lmh)
              ih = ih + 1
c           print *,ibas,ip,itorb,lmh,ias,ih,lh
c           print *,lmh,lh
              denom = eh(itorb,ibas)-absqg1**2
              gan   = 1d0/4d0*rsmh(itorb,ibas)**2
              ppsh(ig1,ih) = ppsh(ig1,ih) - pi4/denom
c     &         * (-img*absqg1)**lh * cy(lmh)* yl(lmh)* exp( gan*denom)
     &         * (-img)**lh * facl(absqg1,lh)
     &         * cy(lmh)* yl(lmh)* exp( gan*denom)
     &         * phase
     &         - phase * plhd(ig1, ih)  ! integral of plane \times hd
            enddo
          enddo
c from tail
          do lma= 1, (lmxa(ibas)+1)**2
            do k  = 0, kmax(ibas)
              la  = ll(lma)
c           print *, ibas,ip,lma,k
              ppsh(ig1,1:ndimh) = ppsh(ig1,1:ndimh)  
     &        - phase * plpkl(ig1, lma, k,ibas)
     &                * bmat(k, lma, ibas, 1:ndimh) ! integral of plane \times hd
            enddo
          enddo
        enddo
      enddo
      if (allocated(cy)) deallocate(cy)
      if (allocated(yl)) deallocate(yl)
      if (allocated(ngvecx)) deallocate(ngvecx)
      if (allocated(ppshx)) deallocate(ppshx)
      end


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
!! ---------------------------------------------
!! ppovl = \int_{cell} \bfr IPW(G) d^3r ==
!! G = ggvec = g2 -g1 = 2pi/alat * \sum_n qlat(:,n)*ngvec(n)
!! 
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
          grmx  = absg* rmax(ibas)
          ppovl = ppovl
     &      - exp( img* sum(ggvec*bas(1:3,ibas))*alat )  
     &        * pi4/absg**3 * ( -grmx * cos(grmx) + sin(grmx) )
        endif
      enddo
      end

c---------------------
      integer function njudge000(ngvec)
      implicit none
      integer:: ngvec(3),ix,i,j,k
      if(sum(abs(ngvec))==0) then
         njudge000=0
         return
      endif
      do i=1,3
         if(i==1) then;  j=2; k=3; endif
         if(i==2) then;  j=3; k=1; endif
         if(i==3) then;  j=1; k=2; endif
         if(ngvec(i)/=0.and.ngvec(j)==0.and.ngvec(k)==0) then
           njudge000=i
           return
         endif
      enddo
      njudge000=-1
      end
      subroutine matgg2bfr(alat,bas,rmax,nbas,voltot, plat,qlat,
     i   ngvec,
     o   bfrppovl)
!! ---------------------------------------------
!! bfrppovl = \int_{cell} \bfr IPW(G) d^3r ==
!! G = ggvec = g2 -g1 = 2pi/alat * \sum_n qlat(:,n)*ngvec(n)
      implicit none
      integer(4) :: nbas,ibas,ngvec(3)
      real(8) :: ggvec(3), alat,bas(3,nbas),rmax(nbas),voltot,absg,grmx
     &           ,tpibaqlat(3,3)
      complex(8) :: bfrppovl(3)
      real(8),parameter :: pi =3.1415926535897932D0,pi4=pi*4d0
      complex(8) :: img =(0d0,1d0)
      real(8):: plat(3,3),qlat(3,3),dabsGdG(3),rcenter(3)
      integer:: ix,njudge000
c      print *,' matgg2bfr:'
      tpibaqlat =  2d0*pi/alat *qlat
      ggvec(1:3) =  matmul(tpibaqlat, ngvec(1:3))
      absg       =  sqrt(sum(ggvec(1:3)**2))
!! full space part
      ix= njudge000(ngvec)
      bfrppovl = 0d0
      if(ix==0) then
        bfrppovl = alat*voltot*0.5d0*(plat(:,1)+plat(:,2)+plat(:,3))
      elseif(ix>0) then
        bfrppovl = voltot*alat*plat(:,ix)/img/ngvec(ix)
      endif
!! subtract MT part
      if(absg/=0d0) dabsGdG = ggvec/absg
      do ibas = 1,nbas
        if(absg==0d0) then
          bfrppovl = bfrppovl - pi4*rmax(ibas)**3/3d0*alat*bas(1:3,ibas) !alat*bas is the center
        else
          grmx  = absg* rmax(ibas)
          bfrppovl = bfrppovl
     &      - img*bas(1:3,ibas)*alat * exp( img* sum(ggvec*bas(1:3,ibas))*alat )  
     &        * pi4/absg**3 * ( -grmx * cos(grmx) + sin(grmx) )
     &      -  exp( img* sum(ggvec*bas(1:3,ibas))*alat )  
     &        * pi4/absg**3 * rmax(ibas)*( grmx * sin(grmx) + 2*cos(grmx) )*dabsGdG(:)
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
#ifdef COMMONLL
      integer(4) ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif
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
