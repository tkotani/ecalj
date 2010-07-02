      subroutine shortn3(pin,noutmx, nout,nlatout)
C- Find shortest pin + nlatout \in  [ pin + any integer linear-combination of (rlat(:,1),rlat(:,2),rlat(:,3)) ].
Cr You first have to call shorn3_initialize, then you can call shortn3 for given rlat.
C
Cio init
Ci pin is on the rlat coodinate.
Cr rlat(3,i) i-th vertor for modulo
Ci rlatp(i,j)= sum( rlat(:,i)*rlat(:,j) )
Ci noutmax: upper limit of nlatout
Co
Co nlatout
C
Cr  length = sum_i sum_j pin(i)*rlatp(i,j)*pin(j)
Cr  pin+nlatout(:,ix), where ix=1:nout, is the shortest vectors. We may have multiple nlatout (# is nout).
c
c Takao think shorbz will be almost OK, but not perfect as an algolism.
c Takao think all shorbz should be replaced by shortn3 in future.
      implicit none
      integer:: nmax(3),nknknk,ik1,ik2,ik3,nout,nk,ik,i,j
      real(8):: rmax2,pin(3),eps=1d-8,rlat(3,3),xmax2(3),rr(3),rmin,nrmax(3)
      integer,allocatable:: nlat0(:,:)
      real(8),allocatable:: rnorm(:)
      integer:: noutmx
      integer:: nlatout(3,noutmx)
      real(8),save:: rlatp(3,3),xmx2(3)
c
      rmax2 = sum(pin*matmul(rlatp,pin)) + eps  ! eps is to make degeneracy safe.
      nrmax(:) =  sqrt(rmax2*xmx2(:))+abs(pin(:)) ! range of ix
      nmax =  nrmax

C we are looking for shortest vectors
      ik=0
      rmin=1d9
      nknknk= (2*nmax(1)+1)*(2*nmax(2)+1)*(2*nmax(3)+1)
      allocate( nlat0(3, nknknk), rnorm(nknknk) )
      do ik1=-nmax(1),nmax(1)
        do ik2=-nmax(2),nmax(2)
          do ik3=-nmax(3),nmax(3)
            ik=ik+1
            nlat0(:,ik) = (/ik1,ik2,ik3/)
            rr= pin + nlat0(:,ik)
            rnorm(ik) = sum(rr*matmul(rlatp,rr))
            if(rnorm(ik)<rmin) rmin=rnorm(ik)
          enddo
        enddo
      enddo
      nk=ik
      !print *,'nk rmin=',nk,rmin
      nout=0
      do ik=1,nk
        rr= pin + nlat0(:,ik)
        !print *,'ik rr   =',ik,rr
        !print *,'ik rnorm=',ik,rnorm(ik)
        if(rnorm(ik)<rmin+eps) then
          nout=nout+1
          if(nout>noutmx) stop 'shortn3: enlarge noutmx'
          nlatout(:,nout)=nlat0(:,ik)
          !print *,'ik nlat0',nlat0(:,ik)
        endif
      enddo
      !write(6,"('pin=',3f8.3,' nmax=',3i4,' nout=',i3)")pin, nmax(1:3),nout
      deallocate(rnorm,nlat0)
      return

      entry shortn3_initialize(rlat)
C- Set translation vactors rlat(:,i),i=1,3
      do i=1,3
        do j=1,3
          rlatp(i,j) = sum(rlat(:,i)*rlat(:,j))
        enddo
      enddo
      call ellipsoidxmax(rlatp,xmx2)
      end

      subroutine ellipsoidxmax(nn, xmx2)
C- Maximum value for x_i for ellipsoid
C  Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.
Ci nn(3,3)
Co xmx2(i)  Maximum of x_i**2
      implicit none
      real(8),target:: nn(3,3)
      real(8):: v2(2),ainv(2,2), rmax2, xmx2(3),det,fac,nv2(2)
      real(8),pointer::n11,n12,n13,n21,n22,n23,n31,n32,n33
      n11=>nn(1,1)
      n12=>nn(1,2)
      n13=>nn(1,3)
      n21=>nn(2,1)
      n22=>nn(2,2)
      n23=>nn(2,3)
      n31=>nn(3,1)
      n32=>nn(3,2)
      n33=>nn(3,3)
      !print *,'nn=',nn
c
      det= n22*n33-n23*n32
      ainv(1,1)=  n33/det
      ainv(2,2)=  n22/det
      ainv(1,2)= -n23/det
      ainv(2,1)= -n32/det
      nv2  = (/n12,n13/)
      fac = n11-sum(nv2 *matmul(ainv,nv2))
      !print *,'ainv=',ainv
      !print *,'ainv*nv=',matmul(ainv,nv2)
      xmx2(1) = 1d0/fac
c
      det= n33*n11-n31*n13
      ainv(1,1)=  n11/det
      ainv(2,2)=  n33/det
      ainv(1,2)= -n31/det
      ainv(2,1)= -n13/det
      nv2  = (/n23,n21/)
      fac = n22-sum(nv2 *matmul(ainv,nv2))
      xmx2(2) = 1d0/fac
c
      det= n11*n22-n12*n21
      ainv(1,1)=  n22/det
      ainv(2,2)=  n11/det
      ainv(1,2)= -n12/det
      ainv(2,1)= -n21/det
      nv2  = (/n31,n32/)
      fac = n33-sum(nv2 *matmul(ainv,nv2))
      xmx2(3) = 1d0/fac
      end



      subroutine shorbz(p,pout,plat,qlat)
C- Shortens vector to equivalent in first Brillouin zone.
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat,qlat lattice vectors and inverse
Ci   p         vector to shorten
Co Outputs:
Co   pout      shortened p.   pout and p can point to the same address.
Cr Remarks
Cr   Switch around plat,qlat to shorten reciprocal space vectors.
Cr   Jan 1997 Adapted from shorps to fix bug:  Example:
Cr   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
Cr   p= 0.0 -0.5 -1.26384
Cr   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
C ----------------------------------------------------------------
C     implicit none
      double precision p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
      double precision tol
      parameter (tol=-1d-10)
      integer i,j,m,j2min,j3min,j1,j2,j3

C --- Reduce to unit cell centered at origin ---
      do  1  i = 1, 3
C ... x is projection of pin along plat(i), with multiples of p removed
        x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
        xx = idnint(x0)
        x(i) = x0-xx
    1 continue
C ... pout is x rotated back to Cartesian coordinates
      do  2  m = 1, 3
        pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)
    2 continue

C --- Try shortening by adding +/- basis vectors ---
   15 continue
      do  10  j1 =  0, 1
        j2min = -1
        if (j1 .eq. 0) j2min = 0
        do  10  j2 = j2min, 1
          j3min = -1
          if (j1 .eq. 0 .and. j2 .eq. 0) j3min = 0
          do  10  j3 = j3min, 1

C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
c takao think this alglorism cause problems (=wrong) for very anisotropic cases.
c
            do  17  i = 1, 3
              x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
   17       continue
            a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
            ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
            j = 0
            if (a2 + 2*ap .lt. tol) j = 1
            if (a2 - 2*ap .lt. tol) j = -1
            if (j .ne. 0) then
              pout(1) = pout(1) + j*x(1)
              pout(2) = pout(2) + j*x(2)
              pout(3) = pout(3) + j*x(3)
              goto 15
            endif
   10 continue

      end
C      subroutine fmain
C      implicit none
C      double precision plat(9),qlat(9),p(3),p1(3),xx
C      integer mode(3)
C
C      data plat /-0.5d0,0.5d0,1.7517d0,
C     .            0.5d0,-.5d0,1.7517d0,
C     .            0.5d0,0.5d0,-1.7517d0/
C
C      integer w(10000)
C      common /w/ w
C
C      data p /0.0d0,-0.5d0,-1.2638400000000001d0/
C
C      call wkinit(10000)
C
CC ... qlat = (plat^-1)^T so that qlat^T . plat = 1
C      call mkqlat(plat,qlat,xx)
C
C      call shorbz(p,p1,plat,qlat)
C      call prmx('p1 from shorbz',p1,1,1,3)
C
C      mode(1) = 2
C      mode(2) = 2
C      mode(3) = 3
C      call shorps(1,plat,mode,p,p1)
C      call prmx('p1 from shorps',p1,1,1,3)
C      end


      subroutine shorbzm(flags,p,pout,plat,qlat)
C- Shortens vector to equivalent in first Brillouin zone.
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat,qlat lattice vectors and inverse
Ci   p         vector to shorten
Co Outputs:
Co   pout      shortened p.   pout and p can point to the same address.
Cr Remarks
Cr   Switch around plat,qlat to shorten reciprocal space vectors.
Cr   Jan 1997 Adapted from shorps to fix bug:  Example:
Cr   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
Cr   p= 0.0 -0.5 -1.26384
Cr   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
C ----------------------------------------------------------------
C     implicit none
      integer,intent(in)::flags(3)
      double precision p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
      double precision tol
      parameter (tol=-1d-10)
      integer i,j,m,j2min,j3min,j1,j2,j3
      real(8),parameter:: eps=1.0d-10

C --- Reduce to unit cell centered at origin ---
      do  1  i = 1, 3
C ... x is projection of pin along plat(i), with multiples of p removed
        x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
        if (flags(i).ne.0) then
          xx = idnint(x0-eps)
        else
          xx = idnint(x0)
        endif
        x(i) = x0-xx
    1 continue
C ... pout is x rotated back to Cartesian coordinates
      do  2  m = 1, 3
        pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)
    2 continue

C --- Try shortening by adding +/- basis vectors ---
   15 continue
      do  10  j1 =  0, 1
        j2min = -1
        if (j1 .eq. 0) j2min = 0
        do  10  j2 = j2min, 1
          j3min = -1
          if (j1 .eq. 0 .and. j2 .eq. 0) j3min = 0
          do  10  j3 = j3min, 1

C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
c takao think this alglorism cause problems (=wrong) for very anisotropic cases.
c
            do  17  i = 1, 3
              x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
   17       continue
            a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
            ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
            j = 0
            if (a2 + 2*ap .lt. tol) j = 1
            if (a2 - 2*ap .lt. tol) j = -1
            if (j .ne. 0) then
              pout(1) = pout(1) + j*x(1)
              pout(2) = pout(2) + j*x(2)
              pout(3) = pout(3) + j*x(3)
              goto 15
            endif
   10 continue

      end

