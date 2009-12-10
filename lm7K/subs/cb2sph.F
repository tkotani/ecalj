      subroutine cb2sph(lhaveu,u,nu,cub,nlma,nlmb,wk,nd,sph)
C- Rotate one block of strux from cubic to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   cub,nlma,nlmb: strux and dimensions in cubic harmonics
Ci   lhaveu: T, rotation matrix u already made (F, cb2sph makes u)
Ci   wk:     work array of length nlma*nlmb*2
Ci   nd:     dimensioning parameter for sph
Cio Inputs/Outputs
Ci   u,nu:   Rotation matrix and dimension to do u S_real u+ = S_sph
Ci           If lhaveu=T, u is input
Ci           If lhaveu=F, u is generated and output
Co Outputs
Co   sph:    strx in spherical harmonics
Co   lhaveu: set to T
Cr Remarks
Cr   This routine makes and uses a rotation matrix u
Cr   that does the following rotation.
Cr   Let R_lm = real harmonic; Y_lm = spherical harmonic.  Then
Cr     Y_l,m = 1/sqrt(2)      (R_l,-m + i * R_l,m),  m<0
Cr     Y_l,m = R_l,m,                                m=0
Cr     Y_l,m = (-1)^m/sqrt(2) (R_l,-m - i * R_l,m),  m>0
Cr
Cr   In particular
Cr    l  m        R_l,m                 Y_l,m
Cr    1  -1   sqrt(3/4/pi)*y            sqrt(3/4/pi/2)(x + i*y)
Cr    1   1   sqrt(3/4/pi)*x            sqrt(3/4/pi/2)(-x + i*y)
Cr    2  -2   2*sqrt(15/16/pi)*x*y      sqrt(15/32/pi)*(x + i*y)^2
Cr    2   2   sqrt(15/16/pi)*(x*x-y*y)
Cr
Cr   The Y_l,m here are the same as those in Jackson, except that
Cr       Y_l,m = Y_l,-m (Jackson)
Cr
Cr   Other remarks:
Cr   lhaveu may be T on subsequent calls if u is preserved.
Cr   If nlma or nlmb=0, only u is made; cub,wk,nd,sph are unused.
Cu Updates
Cu   28 Apr 98 Changed u->u*; definition compatible with e.g. Jackson
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      logical lhaveu
      integer nu,nlma,nlmb,nd
      double precision u(nu,nu,2),cub(nlma,nlmb),wk(nlma,nlmb,2),
     .  sph(nd,nd,2)
C Local variables
      integer na,nb,nl,l,m,la,la0,nla,lb,lb0,nlb,ll,lc
      double precision a,sm

C --- Make rotation matrix u ---
      nl = ll(nu)+1
      if (.not. lhaveu) then
        a = 1/dsqrt(2d0)
        call dpzero(u,nu*nu*2)
        do  10  l = 0, nl-1
          lc = (l+1)**2-l
          u(lc,lc,1) = 1
          sm = -1
          do  12  m = 1, l
            u(lc+m,lc+m,1) = sm*a
            u(lc-m,lc+m,1) = a
            u(lc-m,lc-m,2) = a
            u(lc+m,lc-m,2) = -sm*a
            sm = -sm
   12     continue
   10   continue
        lhaveu = .true.
      endif

C      call yprm('u in cb2sph',2,u,nu*nu,nu,nu,nu)
      if (nlma .le. 0 .or. nlmb .le. 0) return
      na = ll(nlma)+1
      nb = ll(nlmb)+1

C --- make cub u+ by l-blocks ---
      do  20  lb = 0, nb-1
        lb0 = lb**2+1
        nlb = 2*lb+1
        call dgemm('N','T',nlma,nlb,nlb,1d0,cub(1,lb0),nlma,
     .    u(lb0,lb0,1),nu,0d0,wk(1,lb0,1),nlma)
        call dgemm('N','T',nlma,nlb,nlb,-1d0,cub(1,lb0),nlma,
     .    u(lb0,lb0,2),nu,0d0,wk(1,lb0,2),nlma)
   20 continue

C --- make u (cub u+) by la-blocks ---
      do  30  la = 0, na-1
        la0 = la**2+1
        nla = 2*la+1
        call yygemm('N','N',nla,nlmb,nla,1d0,u(la0,la0,1),u(la0,la0,2),
     .    nu,wk(la0,1,1),wk(la0,1,2),nlma,0d0,
     .    sph(la0,1,1),sph(la0,1,2),nd)
   30 continue

C     call yprm('sph in cb2sph',2,sph,nd*nd,nd,nlma,nlmb)

      end
C      subroutine fmain
C      implicit none
C      integer i,l,m
C      double precision u(16,16,2),xx,theta,phi,ylm(16,2),ylms(16,2),sr2
C      double precision scub(16,16), swk(16,2,16), ssph(16,2,16),x,y,z,rr
C      logical lf
C      integer w(1000)
C      common /w/ w
C
C      call pshpr(0)
C      call wkinit(1000)
C      call poppr
C
C      call dpzero(scub,16*16)
C      lf = .false.
C      call cb2sph(lf,u,16,scub,0,0,swk,1,ssph)
C
C      sr2 = 1/dsqrt(2d0)
C      theta = 3
C      phi = 4
C      print *, 'theta,phi=?'
C      read(*,*) theta,phi
C      x = cos(theta)*sin(phi)
C      y = sin(theta)*sin(phi)
C      z = cos(phi)
C      call dpzero(ylm,2)
C      call ropyln(1,x,y,z,3,1,ylm,rr)
C
C      call ygemm('N','N',
C     .  16,1,16,1d0,u,16*16,16,ylm,16,16,0d0,ylms,16,16)
C      print 332, x,y,z,
C     .  exp((0d0,1d0)*phi), exp((0d0,2d0)*phi), exp((0d0,3d0)*phi)
C  332 format(' x,y,z=',3f12.6/
C     .       '  exp(  i phi)=',2f12.6/
C     .       '  exp(2 i phi)=',2f12.6/
C     .       '  exp(3 i phi)=',2f12.6)
C      print *, ' l  m  L       Y_real        Re Y_L      Im Y_L'
C
CC     for this table check that these cols are zero:
CC     mc dat -e2 x5-x7 x6-x8
C      i = 0
C      do  l = 0, 3
C      do  m = -l, l
C        i = i+1
C        if (m .lt. 0) then
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
C     .      sr2*(ylm(i-2*m,1)),sr2*(ylm(i,1))
C  334     format(3i3,2x,f12.6,2x,2f12.6,2x,2f12.6)
C        elseif (m .gt. 0) then
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
C     .      sr2*(-1)**m*ylm(i,1), sr2*(-1)**m*(-1)*ylm(i-2*m,1)
C        else
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2)
C        endif
C
C      enddo
C      enddo
C
C
C      end

