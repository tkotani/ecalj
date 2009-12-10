      subroutine chgmsh(iopt,plat,n,m1,m2,m3,l1,l2,l3,f0,
     .                              n1,n2,n3,k1,k2,k3,f)
C- Retabulate a function on a different real-space mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :0 Use default (smaller of iopt=1,2)
Ci         :1 use Nyquist cutoff
Ci         :2 use cutoff for largest sphere in BZ
Ci         :3 exactly double the mesh
Ci   plat  :primitive lattice vectors, in units of alat
Ci   n     :number of functions to change
Ci   m1..m3:number of divisions for the original mesh
Ci   l1..l3:dimensions of f0
Ci   f0    :function on m1,m2,m3 mesh
Ci   k1..k3:dimensions of f
Cio Inputs/Outputs
Cio  n1..n3:Destination mesh
Cio        :n1,n2,n3 are input, unless iopt eq 3
Cio        :If iopt eq 3  n1..n3 are output as twice (m1..m3)
Co Outputs
Co    f0        FT of f0 is returned in f0
Co    f         function on n1,n2,n3 mesh
Cr Remarks
Cr   f0 and f may occupy the same address space
Cu Updates
Cu   25 Jun 00 added argument n
C ----------------------------------------------------------------------
C     implicit none
      integer iopt,n,m1,m2,m3,n1,n2,n3,l1,l2,l3,k1,k2,k3
      double precision plat(3,3)
      double complex f0(l1,l2,l3,n),f(k1,k2,k3,n)
C Local variables
      integer ng1,ng2,ngmx,iprint,lgunit
      integer ocv1,ocv2,ogv1,ogv2,okv1,okv2
      double precision gmax,gmax1,gmax2,tau(3)
C ... External calls
      external defcc,defi,defrr,dpzero,fftz3,gvctof,gvgetf,
     .         gvlist,gvmtch,gvputf,pchms2,pchmsh,poppr,pshpr,rlse,rx

      integer w(1)
      common /w/ w

      call dpzero(tau,3)
      if (iopt .eq. 3) then
        n1 = 2*m1
        n2 = 2*m2
        n3 = 2*m3
        if (m3 .eq. 1) n3 = 1
      endif
#if AWRITE
      if (iprint() .ge. 30) then
        call awrit7('%N CHGMSH: remake f (%i * %i * %i)  to '//
     .    'f (%i * %i * %i),  opt=%i',' ',80,
     .    lgunit(1),m1,m2,m3,n1,n2,n3,iopt)
      endif
#endif

      if (iopt .eq. 3) goto 100

C ... Lists of vectors for old and target mesh
      call pshpr(iprint()-30)
      call gvctof(iopt,1d0,plat,tau,m1,m2,m3,gmax1,ng1)
      call gvctof(iopt,1d0,plat,tau,n1,n2,n3,gmax2,ng2)
      gmax = dmin1(gmax1,gmax2)
      ngmx = min0(ng1,ng2)

C ... Create the two lists and align them
      call defrr(ogv1, ngmx*3)
      call defrr(ogv2, ngmx*3)
      call defi (okv1, ngmx*3)
      call defi (okv2, ngmx*3)
C      call ogvlst(1d0,plat,tau,m1,m2,m3,gmax,ngmx,ng1,w(ogv1),w(okv1))
C      call ogvlst(1d0,plat,tau,n1,n2,n3,gmax,ngmx,ng2,w(ogv2),w(okv2))
      call gvlist(1d0,plat,w,m1,m2,m3,gmax,8,ngmx,ng1,w(okv1),w(ogv1),
     .  w,w)
      call gvlist(1d0,plat,w,n1,n2,n3,gmax,8,ngmx,ng2,w(okv2),w(ogv2),
     .  w,w)

      if (ng1.ne.ng2) call rx('chgmsh: ng1.ne.ng2')
      call gvmtch(ng1,ogv1,okv1,ng2,ogv2,okv2)
      call poppr

C ... Copy FT to first list, copy back from second list
      call defcc(ocv1,  -ng1*n)
      call defcc(ocv2,  -ng1*n)
      call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
      call gvgetf(ng1,n,w(okv1),l1,l2,l3,f0,w(ocv1))
      call gvputf(ng1,n,w(okv2),k1,k2,k3,w(ocv1),f)
      call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)

      call rlse(ogv1)
      return

  100 continue
      call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
      if (n3 .ne. 1) then
        call pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
      else
        call pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
      endif
      call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)

      end
      subroutine pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
C- Copies Fourier transform on one mesh to a doubled mesh
C     implicit none
      integer m1,m2,m3,l1,l2,l3,k1,k2,k3,n
      double complex f0(l1,l2,l3,n),f(k1,k2,k3,n)
      integer i,i1,i2,i3,i1m,i2m,i3m,j1m,j2m,j3m

C     call zprm3('initial mesh',0,f0,m1,m2,m3)

      call dpzero(f,2*k1*k2*k3*n)
      do  10  i = 1, n
      do  10  i3 = 1, (m3+1)/2
      i3m = m3+1-i3
      j3m = 2*m3+1-i3

        do  20  i2 = 1, (m2+1)/2
        i2m = m2+1-i2
        j2m = 2*m2+1-i2

          do  30  i1 = 1, (m1+1)/2
          i1m = m1+1-i1
          j1m = 2*m1+1-i1

            f(i1,i2,i3,i)   = f0(i1,i2,i3,i)
            f(i1,i2,j3m,i)  = f0(i1,i2,i3m,i)
            f(i1,j2m,i3,i)  = f0(i1,i2m,i3,i)
            f(i1,j2m,j3m,i) = f0(i1,i2m,i3m,i)
            f(j1m,i2,i3,i)  = f0(i1m,i2,i3,i)
            f(j1m,i2,j3m,i) = f0(i1m,i2,i3m,i)
            f(j1m,j2m,i3,i) = f0(i1m,i2m,i3,i)
            f(j1m,j2m,j3m,i)= f0(i1m,i2m,i3m,i)


   30     continue
   20   continue
   10 continue

C     call zprm3('final mesh',0,f,2*m1,2*m2,2*m3)
      end
      subroutine pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
C- 2D analog of pchmsh
C     implicit none
      integer m1,m2,l1,l2,k1,k2,n
      double complex f0(l1,l2,n),f(k1,k2,n)
      integer i,i1,i2,i1m,i2m,j1m,j2m

C     call zprm3('initial mesh',0,f0,m1,1,m2)

      call dpzero(f,2*k1*k2*n)
      do  10  i = 1, n
        do  20  i2 = 1, (m2+1)/2
        i2m = m2+1-i2
        j2m = 2*m2+1-i2

          do  30  i1 = 1, (m1+1)/2
          i1m = m1+1-i1
          j1m = 2*m1+1-i1

            f(i1,i2,i)   = f0(i1,i2,i)
            f(i1,i2,i)   = f0(i1,i2,i)
            f(i1,j2m,i)  = f0(i1,i2m,i)
            f(i1,j2m,i)  = f0(i1,i2m,i)
            f(j1m,i2,i)  = f0(i1m,i2,i)
            f(j1m,i2,i)  = f0(i1m,i2,i)
            f(j1m,j2m,i) = f0(i1m,i2m,i)
            f(j1m,j2m,i) = f0(i1m,i2m,i)


   30     continue
   20   continue
   10 continue

C     call zprm3('final mesh',0,f,2*m1,1,2*m2)
      end

