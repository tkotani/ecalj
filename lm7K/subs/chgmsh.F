Cdelw1 var ['ogv1', '0', 'ogv1', 'ngmx*3', 'real(8)', 'rv_w_', 'chgmsh']
Cdelw1 var ['ogv2', '0', 'ogv2', 'ngmx*3', 'real(8)', 'rv_w_', 'chgmsh']
Cdelw1 var ['okv1', '0', 'okv1', 'ngmx*3', 'integer', 'iv_w_', 'chgmsh']
Cdelw1 var ['okv2', '0', 'okv2', 'ngmx*3', 'integer', 'iv_w_', 'chgmsh']
Cdelw1 var ['ocv1', '0', 'ocv1', '-ng1*n', 'complex(8)', 'zv_w_', 'chgmsh']
Cdelw1 var ['ocv2', '0', 'ocv2', '-ng1*n', 'complex(8)', 'zv_w_', 'chgmsh']
Cdelw1 undel [['wref'], 'ogv1', 'chgmsh']
Cdelw1 undel [['wref'], 'ogv2', 'chgmsh']
Cdelw1 undel [['wref'], 'okv1', 'chgmsh']
Cdelw1 undel [['wref'], 'okv2', 'chgmsh']
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
Cdelw1 do not change  [['wref'], 'ogv1', 'chgmsh']
Cdelw1 do not change  [['wref'], 'ogv2', 'chgmsh']
Cdelw1 do not change  [['wref'], 'okv1', 'chgmsh']
Cdelw1 do not change  [['wref'], 'okv2', 'chgmsh']
Cdelw1       integer ocv1,ocv2,ogv1,ogv2,okv1,okv2
       integer:: ogv1 , ogv2 , okv1 , okv2 
       complex(8) ,allocatable :: zv_w_ocv1(:)
       complex(8) ,allocatable :: zv_w_ocv2(:)

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
Cdelw1 do not change ogv1 because of ['wref']
      call defrr(ogv1, ngmx*3)
Cdelw1 do not change ogv2 because of ['wref']
      call defrr(ogv2, ngmx*3)
Cdelw1 do not change okv1 because of ['wref']
      call defi (okv1, ngmx*3)
Cdelw1 do not change okv2 because of ['wref']
      call defi (okv2, ngmx*3)
C      call ogvlst(1d0,plat,tau,m1,m2,m3,gmax,ngmx,ng1,w(ogv1),w(okv1))
C      call ogvlst(1d0,plat,tau,n1,n2,n3,gmax,ngmx,ng2,w(ogv2),w(okv2))
Cdelw1 do not change ogv1 because of ['wref']
Cdelw1 do not change okv1 because of ['wref']
      call gvlist(1d0,plat,w,m1,m2,m3,gmax,8,ngmx,ng1,w(okv1),w(ogv1),
     .  w,w)
Cdelw1 do not change ogv2 because of ['wref']
Cdelw1 do not change okv2 because of ['wref']
      call gvlist(1d0,plat,w,n1,n2,n3,gmax,8,ngmx,ng2,w(okv2),w(ogv2),
     .  w,w)

      if (ng1.ne.ng2) call rx('chgmsh: ng1.ne.ng2')
Cdelw1 do not change ogv1 because of ['wref']
Cdelw1 do not change ogv2 because of ['wref']
Cdelw1 do not change okv1 because of ['wref']
Cdelw1 do not change okv2 because of ['wref']
      call gvmtch(ng1,ogv1,okv1,ng2,ogv2,okv2)
      call poppr

C ... Copy FT to first list, copy back from second list
Cdelw1       call defcc(ocv1,  -ng1*n)
       allocate(zv_w_ocv1(abs(-ng1*n)))
       if (-ng1*n<0) zv_w_ocv1(:)=0.0d0

Cdelw1       call defcc(ocv2,  -ng1*n)
       allocate(zv_w_ocv2(abs(-ng1*n)))
       if (-ng1*n<0) zv_w_ocv2(:)=0.0d0

      call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
Cdelw1 do not change okv1 because of ['wref']
Cdelw1 do not change okv1 because of ['wref']
Cdelw1       call gvgetf(ng1,n,w(okv1),l1,l2,l3,f0,w(ocv1))
       call gvgetf ( ng1 , n , w ( okv1 ) , l1 , l2 , l3 , f0 , zv_w_ocv1 
     . ) 

Cdelw1 do not change okv2 because of ['wref']
Cdelw1 do not change okv2 because of ['wref']
Cdelw1       call gvputf(ng1,n,w(okv2),k1,k2,k3,w(ocv1),f)
       call gvputf ( ng1 , n , w ( okv2 ) , k1 , k2 , k3 , zv_w_ocv1 
     . , f ) 

      call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)

Cdelw1 rlse name= ogv1 old_list= ogv1 ogv2 okv1 okv2 ocv1 ocv2 
Cdelw1 rlse name= ogv1 new_list= (None)
Cdelw1 not deallocate okv2 because of [wref]
Cdelw1 not deallocate okv1 because of [wref]
Cdelw1 not deallocate ogv2 because of [wref]
Cdelw1 not deallocate ogv1 because of [wref]
Cdelw1       call rlse(ogv1)
       if (allocated(zv_w_ocv2)) deallocate(zv_w_ocv2)
       if (allocated(zv_w_ocv1)) deallocate(zv_w_ocv1)
       call rlse(ogv1)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ogv1 ogv2 okv1 okv2]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

  100 continue
      call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
      if (n3 .ne. 1) then
        call pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
      else
        call pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
      endif
      call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ogv1 ogv2 okv1 okv2]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

