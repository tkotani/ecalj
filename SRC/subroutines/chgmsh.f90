subroutine chgmsh(iopt,plat,n,m1,m2,m3,l1,l2,l3,f0, &
     n1,n2,n3,k1,k2,k3,f)
  !- Retabulate a function on a different real-space mesh
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   iopt  :0 Use default (smaller of iopt=1,2)
  !i         :1 use Nyquist cutoff
  !i         :2 use cutoff for largest sphere in BZ
  !i         :3 exactly double the mesh
  !i   plat  :primitive lattice vectors, in units of alat
  !i   n     :number of functions to change
  !i   m1..m3:number of divisions for the original mesh
  !i   l1..l3:dimensions of f0
  !i   f0    :function on m1,m2,m3 mesh
  !i   k1..k3:dimensions of f
  ! o Inputs/Outputs
  ! o  n1..n3:Destination mesh
  ! o        :n1,n2,n3 are input, unless iopt eq 3
  ! o        :If iopt eq 3  n1..n3 are output as twice (m1..m3)
  !o Outputs
  !o    f0        FT of f0 is returned in f0
  !o    f         function on n1,n2,n3 mesh
  !r Remarks
  !r   f0 and f may occupy the same address space
  !u Updates
  !u   25 Jun 00 added argument n
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: iopt,n,m1,m2,m3,n1,n2,n3,l1,l2,l3,k1,k2,k3
  double precision :: plat(3,3)
  double complex f0(l1,l2,l3,n),f(k1,k2,k3,n)
  integer :: ng1,ng2,ngmx,iprint,lgunit
  complex(8) ,allocatable :: cv1(:)
  double precision :: gmax,gmax1,gmax2,tau(3)
  external dpzero,fftz3,gvctof,gvgetf, &
       gvlist,gvputf,pchms2,pchmsh,poppr,pshpr,rx

  !      integer w(1)
  !      common /w/ w
  real(8),allocatable:: gv1(:),gv2(:)
  integer,allocatable:: kv1(:),kv2(:)
  integer:: wdummy
  !      logical:: fullmesh !not yet for full mesh mode.
  print *,' chgmsh: iopt=',iopt

  tau=0d0 !call dpzero(tau,3)
  if (iopt == 3) then
     n1 = 2*m1
     n2 = 2*m2
     n3 = 2*m3
     if (m3 == 1) n3 = 1
  endif
  if (iprint() >= 30) &
       write(6,"('CHGMSH: remake from ',3i4,' to ',3i4, ' iopt=',i3)") m1,m2,m3,n1,n2,n3,iopt

  if (iopt == 3) goto 100

  ! ... Lists of vectors for old and target mesh
  call pshpr(iprint()-30)
  call gvctof(iopt,1d0,plat,tau,m1,m2,m3,gmax1,ng1)
  call gvctof(iopt,1d0,plat,tau,n1,n2,n3,gmax2,ng2)
  gmax = dmin1(gmax1,gmax2)
  ngmx = min0(ng1,ng2)
  !      if(fullmesh()) then
  !         print *,'full mesh mode'
  !         gmax=1d10
  !         ngmx=min0(m1*m2*m3,n1*n2*n3)
  !      endif
  allocate(gv1(ngmx*3), gv2(ngmx*3),kv1(ngmx*3),kv2(ngmx*3))
  !      call gvlist(1d0,plat,wdummy,m1,m2,m3,gmax,8,ngmx,ng1,kv1,gv1,wdummy,wdummy)
  !      call gvlist(1d0,plat,wdummy,n1,n2,n3,gmax,8,ngmx,ng2,kv2,gv2,wdummy,wdummy)
  call gvlst2(1d0,plat,wdummy,m1,m2,m3,0d0,gmax,0,8,ngmx,ng1,kv1,gv1,wdummy,wdummy)
  call gvlst2(1d0,plat,wdummy,n1,n2,n3,0d0,gmax,0,8,ngmx,ng2,kv2,gv2,wdummy,wdummy)
  if (ng1 /= ng2) call rx('chgmsh: ng1 /= ng2')

  !      call gvmtch(ng1,gv1,kv1,ng2,gv2,kv2)  <--- bug. I had to pass ogv1 and so on in original code. Bug fix on 12Jun2010.
  call pgvmat2 (ng1,gv1, ng2,gv2,kv2 )

  call poppr
  allocate(cv1(ng1*n))!,cv2(ng1*n))
  cv1(:)=0d0
  !      cv2(:)=0d0
  call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
  call gvgetf ( ng1 , n , kv1, l1 , l2 , l3 , f0, cv1 )
  call gvputf ( ng1 , n , kv2, k1 , k2 , k3 , cv1, f )
  call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)
  deallocate(cv1,gv1,gv2,kv1,kv2)
  return

100 continue
  call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
  if (n3 /= 1) then
     call pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
  else
     call pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
  endif
  call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)

end subroutine chgmsh
subroutine pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
  !- Copies Fourier transform on one mesh to a doubled mesh
  !     implicit none
  integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,n
  double complex f0(l1,l2,l3,n),f(k1,k2,k3,n)
  integer :: i,i1,i2,i3,i1m,i2m,i3m,j1m,j2m,j3m
  !     call zprm3('initial mesh',0,f0,m1,m2,m3)
  f=0d0 !call dpzero(f,2*k1*k2*k3*n)
  do  101  i = 1, n
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


30         enddo
20      enddo
10   enddo
101 enddo
end subroutine pchmsh

subroutine pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
  !- 2D analog of pchmsh
  !     implicit none
  integer :: m1,m2,l1,l2,k1,k2,n
  double complex f0(l1,l2,n),f(k1,k2,n)
  integer :: i,i1,i2,i1m,i2m,j1m,j2m

  !     call zprm3('initial mesh',0,f0,m1,1,m2)

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


30      enddo
20   enddo
10 enddo
end subroutine pchms2

