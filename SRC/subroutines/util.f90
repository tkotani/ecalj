!> A Collection of utility routines and small math routines.
module m_lgunit ! file handles for standard output log, and mpilog
  !  stdo: file handle for standard output
  !  stdl: handle for log
  !  stml: mpilog
  public:: m_lgunit_init,ipr
  integer,public :: stdl,stdo=6,stml
  logical:: ipr
  private
contains
  subroutine M_lgunit_init()
!    use m_mpi,only: mpi__root
    logical:: cmdopt0
    stdo= lgunit(1)
    stdl= lgunit(2)
!    if(cmdopt0('--mlog')) stml=lgunit(3)
  end subroutine M_lgunit_init
  integer function lgunit(i)
    ! Returns stdout for i=1, log for i=2, mlog for i=3 (MPI logfile)
    use m_ext,only: sname
    implicit none
    character(10):: i2char
    character*100 ext
    integer :: i, fopn, i1mach, fhndl,ierr,procid
    integer,save:: lgunit1=0,lgunit2=0,lgunit3=0
    include 'mpif.h'
    lgunit = 6
    if (i .eq. 1) return
    if (i .eq. 2) then
       if(lgunit2==0) then
          open(newunit=lgunit2,file='log.'//trim(sname),position='append')
       endif
       lgunit = lgunit2
    elseif (i .eq. 3) then
      return
!       if(lgunit3==0) then
!        call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
!        open(newunit=lgunit3,file='mlog.'//trim(sname)//'_'//trim(i2char(procid)))
!       endif
!       lgunit = lgunit3
    endif
  end function lgunit
end module m_lgunit
real(8) function rydberg()
  rydberg=13.6058d0
END function rydberg
real(8) function bohr()
  bohr=.529177d0
END function bohr
pure subroutine rangedq(qin, qout) ! qout is in [-0.5d0,0.5d0)
  implicit none
  intent(in)::     qin
  intent(out)::         qout
  integer :: ix
  real(8):: qin(3),qout(3)
  real(8),parameter::tol=1d-6
  qout= qin-nint(qin)
  qout= merge(-0.5d0,qout,mask=qout>.5d0-tol)
end subroutine rangedq
character(8) function xt(num)
  integer :: num
  if(num==0) then
     xt=''
     return
  endif
  xt = char(48+mod(num,10))
  if(num>9)     xt = char(48+mod(num/10,10))//xt
  if(num>99)    xt = char(48+mod(num/100,10))//xt
  if(num>999)   xt = char(48+mod(num/1000,10))//xt
  if(num>9999)  xt = char(48+mod(num/10000,10))//xt
  if(num>99999) xt = char(48+mod(num/100000,10))//xt
  if(num>999999) call rx( ' xt:can not produce')
  xt='.'//xt
END function xt
character(10) function i2char(numin) ! convert num to char. See charnum3 to understand this.
  implicit none
  integer(4) ::num,itens,iten,numin
  num=abs(numin)
  if(num>=10**9) call rx('i2char:num>10**9')
  i2char=''
  do itens=8,1,-1 !itens mean itens+1's digit
     iten=10**itens
     if(num>=iten) i2char=trim(i2char)//char(48+mod(num/iten,10))
  enddo
  i2char = trim(i2char)//char(48+mod(num,10)) !1st digit
  if(numin<0) i2char = '-'//trim(i2char)//char(48+mod(num,10)) !1st digit
END function i2char
character(8) function charext(numin)
  integer(4) ::num,numin
  character(1)::sss
  num=numin
  sss=''
  if(num<0) sss='-'
  if(num<0) num=abs(num)
  charext = char(48+mod(num,10))
  if(num>9)   charext= char(48+mod(num/10,10))//charext
  if(num>99)  charext= char(48+mod(num/100,10))//charext
  if(num>999) charext= char(48+mod(num/1000,10))//charext
  if(num>9999)charext= char(48+mod(num/10000,10))//charext
  if(num >99999) charext='LargeInt'
  charext=trim(sss)//charext
END function charext

character(64) function charli(iline,ndat)
  integer:: iline(*),i,ndat
  character(8):: charext
  charli=''
  do i=1, ndat
     charli =trim(charli)//' '//trim(charext(iline(i)))
  enddo
  charli=adjustl(charli)
END function charli

character(3) function charnum3(num)
  integer(4) ::num
  charnum3 = &
       char(48+mod(num/100,10))// &
       char(48+mod(num/10,10))// &
       char(48+mod(num,10))
END function charnum3
character(4) function charnum4(num)
  integer(4) ::num
  charnum4 = &
       char(48+mod(num/1000,10))// &
       char(48+mod(num/100,10))// &
       char(48+mod(num/10,10))// &
       char(48+mod(num,10))
END function charnum4
character(9) function ftoa9(arg)
  real(8):: arg
  write(ftoa9,"(1x,f8.3)") arg
END function ftoa9
character(15) function f2a(arg)
  real(8):: arg
  write(f2a,"(d13.5)") arg
END function f2a
character(20) function xts(num1,num2)
  integer(4) :: num,num1,num2
  num = num2
  xts = char(48+mod(num,10))
  if(num>9)     xts = char(48+mod(num/10,10))//xts
  if(num>99)    xts = char(48+mod(num/100,10))//xts
  if(num>999)   xts = char(48+mod(num/1000,10))//xts
  if(num>9999)  xts = char(48+mod(num/10000,10))//xts
  if(num>99999) xts = char(48+mod(num/100000,10))//xts
  if(num>999999) call rx( ' xts:can not produce')
  xts ='.L'//xts
  num = num1
  xts = char(48+mod(num,10))//xts
  if(num>9)     xts = char(48+mod(num/10,10))//xts
  if(num>99)    xts = char(48+mod(num/100,10))//xts
  if(num>999)   xts = char(48+mod(num/1000,10))//xts
  if(num>9999)  xts = char(48+mod(num/10000,10))//xts
  if(num>99999) xts = char(48+mod(num/100000,10))//xts
  if(num>999999) call rx( ' xts:can not produce')
END function xts
character(20) function xxt(num1,num2)
  integer :: num,num1,num2
  num = num2
  xxt = char(48+mod(num,10))
  if(num>9)     xxt = char(48+mod(num/10,10))//xxt
  if(num>99)    xxt = char(48+mod(num/100,10))//xxt
  if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
  if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
  if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
  if(num>999999) call rx( ' xxt:can not produce')
  xxt ='to'//xxt
  num = num1
  xxt = char(48+mod(num,10))//xxt
  if(num>9)     xxt = char(48+mod(num/10,10))//xxt
  if(num>99)    xxt = char(48+mod(num/100,10))//xxt
  if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
  if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
  if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
  if(num>999999) call rx( ' xxt:can not produce')
END function xxt
real(8) function tripl(a,b,c)  ! tripl (determinant of 3x3 matrix) ==
  double precision :: a(3),b(3),c(3)
  tripl = a(1)*b(2)*c(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2) -a(3)*b(2)*c(1) - a(2)*b(1)*c(3) - a(1)*b(3)*c(2)
END function tripl
pure subroutine cross(a,b,c)
  implicit none
  intent(in)  ::   a,b
  intent(out) ::       c
  real(8):: a(3),b(3),c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross
subroutine minv33tp(plat,qlat)
  implicit none
  real(8),intent(in)::  plat(3,3)
  real(8),intent(out):: qlat(3,3)
  real(8):: det
  call cross(plat(1,2),plat(1,3), qlat     )
  call cross(plat(1,3),plat     , qlat(1,2))
  call cross(plat     ,plat(1,2), qlat(1,3))
  det  = sum( plat(1:3,1)*qlat(1:3,1) )
  qlat = qlat/det
end subroutine minv33tp
subroutine minv33(matrix,inverse) !Inverts 3X3 matrix
  implicit none
  real(8), intent(in) :: matrix(3,3)
  real(8), intent(out) :: inverse(3,3)
  real(8) :: det,ddot,crossf(3)
  call cross(matrix(1,2),matrix(1,3),inverse     )
  call cross(matrix(1,3),matrix     ,inverse(1,2))
  call cross(matrix     ,matrix(1,2),inverse(1,3))
  det = sum(matrix(:,1)*inverse(:,1)) !ddot(3,matrix,1,inverse,1)
  inverse = transpose(inverse)
  inverse = inverse/det
end subroutine minv33
subroutine dinv33(matrix,iopt,invrse,det)  !- Inverts 3x3 matrix
  ! ----------------------------------------------------------------
  !i Inputs
  !i   matrix:  matrix to be inverted
  !i   iopt:  if 0, usual inverse
  !i             1, transpose of inverse
  !i             2, 2*pi*inverse
  !i             3, 2*pi*transpose of inverse
  !o Outputs
  !o   invrse    see iopt
  !o   det:      determinant, or det/2*pi (sign ok ??)
  !r Remarks
  !r   To generate reciprocal lattice vectors, call dinv33(plat,3,rlat)
  ! ----------------------------------------------------------------
  implicit none
  integer :: iopt,i,j
  double precision :: matrix(3,3),invrse(3,3),det,ddot
  double precision :: xx
  call cross(matrix(1,2),matrix(1,3),invrse     )
  call cross(matrix(1,3),matrix     ,invrse(1,2))
  call cross(matrix     ,matrix(1,2),invrse(1,3))
  det = ddot(3,matrix,1,invrse,1)
  if (det == 0d0) call rx('INV33: vanishing determinant')
  if (iopt >= 2) det = det/(8*datan(1d0))
  if (mod(iopt,2) == 0) then
     do    i = 1, 3
        do  j = i+1, 3
           xx = invrse(i,j)
           invrse(i,j) = invrse(j,i)
           invrse(j,i) = xx
        enddo
     enddo
  endif
  call dscal(9,1/det,invrse,1)
end subroutine dinv33
subroutine dpcopy(afrom,ato,n1,n2,fac)
  implicit none
  integer :: n1,n2
  real(8) :: afrom(*),ato(*),fac
  ato(n1:n2) = fac*afrom(n1:n2)
end subroutine dpcopy
subroutine dpzero(array,leng)
  integer :: leng
  real(8) :: array(leng)
  array=0d0
end subroutine dpzero
subroutine dpscop(afrom,ato,nel,n1,n2,fac)
  integer :: n1,n2,nel
  real(8) :: afrom(*),ato(*),fac
  ato(n2:n2+nel-1)= fac*afrom(n1:n1+nel-1) 
end subroutine dpscop
!> taken from https://community.intel.com/t5/Intel-Fortran-Compiler/Weird-Fortran/td-p/1185072?
!> for f90
!$$$Cr Remarks
!$$$C   dmach(1-3) are as returned by the BLAS subroutine dmach and are
!$$$C   defined as follows.
!$$$C        b = base of arithmetic
!$$$C        t = number of base b digits
!$$$C        l = smallest possible exponent
!$$$C        u = largest possible exponent
!$$$C   dmach(1): eps = b**(1-t)
!$$$C   dmach(2): tiny = 100.0*b**(-l+t)
!$$$C   dmach(3): huge = 0.01*b**(u-t)
!$$$C
!$$$C   d1mach(1-5) are as returned by the BLAS subroutine d1mach and are
!$$$C   defined as follows.
!$$$C   d1mach(1) = b**(l-1), the smallest positive magnitude.
!$$$C   d1mach(2) = b**(u*(1 - b**(-t))), the largest magnitude.
!$$$C   d1mach(3) = b**(-t), the smallest relative spacing.
!$$$C   d1mach(4) = b**(1-t), the largest relative spacing.
!$$$C   d1mach(5) = log10(b)
!$$$C   d1mach and dmach call the C segment mkcon found in fsubs.c
!$$$C ----------------------------------------------------------------------
function i1mach(i) result(s)
  implicit none
  integer :: i,s
  s=99999
  if(i==2 .OR. i==4) then
     s=6
  elseif(i==9) then
     s=huge(0)
  else
     call rx('i1mach not defined')
  endif
end function i1mach
function d1mach(i) result(s)
  implicit none
  integer:: i
  double precision :: s,dm(5)
  logical :: beg = .true.
  save dm
  s=1d99
  if(i < 1 .OR. i > 5)stop 'D1MACH(arg < 1 or arg > 5)'
  if(beg)then
     beg=.false.
     dm(1) = tiny(0d0)
     dm(2) = huge(0d0)
     dm(3) = epsilon(0d0)/radix(0d0)
     dm(4) = epsilon(0d0)
     dm(5) = log10(2d0)
  end if
  s = dm(i)
  return
end function d1mach
function dmach(i) result(s)
  implicit none
  integer:: i
  double precision :: s,dm(5)
  logical :: beg = .true.
  real(8):: t,b,l,u,dlamch,eps
  save dm
  if(i < 1 .OR. i > 3)stop 'DMACH(arg < 1 or arg > 3)'
  if(beg)then
     beg=.false.
     !     machine constant from lapack routine
     b = dlamch('b')        !base
     eps = dlamch('p')      !eps*base
     l = dlamch('m')        !emin
     u = dlamch('l')        !emax
     t = int(1d0-(log(eps)/log(b)))
     dm(1) = b**(1d0-t)
     dm(2) = 100d0*b**(l+t)
     dm(3) = (b**(u-t))/100d0
  endif
  s = dm(i)
  return
end function dmach
integer function isw(sw) ! Returns integer 0 logical sw is false, 1 if true
  logical :: sw
  isw = 0
  if (sw) isw = 1
end function isw
subroutine gettime(datim)
  character datim*(*)
  call ftime(datim)
end subroutine gettime
integer function nwordr() !  nword=NWORD_RECORDSIZE
   real(8):: r
   integer:: len   
   !inquire(iolength=len) r
   nwordr = 1 !8/len
!   write(6,*)' nword=',nwordr
!   call flush(6)
! !  stop
END function nwordr
logical function isanrg(i,i1,i2,t1,t2,lreqd)
  logical :: lreqd
  integer :: i,i1,i2,lgunit,iprint,k1,k2,it1
  character*(*) t1,t2
  character strn*80,strn2*80,t3*30
  isanrg=.false.
  if (i>=i1 .AND. i<=i2) return
  isanrg=.true.
  call rx(trim(t1)//' '//trim(t2))
end function isanrg
subroutine fsanrg(f,f1,f2,tol,t1,t2,lreqd)
  logical :: lreqd
  double precision :: f,f1,f2,tol
  character*(*) t1,t2
  character strn*100,strn2*100,t3*30
  if (f>=f1 .AND. f<=f2) return
  if (f1==f2 .AND. f>=f1-tol/2d0 .AND. f<=f2+tol/2d0) return
  call rx(trim(t1)//' '//trim(t2))
end subroutine fsanrg
subroutine setfac(n,fac) !- set up array of factorials.
  integer :: n,i,ik,m
  real(8):: fac(0:n)
  fac(0)=1d0
  do m=1,n
     fac(m)= product([(dble(ik),ik=1,m)])
  enddo   
!  fac=[(product([(dble(ik),ik=1,m)]),m=0,n)] !not working in nvfortran24.1
end subroutine setfac
subroutine stdfac(n,df) !- Set up array of double factorials.
  !  for odd numbers,  makes 1*3*5*..*n
  !  for even numbers, makes 2*4*6*..*n
  integer :: n,ik,m
  real(8):: df(0:n)
  df(0)=1d0
  do m=1,n
     df(m)=product([(dble(ik),ik=m,1,-2)])
  enddo
  !df=[(product([(dble(ik),ik=m,1,-2)]),m=0,n)]
end subroutine stdfac
subroutine ftime(datim)!fortran-callable date and time
  character datim*(*)
  call fdate(datim)!datim=datim(1:24) !takao. If this is not, write(6,*) gives CR at the ene of datim*26.
end subroutine ftime
subroutine readx(ifil,n)
  integer::ifil,n,i,j
  character(72) :: rchar
  do i = 1,n
     read(ifil,*)rchar
     j       = 0
     rchar=trim(adjustl(rchar))
     if(rchar(1:3) == '***')return
     if(rchar(1:3) == '###')return
  enddo
  call rx( 'readx: cannot find the string(rw.f)')
end subroutine readx
!real(8) function derfc(x)
!  real(8)::x
!  derfc= 1d0 - erf(x)
!end function derfc
subroutine getqkey(qx,nqtt,epsd,  nkey,key) !qx is digitized by epsd
!!NOTE: use this with ik=findloc( int(qinput+0.5d0*epsd) - key,value=0)
  intent(in)::    qx,nqtt,epsd
  intent(out)::                  nkey,key
  integer:: nqtt, key(nqtt),isig,ieord(nqtt),nkey,i,ik,kkk,mm
  real(8):: qx(nqtt),epsd
  call sortea(qx,ieord,nqtt,isig)
  ik=0
  do i=1,nqtt
    kkk=(qx(ieord(i))+0.5d0*epsd)/epsd  !kkk is digitized by 1/epsd
    if(i==1)then
       ik=ik+1
       key(ik) = kkk
    elseif( key(ik)<kkk) then   
       ik=ik+1
       key(ik) = kkk
    elseif (key(ik)>kkk) then
       write(6,*) ik,i, key(ik), qx(ieord(i))
       call rx( 'iqindx: bug not sorted well')
    endif
  enddo
  nkey=ik
end subroutine getqkey
real(8) function avwsr(plat,alat,vol,nbas)
  !- Calculate the average ws radius
  !     implicit none
  integer :: nbas
  double precision :: plat(3,3),alat,vol
  vol = alat**3*dabs( &
       plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) + &
       plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) + &
       plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))
  avwsr = (3/(16*datan(1.d0)*nbas)*vol)**(1.d0/3.d0)
END function avwsr
subroutine locase(ps)
  character(*) ps
  integer::i,n,shift
  n=len_trim(ps)
  shift=-ichar('A')+ichar('a')
  do i=1,n
     if ( ichar(ps(i:i)) >= ichar('A') &
          .AND. ichar(ps(i:i)) <= ichar('Z') ) then
        ps(i:i) = char( ichar(ps(i:i))+ shift )
     endif
  enddo
end subroutine locase
module m_factorial
  real(8),allocatable,protected:: factorial(:),factorial2(:)
contains
  subroutine factorial_init(kmaxf,lmaxf)
    integer:: kmaxf,lmaxf
    if(allocated(factorial)) deallocate(factorial,factorial2)
    allocate(factorial(0:kmaxf),factorial2(0:lmaxf))
    call setfac(kmaxf,factorial)
    call stdfac(lmaxf,factorial2)
  end subroutine factorial_init
end module m_factorial
! function zxx(a,b) result(ab)
!   integer:: i,j
!   complex(8) :: a(:),b(:),ab(size(a),size(b))
!   do i=1,size(a)
!      do j=1,size(b)
!         ab(i,j)=dconjg(a(i))*b(j)
!      enddo
!   enddo
! end function zxx
module NaNum
  real(8),parameter :: NaN = -9999999
!  contains
!    subroutine naninit()
!    use, intrinsic :: ieee_arithmetic
!    NaN=ieee_value(1d0,IEEE_QUIET_NAN)
!  end subroutine naninit
end module NaNum
subroutine getdval(dddin, ncount,arr) !Read undefinit number of real(8) array
  integer:: ncount, i,ix,iy
  real(8):: arr(*)
  character(*):: dddin
  character(500):: ddd !  print *,'getdval'
  logical:: debug=.false.
  ddd = trim(dddin)//' '
  if(debug) write(6,*)'getdval:'//trim(ddd)//'###'
  ncount=0
  do i=1,100
     ddd=adjustl(ddd)
     if(debug) write(6,*)'getdval##'//trim(ddd)//'###'
     if(len(trim(ddd))==0) goto 1012  !print *,'ddd:',trim(ddd)
     read(ddd,*,err=1012) arr(i)
     ncount=ncount+1
     ix = scan(ddd,' ')
     iy = scan(ddd,',')
     if(ix==0) ix=99999
     if(iy==0) iy=99999
     ddd=ddd(min(ix,iy)+1:)
  enddo
  call rx('error: getdval')
1012 continue   !write(6,*) 'getdval=',arr(1:ncount)
end subroutine getdval
character(8) function xn(num)
  integer :: num
  if(num==0) then
     xn=''
     return
  endif
  xn = char(48+mod(num,10))
  if(num>9)     xn = char(48+mod(num/10,10))//xn
  if(num>99)    xn = char(48+mod(num/100,10))//xn
  if(num>999)   xn = char(48+mod(num/1000,10))//xn
  if(num>9999)  xn = char(48+mod(num/10000,10))//xn
  if(num>99999) xn = char(48+mod(num/100000,10))//xn
  if(num>999999) call rx( ' xn:can not produce')
!  xn='.'//xn !no dot. see function xt
END function xn
character(8) function xtxx(num) !taken from xt in extension.F
  integer(4) :: num
  xtxx=''
  if(num>0)     xtxx = char(48+mod(num,10))
  if(num>9)     xtxx = char(48+mod(num/10,10))//xtxx
  if(num>99)    xtxx = char(48+mod(num/100,10))//xtxx
  if(num>999)   xtxx = char(48+mod(num/1000,10))//xtxx
  if(num>9999)  xtxx = char(48+mod(num/10000,10))//xtxx
  if(num>99999) xtxx = char(48+mod(num/100000,10))//xtxx
  if(num>999999) call rx( ' xtxx:can not produce')
END function xtxx
double precision function plegn(n,x) ! Legendre polynomical using a recursion relation
  !i Inputs
  !i   n,x
  !o Outputs
  !o   plegn: P_n(x)
  !r   Recursion relation is P_n = [(2*n-1)*x*P_(n-1) - (n-1)*P_(n-2)]/n
  ! ----------------------------------------------------------------------
  implicit none
  integer :: n,j
  double precision :: x,jpjm1,cj,pjp1
  ! jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
  jpjm1 = 0
  plegn = 1
  cj = 1
  do j = 1, n
     pjp1 = (cj*x*plegn - jpjm1)/j
     jpjm1 = j*plegn
     cj = cj + 2
     plegn = pjp1
  enddo
END function plegn
subroutine sortea(ea,ieaord,n,isig)  ! mini-sort routine.
  implicit real*8(a-h,o-z)
  real(8)::        ea(n)
  integer:: ieaord(n),n,isig,itmp,i,ix
  isig = 1
  do i = 1,n
     ieaord(i) = i
  enddo
  do ix= 2,n
     do i=ix,2,-1
        if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
           itmp = ieaord(i-1)
           ieaord(i-1) = ieaord(i)
           ieaord(i) = itmp
           isig= -isig
           cycle
        endif
        exit
     enddo
  enddo
end subroutine sortea
subroutine gvgetf(ng,n,kv,k1,k2,k3,c,c0)!- Gathers Fourier coefficients from 3D array c into list c0.
  implicit none
  integer :: ng,n,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng,n),c(k1,k2,k3,n)
  integer :: ig,i,j1,j2,j3
  do i=1,n
     c0(:,i) = [(c(kv(ig,1),kv(ig,2),kv(ig,3),i), ig=1,ng)]
  enddo   
end subroutine gvgetf
subroutine gvputf(ng,n,kv,k1,k2,k3,c0,c)!- Pokes Fourier coefficients from list c0 into 3D array c.
  implicit none
  integer :: ng,n,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng,n),c(k1,k2,k3,n)
  integer :: ig,i,j1,j2,j3
  c=0d0
  do ig=1,ng
     c(kv(ig,1),kv(ig,2),kv(ig,3),:) = c0(ig,:)
  enddo   
end subroutine gvputf
subroutine gvaddf(ng,kv,k1,k2,k3,c0,c)! Adds Fourier coefficients from list c0 into 3D array c.
  implicit none
  integer :: ng,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng),c(k1,k2,k3)
  integer :: ig,j1,j2,j3
  do  ig = 1, ng
     j1 = kv(ig,1)
     j2 = kv(ig,2)
     j3 = kv(ig,3)
     c(j1,j2,j3) = c(j1,j2,j3) + c0(ig)
  enddo
end subroutine gvaddf

pure real(8) function erfcee(ra) ! erfcee(x)= erfc(|x|)/y0/exp(-x*x)  ! We still use erfcee because gfortran erfc gives NaN for large x at 2022-7.
  ! For example, I observed erfc(25.3526)=NaN (gfortran compilar bug. it must be zero).
  ! ... erfc(x) is evaluated as a ratio of polynomials, to a relative precision of <10^-15 for x<5.
  !     Different polynomials are used for x<1.3 and x>1.3. Numerators and denominators are t,b respectively.
  implicit none
  intent(in):: ra
  real(8)::w, f1,f2, ra
  real(8),parameter:: &
       t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0, &
       t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0, &
       t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0, &
       t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0, &
       t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0, &
       t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0, &
       t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0, &
       t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0, &
       b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0, &
       b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0, &
       b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0, &
       b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0, &
       b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0, &
       b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0, &
       b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0, &
       b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0
  ! ... y0*dexp(-x*x)*f1(w=x-1/2) is erfc(x) for x<1.3
  f1(w) = (((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)*w+t11)*w+t10)&
       /  ((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)*w+b13)*w+b12)*w+b11)*w+1)
  ! ... y0*dexp(-x*x)*f2(w=x-2) is erfc(x) for x>1.3
  f2(w) = (((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)*w+t21)*w+t20)&
       /  ((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)*w+b23)*w+b22)*w+b21)*w+1)
  if (ra > 1.3d0) then
     erfcee = f2(ra - 2d0)
  elseif (ra > 0) then
     erfcee = f1(ra -.5d0)
  elseif (ra > -1.3d0) then
     erfcee = f1(-ra-.5d0)
  else
     erfcee = f2(-ra- 2d0)
  endif
end function erfcee

real(8) function fac2l(i)
  !C A table of (2l-1)!!
  !     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
  integer:: i,l
  integer,parameter:: lmx=100
  logical,save::  init=.true.
  real(8),save:: fac2lx(0:lmx+1)
  if(init) then
     fac2lx(0)=1d0
     do l=1,lmx
        fac2lx(l)=fac2lx(l-1)*(2*l-1)
     enddo
     init=.false.
  endif
  fac2l=fac2lx(i)
END function fac2l
! real(8) function fac2l(l)
!   integer:: l
!   integer,parameter:: nlmax=20,  lmx=  nlmax*2+2, lmin=-nlmax
!   integer:: fac2lx(-nlmax:nlmax*2+3)
!   logical,save:: init=.true.
!   ! --- A table of fac2l(l)=(2l-1)!!
!   !     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
!   if(init) then
!      fac2lx(0) = 1d0
!      do  10  l = 1, lmx+1
!         fac2lx(l) = fac2lx(l-1) * (l+l-1)
! 10   enddo
!      do  11  l = -1, lmin, -1
!         fac2lx(l) = fac2lx(l+1) / (l+l+1)
! 11   enddo
!      init=.false.
!   endif
!   fac2l =fac2lx(l)
! end function fac2l
