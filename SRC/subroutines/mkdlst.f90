integer function mkdlst(strn,fuzz,npmx,xp)
  !- Resolve list (ascii string) into a vector of double-precision numbers
  ! ----------------------------------------------------------------
  !i Inputs
  !i   strn  :string holding list of numbers see Remarks for syntax
  !i Inputs
  !i   fuzz  :uncertainty in the precision.  In a number sequence,
  !i         :the upper bound is known only to precision fuzz.
  !i         :Upper limit is replaced by upper limit + fuzz.
  !i         :fuzz<0 => use |fuzz|*|max element from strn|
  !i   npmx  :npmx>0  => maximum number of points allowed.
  !i          npmx<=0 => mkdlst returns np without filling xp
  !o Outputs:
  !o   xp(1..) the vector of numbers (read only if npmx>0)
  !o   mkdlst: number of points read, or would try to read if npmx>0.
  !o           If failed to parse ascii input, returns -1.
  !o           If np>npmx and npmx>0, returns -npmx
  !r Remarks
  !r   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
  !r   low:high:step
  !r   low, high, and step are real expressions specifying the sequence
  !r     low, low+step, low+2*step, ... high.
  !r   If :step is missing, the step size defaults to 1.  If also :high
  !r   is missing,  the sequence reduces to a single integer. Thus,
  !r     '5+1'       becomes a single number, 6.
  !r     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
  !r     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
  !r   Sequences may be strung together separated by commas, eg
  !r     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
  !u Updates
  !u   25 Aug 04 Added fuzz.  New argument list
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: npmx
  character*(*) strn
  double precision :: fuzz,xp(*)
  ! ... Local parameters
  double precision :: dv(1000),d1mach,dx,tfuzz
  integer :: it(1000),a2vec,ip,i,j,k,n,iprint,i1mach,np
  ip = 0
  call skipbl(strn,len(strn),ip)
  k = a2vec(strn,len(strn),ip,4,',: ',3,3,30,it,dv)
  tfuzz = fuzz
  if (tfuzz < 0) then
     call idmxmn(k,dv,1,i,j)
     tfuzz = max(dabs(dv(i)),dabs(dv(j))) * abs(fuzz)
  endif
  mkdlst = -1
  if (k < 1) return
  mkdlst = -npmx
  ! ... loop over all iv
  np = 0
  i = 0
14 i = i+1
  ! ... Case dv => a single number
  if (it(i) /= 2) then
     np = np+1
     if (npmx > 0) then
        if (np > npmx) return
        xp(np) = dv(i)
     endif
     ! ... Case dv => n1:n2[:n3]
  else
     dx = 1
     if (it(i+1) == 2 .AND. dv(i+2) /= 0) dx = dv(i+2)
     n = int((dv(i+1)+tfuzz-dv(i))/dx)
     do  12  j = 0, n
        np = np+1
        if (npmx > 0) then
           if (np > npmx) return
           xp(np) = dv(i) + j*dx
        endif
12   enddo
     i = i+1
     if (it(i) == 2) i = i+1
  endif
  if (i < k) goto 14
  ! --- Entry for np = npmax
  !  20 continue
  mkdlst = np
end function mkdlst

! --- Test of mkdlst ---
!      subroutine fmain
!C      implicit none
!      integer np,npmx,mkdlst
!      parameter (npmx=20)
!      double precision xp(npmx)

!      call pshpr(101)
!      np = mkdlst('9,3:5:.2 ',npmx,xp)
!      print *, np
!      end

subroutine idmxmn(n,a,incx,imin,imax)
  !- Finds maximum and minimum value of an double precision array
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n     :number of points to evaluate
  !i   a     :array of points
  !i   incx  :spacing between points to consider
  !o Outputs
  !o   imin  :points to lowest
  !o   imax
  !l Local variables
  !l         :
  !r Remarks
  !r   Points a(1),a(1+incx),...,a(1+incx*n) are compared
  !r
  !u Updates
  !u   25 Aug 04 First created from iyamax
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n,incx
  double precision :: a(n)
  ! ... Local parameters
  integer :: i,ix,imin,imax
  double precision :: smin,smax
  imin = 0
  imax = 0
  if (n < 1 .OR. incx <= 0) return
  imin = 1
  imax = 1
  if (n == 1) return
  ix = 1
  smin = dabs(a(1))
  smax = dabs(a(1))
  ix = ix + incx
  do  10  i = 2, n
     if (dabs(a(ix)) <= smax) go to 5
     imax = i
     smax = dabs(a(ix))
5    continue
     if (dabs(a(ix)) >= smin) go to 6
     imin = i
     smin = dabs(a(ix))
6    continue
     ix = ix + incx
10 enddo
end subroutine idmxmn

