! Manage commandline inputs. -vfoobar.
! It is better to replace this by code in python.
subroutine finits()  
  use m_ext,only:sname
  implicit none
  !! Set the numbr of iprint(), extension, and -vnam=val
  integer :: iarg
  double precision :: fcargs(1)
  logical :: lsequ,lext
  integer :: i,fext,nargf,n,it(5),iv(5),k
  logical:: cmdopt
  character strn*256
  character(100) :: extns
  !! Command line arguments and extension ---
  lext = .false.
  do iarg = 1,nargf()-1
     call getarg(iarg,strn)
     extns = strn
     !  ... v encountered ... parse variables
     if (lsequ(strn,'-v',2,' ',n)) then
        i = 2
        call parsyv(strn,len(strn),999,0,i)
     endif
     if (lsequ(strn,'-',1,' ',n)) cycle
  enddo
end subroutine finits

subroutine addsyv(nam,val,ival)
  !- Add a symbolic variable to list
  ! ----------------------------------------------------------------
  !i Inputs
  !i   nam:  name of variable
  !i   val:  value of variable (double precision)
  !o Outputs
  !o   ival  index to which variable is declared or accessed
  !r Remarks
  !r   addsyv  adds a symbolic name and value to the internal table;
  !r   lodsyv  like addsyv, except when symbolic name already exists
  !r           in the table.  In that case, lodsyv does the following:
  !r           if iopt=0, table is not altered and lodsyv returns ival=0
  !r           if iopt>0, lodsyv updates the value for the variable.
  !r   getsyv  retrieves value associated with a name
  !r   watsyv  retrieves name and value associated with index
  !r   chgsyv  changes value associated with an index
  !r   chsyv   changes value associated with a name
  !r   shosyv  displays symbolic variables and associated values
  !r   numsyv  returns the number of variables now declared
  !r   clrsyv  clears all symbolic variables nvar and beyond
  !r   togsyv  toggles the positions of variable nvar and nvar2
  !u Updates
  !u   3 Apr 00 handles `shell-command` as an argument
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character*(*) nam
  double precision :: val
  integer :: ival,first,last,ifmt,nvar,nvar2,ifi,iopt
  ! Local parameters
  integer :: mxnam,namlen
  parameter (mxnam=256,namlen=16)
  character*(namlen) symnam(mxnam), tmpnam
  double precision :: symval(mxnam),tmpval
  integer :: nnam,i,j
  save symnam, symval, nnam
  data nnam /0/

  nnam = nnam+1
  if (nnam > mxnam) stop 'addsyv: too many names'
  symnam(nnam) = nam
  call locase(symnam(nnam))
  symval(nnam) = val
  ival = nnam
  return

  entry getsyv(nam,val,ival)

  tmpnam = nam
  call locase(tmpnam)
  do  10  i = 1, nnam
     if (tmpnam /= symnam(i)) goto 10
     val = symval(i)
     ival = i
     return
10 enddo
  ival = 0
  return

  entry watsyv(nam,val,ival)
  nam = ' '
  val = 0
  if (ival <= nnam) then
     nam = symnam(ival)
     val = symval(ival)
  endif
  return

  entry chgsyv(val,ival)
  if (ival <= nnam) symval(ival) = val
  return

  entry chsyv(nam,val,ival)
  tmpnam = nam
  call locase(tmpnam)
  do  15  i = 1, nnam
     if (tmpnam /= symnam(i)) goto 15
     symval(i) = val
     ival = i
     return
15 enddo
  ival = 0
  return

  entry lodsyv(nam,iopt,val,ival)
  tmpnam = nam
  call locase(tmpnam)
  ! ... Find the name, set value it it exists
  do  16  i = 1, nnam
     if (tmpnam /= symnam(i)) goto 16
     ! ...   Update the table's value if iopt nonzero
     if (iopt == 0) then
        ival = 0
     else
        symval(i) = val
        ival = i
     endif
     return
16 enddo
  ! ... Name not in table: append it and set value
  nnam = nnam+1
  if (nnam > mxnam) stop 'addsyv: too many names'
  symnam(nnam) = tmpnam
  symval(nnam) = val
  ival = nnam
  return

  entry shosyv(first,last,ifmt,ifi)
  j = last
  if (j <= 0 .OR. j > nnam) j = nnam
  if (first > j) return
  if (first == 0) write(ifi,332)
332 format('  Var       Name                 Val')
  do  20  i = max(first,1), j
     if (ifmt == 0) write(ifi,333) i, symnam(i), symval(i)
20 enddo
333 format(i4, 4x, a20, g14.5)
  return

  entry numsyv(nvar)
  nvar = nnam
  return

  entry clrsyv(nvar)
  nnam = nvar
  return

  entry togsyv(nvar,nvar2)
  if (nvar == nvar2) return
  tmpnam = symnam(nvar2)
  symnam(nvar2) = symnam(nvar)
  symnam(nvar) = tmpnam
  tmpval = symval(nvar2)
  symval(nvar2) = symval(nvar)
  symval(nvar) = tmpval
end subroutine addsyv

!$$$#if TEST
!$$$      subroutine fmain
!$$$      implicit none
!$$$      double precision res,rev(10)
!$$$      integer ival,i1mach,a2vec,ip,nvec,nsep,iterm
!$$$      character*1 sep(2), strn*16
!$$$      data sep/':',' '/, strn /'abc=1+1:2+2:3 4 '/
!$$$
!$$$      ip = 4
!$$$      nvec = 4
!$$$      nsep = 2
!$$$      iterm = 2
!$$$C      ival = a2vec(strn,16-1,ip,4,sep,nsep,iterm,nvec,rev)
!$$$      ival = a2vec(strn,len(strn),ip,4,sep,nsep,iterm,nvec,rev)
!$$$      print *, ival,ip,': ',strn(1:ip)
!$$$      print *, (rev(ip), ip=1,ival)
!$$$      call rx('done')
!$$$
!$$$
!$$$      call lodsyv('abc',0,1d0,ival)
!$$$      print *,ival
!$$$      call shosyv(0,0,0,i1mach(2))
!$$$
!$$$      call lodsyv('xyz',0,12d0,ival)
!$$$      print *,ival
!$$$      call shosyv(0,0,0,i1mach(2))
!$$$
!$$$      call lodsyv('pqr',0,12d0,ival)
!$$$      print *,ival
!$$$      call shosyv(0,0,0,i1mach(2))
!$$$
!$$$
!$$$      call lodsyv('xyz',0,2d0,ival)
!$$$      print *,ival
!$$$      call shosyv(0,0,0,i1mach(2))
!$$$
!$$$      call lodsyv('xyz',1,3d0,ival)
!$$$      print *,ival
!$$$      call shosyv(0,0,0,i1mach(2))
!$$$      end
!$$$#endif
!$$$
