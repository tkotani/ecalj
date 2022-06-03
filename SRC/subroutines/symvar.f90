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
  return

end subroutine addsyv
subroutine parsyv(recrd,size,mxdecl,iopt,i)
  !- Parses a string for one or more variable declarations
  ! ----------------------------------------------------------------
  !i Inputs
  !i   recrd(0:*): string recrd is parsed from i to size-1
  !i   iopt:   0: var declaration of existing variable is supressed
  !i           1: var declaration of existing variable supersedes
  !i   mxdecl: parsing ends if number of declarations exceeds mxdecl
  !i   i:      parsing starts at char i (origin at 0)
  !o Outputs
  !o   i:      last character parsed
  !r Remarks
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: size,i,mxdecl,iopt
  character(1) :: recrd(0:*)
  ! Local parameters
  double precision :: dum,dum2
  integer :: j,nextop,i0,ndecl,k,ifi,m
  logical :: a2bin,lrd=.false.
  character(1) :: aops(7)
  character(255) :: a,aa,tmpdir
  data aops/'=','*','/','+','-','^',' '/

  ndecl = 0
10 continue
  call skipbl(recrd,size,i)
  if (i >= size .OR. ndecl >= mxdecl) return
  j=i
  !      call chrps2(recrd,aops,7,size,j,nextop)
  !      if (j .ge. size .or. nextop .eq. 7) return
  call chrps2(recrd,aops,6,size,j,nextop)
  if (j >= size) return
  i0 = j+1
  if (nextop > 1) then
     if (recrd(i0) /= '=') goto 999
     i0 = i0+1
  endif
  j = j-i
  if (j > 15) call fexit(-1,9,'parsyv: var def too long',0)
  call skipbl(recrd,size,i0)
  ! ... This branch writes echo `...` to file tmpdir, reads file contents
  !     and parses result
  if (recrd(i0) == '`') then
     a = 'echo `'
     k = 1
     call strcop(a(7:),recrd(i0+1),min(size-i0,len(a)),'`',k)
     if (a(k+6:k+6) /= '`') goto 999
     call rx('I think no way to go thru here')
     !        tmpdir = ' '
     !        call get_environment_variable('TMPDIR',tmpdir)
     !        if (tmpdir .eq. ' ') tmpdir = '.'
     !        call strcat(tmpdir,len(tmpdir),' ','/parsyv.xxx',11,' ',m)
     !        a(k+7:) = ' > '
     !        call strcat(a(k+10:),len(a)-10,' ',tmpdir,len(tmpdir),' ',m)
     !        call fsystm(a,m)
     !        ifi = fopnx(tmpdir,72,1,-1)
     !        rewind ifi
     !        aa = ' '
     !        read(ifi,'(a)') aa
     !        call fclose(ifi)
     !        m = 0
     !        lrd = a2bin(aa,dum,4,0,' ',m,-1)
     !        i0 = i0+k+1
  else
     lrd = a2bin(recrd,dum,4,0,' ',i0,-1)
  endif
  if (lrd) then
     ndecl = ndecl+1
     call strcop(a,recrd(i),j,'=',i)
     if (nextop == 1) then
        call getsyv(a(1:j),dum2,i)
        if (i == 0) call addsyv(a(1:j),dum,i)
        if (i /= 0 .AND. iopt /= 0) call chsyv(a(1:j),dum,i)
     else
        dum2 = 0
        call getsyv(a(1:j),dum2,i)
        if (i == 0) call addsyv(a(1:j),0d0,i)
        if (nextop == 2) dum = dum2*dum
        if (nextop == 3) dum = dum2/dum
        if (nextop == 4) dum = dum2+dum
        if (nextop == 5) dum = dum2-dum
        if (nextop == 6) dum = dum2**dum
        call chsyv(a(1:j),dum,i)
     endif
  endif
  i = i0
  goto 10
999 continue
  print *, 'parsyv: bad syntax in var def, parsing at: ', &
       (recrd(i), i=i0-1,min(i0+15,size)),' ...'
  call fexit(-1,9,'parsyv aborting',0)
end subroutine parsyv
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
