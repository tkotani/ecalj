subroutine addsvv(nam,nelt,ival)
  !- Add a symbolic vector to list
  ! ----------------------------------------------------------------
  !i Inputs
  !i   nam:  name of variable
  !i   nelt: number of elements of the vector
  !o Outputs
  !o   ival  index to which variable is declared or accessed
  !r Remarks
  !r   addsvv  adds a symbolic name and value to the internal table,
  !r           and allocates memory for the vector.
  !r   lodsvv  sets a range of elements of a vector associated with
  !r           a name or an index, depending on iopt.
  !r           iopt=0: index associated with name
  !r           iopt=1: name associated with index
  !r   getsvv  gets a range of elements of a vector associated with
  !r           a name or an index, depending on iopt.
  !r   sizsvv  returns the length of a vector associated with
  !r           a name or an index, depending on iopt.
  !r   numsvv  returns the number of variables now declared
  !r   watsvv  returns name associated with index
  !r
  !r   Compiler must have either F90 or POINTER capability
  !u Updates
  !u   18 Jan 06 works with F90 compiler
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character*(*) nam
  double precision :: vec(1)
  integer :: ival,first,last,nelt,nvar,ifi,iopt
  ! Local parameters
  integer :: mxnam,namlen
  parameter (mxnam=24,namlen=16)
  character*(namlen) symnam(mxnam), tmpnam
  integer :: size(mxnam)
  integer :: nnam,i,io,i1,i2,i2x
  double precision :: x1,xn
  type row
     real(8), allocatable :: p(:)
  end type row
  type(row) :: symptr(mxnam)
  save symptr
  save symnam, size, nnam
  data nnam /0/
  ! --- Start of addsvv ---
  nnam = nnam+1
  if (nnam > mxnam) call rx('addsvv: too many names')
  symnam(nnam) = nam
  ival = nnam
  call locase(symnam(nnam))
  allocate (symptr(nnam)%p(1:nelt))
  call dvset(symptr(nnam)%p,1,nelt,0d0)
  size(nnam) = nelt
  return

  ! --- lodsvv, getsvv ---
  entry lodsvv(nam,ival,iopt,i1,i2,vec)

  io = -1
  goto  10

  entry getsvv(nam,ival,iopt,i1,i2,vec)

  io = 1
  goto  10

  entry sizsvv(nam,ival,iopt,i1)

  io = -2
  goto  10

  ! --- lodsvv, getsvv ---
  entry numsvv(nvar)
  nvar = nnam
  return

  ! --- watsvv ---
  entry watsvv(nam,ival)
  nam = ' '
  if (ival <= nnam) nam = symnam(ival)
  return

  ! --- Print out table ---
  entry shosvv(first,last,ifi)
  write(ifi,332)
332 format('  Vec       Name            Size   Val[1..n]')
  do  60  i = max(first,1), min(last,nnam)
     call dpscop(symptr(i)%p,x1,1,1,1,1d0)
     call dpscop(symptr(i)%p,xn,1,size(i),1,1d0)
     write(ifi,333) i, symnam(i), size(i), x1, xn
60 enddo
333 format(i4, 4x, a20, i4, 2g14.5)
  return

  ! --- Find an index associated with a name ---
10 continue
  ! ... If iopt=0, find the index associated with this name
  if (iopt == 0) then
     ival = 0
     tmpnam = nam
     call locase(tmpnam)
     do  16  i = 1, nnam
        if (tmpnam /= symnam(i)) goto 16
        ival = i
        goto 20
16   enddo
  endif
  ! --- Set/Retrieve portions of an array[index], depending on io ---
20 continue
  if (io == 0) return
  if (io == -2) then
     i1 = size(ival)
     return
  endif
  ! ... Return unless ival in range
  if (ival <= 0 .OR. ival > nnam) return
  i2x = min(i2,size(ival))
  if (i2x < i1) return
  if (io == -1) call dpscop(vec,symptr(ival)%p,i2x-i1+1,1,i1,1d0)
  if (io ==  1) call dpscop(symptr(ival)%p,vec,i2x-i1+1,i1,1,1d0)
  return
end subroutine addsvv

subroutine parsvv(recrd,recl,indx,mxelt,i1,ip)
  !- Parses a string for one or more elements of a vector variable
  !     implicit none
  ! Passed parameters
  integer :: recl,ip,mxelt,indx,i1
  character recrd*(100)
  ! Local parameters
  double precision :: res
  integer :: nelt,i,k,ix,a2vec
  nelt = 0
  do  33  i = 1, 999
     call skipbl(recrd,recl,ip)
     if (ip >= recl .OR. nelt >= mxelt) goto 99
     k = a2vec(recrd,recl,ip,4,' ',1,1,1,ix,res)
     if (k == -1) return
     call lodsvv(' ',indx,1,i1+nelt,i1+nelt,res)
     nelt = nelt+k
33 enddo
99 continue
end subroutine parsvv
