#define F90 1
#define unix 1
      subroutine addsvv(nam,nelt,ival)
C- Add a symbolic vector to list
C ----------------------------------------------------------------
Ci Inputs
Ci   nam:  name of variable
Ci   nelt: number of elements of the vector
Co Outputs
Co   ival  index to which variable is declared or accessed
Cr Remarks
Cr   addsvv  adds a symbolic name and value to the internal table,
Cr           and allocates memory for the vector.
Cr   lodsvv  sets a range of elements of a vector associated with
Cr           a name or an index, depending on iopt.
Cr           iopt=0: index associated with name
Cr           iopt=1: name associated with index
Cr   getsvv  gets a range of elements of a vector associated with
Cr           a name or an index, depending on iopt.
Cr   sizsvv  returns the length of a vector associated with
Cr           a name or an index, depending on iopt.
Cr   numsvv  returns the number of variables now declared
Cr   watsvv  returns name associated with index
Cr 
Cr   Compiler must have either F90 or POINTER capability
Cu Updates
Cu   18 Jan 06 works with F90 compiler
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) nam
      double precision vec(1)
      integer ival,first,last,nelt,nvar,ifi,iopt
C Local parameters
      integer mxnam,namlen
      parameter (mxnam=24,namlen=16)
      character*(namlen) symnam(mxnam), tmpnam
      integer size(mxnam)
      integer nnam,i,io,i1,i2,i2x
      double precision x1,xn
#if F90
      type row
        real(8), pointer :: p(:)
      end type row
      type(row) :: symptr(mxnam)
      save symptr
#elif POINTER
      double precision arr(20)
      pointer (iarr, arr)
      integer symptr(mxnam)
      save symptr
#endif
      save symnam, size, nnam
      data nnam /0/

C --- Start of addsvv ---
      nnam = nnam+1
      if (nnam .gt. mxnam) call rx('addsvv: too many names')
      symnam(nnam) = nam
      ival = nnam
#if unix
      call locase(symnam(nnam))
#endif
#if F90
      allocate (symptr(nnam)%p(1:nelt))
      call dvset(symptr(nnam)%p,1,nelt,0d0)
      size(nnam) = nelt
#elif POINTER
      call faloc(iarr,4,nelt)
      symptr(nnam) = iarr
      size(nnam) = nelt
#endif
      return

C --- lodsvv, getsvv ---
      entry lodsvv(nam,ival,iopt,i1,i2,vec)

      io = -1
      goto  10

      entry getsvv(nam,ival,iopt,i1,i2,vec)

      io = 1
      goto  10

      entry sizsvv(nam,ival,iopt,i1)

      io = -2
      goto  10

C --- lodsvv, getsvv ---
      entry numsvv(nvar)
      nvar = nnam
      return

C --- watsvv ---
      entry watsvv(nam,ival)
      nam = ' '
      if (ival .le. nnam) nam = symnam(ival)
      return

C --- Print out table ---
      entry shosvv(first,last,ifi)
      write(ifi,332)
  332 format('  Vec       Name            Size   Val[1..n]')
#if F90
      do  60  i = max(first,1), min(last,nnam)
C       should extract element from symptr directly; don't know how
        call dpscop(symptr(i)%p,x1,1,1,1,1d0)
        call dpscop(symptr(i)%p,xn,1,size(i),1,1d0)
        write(ifi,333) i, symnam(i), size(i), x1, xn
C        write(ifi,334) i, symnam(i), size(i), symptr(i)%p
C  334   format(i4, 4x, a20, i4, 20g10.3)
   60 continue
#elif POINTER
      do  60  i = max(first,1), min(last,nnam)
      iarr = symptr(i)
   60 write(ifi,333) i, symnam(i), size(i), arr(1), arr(size(i))
#endif
  333 format(i4, 4x, a20, i4, 2g14.5)
      return

C --- Find an index associated with a name ---
   10 continue
C ... If iopt=0, find the index associated with this name
      if (iopt .eq. 0) then
        ival = 0
#if ! (F90 | POINTER)
        return
#endif
        tmpnam = nam
#if unix
        call locase(tmpnam)
#endif
      do  16  i = 1, nnam
        if (tmpnam .ne. symnam(i)) goto 16
        ival = i
        goto 20
   16 continue
      endif

C --- Set/Retrieve portions of an array[index], depending on io ---
   20 continue
      if (io .eq. 0) return
      if (io .eq. -2) then
        i1 = size(ival)
        return
      endif

C ... Return unless ival in range
      if (ival .le. 0 .or. ival .gt. nnam) return
#if F90
      i2x = min(i2,size(ival))
      if (i2x .lt. i1) return
      if (io .eq. -1) call dpscop(vec,symptr(ival)%p,i2x-i1+1,1,i1,1d0)
      if (io .eq.  1) call dpscop(symptr(ival)%p,vec,i2x-i1+1,i1,1,1d0)
#elif POINTER
      iarr = symptr(ival)
      i2x = min(i2,size(ival))
      if (i2x .lt. i1) return
      if (io .eq. -1) call dpscop(vec,arr,i2x-i1+1,1,i1,1d0)
      if (io .eq.  1) call dpscop(arr,vec,i2x-i1+1,i1,1,1d0)
#endif

C      do  i = 1, nnam
C        call xxx(i,symptr(i)%p,size(i))
C      enddo
      return

      end
      subroutine parsvv(recrd,recl,indx,mxelt,i1,ip)
C- Parses a string for one or more elements of a vector variable
C     implicit none
C Passed parameters
      integer recl,ip,mxelt,indx,i1
      character recrd*(100)
C Local parameters
      double precision res
      integer nelt,i,k,ix,a2vec

      nelt = 0
      do  33  i = 1, 999
        call skipbl(recrd,recl,ip)
        if (ip .ge. recl .or. nelt .ge. mxelt) goto 99
        k = a2vec(recrd,recl,ip,4,' ',1,1,1,ix,res)
        if (k .eq. -1) return
        call lodsvv(' ',indx,1,i1+nelt,i1+nelt,res)
        nelt = nelt+k
   33 continue

   99 continue
      end

C      subroutine xxx(i,vec,n)
C      integer n,i
C      double precision vec(n)
C      
C      print 333, i,vec
C  333 format(i2,':',20f8.3)
C      end

