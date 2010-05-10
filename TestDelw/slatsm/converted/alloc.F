Cdelw1 var ['ipname', '-1', 'undef', '0', 'undef', 'x', 'wkinit']
Cdelw1 var ['ipploc', '-1', 'undef', '0', 'undef', 'x', 'wkinit']
Cdelw1 undel [['undef', 'wref'], 'ipname', 'wkinit']
Cdelw1 undel [['undef', 'wref'], 'ipploc', 'wkinit']
Cdelw1 undel [['call', 'wref'], 'oname', 'redfrr']
Cdelw1 undel [['call', 'wref'], 'oname', 'redfi']
C --- Simple FORTRAN dynamic memory allocation ---
C Memory is allocated from a single heap.  The heap is declared as
C an integer array in common block W by the main program.  Its size
C is fixed by the size of W as declared by the main program.  It
C is possible to allocated memory from the top of the heap, and to
C free memory already allocated from the top of the heap.  It is
C not possible to free memory allocated from the middle of the heap.
C 
C To use these routines, first call WKINIT with NSIZE = size of W.
C 
C Character, integer, real, double precision, complex or double
C complex arrays can then be allocated by calling DEFCH, DEFI,
C DEFR, DEFDR, DEFC or DEFDC with an index IPNAME and LENGTH the
C desired length.  IPNAME is returned as the effective offset of
C the integer array W.  If it is desired that the array be
C initialized to zero, pass LENGTH as the negative of the desired
C length.
C 
C Allocated memory can be freed by calling RLSE(IPNAME) where
C IPNAME must one of the indices from a previous allocation.  All
C blocks of memory allocated after IPNAME are also released.
C
C WKPRNT turns on or off a trace that prints a message each time
C        memory is allocated or released.  Also when turned on,
C        runs through links each time RLSE is called.
C DEFASK returns the amount of memory remaining in the pool
C WKSAV  sets an internal switch so that a heap allocation request
C        returns the the negative of length of the request in IPNAME
C        when not enough memory is available, rather than aborting.
C        In this case, no memory is allocated and no other internal
C        variables are changed.
C WKINFO prints information about the arrays allocated
C WKCHK  runs through links to see if any array was overwritten
C        call wkchk(string)
C ----------------------------------------------------------------
      subroutine wkinit(nsize)
C- Initialize conditions for the 'common pool' data W
C     implicit none
Cdelw1 do not change  [['undef', 'wref'], 'ipname', 'wkinit']
Cdelw1 do not change  [['undef', 'wref'], 'ipploc', 'wkinit']
      integer w(2),irnd,i,j,qqq,ndefmx,limit,nsize,ipmax,ndfdmx,ncall
     .  ,jpr,ipfree,iprint,jprint,length,leng,jopt,i1mach,iopt,imod
     .  ,ipname,ndefd,ixx,lrest,lused,lgunit,ipx,ipy,ipploc
     .  ,i0,i1,i2,i3,ippplc(-1:4),ipppln(-1:4),size,nstk
      integer mfree,mnow,mmax
      character*(*) string
      logical lopt,lerr,lfast
      save
      common /w/ w

      irnd(i) = (i+499)/1000
      size(i) = ippplc(i-1) - ippplc(i) - 1

C --- define storage size ---
C ... Start of first array and max number to be defined:
      qqq = 5
      ndefmx = 300
      limit = nsize
      ipmax = 0
      ndfdmx = 0
      ncall = 0
      jpr = 0
      ipfree = 5
      lerr = .false.
      lfast = .false.
      if (iprint() .gt. 0) then
        write(*,347) irnd(nsize)
  347   format(' WKINIT:  size=',i6,'K')
        write(*,*) ' '
      endif
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return
      entry wkprnt(jprint)
C- Set debug switch for heap management
      if (jprint .eq. 2) then
        jpr = 1-jpr
      else
        jpr = jprint
      endif
      if (jpr .eq. 1) lfast = .false.
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return
      entry wksav(lopt)
      lerr = lopt
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return
      entry wkfast(lopt)
      lfast = lopt
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- Subroutines to define arrays of various types ---
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 87
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defch(ipname,leng)
C- Allocate character array
      length = (leng+3)/4
      goto 10
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 91
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defi(ipname,leng)
C- Allocate integer array
      length = leng
      jopt = 1
      goto 10
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 96
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defr(ipname,leng)
C- Allocate single precision real array
      length = leng*i1mach(17)
      jopt = 2
      goto 10
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 101
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defrr(ipname,leng)
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 102
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defdr(ipname,leng)
C- Allocate double precision real array
      length = leng*i1mach(18)
      jopt = 3
      goto 10
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 107
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defc(ipname,leng)
C- Allocate single precision complex array
      length = leng*2*i1mach(17)
      jopt = 4
      goto 10
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 112
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defcc(ipname,leng)
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 113
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry defdc(ipname,leng)
C- Allocate double precision complex array
      length = leng*2*i1mach(18)
      jopt = 5

C --- Allocate next array from top of heap, of size length ---
   10 iopt = 0
      if (length .lt. 0) then
        iopt = 1
        length = -length
      endif
      if (length .eq. 0) length = 1
C ... Check links, come back here (84); skip if lfast on
      if (lfast) goto 84
      imod = 1
      goto 83
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 129
Cdelw1 do not change ipname because of ['undef', 'wref']
   84 ipname = ipfree
      if (lerr .and. ipfree+length+2 .gt. limit) then
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 131
Cdelw1 do not change ipname because of ['undef', 'wref']
        ipname = -LENGTH
        if (jpr .gt. 0)  print *,
     .    'ALLOC: heap storage exceeded; returning -LENGTH=',-length
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
      endif
C ... ipfree is start of next array, should it be allocated
      ipfree = ipfree + length + 1
      ipfree = 4*((ipfree+2)/4)+1
C ... For record keeping only
      ipmax = max0(ipmax,ipfree)
C ... Save as one of the links
Cdelw1 1 is digit, dropped
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 142
Cdelw1 do not change ipname because of ['undef', 'wref']
      w(ipname-1) = ipfree
      ndefd = ndefd + 1
      ndfdmx = max0(ndfdmx,ndefd)
      ncall = ncall + 1
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 146
Cdelw1 do not change ipname because of ['undef', 'wref']
      if (jpr .gt. 0) write(*,100) ndefd,leng,length,ipname,ipfree-1
  100 format(' define array',I4,':   els=',i9,'   length=',i9,',',
     .   i9,'  to',i9)
C ... Check to ensure array not too big
      if (ipfree .le. limit) then
        if (iopt .ne. 0) goto (201,202,203,204,205) jopt
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
C ---   Initialize arrays to zero ---
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 154
Cdelw1 do not change ipname because of ['undef', 'wref']
  201   call iinit(w(ipname),-leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 156
Cdelw1 do not change ipname because of ['undef', 'wref']
  202   call ainit(w(ipname),-leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 158
Cdelw1 do not change ipname because of ['undef', 'wref']
  203   call dpzero(w(ipname),-leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 160
Cdelw1 do not change ipname because of ['undef', 'wref']
  204   call cinit(w(ipname),-leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 162
Cdelw1 do not change ipname because of ['undef', 'wref']
  205   call zinit(w(ipname),-leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
      endif
      call fexit(-1,11,' ALLOC:  workspace overflow, need %i',ipfree)

C --- Subroutines to rotate top arrays on stack ---
C Initially i1 is top (0th) element; becomes first element, etc
      entry defps4(i1,i2,i3,i0)
      nstk = 3
      goto 60

      entry defps3(i1,i2,i0)
      nstk = 2
      goto 60

      entry defps2(i0,i1)
      nstk = 1
      goto 60

      entry defsw(i0,i1)
      nstk = 1
      goto 60

   60 continue
C ... Check links, save ippplc
      ipx = qqq
      ippplc(0) = -999
      do  62  i = 1, ndefmx
C   ... Exit if links messed up
        if (ipx .lt. qqq .or. ipx .gt. limit) goto 83
        if (ipx .eq. ipfree) goto 66
        do  63  j = 4, 1, -1
   63   ippplc(j) = ippplc(j-1)
        ippplc(0) = ipx - 1
Cdelw1 1 is digit, dropped
        ipx = w(ipx-1)
   62 continue
   66 continue
      ippplc(-1) = ipfree-1
      ipppln(-1) = ipfree-1
      if (ndefd.lt.nstk+1) call rx('defsw: not enough arrays allocated')
C ... Copy bottom array to temp above the top array
      if (ipfree+size(nstk) .gt. limit)
     .  call rx('defsw: no space for swap')
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
      call icopy(size(nstk),w(1+ippplc(nstk)),1,w(1+ipppln(-1)),1)
C ... Create list of the new link locations; shift w(i) to w(i+1)
      ipppln(nstk) = ippplc(nstk)
      do  65  i = nstk-1, -1, -1
        j = i
        if (j .lt. 0) j = nstk
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
        call icopy(size(j),w(1+ippplc(i)),1,w(1+ipppln(i+1)),1)
        ipppln(i) = ipppln(i+1) + size(j) + 1
*       print *, 'copy', size(j), (1+ippplc(i)), (1+ipppln(i+1))
Cdelw1 1 is digit, dropped
        w(ipppln(i+1)) = ipppln(i) + 1
   65 continue
*      print *, (ipppln(i), i=0,nstk)
*      print *, (w(ipppln(i)), i=0,nstk)
C ... Set i0,i1...i3 to elements 0..3
      i0 = ipppln(0)+1
      i1 = ipppln(1)+1
      if (nstk .ge. 2) i2 = ipppln(2)+1
      if (nstk .ge. 3) i3 = ipppln(3)+1
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- Release data up to pointer ---
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 226
Cdelw1 do not change ipname because of ['undef', 'wref']
      entry rlse(ipname)
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 227
Cdelw1 do not change ipname because of ['undef', 'wref']
      if (ipname .gt. limit) call rx('RLSE: pointer exceeds limit')
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 228
Cdelw1 do not change ipname because of ['undef', 'wref']
      if (ipname .lt. 3) call rx('RLSE: pointer less than 3')
      if (jpr .eq. 0) goto 82
C ... Check links, come back here (87)
      imod = 3      
      goto 83
   87 ixx = ipfree
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 234
Cdelw1 do not change ipname because of ['undef', 'wref']
      ipfree = ipname
C ... Check links, come back here (91)
      imod = 4
      goto 83
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 238
Cdelw1 do not change ipname because of ['undef', 'wref']
   91 print 333 ,ndefd+1,ixx-ipname,ipname,ixx
  333 format(' RLSE array',i6,',  size=',i9,':',i9,' to',i9)
Cdelw1 warning(1) , probably  ipname  is not defined yet at linenumber= 240
Cdelw1 do not change ipname because of ['undef', 'wref']
   82 ipfree = ipname
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- Return number of words left in common pool ---
      entry defask(lrest)
      lrest = limit - ipfree - 2
      if (jpr .gt. 0) print 382, lrest
  382 format(' DEFASK:  space left=',i9,'  words')
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return
      entry wkused(lused)
      lused = ipfree
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

      entry wquery(mfree,mnow,mmax)
      mfree = limit-ipfree-2
      mnow = ipfree-1
      mmax = ipmax-1
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- Output workspace information ---
      entry wkinfo
C ... Run through links, come back here (81)
      imod = 2
      goto 83
  81  write(lgunit(1),600) irnd(ipmax-1),irnd(limit),irnd(ncall)
      if (lgunit(2) .gt. 0)
     .write(lgunit(2),600) irnd(ipmax-1),irnd(limit),irnd(ncall)
  600 format(' wkinfo:  used ',i8,' K  workspace of ',i8,
     .  ' K   in',i6,' K calls')
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (ipfree .eq. qqq) return
      if (jpr .gt. 0) write(*,602)
  602 FORMAT(/'    array',5X,'begin',7X,'end',6X,'length')
      ipx = qqq
      do  30  i = 1, ndefmx
Cdelw1 1 is digit, dropped
        ipy = w(ipx-1)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        if (ipx .eq. ipfree) return
        if (ipy .lt. qqq .or. ipy .gt. limit) write(*,*) '   . . . . . '
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        if (ipy .lt. qqq .or. ipy .gt. limit) return
        if (jpr .gt. 0) write(*,603) i,ipx,ipy,ipy-ipx
  603   format(4(i9,3x))
        ipx = ipy
   30 continue
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- Check links ---
      entry wkchk(string)

C --- Run through links to see if any dynamic array was overwritten ---
      imod = 0
      print *, 'WKCHK:  ', string
   83 ndefd = 0
      ipx = qqq
Cdelw1 warning(1) , probably  ipploc  is not defined yet at linenumber= 292
Cdelw1 do not change ipploc because of ['undef', 'wref']
      ipploc = -999
      do  35  i = 1, ndefmx
        if (ipx .lt. qqq .or. ipx .gt. limit) then
Cdelw1 warning(1) , probably  ipploc  is not defined yet at linenumber= 295
Cdelw1 do not change ipploc because of ['undef', 'wref']
          write(*,888) ndefd,ipx,ipploc
  888     format(' ALLOC: link destroyed at start of array',i3,
     .     ',  PTR=',i9,' AT',i9)
          call fexit(-1,19,' ',0)
        endif
        if (ipx .eq. ipfree) goto 86
        ndefd = ndefd + 1
Cdelw1 warning(1) , probably  ipploc  is not defined yet at linenumber= 302
Cdelw1 do not change ipploc because of ['undef', 'wref']
        ipploc = ipx - 1
Cdelw1 warning(1) , probably  ipploc  is not defined yet at linenumber= 303
Cdelw1 do not change ipploc because of ['undef', 'wref']
  35  ipx = w(ipploc)
  86  continue
      goto (84,81,87,91,66), imod
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ipname ipploc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine cinit(array,leng)
C- Initializes complex array to zero
      integer leng
      real array(2*leng)
      call ainit(array,leng+leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine zinit(array,leng)
C- Initializes complex array to zero
      integer leng
      double precision array(2*leng)
      call dpzero(array,leng+leng)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine iinit(array,leng)
C- Initializes integer array to zero
      integer leng,i
      integer array(leng)
      do  10  i=1, leng
   10 array(i) = 0
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine ainit(array,leng)
C- Initializes real array to zero
      integer leng,i
      real array(leng)
      do  10  i=1, leng
   10 array(i) = 0
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine dpzero(array,leng)
C- Initializes double precision array to zero
      integer leng,i
      double precision array(leng)
      do  10  i=1, leng
   10 array(i) = 0
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      real function rval(array,index)
      integer index
C- Returns the real value of ARRAY(INDEX)
      real array(index)
      rval = array(index)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      double precision function dval(array,index)
      integer index
C- Returns the double precision value of ARRAY(INDEX)
      double precision array(index)
      dval = array(index)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      integer function ival(array,index)
C- Returns the integer value of ARRAY(INDEX)
      integer index
      integer array(index)
      ival = array(index)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      integer function ival2(array,nda,i1,i2)
C- Returns the integer value of ARRAY(i1,i2)
C     implicit none
      integer nda,i1,i2,array(nda,1)
      ival2 = array(i1,i2)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      logical function logval(array,index)
C- Returns the integer value of ARRAY(INDEX)
      integer index
      logical array(index)
      logval = array(index)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      complex function cval(array,index)
C- Returns the complex value of ARRAY(INDEX)
      integer index
      complex array(index)
      cval = array(index)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine dvset(array,i1,i2,val)
C- Sets some elements of double precision array to value
      integer i1,i2
      double precision array(i2),val
      integer i
      do  10  i = i1, i2
   10 array(i) = val
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine ivset(array,i1,i2,val)
C- Sets some elements of integer array to value
      integer i1,i2,array(1),val,i
      do  10  i = i1, i2
   10 array(i) = val
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine lvset(array,i1,i2,val)
C- Sets some elements of logical array to value
      integer i1,i2,i
      logical array(1),val
      do  10  i = i1, i2
   10 array(i) = val
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine redfrr(oname,leng)
C- Release to pointer oname, reallocate double oname of length leng
C     implicit none
Cdelw1 do not change  [['call', 'wref'], 'oname', 'redfrr']
      integer oname,leng
Cdelw1 rlse name= oname old_list= (None)

#error ERROR, try to release name= oname ,but list=null at linenumber= 401 list= (None)
Cdelw1       call rlse(oname)

Cdelw1 do not change oname because of ['call', 'wref']
      call defrr(oname,leng)
Cdelw1 w_varlist remains: oname 
Cdelw1 w_varlistundel: [oname]
Cdelw1 w_varlist (undel), remains: [oname]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oname because of [call wref]
      end
      subroutine redfi(oname,leng)
C- Release to pointer oname, reallocate double oname of length leng
C     implicit none
Cdelw1 do not change  [['call', 'wref'], 'oname', 'redfi']
      integer oname,leng
Cdelw1 rlse name= oname old_list= (None)

#error ERROR, try to release name= oname ,but list=null at linenumber= 408 list= (None)
Cdelw1       call rlse(oname)

Cdelw1 do not change oname because of ['call', 'wref']
      call defi(oname,leng)
Cdelw1 w_varlist remains: oname 
Cdelw1 w_varlistundel: [oname]
Cdelw1 w_varlist (undel), remains: [oname]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oname because of [call wref]
      end

