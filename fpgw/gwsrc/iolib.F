C#define FORMATTED_READ
C --- FORTRAN INPUT LIBRARY ---
C
C The i/o library provides fortran support for unformatted input
C resembling a NAMELIST statement, though it is more general.  As a
C guiding philosophy, the execution a set of fortran programs is
C intended to be governed by a single, small control file, that serves
C the dual purpose of governing program execution and document what is
C intended to calculate.  It is intended that this file is to be only
C read.  Only portions of a control file might be used by a single
C program; several programs can share a common control file.
C
C As with the fortran NAMELIST statement a quantity is identified by a
C token, such as 'NSPIN='.  The "contents" of the token are a vector (of
C length 0,1 or more) elements of cast char, logical, integer, real or
C double: 'NSPIN= 1' is a typical input for the control file
C
C Input tokens are actually separated into "categories" and data is read
C category by category.  A token always belongs to a category and tokens
C within different categories are inequivalent.  A new category begins
C (and any previous category ends) as soon as a non-blank character
C appears in the first column of a control file.  In the control file
C
C HEADER  Example of the beginning of a typical
C         control file
C VERS    1
C IO      SHOW=F HELP=F VERBOS=3 WKP=F
C
C the category 'HEADER' begins with the 'H' and continues until the 'V'
C in 'VERS', and so on.  Apart from this, input is essentially
C free-format.  Categories need not be in any particular order, but only
C the first occurence of two categories with equivalent names will be
C read.  This is handy for keeping alternative input within one file.
C Categories can either be optional or mandatory; similarly tokens
C within categories can be either optional or mandatory.  It is up to
C the calling routine to handle optional input that is not found.
C
C Names of both tokens and categories end by specifying the last
C character of the name, such as '=' in 'NSPIN='.
C
C
C --- Use of the routines ---
C
C In practice the routines are called as followed: (1) routine rdfile is
C called to load the entire control file into RAM; (2) a category is
C sought and delimited by calling routine getcat; (3) a token is sought
C and contents loaded by routine partok.  Step 3 is repeated for as many
C elements within a category as are desired, and a new category may be
C sought after the first is finished.
C
C Routines getcat and partok can be called in one of three "modes" (the
C mode in effect is governed by variable optio).  In the first mode
C getcat and partok do not attempt to search for anything; instead they
C display what they would have searched for had the next mode been
C selected.  In the next mode, they carry out the requested searches; in
C the next they print out the contents of the variables corresponding to
C the tokens (essentially an all-purpose output to verify what was
C input).  Again, it is best to look at a practical implementation to
C see how this, and many other small options can be used.
C
C
C --- Subcategories ---
C
C There is a facility to further restrict the range of search for tokens
C within a category.  Routine partok search between offsets iostart and
C subsiz from the beginning of a category.  Unless otherwise specified,
C subsiz defaults to the size of the whole category, but it is possible
C to shrink subsiz by a call to routine subcat.  This is particularly
C useful when looking for several occurences of the same token (and
C related information) within a category
C
C
C --- Other routines ---
C
c iprint is moved to iprint.f
C
C Routine fopna is the file opening routine; there are associated
C entries fopn, fopnn, fopno, fxst, fhndl which when used in conjunction
C with entry fadd, that file logical units with names, so that files
C can be subsequently referenced by name only.
C Routine fclose closes files; dfclose closes and deletes a file.
C All files should be closed by a call to fclose.
C
C Subroutine Query provides a simple "interactive mode".  When query is
C called, nothing happens unless the "query" switch has been set by a
C call to entry initqu.  But if this switch has been set, the user has
C an option to abort program execution, change the verbosity, turn on
C the work array debug switch, or to change a value of a single number
C passed to query.  There is also the option to reset the query switch
C
C There are several string-manipulation routines (copy, concatenate,
C equality, display).
C
C Routine finits performs machine-dependent initialization.
C Routine fexit performs machine-dependent program termination.
C
C --- Variables with common meaning ---
C In common block iolib:
C   recoff: offset to recrd corresponding to first character in file
C           (set by calling program before rdfile is called)
C   reclen: length of file records
C           (can be set in function recln before rdfile is called)
C   nrecs:  number of file records (set by rdfile)
C   maxlen: the maximum allowed number of characters in a category or
C           token (set by calling program before rdfile is called)
C   catbeg: offset to recrd demarcating beginning of a category
C           (set whenever getcat is called)
C   catsiz: the size of a category
C           (set whenever getcat is called)
C   subsiz: the size of a collection of tokens within a subcategory
C           (defaults to catsiz when getcat is called;
C            shrunk when subcat called)
C   iend:   points to the first character after the last match
C           (set whenever getcat or partok is called)
C   ichoos,nchoos: index to current choice of possible inputs and
C           total number of choices, respectively (see description)
C   noerr:  error on return from getcat or partok; see getcat and partok
C           (set whenever getcat or partok is called)
C   optio:  0, 1, or 2 (see description)
C ==============
      subroutine finits
C- Machine and compiler-dependent initialization
C#ifdefC SVS | unix
C      logical lsequ
C      integer j,fext,iargc,iarg,n
C      character*4 extns
C#endif

c      external diolib

c      call initqu(.false.)

C --- For Lahey F77L, open standard output as list-directed ---
C#ifdefC F77LAHEY
C      open(unit=*,carriage control='LIST')
C#endif

C --- Handle floating point exceptions in the IBM VM environment ---
C#ifdefC IBM_VM
C      call errset(208,999,-1)
C#endif

C --- Read in extension from command line if SVS or unix compiler ---
C#ifdefC SVS | unix
C      iarg = 1
C   10 iarg = iarg+1
C        if (iargc() .ge. iarg) then
C          call getarg(iarg,extns)
C          if (lsequ(extns,'-',1,' ',n)) goto 10
C          j = fext(extns)
C        endif
C#endif
      end
      subroutine fexit(retval,strng)
C- Machine and compiler-dependent program termination
      integer retval
      integer fopn,fhndl,lgunit
      double precision cpusec
C#ifdefC APOLLO_BUG
C      character*(10) strng
C#else
      character*(*) strng
C#endif
      print *
      print *, 'Stop:  ', strng
      if (iprint() .ge. 10 .and. cpusec() .ne. 0) then
        write (lgunit(1),10) cpusec()
        write (lgunit(2),10) cpusec()
   10   format(' CPU time, sec:', f12.3)
      endif
      call wkinfo()
      if (fhndl('TMP') .ge. 0) call dfclos(fopn('TMP'))
C#ifdefC unix
C      call set_ret_val(retval,1)
C#endif
Cstop2rx 2013.08.09 kino      stop
      call rx( '')
      end
      double precision function cpusec()
C- process cputime, in seconds
C ----------------------------------------------------------------------
Ci Inputs:
Ci   none
Co Outputs:
Co   returns cpu time, in seconds
Cr  Remarks
Cr    On the Apollo: time (in 4 microsecond units) is
Cr    (time(1)*65536 + time(2))*65536 + time(3)
C ----------------------------------------------------------------------
C#ifdefC APOLLO
C      integer time(3)
C      integer*2 t2(3)
C      call proc1_$get_cput(t2)
C      time(1) = t2(1)
C      time(2) = t2(2)
C      time(3) = t2(3)
C      if (time(1) .lt. 0) time(1) = time(1) + 65536
C      if (time(2) .lt. 0) time(2) = time(2) + 65536
C      if (time(3) .lt. 0) time(3) = time(3) + 65536
C      cpusec = 4d-6 *
C     .  ((dble(time(1))*65536 + dble(time(2)))*65536 + dble(time(3)))
C#else
      cpusec = 0
C#endif
      end
      integer function recln(i)
C- Returns (and optionally sets) the record length of ascii input files
      integer i
C For io routines ...
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio
      if (i .gt. 0) reclen = i
      recln = reclen
      end
C --- INTERACTIVE MODE ROUTINES ---
      subroutine query(instr,cast,var)
      call rx('not support now')
c$$$      
c$$$C- interactive flow control
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   strng: prompt string
c$$$Ci   cast:  <0, if nothing to change, otherwise
c$$$Ci          cast is  0=,logical, 2=int, 3=real 4=double
c$$$Co Outputs
c$$$Co   var:   query will change if requested
c$$$Cr Remarks
c$$$Cr   At the prompt, user enters either nothing, or one of
c$$$Cr     'Snnn', where nnn is number (or T or F for logical variable);
c$$$Cr     'Vnnn', where nnn is the new verbosity;
c$$$Cr     'W' to toggle printing of work array;
c$$$Cr     'I' to turn off interactive prompt;
c$$$Cr     'A' to abort program execution
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed parameters
c$$$      character*(*) instr
c$$$      integer cast,var
c$$$      logical lset
c$$$C Local parameters
c$$$      integer lprmpt
c$$$      parameter (lprmpt=62)
c$$$      character*(lprmpt) prmpt, rdstr
c$$$      logical lquery,lsequ,cnvt
c$$$      integer i,imin,ivbs,iprint,j
c$$$      save lquery
c$$$      character*(1) prmpt2(1),rdstr2(2)
c$$$      equivalence (prmpt,prmpt2), (rdstr,rdstr2)
c$$$      external lsequ,cnvt,iprint
c$$$
c$$$      data prmpt
c$$$     ./'(S)et parm   (V)erbos   (W)kp toggle   (I)act toggle   (A)bort'/
c$$$
c$$$C#ifndef CRAY | BATCH
c$$$      if (.not. lquery) return
c$$$
c$$$      print *
c$$$      imin = 1
c$$$c apr9 2008 takao
c$$$c      goto (1,2,1,1,1), cast+1
c$$$      if( cast+1==2 ) then
c$$$        goto 2
c$$$      elseif( cast+1 <= 5) then
c$$$        goto 1
c$$$      endif
c$$$
c$$$C Case nothing to change
c$$$C#ifdefC APOLLO_BUG
c$$$C      print *, 'QUERY:  '
c$$$C#else
c$$$      print *, 'QUERY:  ',instr
c$$$C#endif
c$$$      imin=14
c$$$      goto 10
c$$$Cstop2rx 2013.08.09 kino    2 stop 'QUERY:  cannot change char or real'
c$$$    2 call rx( 'QUERY:  cannot change char or real')
c$$$    1 continue
c$$$      call disply(instr,'=',var,var,var,var,' ',1,cast)
c$$$   10 continue
c$$$      print 333, (prmpt2(i), i=imin, lprmpt)
c$$$  333 format(' ',62a1)
c$$$  334 format(a62)
c$$$      print 335, iprint()
c$$$C#ifdefC SVS
c$$$C  335 format(' V=',i3,' ? '\)
c$$$C#elseifC APOLLO
c$$$C  335 format(' V=',i3,' ? '$)
c$$$C#else
c$$$  335 format(' V=',i3,' ?')
c$$$C#endif
c$$$      rdstr = ' '
c$$$      read(*,334) rdstr
c$$$
c$$$      i = 0
c$$$      if (lsequ(rdstr,'S',1,' ',j) .or. lsequ(rdstr,'s',1,' ',j)) then
c$$$        if (.not. cnvt(rdstr2(2),var,var,var,var,cast,0,' ',i))
c$$$     .    print *, 'conversion error'
c$$$        goto 10
c$$$      else
c$$$     .if (lsequ(rdstr,'V',1,' ',j) .or. lsequ(rdstr,'v',1,' ',j)) then
c$$$        if (.not. cnvt(rdstr2(2),0,ivbs,0,0,2,0,' ',i))  then
c$$$          print *, 'conversion error'
c$$$        else
c$$$          call popprt
c$$$          call pshprt(ivbs)
c$$$        endif
c$$$        goto 10
c$$$      else
c$$$     .if (lsequ(rdstr,'W',1,' ',j) .or. lsequ(rdstr,'w',1,' ',j)) then
c$$$        call wkprnt(2)
c$$$        goto 10
c$$$      else
c$$$     .if (lsequ(rdstr,'I',1,' ',j) .or. lsequ(rdstr,'i',1,' ',j)) then
c$$$        lquery = .false.
c$$$      else
c$$$     .if (lsequ(rdstr,'A',1,' ',j) .or. lsequ(rdstr,'a',1,' ',j)) then
c$$$        call fexit(1,'QUERY')
c$$$      endif
c$$$
c$$$C#endif
c$$$      return
c$$$
c$$$      entry initqu(lset)
c$$$      lquery = lset
c$$$
      end
C --- FREE-FORMAT INPUT PROCEDURES ---
      subroutine rdfile(unit,recrd,mxrecs,a,reclen,nrecs)
C- Read entire file into RAM, up to a maximum allowed number of records
C ----------------------------------------------------------------
Ci Inputs
Ci   unit,recrd,mxrecs,a,reclen,nrecs
Co Outputs
Co   nrecs
Cr Remarks
Cr   Reads in an entire file from unit.
Cr   Error if nrecs exceeds mxrecs
C ----------------------------------------------------------------
C     implicit none
      integer unit,mxrecs,reclen,nrecs
      logical rdstrn
      integer j
      character*1 recrd(0:*),a(*)

      nrecs = 0
   10 continue
      if (.not. rdstrn(unit,a,reclen,.false.)) return
      call strcop(recrd(reclen*nrecs),a,reclen,'',j)
      nrecs = nrecs + 1
      if (nrecs .ge. mxrecs) return

c      print *, '------------------------'
c      print 20, (recrd(j), j= 0, reclen*nrecs-1)
c   20 format(1x,72a1)

      goto 10
      end
      subroutine getcat(recrd,cat,term,lopt)
C- find a category in a record
C ----------------------------------------------------------------
Ci Inputs
Ci   recrd,cat,term
Ci   lopt: if error, aborts if unmatched category
Ci   optio: 0: show category that would have been sought
Ci          1: attempt to find category
Ci          2: print out category sought
Co Outputs
Co   catbeg,catsiz,noerr (see remarks)
Cr Remarks
Cr   on return, noerr is always true unless a category was actually
Cr   sought (optio=1) and none was found.
C ----------------------------------------------------------------
C     implicit none
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz
      character*1 recrd(0:*),term
      character*(*) cat
      integer irecs,i,j,iend,ichoos,nchoos,optio
      logical lopt,noerr
      logical lsequ
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      noerr = .true.
55    format(70a1)

capr9 2008
c      goto (1,2,3), 1+optio
      if(1+optio==1) goto 1
      if(1+optio==2) goto 2
      if(1+optio==3) goto 3

C --- print message showing category sought ---
    1 continue
      j = 5
      if (lopt) j = 1
      call msgio(recrd(nrecs*reclen),cat,' ',' ',1,j,.true.,i)
      return

C --- print out category ---
    3 continue
C#ifdefC APOLLO_BUG
C      j = 0
C      call chrpos(cat,term,maxlen-1,j)
C      print 11, cat(1:j)
C#else
      print 11, cat
C#endif
   11 format(1x,60a)
      return

C --- find a category ---
    2 continue
      catbeg = recoff
      catsiz = 0
      subsiz = 0
      irecs = nrecs
   10 continue
      irecs = irecs-1
      if (.not. lsequ(recrd(catbeg),cat,maxlen,term,iend)) then
        catbeg = catbeg + reclen
        if (irecs .eq. 0)  then
          noerr = .false.
          if (lopt) then
            call msgio(recrd(nrecs*reclen),cat,' ',' ',3,1,.true.,i)
Cstop2rx 2013.08.09 kino            stop
            call rx( '')
          endif
          return
        endif
        goto 10
      endif
C Found a category.  Next determine the length of the input string
      catsiz = catbeg
   20 continue
      irecs = irecs-1
      catsiz = catsiz + reclen
      if (recrd(catsiz) .eq. term .and. irecs .ge. 0) goto 20
C Clean up and quit
      catsiz = catsiz - catbeg
      subsiz = catsiz
      return
      end
      integer function partok(categ,token,term,result,chr,
     .                  count,cast,istart,lopt)
C- reads a vector of numbers found after a token
C ----------------------------------------------------------------
Ci Inputs
Ci   categ,token,term
Ci   count: maximum number of elements to read (if count passed as <0
Ci          partok will read exactly -count elements or abort w/ error)
Ci   cast:  0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   istart:offset to first character of category where search begins
Ci   In the case of a character field, result should be the max size
Ci   optio: 0: show token that would have been sought
Ci          1: attempt to find token and read contents
Ci          2: print out token and contents
Co Outputs
Co   result:array into which elements are read
Co   noerr: see remarks
Co    iend: offset to first character past match
Co  partok: number of elements actually read (-1 if no attempt to read
Co          an element -- if optio is 0 or 2.)
Cr Remarks
Cr   on return, noerr is always false unless a token was actually
Cr   sought (optio=1) and found (note that this differs from getcat).
Cr   Passing common-block variable iend in place of istart has the
Cr   effect of updating istart to point beyond the end of the match.
Cr
Cr   ichoos is automatically incremented on each call if nchoos .ne. 0
Cr   ichoos and nchoos are automatically reset to 0 if a match is found
Cr   or if ichoos=nchoos
C ----------------------------------------------------------------
C     implicit none
      integer cast
      character*1 categ(0:*),term
      character*(*) chr,token
      integer count,result,istart
      logical parstr,cnvt,errflg,lopt

C local variables
      integer j,k
      character*1 strn(72)
      character*8 nmcast(0:4)
      character*11 option

C common block
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      data nmcast /'logical$','char$','integer$','real$','double$'/

      errflg = (count .lt. 0)
      partok = -1
      noerr = .false.
      if (nchoos .ne. 0) ichoos = ichoos+1
c apr9 2008
c      goto (1,2,3), 1+optio
      if(1+optio==1) goto 1
      if(1+optio==2) goto 2
      if(1+optio==3) goto 3

C --- Print message showing token sought, its cast and length ---
    1 continue
      option = ' (optional)'
      if (lopt) option = ' '
      if (nchoos .ne. 0 .and. ichoos .ne. nchoos) option = '    --- OR:'
      j = 6
      if (count .eq. 0) j=1
      call msgio(strn,' ',token,term,1,j,.false.,k)
      if (count .eq. 0) then
        j = k
      else
        call strcop(strn(k+1),nmcast(cast),9,'$',j)
        j = j+k-1
      endif
      if (iabs(count) .le. 1) then
        print *, (strn(k), k=1,j), option
      else
        call strcop(strn(j+1),' and length $',100,'$',k)
        j = j+k-1
        print *, (strn(k), k=1,j), iabs(count), option
      endif
      goto 42

C --- print out the contents of token ---
    3 if (nchoos .eq. 0 .or. ichoos .eq. 1 .and. count.ne.0) call disply
     .  (token,term,result,result,result,result,chr,iabs(count),cast)
      goto 42

C --- seek token within a category ---
    2 iend = istart
      partok = 0
      if (.not. parstr(categ,token,subsiz,maxlen,term,iend,j)) then
        if (ichoos .lt. nchoos .or. .not. lopt) goto 40
        call msgio(strn,categ,token,term,2,4,.true.,j)
Cstop2rx 2013.08.09 kino        stop
        call rx( '')
      endif
      noerr = .true.

      partok = -1
   10 partok = partok + 1
      if (partok .eq. iabs(count)) goto 40
      call skipbl(categ,subsiz,j)
      if (cast .eq. 1) then
        chr = ' '
        call strcop(chr,categ(j),result,' ',k)
        j = j + k
        goto 10
      else
        if (cnvt(categ,result,result,result,result,
     .             cast,partok,' ',j)) goto 10
      endif

      if (partok .lt. iabs(count) .and. errflg) then
        call msgio(strn,categ,token,term,2,1,.true.,k)
        print 22, iabs(count), partok
   22   format(' sought',i3,' elements but found only',i3)
Cstop2rx 2013.08.09 kino        stop
        call rx( '')
      endif

   40 iend = j
   42 if (noerr .or. ichoos .eq. nchoos) then
        nchoos = 0
        ichoos = 0
      endif
      return
      end
      logical function scat(ifi,categ,term,lrewnd)
C- scan file for a category
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi,categ,term
Ci   lrewnd: true, rewind file before searching; false, do not
Co Outputs
Co   true if category found, false if not
Cr Remarks
Cr   This is a file version of subroutine getcat
C ----------------------------------------------------------------
      integer ifi
      logical lrewnd
      character categ,term
      character*72 a
C Local variables:
      integer catl0,recl0
      parameter (catl0=7,recl0=72)
      integer i
      logical rdstrn,lsequ
      external rdstrn,lsequ
      if (lrewnd) rewind ifi
      scat = .false.
   10 if (.not. rdstrn(ifi,a,recl0,.false.)) return
      if (.not. lsequ(categ,a,catl0,term,i)) goto 10
      scat = .true.
      return
      end
      subroutine subcat(categ,token,term,i)
C- find a subcategory within a category
C ----------------------------------------------------------------
Ci Inputs
Ci   categ,token,term
Ci   i:     offset to first character of category where search begins
Co Outputs
Co   subsiz set to less than range of new token
Cr Remarks
Cr   noerr always returns unchanged.  Does nothing unless optio=1.
C ----------------------------------------------------------------
C passed variables
      integer i
      character*1 categ(0:*),term
      character*(*) token
      logical parstr

C local variables
      integer i2,j

C Common block
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

c apr9 2008
c      goto (1,2,1), 1+optio
      if(1+optio==1) goto 1
      if(1+optio==2) goto 2
      if(1+optio==3) goto 1

C --- seek token within a category ---
    2 subsiz = catsiz
      i2 = i
      if (parstr(categ,token,subsiz,maxlen,term,i2,j)) subsiz = i2

    1 return
      end
      logical function parstr(s1,s2,reclen,len,term,i,j)
C- find a substring within a given string
C ----------------------------------------------------------------
Ci Inputs
Ci   i: character where search within string should begin
Co Outputs
Co   i: index to first position of token
Co   j: index to first position after token
Cr Remarks
Cr    seeks match at string(i), string(i+1), ... string(i+reclen) until
Cr    match is found.  returns false if no match found.
C ----------------------------------------------------------------
      integer reclen,len,i,j
      character*1 s1(0:*),s2(0:*),term
      logical lsequ

      parstr = .false.
      i = i-1
   10 i = i+1
      if (i .eq. reclen) return
      if (.not. lsequ(s1(i),s2,len,term,j)) goto 10
      parstr = .true.
      j = j + i
      return
      end
      subroutine chrpos(s,ch,maxch,ich)
C- Finds position of character in string
C ----------------------------------------------------------------
Ci Inputs
Ci   s:   string (declared as s(0:*)
Ci   ch:  character sought
Ci   ich: start search at s(ich)
ci   maxch: see ich
Co Outputs
Co   ich: position of character ch, not to exceed maxch
Cr Remarks
Cr    seeks match at string(i0), string(i0+1) until ch is found or until
Cr    ich = maxch.
C ----------------------------------------------------------------
      integer ich,maxch
      character*1 ch,s(0:*)

   10 if (ich .eq. maxch  .or.  s(ich) .eq. ch) return
      ich = ich+1
      goto 10
      end
      subroutine skipbl(t,nt,i)
C- Parses string T(I) for blanks
      integer nt,i
      character*1 t(0:nt)
   99 if (t(i) .ne. ' ') return
      i = i + 1
      if (i .ge. nt) return
      goto 99
      end
      subroutine tokmat(string,token,n,len,term,itoken,lopt)
C- compare a string to a list of strings
C ----------------------------------------------------------------
Ci Inputs
Ci   string: test string
Ci   token: vector of strings to compare
Ci   n,len: number and length of strings in token
Ci   lopt:  if true, tokmat stops with error message when no match found
Co Outputs
Co   itoken: index to string tokened, -1 if none
C ----------------------------------------------------------------
      integer n,len,itoken
      logical lsequ,lopt
      character*1 string(*), term
      character*(*) token(0:*)
      integer i
      do  10  itoken = 0, n-1
        if (lsequ(string,token(itoken),len,term,i)) return
   10 continue
      itoken = -1
      if (lopt) then
        print *, 'TOKMAT: unmatched ', (string(i), i=1,len)
Cstop2rx 2013.08.09 kino        stop
        call rx( '')
      endif
      return
      end
      subroutine msgio(strng,categ,token,term,m1,m2,lopt,i)
C- Messages for getcat and partok
C ----------------------------------------------------------------
Ci Inputs
Ci   strng,categ,token,term,m1,m2,lopt,i
Co Outputs
Co   string is printed to std out
Co   i: length of printed string
Cr Remarks
Cr   String is composed as: msg(m1) msg(m2)
Cr   If the first character in <categ> is not blank, the string
Cr   'category <categ>' is appended to msg(m1)
Cr   If the first character in <token> is not blank, the string
Cr   '  token <token>' is prepended to msg(m2)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*1 strng(72)
      character*1 categ(*), token(*), term
      integer i,m1,m2
      logical lopt
C Local parameters
      integer imsg(11),k
      character*1 msgs(1)
      character*70 msg2
      equivalence (msg2,msgs)
      data imsg /1,2,17,28,41,52,1,1,1,1,1/
      data msg2
     ./'$read error in $unmatched $  is missing$(optional)$  of cast $'/
      call strcop(strng,msgs(imsg(m1)),72,'$',i)

      if (categ(1) .ne. ' ') then
        call strcat(strng,100,'$','category: $',100,'$',i)
        call strcat(strng,100,'$',categ,100,' ',i)
      endif
      if (token(1) .ne. ' ') then
        strng(i) = '$'
        call strcat(strng,100,'$','  token  $',100,'$',i)
        call strcat(strng,100,'$',token,100,term,i)
      endif
      call strcop(strng(i+1),msgs(imsg(m2)),72,'$',k)
      i = i+k-1
      if (lopt) print 10, (strng(k), k=1,i)
   10 format(' ',72a1)
      return
      end
      logical function cnvt(instr,resL,resI,resR,resD,cast,count,term,j)
C- ASCII conversion to logical, integer, real, double
C ----------------------------------------------------------------
Ci Inputs
Ci   cast:  0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   j:     offset to first character in string to read
Co Outputs
Co   count'th element of resL,resI,resR, or resD is converted
Co   j is put past last character read
Cr Remarks
Cr   An ascii string is converted to logical, integer, real or double
Cr   using the FORTRAN internal read facility.  An unformatted read
Cr   is preferable, but does not conform to the ANSII 77 standard.
Cr
Cr   IOLIB only uses cnvt as a nucleus for partok.
C ----------------------------------------------------------------
      character*1 instr(0:1),term
      integer cast,count
C these are all equivalent address spaces:
      double precision resD(0:1)
      real resR(0:1)
      integer resI(0:1)
      logical resL(0:1)

      integer k,maxsiz,j
      parameter (maxsiz=15)
      character*(maxsiz) strn
      character*1 strn2(maxsiz),cj
      equivalence (strn,strn2)

C --- Early error checking to avoid problems on some compilers ---
      cnvt = .false.
      if (instr(j) .eq. ' ') call skipbl(instr(j),100,j)
      cj = instr(j)
      if (cast .ne. 0 .and.
     .    .not. (cj .eq. '+' .or.  cj .eq. '-' .or. cj .eq. '.' .or.
     .           cj .ge. '0' .and. cj .le. '9')) return

      strn = ' '
      call strcop(strn,instr(j),maxsiz,term,k)
C Right justify input (required for some compilers)
      if (k .le. maxsiz) then
        strn = ' '
        call strcop(strn2(maxsiz+2-k),instr(j),k-1,term,k)
      endif
      j = j + k + 1

c apr9 2008
c      goto (1,2,3,4,5), cast+1
      if(cast+1==1) goto 1
      if(cast+1==2) goto 2
      if(cast+1==3) goto 3
      if(cast+1==4) goto 4
      if(cast+1==5) goto 5


C#ifdef FORMATTED_READ
    1 read(strn,'(L15)',err=20) resL(count)
C     print *, 'logical string ', strn, 'converted to ', resL(count)
    2 goto 10
    3 read(strn,'(I15)',err=20) resI(count)
C     print *, 'integer string ', strn, 'converted to ', resI(count)
      goto 10
    4 read(strn,'(E15.0)',err=20) resR(count)
C     print *, 'real string ', strn, 'converted to ', resR(count)
      goto 10
    5 continue
      read(strn,'(E15.0)',err=20) resD(count)
C     print *, 'real*8 string ', strn, 'converted to ', resD(count)
      goto 10
C#elseC
C    1 read(strn,*,err=20) resL(count)
CC     print *, 'logical string ', strn, 'converted to ', resL(count)
C    2 goto 10
C    3 read(strn,*,err=20) resI(count)
CC     print *, 'integer string ', strn, 'converted to ', resI(count)
C      goto 10
C    4 read(strn,*,err=20) resR(count)
CC     print *, 'real string ', strn, 'converted to ', resR(count)
C      goto 10
C    5 continue
C      read(strn,*,err=20) resD(count)
CC     print *, 'real*8 string ', strn, 'converted to ', resD(count)
C      goto 10
C#endif
   10 cnvt = .true.
   20 return
      end
      subroutine disply(token,term,resL,resI,resR,resD,chr,count,cast)
C- Called as a nucleus only by partok
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   Used because fortran cannot treat same space as different casts.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer count,cast
C#ifdefC APOLLO_BUG
C      character*(*) token
C      character*4 chr(0:*)
C#else
      character*(*) token,chr(0:*)
C#endif
      character term
C these are all equivalent address spaces:
      double precision resD(0:count)
      real resR(0:count)
      integer resI(0:count)
      logical resL(0:count)
C Local variables:
      integer catl0
      parameter (catl0=7)
      integer i
      character*(catl0) strn

      strn = ' '
      call strcop(strn,token,catl0,term,i)

c apr9 2008
c      goto (1,2,3,4,5), cast+1
      if(cast+1==1) goto 1
      if(cast+1==2) goto 2
      if(cast+1==3) goto 3
      if(cast+1==4) goto 4
      if(cast+1==5) goto 5

      print 10,strn
      return
    1 print 10,strn,(resL(i), i=0,count-1)
   10 format(3X,1P,a8,1x,20(20(L1,1x)/12x))
      return
    2 print *,strn,(chr(i), i=0,count-1)
      return
    3 print 30,strn,(resI(i), i=0,count-1)
   30 format(3X,1P,a8,1x,20(10I5/12x))
      return
    4 print 40,strn,(resR(i), i=0,count-1)
   40 format(3X,1P,a8,1x,20(5G14.7/12x))
      return
    5 print 50,strn,(resD(i), i=0,count-1)
   50 format(3X,1P,a8,1x,20(3G15.8/12x))
      return
      end
      logical function rdstrn(unit,a,len,lopt)
      integer unit,len
      logical lopt
      character*1 a(len)

      rdstrn = .true.
      read(unit,10,end=20) a
      if (lopt) print 10, a
      return
   10 format(72a1)
   20 rdstrn = .false.
      return
      end
      logical function lsequ(s1,s2,len,term,i)
C- Determine whether two strings are equal or not
C ----------------------------------------------------------------
Ci Inputs
Ci   s1,s2: strings to compare
Ci   len:   maximum length of string
Ci   term:  terminator
Co Outputs
Co   lsequ: returned true or false
Co   i:     number of characters tokened (including terminator)
Cr Remarks
Cr   string comparison continues until terminator encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
      integer i,len
      character*1 s1(len),s2(len),term
c     lsequ = .true.
c     do  10  i = 1, len
c       lsequ = (lsequ .and. s1(i) .eq. s2(i))
c       if (s1(i) .eq. term) return
c  10 continue
      lsequ = .false.
      do  10  i = 1, len
        if (s1(i) .ne. s2(i)) return
        if (s1(i) .eq. term) goto 15
   10 continue
   15 lsequ = .true.
      return
      end
      subroutine strcop(dest,source,len,term,i)
C- copy one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   dest,source: source and destination strings
Ci   len:   maximum number of characters to copy
Ci   term:  terminator
Co   i:     number of characters copied (including terminator)
Cr Remarks
Cr   string copy continues until term encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
      integer len
      character*1 dest(len),source(len),term
      integer i
      if (len .eq. 0) return
      i = 0
   10 i = i+1
      dest(i) = source(i)
      if (dest(i) .ne. term .and. i .lt. len) goto 10
      end
      subroutine strcat(s1,len1,term1,s2,len2,term2,i)
C- concatenate one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   s1,s2: source string and string to concatenate
Ci   len1:  maximum length of s1
Ci   term1: terminator for s1
Ci   len2:  maximum length of s2
Ci   term2: terminator for s2
Co Outputs
Co   i:     number of characters in s1 (including terminator)
Cr Remarks
Cr   concatenation continues until term encountered or
Cr   len2 characters are concatenated.
C ----------------------------------------------------------------
      integer len1,len2
      character*1 s1(len1),s2(len2),term1,term2
      integer i,j

      i = 0
   10 i = i+1
      if (s1(i) .ne. term1 .and. i .lt. len1) goto 10
      j = 0
      if (s1(i) .eq. term1) i = i-1
   20 j = j+1
      i = i+1
      s1(i) = s2(j)
      if (s2(j) .ne. term2 .and. j .lt. len2) goto 20
      end
      subroutine strcat2(s1,s2,len,i)
C- concatenate one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   s2,s1: source string and string to concatenate
Ci   len :  number of characters in s1 and s2
Ci   i    : last len - i cahracters of s2 copied
Ci          into first len - i char. in s1
Co Outputs
Cr Remarks
C ----------------------------------------------------------------
      integer len
      character*1 s1(len),s2(len)
      integer i,j

      do 10 j = 1,len
        if( j .gt. len-i ) then
          s1(j) = ' '
        else
          s1(j) = s2(j+i)
        endif
  10  continue
      return
      end
C --- FILE HANDLING ---
      integer function fopna(nam,unit,switch)
C- File opening
C ----------------------------------------------------------------
Ci Inputs
Ci    nam:  LMTO file name
Ci   unit:  logical unit
Ci   switch:switch governing mode of file opening
Co Outputs
Co   fopna returns logical unit number for file name
Co     (file is opened if it is not already open)
Cr Remarks
Cr   fopna(name,integer unit,integer switch): general file opening:
Cr      logical unit is passed as part of the call.
Cr      Returns logical unit.
Cr   fadd(name,unit,switch): sends a new name to append to a
Cr     list of known file names and associated logical units.
Cr     Initializer for fopn, fopno and fopno, routines which return
Cr     logical units for files from the name.
Cr     Returns 1 + logical unit.
Cr   fopn(name) file opening, for which a logical unit and switch has
Cr     already been given by a call to fadd.  Same as fopna, except
Cr     that logical unit and switches have been previously stored.
Cr     Returns logical unit.
Cr   fopnn(name) same as fopn, but open with status='NEW'
Cr   fopno(name) same as fopn, but open with status='OLD'.
Cr   fxst(name) same as fopn, but inquires as to whether file exists
Cr   fext(ext) changes the extension from the default ".dat".  Valid
Cr     only for operating systems which permit extensions.
Cr   fhndl(name) returns with logical unit associated with name,
Cr     and -1 if none exists.
Cr
Cr   Switch is a composite of integers a, b, c, d, stored as digits
Cr   abcd, base 2.
Cr     Bits 0,1: = 0, open the file as 'UNKNOWN'
Cr               = 1, open the file as 'OLD'
Cr               = 2, open the file as 'NEW'
Cr     Bit 2:    if set, open file as unformatted.
Cr     Thus switch=5 opens file unformatted, status='old'
C ----------------------------------------------------------------
C Passed parameters
      character*(*) nam
      integer switch,unit
C Local parameters
      integer fext,fopn,fopno,fopnn,fhndl,fadd,fxst,
     . i,ista,iunit,isw,n,bit,iprint,i1mach
      character*11 ftnfmt
      character*8 fnam,ftnsta
      logical bittst,jsopen,ldum
      integer mxnam
      parameter (mxnam=12)
      character*4 pasnam(mxnam)
      integer passw(mxnam), pasuni(mxnam), nnam
      save pasnam, passw, pasuni, nnam

C --- If an extension is allowed or desired ... ---
C#ifdefC EXTDOTDAT
C      integer extlen
C      parameter (extlen=3)
C      character*(extlen+1) ext
C      save ext
C#endif

      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)

      data nnam /0/
C --- Extension for all file names, if desired ---
C#ifdefC EXTDOTDAT
C      data ext /'.dat'/
C#endif

C --- Purge any leading blanks ---
      fnam  = nam
      i = 0
      call skipbl(fnam,4,i)
C#ifdefC CRAY
C        if ( i .ne. 0 ) then
C        call strcat2(fnam,fnam,4,i)
C        endif
C#else
      fnam = fnam(1+i:4)
C#endif

C#ifdefC unix
C      call str_lo_case(fnam)
C#endif
      iunit = unit
      isw   = switch
      ista = 0
      goto 20

      entry fopn(nam)
      ista = 0
      goto 2

      entry fopno(nam)
      ista = 1
      goto 2

      entry fopnn(nam)
      ista = 2
      goto 2

      entry fxst(nam)
      ista = -1
      goto 2

      entry fhndl(nam)
      ista = -2
      fhndl = -1
      goto 2

      entry fadd(nam,unit,switch)
      fadd = unit+1
      nnam = nnam+1
Cstop2rx 2013.08.09 kino      if (nnam .gt. mxnam) stop 'iolib: too many file names'
      if (nnam .gt. mxnam) call rx( 'iolib: too many file names')
      pasnam(nnam) = nam
C#ifdefC unix
C      call str_lo_case(pasnam(nnam))
C#endif
      pasuni(nnam) = unit
      passw(nnam)  = switch
      return

      entry fext(nam)
C#ifdefC EXTDOTDAT
C      ext = '.'//nam(1:extlen)
C      fext = 0
C#endif
      return

    2 continue
      fnam = nam
C#ifdefC unix
C      call str_lo_case(fnam)
C#endif
      do  10  i = 1, nnam
        if (fnam .eq. pasnam(i)) then
c          fnam = pasnam(i)
          iunit = pasuni(i)
          fhndl = iunit
          isw = passw(i) - mod(passw(i),4) + ista
          goto 20
        endif
c        PRINT *, 'dbg fopn name check', i, nam, pasnam(i)
   10 continue
      if (ista .eq. -2) return
      write(i1mach(4),*) 'fopn: name mismatch: file ', fnam
Cstop2rx 2013.08.09 kino      stop
      call rx( '')

   20 continue
      if (ista .eq. -2) return

      ftnfmt = 'FORMATTED'
      if (bittst(isw,4)) ftnfmt = 'UNFORMATTED'

C --- Attach an extension ---
C#ifdefC EXTDOTDAT
C      call strcat(fnam,4,' ',ext,4,' ',i)
C#endif

C --- Handle INQUIRE statements ---
      if (ista .lt. 0) then
        fxst = 0
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        INQUIRE(UNIT=iunit,EXIST=ldum)
C#else
        INQUIRE(FILE=fnam,EXIST=ldum)
C#endif
        if (ldum) fxst = 1
        return
      endif

      ftnsta = 'UNKNOWN'
      if (mod(isw,4) .eq. 1) ftnsta = 'OLD'
      if (mod(isw,4) .eq. 2) then
        ftnsta = 'NEW'
        call fclose(iunit)
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        open(iunit,FORM=ftnfmt,STATUS='UNKNOWN')
C#else
        open(iunit,FILE=fnam,FORM=ftnfmt,STATUS='UNKNOWN')
C#endif
        close(unit=iunit,status='DELETE')
      endif

c      print *, 'dbg: fopn',fnam,iunit,isw,ista,ftnfmt,ftnsta

      fopna = iunit
      fopn = iunit
      fopno = iunit
      fopnn = iunit

      if (.not. jsopen(iunit)) then
        if (iprint() .ge. 100) print 300, ftnfmt,fnam,ftnsta,iunit
  300   format(/' FOPEN: opening',a12,' file ',a9,' status=',a7,
     .          ', unit=',i2)
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        open(iunit,FORM=ftnfmt,STATUS=ftnsta)
C#else
C#ifdefC APOLLO_PRE_OS10.1
C        if (ftnfmt .eq. 'FORMATTED') then
C          open(iunit,FILE=fnam,FORM=ftnfmt,STATUS=ftnsta)
C        else
C          open(iunit,FILE=fnam,FORM=ftnfmt,STATUS=ftnsta,RECL=9999999)
C        endif
C#else
        open(iunit,FILE=fnam,FORM=ftnfmt,STATUS=ftnsta)
C#endif
C#endif
      endif
      end
      logical function jsopen(unit)
C- adds unit to list of open files, returns whether or not aleady open
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   jsopen
Cr Remarks
Cr   maxfil is maximum number of files that may be open at one time
Cr   minuni is the lowest logical unit number
Cr   All file closings should be made through entry fclose
C ----------------------------------------------------------------
C     implicit none
      integer unit
      integer i,iprint
      integer maxfil,minuni
      parameter (maxfil=15,minuni=10)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      if (iprint() .ge. 102) print 20, unit, (unitab(i), i=0, nopen-1)
   20 format(/' JSOPEN: check logical unit',i3,
     .        ' among open units: ',15i3)

      jsopen = .true.
      do  10  i = nopen-1, 0, -1
        if (unitab(i) .eq. unit) return
   10 continue

      jsopen = .false.
      unitab(nopen) = unit
      nopen = nopen+1
Cstop2rx 2013.08.09 kino      if (nopen .gt. maxfil) stop 'JSOPEN: too many files'
      if (nopen .gt. maxfil) call rx( 'JSOPEN: too many files')
      return
      end
      subroutine fclose(unit)
C- closes an open file, removing unit from the stack
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   none
Cr Remarks
Cr   use in conjunction with jsopen
C ----------------------------------------------------------------
c      implicit none
      integer unit
      integer i,iprint
      character*6 clstat
      integer maxfil,minuni
      parameter (maxfil=15,minuni=10)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      clstat = 'KEEP'
      goto 10

   10 continue
      if (iprint() .ge. 100) print 20, unit
   20 format(/' FCLOSE: closing',i3)

      do  30  i = nopen-1, 0, -1
        if (unitab(i) .eq. unit) unitab(i) = 999
   30 continue
      call ishell(nopen,unitab)
      if (nopen .eq. 0 .or. unitab(max(nopen-1,0)) .ne. 999) then
        if (iprint() .ge. 100)
     .  print *, 'FCLOSE: attempt to close unopened file'
      else
        close(unit=unit,status=clstat)
        nopen = nopen-1
      endif
      return
      entry dfclos(unit)
      clstat = 'DELETE'
      goto 10
      end
      integer function lgunit(i)
C- Returns stdout for i=1, log for i=2
      integer i, fopno, i1mach
      external fopno,i1mach

      lgunit = i1mach(2)
      if (i .eq. 2) lgunit = fopno('LOG')
      end
      subroutine headl2(name,ifi)
C-  Puts a heading line into file, unit ifi
      integer ifi
      character*8 name
      write(ifi,300) name
  300 format(' -------------------------  START ',a8,
     .       ' -------------------------')
      end
      subroutine poseof(iunit)
C- Positions file at end-of-file
C Passed parameters
      integer iunit
C Local parameters
      integer i,nrec

      nrec=0
      rewind iunit
      do 10 i=1,10000
        read(iunit,100,end=90,err=91)
        nrec=i
  100   format(a1)
   10 continue
      write(*,200) iunit
  200 format(' POSEOF: no EOF found for file',i3)
      return
   90 continue
      rewind iunit
      do 11 i=1,nrec
        read(iunit,100)
   11 continue
   91 continue
      end
C --- Virtual memory routines ---
c$$$      subroutine vmem(o1,o2)
c$$$C- Virtual memory routines
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   o1,o2
c$$$Co Outputs
c$$$Co   Nothing
c$$$Cr Remarks
c$$$Cr
c$$$C ----------------------------------------------------------------
c$$$C Passed parameters
c$$$      integer o1,o2
c$$$C Local parameters
c$$$      logical lsave
c$$$      integer len,offst,fopn,fopno,ifi,i,ii
c$$$      save len,lsave,offst
c$$$C heap:
c$$$      integer w(1)
c$$$      common /w/ w
c$$$
c$$$      offst = o1
c$$$      len = -o2
c$$$      lsave = (o2 .lt. 0)
c$$$      if (o2 .lt. 0) o2 = o1
c$$$      return
c$$$
c$$$      entry vmems(i)
c$$$      if (.not. lsave) return
c$$$      ifi = fopn('TMP')
c$$$      rewind ifi
c$$$      do  10  ii = 1, i
c$$$   10 call vmem2(-ifi,w(offst),len)
c$$$      return
c$$$
c$$$      entry vmemg(i)
c$$$      if (.not. lsave) return
c$$$      ifi = fopno('TMP')
c$$$      rewind ifi
c$$$      do  20  ii = 1, i
c$$$   20 call vmem2(ifi,w(offst),len)
c$$$      call fclose(ifi)
c$$$
c$$$      end
c$$$      subroutine vmem2(ifi,w,len)
c$$$      integer ifi,len,w(len)
c$$$
c$$$      if (ifi .gt. 0) read(ifi) w
c$$$      if (ifi .lt. 0) write(-ifi) w
c$$$      end
C#ifdefC FUNIT
C      subroutine funit(unit)
CC- returns next available logical unit for file opening
CC ----------------------------------------------------------------
CCo Outputs
CCo   unit
CCr Remarks
CCr   maxfil is maximum number of files that may be open at one time
CCr   minuni is the lowest logical unit number
CCr   All file closings should be made through entry fclose
CC ----------------------------------------------------------------
CC      implicit none
C      integer unit
C      integer maxfil,minuni
C      parameter (maxfil=15,minuni=10)
C      integer unitab(0:maxfil-1),nfiles
C      integer i
C      save
C      data nfiles /0/
C
C      do  10  i = 0, nfiles-1
C   10 if (unitab(i) .ne. i+minuni) goto 20
C      i = nfiles
C   20 unit = i + minuni
C      unitab(nfiles) = unit
C      nfiles = nfiles+1
C      call ishell(nfiles,unitab)
C      return
C
C      entry fclose(unit)
C      do  30  i = 0, nfiles-1
C   30 if (unitab(i) .eq. unit) goto 40
C      i = nfiles
C   40 continue
C      unitab(i) = 1000
C      call ishell(nfiles,unitab)
C      nfiles = nfiles-1
C
C      endfile unit
C      close(unit=unit)
C
C      return
C      end
C#endif  FUNIT
C      subroutine prtstr(unit,string,term,len)
C      integer unit,len
C      character*1 string(*),term
C      i = 0
C   88 i = i+1
C      if (string(i) .ne. term .and. i .le. len) goto 88
C      write(unit,891) (string(j), j=1,i-1)
C  891 format(72a1)
C      return
C      end
C      character*1 function chrint(i)
C      integer i
C      character*1 csym(10)
C      data csym/'0','1','2','3','4','5','6','7','8','9'/
C      chrint = csym(i)
C      return
C      end
C      logical function ischar(ch)
C      character*1 ch
C      ischar = (ch .ge. 'a' .and. ch .lt. 'z' .or.
C     .          ch .ge. 'A' .and. ch .lt. 'Z' .or.
C     .          ch .ge. '0' .and. ch .le. '9')
C      return
C      end
C --- Default initialization ---
      block data diolib
C     implicit none
      integer catl0,recl0
      parameter (catl0=7,recl0=72)

C for iprint...
c      integer nstack
c      parameter (nstack=5)
c      integer vstack(0:nstack-1),stackp
c      common /iprnt/ vstack,stackp

C For io routines ...
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

C for logical unit routines ...
      integer maxfil,minuni
      parameter (maxfil=15,minuni=10)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      data recoff /0/, reclen /recl0/, maxlen /catl0/,
     .     ichoos /0/, nchoos /0/
c      data vstack /30,30,30,30,30/, stackp /0/
      data nopen /0/

      end





