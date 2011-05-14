#define unix 1
#define PATH 1
c #define FEXTENS 1

ctakao's new file opener and closer for f90. This can be consistent with fopna
c      integer function open90(fnam,options)
c      character*(*)::fnam, options
c        inquire(fnum,number
c        ifio=9038
c        inquire(unit=ifio,opened=nexist)
c        if(nexist) 
c     &  stop 'bndfp: I assumed ifio=9038 is not used(too sloppy)'


      integer function fopna(fnam,lunit,ista)
C- File opening
C ----------------------------------------------------------------
Ci Inputs
Ci    fnam:  file name
Ci   lunit:  file logical lunit.  If passed lunit is -1, program finds
Ci           the next available logical unit; see Remarks.
Ci   ista:   status switches governing mode of file opening; see Remarks
Co Outputs
Co   fopna returns logical lunit number for file name
Cr Remarks
Cr   fopnx, below, is an all-purpose file-handling routine that
Cr   keeps an internal table that maintains an association of filenames
Cr   with their logical unit and status.  It is intended
Cr   to facilitate reference to files without carrying around
Cr   the logical unit number or file status.  Several 'front end'
Cr   routines listed here call fopnx with special arguments.
Cr
Cr   fopng(fnam,lunit,ista):
Cr     opens file 'fnam' with attributes set by 'ista',
Cr     assigns logical unit 'lunit', and returns the logical unit.
Cr     If lunit<0, the logical unit will be assigned if there is no
Cr     record of this file name or logical unit in the internal table
Cr     maintained by fopnx. 
Cr     In this case the file name and lunit are added to the table.
Cr     If fopnx has a record of the logical unit, that unit and
Cr     the corresponding filename are used for file handling.  Else:
Cr     If fopnx has a record of the file name, that name and
Cr     corresponding lunit are used for file handling.
Cr
Cr     ista is a composite of integers a, b, c, d, stored as digits
Cr     abcd, base 2.
Cr     Bits 01: 0, open the file as 'UNKNOWN'
Cr              1, open the file as 'OLD'
Cr              2, open the file as 'NEW'
Cr     Bit 2:   1, open file as unformatted.
Cr     Bit 3:   1, do not add this file to the file table
Cr     Thus ista=5 opens file unformatted, status='old'
Cr
Cr   fopna(fnam,lunit,ista):
Cr     is identical to fopng, except that 'fnam' a path is prepended
Cr     to, and extension is appended to 'fnam'.  'path' and 'ext'
Cr     are set and recovered with:
Cr
Cr   fext(ext)     sets extension
Cr   fpath(path)   sets path.  NB: path should end with '/'
Cr   fextg(ext)    recovers extension
Cr   fpathg(ext)   recovers path
Cr
Cr   fadd(fnam,lunit,ista) adds 'fnam' to the internal table.
Cr
Cr   fopn(fnam)  is equivalent to fopna(fnam,0,-1)
Cr
Cr   fopnn(fnam)   same as fopn, but open with status='NEW'
Cr   fopno(fnam)   same as fopn, but open with status='OLD'.
Cr
Cr   fxst(fnam)    inquires as to whether file exists
Cr
Cr   ffnam(fnam,fnam) returns full file name
Cr
Cr   fhndl(fnam) returns with logical lunit associated with name,
Cr     and -1 if none exists.
Cr
Cr   fclr(fnam,unit) clears the file name from the internal table,
Cr     and closes the file if open.
Cr
Cr   ftflsh(unit) flushes output buffer of logical unit.
Cr     This routine is nonstandard, and works only some some machines.
Cr     unit = stdout calls a C-language routine that flushes stdout
Cr     unit < 0 calls a C-language routine that flushes all buffers
Cr
Cr   fshow lists the internal table
Cu Updates
Cu   15 Feb 02 (ATP) added MPI log for lgunit(3)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) fnam
      integer ista,lunit
C Local parameters
      integer fext,fextg,fopn,fopno,fopnn,fhndl,
     .fadd,fxst,fpath,fpathg,fopng,fopnx
      integer mode
#ifndef AIX
      mode = 100 + 70*mod(ista/8,2)
      fopna = fopnx(fnam,mode,mod(ista,8),lunit)
      return

      entry fopng(fnam,lunit,ista)
      mode =   2 + 70*mod(ista/8,2)
      fopng = fopnx(fnam,mode,mod(ista,8),lunit)
      return

      entry fadd(fnam,lunit,ista)
      mode = 114 + 70*mod(ista/8,2)
      if (lunit .lt. 0) mode = mode-4
      fadd = fopnx(fnam,mode,mod(ista,8),lunit)
      return

      entry fopn(fnam)
      fopn = fopnx(fnam,100,0,-1)
      return

      entry fopno(fnam)
      fopno = fopnx(fnam,100,8+1,-1)
      return

      entry fopnn(fnam)
      fopnn = fopnx(fnam,100,8+2,-1)
      return

      entry fxst(fnam)
      fxst = fopnx(fnam,170,-1,-1)
      return

      entry fhndl(fnam)
      fhndl = fopnx(fnam,184,-1,-1)
      return

      entry fext(fnam)
      fext = fopnx(fnam,120,0,-1)
      return

      entry fpath(fnam)
      fpath = fopnx(fnam,140,0,-1)
      return

      entry fextg(fnam)
      fextg = fopnx(fnam,130,0,-1)
      return

      entry fpathg(fnam)
      fpathg = fopnx(fnam,150,0,-1)
      return
      end

      subroutine fclr(fnam,lunit)
C     implicit none
C Passed parameters
      character*(*) fnam
      integer lunit,fopnx,i

      i = fopnx(fnam,111,0,lunit)
      return

      entry fshow
      i = fopnx(' ',190,0,-1)
      end
      subroutine ffnam(fnam,filnam)
C- Appends extension to fnam
C     implicit none
      character*(*) fnam,filnam
      integer i,fopnx
      filnam = fnam
      i = fopnx(filnam,160,0,-1)
#else
      mode = 100 + 70*mod(ista/8,2)
      fopna = fopnx(fnam,mode,mod(ista,8),lunit)
#endif
      end

#ifdef AIX
      integer function fopng(fnam,lunit,ista)
C      implicit none
      character*(*) fnam
      integer mode,lunit,fopnx,ista
      mode =   2 + 70*mod(ista/8,2)
      fopng = fopnx(fnam,mode,mod(ista,8),lunit)
      end
      integer function fadd(fnam,lunit,ista)
C      implicit none
      character*(*) fnam
      integer ista,lunit
      integer fopnx,mode
      mode = 114 + 70*mod(ista/8,2)
      if (lunit .lt. 0) mode = mode-4
      fadd = fopnx(fnam,mode,mod(ista,8),lunit)
      end
      integer function fopn(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fopn = fopnx(fnam,100,0,-1)
      end
      integer function fopno(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fopno = fopnx(fnam,100,8+1,-1)
      end
      integer function fopnn(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fopnn = fopnx(fnam,100,8+2,-1)
      end
      integer function fxst(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fxst = fopnx(fnam,170,-1,-1)
      end
      integer function fhndl(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fhndl = fopnx(fnam,184,-1,-1)
      end
      integer function fext(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fext = fopnx(fnam,120,0,-1)
      end
      integer function fpath(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fpath = fopnx(fnam,140,0,-1)
      end
      integer function fextg(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fextg = fopnx(fnam,130,0,-1)
      end
      integer function fpathg(fnam)
C      implicit none
      character*(*) fnam
      integer fopnx
      fpathg = fopnx(fnam,150,0,-1)
      end
      subroutine fclr(fnam,lunit)
C      implicit none
C Passed parameters
      character*(*) fnam
      integer lunit,fopnx,i

      i = fopnx(fnam,111,0,lunit)
      return

      entry fshow
      i = fopnx(' ',190,0,-1)
      end
      subroutine ffnam(fnam,filnam)
C- Appends extension to fnam
C      implicit none
      character*(*) fnam,filnam
      integer i,fopnx
      filnam = fnam
      i = fopnx(filnam,160,0,-1)
      end
#endif
      integer function fopnT(fnam,lunit,ista,mode)
C- File opening, using a separate directory
C ----------------------------------------------------------------------
Ci Inputs
Ci   fnam  :file name WITHOUT path or extension
Ci   lunit :file logical lunit.  If passed lunit is -1, program finds
Ci         :the next available logical unit; see Remarks in fopna
Ci   ista  :status switches governing mode of file opening;
Ci         :see Remarks in fopna
Ci   mode  :selects which path is used for file opening
Ci         :0 fopnT opens file using tmpdir for path.
Ci         :  The extension used by fopna is appended
Ci         :1 fopnT opens file using tmpdir for path.
Ci         :  No extension is appended.
Ci         :..Add 2 to cause fopnT returns logical unit associated 
Ci         :  with name if it exists in table, but do not touch file
Ci         :11 fopnT opens no files, but copies fnam into tmpdir
Ci         :  In this case, lunit and ista are not used.
Ci         :12 fopnT opens no files, but gets environment variable
Ci         :  TMPDIR and copies it into tmpdir
Ci         :  In this case, fnam, lunit and ista are not used.
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   30 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character fnam*(*)
      integer ista,lunit,mode
C ... Local parameters
      integer j1,j2,fextg,fopng,fopnx
      character tmpdir*100,pfnam*100,ext*100

      logical:: isanrg, l_dummy_isanrg

      save tmpdir
      data tmpdir /' '/

      if (mode .eq. 12) then
        call gtenv('TMPDIR',tmpdir)
      elseif (mode .eq. 11) then
        tmpdir = fnam
      endif
      if (mode .eq. 11 .or. mode .eq. 12) then
        call word(tmpdir,1,j1,j2)
        if (j2 .ge. j1) then
          if (tmpdir(j2:j2) .ne. '/') tmpdir(j2+1:j2+1) = '/'
        endif
        fopnT = 0
        return
      endif

Ckino isanrg is logical function,       call isanrg(mode,0,4,'fopnt','mode',.true.)
      l_dummy_isanrg=isanrg(mode,0,4,'fopnt','mode',.true.)

      if (mod(mode,2) .eq. 0) then
        j1 = fextg(ext)
        call namcat(pfnam,fnam,ext,.false.)
        ext = pfnam
        call namcat(pfnam,tmpdir,ext,.true.)
      else
        call namcat(pfnam,tmpdir,fnam,.true.)
      endif

      if (mode .le. 1) then
        fopnT = fopng(pfnam,lunit,ista)
      else
        fopnT = fopnx(pfnam,186,-1,-1)
      endif
      end
      integer function fopnx(pnam,mode,ista,lunit)
C- Kernel called by fopna, etc to open files
C ----------------------------------------------------------------
Ci Inputs
Ci   pnam:   file name, or portion of file name, depending on mode
Ci   mode:   1s digit (for 10s digit 0)
Ci           0  Find lunit, ista from name.
Ci              pnam is name sans path, extension
Ci              If name does not exist, add to table, assign lunit.
Ci              Subsequent file handling depends on ista.
Ci           1  Remove name from internal table.  Close file if open.
Ci              Error if name is missing from table.
Ci          ... Add 2 if pnam is already full name
Ci          ... Add 4 to require lunit match internal table
Ci                    returns -1 if lgunit<0 and pnam missing from table
Ci          10s digit: for modifying or retrieving file's path and extension
Ci           ... 1..6 If resulting file name (after addition of path and
Ci                     extension) is new, it is added to internal table.
Ci           ... 1..6,8,9 suppress file handling
Ci           1  Suppress any file handling
Ci           2  copy pnam to extension
Ci           3  copy extension into pnam
Ci           4  copy pnam to path
Ci           5  copy path into pnam
Ci           6  copy full file name into pnam.
Ci           7  Do not add file name to internal table, but allow
Ci              removal of name from table
Ci           8  Do not add file name to internal table
Ci           9  Show internal tables; no file handling
Ci          100s digit
Ci           1  convert pnam to lower case
Ci   ista:  -1  return 1 if file exists, 0 otherwise
Ci         >=0  sets flags file file status on opening.  This status
Ci              is kept and used along with the file table.  Internal
Ci              values are used except when missing or preempted:
Ci           ...The following bits work as independent switches:
Ci           ...bits 0,1:
Ci           0  open the file as 'UNKNOWN', if not open already
Ci           1  open the file as 'OLD'
Ci           2  open the file as 'NEW'
Ci           ...bit 2:
Ci           4  open the file as 'UNFORMATTED'
Ci           ...bit 3:
Ci           8  Override internal ista for bits 0,1
Ci          16  Override internal ista for bit 2
Ci   lunit:  (lunit>0) if in table, uses lunit to match instead of name
Ci              UNLESS mode, 10s digit set.
Cb Bugs
Cb   This routine mixes up many different functions.  
Cb   Operations handling internal tables of file names and logical units
Cb   should be separately packaged from those accessing (opening or
Cb   inquiring about) files.
Cr Remarks
Cr   mxnam is maximum number of files that may be open at one time.
Cr   All file closings should be made through entry fclose
Cr   17 June 97 internal ista bits 0,1 ALWAYS 0,
Cr              making default file opening 'unknown'
C ----------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer ista,lunit,mode
      character*(*) pnam
      character*11 ftnfmt
C Local variables
      logical bittst,isopen,ldum,mod100
      integer bit,bitand,getdig,i,i1,i2,iprint,locuni,lsta,matchn,
     .matchu,mode1,mode10,n,IDBG,isw
      character ftnsta*7, outs*80
C Parameters concerning limits to string size in filenames:
C plen: path length; extlen: ext length; maxsiz: name (including ext)
      integer maxsiz,extlen,plen
c #if unix
      parameter (extlen=20, plen=80, maxsiz=100)
c #elif MSDOS
c      parameter (extlen=4, plen=1, maxsiz=12)
c #endif
c #if DEBUG
c      parameter (IDBG=10)
c #else
      parameter (IDBG=110)
c #endif
C Save list of filenames
      integer mxnam
      parameter (mxnam=30)
      character*(maxsiz) tabnam(mxnam),pfnam,fnam,spfnam
      character*(extlen) ext
      character*(plen) path
      integer tabsw(mxnam), tabuni(mxnam),nnam,lgunit
      save tabnam, tabsw, tabuni, nnam, ext, path
      data nnam /0/
c #if FEXTENS
      data ext /'.dat'/
c #else
c      data ext /' '/
c #endif
      data path /' '/


      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)

      mode1   = getdig(mode,0,10)
      mode10  = getdig(mode,1,10)
      mod100 = getdig(mode,2,10) .ne. 0

C --- Cases handling partial file names ---
      fopnx = 0
C ... Change extension
      if (mode10 .eq. 2) then
        call namcat(ext,pnam,' ',mod100)
        if (iprint() .ge. IDBG) call awrit0(
     .  ' FOPNX: set ext to "'//ext//'%a"',' ',80,lgunit(1))
        return
C ... Copy extension into pnam
      elseif (mode10 .eq. 3) then
        call namcat(pnam,ext,' ',.false.)
        return
C ... Change internal path name
      elseif (mode10 .eq. 4) then
        call word(pnam,1,i1,i2)
        if (i2 .le. 0) then
          path = ' '
C        elseif (pnam(i2:i2) .eq. '/') then
C          call namcat(path,pnam,' ',mod100)
C        else
C          call namcat(path,pnam,'/',mod100)
        else
          call namcat(path,pnam,' ',mod100)
        endif
        if (iprint() .ge. IDBG) call awrit0(
     .  ' FOPNX: set path to "'//path//'%a"',' ',80,lgunit(1))
        return
C ... Copy path into pnam
      elseif (mode10 .eq. 5) then
        call namcat(pnam,path,' ',.false.)
        return
C ... Show internal table
      elseif (mode10 .eq. 9) then
        call awrit3('%N fopen:  %i files in table.'//
     .  '%?;n;  path='''//path//'%a'';;'//
     .  '%?;n;  ext='''//ext//'%a'';;',' ',-80,lgunit(1),nnam,
     .  isw(path.ne.' '),isw(ext.ne.' '))
        n = 0
        do  10  i = 1, nnam
          call word(tabnam(i),1,i1,i2)
          n = max(n,i2)
   10   continue
        n = n+2
        call awrit1(' file%npunit   mode  open',' ',80,lgunit(1),n)
        call pshpr(0)
        do  11  i = 1, nnam
          call awrit4(' '//tabnam(i)//
     .    '%np%,2i    %?;n==0;  fmt;unfmt;   %l',' ',80,lgunit(1),n,
     .    tabuni(i),bitand(tabsw(i),4),isopen(tabuni(i),.false.))
   11   continue
        call poppr
        print *, ' '
        return

      endif

C --- Make full file name ---
      if (bittst(mode1,2)) then
        call namcat(pfnam,pnam,' ',mod100)
      else
        call namcat(pfnam,pnam,' ',mod100)
        call namcat(fnam,pfnam,ext,.false.)
        call namcat(pfnam,path,fnam,.false.)
      endif

C ... Return full file name in pnam
      if (mode10 .eq. 6) then
        call namcat(pnam,pfnam,' ',.false.)
        return
      endif

C --- Match lunit with table.  matchu=-1 unless a match exists ---
      matchu = -1
      if (lunit .gt. 0) then
        do  110  i = 1, nnam
          matchu = i
  110   if (tabuni(i) .eq. lunit) goto 111
        matchu = -1
  111   continue
      endif

C --- If no match, get next available unit number => locuni ---
      locuni = lunit
      if (matchu .eq. -1 .and. lunit .lt. 0) then
        locuni = 7
   12   locuni = locuni+1
        if (isopen(locuni,.false.)) goto 12
        do  14  i = 1, nnam
   14   if (tabuni(i) .eq. locuni) goto 12
      endif

C --- Match pfnam with table, unless not needed ---
      matchn = matchu
      if (matchu .le. 0 .or. bittst(mode1,4)) then
        do  30  i = 1, nnam
          matchn = i
   30   if (tabnam(i) .eq. pfnam) goto 31
        matchn = -1
   31   continue
      endif
      if (lunit .lt. 0 .and. matchn .gt. 0) then
        matchu = matchn
        locuni = tabuni(matchn)
      elseif (lunit .lt. 0 .and. bittst(mode1,4)) then
C        call rxs('fopnx: no unit or name specified, file ',pnam)
        fopnx = -1
        return
      endif

      fopnx = locuni

C ... Force match between name and unit, or error
      call awrit1('%xFOPNX: file '//pfnam//'%a, unit %i, '//
     .'conflicts with prior use',outs,80,0,locuni)
      if (matchu .ne. matchn .and. bittst(mode1,4)) then
        if (matchu .gt. 0)
     .  call awrit0('%a%10b file '//tabnam(matchu),outs,80,0)
        call rx(outs)
      endif

C --- Add a file to internal table ---
      if (matchu .lt. 0 .and. mode10 .le. 6) then
        if (mod(mode1,2) .eq. 1) then
          if (lunit .gt. 0)
     .    call rxi('FOPNX: attempt to remove nonexistent unit',lunit)
          call rxs('FOPNX: attempt to remove nonexistent name ',pfnam)
        endif
        nnam = nnam+1
        if (nnam .gt. mxnam) call
     .  rxi('FOPNX: too many file names.  max=',mxnam)
        tabnam(nnam) = pfnam
        tabuni(nnam) = locuni
        tabsw(nnam)  = mod(ista,8) - mod(ista,4)
        matchu = nnam
C --- Remove a file entry from the internal table ---
      elseif (mod(mode1,2) .eq. 1 .and. mode10 .le. 7) then
        i = tabuni(matchu)
        if (isopen(i,.false.)) call fclose(i)
        do  130  i = matchu+1, nnam
          tabnam(i-1) = tabnam(i)
          tabuni(i-1) = tabuni(i)
          tabsw(i-1)  = tabsw(i)
  130   continue
        nnam = nnam-1
        return
      endif

#if DEBUG
      if (iprint() .ge. IDBG) then
        call pshpr(0)
        call awrit3(' FOPEN: file '//pfnam//'%a unit %i  status %i'//
     .  ' open=%l',' ',80,lgunit(1),locuni,ista,isopen(locuni,.false.))
        call poppr
      endif
#endif

C --- Cases ista=-2,-1 ---
      if (ista .eq. -2 .or. mode10 .ne. 0 .and. mode10 .ne. 7) return
      if (ista .eq. -1) then
#if PRECONNECTED_UNITS | IBM_VM
        INQUIRE(UNIT=lunit,EXIST=ldum)
#else
        INQUIRE(FILE=pfnam,EXIST=ldum)
#endif
        fopnx = 0
        if (ldum) fopnx = 1
        return
      else
      endif

C ... Format status
      lsta = ista
      if (matchu .gt. 0) then
C        lsta = tabsw(matchu)
C        if (bittst(ista,8))  lsta = bitand(lsta,1023-3) + bitand(ista,3)
C        if (bittst(ista,16)) lsta = bitand(lsta,1023-4) + bitand(ista,4)
C        tabsw(matchu) = lsta - mod(lsta,4)
        lsta = bitand(lsta,1023-4) + bitand(tabsw(matchu),4)
        if (bittst(ista,16)) then
          lsta = bitand(lsta,1023-4) + bitand(ista,4)
          tabsw(matchu) = bitand(tabsw(matchu),1023-4) + bitand(ista,4)
        endif
      endif
      ftnfmt = 'FORMATTED'
      if (bittst(lsta,4)) ftnfmt = 'UNFORMATTED'

C ... Handle fortran status (slow, but apparently stable)
      ftnsta = 'UNKNOWN'
      if (mod(lsta,4) .eq. 1) ftnsta = 'OLD'
      if (mod(lsta,4) .eq. 2) then
        ftnsta = 'NEW'
        call fclose(locuni)
#if PRECONNECTED_UNITS | IBM_VM
        open(locuni,FORM=ftnfmt,STATUS='UNKNOWN')
#else
        open(locuni,FILE=pfnam,FORM=ftnfmt,STATUS='UNKNOWN')
#endif
        close(unit=locuni,status='DELETE')
      endif

c      print *, 'dbg: fopn',pfnam,locuni,lsta,ftnfmt,ftnsta

      if (.not. isopen(locuni,.true.)) then
        call word(pfnam,1,i1,i2)
        call chrpos(pfnam,' ',plen+maxsiz,i)
C       if (iprint() .ge. IDBG) print *, 'FOPEN: opening ',ftnfmt,
C    .    ' file ''',pfnam(1:i2),''', status=',ftnsta,', unit=',locuni
        if (iprint() .ge. IDBG) print '(a,a,a,a,a,a,a,i3)',
     .  ' FOPEN: opening ',ftnfmt,' file ''',pfnam(1:i2),
     .  ''', status=',ftnsta,', unit=',locuni

#if PRECONNECTED_UNITS | IBM_VM
        open(locuni,FORM=ftnfmt,STATUS=ftnsta)
#elif SUN_ULTRA
        spfnam = pfnam(i1:i2)
        open(locuni,FILE=spfnam,FORM=ftnfmt,STATUS=ftnsta)
        rewind locuni
#else
        open(locuni,FILE=pfnam(i1:i2),FORM=ftnfmt,STATUS=ftnsta)
        rewind locuni
#endif
      endif

      end
      logical function isopen(unit,sw)
C- Returns whether unit aleady open; optionally adds unit to list if not
C ----------------------------------------------------------------
Ci Inputs
Ci   unit: if unit not among list of open files, add to list
Ci   sw:   if false, return isopen but do not add to list
Co Outputs
Co   isopen
Cr Remarks
Cr   maxfil is maximum number of files that may be open at one time
Cr   All file closings should be made through entry fclose
C ----------------------------------------------------------------
C     implicit none
      integer unit
      logical sw
      integer i,iprint
      integer maxfil
      parameter (maxfil=15)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      if (iprint() .gt. 110) print 20, unit, (unitab(i), i=0, nopen-1)
   20 format(/' ISOPEN: check logical unit',i3,
     .' among open units: ',15i3)

      isopen = .true.
      do  10  i = nopen-1, 0, -1
   10 if (unitab(i) .eq. unit) return

      isopen = .false.
      if (.not. sw) return
      unitab(nopen) = unit
      nopen = nopen+1
      if (nopen .gt. maxfil) call rx('ISOPEN: too many files')
      end
      subroutine fclose(unit)
C- Closes an open file, removing unit from the stack
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   none
Cr Remarks
Cr   use in conjunction with isopen
C ----------------------------------------------------------------
C     implicit none
      integer unit
      integer i,iprint
      character*6 clstat
      integer maxfil
      parameter (maxfil=15)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      clstat = 'KEEP'
      goto 10

   10 continue
      if (iprint() .ge. 100) print 20, unit
   20 format(' FCLOSE: closing',i3)

      do  30  i = nopen-1, 0, -1
   30 if (unitab(i) .eq. unit) unitab(i) = 999
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

C --- Close all open files ---
      entry closea
      do  50  i = 0, nopen-1
        if (iprint() .ge. 100) print *,'CLOSEA: closing file ',unitab(i)
        close(unit=unitab(i),status='KEEP')
        unitab(i) = 999
   50 continue
      nopen = 0
      return
      end
      integer function lgunit(i)
C- Returns stdout for i=1, log for i=2, mlog for i=3 (MPI logfile)
C     implicit none
      integer i, fopn, i1mach, fhndl
      lgunit = i1mach(2)
      if (i .eq. 1) return
      if (i .eq. 2) then
        lgunit = -1
        if (fhndl('LOG') .lt. 0) return
        lgunit = fopn('LOG')
      elseif (i .eq. 3) then
        lgunit = -1
        if (fhndl('MLOG') .lt. 0) return
        lgunit = fopn('MLOG')
      endif
      end
      subroutine namcat(dest,fnam1,fnam2,lcase)
C- Concatenates two names, optionally lowering the case
C     implicit none
      logical lcase
      character*(*) fnam1,fnam2,dest
      integer iprint,i1,i2,j1,j2

      call word(fnam1,1,i1,i2)
      call word(fnam2,1,j1,j2)
      i2 = max(i1-1,i2)
      j2 = max(j1,j2)
      if (i2 .ge. i1) then
        dest = fnam1(i1:i2)//fnam2(j1:j2)
      else
        dest = fnam2(j1:j2)
      endif
      if (i2-i1+j2-j1+2 .gt. len(dest) .and. iprint() .gt. 0) then
        print *, 'namcat: string "',fnam1(i1:i2),fnam2(j1:j2),
     .  '" was truncated to "',dest,'"'
      endif

#if unix
      if (lcase) call locase(dest)
#endif

      end
      subroutine ftflsh(ifi)
C- Flushes output buffer for file ifi
C  NB: this routine is nonstandard, and machine-dependent.
C     implicit none
      integer ifi
      character fname*101,fmt*33
      integer i1mach,j1,j2
      logical open

      if (ifi .eq. i1mach(2)) then
        call flushs(1)
      elseif (ifi .lt. 0) then
        call flushs(-1)
#if AIX
      else
        inquire (UNIT=ifi,OPENED=open,NAME=fname,FORM=fmt)
        call word(fname,1,j1,j2)
C        write (6,220) ifi,open,fname(j1:j2),fmt
C  220   format(' unit',i3,' open=',l1,'  name=',a,'  fmt=',a)
        if (open) then
          close(ifi)
          open(ifi, file=fname, form=fmt, status='old')
        endif
#endif
      endif
      end
#if TEST
      subroutine fmain
C      implicit none
      integer fopna,fext,fextg,fopn,fopno,fopnn,fhndl,fadd,fxst,fpath,
     .fopng,fopnx,fopnT
      integer ifi,jfi,kfi,lgunit,i
      double precision s(2,4)
      character fnam*20
*      call pshpr(110)

C     Check fopnT
C      i = fopnT(' ' ,-1,0,12)
C      i = fopnT('tmp' ,-1,0,0)
C      ifi = fopn('tmp')      
C      call fshow
C      i = fhndl('tmp')
C      print *, i
C      call dfclos(fopn('TMP'))
C      i = fopnT('tmp' ,-1,0,2)
C      print *, i
C      call fshow
C      i = fhndl('tmp')
C      print *, i
C      call rx('done')

      s(1,1) = 11
      s(1,2) = 12
      s(2,1) = 21
      s(2,2) = 22

      fnam = 'xyz.123'
      ifi = -1
      print *, 'test fopng(fnam=xyz.123,-1,4)'
      ifi = fopng(fnam,-1,4)
      call ywrm(1,'testing fopng',1,ifi,' ',s,2,2,2,2)
      call fshow
      call awrit0('%N test fclose:',' ',80,lgunit(1))
      call fclose(ifi)
*     print *, 'test close(fnam=xyz.123,-1,4)'
      call fshow

      call awrit0('%N test fopna(fnam=pqr,ifi,0)',' ',80,lgunit(1))
      ifi = 10
      ifi = fopna('pqr' ,ifi,0)
      call fshow
      call dscal(4,2d0,s,1)
      call ywrm(0,' testing fopng',1,ifi,'(5f12.2)',s,2,2,2,2)

      call awrit0('%N test fhndl(xyz),fhndl(pqr)',' ',80,lgunit(1))
      ifi = fhndl('xyz')
      jfi = fhndl('pqr')
      print '(2i4)', ifi,jfi

      call awrit0('%N retest fopng(fnam=xyz.123,-1,4)',' ',80,lgunit(1))
      ifi = fopng(fnam,-1,4)
      call fshow

      call rdm(ifi,0002,8,' ',s,2,2)
      call ywrm(0,' testing fopng',1,lgunit(1),'(5f12.2)',s,2,2,2,2)

      call awrit0('%N retest fopng, file already open',' ',80,lgunit(1))
      ifi = fopng(fnam,-1,4)
      call fshow

      call awrit0('%N test fopng, file already open',' ',80,lgunit(1))
      ifi = fopng(fnam,-1,4)
      call fshow

      call awrit0('%N test fadd, file=abc',' ',80,lgunit(1))
      ifi = 14
      ifi = fadd('abc',ifi,0)
      call fshow

      call awrit0('%N call fext ext=dat2',' ',80,lgunit(1))
      ifi = fext('.dat2')
      call fshow

      call awrit0('%N test fadd, file=abc with diff ext',' ',80,
     .lgunit(1))
      ifi = fadd('abc',-1,0)
      call fshow

      call awrit0('%N test fopn (abc)',' ',80,lgunit(1))
      ifi = fopn('abc')
      call ywrm(0,' testing fopn',1,ifi,'(5f12.2)',s,2,2,2,2)
      ifi = fopn('abc')
      call ywrm(0,' testing fopn',1,ifi,'(5f12.2)',s,2,2,2,2)
      call fshow
      jfi = fopna('abc',-1,0)
      if (ifi .ne. jfi) stop 'bug in fopna'

      call awrit0('%N test ftflsh (stdout, abc.dat2)',' ',80,lgunit(1))
      call cwrite(' ... testing ftflsh(6) ... ',0,25,0)
      call ftflsh(6)
      write(*,*) 'after ftflsh(6)'
      write(*,*) '... calling ftflsh(-1)'
      write(ifi,*) '... before ftflsh(-1)'
      call ftflsh(-1)
      write(*,*) '... calling ftflsh(ifi)'
      write(ifi,*) '... before ftflsh(ifi)'
      call ftflsh(ifi)
      write(ifi,*) '... after ftflsh(ifi)'

      print *, '... invoking system  call ''cat abc.dat2'' ...'
      call ftflsh(6)
      call fsystm('cat abc.dat2',i)
      call awrit1('%N ... system  call ''cat abc.dat2'' returned... %i',
     .' ',80,lgunit(1),i)

      call awrit0('%N restore ext to dat',' ',80,lgunit(1))
      ifi = fext('.dat')
      call awrit0('%N test fclr (abc) by name',' ',80,lgunit(1))
      call fclr('abc',-1)
      call fshow

      call awrit0('%N test fclr (pqr) by number',' ',80,lgunit(1))
      call fclr(' ',10)
      call fshow

      call awrit0('%N test fopna (xyz) to show suppression of list',
     .' ',80,lgunit(1))
      ifi = fopna('xyz',-1,8)
      write (ifi,*) 'test'
      call fclose(ifi)
      call fshow
      ifi = fopna('xyz',-1,0)
      read(ifi,'(a10)') fnam
      print *, 'after opening, read from xyz: ',fnam

      call awrit0('%N ... test fxst (xyz) to show suppression of list',
     .' ',80,lgunit(1))
      ifi = fxst('xyz')
      jfi = fxst('zzz')
      kfi = fopnx('xyz.123',172,-1,-1)
      call awrit3(' ifi for xyz is %i, for zzz is %i, for xyz.123 %i',
     .' ',80,lgunit(1),ifi,jfi,kfi)
      call fshow

      call awrit0('%N ... show fopn works without fadd',
     .' ',80,lgunit(1))
      ifi = fopn('TMP')
      call fshow

      call awrit0('%N ... call fopnx to test locase, "NEW","OLD"',
     .' ',80,lgunit(1))
      ifi = fopnx('TMP',2,2,-1)
      call fclose(ifi)
      call fshow
      ifi = fopnx('TMP',2,8+1,-1)
      call fshow

      call awrit0('%N ... use path as pre-extension',' ',80,lgunit(1))
      ifi = fext(' ')
      ifi = fpath('dat.')

      ifi = fopn('tmp')
      call fshow
      call fclr('TMP',-1)
      call fshow

      call awrit0('%N ... check opening as NEW',' ',80,lgunit(1))
      ifi = fopnx('STR',100,16+8+4+2,-1)
      call fshow
      call fclose(ifi)

      call awrit0('%N ... check fopn, UNFORMATTED',' ',80,lgunit(1))
      ifi = fopn('STR')
      call fshow

      call closea
      call fshow

      end
#endif

