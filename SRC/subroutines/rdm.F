      integer function rdm(ifi,sw,mxelt,filel,s,nr,nc)
C- Matrix input, parsing file with rdfiln.
C ----------------------------------------------------------------
Ci  nr,nc: number of rows and columns (also output; see Remarks)
Ci  sw:    one's digit: 1 'quick read'--- read array with fortran read
Ci                      2 'binary read'--- read a binary file
Ci                      3 'binary read'--- nr,nc,ic not read
Ci                        (nr and nc must be passed; array assumed real)
Ci                      4 'binary read'--- same as 3 except that
Ci                        array is assumed to be c16 format
Ci                        (NOT IMPLEMENTED)
Ci         ten's digit: 0, accept real matrix s(nr,nc) only
Ci                      1, accept real s(nr,nc) or cmplx s(nr,nc,2)
Ci                      2, accept real s(nr,nc) or cmplx s(2,nr,nc)
Ci                         NB: implemented only for binary read
Ci                      3, accept cmplx s(nr,nc,2) only (See remarks)
Ci                      4, accept cmplx s(2,nr,nc) only
Ci        1000's digit: 0, ignore label if any exists
Ci                      1, return label in string filel
Ci       10000's digit: 0, matrix elements in normal order, ie
Ci                         s(1,1), s(1,2), ..., s(1,nc),
Ci                         s(2,1), s(2,2), ..., s(2,nc)...
Ci                      1, matrix elements in transpose order ie
Ci                         s(1,1), s(2,1), ..., s(nr,1),
Ci                         s(1,2), s(2,2), ..., s(nr,2)...
Ci                       NB: for this case, nr and nc must be spec'd
Ci                       explicitly, ie by priority 1 or 2 in Remarks
Ci                       NB: can also be set in the first line of the
Ci                       input; see Remarks
Ci  mxelt: maximum dimension of s (considered as a 1-dimensional array)
Ci         mxelt=0 => determine nr,nc but make no attempt to read s
Co  filel :text label describing file contents (see 1000s digit sw)
Co  s:    Data is read into this array
Co  rdm:  1 for real, 2 for complex, 3 for complex in complex*16 format,
Co        provided it evaluates all expressions with no errors and reads
Co        enough expressions to be consistent with nr and nc specified.
Co        If not, returns rdm = -1.
Cr Remarks
Cr  rdm attempts to read in matrix s by parsing file ifi, using routine
Cr  rdfiln for conditional reading of lines, looping syntax and
Cr  expression substitutions.  The number of rows nr and columns nc of s
Cr  need not be prescribed in advance.  rdm determines nr and nc
Cr  according to the following priority.
Cr  1. If on entry nc>0, nc always takes that value.  rdm will not set
Cr     nc.  If file also specifies nc (see below) not compatible
Cr     with passed nc, rdm fails and returns -1.  Similarly for nr.
Cr  ASCII read:
Cr  2. rdm reads the first line of the file.  If it has the syntax
Cr     % rows expr1 cols expr2
Cr     nr is assigned value expr1 and nc value expr2,
Cr     provided neither has been set in step 1
Cr     (NB: either rows declaration or cols declaration may be absent).
Cr  3. rdm parses the first line, and sets nc to the number of
Cr     expressions it finds on that line.  At this point, rdm has an
Cr     unambiguous value for nc (or else it has aborted)
Cr  4. If nr has not yet been assigned, the file is parsed until EOF is
Cr     reached or it encounters a line beginning with '%'.
Cr     nr is assigned the (total no expressions evaluated) / nc.
Cr     When nr is determined in this way, mxelt must be double the
Cr     number of elements in s, since rdm must read in s-transpose,
Cr     and copy s-transpose back to s.
Cr     NB: nr may not be determined this way using 'quick-read'
Cr  BINARY read:
Cr  2. rdm reads the first record in the file.  It must contain
Cr     nr nc ic
Cr     It MAY contain additionally
Cr     nr nc ic iswitch
Cr       where:
Cr     If nr is >0 on input, file nr must match passed nr; ditto for nc.
Cr     ic is the cast (0=int, 1=dble, 2=cmplx, split real, imag, 3=c*16)
Cr     iswitch indicates next recrd is label
Cr   NB: this step is missing if sw = 1's digit of ic is 3
Cr  Complex matrices: a file is unambiguously
Cr     specified to be real or complex if the first line begins
Cr     with % and contain 'real' or 'complex'  (eg, '% rows 7 complex').
Cr     Otherwise, it is assumed to real unless a complex matrix was
Cr     requested (10's digit of sw=3,4).  rdm will convert if the cast
Cr     sought conflicts with the cast specified.  The resultant cast
Cr     depends on the ten's digit of sw, as follows:
Cr             |   sw=0         sw=3,4       sw=1,2
Cr             |  (want R)     (want C)     (unspec)
Cr      -------|---------------------------------------
Cr     Have R  |      R          R->C          R
Cr     Have C  |   C->R             C          C
Cr     unspec  |      R             C          R
Cr  Tranpose of matrices: rdm will read the matrix elements in transpose
Cr    order if the appropriate switch is set (see sw, above) or if the
Cr    first line contains
Cr    % .... trans   ...
Cr
Cb  Bugs: no check to see whether matrix asserted to be symmetric
Cb        (hermitian) actually is.
Cu  Updates
Cu  09 Nov 05  rdm can read complex arrays into c*16 format
Cu  16 jan 01  rdm can read complex arrays, sparse format
Cu  17 Mar 99  rdm can read integer array (converting it to real)
Cu  17 Nov 98  rdm can read real sparse matrices
C ----------------------------------------------------------------
C     implicit none
      integer ifi,mxelt,nr,nc,sw
      character*(*) filel
      double precision s(1)
C Local variables
      integer recl,nl,wantC,haveC,sw10,ncx,ioff,iprint
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl), a*(recl)
      character ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev),a2bin,qread,ltrnsi,ltrns,pvrdm1,
     .lddump,lidump,lpr,lsprse
#if CCAR
      logical cmdopt
#endif
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .nlist(0:mxlev),i1mach,i,ip,k,nelt,ic,ir,retval,ii,ilbl,nrx,ipr,
     .ndir,a2vec,ix(4),sw1
      double precision xwk(4),xx
      parameter (ndir=8)
      character vnam(mxlev)*16,dir(ndir)*7,llbl*1000,rdarg*7
      data rdarg /'#{}% ct'/
      data dir /'rows','cols','real','symm','complex','herm',
     .'trans','sparse'/

C --- Initialization ---
      ipr=iprint()
      lpr = ipr .gt. 0
      ltrnsi = .false.
      ltrns = mod(sw/10000,10) .eq. 1
      lsprse = .false.
      sw1 = mod(sw,10)
      qread = sw1 .eq. 1
      sw10 = mod(sw/10,10)
      call dpzero(xwk,4)
C     wantC and haveC:  0, unspecified, 1 real, 2 complex
      wantC = 1
      if (sw10 .eq. 1 .or. sw10 .eq. 2) wantC = 0
      if (sw10 .eq. 3 .or. sw10 .eq. 4) wantC = 2
      haveC = 0
      nelt = 0
      ir = 0
      ic = 0
      ioff = 0
      rdm = -max(wantC,1)
      nl = 0
      call dpzero(s,mxelt)
      retval = 0

C --- Binary read ---
      if (sw1 .eq. 2 .or. sw1 .eq. 3) then
        if (sw1 .eq. 2) then
          read(ifi,err=198,end=198) ir,ii,ic
          if (nr .eq. 0) nr = ir
          if (nc .eq. 0) nc = ii
          if (ir .ne. nr .or. ii .ne. nc) goto 196
C     ... look for a label
          backspace ifi
          read(ifi,err=197,end=197) ii,ii,ii,ilbl
          if (ilbl .eq. 0) goto 195
C     ... there is a label; read into llbl
          llbl = ' '
          read(ifi,err=198,end=198) llbl(1:min(len(llbl),ilbl))
          if (mod(sw/1000,10) .ne. 0) filel = llbl
          goto 195
C     ... recover from missing ilbl
  197     continue
          backspace ifi
          read(ifi,err=198,end=198) ii
C         read(ifi) s(1)
  195     continue
        else
          ic = 1
          if (nr .le. 0 .or. nc .le. 0) goto 199
        endif
        haveC = mod(ic,10)
        if (haveC .eq. 3) haveC = 2
C       if (ic .eq. 0 .and. wantC .gt. 0) then
        if (ic .eq. 0) then
          if (mxelt .lt. 2*nr*nc) goto 99
          if (.not. lidump(s(1+nr*nc),nr*nc,ifi)) goto 198
          call idscop(s(1+nr*nc),s,nr*nc,1,1)
          haveC = 1
          ic = 1
C       elseif (ic .eq. 0) then
C          if (mxelt .lt. nr*nc) goto 99
C         if (.not. lidump(s,nr*nc,ifi)) goto 198
C         rdm = ic
C         return
        else
          if (nr*nc*haveC .gt. mxelt) goto 99
          if (.not. lddump(s,nr*nc*haveC,ifi)) goto 198
        endif
        if (haveC .eq. 1 .and. wantC .eq. 2 .and. lpr) print '(a)',
     .  '#rdm (warning): sought complex matrix but file is real'
        if (haveC .ge. 2 .and. wantC .eq. 1 .and. lpr) print '(a)',
     .  '#rdm (warning): sought real matrix but file is complex'
        rdm = ic
        ltrnsi = ltrns
        goto 199
      endif

      if (sw10 .eq. 2)
     .call rx('rdm not ready for complex*16 s')

C --- Determine nr,nc from % rows and % cols ---
      recrd = ' '
#if CCAR
      if (cmdopt('-C',2,0,a)) rdarg(1:1) = a(3:3)
#endif
   21 call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
C ... Empty file kosher unless nr or nc nonzero
      if (nl .le. 0 .and. nr .eq. 0 .and. nc .eq. 0) rdm = max(wantC,1)
      if (nl .le. 0) return
C ... Skip blank record
      if (recrd .eq. ' ') goto 21
C ... First record: check for % rows and % cols
      if (recrd(1:1) .eq. '%') then
        ip = 1
   20   call skipbl(recrd,recl,ip)
        if (ip .ge. recl) goto 30
        k = ip-1
        call tokmat(recrd(ip+1:recl),dir,ndir,7,' ',i,ip,.false.)
        ip = ip+k
C       print *, 'after tokmat', i,ip
        if (i .lt. 0) then
          call skp2bl(recrd,recl,ip)
          goto 20
        endif
C   ... Read matrix elements in transpose order?
        if (i .eq. 6) then
          ltrns = .true.
          if (sw10 .eq. 4) call rx('rdm not ready for complex*16 s')
          goto 20
C   ... Read sparse matrix?
        elseif (i .eq. 7) then
          lsprse = .true.
          s(1) = 0
          if (sw10 .eq. 4) call rx('rdm not ready for complex*16 s')
          goto 20
        endif
C   ... Matrix is specified real or complex
        if (i .ge. 2 .and. i .le. 5) then
          haveC = 1
          if (i .gt. 3) haveC = 2
          if (i .eq. 3 .or. i .eq. 5) retval = 10
          goto 20
        endif
        if (ip .ge. recl) goto 30
        if (.not. (a2bin(recrd,xx,4,0,' ',ip,recl))) return
        if (i .eq. 0 .and. nr .ne. 0 .and. nr .ne. nint(xx)) goto 196
        if (i .eq. 1 .and. nc .ne. 0 .and. nc .ne. nint(xx)) goto 196
        if (i .eq. 0 .and. nr .eq. 0) nr = nint(xx)
        if (i .eq. 1 .and. nc .eq. 0) nc = nint(xx)
        goto 20
      else
C   ... First record not '%': skip reading this first record
        if (nr .ne. 0 .and. nc .ne. 0 .and. mxelt .eq. 0) goto 98
        goto  31
      endif

C --- Read next record from file ---
   30 continue
      if (nr .ne. 0 .and. nc .ne. 0 .and. mxelt .eq. 0) goto 98
      recrd = ' '
      if (.not. (qread .and. lsprse))
     .call rdfiln(ifi,'#{}%',mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
   31 continue

C --- Before reading, decide whether file is real or complex ---
      if (nelt .eq. 0) then
        if (haveC .ne. 0) then
          if (wantC .eq. 0) wantC = haveC
        elseif (wantC .eq. 1 .or. wantC .eq. 2) then
          haveC = wantC
        else
          wantC = 1
          haveC = 1
        endif
        if (haveC .eq. 1 .and. wantC .eq. 2 .and. lpr) print '(a)',
     .  '#rdm (warning): sought complex matrix but file is real'
        if (haveC .eq. 2 .and. wantC .eq. 1 .and. lpr) print '(a)',
     .  '#rdm (warning): sought real matrix but file is complex'
        if (wantC .eq. 0 .or. haveC .eq. 0) stop 'bug in rdm'
      endif

C --- nr and nc must be spec'd if lsprse ---
      if ((nr .eq. 0 .or. nc .eq. 0) .and. lsprse)
     .call rx('rdm not set up to read sparse matrix'//
     .' without explicit nr,nc')

C --- if nc zero, set nc to number of expressions in 1st recrd ---
C ... Case ltrns: nc and nr exchange roles
      k = 0
      if ((nc .eq. 0 .and. .not. ltrns .or.
     .nr .eq. 0 .and. ltrns .or. qread) .and.
     ..not. lsprse) then
        if (ltrns) then
          ip = nr
          nr = nc
          nc = ip
        endif
        ip = 0
   41   call skipbl(recrd,recl,ip)
        if (ip .ge. recl) then
          if (nc .eq. 0) nc = k
          goto 42
        endif
C   ... Read into scratch unless qread (then load as if s a vector)
        i = 0
        if (qread .and. mxelt .gt. 0) i=k
        if (i .gt. mxelt) goto 99
        if (.not. a2bin(recrd,s,4,i,' ',ip,recl)) return
        k = k+1
        goto 41
   42   continue
C ... Restore nr, nc to proper roles
        if (ltrns) then
          ip = nr
          nr = nc
          nc = ip
        endif
      endif

C ... nrx,  ncx are nr,nc or if ltrns,  nc,nr
      nrx = nr
      ncx = nc
      if (ltrns) then
        nrx = nc
        ncx = nr
      endif

C ... Now ncx's value is known.  ncx=0 => empty file.
      if (ncx .eq. 0 .or. (nrx .ne. 0 .and. mxelt .eq. 0)) then
        rdm = max(wantC,1) + retval
        return
      endif

C --- Quick-read: read remaining elements with unformatted read ---
   55 continue
      if (qread) then
C   ... Quick-read option requires that nr,nc are known by now.
        if (nr*nc .eq. 0) then
          print
     .    '(a)', '#rdm: quick-read set but nr or nc unspecified ...'
          return
        endif
        if (lsprse) then
          ltrns = .false.
          ltrnsi = .false.
          i = 3
          if (haveC .eq. 2) i = 4
          if (.not. pvrdm1(ifi,xwk,i)) goto 70
          i = xwk(1)
          ii = xwk(2)
          if (abs(i-xwk(1))+abs(ii-xwk(2)) .ne. 0 .or.
     .    i .gt. nr .or. ii .gt. nc) return
          s(i+nr*(ii-1)) = xwk(3)
          if (wantC .eq. 2) s(i+nr*(ii-1)+nr*nc) = xwk(4)
          goto 55
        else
          if (.not. pvrdm1(ifi,s(k+1),min(haveC,wantC)*nr*nc-k)) return
          ltrnsi = .not. ltrns
          goto 70
        endif
      endif

C --- Parse this record for array elements; load into s(ir,ic) ---
C     NB: if nr unknown, read in s-transpose, since nc is known.
C ... Exit if we are done reading this file
      if (recrd(1:1) .eq. '%') goto 70
      if (nl .le. 0) goto 70
      ip = 0
   61 call skipbl(recrd,recl,ip)
C   ... Pick up next record when this one parsed to the end
      if (ip .ge. recl) goto 30
C   ... Read into s(ir,ic) (or scratch if mxelt is zero)
      i = max(nr,1)*ic + ir + ioff
C       s(2,i,j) --- complex --- format
      if (sw10 .eq. 4) then
        i = 2*max(nr,1)*ic + 2*ir + ioff          
      endif
      if (mxelt .eq. 0) i = 0
      if (i .gt. mxelt) goto 99
      if (lsprse) then
        i = 3
        if (haveC .eq. 2) i = 4
        if (a2vec(recrd,recl,ip,4,' ',1,2,i,ix,xwk) .ne. i) return
        i = xwk(1)
        ii = xwk(2)
        if (abs(i-xwk(1))+abs(ii-xwk(2)) .ne. 0 .or.
     .  i .gt. nr .or. ii .gt. nc) return
        s(i+nr*(ii-1)) = xwk(3)
        if (wantC .eq. 2) s(i+nr*(ii-1)+nr*nc) = xwk(4)
        goto 30
      else
        if (.not. a2bin(recrd,s,4,i,' ',ip,recl)) return
      endif
C        print '('' nelt,i'',2i4,''  ir,ic,ioff='',3i4,f12.6)',
C     .    nelt,i,ir,ic,ioff,s(i+1)
      nelt = nelt+1
      ic = ic+1
      if (ltrns) then
        ic = ic-1
        ir = ir+1
      endif
      if (nc .gt. 0 .and. ltrns) then
        if (ir .gt. nr) call rx('bug in rdm')
        if (ir .eq. nr) then
          ic = ic+1
C       ... patch exit for now ... need to read to EOF too!
          if (nc .gt. 0 .and. ic .eq. nc) then
            if (min(wantC,haveC) .eq. 1 .or. ioff .ne. 0) goto 70
            ioff = nc*nr
            ic = 0
C             s(2,i,j) --- complex --- format
            if (sw10 .eq. 4) then
              ioff = 1
            endif
          endif
          ir = 0
        endif
      elseif (nr .gt. 0 .and. .not. ltrns) then
        if (ic .gt. nc) call rx('bug in rdm')
        if (ic .eq. nc) then
          ir = ir+1
C       ... patch exit for now ... need to read to EOF too!
          if (nr .gt. 0 .and. ir .eq. nr) then
            if (min(wantC,haveC) .eq. 1 .or. ioff .ne. 0) goto 70
            ioff = nr*nc
            ir = 0
C             s(2,i,j) --- complex --- format
            if (sw10 .eq. 4) then
              ioff = 1
            endif
          endif
          ic = 0
        endif
      endif
      goto 61

C --- Cleanup.  If nr is zero, set nr and copy s-transpose to s ---
   70 continue
      if (nc .eq. 0 .and. ltrns) then
        nc = nelt/nr/haveC
      elseif (nr .eq. 0 .and. .not. ltrns) then
        nr = nelt/nc/haveC
        ltrnsi = mxelt .ne. 0
      elseif (nelt .lt. nr*nc*haveC .and. .not. (qread .or. lsprse)
     ..and. lpr) then
        call awrit4('# rdm (warning) expected %i elements (%i rows,'//
     .  ' %i cols) but read %i',' ',80,i1mach(2),nr*nc*haveC,nr,nc,
     .  nelt)
      endif
      ncx = nc
      if (ltrnsi .and. wantC .eq. 2) ncx = 2*nc
      if (ltrnsi) then
        if (nelt*2 .gt. mxelt .and. mxelt .gt. 0) goto 99
        ioff = nr*ncx
        call dcopy(nr*nc,s,1,s(1+ioff),1)
        do    ir = 1, nr
          do  ic = 1, nc
             s(ir+nr*(ic-1)) = s(ioff+ic+nc*(ir-1))
          enddo
        enddo
        if (ncx .eq. 2*nc) then
          k = nr*nc
          call dcopy(k,s(1+k),1,s(1+ioff),1)
          do    ir = 1, nr
            do  ic = 1, nc
               s(ir+nr*(ic-1)+k) = s(ioff+ic+nc*(ir-1))
            enddo
          enddo
        endif
      endif

C --- Force matrix symmetric or hermitian if so stated ---
      if (retval .eq. 10) call dosymm(max(wantC,1),s,nr,nc)

C --- Normal exit ---
   98 rdm = max(wantC,1) + retval
      return

C --- Exit when not enough space to read s ---
   99 continue
      print '(a)', '#rdm: insufficent space ...'
      return

C --- Exit for incompatible input ---
  196 continue
      print '(a)', '# rdm: input dimensions incompatible with file'
      rdm = -1
      return

C --- Exit for binary read ---
  198 continue
      if (lpr) print '(a)', '#rdm: binary read failed...'
      rdm = -1
  199 continue
      end
      logical function pvrdm1(ifi,s,ns)
C- Kernel to read array s with unformatted read
C     implicit none
      integer ifi,ns
      double precision s(ns)

      pvrdm1 = .false.
      if (ns .gt. 0) read(ifi,*,err=99,end=99) s
      pvrdm1 = .true.
   99 continue
      end
      subroutine dosymm(icast,s,nr,nc)
C- Kernel to symmetrize array s
      integer nr,nc,nn
      double precision s(nr,nc,2),xx
      integer i,j,icast

      nn = min(nr,nc)
      do   i = 1, nn
        do  j = 1, i
          xx = (s(i,j,1) + s(j,i,1))/2
          s(i,j,1) = xx
          s(j,i,1) = xx
      enddo
      enddo

      if (icast .eq. 2) then
        do    i = 1, nn
          do  j = 1, i
            xx = (s(i,j,2) - s(j,i,2))/2
            s(i,j,2) = xx
            s(j,i,2) = -xx
         enddo
         enddo
      endif

      end

