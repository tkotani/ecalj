      integer function iosite(lio,vn,errh,filnam,ifi,slabl,alat,plat,
     .nbas,nspec,
     .pos,vel,eula,vshft,ips,ipl,irlx)
C- File i/o of site data
C ----------------------------------------------------------------
Ci Inputs
Ci   lio  :controls file handling and header I/O
Ci         1s digit: file handling
Ci           0 for read
Ci           1 for write
Ci           2 add to use given ifi instead of obtaining it from fnam
Ci         10s-100s digit: bits set conditions for file read
Ci           1    file alat must match passed array
Ci           2    file plat must match passed array
Ci           4    file nbas must match passed value
Ci           8    Must read nbas data elements from file
Ci          16    version-id must match file
Ci        1000s to 10000s digit: bits indicate what info is read in
Ci           1    I/O alat
Ci                ... alat is optional, if sought; no value of alat
Ci                    is assigned if missing (standard format)
Ci           2    I/O plat required, if sought
Ci           4    I/O nbas required, if sought
Ci           8    I/O species, pos information
Ci                    If bit is not set, pos array is untouched
Ci                ... these data are optional, if sought, that is
Ci                    on file read, arrays are assigned to zero if
Ci                    missing from file data.
Ci                (version 3)
Ci          16    I   return local io as functional value
Ci          32    I/O vel,eula,PL,rlx,vshft (standard format only)
Ci                If data sought, but not available, parameters
Ci                are assigned to zero.
Ci                If data not sought, parameters are not assigned
Ci          64    I/O (empirical TB only)
Ci                (version 2)
Ci          16    I/O vel,eula,PL,rlx,vshft (standard format only)
Ci          32    I/O (empirical TB only)
Ci
Ci   errh :controls error handling when data is missing or a required
Ci         match (described in lio, above) is not satisfied.
Ci         0 flags iosite to abort
Ci        >0 iosite returns <0 and prints a warning if verbosity > errh
Ci
Ci   nspec :number of species
Ci  filnam :file name.  See description of argument ifi.
Ci
Ci   vn    : version id.
Ci   slabl :
Ci
Cio  alat  :Scaling of lattice and position vectors
Cio  plat  :Lattice vectors (dimensionless)
Ci
Cio  ifi   :File logical unit.
Cio       *If 1s digit of lio contains 2's bit,
Cio        iosite uses ifi as the file handle and fnam is not used
Cio       *Otherwise, the file 'filnam' is opened and assigned to ifi.
Ci
Cio nbas:   size of basis
Cio
Cio sspec
Cr Remarks
Cr  *File consists of 1st record specifying version number
Cr     % version-id ...
Cr   The file version-id consists of a floating-point number,
Cr   followed by an optional string which marks file types different
Cr   from the standard. See below for possible formats.
Cr   A check is made that caller version no matches file.
Cr
Cr  *Standard format
Cr   First line has nbas and plat data, in this format:
Cr     % version-id [options] io=# nbas=# plat=# # ..
Cr   where the tokens specify:
Cr     version-id is the current version number (2.0)
Cr     io=  indicates file's contents (1000-10000s digits of lio)
Cr     nbas= is the size of the basis
Cr     plat= are the lattice vectors
Cr     ... and there are the following optional tokens:
Cr     fast => data is not parsed through preprocessor and
Cr             may be read using fortran read (no algebra)
Cr     xpos => site positions are expressed as multiples of plat
Cr             as opposed to Cartesian coordinates
Cr
Cr   Following the header line is one single line for each site.
Cr   Its structure is (version 2.0)
Cr     spid   x y z   vx vy vz   euler-angles  PL  rlx(xyz)
Cr   Example site file:
Cr     % vn=3.0 fast io=14 nbas=1 plat=-.5 .5 .5 .5 -.5 .5 .5 .5 -.5
Cr      FE    0 0 1
Cr
Cr  *Kotani format
Cr   Any line beginning with '#' is ignored (comment line)
Cr   First line has
Cr     % [#:]kotani   where `#' is version # (2.0)
Cr   If no `#' is specified, no check is made on version number
Cr   Then follows these lines:
Cr       alat                     (lattice constant, in a.u.)
Cr       plat(1:3,1)              (first lattice vector)
Cr       plat(1:3,2)              (second lattice vector)
Cr       plat(1:3,3)              (third lattice vector)
Cr   Next follow site data, in this format, one for each site
Cr       ibas, iclass, species, pos(ibas)
Cr   Example: CaO
Cr      % vn=kotani
Cr      10.26                ! alat
Cr       0   .5   .5         ! plat(1:3,1)
Cr      .5   .0   .5         ! plat(1:3,2)
Cr      .5   .5   .0         ! plat(1:3,3)
Cr      1 1 Ca  0.0 0.0 0.0  ! ibas, iclass, species, pos
Cr      2 2 Ca  0.3 0.3 0.3  ! etc
Cr      3 3 O   0.5 0.5 0.5  !
Cr      4 3 O   1.0 1.0 1.0  !
Cr   NB: iosite ignores 'class' column.
Cb Bugs
Cb   No check is made on version number, Kotani style input.
Cm MPI
Cm   Note: this routine should only be called by the master node!
Cu Updates
Cu   30 Oct 07 New file, adapted from iosits.f
Cu   20 Apr 07 Bug fix, blank lines at EOF
Cu   14 Apr 03 MPI enabled
Cu   24 May 02 standard format will read alat (optional)
Cu   11 Jan 02 xpos option implemented
Cu   11 Jan 02 iosite can read input using a2vec
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lio,errh,nbas,nspec,ifi
      character*8 slabl(1)
      character*(*) filnam
      double precision sspec(1),alat,plat(3,3),vn
      double precision pos(3,*),vel(3,*),eula(3,*),vshft(*)
      integer ips(*),ipl(*),irlx(3,*)
C ... Local parameters
C ... These are for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),aa*(recl),
     .vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .nlist(0:mxlev)
C Local variables
      character locnam*8, spid*8, readin*6
      double precision alatx,posl(3),vell(3),eulal(3),plx(9),ddot,
     .vshftl,delta(6),vnx
      logical bittst,ltmp,parstr,a2bin,lfast,lxpos,rdstrn
      integer a2vec,ipr,lio0,lio12,lio345,llio345,fopnx,ib,is,irlxl(3),
     .ipll,lgunit,i,j,nbasl,ix(9),ndel,iirlx,j1,j2,getdig,bitand
      integer procid,master
      data procid /0/ master /0/
      data rdarg /'#{}% c'/

C ... Setup
      iosite = -1
      call getpr(ipr)
      locnam = filnam
      lio0   = mod(lio,10)
      lio12 = mod(lio/10,100)
      lio345 = mod(lio/1000,1000)
      if (procid .eq. master) then
        if (.not. bittst(lio0,2)) ifi = fopnx(filnam,100,16+8+0,-1)
        if (lio0 .lt. 4) rewind ifi
      endif

C --- File write ---
      if (mod(lio0,2) .eq. 1 .and. procid .eq. master) then

C   ... Write header
        if (alat .gt. 0) then
          call awrit4('%% site-data vn=2.0 fast io=%i nbas=%i'//
     .    ' alat=%;8,1d plat=%9:1;8,1d',a,len(a),ifi,lio345,nbas,alat,
     .    plat)
        else
          call awrit3(
     .    '%% site-data vn=2.0 fast io=%i nbas=%i plat=%9:1;8,1d',a,
     .    len(a),ifi,lio345,nbas,plat)
        endif
        write(ifi,344)
  344   format('#',t26,'pos',t61,'vel',t85,'eula',t103,'vshft',t111,
     .  'PL rlx')
        do  10  ib = 1, nbas
C          call upack('site spec posl vell eulal',ssite,ib,is,posl,vell,eulal)
C          call upack('site pl relax vshft',ssite,ib,ipll,irlxl,vshftl,0)
          stop 'fix'
          spid = slabl(is)
          a = ' '//spid
C          call awrit6('%9p%3:1;10F %3:1;8F %3:1;8F  %;8F%,3i %3,1i',
C     .      a,len(a),-ifi,posl,vell,eulal,vshftl,ipll,irlxl)
          call awrit6('%8p%3;12,7D%3:1;12,7D%3:1;12,7D %;8F%,3i %3,1i',
     .    a,len(a),-ifi,posl,vell,eulal,vshftl,ipll,irlxl)
C         Write tbe data
          if (bittst(lio345,64)) then
            call dpzero(delta,6)
C            call upack('site ndelta delta',ssite,ib,ndel,delta,0,0)
            call awrit1(' tbe:%9p%6:1;9F',' ',80,ifi,delta)
          endif
   10   continue

        if (ipr .ge. 30) call awrit1(' IOSITE: wrote to file '''//
     .  locnam//'%a'', %i sites',' ',80,lgunit(1),nbas)

C --- File read ---
      else
        if (procid .eq. master) then
          ib = 0
          readin = 'header'
          call getpr(ipr)
          nr = 0
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          if (nr .eq. 0) goto 99

C   ... Check for alternative formats
          j1 = 0
          if (.not. parstr(a,'vn=',len(a)-3,3,'=',j1,j2)) goto 98
          if (a(j1:j1) .ne. ' ') goto 40
          call nword(a,1,j1,j2)
          i = 0
          if (.not. parstr(a(j1+3:),'kotani',j2-j1-6,6,'=',i,j)) goto 40

C   --- File read, Kotani's format : no algebraic expressions ---
          llio345 = 15
          read(ifi,*) alatx
          if (bittst(lio345,1)) alat = alatx
          if (bittst(lio12,1)) then
            call fsanrg(alatx,alat,alat,1d-6,'iosite:','alat',.true.)
          endif

C   ... Read plat
          read(ifi,*,err=98,end=98) plx(1),plx(2),plx(3)
          read(ifi,*,err=98,end=98) plx(4),plx(5),plx(6)
          read(ifi,*,err=98,end=98) plx(7),plx(8),plx(9)
          if (bittst(lio345,2)) call dcopy(9,plx,1,plat,1)
          if (bittst(lio12,2)) then
            call daxpy(9,-1d0,plat,1,plx,1)
            if (ddot(9,plx,1,plx,1) .gt. 1d-14) goto 99
          endif

C   ... Read nbas.  Find indirectly by counting lines.
C       First char = `#' => comment line
          if (bittst(lio345,4)) then
C         i = 0
            nbasl = 0
   31       continue
            if (rdstrn(ifi,a,recl,.false.)) then
C           Skip comment lines
              if (a(1:1) .eq. '#') goto 31
C           Formatted read ibas as a float, for portability
              call word(a,1,j1,j2)
              read(a(j1:j2),'(E72.0)',err=98) alatx
C           Error unless ibas is an integer
              ib = alatx
              if (dble(ib)-alatx .ne. 0) goto 98
C           Keep track of largest value of ib
C           i = max(i,ib)
              nbasl = nbasl+1
C           For now, sites must be ordered
              if (ib .ne. nbasl) goto 98
              goto 31
            endif
            if (bittst(lio345,4)) nbas = nbasl
            if (bittst(lio12,4))
     .      call isanrg(nbasl,nbas,nbas,'iosite:','file''s nbas',.true.)
C         file pointer now at EOF.  Restore position
            rewind ifi
            nr = 0
            call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .      ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
            do  i = 1, 4
              ltmp = rdstrn(ifi,a,len(a),.false.)
            enddo
          endif

C  ...  Read site data, Kotani format
          if (lio345 .ge. 8) then
            do  ib = 1, nbas
              if (.not. rdstrn(ifi,a,len(a),.false.)) goto 98
C         First and second columns (site and class index) ignored
              call word(a,3,j1,j2)
C         Column 3 is species index
              spid = ' '
              read(a(j1:j2),'(a)',err=98,end=98) spid
              call tokmat(spid,slabl,nspec,8,' ',is,j,.false.)
              is = is+1
C         Error if species not found in list
              if (is .le. 0) goto 98
C         Columns 4-6 are site position
C         call pack5('site spec posl vell eulal',ssite,ib,is,posl,vell,eulal)
C         call pack5('site pl relax vshft',ssite,ib,ipll,irlxl,vshftl,0)
              read(a(j2+1:),*,err=98,end=98) posl
              call dpzero(vel(1,ib),3)
              call dpzero(eula(1,ib),3)
              call ivset(irlx(1,ib),1,3,1)
              ipl(ib) = 0
              vshft(ib) = 0
              call dcopy(3,posl,1,pos(1,ib),1)
              ips(ib) = is
            enddo
          endif
          goto 999

C   ... Read, std format.  Return here to resume parsing for arguments
   40     continue
C   ... Read version and check compatibility with calling program
          i = j1+2
          if (.not. a2bin(a,vnx,4,0,' ',i,len(a)-5)) goto 99
          call isanrg(int(vnx),2,3,'IOSITE:',
     .    'file''s version number',.true.)
C   ... Read lio345
          i = 0
          readin = ' "io" '
          if (.not. parstr(a,'io=',len(a)-3,3,' ',i,j)) goto 99
          i = j-1
          if (.not. a2bin(a,llio345,2,0,' ',i,len(a)-3)) goto 99
C       Remake convention for llio345 if file is v2 but ask v3
          if (vnx .lt. 3 .and. vn .ge. 3) then
            llio345 = mod(llio345,16) + 32*(llio345/16)
          endif
C   ... Read nbas
          if (bittst(lio12,4) .or. bittst(lio345,4)) then
            i = 0
            readin = '"nbas" '
            if (.not. parstr(a,'nbas=',len(a)-5,5,' ',i,j)) goto 99
            i = j-1
            if (.not. a2bin(a,nbasl,2,0,' ',i,len(a)-5)) goto 99
            if (bittst(lio345,4)) nbas = nbasl
            if (bittst(lio12,4))
     .      call isanrg(nbasl,nbas,nbas,'iosite:','file''s nbas',.true.)
          endif
C   ... Read alat (optional if sought), plat (reqd if sought)
          if (bittst(lio12,1) .or. bittst(lio345,1)) then
            i = 0
            readin = '"alat" '
            if (parstr(a,'alat=',len(a)-5,5,' ',i,j)) then
              i = j-1
              if (a2vec(a,len(a)-5,i,4,', ',2,-3,1,ix,plx) .eq. 1) then
                if (bittst(lio345,1)) alat = plx(1)
              endif
            endif
          endif
          if (bittst(lio12,2) .or. bittst(lio345,2)) then
            i = 0
            readin = '"plat" '
            if (.not. parstr(a,'plat=',len(a)-5,5,' ',i,j)) goto 99
            i = j-1
            if (a2vec(a,len(a)-5,i,4,', ',2,-3,9,ix,plx) .ne. 9) goto 99
            if (bittst(lio345,2)) call dcopy(9,plx,1,plat,1)
            if (bittst(lio12,2)) then
              call daxpy(9,-1d0,plat,1,plx,1)
              if (ddot(9,plx,1,plx,1) .gt. 1d-14) then
                call info(0,0,0,
     .          ' IOSITE: file plat does not match passed data',0,0)
                call rxs('IOSITE: file mismatch, file ',filnam)
              endif
            endif
          endif
C       See whether input file can be read with fortran read
          i = 0
          lfast = parstr(a,'fast ',len(a)-5,5,' ',i,j)
          i = 0
          lxpos = parstr(a,'xpos ',len(a)-5,5,' ',i,j)

          readin = 'data'
          if (lio345-bitand(lio345,16) .ge. 8) then
C         Exit if file doesn't contain enough data
            if (mod(lio345,16) .gt. mod(llio345,16)) goto 99
   41       a(1:1) = ' '
            aa(1:1) = ' '
            if (lfast) then
   43         ltmp = rdstrn(ifi,a,len(a),.false.)
              nr = nr+1
              if (ltmp .and. a(1:1) .eq. '#' .or.
     .        a .eq. ' ' .and. ib .lt. nbas) goto 43
              if (lio345 .ge. 64) then
   44           ltmp = rdstrn(ifi,aa,len(aa),.false.)
                nr = nr+1
                if (ltmp .and. aa(1:1) .eq. '#' .or. aa .eq. ' ') goto 44
              endif
            else
              call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .        ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
              ltmp =  nr .gt. 0
              if (ltmp .and. lio345 .ge. 64)
     .        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .        ilist,nlist,vnam,ctbl,mxchr,aa,recrd,recl,nr)
              ltmp =  nr .gt. 0
            endif

C         If all sites have been read, exit
            if (ib .ge. nbas) then
              call info2(40,0,0,' IOSITE: read pos'//
     .        '%?!n>15!,vel,eula,vshft,ipl,irlx!! from file unit %i',
     .        llio345,ifi)
              goto 999
            endif

            if (.not. ltmp) then
              if (ib .lt. nbas .and. bittst(lio12,8)) then
                if (errh.le.ipr .or. errh.eq.0) print *,
     .          'IOSITE: file does not contain',nbas,' sites'
                goto 99
              endif
              goto 999
            endif
C         Read data for this site
            ltmp = .false.
            call word(a,1,j1,j2)
            spid = ' '
            read(a(j1:j2),'(a)',err=98,end=98) spid
            call tokmat(spid,slabl,nspec,8,' ',is,j,.false.)
            is = is+1
            if (is .le. 0) goto 98
C         No vel,eula,vshft,ipl,irlx available on disk
            if (llio345 .lt. 32) then
              call dpzero(posl,3)
              if (lio345 .ge. 32) then  ! Default values if sought but not available
                call dpzero(vell,3)
                call dpzero(eulal,3)
                vshftl = 0
                ipll = 0
                iirlx = 111
              endif
              if (lfast) then
                read(a(j2+1:),*,err=98,end=98) posl
              else
                i = 0
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,posl)
     .          .ne. 3) goto 98
              endif
C         vel,eula,vshft,ipl,irlx are available on disk
            else
              if (lfast) then
                read(a(j2+1:),*,err=98,end=98)
     .          posl,vell,eulal,vshftl,ipll,iirlx
              else
                i = 0
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,posl)
     .          .ne. 3) goto 98
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,vell)
     .          .ne. 3) goto 98
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,eulal)
     .          .ne. 3) goto 98
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,1,ix,vshftl)
     .          .ne. 1) goto 98
                if (a2vec(a(j2+1:),len(a)-j2,i,2,', ',2,3,1,ix,ipll)
     .          .ne. 1) goto 98
                if (a2vec(a(j2+1:),len(a)-j2,i,2,', ',2,3,1,ix,iirlx)
     .          .ne. 1) goto 98
              endif
            endif
            if (lxpos) then
              call dcopy(3,posl,1,plx,1)
              call dgemm('N','N',3,1,3,1d0,plat,3,plx,3,0d0,posl,3)
            endif

            ib = ib+1
            call dcopy(3,posl,1,pos(1,ib),1)
            ips(ib) = is

C          Copy vel,eula,vshft,ipl,irlx if sought
            if (lio345 .ge. 32) then
              do  46  i = 1, 3
                irlxl(i) = getdig(iirlx,3-i,10)
   46         continue
              call dcopy(3,vell,1,vel(1,ib),1)
              call dcopy(3,eulal,1,eula(1,ib),1)
              vshft(ib) = vshftl
              ipl(ib) = ipll
              call icopy(3,irlxl,1,irlx(1,ib),1)
            endif

            if (lio345 .ge. 64) then
              call word(a,1,j1,j2)
              locnam = ' '
              read(a(j1:j2),'(a)',err=98,end=98) locnam
              read(aa(j2+1:),*,err=98,end=98) delta
            endif
C         We can quit now if ib eq nbas but read one more line in case
C         site file has comments at the end of file
C         if (ib .lt. nbas) goto 41
            goto 41
          endif

C   ... Exit point for read
  999     continue
        endif
      endif

      iosite = 0
C     Copy llio345 if sought
      if (bittst(lio345,16)) iosite = llio345
      return

C ... Error handling
   98 if (errh.le.ipr .or. errh.eq.0) then
        call awrit1(' IOSITE line %i: missing or incompatible data:',
     .  ' ',80,lgunit(1),nr)
        a(70:) = ' ...'
        call awrit0('%a',a,len(a),-lgunit(1))
        call rx('failed to read site data')
      endif

   99 continue
C     Copy llio345
      if (bittst(lio345,16)) lio = llio345

      if (errh .eq. 0) then
        print 369, ':', readin, filnam
      elseif (errh .le. ipr) then
        print 369, '(warning)', readin, filnam
      endif
  369 format(1x,'IOSITE ',a,' failed to read ',a,' from file ',a)
      if (errh .eq. 0) call rxs('IOSITE: file mismatch, file ',filnam)

      end

