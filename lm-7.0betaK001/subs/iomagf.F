      subroutine iomagf(nbas,nl,bfield,nbf,ifi)
C- I/O of external magnetic field
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ifi   :file logical unit, <0 for write, >0 for read
Co Outputs
Co   B-field vectors are read in or written out
Co   nbf  is read in or written out (either 1, nl or nl**2)
Cr Remarks
Cr   Aborts on read when nbas does not match file
Cr   19 Nov 97 use rdfiln for reading first line
Cu Updates
Cu   27 Mar 04 
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nbf,nl,ifi
      double precision bfield(nbas,nbf,3)
C ... Local parameters
      integer i,j,k,ipr,nc,rdm
      logical parstr,a2bin
      character*80 ss,ss2*4
      double precision tmp(125)
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
      equivalence (a,ss)
      data rdarg /'#{}% c'/

C --- Write ---
      if (ifi .lt. 0) then
C       rewind (-ifi)
        if (nbf .ne. 1 .and. nbf .ne. nl .and. nbf .ne. nl*nl)
     .    call rxi('iomagf: bad dim for bfield: nbfield=',nbf)
        call awrit3('%% rows %i cols %i nbf %i',ss,80,0,nbas*nbf,3,
     .    nbf)
        call awrit0('%a',ss,80,ifi)
        do  i = 1, nbas
          if (nbf .gt. 1) call awrit1('# ib %,4i',' ',80,-ifi,i)
          write(-ifi,'(3f16.12)') ((bfield(i,j,k),k=1,3),j=1,nbf)
        enddo
C --- Read ---
      else
        call getpr(ipr)
        rewind ifi
        nr = 0
        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,ss,recrd,recl,nr)
C   ... No error if file is empty, but don't bother reading
        if (nr .eq. 0) then
          if (ipr .gt. 30) print '('' IOMAGF:   empty file'')'
          return
        endif
        if (ss(1:1) .eq. '%') then
C   ... Read nbf if it is there
        nbf = 1
        i = 0
        if (parstr(ss,'nbf ',len(ss)-4,4,' ',i,j)) then
          if (.not. a2bin(ss,nbf,2,0,' ',j,len(ss)))
     .      call rx('iomagf: failed to parse for nbf')
        endif
        else
          nbf = 1
          rewind ifi
        endif

C   ... Expect nbas*nbf*3 rows in file ... skip reading if not
        if (nbf .ne. 1 .and. nbf .ne. nl .and. nbf .ne. nl*nl)
     .    call rxi('iomagf: illegal dimension: nbf =',nbf)
        nc = 3
        ss2 = 'in'
        if (nbf .eq. nl) ss2 = ' '
        if (nbf .eq. nl*nl) ss2 = '%bm-'
        call info2(31,0,0,' IOMAGF:  reading l-'//ss2//
     .    '%adependent magnetic field (file %i)',ifi,0)

C   ... Read the data, 1 atom at a time ...
        call dpzero(bfield,nbas*nbf*3)
        do  20  i = 1, nbas
          call rxx(rdm(ifi,0,nbf*nc,' ',tmp,nbf,nc).ne.1,
     .      'iomagf:  bad file')
          call dmcpy(tmp,nbf,1,bfield(i,1,1),nbas*nbf,nbas,nbf,nc)
   20   continue
        if (ipr .ge. 50) call pvioeu(11,0d0,bfield,nbas,nbf)
      endif
      end

