      logical function iosg(lio,name,ivl,nl,nbas,ckbas,nttabg,
     .  oiax,onpr,onprg,osg)
C- Real space value-laplacian structure constants file read or write
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :specifies conditions for file i/o (condensed into digits)
Ci         :digit  specification
Ci         :1-10s   bits in 1-10s digit control program flow
Ci                  0 for file read
Ci                  1 for file write
Ci                  2 (read only) No header matching requirement;
Ci                    header parms returned into passed arguments
Ci                  4 read/write file header only and exit.
Ci                  8 Add this number to allocate arrays before reading
Ci                    or to release arrays after writing
Ci                 16 file is formatted (not implemented)
Ci   name  :file name holding strux (eg STRG)
Ci   ivl   :specifies functions used to built the value-laplacian strux
Ci   nl    :1 + global-maximum-l 
Ci          specifies the leading dimensions of sg 
Ci   nbas  :number of atoms in the basis
Ci   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Cio Inputs/Outputs:
Cio  nttabg:number of strux for the value-laplacian basis 
Cio  oiax  :offset for iax array
Cio  onpr  :offset for the total number of sites in each cluster ib, ntab(ib+1)
Cio  onprg :offset for the number of sites in each value-laplacian cluster, ntabg(ib)
Cio  osg   :offset for the value-laplacian structure constant matrix
Cb Bugs
Cb   MPI is not implemented
Cb   All sites are assumed non equivalent
Cb   Program does not check if ehvl used to build struxs are same as those in ctrl.*
Cb   in the reading mode
Cr Remarks
Cr   Structure constant matrix sg(nl,nl,2,2,nttabg) is treated as as 1D
Cr   array sg(ndim) throughout.
Cr
Cr   If file read and 8th bit of lio is on, allocates memory for and reads 
Cr   iax, ntab, ntabg, and sg
Cr
Cr   Input parameters nl,nbas, and bas must match file
Cr   for iosg to complete file read.
Cr
Cr   Content of file STRG (including header), version 2
Cr     record  content
Cr       1     -77 vsn ihold 
Cr       2     ivl,nttab,nttabg,nl,nbas,ckbas
Cr     .... end of header                      
Cr       3     ntab,ntabg,iax,sg,-999
Cr                                    
Cr   In the reading mode program stops if either of the following happens:                                    
Cr       first number is not -77
Cr       version number (vsn) is not 2
Cr       if ivl, nl, nbas or chbas do not match those in the input file
Cr       last number is not -999 (to verify if size of sg is correct)
Cu Updates
Cu  24 Jan 08 (S. Lozovoi) adapted from iostr.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character name*(*)
      integer lio,nttabg,ivl,nl,nbas,oiax,onpr,onprg,osg
      double precision ckbas
C ... Local parameters
      character locnam*8
      logical iosg1,iosg2,lerr
      integer nlr,nbr,ivlr,ifi,ipr,nttab,lgunit,niax,
     .  lio01,getdig,nds,lior,ndim,ng,ig          
      integer vsn
      double precision tiny,ckbasr
      parameter (tiny=1d-8, niax=10, vsn=1)
C ... Heap
      integer w(1)
      common /w/ w

C --- Setup information ---
      call getpr(ipr)
      locnam = name
      lio01 = mod(lio/1,100)

C --- File write ---
      if (mod(lio,2) .ne. 0) then

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        ndim = (2*nds)**2*nttabg

C   ... Write header
        nttab = w(onpr+nbas)
        if (.not. iosg1(lio,name,ivl,nttab,nl,nbas,
     .    nttabg,ckbas,ifi)) goto 99
C   ... Write strux
          if (.not. iosg2(-ifi,nbas,nttab,ndim,
     .        w(onpr),w(onprg),w(oiax),w(osg))) goto 99

        if (ipr .ge. 40) call awrit3(' IOSG: wrote to file '
     .    //locnam//'%a %i sites'//
     .    ' out of %i total',
     .    ' ',80,lgunit(1),nttabg,nttab,0)
        call fclose(ifi)

C --- File read ---
      else
C   ... Open file, and read header
        lior = 0
        lerr = .not. iosg1(lior,name,ivlr,nttab,nlr,nbr,
     .    nttabg,ckbasr,ifi)

C   ... If failed to read header, give up
        if (lerr) goto 99

C       Copy file parameters to those in argument list
        if (getdig(lio01,1,2) .eq. 1) then
          nl = nlr
          nbas = nbr
          ckbas = ckbasr
        endif

        call isanrg(ivlr,ivl,ivl,'file: IOSG:','ivl',.true.)
        call isanrg(nlr,nl,nl,'file: IOSG:','nl',.true.)
        call isanrg(nbr,nbas,nbas,'file: IOSG:','nbas',.true.)
        call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSG:','cksum',
     .    .true.)

C   ... If only header is to be read, exit
        if (getdig(lio01,2,2) .eq. 1) then
          iosg = .true.
          return
        endif

C   ... Allocate memory for arrays
        nds = nl*nl 
        ndim = (2*nds)**2*nttabg 
        if (getdig(lio01,3,2) .eq. 1) then
          call defi(oiax,niax*nttab)
          call defi(onpr,nbas+1)
          call defi(onprg,nbas)
          call defdr(osg,ndim)
        endif

C   ... Read strux
        lerr = iosg2(ifi,nbas,nttab,ndim,
     .                 w(onpr),w(onprg),w(oiax),w(osg))
        if (.not. lerr) goto 99

C ... Sanity check
        ng = 0
        do ig = 1, nbas
          ng = ng + w(onprg+ig-1)
        enddo
        call isanrg(ng,nttabg,nttabg,'file: IOSG:','nttabg',.true.)

C ... Normal exit
        if (ipr .ge. 10) call awrit3(' IOSG: read %i sites'//
     .    ' out of %i total',' ',80,
     .    lgunit(1),nttabg,nttab,0)
      endif
      iosg = .true.
      return

C ... Error handling
   99 continue
      if ((lio/10) .eq. 0) call rxs('IOSG: file mismatch, file ',name)
      if (ipr .ge. 10)
     .  print *, 'IOSG (warning): failed to read file ',name
      iosg = .false.

      end

      logical function iosg1(lio,name,ivl,nttab,nl,nbas,nttabg,
     .  ckbas,ifi)
C- Low level structure constants file open and header I/O
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  lio   :1s digit 0 for read, otherwise write
Cio        :On file READ, returns with file contents of lio
Cio  name  :file name
Ci   ivl   :specifies functions used to built the value-laplacian strux
Cio  nttab :total number of pairs in neighbor table iax (pairc.f)
Cio  nl    :global maximum l + 1
Cio  nbas  :number of atoms in the basis
Cio  nttabg:total number of pairs in the value-laplacian cluster (pairg.f)
Cio  ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Co Outputs
Co   ifi   :file logical unit number
Co   iosg1:T if read successful, false otherwise
Cr Remarks
Cr   Reads strux header data for the value-laplacian basis from file, 
Cr   corresponds to iostr version 2.
Cu Updates
Cu   24 Jan 08 Adapted from iostr1 (iostr.f) 
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character name*(*)
      integer lio,ivl,nttab,nttabg,nl,nbas
      double precision ckbas
C ... Local parameters
      integer ifi,fopnx
      integer i,vsn,rdvsn
      character strn*120
      parameter (vsn=2)
      
      ifi = fopnx(name,100,16+8+4+0,-1)
      rewind ifi

C ... File READ
      if (mod(lio,2) .eq. 0) then

C ... Check for version number.  If nonexistent, stop
        read(ifi,end=99,err=99) i,rdvsn,lio
        if (i .ne. -77) call rxs('IOSG: file mismatch, file ',name)
        if (rdvsn .eq. 2) then
          read(ifi,end=99,err=99) ivl,nttab,nttabg,nl,nbas,ckbas
        else
C ... Only ver=2 is currently permitted. Stop
          call rxi('iosg: version is not supported. ver = ',rdvsn)
        endif

        strn = ' IOSG : read file='//name//
     .    '%a lio=%i nl=%i nbas=%i nttab=%i nttabg=%i'
        i = lio

C ... File WRITE
      else
        i = lio-1
        write(ifi) -77,vsn,i
        write(ifi) ivl,nttab,nttabg,nl,nbas,ckbas

        strn = ' IOSG : write file='//name//
     .    '%a lio=%i nl=%i nbas=%i nttab=%i nttabg=%i'

      endif

      call info5(110,1,0,strn,i,nl,nbas,nttab,nttabg)

      iosg1 = .true.
      return

   99 iosg1 = .false.
      end

      logical function iosg2(ifi,nbas,nttab,ndim,ntab,ntabg,
     .  iax,sg)
C- Low level structure constants I/O
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   nbas  :number of atoms in the basis
Ci   nttab :total number of pairs in neighbor table iax (pairc.f)
Ci   ndim  :size of sg
Cio Inputs/Outputs
Cio  ntab  :ntab(ib) no. pairs in neighbor table iax preceding ib (pairc.f)
Cio  ntabg :ntabg(ib) no. pairs for the value-laplacian basis in cluster ib (pairg.f)
Cio  iax   :array of parameters containing info about each pair
Cio  sg    :real-space structure constant matrix for the value-laplacian basis
Co Outputs
Co   iosg2:T if read successful, false otherwise
Cr Remarks
Cr   nctrl is placed at the end of the file to ensure that array dimensions match
Cu Updates
Cu   24 Jan 08 Adapted from iostr4 (iostr.f) 
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,nbas,nttab,ndim
      integer niax
      parameter (niax=10)
      integer iax(niax,nttab),ntab(nbas+1),ntabg(nbas)
      double precision sg(ndim)
C ... Local parameters
      integer nctrl,ncr
      parameter (nctrl=-999)

      if (ifi .gt. 0) then
        read(ifi,end=99,err=99) ntab,ntabg,iax,sg,ncr
        if (ncr .ne. nctrl) 
     .    call rx('iosg2(r): array dimension mismatch')
      else
        write(-ifi) ntab,ntabg,iax,sg,nctrl
      endif

      iosg2 = .true.
      return

   99 continue
      iosg2 = .false.
      end

