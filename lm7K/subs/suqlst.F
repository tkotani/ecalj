c$$$#define MULLIKEN 1
c$$$      subroutine suqlst(sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,
c$$$     .qp, onespdummy) !onesp removed 
c$$$      use m_globalvariables
c$$$C- Set up a list of q-points in various modes for energy bands
c$$$C  See also entry suqlsw(nband,jsp,nsp,evl)
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   sopts :character string describing special options that can
c$$$Ci         :selects the mode for which qp are read.
c$$$Ci         :*default mode: the qp file consists of a list of
c$$$Ci         :               lines and the number of qp in each line.
c$$$Ci         :*list mode   : the qp file consists of a list of qp
c$$$Ci         :*contour mode: input for specifying qp on a uniform 2D
c$$$Ci                       : mesh, for contour plots.
c$$$Ci         :See Remarks for the syntax of each mode.
c$$$Ci
c$$$Ci         :Options are separated by delimiters; the first character is
c$$$Ci         :the delimiter.  The following list the strings declaring
c$$$Ci         :options assuming the delimiter is '/'. (Space as the first
c$$$Ci         :implies that there are no options)
c$$$Ci         :/fn=fnam     read qp from file 'fnam'
c$$$Ci         :/con         input file for contour plot mode (mode 3)
c$$$Ci         :/qp          input file specifies a list of qp (mode 2)
c$$$Ci         :             (default assumes a list of symmetry lines.)
c$$$Ci         :/long        write bands with extra digits precision
c$$$Ci         :             (has no effect for symmetry line mode)
c$$$Ci         :/spin1       generate bands only for first spin
c$$$Ci         :/lst=list    write only those bands specified in a list.
c$$$Ci         :             For syntax of list, see slatsm/mkilst.f
c$$$Ci         :/evn=#       keep track of smallest, largest eval for
c$$$Ci                       #th band, and printout at close.
c$$$Ci         :/ef=#        change efermi to #.
c$$$Ci         :Example: --band~long~qp~lst=2:5
c$$$Ci   iop   :options passed by the caller
c$$$Ci         :1s digit is suppress actions for parallel mode
c$$$Ci         : mode 1:
c$$$Ci         : 1: suppress writing line header info file
c$$$Ci         : 2: Return total number of qp to be generated in all lines
c$$$Ci         :    and also suppress writing line header info file
c$$$Ci         :    Thus iop=2 => all modes return total no qp to be generated
c$$$Ci         : mode 3:
c$$$Ci         : nonzero: do not allocate evsav until iq>nq
c$$$Ci
c$$$Ci   nband :(suqlst) maximum number of energy bands to write
c$$$Ci         :(suqlsw) actual number of energy bands to write
c$$$Ci
c$$$Ci   efermi:Fermi energy (written to bnds file)
c$$$Ci
c$$$Ci   nsp   :2 for spin-polarized case, otherwise 1
c$$$Ci         :NB used only in file write to indicate how many times
c$$$Ci         :nsp is ALTERED to nsp=1 if spin1 option is set
c$$$Ci
c$$$Ci   evl   :eigenvalues for current qp (used only in suqlsw)
c$$$Cio Inputs/Outputs (see Remarks)
c$$$Cio  nq    :Input nq=0 :if zero, flags suqlst to set up bands mode:
c$$$Cio                    :some initializations; then returns ready
c$$$Cio                    :for first block of qp
c$$$Cio        :           :if nonzero, setup for for next qp
c$$$Cio        :Output nq  :if zero, suqlst has no more qp to calculate.
c$$$Cio                    :if nonzero, nq is the number of k-points in the
c$$$Cio                    :current block.  Caller then calls suqlst nq
c$$$Cio                    :times (one for each k-point), which will return
c$$$Cio                    :qp for each k-point.  After exactly nq calls,
c$$$Cio                    :suqlst will start another block, if it exists.
c$$$Cio                    :See Remarks for schematic of calling sequence.
c$$$Co Outputs
c$$$Co   nfbn  :(color weights) number of elements in iblst
c$$$Co   iblst :(color weights) list of orbital indices for color weights:
c$$$Co         :color weight is the sum of contributions from eigenvector
c$$$Co         :components in iblst
c$$$Co   qp    :not set, if suqlst is called starting a new block.
c$$$Co         :(see Remarks).  Else, k-point at which to generate bands
c$$$Co    onesp(removed by takao) :
c$$$Co    if spin1 flag is encountered, onesp is set to 1 and nsp is set to 1
c$$$Cl Local variables
c$$$Cl   iq    :current qp in this block
c$$$Cl   mode  :1 symmetry-line mode
c$$$Cl         :2 list-of-qp mode
c$$$Cl         :3 contour mode
c$$$Cl   q1    :starting qp for symmetry mode; only meaningful in that mode
c$$$Cl   q2    :ending   qp for symmetry mode; only meaningful in that mode
c$$$Cl   nevn  :(optional) band index; routine monitors largest, smallest value
c$$$Cl         :for that index.
c$$$Cl   evmnn :smallest value found for qp of specified band index
c$$$Cl   evmxn :largest  value found for qp of specified band index
c$$$Cl   ifiq  :file logical unit for input qp file
c$$$Cl   ifib  :file logical unit for output bands file
c$$$Cl   ql    :local copy of current qp
c$$$Cl   oqp
c$$$Cr Remarks
c$$$Cr   suqlst is designed to be called to generate qp in groups or
c$$$Cr   blocks.  The calling sequence is:
c$$$Cr
c$$$Cr     nq = 0  <- flags that first call to suqlst, to set up mode
c$$$Cr     do  iblock = 1, forever
c$$$Cr       This call generates nq, the number of points in this block
c$$$Cr       call suqlst(nband,efermi,nsp,evl,nq,qp,onesp)
c$$$Cr       if (nq .eq. 0) stop
c$$$Cr       do  iq = 1, nq
c$$$Cr         This call generates qp for current block
c$$$Cr         call suqlst(ndimh,ef0,nsp,w,nkp,qp,onesp) <- returns qp
c$$$Cr         do  isp = 1, nsp
c$$$Cr         call suqlsw(ndimh,qp,evl(1,isp)) <- saves evl for this qp
c$$$Cr                                             (call is optional)
c$$$Cr         enddo
c$$$Cr       enddo
c$$$Cr     enddo
c$$$Cr
c$$$Cr   The following modes are implemented:
c$$$Cr     mode=1 reads qp from syml file, and generates qp along each
c$$$Cr            specified symmetry line.  Structure of qp file:
c$$$Cr            file has one line for each symmetry line as follows:
c$$$Cr               nq      q1x   q1y   q1z      q2x   q2y   q2z
c$$$Cr               ...
c$$$Cr            the line entries have meanings:
c$$$Cr            --# qp-   ---starting qp---    --- ending qp ---
c$$$Cr            Any line with nq=0 implies no more lines.
c$$$Cr     mode=2 reads qp from specified file and generates qp for each
c$$$Cr            specified qp.  File consists sets of qpx,qpy,qpz for each
c$$$Cr            qp sought.  Typically a file would consists of lines like
c$$$Cr               q1x   q1y   q1z
c$$$Cr               q2x   q2y   q2z
c$$$Cr               ...
c$$$Cr     mode=3 generates qp for a uniform mesh in a plane (contour plot)
c$$$Cr            The file supplies information describing a rectangle in
c$$$Cr            the Brillouin zone.  It consists of a single line,
c$$$Cr            which contains the following:
c$$$Cr             v1    range  n     v2    range  n   height  list-of-bands
c$$$Cr
c$$$Cr            v1 and v2 are two vectors specifying the plane of the
c$$$Cr            contour.  range and n (one each for v1 and v2) are the
c$$$Cr            starting and final amplitudes of those vectors, and the
c$$$Cr            the number of points within the vector.  list-of-bands
c$$$Cr            is a list of integers which specify which bands are to
c$$$Cr            be written to the output file.  'height' is the 'z' axis.
c$$$Cr            For example,
c$$$Cr             v1    range  n     v2    range  n   height  list-of-bands
c$$$Cr            1 0 0  -1 1   51   0 1 0  -1 1   51   0.00    4,5
c$$$Cr            creates a file of 51x51 points, with the four corners
c$$$Cr            (-1,-1,0),  (1,-1,0),  (-1,1,0),  (1,1,0)
c$$$Cr            going through the gamma-point. Two bands (4,5) are stored.
c$$$Cu Updates
c$$$Cu   08 Jul 08 Extend to case where number of bands can be q dependent
c$$$Cu             modes 1,2: suqlsw writes out number of bands with qp
c$$$Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1 
c$$$c                                                          -->forbidden 2014takao
c$$$Cu   02 Jul 06 Color mode extended to two colors
c$$$Cu   02 Jul 06 New color mode (one color weight only)
c$$$Cu   14 Feb 05 contour mode saves both spins in spin-polarized case
c$$$Cu   20 Oct 03 suqlst works properly in contour mode
c$$$Cu   28 Aug 01 added ef switch
c$$$Cu   23 Jan 01 first written
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      integer nq,nband,nsp,nfbn(2),ifblst(nband,*),onespdummy,iop!onesp
c$$$      double precision efermi,qp(3),evl(nband)
c$$$C ... Local variables
c$$$      character strn*120, strn2*120, dc*1, fn*120, prfmt*40
c$$$      logical rdstrn,fmlong
c$$$      integer fopno,i,iprint,j,j1,j2,k,nglob,fopnn,rdm,
c$$$     .parg,a2vec,iv(7),stdo,jsp,op1
c$$$      double precision xx,xx1,xx2,sclp,plndst,xxv(7)
c$$$      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c$$$      integer oqp,nqall,imode
c$$$      double precision x1,x2,y1,y2,evmxn,evmnn,
c$$$     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
c$$$      real(8),allocatable :: evsav(:,:,:,:)
c$$$C ... MPI
c$$$      integer procid,master,mpipid
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$ctakao
c$$$      real(8)::qlx(3)
c$$$C ... Heap
c$$$c      integer w(1)
c$$$c      common /w/ w
c$$$      common /suqlsd/
c$$$     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
c$$$     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c$$$c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
c$$$      save evsav
c$$$      real(8),allocatable:: w_oqp(:,:)
c$$$
c$$$      procid = mpipid(1)
c$$$      master = 0
c$$$Changenglob      stdo = nglob('stdo')
c$$$      stdo = globalvariables%stdo
c$$$      op1 = mod(iop,10)
c$$$      nqall = 0
c$$$
c$$$C --- First call ... setup and parse options ---
c$$$      if (nq .eq. 0) then
c$$$C   ... Defaults
c$$$        fmlong = .false.
c$$$        mode = 1
c$$$        fn = 'syml'
c$$$        nblst = 0
c$$$        nevn = 0
c$$$        evmxn = -99d9
c$$$        evmnn = 99d9
c$$$        nfbn(1) = 0
c$$$        nfbn(2) = 0
c$$$
c$$$        dc = sopts(1:1)
c$$$        if (dc .ne. ' ') then
c$$$C   ... Return here to resume parsing for arguments
c$$$          j2 = 0
c$$$   10     continue
c$$$          j2 = j2+1
c$$$          if (sopts(j2:j2) .eq. dc) goto 10
c$$$          j1 = min(len(sopts),j2)
c$$$          call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$          if (j2 .ge. j1) then
c$$$            if (.false.) then
c$$$            elseif (sopts(j1:j1+2) .eq. 'fn=')  then
c$$$              if (j1+3 .gt. j2) call rx('suqlst: bad file name')
c$$$              fn = sopts(j1+3:j2)
c$$$            elseif (sopts(j1:j2) .eq. 'qp')  then
c$$$              mode = 2
c$$$            elseif (sopts(j1:j2) .eq. 'spin1')  then
c$$$               call rx('suqlst:spin removed! takao')
c$$$c              onesp = 1
c$$$c              nsp = 1
c$$$            elseif (sopts(j1:j2) .eq. 'con')  then
c$$$              mode = 3
c$$$            elseif (sopts(j1:j2) .eq. 'long')  then
c$$$              fmlong = .true.
c$$$            elseif (sopts(j1:j1+3) .eq. 'col=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, col=..')
c$$$              call mkils0(sopts(j1+4:j2),nfbn,ifblst)
c$$$              call mkilst(sopts(j1+4:j2),nfbn,ifblst)
c$$$              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
c$$$            elseif (sopts(j1:j1+4) .eq. 'col2=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, col2=..')
c$$$              call mkils0(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
c$$$              call mkilst(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
c$$$              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
c$$$            elseif (sopts(j1:j1+3) .eq. 'lst=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, lst=..')
c$$$              call mkils0(sopts(j1+4:j2),nblst,iblst)
c$$$              if (nblst .gt. 100) call rx('increase size of iblst')
c$$$              call mkilst(sopts(j1+4:j2),nblst,iblst)
c$$$            elseif (sopts(j1:j1+2) .eq. 'ef=')  then
c$$$              j = 0
c$$$              i = parg('ef=',4,sopts(j1:),j,len(sopts(j1:)),
c$$$     .        dc//' ',1,1,i,efermi)
c$$$              if (i .le. 0) call rxs(
c$$$     .        'suqlst: failed to parse string for fermi level:  ',
c$$$     .        sopts(j1:))
c$$$            elseif (sopts(j1:j1+3) .eq. 'evn=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list')
c$$$              i = j1+3
c$$$              xxv(1) = 0
c$$$              j = a2vec(sopts,j2,i,2,dc//' ',2,3,1,iv,nevn)
c$$$              if (j .ne. 1 .or. nevn .gt. nband)
c$$$     .        call rx('suqlst: bad value for evn')
c$$$            else
c$$$              call rxs('suqlst: failed to parse argument, ',sopts(j1:j2))
c$$$            endif
c$$$            goto 10
c$$$          endif
c$$$        endif
c$$$
c$$$        if (nfbn(1) .eq. 0) then
c$$$          call info2(20,0,0,' suqlst:  generate bands, mode %i',mode,0)
c$$$        elseif (nfbn(2) .eq. 0) then
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,
c$$$     .    ' suqlst:  generate bands with color weights, mode %i.'//
c$$$     .    ' %N%10f%i components: '//strn//'%a',mode,nfbn)
c$$$        else
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,' suqlst:  generate bands '//
c$$$     .    'with two color weights, mode %i.',mode,0)
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,
c$$$     .    '%10f%i components, first color:  '//strn//'%a',nfbn,0)
c$$$          call ilst2a(ifblst(1,2),nfbn(2),strn)
c$$$          call info2(20,0,0,
c$$$     .    '%10f%i components, second color: '//strn//'%a',nfbn(2),0)
c$$$        endif
c$$$        if (nfbn(1) .eq. 0 .and. nfbn(2) .ne. 0) call info0(20,0,0,
c$$$     .  ' suqlst:  (warning) col2 specified, but col not')
c$$$
c$$$C   ... Open qp file
c$$$        if (procid .eq. master) then
c$$$          ifiq = fopno(fn)
c$$$          rewind ifiq
c$$$
c$$$C   ... open bands file
c$$$          ifib = fopnn('BNDS')
c$$$          rewind ifib
c$$$
c$$$C   ... Write header
c$$$          if (mode .eq. 1) then
c$$$            i = nblst
c$$$            if (nblst .eq. 0) i = nband
c$$$            if (nfbn(1) .eq. 0) then
c$$$C           Use separate format statment to circumvent gfortran bug
c$$$              write(ifib,335) i,efermi,0
c$$$  335         format(i5,f10.5,i6)
c$$$            elseif (nfbn(2) .eq. 0) then
c$$$              call ilst2a(ifblst,nfbn,strn)
c$$$              call strip(strn,j1,j2)
c$$$              write(ifib,336) i,efermi,1,strn(j1:j2)
c$$$  336         format(i5,f10.5,i6:'  col= ',a:'  col2= ',a)
c$$$            else
c$$$              call ilst2a(ifblst,nfbn,strn)
c$$$              call strip(strn,j1,j2)
c$$$              call ilst2a(ifblst(1,2),nfbn(2),strn2)
c$$$              call strip(strn2,j,k)
c$$$              write(ifib,336) i,efermi,2,strn(j1:j2),strn2(j:k)
c$$$            endif
c$$$          endif
c$$$        endif
c$$$
c$$$C   ... Other initializations
c$$$        iq = 0
c$$$        nq = -1
c$$$      endif
c$$$
c$$$C --- Setup for a new block of k-points, depending on mode ---
c$$$C     This branch occurs on completion of the last qp of the current block
c$$$C     which is marked by iq>nq
c$$$C     At the completion of this block: nq must be computed and:
c$$$C     (mode=1) q1,q2 set up.
c$$$C              Note: if 1s digit of iop is set in this mode,
c$$$C              this branch returns sum of all qp in all lines.
c$$$C              No setup for q1,q2; no
c$$$C     (mode=2) oqp allocated and loaded
c$$$C     (mode=3) nblst = number of bands to save
c$$$      if (iq .gt. nq) then
c$$$        iq = 1
c$$$C       Bands along specified symmetry lines
c$$$        if (mode .eq. 1) then
c$$$  725     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
c$$$          if (strn(1:1) .eq. '#') goto 725
c$$$C         Exit if first entry in line is zero
c$$$          i = 0
c$$$          xxv(1) = 0
c$$$          i = a2vec(strn,len(strn),i,4,', ',2,3,1,iv,xxv)
c$$$          if (i .eq. 1 .and. xxv(1) .eq. 0) goto 999
c$$$          i = 0
c$$$          i = a2vec(strn,len(strn),i,4,', ',2,3,7,iv,xxv)
c$$$          if (i .ne. 7 .and. iprint().ge.10) then
c$$$            write(stdo,
c$$$     .      '(/'' suqlst (warning) skipping line:''/''  '',a)') strn
c$$$            goto 725
c$$$          endif
c$$$          nq = xxv(1)
c$$$C         1 qp is nonsensical for a line
c$$$          if (nq .eq. 1) nq = 0
c$$$C         No qp: exit
c$$$          if (nq .le. 0) goto 999
c$$$C         setup q1,q2
c$$$          call dcopy(3,xxv(2),1,q1,1)
c$$$          call dcopy(3,xxv(5),1,q2,1)
c$$$          if (iprint().ge.10) write(stdo,785) nq,q1,q2
c$$$  785     format(/' suqlst:  nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
c$$$C         Write line header information to disk
c$$$          if (op1 .eq. 0) then
c$$$            write(ifib,337) nq*nsp
c$$$  337       format(2i5)
c$$$C         Accumulate all qp and cycle until all lines are read
c$$$          elseif (op1 .eq. 2) then
c$$$            nqall = nqall + nq
c$$$            goto 725
c$$$          endif
c$$$C       Bands for a list of specified qp
c$$$        elseif (mode .eq. 2) then
c$$$C         Only one block for this mode.
c$$$C         Flag that prior block already completed: nq>0
c$$$          if (nq .gt. 0) goto 999
c$$$          nq = 0
c$$$          j = 3
c$$$          if (rdm(ifiq,10000,0,' ',xx,j,nq) .ne. 1) call
c$$$     .    rx('suqlst:  bad qp file')
c$$$
c$$$c          call defrr(oqp,3*nq)
c$$$          allocate(w_oqp(3,nq))
c$$$
c$$$          rewind ifiq
c$$$          j = rdm(ifiq,10000,3*nq,' ',w_oqp,j,nq)
c$$$          call awrit1('%N suqlst: read %i qp from file '//fn//
c$$$     .    '%a',' ',80,stdo,nq)
c$$$          if (j .ne. 1) call rx('suqlst: failed to read qp')
c$$$          i = nblst
c$$$          if (nblst .eq. 0) i = nband
c$$$          if (nfbn(1) .ne. 0) i = i*2
c$$$          call awrit2('%% rows %i cols %i',' ',80,ifib,nq,i+4)
c$$$C       Bands on a uniform mesh in a specified plane (contour plot)
c$$$        elseif (mode .eq. 3 .and. nq .eq. -1) then
c$$$  825     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 998
c$$$          if (strn(1:1) .eq. '#') goto 825
c$$$          call words(strn,i)
c$$$          if (i .ne. 14) call rxi('suqlst con mode: expected 14 '//
c$$$     .    'arguments from input file but read',i)
c$$$          iblst(1) = -1
c$$$          nblst = 1
c$$$          backspace ifiq
c$$$          read(ifiq,*) vectx,x1,x2,nqx,vecty,y1,y2,nqy,plndst
c$$$          call word(strn,14,j1,j2)
c$$$          call mkilss(11,strn(j1:j2),nblst,iblst)
c$$$          if (nblst .le. 0) call rx('suqlst: no bands in list')
c$$$          call dscal(3,1/dsqrt(sclp(vectx,vectx)),vectx,1)
c$$$          call dscal(3,1/dsqrt(sclp(vecty,vecty)),vecty,1)
c$$$C         Subtract from vecty projection onto vectx
c$$$C         call daxpy(3,-sclp(vectx,vecty),vectx,1,vecty,1)
c$$$          call cross(vectx,vecty,vect3)
c$$$          call dscal(3,plndst/dsqrt(sclp(vect3,vect3)),vect3,1)
c$$$          nq = nqx*nqy
c$$$          if (iprint() .ge. 10) then
c$$$            write(stdo,717) vectx,x1,x2,nqx,plndst,vecty,y1,y2,nqy,nq
c$$$  717       format(' vx=',3f9.6,'  x1,x2=',2f9.6,'  nx=',i3,'   h=',f9.6/
c$$$     .      ' vy=',3f9.6,'  y1,y2=',2f9.6,'  ny=',i3,'  np=',i5)
c$$$            write(strn,'('' save %i bands: %'',i2,'':1i'')') nblst
c$$$            if (strn(17:18) .eq. '% ') strn(17:18) = ' %'
c$$$            call awrit2(strn,strn,80,stdo,nblst,iblst)
c$$$            if (dabs(sclp(vectx,vecty)) .gt. 1d-6 .and. iprint().gt.20)
c$$$     .      print *,'suqlst (warning): input plane vectors not orthogonal'
c$$$          endif
c$$$
c$$$          if (iop .eq. 0) then
c$$$            allocate(evsav(nqx,nqy,nblst,nsp))
c$$$          endif
c$$$
c$$$C       Contour plot, cleanup.
c$$$C       Note: cleanup handled by suqlsw when last qp is called.
c$$$        elseif (mode .eq. 3) then
c$$$C          if (op1 .ne. 0) then
c$$$C            allocate(evsav(nqx,nqy,nblst,nsp))
c$$$C            return
c$$$C          endif
c$$$          call rx('suqlst: caller should never reach this branch')
c$$$        endif
c$$$
c$$$C --- Generate qp for this iq, depending on mode ---
c$$$      else
c$$$        if (mode .eq. 1) then
c$$$          xx = dble(iq-1)/dble(nq-1)
c$$$          qp(1) = xx*q2(1) + (1-xx)*q1(1)
c$$$          qp(2) = xx*q2(2) + (1-xx)*q1(2)
c$$$          qp(3) = xx*q2(3) + (1-xx)*q1(3)
c$$$        elseif (mode .eq. 2) then
c$$$          call dpscop(w_oqp,qp,3,iq*3-2,1,1d0)
c$$$        elseif (mode .eq. 3) then
c$$$C         Inner Loop:  excursions in y; outer loop: excursions in x
c$$$          j = mod(iq-1,nqy)
c$$$          i = (iq-1-j)/nqy
c$$$          if (nqx .le. 1) then
c$$$            xx1 = x1
c$$$          else
c$$$            xx1 =i*(x2-x1)/(nqx-1)+x1
c$$$          endif
c$$$          if (nqy .le. 1) then
c$$$            xx2 = y1
c$$$          else
c$$$            xx2 =j*(y2-y1)/(nqy-1)+y1
c$$$          endif
c$$$          do  93  k = 1, 3
c$$$            qp(k) = xx1*vectx(k) + xx2*vecty(k) + vect3(k)
c$$$   93     continue
c$$$          if (j .eq. 0 .and. iprint().ge.20) write(stdo,718) i+1,nqx,qp
c$$$  718     format(' line',i3,' of',i3,'   q(1)=',3f10.6)
c$$$        else
c$$$          call rx('suqlst: bad mode')
c$$$        endif
c$$$        iq = iq+1
c$$$C       Hold onto local copy of qp
c$$$        call dcopy(3,qp,1,ql,1)
c$$$      endif
c$$$      if(allocated(w_oqp)) deallocate(w_oqp)
c$$$      return
c$$$
c$$$C --- No more qp blocks: cleanup ---
c$$$  999 continue
c$$$      nq = 0
c$$$      if (nevn .ne. 0) call awrit3(' eval no. %i:  minimum eval'//
c$$$     .' = %;8F  maximum eval = %;8F',' ',80,stdo,nevn,evmnn,evmxn)
c$$$
c$$$      if (mode .eq. 1) then
c$$$        if (op1 .eq. 0) then
c$$$          write(ifib,337) 0
c$$$        elseif (op1 .eq. 2) then
c$$$          nq = nqall
c$$$        endif
c$$$      elseif (mode .eq. 3) then
c$$$        call rx('not ready')
c$$$C       call xxxbnd(w(oev),nblst,nqx,nqy,ifib)
c$$$C       return
c$$$      endif
c$$$      if(allocated(w_oqp)) deallocate(w_oqp)
c$$$      return
c$$$
c$$$C --- Error exit ---
c$$$  998 call rxs('suqlst: failed to read file contents, file ',fn)
c$$$
c$$$      entry suqlsm(imode)
c$$$C- Return qlist mode
c$$$      imode = mode
c$$$      return
c$$$
c$$$      entry suqlsw(nband,jsp,nsp,evl)
c$$$C- Write or store the energy bands to file for this qp
c$$$
c$$$      if (mode .eq. 1) then
c$$$        if (nblst .eq. 0) then
c$$$          prfmt = '(3f10.5,i6/(10f8.4))'
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.5/(10f8.4))'
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong .and. nblst .eq. 0) then
c$$$          prfmt = '(3f15.10,i6/(5f15.10))'
c$$$        elseif (nblst .eq. 0) then
c$$$          prfmt = '(3f10.6,i6/(8f10.6))'
c$$$        elseif (fmlong .and. nblst .gt. 0) then
c$$$          prfmt = '(3f15.10/(5f15.10))'
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.6/(8f10.6))'
c$$$C         if (nblst .le. 5) prfmt = '(8f10.6)'
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$C ... Keep running tab on smallest,largest eval
c$$$      if (nevn .ne. 0) then
c$$$        evmxn = max(evmxn,evl(nevn))
c$$$        evmnn = min(evmnn,evl(nevn))
c$$$      endif
c$$$
c$$$
c$$$      entry suqlsw2(nband,jsp,nsp,evl,qlx) !takao added ql
c$$$C- Write or store the energy bands to file for this qp
c$$$      ql=qlx
c$$$      if (mode .eq. 1) then
c$$$        if (nblst .eq. 0) then
c$$$          prfmt = '(3f10.5,i6/(10f8.4))'
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.5/(10f8.4))'
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong .and. nblst .eq. 0) then
c$$$          prfmt = '(3f15.10,i6/(5f15.10))'
c$$$        elseif (nblst .eq. 0) then
c$$$          prfmt = '(3f10.6,i6/(8f10.6))'
c$$$        elseif (fmlong .and. nblst .gt. 0) then
c$$$          prfmt = '(3f15.10/(5f15.10))'
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.6/(8f10.6))'
c$$$C         if (nblst .le. 5) prfmt = '(8f10.6)'
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$c  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$C ... Keep running tab on smallest,largest eval
c$$$      if (nevn .ne. 0) then
c$$$        evmxn = max(evmxn,evl(nevn))
c$$$        evmnn = min(evmnn,evl(nevn))
c$$$      endif
c$$$      end
c$$$
c$$$      subroutine suqlse(nband,jsp,nsp,ndimhx,ifbn,nfbn,ifblst,ndlst,
c$$$     .evec,wk)
c$$$      use m_globalvariables
c$$$C- Write to file the projection of eigenvector subblock for this qp
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nband :number of energy bands to write
c$$$Ci   jsp   :current spin index (not used now)
c$$$Ci   nsp   :number of spins (not used now)
c$$$Ci   ndimhx:dimensions evec
c$$$Ci   ifbn  :index to color list (nfbn,ifblst)
c$$$Ci   nfbn  :number of elements for color weights projection
c$$$Ci   ifblst:list of elements for  color weights projection
c$$$Ci   ldlst :leading dimension of ifblst
c$$$Ci   evec  :eigenvectors
c$$$Ci   wk    :work array of same dimension as evl
c$$$Co Outputs
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr   Inefficient, but it works
c$$$Cu Updates
c$$$Cu   08 Jul 08 New argument ndlst so ifblst can be dimensioned
c$$$Cu             independently from nband
c$$$Cu   05 Jun 06 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer jsp,nsp,nband,ifbn,nfbn(2),ndlst,ifblst(ndlst,*)
c$$$      integer ndimhx
c$$$      double precision wk(1)
c$$$      double complex evec(ndimhx,ndimhx)
c$$$C ... Local parameters
c$$$      double complex zdotc
c$$$      character prfmt*40
c$$$      logical fmlong
c$$$      integer i,iprint,j,j1,k,nglob,stdo
c$$$      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c$$$c      integer oqp
c$$$      double precision x1,x2,y1,y2,evmxn,evmnn,
c$$$     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
c$$$
c$$$      common /suqlsd/
c$$$     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
c$$$     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c$$$c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
c$$$
c$$$      integer ipiv(ndimhx)
c$$$      real(8),allocatable :: evsav(:,:,:,:)
c$$$      complex(8),allocatable:: evecc(:,:),work(:,:)
c$$$
c$$$      allocate(evecc(ndimhx,ndimhx),work(ndimhx,ndimhx))
c$$$
c$$$#if MULLIKEN
c$$$      call zcopy(ndimhx**2,evec,1,evecc,1)
c$$$      call zgetrf(ndimhx,ndimhx,evecc,ndimhx,ipiv,j)
c$$$      if (j .ne. 0) call rx('mullmf: failed to generate overlap')
c$$$      call zgetri(ndimhx,evecc,ndimhx,ipiv,work,ndimhx**2,j)
c$$$#endif
c$$$      do  i = 1, ndimhx
c$$$        wk(i) = 0
c$$$        do  j = 1, nfbn(ifbn)
c$$$          k = ifblst(j,ifbn)
c$$$          if (k .le. 0 .or. k .gt. ndimhx) then
c$$$            call fexit2(-1,111,' Exit -1 : suqlst: component %i'//
c$$$     .      ' outside range (1:%i)',k,ndimhx)
c$$$          endif
c$$$#if MULLIKEN
c$$$          wk(i) = wk(i) + evecc(i,k)*evec(k,i)
c$$$#else
c$$$          wk(i) = wk(i) + dconjg(evec(k,i))*evec(k,i)
c$$$#endif
c$$$        enddo
c$$$
c$$$#ifndef MULLIKEN
c$$$        wk(i) = wk(i) / zdotc(ndimhx,evec(1,i),1,evec(1,i),1)
c$$$        if (wk(i) .lt. 0 .or. wk(i) .gt. 1+1d-14) then
c$$$          call rx('bug in suqlse')
c$$$        endif
c$$$#endif
c$$$
c$$$      enddo
c$$$
c$$$      if (mode .eq. 1) then
c$$$        prfmt = '(3f10.5/(10f8.4))'
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, (wk(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (wk(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong) then
c$$$          prfmt = '(5f15.10)'
c$$$        else
c$$$          prfmt = '(8f10.6)'
c$$$C         if (nblst .gt. 0 .and. nblst .le. 5) prfmt = '(8f10.6)'
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) (wk(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) (wk(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        call rx('need copy weights into different place, mode=3')
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = wk(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$      deallocate(evecc,work)
c$$$
c$$$      end
c$$$
c$$$      subroutine suqlsr(mode,ifi,nsp,nbf,lde,nb1,nb2,ls1,ls2,nq,qp,eb)
c$$$C- Count qp, optionally read them and energy bands from file
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :Any combination of the following is allowed:
c$$$Ci         :1 count number of qp in file, return as nq
c$$$Ci         :  if 1s bit mode is 0, nq is an input,
c$$$Ci         :  the file nq cannot exceed input nq
c$$$Ci         :2 return qp read in file into array qp
c$$$Ci         :4 return bands read in file into array eb
c$$$Ci   ifi   :read from logical unit ifi
c$$$Ci   nsp   :number of spins in band file
c$$$Ci   nbf   :number of bands in band file
c$$$Ci   lde   :leading dimension of eb
c$$$Ci   nb1,nb2: read bands nb1..nb2 into eb
c$$$Ci   ls1   :(nsp=2): read first spin only
c$$$Ci   ls2   :(nsp=2): read second spin only
c$$$Cio Inputs/Outputs
c$$$Cio  nq    :numbefr of k-points (returned if 1s bit of mode is set)
c$$$Co Outputs
c$$$Co   qp    :k-points, returned if 4s bit of mode is set
c$$$Co   eb    :energy bands, returned if 4s bit of mode is set
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr
c$$$Cu Updates
c$$$Cu   06 Jun 07  First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer mode,ifi,nsp,nq,nbf,lde,nb1,nb2
c$$$      logical ls1,ls2
c$$$      double precision qp(3,nq),eb(lde,nsp,nq)
c$$$C ... Local parameters
c$$$      integer nq1,nql,mode1,mode2,mode4,iq1,i,nline,isp
c$$$      double precision ql(3),evl(3000)
c$$$      logical :: l_dummy_isanrg,isanrg
c$$$
c$$$      if (mode .eq. 0) return
c$$$      mode1 = mod(mode,2)
c$$$      mode2 = mod(mode/2,2)
c$$$      mode4 = mod(mode/4,2)
c$$$      if (mode4 .ne. 0) then
c$$$Ckino isanrg is logical function,         call isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
c$$$        l_dummy_isanrg=isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
c$$$     .  .true.)
c$$$Ckino isanrg is logical function,         call isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
c$$$        l_dummy_isanrg=isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
c$$$      endif
c$$$
c$$$C --- For each panel, do ---
c$$$      nql = 0
c$$$      nline = 0
c$$$   91 continue
c$$$      read(ifi,*) nq1
c$$$      if (nq1 .le. 0) goto 90
c$$$      isp = 0
c$$$      do  iq1 = 1, nq1
c$$$        isp = mod(isp,2)+1
c$$$        nql = nql+1
c$$$        if (mode1 .eq. 0 .and. nql .gt. nq) call rxi(
c$$$     .  'suqlsr: file more q-points than allocated: nqmx=',nq)
c$$$        read(ifi,*,END=999,ERR=999) ql(1),ql(2),ql(3)
c$$$        if (mode2 .ne. 0) call dcopy(3,ql,1,qp(1,nql),1)
c$$$        read(ifi,*,END=999,ERR=999) (evl(i),i=1,nbf)
c$$$        if (mode4 .ne. 0) then
c$$$C         Copy only if appropriate spin
c$$$          if (isp .eq. 1 .and. ls2) then
c$$$          elseif (isp .eq. 2 .and. ls1) then
c$$$          else
c$$$            call dcopy(nb2-nb1+1,evl(nb1),1,eb(1,1,nql),1)
c$$$          endif
c$$$        endif
c$$$      enddo
c$$$      nline = nline+1
c$$$      goto 91
c$$$C     End of loop over lines
c$$$   90 continue
c$$$
c$$$
c$$$      if (mode1 .eq. 0) then
c$$$        call info2(30,1,0,' suqlsr: found %i qp in %i lines from file',
c$$$     .  nql,nline)
c$$$      else
c$$$        nq = nql/nsp
c$$$        call info2(30,1,0,' suqlsr: read %i qp in %i lines from file',
c$$$     .  nql,nline)
c$$$      endif
c$$$
c$$$      if (mod(nql,2) .ne. 0 .and. nsp .eq. 2) then
c$$$        call info0(20,0,0, ' suqlsr (warning): '//
c$$$     .  'odd number of bands encountered but nsp=2')
c$$$      endif
c$$$
c$$$      return
c$$$  999 continue
c$$$      call rxi('suqlsr: failed to read bands file, nq=',nql)
c$$$
c$$$      end
c$$$
