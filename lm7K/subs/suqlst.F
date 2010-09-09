#define MULLIKEN 1
      subroutine suqlst(sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,
     .qp,onesp)
      use m_globalvariables
C- Set up a list of q-points in various modes for energy bands
C  See also entry suqlsw(nband,jsp,nsp,evl)
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :character string describing special options that can
Ci         :selects the mode for which qp are read.
Ci         :*default mode: the qp file consists of a list of
Ci         :               lines and the number of qp in each line.
Ci         :*list mode   : the qp file consists of a list of qp
Ci         :*contour mode: input for specifying qp on a uniform 2D
Ci                       : mesh, for contour plots.
Ci         :See Remarks for the syntax of each mode.
Ci
Ci         :Options are separated by delimiters; the first character is
Ci         :the delimiter.  The following list the strings declaring
Ci         :options assuming the delimiter is '/'. (Space as the first
Ci         :implies that there are no options)
Ci         :/fn=fnam     read qp from file 'fnam'
Ci         :/con         input file for contour plot mode (mode 3)
Ci         :/qp          input file specifies a list of qp (mode 2)
Ci         :             (default assumes a list of symmetry lines.)
Ci         :/long        write bands with extra digits precision
Ci         :             (has no effect for symmetry line mode)
Ci         :/spin1       generate bands only for first spin
Ci         :/lst=list    write only those bands specified in a list.
Ci         :             For syntax of list, see slatsm/mkilst.f
Ci         :/evn=#       keep track of smallest, largest eval for
Ci                       #th band, and printout at close.
Ci         :/ef=#        change efermi to #.
Ci         :Example: --band~long~qp~lst=2:5
Ci   iop   :options passed by the caller
Ci         :1s digit is suppress actions for parallel mode
Ci         : mode 1:
Ci         : 1: suppress writing line header info file
Ci         : 2: Return total number of qp to be generated in all lines
Ci         :    and also suppress writing line header info file
Ci         :    Thus iop=2 => all modes return total no qp to be generated
Ci         : mode 3:
Ci         : nonzero: do not allocate evsav until iq>nq
Ci
Ci   nband :(suqlst) maximum number of energy bands to write
Ci         :(suqlsw) actual number of energy bands to write
Ci
Ci   efermi:Fermi energy (written to bnds file)
Ci
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :NB used only in file write to indicate how many times
Ci         :nsp is ALTERED to nsp=1 if spin1 option is set
Ci
Ci   evl   :eigenvalues for current qp (used only in suqlsw)
Cio Inputs/Outputs (see Remarks)
Cio  nq    :Input nq=0 :if zero, flags suqlst to set up bands mode:
Cio                    :some initializations; then returns ready
Cio                    :for first block of qp
Cio        :           :if nonzero, setup for for next qp
Cio        :Output nq  :if zero, suqlst has no more qp to calculate.
Cio                    :if nonzero, nq is the number of k-points in the
Cio                    :current block.  Caller then calls suqlst nq
Cio                    :times (one for each k-point), which will return
Cio                    :qp for each k-point.  After exactly nq calls,
Cio                    :suqlst will start another block, if it exists.
Cio                    :See Remarks for schematic of calling sequence.
Co Outputs
Co   nfbn  :(color weights) number of elements in iblst
Co   iblst :(color weights) list of orbital indices for color weights:
Co         :color weight is the sum of contributions from eigenvector
Co         :components in iblst
Co   qp    :not set, if suqlst is called starting a new block.
Co         :(see Remarks).  Else, k-point at which to generate bands
Co   onesp :if spin1 flag is encountered, onesp is set to 1 and nsp is set to 1
Cl Local variables
Cl   iq    :current qp in this block
Cl   mode  :1 symmetry-line mode
Cl         :2 list-of-qp mode
Cl         :3 contour mode
Cl   q1    :starting qp for symmetry mode; only meaningful in that mode
Cl   q2    :ending   qp for symmetry mode; only meaningful in that mode
Cl   nevn  :(optional) band index; routine monitors largest, smallest value
Cl         :for that index.
Cl   evmnn :smallest value found for qp of specified band index
Cl   evmxn :largest  value found for qp of specified band index
Cl   ifiq  :file logical unit for input qp file
Cl   ifib  :file logical unit for output bands file
Cl   ql    :local copy of current qp
Cl   oqp
Cr Remarks
Cr   suqlst is designed to be called to generate qp in groups or
Cr   blocks.  The calling sequence is:
Cr
Cr     nq = 0  <- flags that first call to suqlst, to set up mode
Cr     do  iblock = 1, forever
Cr       This call generates nq, the number of points in this block
Cr       call suqlst(nband,efermi,nsp,evl,nq,qp,onesp)
Cr       if (nq .eq. 0) stop
Cr       do  iq = 1, nq
Cr         This call generates qp for current block
Cr         call suqlst(ndimh,ef0,nsp,w,nkp,qp,onesp) <- returns qp
Cr         do  isp = 1, nsp
Cr         call suqlsw(ndimh,qp,evl(1,isp)) <- saves evl for this qp
Cr                                             (call is optional)
Cr         enddo
Cr       enddo
Cr     enddo
Cr
Cr   The following modes are implemented:
Cr     mode=1 reads qp from syml file, and generates qp along each
Cr            specified symmetry line.  Structure of qp file:
Cr            file has one line for each symmetry line as follows:
Cr               nq      q1x   q1y   q1z      q2x   q2y   q2z
Cr               ...
Cr            the line entries have meanings:
Cr            --# qp-   ---starting qp---    --- ending qp ---
Cr            Any line with nq=0 implies no more lines.
Cr     mode=2 reads qp from specified file and generates qp for each
Cr            specified qp.  File consists sets of qpx,qpy,qpz for each
Cr            qp sought.  Typically a file would consists of lines like
Cr               q1x   q1y   q1z
Cr               q2x   q2y   q2z
Cr               ...
Cr     mode=3 generates qp for a uniform mesh in a plane (contour plot)
Cr            The file supplies information describing a rectangle in
Cr            the Brillouin zone.  It consists of a single line,
Cr            which contains the following:
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr
Cr            v1 and v2 are two vectors specifying the plane of the
Cr            contour.  range and n (one each for v1 and v2) are the
Cr            starting and final amplitudes of those vectors, and the
Cr            the number of points within the vector.  list-of-bands
Cr            is a list of integers which specify which bands are to
Cr            be written to the output file.  'height' is the 'z' axis.
Cr            For example,
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr            1 0 0  -1 1   51   0 1 0  -1 1   51   0.00    4,5
Cr            creates a file of 51x51 points, with the four corners
Cr            (-1,-1,0),  (1,-1,0),  (-1,1,0),  (1,1,0)
Cr            going through the gamma-point. Two bands (4,5) are stored.
Cu Updates
Cu   08 Jul 08 Extend to case where number of bands can be q dependent
Cu             modes 1,2: suqlsw writes out number of bands with qp
Cu   09 Jul 07 configured to with with MPIK mode
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1
Cu   02 Jul 06 Color mode extended to two colors
Cu   02 Jul 06 New color mode (one color weight only)
Cu   14 Feb 05 contour mode saves both spins in spin-polarized case
Cu   20 Oct 03 suqlst works properly in contour mode
Cu   28 Aug 01 added ef switch
Cu   23 Jan 01 first written
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) sopts
      integer nq,nband,nsp,nfbn(2),ifblst(nband,*),onesp,iop
      double precision efermi,qp(3),evl(nband)
C ... Local variables
      character strn*120, strn2*120, dc*1, fn*120, prfmt*40
      logical rdstrn,fmlong
      integer fopno,i,iprint,j,j1,j2,k,nglob,fopnn,rdm,
     .parg,a2vec,iv(7),stdo,jsp,op1
      double precision xx,xx1,xx2,sclp,plndst,xxv(7)
      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
      integer oqp,nqall,imode
      double precision x1,x2,y1,y2,evmxn,evmnn,
     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
      real(8),allocatable :: evsav(:,:,:,:)
C ... MPI
      integer procid,master,mpipid

ctakao
      real(8)::qlx(3)
C ... Heap
c      integer w(1)
c      common /w/ w
      common /suqlsd/
     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
      save evsav
      real(8),allocatable:: w_oqp(:,:)

      procid = mpipid(1)
      master = 0
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      op1 = mod(iop,10)
      nqall = 0

C --- First call ... setup and parse options ---
      if (nq .eq. 0) then
C   ... Defaults
        fmlong = .false.
        mode = 1
        fn = 'qp'
        nblst = 0
        nevn = 0
        evmxn = -99d9
        evmnn = 99d9
        nfbn(1) = 0
        nfbn(2) = 0

        dc = sopts(1:1)
        if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
          j2 = 0
   10     continue
          j2 = j2+1
          if (sopts(j2:j2) .eq. dc) goto 10
          j1 = min(len(sopts),j2)
          call nwordg(sopts,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (sopts(j1:j1+2) .eq. 'fn=')  then
              if (j1+3 .gt. j2) call rx('suqlst: bad file name')
              fn = sopts(j1+3:j2)
            elseif (sopts(j1:j2) .eq. 'qp')  then
              mode = 2
            elseif (sopts(j1:j2) .eq. 'spin1')  then
              onesp = 1
              nsp = 1
            elseif (sopts(j1:j2) .eq. 'con')  then
              mode = 3
            elseif (sopts(j1:j2) .eq. 'long')  then
              fmlong = .true.
            elseif (sopts(j1:j1+3) .eq. 'col=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, col=..')
              call mkils0(sopts(j1+4:j2),nfbn,ifblst)
              call mkilst(sopts(j1+4:j2),nfbn,ifblst)
              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
            elseif (sopts(j1:j1+4) .eq. 'col2=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, col2=..')
              call mkils0(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
              call mkilst(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
            elseif (sopts(j1:j1+3) .eq. 'lst=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, lst=..')
              call mkils0(sopts(j1+4:j2),nblst,iblst)
              if (nblst .gt. 100) call rx('increase size of iblst')
              call mkilst(sopts(j1+4:j2),nblst,iblst)
            elseif (sopts(j1:j1+2) .eq. 'ef=')  then
              j = 0
              i = parg('ef=',4,sopts(j1:),j,len(sopts(j1:)),
     .        dc//' ',1,1,i,efermi)
              if (i .le. 0) call rxs(
     .        'suqlst: failed to parse string for fermi level:  ',
     .        sopts(j1:))
            elseif (sopts(j1:j1+3) .eq. 'evn=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list')
              i = j1+3
              xxv(1) = 0
              j = a2vec(sopts,j2,i,2,dc//' ',2,3,1,iv,nevn)
              if (j .ne. 1 .or. nevn .gt. nband)
     .        call rx('suqlst: bad value for evn')
            else
              call rxs('suqlst: failed to parse argument, ',sopts(j1:j2))
            endif
            goto 10
          endif
        endif

        if (nfbn(1) .eq. 0) then
          call info2(20,0,0,' suqlst:  generate bands, mode %i',mode,0)
        elseif (nfbn(2) .eq. 0) then
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,
     .    ' suqlst:  generate bands with color weights, mode %i.'//
     .    ' %N%10f%i components: '//strn//'%a',mode,nfbn)
        else
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,' suqlst:  generate bands '//
     .    'with two color weights, mode %i.',mode,0)
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,
     .    '%10f%i components, first color:  '//strn//'%a',nfbn,0)
          call ilst2a(ifblst(1,2),nfbn(2),strn)
          call info2(20,0,0,
     .    '%10f%i components, second color: '//strn//'%a',nfbn(2),0)
        endif
        if (nfbn(1) .eq. 0 .and. nfbn(2) .ne. 0) call info0(20,0,0,
     .  ' suqlst:  (warning) col2 specified, but col not')

C   ... Open qp file
        if (procid .eq. master) then
          ifiq = fopno(fn)
          rewind ifiq

C   ... open bands file
          ifib = fopnn('BNDS')
          rewind ifib

C   ... Write header
          if (mode .eq. 1) then
            i = nblst
            if (nblst .eq. 0) i = nband
            if (nfbn(1) .eq. 0) then
C           Use separate format statment to circumvent gfortran bug
              write(ifib,335) i,efermi,0
  335         format(i5,f10.5,i6)
            elseif (nfbn(2) .eq. 0) then
              call ilst2a(ifblst,nfbn,strn)
              call strip(strn,j1,j2)
              write(ifib,336) i,efermi,1,strn(j1:j2)
  336         format(i5,f10.5,i6:'  col= ',a:'  col2= ',a)
            else
              call ilst2a(ifblst,nfbn,strn)
              call strip(strn,j1,j2)
              call ilst2a(ifblst(1,2),nfbn(2),strn2)
              call strip(strn2,j,k)
              write(ifib,336) i,efermi,2,strn(j1:j2),strn2(j:k)
            endif
          endif
        endif

C   ... Other initializations
        iq = 0
        nq = -1
      endif

C --- Setup for a new block of k-points, depending on mode ---
C     This branch occurs on completion of the last qp of the current block
C     which is marked by iq>nq
C     At the completion of this block: nq must be computed and:
C     (mode=1) q1,q2 set up.
C              Note: if 1s digit of iop is set in this mode,
C              this branch returns sum of all qp in all lines.
C              No setup for q1,q2; no
C     (mode=2) oqp allocated and loaded
C     (mode=3) nblst = number of bands to save
      if (iq .gt. nq) then
        iq = 1
C       Bands along specified symmetry lines
        if (mode .eq. 1) then
  725     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
          if (strn(1:1) .eq. '#') goto 725
C         Exit if first entry in line is zero
          i = 0
          xxv(1) = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,1,iv,xxv)
          if (i .eq. 1 .and. xxv(1) .eq. 0) goto 999
          i = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,7,iv,xxv)
          if (i .ne. 7 .and. iprint().ge.10) then
            write(stdo,
     .      '(/'' suqlst (warning) skipping line:''/''  '',a)') strn
            goto 725
          endif
          nq = xxv(1)
C         1 qp is nonsensical for a line
          if (nq .eq. 1) nq = 0
C         No qp: exit
          if (nq .le. 0) goto 999
C         setup q1,q2
          call dcopy(3,xxv(2),1,q1,1)
          call dcopy(3,xxv(5),1,q2,1)
          if (iprint().ge.10) write(stdo,785) nq,q1,q2
  785     format(/' suqlst:  nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
C         Write line header information to disk
          if (op1 .eq. 0) then
            write(ifib,337) nq*nsp
  337       format(2i5)
C         Accumulate all qp and cycle until all lines are read
          elseif (op1 .eq. 2) then
            nqall = nqall + nq
            goto 725
          endif
C       Bands for a list of specified qp
        elseif (mode .eq. 2) then
C         Only one block for this mode.
C         Flag that prior block already completed: nq>0
          if (nq .gt. 0) goto 999
          nq = 0
          j = 3
          if (rdm(ifiq,10000,0,' ',xx,j,nq) .ne. 1) call
     .    rx('suqlst:  bad qp file')

c          call defrr(oqp,3*nq)
          allocate(w_oqp(3,nq))

          rewind ifiq
          j = rdm(ifiq,10000,3*nq,' ',w_oqp,j,nq)
          call awrit1('%N suqlst: read %i qp from file '//fn//
     .    '%a',' ',80,stdo,nq)
          if (j .ne. 1) call rx('suqlst: failed to read qp')
          i = nblst
          if (nblst .eq. 0) i = nband
          if (nfbn(1) .ne. 0) i = i*2
          call awrit2('%% rows %i cols %i',' ',80,ifib,nq,i+4)
C       Bands on a uniform mesh in a specified plane (contour plot)
        elseif (mode .eq. 3 .and. nq .eq. -1) then
  825     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 998
          if (strn(1:1) .eq. '#') goto 825
          call words(strn,i)
          if (i .ne. 14) call rxi('suqlst con mode: expected 14 '//
     .    'arguments from input file but read',i)
          iblst(1) = -1
          nblst = 1
          backspace ifiq
          read(ifiq,*) vectx,x1,x2,nqx,vecty,y1,y2,nqy,plndst
          call word(strn,14,j1,j2)
          call mkilss(11,strn(j1:j2),nblst,iblst)
          if (nblst .le. 0) call rx('suqlst: no bands in list')
          call dscal(3,1/dsqrt(sclp(vectx,vectx)),vectx,1)
          call dscal(3,1/dsqrt(sclp(vecty,vecty)),vecty,1)
C         Subtract from vecty projection onto vectx
C         call daxpy(3,-sclp(vectx,vecty),vectx,1,vecty,1)
          call cross(vectx,vecty,vect3)
          call dscal(3,plndst/dsqrt(sclp(vect3,vect3)),vect3,1)
          nq = nqx*nqy
          if (iprint() .ge. 10) then
            write(stdo,717) vectx,x1,x2,nqx,plndst,vecty,y1,y2,nqy,nq
  717       format(' vx=',3f9.6,'  x1,x2=',2f9.6,'  nx=',i3,'   h=',f9.6/
     .      ' vy=',3f9.6,'  y1,y2=',2f9.6,'  ny=',i3,'  np=',i5)
            write(strn,'('' save %i bands: %'',i2,'':1i'')') nblst
            if (strn(17:18) .eq. '% ') strn(17:18) = ' %'
            call awrit2(strn,strn,80,stdo,nblst,iblst)
            if (dabs(sclp(vectx,vecty)) .gt. 1d-6 .and. iprint().gt.20)
     .      print *,'suqlst (warning): input plane vectors not orthogonal'
          endif

          if (iop .eq. 0) then
            allocate(evsav(nqx,nqy,nblst,nsp))
          endif

C       Contour plot, cleanup.
C       Note: cleanup handled by suqlsw when last qp is called.
        elseif (mode .eq. 3) then
C          if (op1 .ne. 0) then
C            allocate(evsav(nqx,nqy,nblst,nsp))
C            return
C          endif
          call rx('suqlst: caller should never reach this branch')
        endif

C --- Generate qp for this iq, depending on mode ---
      else
        if (mode .eq. 1) then
          xx = dble(iq-1)/dble(nq-1)
          qp(1) = xx*q2(1) + (1-xx)*q1(1)
          qp(2) = xx*q2(2) + (1-xx)*q1(2)
          qp(3) = xx*q2(3) + (1-xx)*q1(3)
        elseif (mode .eq. 2) then
          call dpscop(w_oqp,qp,3,iq*3-2,1,1d0)
        elseif (mode .eq. 3) then
C         Inner Loop:  excursions in y; outer loop: excursions in x
          j = mod(iq-1,nqy)
          i = (iq-1-j)/nqy
          if (nqx .le. 1) then
            xx1 = x1
          else
            xx1 =i*(x2-x1)/(nqx-1)+x1
          endif
          if (nqy .le. 1) then
            xx2 = y1
          else
            xx2 =j*(y2-y1)/(nqy-1)+y1
          endif
          do  93  k = 1, 3
            qp(k) = xx1*vectx(k) + xx2*vecty(k) + vect3(k)
   93     continue
          if (j .eq. 0 .and. iprint().ge.20) write(stdo,718) i+1,nqx,qp
  718     format(' line',i3,' of',i3,'   q(1)=',3f10.6)
        else
          call rx('suqlst: bad mode')
        endif
        iq = iq+1
C       Hold onto local copy of qp
        call dcopy(3,qp,1,ql,1)
      endif
      if(allocated(w_oqp)) deallocate(w_oqp)
      return

C --- No more qp blocks: cleanup ---
  999 continue
      nq = 0
      if (nevn .ne. 0) call awrit3(' eval no. %i:  minimum eval'//
     .' = %;8F  maximum eval = %;8F',' ',80,stdo,nevn,evmnn,evmxn)

      if (mode .eq. 1) then
        if (op1 .eq. 0) then
          write(ifib,337) 0
        elseif (op1 .eq. 2) then
          nq = nqall
        endif
      elseif (mode .eq. 3) then
        call rx('not ready')
C       call xxxbnd(w(oev),nblst,nqx,nqy,ifib)
C       return
      endif
      if(allocated(w_oqp)) deallocate(w_oqp)
      return

C --- Error exit ---
  998 call rxs('suqlst: failed to read file contents, file ',fn)

      entry suqlsm(imode)
C- Return qlist mode
      imode = mode
      return

      entry suqlsw(nband,jsp,nsp,evl)
C- Write or store the energy bands to file for this qp

      if (mode .eq. 1) then
        if (nblst .eq. 0) then
          prfmt = '(3f10.5,i6/(10f8.4))'
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.5/(10f8.4))'
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        else
          call rx('bug in suqlsw')
        endif

      elseif (mode .eq. 2) then
        if (fmlong .and. nblst .eq. 0) then
          prfmt = '(3f15.10,i6/(5f15.10))'
        elseif (nblst .eq. 0) then
          prfmt = '(3f10.6,i6/(8f10.6))'
        elseif (fmlong .and. nblst .gt. 0) then
          prfmt = '(3f15.10/(5f15.10))'
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.6/(8f10.6))'
C         if (nblst .le. 5) prfmt = '(8f10.6)'
        else
          call rx('bug in suqlsw')
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo
            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

C ... Keep running tab on smallest,largest eval
      if (nevn .ne. 0) then
        evmxn = max(evmxn,evl(nevn))
        evmnn = min(evmnn,evl(nevn))
      endif


      entry suqlsw2(nband,jsp,nsp,evl,qlx) !takao added ql
C- Write or store the energy bands to file for this qp
      ql=qlx
      if (mode .eq. 1) then
        if (nblst .eq. 0) then
          prfmt = '(3f10.5,i6/(10f8.4))'
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.5/(10f8.4))'
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        else
          call rx('bug in suqlsw')
        endif

      elseif (mode .eq. 2) then
        if (fmlong .and. nblst .eq. 0) then
          prfmt = '(3f15.10,i6/(5f15.10))'
        elseif (nblst .eq. 0) then
          prfmt = '(3f10.6,i6/(8f10.6))'
        elseif (fmlong .and. nblst .gt. 0) then
          prfmt = '(3f15.10/(5f15.10))'
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.6/(8f10.6))'
C         if (nblst .le. 5) prfmt = '(8f10.6)'
        else
          call rx('bug in suqlsw')
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo
            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

C ... Keep running tab on smallest,largest eval
      if (nevn .ne. 0) then
        evmxn = max(evmxn,evl(nevn))
        evmnn = min(evmnn,evl(nevn))
      endif
      end

      subroutine suqlse(nband,jsp,nsp,ndimhx,ifbn,nfbn,ifblst,ndlst,
     .evec,wk)
      use m_globalvariables
C- Write to file the projection of eigenvector subblock for this qp
C ----------------------------------------------------------------------
Ci Inputs
Ci   nband :number of energy bands to write
Ci   jsp   :current spin index (not used now)
Ci   nsp   :number of spins (not used now)
Ci   ndimhx:dimensions evec
Ci   ifbn  :index to color list (nfbn,ifblst)
Ci   nfbn  :number of elements for color weights projection
Ci   ifblst:list of elements for  color weights projection
Ci   ldlst :leading dimension of ifblst
Ci   evec  :eigenvectors
Ci   wk    :work array of same dimension as evl
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr   Inefficient, but it works
Cu Updates
Cu   08 Jul 08 New argument ndlst so ifblst can be dimensioned
Cu             independently from nband
Cu   05 Jun 06 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jsp,nsp,nband,ifbn,nfbn(2),ndlst,ifblst(ndlst,*)
      integer ndimhx
      double precision wk(1)
      double complex evec(ndimhx,ndimhx)
C ... Local parameters
      double complex zdotc
      character prfmt*40
      logical fmlong
      integer i,iprint,j,j1,k,nglob,stdo
      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c      integer oqp
      double precision x1,x2,y1,y2,evmxn,evmnn,
     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)

      common /suqlsd/
     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy

      integer ipiv(ndimhx)
      real(8),allocatable :: evsav(:,:,:,:)
      complex(8),allocatable:: evecc(:,:),work(:,:)

      allocate(evecc(ndimhx,ndimhx),work(ndimhx,ndimhx))

#if MULLIKEN
      call zcopy(ndimhx**2,evec,1,evecc,1)
      call zgetrf(ndimhx,ndimhx,evecc,ndimhx,ipiv,j)
      if (j .ne. 0) call rx('mullmf: failed to generate overlap')
      call zgetri(ndimhx,evecc,ndimhx,ipiv,work,ndimhx**2,j)
#endif
      do  i = 1, ndimhx
        wk(i) = 0
        do  j = 1, nfbn(ifbn)
          k = ifblst(j,ifbn)
          if (k .le. 0 .or. k .gt. ndimhx) then
            call fexit2(-1,111,' Exit -1 : suqlst: component %i'//
     .      ' outside range (1:%i)',k,ndimhx)
          endif
#if MULLIKEN
          wk(i) = wk(i) + evecc(i,k)*evec(k,i)
#else
          wk(i) = wk(i) + dconjg(evec(k,i))*evec(k,i)
#endif
        enddo

#ifndef MULLIKEN
        wk(i) = wk(i) / zdotc(ndimhx,evec(1,i),1,evec(1,i),1)
        if (wk(i) .lt. 0 .or. wk(i) .gt. 1+1d-14) then
          call rx('bug in suqlse')
        endif
#endif

      enddo

      if (mode .eq. 1) then
        prfmt = '(3f10.5/(10f8.4))'
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, (wk(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (wk(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 2) then
        if (fmlong) then
          prfmt = '(5f15.10)'
        else
          prfmt = '(8f10.6)'
C         if (nblst .gt. 0 .and. nblst .le. 5) prfmt = '(8f10.6)'
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) (wk(i),i=1,nband)
        else
          write(ifib,prfmt) (wk(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        call rx('need copy weights into different place, mode=3')
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = wk(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo

            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

      deallocate(evecc,work)

      end

      subroutine suqlsr(mode,ifi,nsp,nbf,lde,nb1,nb2,ls1,ls2,nq,qp,eb)
C- Count qp, optionally read them and energy bands from file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :Any combination of the following is allowed:
Ci         :1 count number of qp in file, return as nq
Ci         :  if 1s bit mode is 0, nq is an input,
Ci         :  the file nq cannot exceed input nq
Ci         :2 return qp read in file into array qp
Ci         :4 return bands read in file into array eb
Ci   ifi   :read from logical unit ifi
Ci   nsp   :number of spins in band file
Ci   nbf   :number of bands in band file
Ci   lde   :leading dimension of eb
Ci   nb1,nb2: read bands nb1..nb2 into eb
Ci   ls1   :(nsp=2): read first spin only
Ci   ls2   :(nsp=2): read second spin only
Cio Inputs/Outputs
Cio  nq    :numbefr of k-points (returned if 1s bit of mode is set)
Co Outputs
Co   qp    :k-points, returned if 4s bit of mode is set
Co   eb    :energy bands, returned if 4s bit of mode is set
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 Jun 07  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ifi,nsp,nq,nbf,lde,nb1,nb2
      logical ls1,ls2
      double precision qp(3,nq),eb(lde,nsp,nq)
C ... Local parameters
      integer nq1,nql,mode1,mode2,mode4,iq1,i,nline,isp
      double precision ql(3),evl(3000)

      if (mode .eq. 0) return
      mode1 = mod(mode,2)
      mode2 = mod(mode/2,2)
      mode4 = mod(mode/4,2)
      if (mode4 .ne. 0) then
        call isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
     .  .true.)
        call isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
      endif

C --- For each panel, do ---
      nql = 0
      nline = 0
   91 continue
      read(ifi,*) nq1
      if (nq1 .le. 0) goto 90
      isp = 0
      do  iq1 = 1, nq1
        isp = mod(isp,2)+1
        nql = nql+1
        if (mode1 .eq. 0 .and. nql .gt. nq) call rxi(
     .  'suqlsr: file more q-points than allocated: nqmx=',nq)
        read(ifi,*,END=999,ERR=999) ql(1),ql(2),ql(3)
        if (mode2 .ne. 0) call dcopy(3,ql,1,qp(1,nql),1)
        read(ifi,*,END=999,ERR=999) (evl(i),i=1,nbf)
        if (mode4 .ne. 0) then
C         Copy only if appropriate spin
          if (isp .eq. 1 .and. ls2) then
          elseif (isp .eq. 2 .and. ls1) then
          else
            call dcopy(nb2-nb1+1,evl(nb1),1,eb(1,1,nql),1)
          endif
        endif
      enddo
      nline = nline+1
      goto 91
C     End of loop over lines
   90 continue


      if (mode1 .eq. 0) then
        call info2(30,1,0,' suqlsr: found %i qp in %i lines from file',
     .  nql,nline)
      else
        nq = nql/nsp
        call info2(30,1,0,' suqlsr: read %i qp in %i lines from file',
     .  nql,nline)
      endif

      if (mod(nql,2) .ne. 0 .and. nsp .eq. 2) then
        call info0(20,0,0, ' suqlsr (warning): '//
     .  'odd number of bands encountered but nsp=2')
      endif

      return
  999 continue
      call rxi('suqlsr: failed to read bands file, nq=',nql)

      end

