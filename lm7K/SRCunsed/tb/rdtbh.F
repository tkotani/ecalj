      subroutine rdtbh(ifi,nclass,dclabl,z,nlmesp,ltb,nterm,ipass,
     .  fitpar,iam,npm,tabme,tabcf,tabov,tbocf,decay,deccf,decov,dcocf,
     .  itab,itbcf,itbov,itocf,idec,idcf,idov,idocf,V0,k,memode)
C- Disk-read of tight-binding hamiltonian and pair potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi: file handle for CTRL
Ci   nlmesp: (number of ME/spin channel) * (# spin channels)
Ci           for which ME are input
Ci   nclass,nsp,ltb
Ci   nterm: number of parameters for each matrix element
Ci   ipass: 1, read tabme; 2, generate iam, V0
Ci   fitpar: if true then fitting parameters, read freeze switches
Co Outputs
Co   iam(1..3,kk) are a group of three numbers that associate a pair
Co          with a set of matrix elements.  For a pair of classes 
Co          iam(1,kk) and iam(2,kk), iam(3,kk) is an index specifying which
Co          rule in tabme (and possibly tabcf,tabov,tbocf) that defines the
Co          matrix elements between these classes.  It is possible that
Co          there is no rule connecting two particular classes, in which
Co          case the matrix element is taken to be zero.
Co          iam is ordered according to increasing class.
Co   npm    is a table of offsets to iam so that every rule associated
Co          with a given class can be easily identified.  npm(0,i) is
Co          the number of entries in iam containing class i, and
Co          npm(1,i) is the offset to iam for the first entry.
Co   tabme: a set of parameters that correspond to coefficients of
Co          Slater-Koster, or Slater-Koster-like matrix elements.
Co          The meaning of the coefficients, depends on mode;
Co          see Remarks, and memode, below.  Which pair of atoms
Co          the coefficients apply is defined by a set of rules,
Co          as described in Remarks.
Co          tabme is a list of parameter, for as many rules
Co          are are read in.  
Co          When ipass=1, the number of rules is returned (in k)
Co          and the corresponding coeficients are read in tabme
Co   tabcf: table of crystal field MEs [see D.J. Chadi in ``Atomistic
Co          Simulation of Materials Beyond Pair Potentials'', edited
Co          by V. Vitek and D. Srolovitz (Plenum, 1989), page 309].
Co          tabcf has the same structure as tabme
Co   tabov: table of overlap matrix elements, structured as tabme
Co   tbocf: table of crystal field MEs for overlap, structured as tabme
Co   decay,deccf,decov,dcocf: exponential or power decay parameters
Co      matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Co   itab,itbcf,itbov,itocf,idec,idcf,idov,idocf:  switches for all the
Co      matrix element and decay parameters: if 0 vary param., if 1 fix;
Co      only read if fitpar=.true.
Co   V0: parameters for pair potential
Co       Ordering: a_1 b_1 c_1 a_2 ... c_3
Co       V0 = \sum_i=1,3 a_i d^b_i exp(-c_i d) unless b_1 > 0 and then
Co       V0 = a_1 eps + a_2 eps^2 where eps = (d - c_1) / c_1
Co       If b_1 > 0 and c_1 < 0 then Goodwin-Skinner-Pettifor
Co       V0 = A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)]
Co            ordering: A (b_1 > 0) (c_1 < 0) m mc r0 rc
Co   k (pass 1) number of ME rules (pass 2) number of ME pairs
Co   memode: 0, fixed MEs
Co           1  Harrison universal MEs
Co           2, exponential decay
Co           3, power decay
Co           4, ME = \sum_i=1,3 a_i d^b_i exp(-c_i d), the ordering is:
Co              a_1 b_1 c_1 a_2 ... c_3 for ss-sigma, then sp-sigma, etc
Co           5, Goodwin-Skinner-Pettifor,
Co              v_ij (r0/d)^n exp[n (-{d/rc}^nc + {r0/rc}^nc)]
Co              ordering: v_ij n nc r0 rc for ss-sigma, etc
Co           6, (nl=1) Slater-Koster + Srini's extra term
Co              NB: NOT implemented for spin pol.
Co           7, a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc)
Co              ordering: a b c d0 for ss-sigma, etc
Co           nterm=1 for memode=0-3, nterm=9 for memode=4,
Co           nterm=5 for memode=5, nterm=2 for memode=6,
Co           nterm=4 for memode=7
Co           memode >= 10, use canonical TB Hamiltonian
Cr Remarks
Cr *rdtbh parses category ME, looking for a set of rules defining the
Cr  matrix elements.  For each rule, rdtbh reads a vector of numbers;
Cr  the precise meaning of these numbers depends on memode.  But
Cr  typically, the numbers are coefficients to Slater-Koster matrix
Cr  elements, whose ordering is:
Cr     ss-sigma, sp-sigma, pp-sigma, pp-pi, sd-sigma, pd-sigma,
Cr     pd-pi, dd-sigma, dd-pi, dd-delta
Cr 
Cr *For each rule, rdtbh reads in two lists of classes (precisely how
Cr  the lists are defined is described below), followed by a vector
Cr  of numbers that define the matrix elements (eg ss-sigma, etc), then
Cr  optionally followed by another vector of numbers for the overlap,
Cr  crystal-field terms, etc, and finally another vector of numbers
Cr  for the pairwise potential.
Cr
Cr *A rule takes the following schematic form:
Cr     i-class-list j-class-list | me ! V0
Cr  If crystal field terms and/or overlap the form is
Cr     i-class-list j-class-list | tb me & cf me @ ovl me % ovl cf ! V0
Cr                          or
Cr     i-class-list j-class-list | tb me @ ovl me & cf me ! V0
Cr                          or etc.
Cr
Cr  where me (or tb me, cf me, ovl me, V0) are a vector of numbers.  How
Cr  many numbers there are depends on nl and memode, except for the
Cr  pair potential vector, which is always of length 9.
Cr   
Cr *In the spin-orbit case all the (+,+) MEs for a given atom pair are
Cr  listed first, followed by the (-,-) MEs and then the (+,-) MEs
Cr
Cr *Syntax of class-list (i- and j-class-list): This string specifies
Cr  which groups of classes are to be associated with this rule.
Cr  It can be a simple number, eg for memode=0 and nl=2, a rule can be:
Cr     1 2 | vsss vsps vpps vppp
Cr  More generally,  the i- and j- class-list can be a list of integers
Cr  such as 1,2,3 or 1:5; see mkilst for the syntax of an integer list.
Cr
Cr  It is also possible to define the list in one of two other alternate
Cr  styles.   The second style is to define the list according to
Cr  an expression, that is i- or j-class-list is some expression 
Cr  the variables ic and z; thus class-list may look like ic<6&z==14.
Cr  Every class for which the expression evaluates to true belongs in
Cr  the list.  Thus the rule
Cr     ic<6 ic>=6 | vsss vsps vpps vppp
Cr  defines the rule that will connect every class up to the fifth
Cr  class with every class beyond the fifth class.
Cr
Cr  The third alternative is specifically for unix systems.  i- and j-
Cr  class-list are filenames with the usual unix wildcards, eg
Cr  a[1-6].  Each class has an associated file; rdtbh invokes a 
Cr  shell 'ls class-list | grep classname'.  Any class which ls finds
Cr  is added to the list.  The following example
Cr     a[1-6] b* | vsss vsps vpps vppp
Cr  creates a rule connecting any class a1, a2, a3, ...  a6 to any class
Cr  whose name begins with 'b'
Cr 
Cr *You tell rdtbh to use one of these styles at the start of
Cr  the ME category, directly after memode, eg
Cr    ME  3 CLSTYL=2
Cb Bugs
Cb    In some cases if the pair potential line (!) is missing then rdtbh
Cb    fails. Specifically, at line 296 iend is not properly shifted in
Cb    which case at line 337 the overlap parameters are read again
Cb    into a rule k+1 for which there is no space to store them. 
Cb    The workaround is to always include a pair potential line even if
Cb    it is a series of zeros.
Cu Updates
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
Cu   21 Jul 02 (MvS) Bug fix in partok call
Cu   15 Feb 02 (ATP) Added MPI parallelization
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
C#ifdefC MPI
C      implicit none
C      include "mpif.h"
C      integer numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C#endif
      integer ifi,nclass,nlmesp,nterm,ipass,k,memode,ltb
      integer procid, master, mpipid
      double precision dclabl(nclass)
      integer iam(3,2),npm(0:1,1),itab(nterm*nlmesp,2),
     .  itbcf(nterm*nlmesp,2),itbov(nterm*nlmesp,2),
     .  itocf(nterm*nlmesp,2),idec(nlmesp,2),idcf(nlmesp,2),
     .  idov(nlmesp,2),idocf(nlmesp,2)
      double precision tabme(nterm*nlmesp,*),tabcf(nterm*nlmesp,2),
     .  tabov(nterm*nlmesp,2),tbocf(nterm*nlmesp,2),decay(nlmesp,2),
     .  deccf(nlmesp,2),decov(nlmesp,2),dcocf(nlmesp,2),V0(9,2),z(1)
      logical fitpar
C For iolib
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio
C heap:
      integer w(1)
      common /w/ w
C Local parameters
      integer memx,i,j,inext,owk,partok,iprint,is,js,is0,is1,ijmx,nt,
     .  memod0,ntm0,ival,ii,jj,jjj,kk,recln0,recln
      character*8 clabl
      parameter (ijmx=2000,recln0=120)
      integer ij(ijmx,2),nij(2),clstyl,iv0,ic,fextg
      double precision decay0
      logical scat,lsequ,cryf,ovl,ocryf,a2bin,sw
      character*1 ext*40
      character*24 mode(0:7)
      parameter (memx=500)
      character*(recln0) mfile(memx),strn
C ... for rdfiln
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,ctlen=120)
      character vnam(mxlev)*16,ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev),bittst
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
      data mode /'Fixed','Universal ME','Exp. decay','Power decay',
     .  'a d^b exp(-c d)', 'Goodwin-Skinner-Pettifor',
     .  'Fixed+extension','a d^-b/{1+exp[c(d-d0)]}'/

C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#endif

C --- Initialization ---
      i = recln(recln0)
      optio = 1
      recoff = 0
      decay0 = 0
      clstyl = 0
      ovl   = bittst(ltb,1)
      cryf  = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      procid = mpipid(1)
      master = 0

C --- Read in ME category from CTRL file ---
      if (procid .eq. master) then
      if (.not. scat(ifi,'ME ',' ',.true.))
     .  call rx('RDTBH: missing category ME')
      backspace ifi

      nrecs = 0
      do  10  i = 1, memx
        mfile(i) = 'z'
        call rdfiln(ifi,'#{}% c',mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,strn,mfile(i),recln0,nrecs)
*       if (.not. rdstrn(ifi,mfile(i),recln0,.false.)) goto 20
        if (i .gt. 1 .and. mfile(i)(1:1) .ne. ' ') goto 20
*       print *, mfile(i)
   10 continue
      call rx('RDTBH: increase memx')
   20 continue
      endif
C#ifdefC MPI
C      call MPI_BCAST(mfile,memx*recln0,MPI_CHARACTER,
C     .               master,MPI_COMM_WORLD,ierr)
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#endif
      call getcat(mfile,'ME ',' ',.true.)

C --- Determine ME mode, clstyl ---
      j = partok(mfile,'ME ',' ',memode,' ',-1,2,0,1)
C ... Pick up optional clstyl
      is = iend
      call skipbl(mfile,memx*recln0,is)
      js = is
      is = is+1
      call skp2bl(mfile,memx*recln0,js)
      if (mfile(1)(is:is+6) .eq. 'CLSTYL=') then
        j = partok(mfile,'CLSTYL=','=',clstyl,' ',-1,2,0,1)
      endif

      memod0 = mod(memode,10)
      if (iprint() .ge. 20 .and. ipass .eq. 2 .and. memod0 .ge. 0) print
     .  '(/'' RDTBH:  matrix elements with mode: '', a)', mode(memod0)

      ntm0 = 1
      if (memod0 .eq. 4) ntm0 = 9
      if (memod0 .eq. 5) ntm0 = 5
      if (memod0 .eq. 6) ntm0 = 2
      if (memod0 .eq. 7) ntm0 = 4
      if (ipass .eq. 2)
     .  call rxx(nterm .ne. ntm0,'RDTBH: dimension mismatch')

C --- Read in global exponential or power decay parameter ---
      inext = iend
      if (memod0 .eq. 2 .or. memod0 .eq. 3) then
        j = partok(mfile,'DECAY0=','=',decay0,' ',-1,4,0,0)
      else
        j = partok(mfile,'DECAY0=','=',decay0,' ',1,4,0,0)
      endif
      if (ipass .eq. 2) then
        call dcopy(k*nlmesp,decay0,0,decay,1)
        if (cryf)  call dcopy(k*nlmesp,decay0,0,deccf,1)
        if (ovl)   call dcopy(k*nlmesp,decay0,0,decov,1)
        if (ocryf) call dcopy(k*nlmesp,decay0,0,dcocf,1)
        if (fitpar) then
          call iinit(idec,k*nlmesp)
          if (cryf)  call iinit(idcf, k*nlmesp)
          if (ovl)   call iinit(idov, k*nlmesp)
          if (ocryf) call iinit(idocf,k*nlmesp)
        endif
      endif
      if (j .ne. 0) inext = iend
      k = 0
      kk = 0

C --- Loop over input lines ---
      do 60  i = 1, memx
        k = k+1
        is = inext
        jj = 0

C --- Read V0 parameters ---
        ii = 1
        if (ipass .eq. 2) ii = k
        jj = partok(mfile,'!','!',V0(1,ii),' ',-9,4,is,0)
        if (jj .ne. 0) inext = iend

        if (ipass .eq. 2) then
C --- Read crystal field ME and decay parameters ---
          if (cryf) then
            jjj = partok(mfile,'&','&',tabcf(1,k),' ',
     .        -ntm0*nlmesp,4,is,0)
            if (jjj .ne. 0) then
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(1)(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',deccf(1,k),' ',
     .          -nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZ=',4,'=',is0))
     .            jjj = partok(mfile,'FRZ=','=',itbcf(1,k),' ',
     .            -ntm0*nlmesp,2,is1,0)
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idcf(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif

C --- Read overlap ME and decay parameters ---
          if (ovl) then
            jjj = partok(mfile,'@','@',tabov(1,k),' ',
     .        -ntm0*nlmesp,4,is,0)
            if (jjj .ne. 0) then
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(1)(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',decov(1,k),' ',
     .          -nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZ=',4,'=',is0))
     .            jjj = partok(mfile,'FRZ=','=',itbov(1,k),' ',
     .            -ntm0*nlmesp,2,is1,0)
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idov(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif

C --- Read overlap crystal field ME and decay parameters ---
          if (ocryf) then
            jjj = partok(mfile,'%','%',tbocf(1,k),' ',
     .        -ntm0*nlmesp,4,is,0)
            if (jjj .ne. 0) then
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(1)(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',dcocf(1,k),' ',
     .          -nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZ=',4,'=',is0))
     .            jjj = partok(mfile,'FRZ=','=',itocf(1,k),' ',
     .            -ntm0*nlmesp,2,is1,0)
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(1)(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idocf(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif
        endif

C --- Read tight-binding ME and decay parameters ---
        j = ntm0*nlmesp
        if (memod0 .eq. 6) call rx('check branch')
C       if (memod0 .eq. 6) j = 8
        j = partok(mfile,'|','|',tabme(1,k),' ',-j,4,is,0)
        if (j .ne. 0) then
          nt = recln0*memx
          is1 = iend
          call skipbl(mfile,nt,is1)
          is0 = 1
          if (lsequ(mfile(1)(is1+1:is1+1),'DECAY=',6,'=',is0))
     .      jjj = partok(mfile,'DECAY=','=',decay(1,k),' ',
     .      -nlmesp,4,is1,0)
          if (ipass .eq. 2 .and. fitpar) then
            nt = recln0*memx
            is1 = iend
            call skipbl(mfile,nt,is1)
            is0 = 1
            if (lsequ(mfile(1)(is1+1:is1+1),'FRZ=',4,'=',is0))
     .        jjj = partok(mfile,'FRZ=','=',itab(1,k),' ',
     .        -ntm0*nlmesp,2,is1,0)
            nt = recln0*memx
            is1 = iend
            call skipbl(mfile,nt,is1)
            is0 = 1
            if (lsequ(mfile(1)(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .        jjj = partok(mfile,'FRZDEC=','=',idec(1,k),' ',
     .        -nlmesp,2,is1,0)
          endif
        endif
        if (j .ne. 0) inext = max(inext,iend)
        if (j .eq. 0 .and. jj .eq. 0) goto 70

        if (ipass .eq. 2) then
C --- Assemble {ij} ---
          do  40  ii = 1, 2
            call skipbl(mfile,memx*recln0,is)
            js = is
            call skp2bl(mfile,memx*recln0,js)
            if (clstyl .lt. 2) then
              call mkilst(mfile(1)(is+1:js+1),nij(ii),ij(1,ii))
            elseif (clstyl .eq. 2) then
              call numsyv(iv0)
              nij(ii) = 0
              do  42  ic = 1, nclass
                call lodsyv('ic',1,dble(ic),ival)
                call lodsyv('z',1,z(ic),ival)
                is1 = is
                if (a2bin(mfile,sw,0,0,' ',is1,js)) then
                  if (sw) then
                    nij(ii) = nij(ii)+1
                    ij(nij(ii),ii) = ic
                  endif
C           ... Abort if a2bin can't parse expression
                else
                  call rx('rdtbh: failed to parse'//mfile(1)(is+1:js+1))
                endif
   42         continue
              call clrsyv(iv0)
            elseif (clstyl .eq. 3) then
              nij(ii) = 0
              j = fextg(ext)
              do  44  ic = 1, nclass
                call r8tos8(dclabl(ic),clabl)
C#ifdefC GNU
C                call awrit0('%xls '// mfile(1)(is+1:js) // ext //'%a|'//
C     .            'grep -s ''\\b'//clabl//'%a'//ext//'%a''>/dev/null',
C     .            strn,len(strn),0)
C#else
                call awrit0('%xls '// mfile(1)(is+1:js) // ext //'%a|'//
     .            'grep -s '//clabl//'%a'//ext//'%a>/dev/null',
     .            strn,len(strn),0)
C#endif
                call locase(strn)
                call fsystm(strn,j)
                if (iprint().ge. 70) print '(''found='',i3,1x,a)',j,strn
                if (j .eq. 0) then
                  nij(ii) = nij(ii)+1
                  ij(nij(ii),ii) = ic
                endif
   44         continue
            else
              call rx('rdtbh: bad clstyl')
            endif
            is = js
   40     continue
          if (nij(1) .gt. ijmx .or. nij(2) .gt. ijmx)
     .      call rx('RDTBH: increase ijmx')

C --- Generate iam ---
          do  50  ii = 1, nij(1)
          do  50  jj = 1, nij(2)
            if (ij(ii,1) .gt. nclass .or. ij(ii,1) .lt. 1 .or.
     .          ij(jj,2) .gt. nclass .or. ij(jj,2) .lt. 1) then
              print *, ' rdtbh: no such pair ', ij(ii,1), ij(jj,2)
            else
              kk = kk+1
              iam(1,kk) = ij(ii,1)
              iam(2,kk) = ij(jj,2)
              iam(3,kk) = k
            endif
   50     continue
        endif

   60 continue
      call rx('RDTBH:  increase memx')

   70 continue
      k = k-1

C --- Sort iam, generate npm ---
      if (ipass .eq. 2) then
        k = kk
        call defi(owk,k)
        call ivshel(3,k,iam,w(owk),.false.)
        call rlse(owk)

        if (iprint() .gt. 40) then
          print '(/'' RDTBH:'', i5, '' matrix element pairs found:'')',k
          if (iprint() .ge. 50) print 333, ((iam(i,j), i=1,3), j=1,k)
  333     format(3i5)
        endif

        call iinit(npm,2*nclass)
        do  80  j = 1, k
          if (iam(1,j) .gt. nclass .or. iam(2,j) .gt. nclass .or.
     .        iam(1,j) .lt. 1    .or. iam(2,j) .lt. 1) then
            print *, ' j, iam(1,j), iam(2,j)=', j, iam(1,j), iam(2,j)
            call rx('RDTBH:  bad basis index in category ME')
          endif
          npm(0,iam(1,j)) = npm(0,iam(1,j))+1
   80   continue

        npm(1,1) = 0
        do  90  i = 2, nclass
   90   npm(1,i) = npm(1,i-1) + npm(0,i-1)

      endif
      return

      end

