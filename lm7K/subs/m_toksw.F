      module m_toksw
C- Module to create token switches, and read contents of tokens
C ----------------------------------------------------------------------
Cr First, a list of tokens must be created, one list each for
Cr any program which uses this module.
Cr Suppose a particular program 'foo' reads tokens
Cr   STRUC_ALAT STRUC_PLAT SITE_ATOM_POS HAM_TOL ...
Cr Further suppose that HAM_TOL is "optional", i.e. no error occurs
Cr if the HAM_TOL is not found
Cr Create the token list corresponding to 'foo' as follows:
Cr    call nwsadd()  ! starts list for a new program
Cr    call tkadd(" foo::" ) ! Identifies foo as progn associated w/ this list
Cr    call tkadd(" STRUC_PLAT") ! Append this token to the list
Cr    call tkadd(" HAM_TOL~")   ! ~ => token can be optionally read
Cr    call tkadd(" STRUC_PLAT SITE_ATOM_POS") ! lump several tokens in 1 call
Cr    etc
Cr
Cr Examples see how routine toksw_init in m_rdctrl.f bundles the above
Cr calls for a particular set of programs, e.g. "LMF", "LM", etc.
Cr
Cr After the lists have been created, function tksw(foo,tok) returns
Cr a value 0,1,2 indicating :
Cr   0 tok is sought by foo, with input optional
Cr   1 tok is sought by foo, with input required
Cr   2 tok is is not sought by foo
Cr
Cr The input file reader contains commands to read the contents of all
Cr tokens any program using this input system will seek.
Cr It calls tksw to determine whether a particular token is be read
Cr by the given calling program.
Cr
Cr This module also contains module m_gtv, the module that read and
Cr parses tokens and tokens' contents.  The entry points are:
Cr   gtv_r8,gtv_r8v,gtv_i4,gtv_i4v,gtv_lg,gtv_char,gtv_none
Cr and are usually invoked by the interface gtv, e.g. this call:
Cr     call gtv('HEADER',tksw('LMF','HEADER'),header,note=
Cr    .  'Contents displayed at beginning of program execution')
Cr reads string 'header'
C ----------------------------------------------------------------------
      implicit none
C ... Module variables
      logical,private:: debug=.false.
      integer,parameter,private :: lenmax=3000,nswmx=25
      character(lenmax),private:: swtok(nswmx)
      integer,private ::  nsw=0  !lenc=0,

      contains

      subroutine clear_swtok(debug_in)
C- Resets all tokens in swtok and sets internal debug switche
C ----------------------------------------------------------------------
Ci Inputs
Ci   debug_in : sets internal debug switch for future calls
Co Outputs
Cm Module variables
Cm  swtok : is initialized
Cm  nsw   : is initialized to zero
Cm  debug : is assigned to debug_in
Cr Remarks
Cr
Cu Updates
Cu   16 Jan 09
C ----------------------------------------------------------------------
C ... Passed parameters
      logical:: debug_in
C ... Module variables
      debug = debug_in
      swtok = ' '
      nsw = 0
      end subroutine

      subroutine nswadd()
C- Flag that all tokens have been appended for current main program
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Cm Module variables
Cm  nsw : nsw is incremented, starting a new
Cm      : entry in character array swtok
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none

      nsw = nsw+1
      if (nsw>nswmx) call rx('m_toksw: enlarge nswmx')
      end subroutine

      subroutine tkadd(tok)
C- Append tok to current set in swtok
C ----------------------------------------------------------------------
Ci Inputs
Ci   tok  : Append 'tok' to swtok(nsw)
Ci        : Note: nsw must be set previously by calling nswadd()
Co Outputs
Co  swtok : tok is appended to (module string) swtok(nsw)
Cl Local variables
Cr Remarks
Cr   This routine builds up character array swtok for a particular
Cr   main program.
Cr
Cr   swtok contains a set of lists of tokens, one list
Cr   for each main program, which are read by the main program,  e.g.
Cr     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
Cr     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
Cr   The list for each main program contains a sequence of strings,
Cr   separated by a space, each string corresponding to the full name
Cr   of a token (incl parents).  Optionally the name may be appended
Cr   by a "~" which signifies the token can be optionally read.
Cr
Cr   swtok is built up by the following sequence:
Cr     call nwsadd()  ! starts new list for another program
Cr     call tkadd(" prgn::" )  ! prgn = main program name, e.g. "LMF"
Cr     call tkadd(" tok")      ! tok = token name, e.g. SITE_ATOM_POS
Cr     call tkadd(" tok")      ! etcetera
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*):: tok
C ... Local parameters
      integer:: lensw

      if (nsw > nswmx) call rx('m_toksw: enlarge nswmx')
      lensw = len(trim(swtok(nsw)))
      if ( lensw+30>lenmax ) call rx('m_toksw: enlarge lenmax')
      swtok(nsw) = trim(swtok(nsw))//trim(tok)

C     lenc = max(lenc,lensw)  ! lenc (to check the length of buff)
C     print *, 'toksw for nsw=',nsw,' is now:', trim(swtok(nsw))
      end subroutine

      integer function tksw(prgn,name)
C- Return 0, 1 or 2 for given program and token name
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgn : main program name
Ci        : From a list of tokens 'prgn' uses (see Remarks),
Ci        : determine whether the supplied name is in the list.  Return:
Ci        : 0 if 'name'~ is present => input sought by 'prgn', but input
Ci            is optional
Ci        : 1 if 'name' is present => input is sought by 'prgn'
Ci        : 2 if 'name' is missing => input not sought by 'prgn'
Ci  name  : full name of token (incl parents) e.g. SITE_ATOM_POS
Cr Remarks
Cr   This routine searches character array swtok (which must be
Cr   generated in advance) to determine how a program is to treat an
Cr   input token.
Cr
Cr   swtok contains a set of lists of tokens, one list
Cr   for each main program, which are read by the main program,  e.g.
Cr     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
Cr     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
Cr   The list for each main program contains a sequence of strings,
Cr   separated by a space, each string corresponding to the full name
Cr   of a token (incl parents).  Optionally the name may be appended
Cr   by a "~" which signifies the token can be optionally read.
C------------------------------------------------------
      character*(*):: name,prgn
      integer:: isw,ix,ixe,sw,ix1,ix2,npgsz=30

C ... Find which of swtok corresponds to prgn
      do  isw = 1, nsw
        ix = index(swtok(isw)(1:npgsz),trim(prgn)//'::')
cki        if (ix/=0 .and. swtok(isw)(ix-1:ix-1) .eq. ' ') goto 100
        if (ix/=0) then
         if (swtok(isw)(ix-1:ix-1) .eq. ' ') goto 100
        endif
      enddo
      call rx('tksw: no switches for program '//prgn)
  100 continue

C      ix= index(swtok(isw),name)
C      print *,'tksw: name= '//'###'//name//'###',ix

C     See if found among required tokens
      ix1 = index(swtok(isw),' '//trim(name)//' ')
C     See if found among optional tokens
      ix2 = index(swtok(isw),' '//trim(name)//'~')

C ... Determine whether optional or required
      if (ix1 == 0 .and. ix2 == 0) then
        ix = 0
      elseif (ix1 == 0 .and. ix2 /= 0) then
        ix = ix2+1
      elseif (ix2 == 0 .and. ix1 /= 0) then
        ix = ix1+1
      else
        call rx('tksw, prgn '//trim(prgn)//': '//trim(name)//
     .  ' appears in both optional and required lists')
      endif
      ixe = ix+len(trim(name))
      if (ix == 0) then
        sw=2 !neither optional nor required: ignore
        if (debug) print *,'tksw: '//trim(name)//' ignored'
      elseif (swtok(isw)(ixe:ixe)=='~') then
        sw=0 !optional
        if (debug) print *,'tksw: '//trim(name)//' optional'
      else
        sw=1 !required
        if (debug) print *,'tksw: '//trim(name)//' required'
      endif
c      print *,'tksw: ', trim(prgn),' ',trim(name),sw
      tksw = sw
      end function
      end module


      subroutine find_region(lev,instr,token,toks,tokt,ntk,itrm,eor,
     .is,ie)
C- Finds substring corresponding to ntk-th token.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lev   :token level, reserved for acceleration (not used now).
Ci   instr :input string, where token and contents are embedded
Ci   token :token to match in instr; see Remarks
Ci   toks: :(only used if 1000s digit of itrm is set)
Ci         :list of characters that must precede token to satisfy a
Ci         :match.  Eliminates strings matching token embedded within
Ci         :other string.  Example:
Ci         :Example: given token ATOM, and
Ci         :toks='[ ' matches ' ATOM' and '[ATOM] but not 'CATOM'
Ci   tokt: :terminator(s) to token
Ci         :Example: given token ATOM,
Ci         :tokt='=:' matches matches ATOM= and ATOM:
Ci   ntk   :search for ntk-th occurence of token
Ci   itrm  :(itrm=0) start-of-region begins with second character following
Ci         :         token (1st char is token terminator)
Ci         :         end-of-region is delimited by first occurence of string
Ci         :         contained in eor following start-of-region
Ci         :
Ci         :(itrm=1) the first nonblank character after the token and its
Ci         :         terminator must be '['.  Note: '[' may also server as
Ci         :         start-of-region begins with the character following it.
Ci         :         end-of-region is delimited by the matching ']'
Ci         :         Note: '[...]' pairs may be nested; it is an error for
Ci         :         '[' not to have a matching ']'.
Ci         :
Ci         :(itrm=2) If the first nonblank character after the token and its
Ci         :         terminator is '[', follow the syntax of itrm=1.
Ci         :         Otherwise, follow the syntax of itrm=0.
Ci         :(itrm=11)Identical to itrm=1, with the addition that '[' may
Ci         :         serve both as terminator and start-of-region marker
Ci         :(itrm=12)Identical to itrm=2, with the addition that '[' may
Ci         :         serve both as terminator and start-of-region marker
Ci
Ci         :Adding 100 to itrm causes find_region to move
Ci         :start-of-region marker to 1st nonblank character
Ci         :Adding 1000 to itrm turns on the pre-token matching;
Ci         :see toks above
Ci
Ci   eor   :string demarcating end-of-region.  Its use depends on
Ci         :how start-of-region was determined (see itrm)
Ci         :If start-of-region is NOT determined by '[',
Ci         :the first string after start-of-region that matches the
Ci         :contents of eor demarcates end-of-region
Ci         :If start-of-region IS determined by '[', eor is not used
Co Outputs
Co   is    : start-of-region containing token contents, i.e.
Co         : token(is:is) is the first character after the terminator
Co         : If token is not matched, is=-99999
Co         : If token is matched but no '[' follows and itrm=1, is=-99998
Co         : If token is matched and itrm=1, but the matching ']' terminator
Co         :          cannot be found, is=-99997
Co   ie    : end-of-region containing ntk-th occurence of token contents
Co         : This index is EITHER:
Co         : (itrm=0) index to end-of-region.
Co         : (itrm>0) start of ntk+1 th occurence of token.
Co         : In either case, when marker is not found, ie=end of instr
Co         : -99999  token not found
Co         : -99998  missing start-of-region '['
Co         : -99997  missing end-of-region ']'
Cl Local variables
Cl         :
Cr Remarks
Cr   Find pointers is and ie of instr(is:ie) that demarcate token
Cr   contents.  How start-of region and end-of-region are determined depends
Cr   on itrm; see above.
Cr   Examples:
Cr      token        tokt   itrm  Matches
Cr     '@CAT@HAM'    ' '     0    '@CAT@HAM '
Cr     'SIG'         ' =:'   1    'SIG= [' or  'SIG: [' or 'SIG ['
Cr     'ATOM'        '='     2    'ATOM='  or 'ATOM= ['
Cu Updates
Cu   07 Aug 07 Adapted from region_finder (TK)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character(*),intent(in) :: instr,token,toks,tokt,eor
      integer is,ntk,lev,ie,itrm
C ... Local parameters
      integer j,i0,k,lentok,litrm,nnest,isave,i1,iprint,itrm2
      logical :: debug
      character(1) :: cc

      debug = iprint() .ge. 110
      ie = len(instr)
      lentok = len(token)

C --- Find ntk-th occurence of token ---
C     After this do loop:
C     1.  token+terminator was located
C     2.  If terminator is '[', nnest=1 and litrm=1
      is = 1
      nnest = 0
      itrm2 = mod(itrm,100)
      litrm = itrm2
      do  j = 1, ntk
   10   continue
C       i0 = 0 if token not found
        i0 = index(instr(is:ie),token)
C       No match ; exit
        if ( i0==0 ) then
          is = -99999
          if (debug) write(*,333) token, instr(1:min(10,len(instr)))
  333     format(' find_region: token `',a,
     .    ''' not found in string ',a,' ...')
          return
        endif
        is = is + i0-1 + lentok
C   ... One of toks must precede token; otherwise no match
        if (itrm .ge. 1000 .and. is-lentok .gt. 1) then
          do  k = 1, len(toks)
            if (instr(is-lentok-1:is-lentok-1) == toks(k:k)) goto 15
          enddo
          goto 10
   15     continue
        endif
C   ... Terminator must follow; otherwise no match
        if (itrm2 .gt. 10) then                ! Special case TOKEN[...
          cc = adjustl(instr(is:))
          if (cc .eq. '[') then
            is = is-1+index(instr(is:ie),cc)
            nnest = 1
            litrm = 1
            goto 20
          endif
        endif
        do  k = 1, len(tokt)
          if (instr(is:is) == tokt(k:k)) goto 20
        enddo
        goto 10
   20   continue
      enddo
      is = is+1

C --- Find is = start-of-region ---
C     If itrm=0, token terminator marks start-of-region
C     In this case, this branch is not executed.
C     If nnest>0, terminator was '[' which marks start-of-region
C     In this case, this branch is not executed.
C     In remaining cases, if the next nonblank character is '[',
C     it marks start-of-region
C     litrm is either 0 or 1 after this branch
      if (itrm2 .gt. 0 .and. nnest .eq. 0) then
        cc = adjustl(instr(is:ie))
        if (itrm2 .eq. 1 .or. itrm2 .eq. 11) then
          if (cc .ne. '[') then
            if (debug) write(*,"(a)") ' find_region: missing ''['' '//
     .      'after '//instr(is-lentok:is+1)//'...'
            is = -99998
            return
          else
            i0 = index(instr(is:ie),'[')
C           if (i0 .eq. 0) call rx('bug in find_region')
            is = is + i0
          endif
        elseif (itrm2 .eq. 2 .or. itrm2 .eq. 12) then
          if (cc .ne. '[') then
            litrm = 0
          else
            i0 = index(instr(is:ie),'[')
            is = is + i0
            litrm = 1
          endif
        else
          call rxi('illegal value for itrm:',itrm)
        endif
      endif

C --- Find ie = end-of-region.  Action depends on litrm ---
      if (litrm .eq. 0) then
        i0 = index(instr(is:ie),eor)
C       i0=0 => no eor was found => ie remains (length of instr)
        if (i0 .ne. 0) ie = is-1 + i0-1
C     Require that end-of-region correspond to ']' matching '['
      else
        isave = is
        nnest = 1
        do  while (nnest .gt. 0)
          i0 = index(instr(is:ie),']')
          if (i0 .eq. 0) then
            if (debug) write(*,"(a)") ' find_region: missing '']'' '//
     .      'after '//instr(isave-lentok:min(isave+10,ie))//'...'
            is = -99997
            return
          endif
          i1 = index(instr(is:ie),'[')
          if (i1 .gt. 0 .and. i1 .lt. i0) then
            is = is+i1
            nnest = nnest+1
          else
            is = is+i0
            nnest = nnest-1
          endif
        enddo
        ie = is-2
        is = isave
      endif

C ... Move start-of-region to first nonblank character
      if (ie .gt. is .and. mod(itrm,1000) .ge. 100) then
        cc = adjustl(instr(is:ie))
        if (cc .ne. ' ') then
          i1 = index(instr(is:ie),cc)
C         print *, instr(is:ie)
          is = is+i1-1
C         print *, instr(is:ie)
        endif
      endif

C ... Printout
      if (debug)
     .write(*,"(' find_region: contents of ', a,' : |',a,'|')")
     .token, instr(is:ie)
      end subroutine

