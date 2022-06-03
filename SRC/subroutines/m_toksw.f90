module m_toksw
  !- Module to create token switches, and read contents of tokens
  ! ----------------------------------------------------------------------
  !r First, a list of tokens must be created, one list each for
  !r any program which uses this module.
  !r Suppose a particular program 'foo' reads tokens
  !r   STRUC_ALAT STRUC_PLAT SITE_ATOM_POS HAM_TOL ...
  !r Further suppose that HAM_TOL is "optional", i.e. no error occurs
  !r if the HAM_TOL is not found
  !r Create the token list corresponding to 'foo' as follows:
  !r    call nwsadd()  ! starts list for a new program
  !r    call tkadd(" foo::" ) ! Identifies foo as progn associated w/ this list
  !r    call tkadd(" STRUC_PLAT") ! Append this token to the list
  !r    call tkadd(" HAM_TOL~")   ! ~ => token can be optionally read
  !r    call tkadd(" STRUC_PLAT SITE_ATOM_POS") ! lump several tokens in 1 call
  !r    etc
  !r
  !r Examples see how routine toksw_init in m_rdctrl.f bundles the above
  !r calls for a particular set of programs, e.g. "LMF", "LM", etc.
  !r
  !r After the lists have been created, function tksw(foo,tok) returns
  !r a value 0,1,2 indicating :
  !r   0 tok is sought by foo, with input optional
  !r   1 tok is sought by foo, with input required
  !r   2 tok is is not sought by foo
  !r
  !r The input file reader contains commands to read the contents of all
  !r tokens any program using this input system will seek.
  !r It calls tksw to determine whether a particular token is be read
  !r by the given calling program.
  !r
  !r This module also contains module m_gtv, the module that read and
  !r parses tokens and tokens' contents.  The entry points are:
  !r   gtv_r8,gtv_r8v,gtv_i4,gtv_i4v,gtv_lg,gtv_char,gtv_none
  !r and are usually invoked by the interface gtv, e.g. this call:
  !r     call gtv('HEADER',tksw('LMF','HEADER'),header,note=
  !r    .  'Contents displayed at beginning of program execution')
  !r reads string 'header'
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Module variables
  logical,private:: debug=.false.
  integer,parameter,private :: lenmax=3000,nswmx=25
  character(lenmax),private:: swtok(nswmx)
  integer,private ::  nsw=0  !lenc=0,

contains

  subroutine clear_swtok(debug_in)
    !- Resets all tokens in swtok and sets internal debug switche
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   debug_in : sets internal debug switch for future calls
    !o Outputs
    !m Module variables
    !m  swtok : is initialized
    !m  nsw   : is initialized to zero
    !m  debug : is assigned to debug_in
    !r Remarks
    !r
    !u Updates
    !u   16 Jan 09
    ! ----------------------------------------------------------------------
    ! ... Passed parameters
    logical:: debug_in
    ! ... Module variables
    debug = debug_in
    swtok = ' '
    nsw = 0
  end subroutine clear_swtok

  subroutine nswadd()
    !- Flag that all tokens have been appended for current main program
    ! ----------------------------------------------------------------------
    !i Inputs
    !o Outputs
    !m Module variables
    !m  nsw : nsw is incremented, starting a new
    !m      : entry in character array swtok
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none

    nsw = nsw+1
    if (nsw>nswmx) call rx('m_toksw: enlarge nswmx')
  end subroutine nswadd

  subroutine tkadd(tok)
    !- Append tok to current set in swtok
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   tok  : Append 'tok' to swtok(nsw)
    !i        : Note: nsw must be set previously by calling nswadd()
    !o Outputs
    !o  swtok : tok is appended to (module string) swtok(nsw)
    !l Local variables
    !r Remarks
    !r   This routine builds up character array swtok for a particular
    !r   main program.
    !r
    !r   swtok contains a set of lists of tokens, one list
    !r   for each main program, which are read by the main program,  e.g.
    !r     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
    !r     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
    !r   The list for each main program contains a sequence of strings,
    !r   separated by a space, each string corresponding to the full name
    !r   of a token (incl parents).  Optionally the name may be appended
    !r   by a "~" which signifies the token can be optionally read.
    !r
    !r   swtok is built up by the following sequence:
    !r     call nwsadd()  ! starts new list for another program
    !r     call tkadd(" prgn::" )  ! prgn = main program name, e.g. "LMF"
    !r     call tkadd(" tok")      ! tok = token name, e.g. SITE_ATOM_POS
    !r     call tkadd(" tok")      ! etcetera
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    character*(*):: tok
    ! ... Local parameters
    integer:: lensw

    if (nsw > nswmx) call rx('m_toksw: enlarge nswmx')
    lensw = len(trim(swtok(nsw)))
    if ( lensw+30>lenmax ) call rx('m_toksw: enlarge lenmax')
    swtok(nsw) = trim(swtok(nsw))//trim(tok)

    !     lenc = max(lenc,lensw)  ! lenc (to check the length of buff)
    !     print *, 'toksw for nsw=',nsw,' is now:', trim(swtok(nsw))
  end subroutine tkadd

  integer function tksw(prgn,name)
    !- Return 0, 1 or 2 for given program and token name
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   prgn : main program name
    !i        : From a list of tokens 'prgn' uses (see Remarks),
    !i        : determine whether the supplied name is in the list.  Return:
    !i        : 0 if 'name'~ is present => input sought by 'prgn', but input
    !i            is optional
    !i        : 1 if 'name' is present => input is sought by 'prgn'
    !i        : 2 if 'name' is missing => input not sought by 'prgn'
    !i  name  : full name of token (incl parents) e.g. SITE_ATOM_POS
    !r Remarks
    !r   This routine searches character array swtok (which must be
    !r   generated in advance) to determine how a program is to treat an
    !r   input token.
    !r
    !r   swtok contains a set of lists of tokens, one list
    !r   for each main program, which are read by the main program,  e.g.
    !r     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
    !r     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
    !r   The list for each main program contains a sequence of strings,
    !r   separated by a space, each string corresponding to the full name
    !r   of a token (incl parents).  Optionally the name may be appended
    !r   by a "~" which signifies the token can be optionally read.
    !------------------------------------------------------
    character*(*):: name,prgn
    integer:: isw,ix,ixe,sw,ix1,ix2,npgsz=30

    ! ... Find which of swtok corresponds to prgn
    do  isw = 1, nsw
       ix = index(swtok(isw)(1:npgsz),trim(prgn)//'::')
       ! i        if (ix/=0 .and. swtok(isw)(ix-1:ix-1) .eq. ' ') goto 100
       if (ix/=0) then
          if (swtok(isw)(ix-1:ix-1) == ' ') goto 100
       endif
    enddo
    call rx('tksw: no switches for program '//prgn)
100 continue

    !      ix= index(swtok(isw),name)
    !      print *,'tksw: name= '//'###'//name//'###',ix

    !     See if found among required tokens
    ix1 = index(swtok(isw),' '//trim(name)//' ')
    !     See if found among optional tokens
    ix2 = index(swtok(isw),' '//trim(name)//'~')
    
    ! ... Determine whether optional or required
    ix=-999999
    if (ix1 == 0 .AND. ix2 == 0) then
       ix = 0
    elseif (ix1 == 0 .AND. ix2 /= 0) then
       ix = ix2+1
    elseif (ix2 == 0 .AND. ix1 /= 0) then
       ix = ix1+1
    else
       call rx('tksw, prgn '//trim(prgn)//': '//trim(name)// &
            ' appears in both optional and required lists')
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
    !      print *,'tksw: ', trim(prgn),' ',trim(name),sw
    tksw = sw
  end function tksw
end module m_toksw


subroutine find_region(lev,instr,token,toks,tokt,ntk,itrm,eor, &
     is,ie)
  !- Finds substring corresponding to ntk-th token.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lev   :token level, reserved for acceleration (not used now).
  !i   instr :input string, where token and contents are embedded
  !i   token :token to match in instr; see Remarks
  !i   toks: :(only used if 1000s digit of itrm is set)
  !i         :list of characters that must precede token to satisfy a
  !i         :match.  Eliminates strings matching token embedded within
  !i         :other string.  Example:
  !i         :Example: given token ATOM, and
  !i         :toks='[ ' matches ' ATOM' and '[ATOM] but not 'CATOM'
  !i   tokt: :terminator(s) to token
  !i         :Example: given token ATOM,
  !i         :tokt='=:' matches matches ATOM= and ATOM:
  !i   ntk   :search for ntk-th occurence of token
  !i   itrm  :(itrm=0) start-of-region begins with second character following
  !i         :         token (1st char is token terminator)
  !i         :         end-of-region is delimited by first occurence of string
  !i         :         contained in eor following start-of-region
  !i         :
  !i         :(itrm=1) the first nonblank character after the token and its
  !i         :         terminator must be '['.  Note: '[' may also server as
  !i         :         start-of-region begins with the character following it.
  !i         :         end-of-region is delimited by the matching ']'
  !i         :         Note: '[...]' pairs may be nested; it is an error for
  !i         :         '[' not to have a matching ']'.
  !i         :
  !i         :(itrm=2) If the first nonblank character after the token and its
  !i         :         terminator is '[', follow the syntax of itrm=1.
  !i         :         Otherwise, follow the syntax of itrm=0.
  !i         :(itrm=11)Identical to itrm=1, with the addition that '[' may
  !i         :         serve both as terminator and start-of-region marker
  !i         :(itrm=12)Identical to itrm=2, with the addition that '[' may
  !i         :         serve both as terminator and start-of-region marker
  !i
  !i         :Adding 100 to itrm causes find_region to move
  !i         :start-of-region marker to 1st nonblank character
  !i         :Adding 1000 to itrm turns on the pre-token matching;
  !i         :see toks above
  !i
  !i   eor   :string demarcating end-of-region.  Its use depends on
  !i         :how start-of-region was determined (see itrm)
  !i         :If start-of-region is NOT determined by '[',
  !i         :the first string after start-of-region that matches the
  !i         :contents of eor demarcates end-of-region
  !i         :If start-of-region IS determined by '[', eor is not used
  !o Outputs
  !o   is    : start-of-region containing token contents, i.e.
  !o         : token(is:is) is the first character after the terminator
  !o         : If token is not matched, is=-99999
  !o         : If token is matched but no '[' follows and itrm=1, is=-99998
  !o         : If token is matched and itrm=1, but the matching ']' terminator
  !o         :          cannot be found, is=-99997
  !o   ie    : end-of-region containing ntk-th occurence of token contents
  !o         : This index is EITHER:
  !o         : (itrm=0) index to end-of-region.
  !o         : (itrm>0) start of ntk+1 th occurence of token.
  !o         : In either case, when marker is not found, ie=end of instr
  !o         : -99999  token not found
  !o         : -99998  missing start-of-region '['
  !o         : -99997  missing end-of-region ']'
  !l Local variables
  !l         :
  !r Remarks
  !r   Find pointers is and ie of instr(is:ie) that demarcate token
  !r   contents.  How start-of region and end-of-region are determined depends
  !r   on itrm; see above.
  !r   Examples:
  !r      token        tokt   itrm  Matches
  !r     '@CAT@HAM'    ' '     0    '@CAT@HAM '
  !r     'SIG'         ' =:'   1    'SIG= [' or  'SIG: [' or 'SIG ['
  !r     'ATOM'        '='     2    'ATOM='  or 'ATOM= ['
  !u Updates
  !u   07 Aug 07 Adapted from region_finder (TK)
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  character(*),intent(in) :: instr,token,toks,tokt,eor
  integer :: is,ntk,lev,ie,itrm
  ! ... Local parameters
  integer :: j,i0,k,lentok,litrm,nnest,isave,i1,iprint,itrm2
  logical :: debug
  character(1) :: cc

  debug = iprint() .ge. 110
  ie = len(instr)
  lentok = len(token)

  ! --- Find ntk-th occurence of token ---
  !     After this do loop:
  !     1.  token+terminator was located
  !     2.  If terminator is '[', nnest=1 and litrm=1
  is = 1
  nnest = 0
  itrm2 = mod(itrm,100)
  litrm = itrm2
  do  j = 1, ntk
10   continue
     !       i0 = 0 if token not found
     i0 = index(instr(is:ie),token)
     !       No match ; exit
     if ( i0==0 ) then
        is = -99999
        if (debug) write(*,333) token, instr(1:min(10,len(instr)))
333     format(' find_region: token `',a, &
             ''' not found in string ',a,' ...')
        return
     endif
     is = is + i0-1 + lentok
     !   ... One of toks must precede token; otherwise no match
     if (itrm >= 1000 .AND. is-lentok > 1) then
        do  k = 1, len(toks)
           if (instr(is-lentok-1:is-lentok-1) == toks(k:k)) goto 15
        enddo
        goto 10
15      continue
     endif
     !   ... Terminator must follow; otherwise no match
     if (itrm2 > 10) then                ! Special case TOKEN[...
        cc = adjustl(instr(is:))
        if (cc == '[') then
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

  ! --- Find is = start-of-region ---
  !     If itrm=0, token terminator marks start-of-region
  !     In this case, this branch is not executed.
  !     If nnest>0, terminator was '[' which marks start-of-region
  !     In this case, this branch is not executed.
  !     In remaining cases, if the next nonblank character is '[',
  !     it marks start-of-region
  !     litrm is either 0 or 1 after this branch
  if (itrm2 > 0 .AND. nnest == 0) then
     cc = adjustl(instr(is:ie))
     if (itrm2 == 1 .OR. itrm2 == 11) then
        if (cc /= '[') then
           if (debug) write(*,"(a)") ' find_region: missing ''['' '// &
                'after '//instr(is-lentok:is+1)//'...'
           is = -99998
           return
        else
           i0 = index(instr(is:ie),'[')
           !           if (i0 .eq. 0) call rx('bug in find_region')
           is = is + i0
        endif
     elseif (itrm2 == 2 .OR. itrm2 == 12) then
        if (cc /= '[') then
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

  ! --- Find ie = end-of-region.  Action depends on litrm ---
  if (litrm == 0) then
     i0 = index(instr(is:ie),eor)
     !       i0=0 => no eor was found => ie remains (length of instr)
     if (i0 /= 0) ie = is-1 + i0-1
     !     Require that end-of-region correspond to ']' matching '['
  else
     isave = is
     nnest = 1
     do  while (nnest .gt. 0)
        i0 = index(instr(is:ie),']')
        if (i0 == 0) then
           if (debug) write(*,"(a)") ' find_region: missing '']'' '// &
                'after '//instr(isave-lentok:min(isave+10,ie))//'...'
           is = -99997
           return
        endif
        i1 = index(instr(is:ie),'[')
        if (i1 > 0 .AND. i1 < i0) then
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

  ! ... Move start-of-region to first nonblank character
  if (ie > is .AND. mod(itrm,1000) >= 100) then
     cc = adjustl(instr(is:ie))
     if (cc /= ' ') then
        i1 = index(instr(is:ie),cc)
        !         print *, instr(is:ie)
        is = is+i1-1
        !         print *, instr(is:ie)
     endif
  endif

  ! ... Printout
  if (debug) &
       write(*,"(' find_region: contents of ', a,' : |',a,'|')") &
       token, instr(is:ie)
end subroutine find_region

