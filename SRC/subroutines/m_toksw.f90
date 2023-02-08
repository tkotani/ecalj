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


