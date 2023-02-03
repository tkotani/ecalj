!=========================================================================
module m_gtv
  use m_ftox
  !- module to read ctrl file -------------------------------
  !r gtv is the main rouitne. Its has interface.
  !r
  !r How gtv work, step by step.
  !r   0. call it with 'name' like HAM_RSRNGE and sw. See options below.
  !r   1. Af first, gtv check sw. If sw=2, exit with nout=0.
  !r   2. Then read ctrl. If readin number of data is nmin or more, normal exit.
  !r   3. If readin number of data is less, Error exit (stop) if sw=1
  !r   4.   If sw=0 (optional case) and default exist, read it and notmal exit
  !r   5.   If sw=0 (optional case) and no default, notmal exit with nout=0.
  !r
  !r   gtv( name, sw, data,
  !r       [default,cindx, note, or, nout, nmin, exist])
  !r  * default and below are optional.
  !r    So use f77 like binding for name, sw and data,
  !r    But use name-to-name binding for others.
  !r    Then f90 restriction, use defalut_*,where * is i4,i4vr8,r8v,lg.
  !r  * data can be i4 i4v r8 r8v logical char
  !r    In addition, no default in the case of char.
  !r  * In addition, basic settin in m_gtv is used: these are given by
  !r      subroutine gtv_setio(io_show_in, io_help_in)
  !r      subroutine gtv_setst(debug_in, stdo_in,stdl_in,stde_in)
  !r    These must be called before calling gtv.
  !r
  !i   name: like HAM_PMIN, SPEC_ATOM_Z and so on. _ devide tokens.
  !i         In principle, A1_A2_A3_A4 is readable, but
  !i         the grammer in ctrl file is problematic---so not so safe.
  !i     sw; sw=tksw(prgn,name)
  !o   data: this is first cleared as null=-99999.
  !o         data can be (i4 i4v r8 r8v logical char) types.
  !o         As f90, number of data(maximum number of readin data) is
  !o         automatically passed to.
  !i   def_* : default values. This can be an array for i4f,r8v.
  !i               As f90, number of array is recognized.
  !i   cindx: is now (/1,j/), j is a integer --- this is used for SPEC and SITE.
  !i          cindx specify j-th atom of second token (ATOM for SPEC_ATOM_Z).
  !i   note: note for input help
  !i     or: just for input help mode. Show "OR" if or=.true. default=.false.
  !o   nout: how many data are output. zero can be.
  !i   nmin: gtv fails to read ctrl if gtv can not read nmin or more data.
  !i         nmin=0 is possible. default=1.
  !o  exist: if the given name is found, exist=.true.
  ! ---------------------------------------------------------------
  implicit none
  integer,parameter,private :: nrcd=600000
  character(nrcd),private:: rcd  !ifc henry is strange
  !     character(nrcd),allocatable,private:: rcd(1)

  integer,private:: io_show,io_help,stdo,stdl,stde
  logical,private::  debug=.false.
  interface gtv
     module procedure &
          gtv_r8,gtv_i4,gtv_r8v,gtv_i4v,gtv_char,gtv_lg,gtv_none
  end interface gtv

  interface getinput
     module procedure &
          getinput_r8,getinput_i4,getinput_r8v,getinput_i4v, &
          getinput_char,getinput_none
  end interface getinput

contains
  subroutine gtv_setio(debug_in,io_show_in,io_help_in)
    !- Set io_show and io_help for m_gtv
    implicit none
    logical::  debug_in
    integer:: io_show_in,io_help_in
    debug   = debug_in
    io_show = io_show_in
    io_help = io_help_in
  end subroutine gtv_setio

  subroutine replacetab2space( recrd,recln )
    integer:: i,recln
    character*(*):: recrd
    do i = 1, recln
       if(recrd(i:i)==char(9)) recrd(i:i)=' '
    enddo
  end subroutine replacetab2space

  subroutine gtv_setrcd(recrd,nrecs,recln,stdo_in,stdl_in,stde_in)
    !! Copy recrd to rcd
    implicit none
    intent(in)::          recrd,nrecs,recln,stdo_in,stdl_in,stde_in
    integer:: stdo_in,stdl_in,stde_in
    ! ... Passed parameters
    integer :: nrecs,recln
    character*(recln):: recrd(nrecs),rrr
    ! ... Local parameters
    integer :: i,ioff,ioff2,j,l
    stdo = stdo_in
    stdl = stdl_in
    stde = stde_in
    rcd = ' '
    ioff2 = 1
    do  i = 1, nrecs
       call replacetab2space(recrd(i),recln) !takao replace tab to space. 7thJune2010.
       if(recrd(i)(1:1)=='#' .OR. recrd(i)(1:1)=='!') then
          continue
       elseif (recrd(i)(1:1)==' ') then !C       No character in 1st column
          rcd(ioff2:) = ' '//trim(adjustl(recrd(i)))
          ioff2 = ioff2 + len_trim(adjustl(recrd(i))) + 1
       else
          rrr=recrd(i)
          if (rrr(1:5)=='CLASS') then
             rrr(1:5) = 'SPEC '
          endif
          rcd(ioff2:) = ' @CAT@'//trim(adjustl(rrr))//' '
          ioff2 = len_trim(rcd)+1
       endif
    enddo
    rcd(ioff2:) = ' @CAT@EOF'
    !      write(6,*) 'xxxxxxxxxx rcd=',trim(rcd)
    ioff2 = ioff2+9
    if (ioff2 > nrcd) call rxi('m_gtv: increase size of rcd; need at least',ioff2)
  end subroutine gtv_setrcd

  subroutine gtv_entrance()
    !- General purpose routine to read token contents
    !  True entry points are gtv_r8,gtv_r8v,gtv_i4,gtv_i4v,gtv_lg,gtv_char,gtv_none
    ! ----------------------------------------------------------------------
    !i Inputs (required)
    !i  name   :string defining a token 'tree'.
    !i         :It delimits a region where data is to be read; See Remarks
    !i    sw   :sw=0: Token's presence is optional
    !i         :sw=1: Token's presence is required
    !i         :sw=2: Token's presence is ignored
    !i Inputs (optional)
    !i  cindx  :used when multiple instances of data are used for
    !i         :a particular token e.g. multiple species data.
    !i         :Token region corresponds to cindx-th occurrence of token
    !i         :If cindx exists, the how the end of the token's contents
    !i         :is determined is described in Remarks, case D.
    !i
    !i  note   :string for help printout
    !i
    !i    or   :if present, help printout flags that a later token
    !i         :may be alternately read to supply substitute information
    !i
    !i  nmin   :minimum number of elements that must be input EITHER from
    !i         :the token contents, OR supplied as default values.
    !i         :nmin's function depends on whether default values are supplied.
    !i         :Finally, whether nmin is used at all or not depends on sw
    !i         :See Examples below.
    !i         :nmin is supplied as one of these possibilities:
    !i         :  not present  or  nmin=0   or   nmin>0   or   nmin<0
    !i         :|nmin| corresponds to how many values must be returned.
    !i         :The sign of nmin has a special meaning; see below;
    !i
    !i         :Case defaults ARE supplied:
    !i         :1. nmin=(not present): gtv uses internally nmin=1
    !i                         gtv fills data after last parsed entry with
    !i                         default values.
    !i         :2. nmin=0      Not necessary that any values be parsed
    !i                         gtv fills data after last parsed entry with
    !i                         default values.
    !i         :3. nmin>0      Same as case 1.
    !i                         Note: nonsensical for nmin>nndef where
    !i                         nndef = number of default values passed.
    !i         :4. nmin<0      If a token is present, at least |nmin| values
    !i                         must be parsed, regardless of nndef
    !i
    !i         :Case defaults ARE NOT supplied:
    !i         :1. nmin=(not present): gtv uses internally nmin=1
    !i                         gtv fills data after last parsed entry with NULLI
    !i         :2. nmin=0      Not necessary that any values be parsed
    !i         :3. nmin>0      If a token is present, at least |nmin| values
    !i                         must be parsed
    !i         :4. nmin<0      Same as nmin>0
    !i
    !i         :Help mode: Pass nmin=NULLI together with io_show=1 =>
    !i         :size of input is unknown
    !i
    !i         :Character input: nmin has different meaning.
    !i         :1s digit: 0 => no string need be present
    !i         :        : 1 => a nonblank string must be read
    !i         :10s digit affects how character string is delimited
    !i         :Not used if quote marks delimit string ("'" or `"')
    !i         :Otherwise: 10s digit means:
    !i         :0  string delimited by token's region.  See Remarks below,
    !i             description below of *For character input, option 2.
    !i         :1  end-of-string delimited by first blank after start-of-string
    !i             See Remarks below, description below of *For character input,
    !i             option 3
    !i
    !i def_r8,def_r8v,def_i4,def_i4v,def_lg:
    !i         :Default values for input data.  If fewer than nmin elements
    !i         :are read from input, these default values may substitute,
    !i         :if they exist.  At most nmin elements in the default array
    !i         :are used.
    !i         :NB: Pass def_*(1)=NULL => default value depends on other input
    !i         :Works for real and integer input.
    !i
    !o Outputs
    !o dat,datv,idat,idatv,lg,char: are all the third argument of the
    !o         :generic entry point gtv.  They differ in their cast.
    !o         :Only one is used.  Input is stored into this argument.
    !o         :The dimension of this argument sets how many elements
    !o         :gtv will attempt to read.
    !o Outputs (optional)
    !o    nout :number of elements read into dat (or its equivalent),
    !o         :either from the input file, or assigned by default
    !o         :nout returned -1 if nothing sought (sw=2)
    !o  Texist :Returns information about whether token found or not
    !o         :If sw==2, always return Texist=F.  Otherwise:
    !o         :If io_help>0, always return Texist=T.  Otherwise:
    !o         :If token found, return Texist=T
    !o         :If token not found, return Texist=F
    !o         :Texist is only set if variable is present
    !l Local variables
    !l  ig     :string defining cast of object
    !l  nndef  :size of (number of elements in) default array
    !l  nminl  :nminl=|nmin|, if nmin exists; otherwise nminl=1
    !l  sizez  :number of elements in result array; number of elements to read
    !l  cindx2 :local copy of cindx(2), if it exists.
    !l         :otherwise, -1
    !r Remarks
    !r   'name' is a tree of tokens, with elements separated by underscores, e.g.
    !r     'HAM_TOL'  or  'SPEC_ATOM_Z'
    !r
    !r   The top-level token is called the 'category'.  In the second example
    !r   it is 'SPEC'; the second and third level tokens are 'ATOM' and 'Z'.
    !r   We will call the token of the preceding level the 'parent' token.  In
    !r   the example, 'ATOM' is the parent to 'Z'; 'SPEC' is the parent to 'ATOM'.
    !r
    !r   The input string is held in char variable  rcd .  rcd contains tokens
    !r   and their contents.  Following a token is a single character, the token
    !r   terminator, e.g. '=' (the terminator's value depends a little on the
    !r   context.)  The "data region" follows the token terminator, in a
    !r   substring rcd(i1:ie).  The rules for starting and end delimiters i1
    !r   and ie can depend on the context; see  *Finding a token and
    !r   *End delimiter below.
    !r
    !r   Tokens and their data regions are embedded within the region of the
    !r   parent.  Thus the region of 'SPEC_ATOM_Z' is contained within the
    !r   region of 'SPEC_ATOM,' which is in turn contained within the region
    !r   'SPEC'.  This means that the regions corresponding to tokens of the
    !r   same name but members of different trees, eg TOL in the HAM_TOL and
    !r   MIX_TOL, will be different and not overlap.
    !r
    !r   *Finding a token and the beginning of its data region
    !r   1.  Categories (top-level tokens) have a special syntax.
    !r       From the user's point of view, a new category begins whenever a
    !r       new line in the input file begins with a nonblank character.
    !r       The region of the category ends where the next new category starts.
    !r
    !r       The computer handles this by checking each input line.  Any line
    !r       which begins with a nonblank character has string @CAT@ prepended
    !r       to it.  Thus the start of a category with name 'NAM' is determined
    !r       by the first occurrence of the string '@CAT@NAM ', where the last
    !r       (terminating) character is a space.  The region ends just before
    !r       the first occurence of '@CAT@' following '@CAT@NAM '.
    !r
    !r   2.  Tokens below the top level must be found between the start and end
    !r       region of its parent token. A token must have a space preceding
    !r       it.  Thus, a token with nam 'NAM' is identified by the first occurence
    !r       of the string ' NAMx' within the region of the parent token.  Here x
    !r       is the terminating character; usually x is '=' or ':'.
    !r
    !r       In cases 1 and 2, a token's contents (data region) begin after
    !r       the token terminator.  What defines the start-of-region depends
    !r       on the context.  There are two types:
    !r       a. The first nonblank character following the token is a '['
    !r          In this case, start-of-region is the first character after '['
    !r       b. start-of-region begins immmediately after the token terminator
    !r       c. The calling program may use either of these rules, or both:
    !r          rule (a) applies if the first nonblank char is '['; if not
    !r          rule (b) applies.
    !r       Note: if only rule (a) is allowed, and first nonblank char
    !r       following the token is not '[', the parser exits with an
    !r       'error-match'
    !r
    !r   *End delimiter of token region.  It determined as follows:
    !r   A.  If start-of-region is determined by rule (a) above, the region
    !r       is terminated by ']'.  Since tokens may be nested in a tree
    !r       structure, the end delimiter must be the ']' matching its
    !r       corresponding '[' It is an error for any opening '['
    !r       to be missing its corresponding ']'.
    !r
    !r   B   If the token is top-level (category), end-of-region is defined
    !r       in notes (1) above.
    !r
    !r       Otherwise:
    !r
    !r   C.  With the exception noted in D below, the end of a token's region
    !r       coincides with the end of its parent.
    !r
    !r   D.  Certain tokens are used to tag multiple instances of data,
    !r       such as site positions or species data.  In these special
    !r       cases, the nth occurence of the token must also be specified
    !r       (cindx) and the end delimiter is the SMALLER of:
    !r         the endpoint of the parent token
    !r         the character preceding next occurence of the token
    !r       Example:   Given token 'ATOM'
    !r           ATOM=A  A-contents  ATOM=B  B-contents end  ATOM=C ...
    !r       the region of the 2nd occurence is ' B-contents end  '
    !r
    !r  --- Parsing of tokens, and how missing information is handled ---
    !r  A lot of flexibility is available to handle cases when tokens are
    !r  missing, or what to do when fewer than maximum elements are read.
    !r
    !r  *'sizez' is the maximum number of elements that can be read
    !r   If a token is found, gtv will attempt to read 'sizez' elements
    !r   of data (results of numerical expressions) from token contents.

    !r  *Character input is treated specially (number of elements is meaningless)
    !r   If the first character is a quote ("'" or `"'), the string is delimited
    !r   by pairs of quotation marks; see below.
    !r
    !r  *'nmin' dictates the minimum number of elements that must be input.
    !r    "Input" can be either from token contents or defaults; see below.
    !r    If 'nmin' is not supplied, it locally takes the value 1.
    !r
    !r  If gtv succeeds in reading only k elements, with k<nmin, gtv will fill
    !r  elements k+1..nmin with default values, if caller supplies them.
    !r  If caller supplies fewer than nmin default values, gtv cannot supply
    !r  nmin elements and program aborts.
    !r
    !r  Examples:
    !r  sw nmin  # defaults  Action:
    !r   1  2      none      Token must be found; 2 expressions must be parsed
    !r                       from token contents
    !r   1  1      1         Token must be found; 1 expressions must be found;
    !r                       since a default is available, it will substitute
    !r                       if no expression is found.  Note: it is not common
    !r                       to require a token be present with the expectation
    !r                       that an expression follow, yet to be no error if
    !r                       the expression is missing.
    !r                       Still, such instances may occur.
    !r   1  2      1         Same as first example.  This is a poor combination
    !r                       of parameters.  There are fewer than nmin
    !r                       defaults; so all data must be parsed anyway.
    !r                       Better to use first example, or supply at least
    !r                       nmin default values.
    !r   0  0      none      Token need not be present.  If it is present,
    !r                       no error occurs if no expressions are parsed.
    !r                       (The expression will evaluate to NULL)
    !r   0  1      none      Token need not be present.  If it is present,
    !r                       at least one expression must be parsed.
    !r                       Another unusual combination of conditions.
    !r   0  1      1         Token need not be present.  If it is, result is
    !r                       value of expression, if one is successfully parsed.
    !r                       Otherwise, result is the supplied default value.
    !r
    !r *For character input,  'character string' corresponds to 'expression'
    !r  in the description above.
    !r  A character string may be delimited in one of three ways:
    !r  1. If first nonblank character after start-of-region is a single
    !r     or double quote ("'" or `"'), string starts after quote,
    !r     ends before next occurence of quote, or end-of-region, whichever
    !r     comes first.
    !r  2. start-of-string is first nonblank character after start-of-region
    !r     end-of-string is end-of-region (this is the default)
    !r  3. start-of-string is same as 2.
    !r     end-of-string is first blank character after start-of-string
    !r     For this option, set 10's digit nmin=2; see descr. of nmin above
    !r
    !r  Thus, sw=0, nmin=0 =>Token need not be present.  If it is present,
    !r                       no error occurs if no string is found
    !r                       (The result is empty string)
    !r        sw=0, nmin=1 =>Token need not be present.  If it is present,
    !r                       a nonblank string must be read.
    !r  At present:
    !r    There is no capability to read vectors of character strings.
    !r    No default strings can be supplied
    !r
    !u Updates
    !u   27 Jul 07
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    character*(*),intent(in):: name
    real(8),intent(out)     :: dat,  datv(1:)
    integer,intent(out)  :: idat,idatv(1:)
    integer,intent(in)   :: sw
    logical,optional,intent(in)::  def_lg,or
    integer,optional,intent(in)::  nmin, cindx(1:)
    character*(*),optional,intent(in):: note
    integer,optional,intent(in)::  def_i4, def_i4v(1:)
    real(8),optional,intent(in)::     def_r8, def_r8v(1:)
    logical,intent(out)     :: lg
    integer,optional,intent(out):: nout
    logical,optional,intent(out):: Texist
    ! ... Local parameters
    integer::  cindx2
    integer:: nn,sizez,nparse,isw,iprint
    real(8):: ddat(1000)
    !     character(128*2) aaa
    character*(*):: char
    character(3)::ig
    character(256) sout
    !     character(21+124):: rrrx
    logical :: lgx,lexist,lor,lrqdn
    real(8):: defa(256)
    integer:: nndef
    character(4) opts(3)
    integer:: nminl, NULLI = -99999
    integer(2),intent(in):: nono
    !     A string for help printout
    data opts /'opt','reqd','skip'/

    logical:: isanrg, l_dummy_isanrg

    entry gtv_r8(name,sw,dat,def_r8,cindx,note,or,nmin,nout,Texist)
    nndef = 0; sizez = 1; defa(1) = 0d0
    if ( present(def_r8)  ) then
       nndef = 1; defa(1) = def_r8
    endif
    !       not sure we want to do this ...
    if (sw /= 2) dat = nulli
    ig = 'r8' ;    goto 990

    entry gtv_r8v(name,sw,datv,def_r8v,cindx,note,or,nmin,nout,Texist)
    nndef = 0; sizez = size(datv); defa(1) = 0d0
    if ( present(def_r8v) ) then
       nndef = min(size(def_r8v),sizez); defa= def_r8v
    endif
    if (sw /= 2) datv = nulli
    ig = 'r8v' ;   goto 990

    entry gtv_i4(name,sw,idat,def_i4,cindx,note,or,nmin,nout,Texist)
    nndef = 0; sizez = 1; defa(1) = 0d0
    if ( present(def_i4)  ) then
       nndef = 1; defa(1) = def_i4
    endif
    if (sw /= 2) idat = nulli
    ig = 'i4';    goto 990

    entry gtv_i4v(name,sw,idatv,def_i4v,cindx,note,or,nmin,nout, &
         Texist)
    nndef = 0; sizez = size(idatv); defa(1) = 0d0
    if (present(def_i4v) ) then
       nndef = min(size(def_i4v),sizez); defa= def_i4v
    endif
    if (sw /= 2) idatv = nulli
    ig='i4v';    goto 990

    entry gtv_lg(name,sw,lg,def_lg,cindx,note,or,nmin,nout,Texist)
    sizez = 1; nndef = 0
    if ( present(def_lg)  ) then
       nndef = 1 ; defa(1) = isw(def_lg)
    endif
    ig='lg' ;    goto 990

    entry gtv_char(name,sw,char,cindx,note,or,nmin,nout,Texist)
    nndef = 0
    sizez = 1;   ig='chr';   goto 990

    entry gtv_none(name,sw,nono,cindx,note,or,Texist)
    nndef = 0; defa(1) = 0d0
    sizez = 0;   ig='---';   goto 990

    ! --- Start of statements common to all entry routines ---
990 continue
    if (sw == 2) then
       if (ig /= '---') then
          if (present(nout)) nout=-1
       endif
       if (present(Texist) ) Texist = .FALSE. 
       !        if (debug) then
       !        endif
       return
    endif
    !     Assign cindx2 = cindx(2) for now
    cindx2 = -1
    if (present(cindx)) then
       !       if (sum(cindx)-cindx(2)/=size(cindx)-1) !only cont(2) allowed now
       if (cindx(1) /= 1) call rx('gtv:  cindx(1)>1 not implemented')
       cindx2 = cindx(2)
    endif
    ! ino isanrg is logical function,       call isanrg(sw,0,2,' gtv','sw',.true.)
    l_dummy_isanrg=isanrg(sw,0,2,' gtv','sw',.true.)
    nndef = min(nndef,sizez)  ! Only use defaults up to size of input

    ! ... Local copy of nmin that always exists
    nminl = min(1,sizez)
    !     DEC fortran: nmin can exis for gtv_none, even if not in arg lst
    lrqdn=.true.
    if (ig /= '---') then
       lrqdn = .false.
       if (present(nmin)) then
          !       nmin<0 => nparse MUST be at least nminl
          lrqdn = nmin .lt. 0
          nminl = iabs(nmin)
          if (io_help == 0 .AND. (ig /= 'chr')) &
                                ! ino isanrg is logical function,      .    call isanrg(iabs(nmin),-1,sizez,' gtv','nmin',.true.)
               l_dummy_isanrg=isanrg(iabs(nmin),-1,sizez,' gtv','nmin', .TRUE. )
       elseif (ig == 'chr') then
          nminl = 0
       endif
    endif
    ! ... local copy of present(or)
    lor = present(or)
    if (lor) then
       lor = or
    endif

    ! ... For printout
    if (debug .AND. io_help /= 0) then
       !        call info2(0,0,0,' gtv:  name='//name//'; cast='//ig//
       !     .    ' ... help mode',0,0)
!    elseif (debug) then
!       call info2(0,0,0,' gtv:  name='//name//'; cast='//ig// &
!            '%?!n>0!; contents from occurence #%-1j%i of token !!'// &
!            '%?!n>1!; attempt to read %-1j%i elements', &
!            cindx2,sizez)
    endif
    if (-nminl == NULLI) then
       write(sout, &
            "(1x,a,t20,a,3x,a,21x,'size depends on other input')") &
            name,opts(sw+1),ig
    elseif (sizez == 0) then
       write(sout,"(1x,a,t20,a,3x,a)") name,opts(sw+1),ig
    else
       write(sout,"(1x,a,t20,a,3x,a,i7,',',i3)") &
            name,opts(sw+1),ig,sizez,mod(nminl,10)
    endif

    ! ... Logical case: treat locally as integer
    lgx = .false.
    if (ig == 'lg') then
       ig = 'i4'
       lgx = .true.
    endif

    !  --- Help mode ---
    if (io_help /= 0) then
       if (lgx) then
          if(nndef==0) write(stdo,ftox)trim(sout)
          if(nndef==1) write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'default=',defa(1)/=0
       elseif (nndef >= 1 .AND. defa(1) == NULLI) then
          if (nminl == NULLI) then
             write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'size and defaults depend on other input'
          else
             write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'default depends on other input'
          endif
       elseif (nminl == NULLI .AND. nndef >= 1) then
          write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'... size depends on other input def=',defa(1:nndef)
       elseif (nminl == NULLI) then
          write(stdo,ftox) trim(sout)//repeat(" ",50-len(trim(sout)))//'size depends on other input'
       elseif (nndef >= 1 .AND. nndef <= 4) then
          if(sum(abs(defa(1:nndef)-nint(defa(1:nndef))))<1d-12) then
             write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'default=',nint(defa(1:nndef))
          elseif(sum(abs(defa(1:nndef)))<1d-2) then
             write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'default=',ftod(defa(1:nndef),3)
          else   
             write(stdo,ftox)trim(sout)//repeat(" ",50-len(trim(sout)))//'default=',ftof(defa(1:nndef),3)
          endif   
       else
          write(stdo,*)trim(sout)
       endif

       if (present(note)) then
          write(stdo,*)'   '//note
       endif
       if ( lor ) then
          write(stdo,"(a)") &
               ' * If token is not parsed, attempt to read the following:'
       endif
       if (ig /= '---') then
          if (present(nout)) nout=0
       endif
       if (present(Texist) ) Texist = .TRUE. 
       return
    endif

    ! --- Check for token match; no data input ---
    if (ig == '---') then
       call getinput(name, cindx2, lexist)
       if (present(Texist) ) Texist = lexist
       !   ... Printout
       if ((io_show>0 .OR. debug) .AND. iprint() /= 0) then
          if (lexist) then
             write(sout(1+len_trim(sout):),"(27x,'present')")
          else
             write(sout(1+len_trim(sout):),"(27x,'missing')")
          endif
          write(stdo,'(a)') trim(sout)
       endif
       if (sw == 1 .AND. .NOT. lor .AND. .NOT. lexist) then
          sout = ' gtv (abort): no token '//trim(name)//' found'
          call rx(sout)
       endif
       return
    endif

    ! --- Character input ---
    if (ig == 'chr') then
       char = ' '
       !        call getinput(name, char, 1, cindx2, lexist, nn)
       call getinput(name, char, nminl/10, cindx2, lexist, nn)
       !       print *, sout
       if (( .NOT. lexist .AND. .NOT. lor .AND. sw == 1)) then
          sout = ' gtv (abort): no token '//trim(name)//' found'
          call rx(sout)
       elseif (lexist .AND. mod(nminl,10) > 0 .AND. nn == 0) then
          sout = ' gtv (abort): no string for token '//trim(name)
          call rx(sout)
       endif
       if (present(Texist) ) Texist = lexist
       if (present(nout) ) nout = nn
       !   ... Printout
       if (io_show>0 .OR. debug) then
          if ( .NOT. lexist) then
             write(sout(1+len_trim(sout):),"(',   *')")
          elseif (nn == 0) then
             write(sout(1+len_trim(sout):),"(',   0')")
          else
             write(sout(1+len_trim(sout):),"(',   1')")
             sout(57:) = char
          endif
          write(stdo,'(a)') trim(sout)
       endif
       return
    endif

    ! --- Get token contents; try to read d.p. vector of size sizez ---
    !     nparse = number of elements actually read
    call getinput(name, ddat, sizez, cindx2, lexist, nparse)
    if (present(Texist)) Texist = lexist
    nn = nparse
    ! ... Case fewer values read than sought
    if ( nparse < nminl ) then
       !       No error, if defaults are available to fill nparse+1 ... nmin
       if (nndef >= nminl .AND. (sw == 0 .OR. lexist) .AND. &
            .NOT. lrqdn) then
          nn = nndef
          ddat(nparse+1:nn) = defa(nparse+1:nn)
          !       No error if token's presence not required
       elseif ((sw == 0 .OR. lor) .AND. .NOT. lexist) then
          continue
          !       Otherwise, error exit
       else
          if ( .NOT. lexist) then
             sout = ' gtv (abort): no token '//trim(name)//' found'
          elseif (nminl == 1) then
             sout = ' gtv (abort): no expression read for token '//name
          else
             !call info(0,0,0,' gtv: parsed %i%-1j elements: %n:1g',nn, &
             !     ddat)
             write(opts(1),"(i4)") nn
             write(opts(2),"(i4)") nminl
             sout = ' gtv (abort): only '// trim(adjustl(opts(1))) // &
                  ' expressions(s) read for ' // name // &
                  ' when ' // trim(adjustl(opts(2))) // ' required'
          endif
          call rx(sout)
       endif
    elseif ( nndef > nminl .AND. nparse < nndef ) then
       nn = nndef
       ddat(nparse+1:nn) = defa(nparse+1:nn)
    endif
    ! ... Copy result to one of (dat,idat,datv,idatv,lg)
    if (present(nout)) nout = nn
    if (ig=='r8' .AND. nn==1)      then ;   dat  = ddat(1)
    elseif (lgx .AND. nn==1)       then ;   lg = nint(ddat(1)) /= 0
    elseif (ig=='i4' .AND. nn==1)  then ;   idat = ddat(1)
    elseif (ig=='r8v' .AND. nn > 0) then ; datv(1:nn)  = ddat(1:nn)
    elseif (ig=='i4v' .AND. nn > 0) then ; idatv(1:nn) = ddat(1:nn)
    endif
    ! ... Printout
    if (io_show>0 .OR. debug) then
       if ( .NOT. lexist .AND. nndef == 0) then
          write(sout(1+len_trim(sout):),"(',   *, --')")
       elseif (nn == nparse .AND. nndef /= 0) then
          write(sout(1+len_trim(sout):),"(',',i4,',',i3)")nparse,nn-nparse
       elseif (nn == nparse) then
          write(sout(1+len_trim(sout):),"(',',i4,', --')")nparse
       elseif ( .NOT. lexist) then
          write(sout(1+len_trim(sout):),"(',   *,',i3)")         nn-nparse
       else
          write(sout(1+len_trim(sout):),"(',',i4,',',i3)")nparse,nn-nparse
       endif
       if (nparse > 0 .AND. lgx)       write(stdo,ftox) trim(sout),'     ',lg
       if (nparse > 0 .AND. .NOT. lgx) write(stdo,ftox) trim(sout),'     ',ftom(ddat(1:nn))
       !        if (nn .gt. 0 .and. lgx)       write(stdo,ftox) trim(sout),'     ',lg
       !        if (nn .gt. 0 .and. .not. lgx) write(stdo,ftox) trim(sout),'     ',ftom(ddat(1:nn))
    endif
1012 format(a,d13.5)
2012 format(a,200d13.5)
  end subroutine gtv_entrance

  subroutine getinput_entrance()
    !- Find token and read contents (true entry points are below)
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  name   :Name of token, including parents
    !i  nin    :number of arguments to read from token
    !i         :For character input, nin means the following:
    !i         :Not used if quote marks delimit string ("'" or `"')
    !i         :Otherwise:
    !i         :0  string delimited by token's region.
    !i         :1  start-of-string delimited by first non-blank character in region
    !i         :   end-of-string delimited by first blank after start-of-string.
    !i  cindx2 :used to indicate multiple occurences of a token
    !i         :If cindx2>0, use cindx2-th occurence of token
    !i         :Otherwise, cindx2 should be -1
    !i         :Note: the syntax for delimiting the token's region
    !i         :can be different when cindx2>0; see Remarks in
    !i         :subroutine gtv_entrance above.
    !o Outputs
    !o dat,datv,idat,idatv,char: are all the second argument of the
    !o         :generic entry point getinput.  They differ in their cast.
    !o         :Only one is used.  Input is stored into this argument.
    !o  Texist :.true. if the token could be matched; otherwise .false.
    !o Outputs (optional)
    !o   nout  :number of elements actually read
    !o         :For character input,
    !o         :nout is 0 if no nonblank string found
    !o         :nout is 1 if a nonblank string is found
    !l Local variables
    !l         :
    !r Remarks
    !r
    !u Updates
    !u   10 Aug 07
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer,intent(in) :: nin, cindx2
    character(*),intent(in) :: name
    character(*),intent(out):: char
    real(8),intent(out) :: dat,  datv(1:)
    integer,intent(out) :: idat, idatv(1:)
    integer,optional,intent(out) :: nout
    logical,intent(out) :: Texist
    ! ... Local parameters
    character(1024) :: nameb
    ! ino del spoint
    integer :: i1,iend,i,ilev,iix,iex,itrm, &
         mxlev,ie,ii0,ii
    ! ino      external :: spoint
    integer ::  n
    character(3) :: ig
    character(1) :: keye
    character(5) :: head="@CAT@"
    integer :: iarr(1000)
    real(8) :: arr(1000)

    character(50),save:: tokencut(10)
    integer :: a2vec

    entry getinput_r8 (name,  dat,nin,cindx2,Texist,nout)
    ig = 'r8' ;  goto 990
    entry getinput_r8v(name, datv,nin,cindx2,Texist,nout)
    ig = 'r8v' ; goto 990
    entry getinput_i4 (name, idat,nin,cindx2,Texist,nout)
    ig = 'i4' ;  goto 990
    entry getinput_i4v(name,idatv,nin,cindx2,Texist,nout)
    ig = 'i4v';  goto 990
    entry getinput_char(name,char,nin,cindx2,Texist,nout)
    ig = 'chr'; goto 990
    entry getinput_none(name,cindx2,Texist,nout)
    ig = '---'; goto 990

    ! --- Start of statements common to all entry routines ---
990 continue

    ! --- Split token=string1[_string2[_string3...]]; -> tokencut(1..mxlev) ---
    ilev = 0
    nameb = name
    mxlev = 0
    do  !token cut by '_'
       ilev=ilev+1
       iend  = index(nameb,'_') - 1
       if (iend == -1) then
          iend = len_trim(nameb) + 1
          mxlev = ilev
       endif
       if (ilev/=1) tokencut(ilev) = adjustl(nameb(1:iend))
       if (ilev==1) tokencut(ilev) = head//adjustl( nameb(1:iend) )
       nameb = nameb(iend+2:)
       if (mxlev/=0) exit
    enddo
    if (debug) &
         write(stdo,"(10a)") ' getinput: ',name,' partitioned into: ', &
         (trim(tokencut(i))//'->',i=1,mxlev-1), trim(tokencut(mxlev))

    ! --- Find region of token's contents ---
    !     Start from top-level token in tree structure to define
    !     region of that level.  A region of a level must be contained
    !     in the region of the prior level.
    if (present(nout) ) nout = 0

    i1 = 1
    ie = len_trim(rcd)
    do  ilev = 1, mxlev
       itrm = 12
       if (ig == 'chr') itrm = itrm+100
       if (debug) call pshpr(111)
       !       Categories: Terminator=' '  eor=head
       if (ilev == 1) then
          call find_region(ilev,rcd(i1:ie),trim(tokencut(ilev)),' ', &
               ' ',1,itrm,head,iix,iex)
       elseif (ilev == 2 .AND. cindx2 > 0) then
          call find_region(ilev,rcd(i1:ie),trim(tokencut(ilev)),' ', &
               ':=',cindx2,1000+itrm,trim(tokencut(ilev)),iix,iex)
          !       If below highest nesting level, terminator '[' is required
       elseif (ilev < mxlev) then
          call find_region(ilev,rcd(i1:ie),trim(tokencut(ilev)),' ', &
               ' ',1,1011,head,iix,iex)
          !       All other tokens: eor=head and Terminator=':='
       else
          call find_region(ilev,rcd(i1:ie),trim(tokencut(ilev)),' ', &
               ':=',1,1000+itrm,head,iix,iex)
       endif
       if (debug) call poppr
       !       Error exit if token//'[' found has no matching ']'
       if (iix == -99997) call rx('getinput: token '// &
            trim(tokencut(ilev))//'[ missing its corresponding '']''')
       !       Exit if token not found
       if (iix < -99990) then
          if (present(nout) ) nout = 0
          Texist = .false.
          return
       endif
       !       Narrow region rcd(i1:ie)
       ii0 = i1
       i1 = ii0-1 + iix
       ie = ii0-1 + iex
    enddo
    Texist = .true.

    ! --- Token has no arguments ---
    if ( ig == '---' ) then
!       if (debug) call info0(0,0,0, &
!            ' getinput: found token '//trim(name))
       return
    endif

    ! --- Character input ---
    if ( ig == 'chr' ) then
       if (len(char) == 0) then        ! Input string of null length
          if (present(nout) ) nout = 0
       elseif (ie < i1) then          ! Token region of null length
          char = ' '
          if (present(nout) ) nout = 0
          !   ... Normal case: fix delimiters
       else
          if (present(nout) ) nout = 1
          if (rcd(i1:i1) == '''') then    ! string delimited by '...'
             ii = index(rcd(i1+1:),'''')
             if (ii > 0) then
                char = rcd(i1+1:i1+ii-1)
             else
                char = rcd(i1+1:ie)
             endif
          elseif (rcd(i1:i1) == '"') then ! string delimited by "..."
             ii = index(rcd(i1+1:),'"')
             if (ii > 0) then
                char = rcd(i1+1:i1+ii-1)
             else
                char = rcd(i1+1:ie)
             endif
             !     ... String may be further delimited, depending on nin
          else
             if (nin == 1) then   ! reduce the range of i1:ie
                if (rcd(i1:i1) == ' ') then  ! Shift i1 to 1st nonblank
                   keye = adjustl(rcd(i1:ie))
                   n = index(rcd(i1:ie),keye)
                   if (n > 0) i1 = i1 + n-1
                endif
                n = index(rcd(i1:ie),' ')
                if (n > 0) ie = i1 + n-1
             elseif (nin > 1) then
                call rxi('getinput: illegal value, nin=',nin)
             endif
             char = rcd(i1:ie)
          endif
       endif
       return
    endif

    ! --- ASCII-numerical conversion. Simple matematical conversion.
    ii = 0
    n = a2vec(rcd(i1:ie),ie-i1+1,ii,4,', ',2,-3,nin,iarr,arr)
    if (n < 0) n = -n-1
!  if(debug) write(stdo,ftox)' getinput: sought',nin,'numbers. Read',n,'from'//trim(name)
    if (present(nout) ) nout = n
    ! --- Copy array to data ---
    if (n == 0) then;
    elseif (ig=='r8')  then ;   dat = arr(1)
    elseif (ig=='i4')  then ;   idat = arr(1)
    elseif (ig=='r8v') then;  datv(1:n) = arr(1:n)
    elseif (ig=='i4v') then;  idatv(1:n) = arr(1:n)
    endif
  end subroutine getinput_entrance

end module m_gtv
