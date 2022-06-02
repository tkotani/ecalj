module m_rdfiln
  public M_rdfiln_init, recln,recrd,nrecs,mxrecs

  integer,parameter:: recln=511, mxrecs=10000
  character*(recln),protected:: recrd(mxrecs)
  integer,protected:: nrecs

  private
contains
  subroutine m_rdfiln_init() !Read recln,recrd,nrecs,mxrecs for m_lmfinit_init
    !     -vfoobar replaced simplified contents of ctrl into recrd
    use m_lgunit,only:stdo
    use m_ext,only: sname
    integer:: nfilin,master=0,ierr,i,procid=0,mpipid
    character(8) :: alabl
    character:: strn*1000
    logical:: fileexist,lshowp,lshow,cmdopt0,mlog
    include "mpif.h"
    procid = mpipid(1)
    mlog = cmdopt0('--mlog')  !! set log for --mlog (not maintained well)
    !      stdo = lgunit(1)
    if(procid==master) then
       inquire(file='ctrl.'//trim(sname),exist=fileexist)
       if( .NOT. fileexist) call rx("No ctrl file found! ctrl."//trim(sname))
       open(newunit=nfilin,file='ctrl.'//trim(sname))
       call findctrlstart(nfilin) ! if a tag 'ctrlstart' in ctrl, ctrl is read from the tag.
       alabl = '#{}% ct '
       call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln,nrecs) !read ctrl into recrd
       close(nfilin)
    endif
    call mpibc1( nrecs,1,2,mlog,'main','nrecs')
    call MPI_BCAST( recrd,recln*(nrecs+1),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
    !! Show or not show readin ctrl file.
    lshowp = cmdopt0('--showp')
    lshow  = cmdopt0('--show')
    if(procid==master .AND. (lshow .OR. lshowp) ) then
       write(stdo,"('---- preprocessed ctrl file -------')")
       do i = 1, nrecs
          if(trim(recrd(i))=='') cycle
          write(stdo,"(a)")'%% '//trim(recrd(i))//' %%'
       enddo
       write(stdo,"(a,' preprocesses ctrl file. nrecs=', i5)")trim(sname),nrecs
    endif
    if(lshowp) call Rx0('end of --showp')
    !  recrd already replace foobar with foobar2 if we have an command-line option
    !    -vfoobar=foobar2 option (replacement of %const defined in ctrl file).
    ! Following ctrlrecrd is just for check recrd contains simplified ctrl with -vfoobar replacement.
    !      do i=1,nrecs
    !         if(len_trim(recrd(i))/=0) write(6,"(a)") trim(recrd(i))
    !     enddo
  end subroutine m_rdfiln_init

  subroutine parchv(recrd,recl,mxchr,cex,sep,opts,nchr,ctbl,j)
    !- Parses a string for one or more character variable declarations
    ! ----------------------------------------------------------------
    !i Inputs
    !i   recrd(0:*): string recrd is parsed from j to recl-1
    !i   recl :  record length: characters limited to recrd(0:recl-1)
    !i   mxchr:  maximum number of entries in character table
    !i   cex  :  characters demarcating string substitution, eg '{}'
    !i   sep  :  characters demarcating symbol and string-value; see Remarks
    !i   opts:   0: var declaration of existing variable is supressed
    !i           1: var declaration of existing variable supersedes
    !i        10's digit:
    !i           1  variable's value is substituted for shell environment
    !i              variable of the same name
    !i   nchr :  number of variables defined so far
    !i   ctbl :  table of character variables (col 1) and values (col 2)
    !i   j:      parsing starts at char j (origin at 0)
    !o Outputs
    !o   nchr :  is updated
    !o   j:      last character parsed
    !r Remarks
    !r   Declarations are typically of the form symbol=string-value
    !r   The separator '=' need not be '=', but is defined by 'sep'
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: recl,nchr,j,mxchr,opts,ctlen
    parameter (ctlen=120)
    character(1) :: recrd(0:*),ctbl(mxchr,2)*(ctlen),sep*(*),cex(2)
    ! Local parameters
    integer :: i,k,j0,jr,opt1,opt2
    character ctmp*1,strn*(ctlen)

    opt1 = mod(opts,10)
    opt2 = mod(opts/10,10)

    ! --- Parse for next declaration ---
2   call skipbl(recrd,recl,j)
    if (j >= recl .OR. nchr >= mxchr) return
    i = j
    call chrps2(recrd,sep,len(sep),recl,i,k)
    if (i >= recl) return
    !      print *, (recrd(j0),j0=0,j)
    !      print *, (recrd(j0),j0=0,i)
    ctmp = sep(k:k)
    i = i-j
    if (i > ctlen) call fexit(-1,9,'parchv: name too long',0)
    call strcop(strn,recrd(j),i,ctmp,k)
    strn(i+1:i+1) = ' '
    call tokmat(strn,ctbl,nchr,ctlen,' ',i,k,.false.)
    ! ... Replace existing string definition
    if (i >= 0 .AND. opt1 == 1) then
       i = i+1
       ! ... Prior declaration of variable; ignore this declaration
    else if (i >= 0 .AND. opt1 == 0) then
       j = j+k
       call eostr(recrd,recl,11,' ',j)
       !        call skipbl(recrd,recl,j)
       !        call skp2bl(recrd,recl,j)
       goto 2
       ! ... Create new string definition
    else
       i = nchr+1
       ctbl(i,1) = ' '
       call strcop(ctbl(i,1),strn,ctlen,' ',k)
       nchr = i
    endif

    ! --- Copy variable contents into the table ---
    ctbl(i,2) = ' '
    j = j+k
    call skipbl(recrd,recl,j)
    if (j >= recl) return
    j0 = 1
    ctmp = recrd(j)
    if (ctmp /= '"') ctmp = ' '
    if (ctmp == '"') j=j+1
    jr = j
    call chrpos(recrd,ctmp,recl,jr)

    call pvfil1(jr,ctlen,j,recrd,cex,nchr,mxchr,ctbl,j0,ctbl(i,2),0)
    if (ctmp == '"') j=j+1
    if (opt2 == 1) call get_environment_variable(ctbl(i,2),ctbl(i,2))
    goto 2

  end subroutine parchv


end module m_rdfiln


subroutine rdfiln(unit,cch,mxlev,loop0,nlin,list,lstsiz, &
     ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
  !- Read one line of a file, with variables substitution and loop parsing
  ! ----------------------------------------------------------------
  !i Inputs
  !i   unit: file logical unit number
  !i   recl: maximum record length allowed
  !i   nr:   number of records already read.
  !i         Must set nr=0 for first call to initialize variables.
  !i   cch: comment/command character
  !i   cch(1:1) comment char, eg '#': line ignored if first char=cch(1:1)
  !i   cch(2:3) expression chars, eg '{}'.  See Remarks.
  !i   cch(4:4) directive char, eg '%', signaling rdfiln directives.
  !i   cch(5:5) if 'k', save all declared variables
  !i   cch(6,6) if 'c', try to parse for command-line char variable
  !i            declarations (-c=nam,strn nam,strn ...) when nr=0
  !i   cch(7,7) if 't', replace tabs with spaces
  !i   cch(8,8) if 'p', write contents to stdout
  !i   loop0,nlin,list,ilist,nlist,vnam: arrays to keep track of loop
  !i         nesting levels
  !i   mxlev:   maximum depth allowed for nesting.
  !i   ctbl,mxchr: table of character variables.
  !i   a: character work array of length recl.
  !o Outputs
  !o   recrd: next record input
  !o   nr is incremented by one.  Returns -nr if EOF encounterd.
  !r Remarks
  !r  *Rdfiln reads an ascii string from logical file 'unit' into character
  !r   array 'recrd', with the option of parsing algebraic expressions in
  !r   the line, as described below.  Rdfiln also supplies some facility
  !r   for branching control to skip over reading of certain lines or
  !r   repeatedly reading a block of lines.  The caller calls rdfiln until
  !r   all the lines desired are read in, or until EOF is reached.
  !r
  !r   -------------- Conceptual overview --------------
  !r   For concretness, let's assume the following
  !r     cch(1) = '#' comment character
  !r     cch(2) = '{' character marking start of expression substitution
  !r     cch(3) = '}' character marking end of expression substitution
  !r     cch(4) = '%' rfiln directive character
  !r
  !r  *Expression substitution
  !r   rdfiln parses anything inside `{...}' and substitutes the contents for
  !r   `something else'.  For example, the contents of `{...}' may by an
  !r   algebraic expression.  In that case, rdfiln evaluates the expression
  !r   numerically, and turns it back into a number, in ascii form. Thus
  !r   a string read from 'unit' as
  !r      talk {4/2} me
  !r   becomes
  !r      talk 2 me
  !r   rdfiln evaluates '4/2' as a floating-point expression using a2bin.f,
  !r   turns the result back into an ascii string, using bin2a.f, and
  !r   substitutes the resultant substring for the original.
  !r   Caution: if the substituted string is longer than recl, it is truncated.
  !r
  !r  *Variables
  !r  *Rdfiln permits three kinds of variables, floating point scalar,
  !r   floating-point vector, and character.  It uses a2bin to evaluate
  !r   floating-point expressions and bin2a to recast the result as a string.
  !r   The scalar symbols table is maintained in the standard variables table,
  !r   File symvar.f contains the source code maintaining this table.; file
  !r   symvec.f contains source code maintaining the vector variables table.
  !r   (NB: symvec allocates space for these vectors through the malloc utility;
  !r   your compiler must have pointer capability to use this table.)  The table
  !r   of character variables is maintained in the character array ctbl, which is
  !r   passed as an argument to rdfiln.
  !r
  !r   NB: when EOF is reached, rdfiln eliminates any variables declared within
  !r   the input file using `%' constructs described below (see % save directive
  !r   below for exceptions to this rule).
  !r
  !r  *Comments
  !r   Lines beginning with the comment character are purged.
  !r
  !r  *Rdfiln directives
  !r   Lines beginning with '% directive', where directive is one of:
  !r     const cconst cvar udef var vec char char0 cchar getenv
  !r     if ifdef ifndef iffile else elseif elseifd endif include includo
  !r     while repeat end
  !r     echo show stop exit save trace vfind
  !r   are interpreted by rdfiln not as part of the input, but as a directive to
  !r   do something, such as assign a value to a variable; to conditionally skip
  !r   over a block of lines; to repeatedly read a block of lines using a 'loop'
  !r   construct; and some miscellaneous commands.  Each of these is described
  !r   in the 'Rdfiln directives' section below.
  !r
  !r   -------------- Expression substitution --------------
  !r    A string in `{}', e.g. `{strn}' may contain one of the following.
  !r    These are listed in order of precedence in which rdfiln tries to
  !r    parse the contents of `{...}' .
  !r    NB: The {} can be nested.
  !r
  !r      1. the name of a character variable, say `myvar'
  !r         In that case, rdfiln replaces string `{myvar}' with contents of
  !r         `myvar'.
  !r
  !r      2. a character variable name, say `myvar', followed by a qualifier (...)
  !r         which can be one of the following:
  !r
  !r        *(integer1,integer2) --- returns a substring of `myvar'
  !r         {myvar(n1,n2)} is replaced by the (n1,n2) substring of `myvar'.
  !r
  !r        *('char-list',n) --- marks a position in contents of `myvar'
  !r         {myvar('char-list',n)} is replaced by integer, which is the
  !r         index to the n'th occurence of one element in 'char-list'.
  !r         n is optional.
  !r         Example: If myvar="foo bar", {myvar('abc',2)} evaluates to 6.
  !r
  !r        *(:e) --- returns an integer marking last nonblank character
  !r         Example: If myvar='foo bar', {myvar(:e)} evaluates to 7.
  !r
  !r        *(/'str1'/'str2'/,n1,n2) --- string substitution
  !r         {myvar(/'str1'/'str2'/,n1,n2)} substitutes str2 for str1
  !r         It does it for the n1'th to n2'th occurence.
  !r         n1 and n2 are optional, as are the quotation marks.
  !r         Example: If myvar="foo boor", {myvar(/'oo'/a/,2,2)} = "foo bar"
  !r
  !r      3. The name of a vector variable, say `myvec'
  !r         rdfiln replaces '{myvec}' with a sequence of numbers each separated
  !r         by one space, which are the contents of myvec'  Thus
  !r           % vec myvec[5] 5 4 3 2 1
  !r           {myvec}
  !r         becomes
  !r           5 4 3 2 1
  !r         (The first line declares `myvec' to be vector of length 5 and
  !r         initializes its contents; see description of % vec below)
  !r         Alternatively you can substitute a single element.  Thus
  !r           {myvec(2)}
  !r         is transformed into
  !r           4
  !r
  !r      4. a string consisting an algebraic expression of scalar numbers and
  !r         previously declared variables.  (variables are declared and set with
  !r         '%' directives; see 'Rdfiln directives' section below.)  rdfiln
  !r         parses the expression, turns the result into a string, and
  !r         substitutes the string in place of {expression}.  This is a special
  !r         case of the following:
  !r
  !r      5. A variable assignment, or a sequence of assignments separated by
  !r         commas.  This syntax returns the value of the (last) expression,
  !r         while assigning variables to evaluated expressions.
  !r         NB: the last expression need not have an assignment operator
  !r         {x=3}               ->  is replaced by '3'
  !r         {x=3,y=4}           ->  is replaced by '4'
  !r         {x=3,y=4,x*=y}      ->  is replaced by '4'
  !r         {x=3,y=4,x*=y,x*2}  ->  is replaced by '24'
  !r
  !r         The general syntax is:  {var assignment-op expr [, ... ]} .
  !r         The following are the allowed operators:
  !r         assignment-op         function
  !r           '='            simple assignment
  !r           '*='           replace 'var' by var*expr
  !r           '/='           replace 'var' by var/expr
  !r           '+='           replace 'var' by var+expr
  !r           '-='           replace 'var' by var-expr
  !r           '^-'           replace 'var' by var^expr
  !r         NB: it is permissible to omit the 'var assignment-op' pair;
  !r         may be be useful for the final expression, as in the last example.
  !r
  !r      6. A C-like syntax of the form '{?~expr~strn1~strn2}'
  !r         If expr evaluates to nonzero, the {...} is replaced by strn1
  !r         If expr evaluates to zero, the {...} is replaced by strn2
  !r         NB:  the '~' above can be any character
  !r
  !r   To summarize, emphasizing the order of precedence: rdfiln first looks to
  !r   see if the contents of `{...}' is the name of a character variable, or a
  !r   name followed by qualifiation (...).  If so, `{...}' is replaced by the
  !r   (possibly qualified) value of the variable.  If not, rdfiln sees whether
  !r   the contents of `{...}' is the name of a vector variable.  If so, rdfiln
  !r   substitutes `{...}' into a character representation of the vector as
  !r   described in step 3 above.  If this fails, rdfiln parses `{...}' as a
  !r   scalar expression, or a sequence of expressions, and `{..}' is replaced by
  !r   a character representation of the result (4 and 5 above).
  !r
  !r   Example:  suppose that the variables table looks like:
  !r     Var       Name                 Val
  !r      1        t                   1.0000
  !r      2        f                  0.00000
  !r      3        pi                  3.1416
  !r      4        a                   2.0000
  !r   ...
  !r     Vec       Name            Size   Val[1..n]
  !r      1        firstnums          5    1.0000        5.0000
  !r      2        nextnums           5    6.0000        10.000
  !r   ...
  !r       char symbol                     value
  !r      1 c                               half
  !r      2 a                               whole
  !r      3 blank
  !r
  !r   NB: The scalar variables table always begins with predefined variables
  !r   t=1,f=0 and pi.  It is STRONGLY ADVISED that you never alter any of
  !r   these variables.
  !r
  !r   You can print out the current tables of variables with the 'show'
  !r   command; see below.  (Because the vector variables can have arbitrary
  !r   length, 'show' prints only the size of the vector and the first and
  !r   last entries.  As described in more detail below, you can create such
  !r   a variables table with the following directives:
  !r
  !r   % const a=2
  !r   % char c half a whole blank " "
  !r   % vec firstnums[5] 1 2 3 4 5
  !r   % vec nextnums[5] 6 7 8 9 10
  !r
  !r   Then rdfiln substitutes for the line
  !r    {c} of the {a} {pi} is {pi/2}
  !r   yields the following:
  !r    half of the whole 3.1415926536 is 1.5707963268
  !r
  !r   whereas the line
  !r    one quarter is {1/(nextnums(4)-5)}
  !r   becomes
  !r    one quarter is .25
  !r
  !r   The following illustrates substitution of character substrings:
  !r   % char c half a whole
  !r    To {c(1,3)}ve a cave is to make a {a(2,5)}!
  !r   becomes
  !r    To halve a cave is to make a hole!
  !r
  !r   The following line illustrates substitution of vector name
  !r    {firstnums}, I caught a hare alive, {nextnums} ...
  !r   becomes
  !r    1 2 3 4 5, I caught a hare alive, 6 7 8 9 10 ...
  !r
  !r  *Nesting of {...}.  If the contents of {...} contain an
  !r   inner block of {}, the inner block is subtituted first, as
  !r   the following illustrates.  The following line
  !r      % const xx{1{2+{3+4}1}} = 2
  !r   undergoes substitution in three passes
  !r      % const xx{1{2+71}} = 2
  !r      % const xx{173} = 2
  !r      % const xx173 = 2
  !r
  !r   This line combines nesting and '{?~expr~strn1~strn2}' syntax:
  !r      MODE={?~k~B~C}3
  !r   evaluates to, if k is nonzero
  !r      MODE=B3
  !r   or, if k is zero:
  !r      MODE=C3
  !r
  !r   -------------- Rdfiln directives --------------
  !r  This section describes the syntax for each of the directives rdfiln
  !r  understands.
  !r
  !r  *'const', 'cconst' and 'var' load or alter the variables table.
  !r   A variable 'myvar' is declared eg,  % const  myvar = expr.
  !r   'expr' may be multiplied into, divided into, added into,
  !r   subtracted from or exponentiated into an already-declared variable
  !r   using one of the following C-like syntax:
  !r     myvar*=expr  myvar/=expr  myvar+=expr  myvar-=expr  myvar^=expr
  !r
  !r   'const' and 'var' are equivalent except that, for a variable
  !r   already declared, 'const' ignores a re-declaration of the
  !r   variable (nam=val), thus preserving its original value, while
  !r   'var' alters its value; see example below.
  !r
  !r   'cconst' is a conditional 'const': the first argument following
  !r   'cconst' is an expression; declarations following the expression
  !r   are parsed only if the expression evaluates to true.
  !r
  !r   'cvar' is a conditional 'var': the first argument following
  !r   'cvar' is an expression; declarations following the expression
  !r   are parsed only if the expression evaluates to true.
  !r
  !r   .... Example: for the input file
  !r     % const a = 2 b=3 c=4 d=5
  !r     a={a} b={b} c={c} d={d}
  !r     % const a=3
  !r     % var d=-1
  !r     % const b*=2 c+=3
  !r     a={a} b={b} c={c} d={d}
  !r     % cconst b==6  b+=3 c-=3
  !r     a={a} b={b} c={c} d={d}
  !r     % cconst b==6  b+=3 c-=3
  !r     a={a} b={b} c={c} d={d}
  !r
  !r   generates four lines:
  !r     a=2 b=3 c=4 d=5
  !r     a=2 b=6 c=7 d=-1
  !r     a=2 b=9 c=4 d=-1
  !r     a=2 b=9 c=4 d=-1
  !r   'a' is unchanged from its initial declaration while 'd' changes.
  !r    The two 'cconst' show that 'b' and 'c' are altered in the first
  !r    instance, since then 'b==6' is true, while are unchanged in
  !r    the second instance, since this time 'b==6' is no longer true.
  !r
  !r  *'char' and 'cchar' load or alter the character table. Directive
  !r   % char  c half     a whole      blank
  !r   loads the character table as follows:
  !r       char symbol                     value
  !r      1 c                               half
  !r      2 a                               whole
  !r      3 blank
  !r   The last value may be a blank string. 'cchar' has the syntax
  !r   % cchar nam  expr1 str1 expr2 str2 ...
  !r   expr1 expr2 etc are algebraic expressions and 'nam' takes the
  !r   value 'str1' if expr1 evaluates to true (ie nearest integer is
  !r   nonzero), the value 'str2' if expr2 evaluates to true, etc.
  !r   Re-declaration of any previously defined variable has the effect
  !r   of changing the contents of the variable
  !r  *'char0' is the same as 'char', except re-declaration of existing
  !r   variables is ignored.
  !r  *'getenv' is the same as 'char', except the string char holds
  !r   is used as a name for an environment variable, and its value
  !r   is replaced by the value of the enviroment variable.  Thus
  !r  % getenv myhome HOME
  !r   puts the string of your home directory into variable 'myhome.'
  !r
  !r  *'vec' loads or alters elements in the table of vector variables.
  !r   % vec v[n]                  creates a vector variable of length n
  !r   % vec v[n] n1 n2 n3 ...     ditto, first elements are also set
  !r   NB: Once 'v' is already declared, elements of v may be set with
  !r   the following syntax, which sets all elements bewtween i1..i2
  !r   % vec v(i) n                or
  !r   % vec v(i1:i2)  n1 n2 ... nn
  !r   There must be exactly i2-i1+1 elements n1 ... nn.
  !r   Also, if 'v' is already declared, it is an error to re-declare it.
  !r
  !r  *'vfind' finds the entry in an array that matches a specified value:
  !r   % vfind v(i1:i2)  name match-value
  !r   parses v(i) for i=i1..i2 and sets variable 'name' to i when it
  !r   finds v(i)=match-value.  If no match, 'name' is set to zero.
  !r   .... Example, the lines
  !r   % vec  a[3] 101 2002 30003
  !r   % vfind a(1:3) i 2002    <---- will set  i to 2
  !r   % vfind a(1:3) i 10      <---- will set  i to 0
  !r
  !r  *'save' preserves some variables for future use
  !r   % save              preserves all variables defined to this point
  !r   % save name [name2 ...]                saves only variables named
  !r   NB: only scalar variables may be saved.
  !r
  !r  *'udef' deletes a variable and its definition.
  !r   Only scalar and character variables may be deleted
  !r   rdfiln aborts with error if no variable exists to 'undefine'
  !r  *'udef -f' is equivalent to 'udef' except that 'udef -f' does
  !r   nothing if an attempt is made to udefine a nonexistent variable
  !r
  !r  *'trace' when turned on, chatters about how rdfiln parses the input
  !r           Invoking 'trace' with no argument toggles whether it is
  !r           on or off.
  !r           'trace 0' turns the tracing off (the default)
  !r           'trace 1' turns the tracing to lowest level:
  !r                     all directives having to do with execution flow
  !r                     (if-else-endif, repeat/while-end)
  !r           'trace 2' prints some information about most directives.
  !r
  !r  *'echo' echos the current line to stdout (i1mach(2)).
  !r  *'stop expr msg' aborts with 'msg' if 'expr' evaluates to true
  !r  *'show' prints out various things:
  !r   % show lines       (echos each line generated to the screen until:
  !r   % show stop         is encountered)
  !r   % show vars        (prints out the state of the variables table)
  !r
  !r   Expressions are evaluated for both echo and stop before printout.
  !r
  !r  *'if expr', 'elseif expr', 'else' and 'endif' are conditional read
  !r   blocks.  Lines between these directives are read or not,
  !r   depending on the value of the expression following 'if.'  For
  !r   .... Example, the lines
  !r   % if Quartz
  !r    is clear
  !r   % elseif Ag
  !r    is bright
  !r   % else
  !r    neither is right
  !r   % endif
  !r   generate one line ' is clear', if Quartz is true, ' is bright' if
  !r   Ag is false but Quartz is true, ' neither is right' otherwise.
  !r
  !r  *ifdef is similar to if, but has a more general idea of what
  !r   constitutes an expression.  First, 'if' requires a valid
  !r   expression, while 'ifdef' treats an invalid expression (eg one
  !r   containing an undefined variable) as a valid expression evaluating
  !r   to false.  The syntax of ifdef allows several expressions :
  !r   ifdef expr1 [| expr2 | expr3 ...]
  !r   and if any of expr1, expr2, ... evaluate to true, the result
  !r   is true, whether or not preceding expressions are valid.
  !r   The spaces here are syntatically significant here, since
  !r   expr1|expr2 is only true if both expr1 and expr2 are valid
  !r   expressions, while expr1 | expr2 may be true if either is valid.
  !r   'ifdef'  allows a limited use of character variables in
  !r   expressions. Either of the following are permissible expressions:
  !r     char-variable            (T if char-variable exists, otherwise F)
  !r     char-variable=='string'  (T ifchar-variable equals "string")
  !r   .... Example: ifdef  x1==2 | atom=='Mg'
  !r     is true if scalar 'x1' is 2, or if character variable
  !r     "atom" is equal to "Mg".
  !r   Also, the 'expr' above can be groups of subexpressions of any type
  !r   just mentioned, separated by ' & '.
  !r   .... Example: ifdef  x1==2 & atom=='Mg' | x1===1
  !r     is true if scalar 'x1' is 1, or 'x1' is 2, and character
  !r     variable "atom" is equal to "Mg".
  !r  *'elseifd' is to 'elseif' as 'ifdef' is to 'if'.
  !r  'if' and/or 'ifdef' constructs may be nested to a depth of mxlev.
  !r  *'ifndef' expr ... is equivalent syntatically to
  !r   'ifdef'  expr ... followed immediately by 'else'
  !r
  !r  *iffile file-name  is another conditional read block.  Instead of
  !r   the condition being set by an expression, it it set by whether
  !r   file 'file-name' exists.
  !r
  !r  *'while' and 'end', or 'repeat' and 'end' are looping constructs,
  !r   as illustrated below.  The 'while' construct has the syntax
  !r     % while [assignment assignment ...] expression
  !r      lines here are repeatly read in while expression is true
  !r     % end
  !r   here the (optional) assignments following expression have the
  !r   same syntax and meaning as the 'const' construct.  That is,
  !r   they look like 'nam=expr' or 'nam op= expr'.  As in the 'const'
  !r   case, 'nam=expr' only has effect when the variable 'nam' has
  !r   is not in the variables table.  This is made evident in the
  !r   example below.  'repeat' has the syntax
  !r     % repeat varnam list
  !r       lines here are reread, for each value 'varnam' takes in list;
  !r       list can be an integer, eg '7' or a more complex integer list,
  !r       eg '1:3,6,2' -- see mkilst.f for the syntax of an integer list.
  !r     % end
  !r   .... Example:  note in the 'while' construct the assignment
  !r        'db=-1' is operative only the first time, while 'db+=2'
  !r        is changes db in every loop.
  !r   % const nm=-3 nn=4
  !r   % while db=-1 db+=2 db<=3
  !r   % repeat k= 2,7
  !r   this is k={k} and db={db}
  !r   {db+k+nn+nm} is db + k + nn+nm, where nn+nm={nn+nm}
  !r   % end (loop over k)
  !r   % end (loop over db)
  !r   .... is expanded into
  !r   this is k=2 and db=1
  !r   4 is db + k + nn+nm, where nn+nm=1
  !r   this is k=7 and db=1
  !r   9 is db + k + nn+nm, where nn+nm=1
  !r   this is k=2 and db=3
  !r   6 is db + k + nn+nm, where nn+nm=1
  !r   this is k=7 and db=3
  !r   11 is db + k + nn+nm, where nn+nm=1
  !r
  !r  *include file-name causes rdfiln to open file 'file-name', and
  !r   input is read from the new file until EOF, after which lines are
  !r   read from the calling file.  %include may be nested to a depth
  !r   of 10.  NB:  repeat-end and if-endif constructs MUST reside in the
  !r   same file.  'includo' is identical to 'include', except that the
  !r   rdfiln aborts if the file does not exist.
  !r   Sandwiching include directives inside constructs is permissible.
  !u Updates
  !u   04 Aug 07 Bug fix when trace used together with macros
  !u   21 Dec 05 cch(8) => write contents to stdout
  !u   12 Aug 04 Comment character redefined: characters following comment
  !u             character are replaced by blanks
  !u    3 Aug 04 Changed call to nargc with call to nargf
  !u   19 Dec 02 Added macro declarations
  !u   27 Aug 01 Added option to purge tabs from input lines
  !u   13 Aug 97 {} can be nested
  !u   17 Sep 97 rdfiln allows ifdef expr & expr ...
  !u    6 Oct 97 added % vfind
  !u   24 Nov 98 added {?~expr~string1~string2}
  !u             cchar var ... doesn't alter ctbl unless an expr is true
  !u             udef char-var works
  !u             'stop' with no arguments stops
  !u             first pass at substrings in char variables
  !u   22 Oct 99 rdfiln always loads arg[0] as char variable 'progname'
  !u             added 'exit'
  !u   16 Feb 00 implemented 'iffile' and 'udef -f',
  !u             expression substitution for cchar, and
  !u             {var assignment-op expression [, ...}, {char-var(:e)}
  ! ----------------------------------------------------------------
  implicit none
  integer :: unit,recl,nr,mxchr,mxlev,lstsiz, &
       nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),nlist(0:mxlev)
  logical rdstrn,loop0(0:mxlev)
  integer :: ctlen
  parameter (ctlen=120)
  character*(*) cch, recrd(0:*)*1, a, ctbl(mxchr,2)*(ctlen)
  character vnam(mxlev)*16
  ! Local variables
  integer :: i,j,nlev,ja,jr,k,lvnam,ipr,nrd,nsyv,niff,iecho,j0,i1mach, &
       nchr,lunit,nkeepv,nelt,ival,i1,i2,nargf,ndir,ltrace,isw,nchr0, &
       mxiff
  parameter (ndir=31,mxiff=15)
  character cc1*1,cc2*1,cex*2,dir(ndir)*7,sdir(4)*5,tab*1
  !     liff  is the conditional read truth table for each nesting level
  !     liff0 is like liff, but is evaluated for each nesting level
  !     independently of the others.  It is needed to re-evaluate liff
  !     as the nesting level changes.
  logical liff(0:mxiff),liff0(0:mxiff)
  logical :: lcc,lex,a2bin,noskip,notab
  logical :: lshow,ltmp,ltmp2,lsequ,namchr,eee
  integer :: unit0(10),filstk,nrdstk(10)
  double precision :: val,xx
  character fname*120,ctmp*1,aa*512,asbst*512,strn*512
  character(2) :: aops(6)
  integer :: LIF,LFDEF,LLSEID,LLSEIF,LLSE,LNDIF,LEPEAT,LND,LONST,LAR, &
       LHAR,LCHAR,LHOW,LCHO,LTOP,LNCLUD,LAVE,LEC,LENV,LHAR0, &
       LCONST,LVAR,LDEF,LNCLUO,LHILE,LRACE,LFNDEF,LFIND,LXIT,LFILE, &
       LACRO
  parameter (LIF=0,LFDEF=1,LLSEID=2,LLSEIF=3,LLSE=4,LNDIF=5, &
       LEPEAT=6,LND=7,LONST=8,LAR=9,LHAR=10,LCHAR=11,LHOW=12,LCHO=13, &
       LTOP=14,LNCLUD=15,LAVE=16,LEC=17,LENV=18,LHAR0=19,LCONST=20, &
       LVAR=21,LDEF=22,LNCLUO=23,LHILE=24,LRACE=25,LFNDEF=26,LFIND=27, &
       LXIT=28,LFILE=29,LACRO=30)
  ! #if CSYTLETABCHR
  !      parameter (tab='\t')
  ! #else
  parameter (tab='        ')
  ! #endif
  save nchr,lshow,cc1,cc2,cex,lcc,lex,notab,nkeepv,liff,liff0, &
       nlev,nrd,niff,lunit,unit0,filstk,nrdstk,ltrace
  ! ... Allowed characters in a name
  namchr(ctmp) =ctmp .ge. 'a' .and. ctmp .le. 'z' .or. ctmp .eq. '_' &
       .or. ctmp .ge. 'A' .and. ctmp .le. 'Z' &
       .or. ctmp .ge. '0' .and. ctmp .le. '9'
  data dir /'if','ifdef','elseifd','elseif','else','endif','repeat', &
       'end','const','var','char','cchar','show','echo','stop','include', &
       'save','vec','getenv','char0','cconst','cvar','udef','includo', &
       'while','trace','ifndef','vfind','exit','iffile','macro'/
  data sdir /'stop','all','lines','vars'/ ltrace /0/
  data aops/'= ','*=','/=','+=','-=','^='/



  call getpr(ipr)
  if (recl > 512) call rx('rdfiln: increase length of aa')
  aa = ' '
  !     10 Oct 02 patch to avoid bug in Intel ifc compiler
  !     do i = 1, recl
  !       a(i:i) = ' '
  !     enddo

  ! --- if nr=0, initialization for rdfiln ---
  if (nr == 0) then
     !   ... The following guarantees basic constants loaded:
     lshow = ipr .ge. 110
     lunit = unit
     nchr = 0
     fname = '-cprogname='
     call getarg(0,fname(12:))
     i = 2
     call parchv(fname,len(fname),mxchr,cex,'=',0,nchr,ctbl,i)
     i = 2
     fname = '-cext=EXT'
     call parchv(fname,len(fname),mxchr,cex,'=',11,nchr,ctbl,i)
     j = 0
     lcc = a2bin('1 ',i,2,0,' ',j,2)
     call numsyv(nsyv)
     nkeepv = nsyv
     filstk = 0
     nrd = 0
     nlin(0) = 0
     cc1 = cch(1:1)
     lex  = .false.
     lcc  = .false.
     notab = .false.
     nlev = 0
     !   ... loop0 is true unless within a % repeat or while with null list
     loop0(nlev) = .true.
     !   ... liff: truth table for %if blocks; niff is the nesting depth
     liff(0)  = .true.
     liff0(0)  = .true.
     !   ... noskip is true unless %if evaluates to F or loop0 is F
     noskip = .true.
     niff = 0
     if (len(cch) >= 3) then
        cex = cch(2:3)
        lex = cex .ne. ' '
     endif
     if (len(cch) >= 4) then
        lcc = .true.
        cc2 = cch(4:4)
        if (cc2 == ' ') lcc = .FALSE. 
     endif
     if (len(cch) >= 5) then
        if (cch(5:5) == 'k') nkeepv = -1
     endif
     !   ... Look for char variable definitions
     if (len(cch) >= 6 .AND. cch(6:6) == 'c') then
        do  5  j = 1, nargf()-1
           call getarg(j,fname)
           if (lsequ(fname,'-c',2,' ',k)) then
              k = len(fname)
              jr = 0
              i  = 0
              call chrpos(fname,' ',k,jr)
              call chrpos(fname,'=',jr,i)
              if (i < jr) then
                 i = 2
                 call parchv(fname,k,mxchr,cex,'=',0,nchr,ctbl,i)
              endif
           endif
5       enddo
     endif
     if (len(cch) >= 7) then
        if (cch(7:7) == 't') notab = .TRUE. 
     endif
     if (len(cch) >= 8) then
        if (cch(8:8) == 'p') lshow = .TRUE. 
     endif
  endif

  !      call shosyv(0,0,0,6)
  goto 10

  ! --- Entry point for new value of loop variable ---
12 continue
  if (nlev <= 0) call rx('bug in rdfiln ... aborting')
  ilist(nlev) = ilist(nlev)+1
  if (ilist(nlev) <= nlist(nlev) .AND. liff(niff)) then
     val = list(ilist(nlev),nlev)
     call lodsyv(vnam(nlev),1,val,i)
  endif
  !      call shosyv(0,0,0,6)
  ! --- Entry point for new line read ---
10 continue
  noskip = loop0(nlev) .and. liff(niff)
  if ( .NOT. rdstrn(lunit,a,recl, .FALSE. )) goto 99
  !       12 Aug 04 Blank out remining line following comment char
  if (cc1 /= ' ') then
     ltmp = .false.
     do  i = 2, recl
        if (a(i:i) == cc1) ltmp = .TRUE. 
        if (ltmp) a(i:i) = ' '
     enddo
  endif
  if (notab) then
     do  i = 1, recl
        if (a(i:i) == tab) a(i:i) = ' '
     enddo
  endif
  !       call strncp(aa,a,1,1,recl)
  !       Keep track of line number
  iecho = 0
  nrd = nrd+1
  do   j = 0, nlev
     nlin(j) = nlin(j)+1
  enddo
  !       print *, 'nlev=',nlev,' nlin=',(nlin(j),j=1,nlev),
  !     .    '  a=',(a(j),j=1,recl)
  j0 = 0

  !   --- Directives to rdfiln ---
  if (lcc .AND. a(1:1) == cc2) then
     j = 1
     call skipbl(a,recl,j)
     k = j-1
     call tokmat(a(j+1:),dir,ndir,7,' ',i,j,.false.)
     if (ltrace > 0 .AND. i >= 0) then
        aa = ' '
        call awrit1('#rf %i: '''//dir(i+1),aa,80,0,nrd)
     endif

     if (i == LIF .OR. i == LFDEF .OR. i == LFNDEF .OR. &
          i == LFILE) then
        niff = niff+1
        if (niff > mxiff) call rx('increase niff in rdfiln')
        j = j+k
        liff0(niff) = .false.
        liff(niff) = .false.
        if ( .NOT. noskip .AND. ltrace > 2) goto 81
        if ( .NOT. noskip) goto 10
        !       ... make substitutions for eg % ... {}...
        ja = 1
        asbst = ' '
        jr = 0
        call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja, &
             asbst,ltrace)
        !       ... Determine whether line results in T or F
        if (i == LFILE) then
           i1 = j+1
           call nword(asbst,1,i1,i2)
           liff0(niff) = .false.
           if (i2 >= i1) then
              inquire(file=asbst(i1:i2),exist=eee) !2022apr24 right?
              liff0(niff) = eee !i .eq. 1
           endif
        else
           call rdfilx(asbst,ja,i.eq.LIF,j,nrd,ctbl,mxchr,nchr, &
                liff0(niff))
           if (i == LFNDEF) liff0(niff) = .NOT. liff0(niff)
        endif
        liff(niff) = liff(niff-1) .and. liff0(niff)
        goto 81
     elseif (i == LLSEID .OR. i == LLSEIF .OR. i == LLSE) then
        if (niff <= 0) call fexit(-1,1,' Exit -1 rdfiln:'// &
             ' else encountered with matching if, line %i',nrd)
        if ((i == LLSEID .OR. i == LLSEIF) .AND. &
             .NOT. liff0(niff)) then
           j = j+k
           !         ... make substitutions for eg % ... {}...
           ja = 1
           asbst = ' '
           jr = 0
           call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja, &
                asbst,ltrace)
           call rdfilx(asbst,ja,i.eq.LLSEIF,j,nrd,ctbl,mxchr,nchr, &
                liff0(niff))
           liff(niff) = liff(niff-1) .and. liff0(niff)
        else
           liff(niff) = liff(niff-1) .and. .not. liff0(niff)
        endif
        goto 81
     elseif (i == LNDIF) then
        niff = niff-1
        if (niff < 0) call rx &
             ('rdfiln:  endif  encountered before  if')
        goto 81
     elseif (i == LHILE .AND. liff(niff)) then
        nlev = nlev+1
        if (nlev > mxlev) call fexit(-1,1,' Exit -1 RDFILN:  '// &
             'loop directives nested deeper than %i',mxlev)
        !       ... Make %end backspace to this %while line
        nlin(nlev) = 1
        !       ... Flag the matching %end that it matches a %while construct
        nlist(nlev) = -1
        ilist(nlev) = 0
        !       ... skip while block if in a nested null-list
        loop0(nlev) = .false.
        if ( .NOT. loop0(nlev-1)) goto 82
        j = j+k
        !       ... make variable substitutions
        ja = 1
        asbst = ' '
        jr = 0
        call strncp(asbst,a,1,1,recl)
        call pvfil1(recl,recl,jr,asbst,cex,nchr,mxchr,ctbl,ja, &
             a,ltrace)
        !       ... Come back here until next string not an assignment
13      call skipbl(a,recl,j)
        if (a(j+1:j+1) /= cch(2:2)) then
           j0 = j+1
           !         ... Return here until char not one belonging to 'name'
23         j0 = j0+1
           if (namchr(a(j0:j0))) goto 23
           !         ... An assignment?  If so, returns i>0, or else '='
           call tokmat(a(j0:),aops,6,2,' ',i,k,.false.)
           if (i > 0 .OR. a(j0:j0) == '=' .AND. a(j0+1:j0+1) /= '=') then
              call parsyv(a,recl,1,0,j)
              call strncp(aa,a,1,1,recl)
              goto 13
           endif
        endif
        !       ... Done with assignments, next follows expression
        ltmp = .false.
        call rdfilx(a,recl,.false.,j,nrd,ctbl,mxchr,nchr,ltmp)
        loop0(nlev) = ltmp .and. noskip
        goto 82
     elseif (i == LEPEAT .AND. liff(niff)) then
        nlev = nlev+1
        if (nlev > mxlev) call fexit(-1,1,' Exit -1 RDFILN:  '// &
             'loop directives nested deeper than %i',mxlev)
        ilist(nlev) = 0
        nlin(nlev) = 0
        nlist(nlev) = 0
        loop0(nlev) = .false.
        if ( .NOT. loop0(nlev-1)) goto 82
        j = j+k
        !       ... make variable substitutions
        ja = 1
        asbst = ' '
        jr = 0
        call strncp(asbst,a,1,1,recl)
        call pvfil1(recl,recl,jr,asbst,cex,nchr,mxchr,ctbl,ja, &
             a,ltrace)
        !           Pick up variable name
        call skipbl(a,recl,j)
        vnam(nlev) = ' '
        call strcop(vnam(nlev),a(j+1:),16,'=',lvnam)
        vnam(nlev)(lvnam:lvnam) = ' '
        call strncp(strn,a,1,1,recl)
        call mkilst(strn(j+1+lvnam:),nlist(nlev),list(1,nlev))
        if (nlist(nlev) > lstsiz) &
             call fexit2(-1,1,'Exit -1 rdfiln %i: list exceeded'// &
             ' maximum number of values allowed (%i)',nrd,lstsiz)
        if (nlist(nlev) < 0 .AND. ipr >= 40) call shosyv(0,0,0,6)
        if (nlist(nlev) < 0) call rx('rdfiln: bad or null list')
        loop0(nlev) = nlist(nlev) .gt. 0 .and. loop0(nlev-1)
        lvnam = lvnam-1
        if (ltrace > 0) then
           write (fname,348) nrd, vnam(nlev)(1:lvnam)
348        format(' rdfiln', i5, ': ', a,' -> ')
           call awrit2(cc1//'%a%n:1i',fname,len(fname), &
                -i1mach(4),nlist(nlev),list(1,nlev))
        endif
        if (ltrace > 0) &
             call awrit4('%a'' over %i values; read following lines: ' &
             //'%?#n#yes#no#%?#n#%2b(blocked by level %i)',aa,80, &
             -i1mach(2),nlist(nlev),isw(loop0(nlev)), &
             isw( .NOT. loop0(nlev-1) .AND. nlev > 1),nlev-1)
        goto 12
     elseif (i == LND .AND. liff(niff)) then
        !           print *, 'nlev=',nlev,'  ','nlist=',(nlist(i),i=1,nlev)
        !           print *, 'vnam now ', vnam(nlev)
        if (nlev <= 0) call fexit(-1,1,' rdfiln %i (stop):'// &
             ' end encountered before while or repeat',nrd)
        !       ... if incomplete a repeat:end sequence or a while-end
        if (ilist(nlev) < nlist(nlev) .OR. &
             nlist(nlev) == -1 .AND. loop0(nlev)) then
           !             print *, 'nlev=',nlev,' backspacing by',nlin(nlev)
           k = nlin(nlev)
           do  i = 1, nlin(nlev)
              backspace(lunit)
           enddo
           nrd = nrd-nlin(nlev)
           call rxx(nrd.lt.0,'rdfiln: attempt to loop across files')
           do  i = 0, nlev
              nlin(i) = nlin(i)-nlin(nlev)
           enddo
           !         ... Case while-end
           if (nlist(nlev) == -1) then
              nlev = nlev-1
              goto 82
              !         ... Case repeat-end
           else
              if (ltrace > 0) call awrit2( &
                   '%a'', repeat loop over '''//vnam(nlev)// &
                   '%a'' with val=%i; reread %i lines', &
                   aa,80,-i1mach(2),list(ilist(nlev),nlev),k)
              goto 12
           endif
        endif
        !       ... Case loop has finished
        !            nlist(nlev) = 0
        !            loop0(nlev) = .true.
        nlev = nlev-1
        goto 82
     elseif ((i == LCONST .OR. i == LVAR) .AND. noskip) then
        call numsyv(nsyv)
        j = j+k
        ltmp = noskip
        call rdfilx(a,recl,.false.,j,nrd,ctbl,mxchr,nchr,ltmp)
        if (ltmp) then
           !         ... make substitutions for eg % const a=nam{xx}
           ja = 1
           aa = ' '
           call pvfil1(recl,len(aa),j,a,cex,nchr,mxchr,ctbl,ja,aa, &
                ltrace)
           j = 0
           call parsyv(aa,ja+1,999,i-LCONST,j)
        endif
        !           call shosyv(0,0,0,6)
        if (noskip) goto 83
        goto 10
     elseif ((i == LONST .OR. i == LAR) .AND. noskip) then
        call numsyv(nsyv)
        j = j+k
        !       ... make substitutions for eg % const a=nam{xx}
        ja = 1
        aa = ' '
        call pvfil1(recl,len(aa),j,a,cex,nchr,mxchr,ctbl,ja,aa, &
             ltrace)
        j = 0
        call parsyv(aa,ja+1,999,i-LONST,j)
        !            call shosyv(0,0,0,6)
        !            pause
        if (noskip) goto 83
        goto 10
     elseif ((i == LEC .OR. i == LFIND) .AND. noskip) then
        call numsvv(nsyv)
        j = j+k
        !       ... make substitutions for eg % vec a{xx}...
        ja = 1
        aa = ' '
        jr = 0
        call pvfil1(recl,len(aa),jr,a,cex,nchr,mxchr,ctbl,ja,aa, &
             ltrace)
        call strncp(a,aa,1,1,min(ja,recl))
        call skipbl(a,recl,j)
        j0 = j
        call chrps2(a,'([',2,recl,j,k)
        if (k < 1) goto 999
        aa = ' '
        call strcop(aa,a(j0+1:),j-j0,' ',i1)
        !       ... does vector variable aa already exist?
        call getsvv(aa,ival,0,1,1,val)
        !       ... case not already declared:
        if (ival == 0) then
           if (k /= 2) call fexit(-1,1,' Exit -1: rdfiln: attempt' &
                //' to use undeclared vector, '//'line %i',nrd)
           j = j+1
           if ( .NOT. a2bin(a,nelt,2,0,']',j,recl)) goto 999
           call addsvv(aa,nelt,ival)
           call parsvv(a,recl,ival,nelt,1,j)
           !       ... case already declared:
        else
           if (k /= 1) then
              call awrit1('#rf (warning): attempt' &
                   //' to redeclare variable '//aa//'%a, line %i', &
                   ' ',80,i1mach(2),nrd)
              goto 91
              goto 10
           endif
           !       ...   Find indices for vector subblock
           j = j+1
           j0 = j
           ltmp2 = a2bin(a,i1,2,0,':',j0,recl)
           if (ltmp2) j = j0
           if ( .NOT. a2bin(a,nelt,2,0,')',j,recl)) goto 999
           if ( .NOT. ltmp2) i1 = nelt
           if (i == LFIND) then
              call strncp(aa,a,1,1,recl)
              call skipbl(a,recl,j)
              fname = ' '
              call strcop(fname,a(j+1:),ctlen,' ',k)
              j = j+k
              if ( .NOT. a2bin(a,val,4,0,' ',j,recl)) goto 999
              k = 0
              do  45  i = i1, nelt
                 call getsvv(' ',ival,1,i,i,xx)
                 if (xx == val) then
                    k = i
                    goto 46
                 endif
45            enddo
46            continue
              call lodsyv(fname,1,dble(k),k)
              goto 91
           else
              call parsvv(a,recl,ival,nelt-i1+1,i1,j)
           endif
        endif
        if (noskip) goto 91
        goto 10
     elseif ((i == LHAR .OR. i == LHAR0 .OR. i == LENV) &
          .AND. noskip) then
        nchr0 = nchr
        j = j+k
        k = 1
        if (i == LHAR0) k = 0
        if (i == LENV) k = k+10

        !       ... make substitutions for eg % ... {}...
        ja = 1
        asbst = ' '
        jr = 0
        call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja, &
             asbst,ltrace)
        call parchv(asbst,ja,mxchr,cex,'= ',k,nchr,ctbl,j)
        if (noskip) goto 85
        goto 10
     elseif (i == LCHAR .AND. noskip) then
        nchr0 = nchr
        j0 = k+j
        i1 = i

        ja = 1
        asbst = ' '
        jr = 0
        call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja, &
             asbst,ltrace)

        call skipbl(asbst,recl,j0)
        if (j0 >= recl) goto 10
        call tokmat(asbst(j0+1:),ctbl,nchr,ctlen,' ',i1,k,.false.)
        ! ... Replace existing string definition
        if (i1 >= 0) then
           i1 = i1+1
           ! ... Create new string definition
        else
           i1 = nchr+1
           ctbl(i1,1) = ' '
           call strcop(ctbl(i1,1),asbst(j0+1:),ctlen,' ',k)
           nchr = i1
           if (nchr > mxchr) &
                call rx('rdfiln: too many char variables')
        endif
        j0 = j0+k
        ! ...       Loop until expression evaluates to true
73      call skipbl(asbst,recl,j0)
        ltmp2 = .true.
        ltmp = noskip
        j = j0
        call rdfilx(asbst,recl,.true.,j0,nrd,ctbl,mxchr,nchr,ltmp)
        !           print *, aa(1:j0)
        !           Expression evaluated to false; skip over associated expr.
        if ( .NOT. ltmp) then
           !             If at end, give up
           if (j0 >= recl) goto 87
           !             This ensures movement past end of expression
           call chrpos(asbst,' ',recl,j0)
           call skipbl(asbst,recl,j0)
           !             (chrpos works except doesn't allow spaces in string)
           !             call chrpos(asbst,' ',recl,j0)
           !             If at end, give up
           if (j0 >= recl) goto 87
           ctmp = asbst(j0+1:j0+1)
           if (ctmp /= '"') ctmp = ' '
           if (ctmp == '"') j0=j0+1
           call chrpos(asbst,ctmp,recl,j0)
           if (ctmp == '"') j0=j0+1
           !             print *, (asbst(m), m=1,j0), '|'
           !             If at end, give up
           if (j0 >= recl) goto 87
           goto 73
        endif
        call skipbl(asbst,recl,j0)
        !           Expression evaluates to true, but no string: set to blank
        if (j0 == recl) then
           ctbl(i1,2) = ' '
           if (noskip) goto 85
           goto 10
        endif
        !           print *, (asbst(m), m=1,j0), '|'
        !           call strcop(ctbl(i1,2),asbst(j0+1),ctlen,' ',k)
        k = 1
        ctmp = asbst(j0+1:j0+1)
        if (ctmp /= '"') ctmp = ' '
        if (ctmp == '"') j0=j0+1
        ja = j0
        call chrpos(asbst,ctmp,recl,ja)
        aa = ' '
        call pvfil1(ja,ja,j0,asbst,cex,nchr,mxchr,ctbl,k,aa,ltrace)
        ctbl(i1,2) = aa
        if (ctmp == '"') ja=ja+1
        if (noskip) goto 85
        goto 10
     elseif (i == LHOW .AND. noskip) then
        j = j+k
        call skipbl(a,recl,j)
        i = 2
        if (j >= recl) goto 32
        call tokmat(a(j+1:),sdir,4,5,' ',i,j,.false.)
        goto (31,32,32,34), i+1
        call fexit(-1,1,' Exit -1: rdfiln failed to parse show, '// &
             'line %i',nrd)
        !      ...  'show stop'
31      continue
        lshow = ipr .gt. 110
        goto 10
        !      ...  'show all' or 'show lines'
32      continue
        lshow = .true.
        if (i == 2) goto 10
        !      ...  'show vars'
34      continue
        call awrit1('#rf %i:',' ',80,i1mach(2),nrd)
        call shosyv(0,0,0,6)
        call numsvv(ival)
        if (ival > 0) then
           print *, '----'
           call shosvv(1,ival,i1mach(2))
        endif
        if (nchr > 0) print *, '---- character variables:'
        do  35  k = 1, nchr
           print '(i4,2x,a,2x,a)', k, ctbl(k,1)(1:20),ctbl(k,2)
           ! ccccccccccccc
           !           call awrit1('%,4i  '//ctbl(k,1)(1:20)//ctbl(k,2)//'%a',
           !     .          ' ',ctlen+25,i1mach(2),k)
           ! ccccccccccccc
35      enddo
        goto 10
     elseif (i == LRACE) then
        j0 = j+k
        ltmp = noskip
        if ( .NOT. ltmp) goto 10
        !       ... No expression, toggle trace and exit
        call skipbl(a,recl,j0)
        if (j0 >= recl) then
           if (ltrace <= 0) then
              ltrace = 1
           else
              ltrace = 0
           endif
        else
           j = j0
           call rdfilx(a,recl,.true.,j0,nrd,ctbl,mxchr,nchr,ltmp)
           ltrace = 0
           if (ltmp) ltrace = 1
           !         ... If a valid expression, set to numerical value
           if (ltmp) then
              ltmp2 = a2bin(a,k,2,0,' ',j,recl)
              if (ltmp2) ltrace = k
           endif
        endif
        call awrit3(' rdfiln %i: trace %?#n>0#set to %i'// &
             '#turned off#',' ',80,i1mach(2),nrd,ltrace,ltrace)
        goto 10
     elseif (i == LCHO .AND. noskip) then
        iecho = 1
        j0 = j+k
     elseif ((i == LXIT .OR. i == LTOP) .AND. noskip) then
        j0 = j+k
        call skipbl(a,recl,j0)
        if (j0 >= recl) then
           if (i == LXIT) goto 99
           call awrit1('#rf %i: stop encountered', &
                ' ',80,i1mach(2),nrd)
           call cexit(-1,1)
        endif
        ltmp = a2bin(a,ltmp2,0,0,' ',j0,recl)
        j0 = j0-1
        iecho = 4
        if (ltmp .AND. ltmp2) iecho = 3
     elseif ((i == LNCLUD .OR. i == LNCLUO) .AND. noskip) then
        filstk = filstk+1
        nrdstk(filstk) = nrd
        unit0(filstk) = lunit
        lunit = 99-filstk
        nrd = 0
        j = j+k
        call skipbl(a,recl,j)
        ja = 1
        fname = ' '
        call pvfil1(recl,len(fname),j,a,cex,nchr,mxchr,ctbl,ja, &
             fname,ltrace)
        call skpblb(fname,len(fname),k)
        if (ltrace > 0) &
             call awrit1('%x#rf %i: '''//dir(i+1)//'%a'' opening'// &
             ' file '''//fname//'%a''',aa,120,i1mach(2),nrdstk(filstk))
        !       ... Should we use fopna?
        if (i == LNCLUD) then
           open(lunit,FILE=fname(1:k+1),STATUS='UNKNOWN',ERR=98)
        elseif (i == LNCLUO) then
           open(lunit,FILE=fname(1:k+1),STATUS='OLD',ERR=98)
        endif
        !       ... Read new line with different logical unit
        goto 10
        !       ... Error if cannot open file
98      call rx('rdfiln: error opening include file "'// &
             fname(1:k+1)//'"')
     elseif (i == LDEF .AND. noskip) then
        !       ... For each name, swap with top and remove top variable
        j0 = k+j

        ja = 1
        asbst = ' '
        jr = 0
        call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja, &
             asbst,ltrace)

        call skipbl(asbst,recl,j0)
        !           set ltmp = .false. unless udef -f , then ltmp = .true.
        ltmp = .false.
        if (j0 < recl-3) then
           if (asbst(j0+1:j0+3) == '-f ') then
              j0 = j0+2
              ltmp = .true.
           endif
        endif
43      continue
        call skipbl(asbst,recl,j0)
        if (j0 >= recl) goto 10
        fname = ' '
        k = 0
        call strcop(fname,asbst(j0+1:),ctlen,' ',k)
        j0 = j0+k
        !       ... First try to udef a char variable
        call tokmat(fname,ctbl,nchr,ctlen,' ',i1,k,.false.)
        i1 = i1+1
        if (i1 > 0) then
           do  k = i1+1, nchr
              ctbl(k-1,1) = ctbl(k,1)
              ctbl(k-1,2) = ctbl(k,2)
           enddo
           nchr = nchr-1
           goto 43
        endif
        !       ... Otherwise, try to udef a scalar variable
        call getsyv(fname,val,k)
        if (k == 0 .AND. ltmp) goto 43
        if (k == 0) call fexit(-1,1,' Exit -1 rdfiln line %i:'// &
             '  no variable "'//fname//'%a" to undefine',nrd)
        call numsyv(nsyv)
        call togsyv(k,nsyv)
        nsyv = nsyv-1
        call clrsyv(nsyv)
        goto 43

     elseif (i == LAVE .AND. noskip) then
        !       ... If keeping all anyway, forget this command
        if (nkeepv < 0) goto 10
        j0 = k+j
        !       ... if no arguments, save everything up to this point
        call skipbl(a,recl,j0)
        if (j0 >= recl) then
           call numsyv(nkeepv)
           goto 10
        endif
        !       ... For each name, swap with name at nkeepv
42      continue
        call skipbl(a,recl,j0)
        if (j0 >= recl) goto 10
        fname = ' '
        k = 0
        call strcop(fname,a(j0+1:),ctlen,' ',k)
        j0 = j0+k
        call getsyv(fname,val,k)
        if (k == 0) call fexit(-1,1,' Exit -1 rdfiln line %i:'// &
             '  no variable "'//fname//'%a" to save',nrd)
        if (k <= nkeepv) goto 42
        nkeepv = nkeepv+1
        call togsyv(k,nkeepv)
        goto 42
     elseif (i == LACRO .AND. noskip) then
        call strncp(aa,a,1,1,recl)
        j = j+k
        call macset(aa(j+1:),k)
        if (k < 0) &
             call rxs('rdfiln failed to parse macro',aa(j+1:))
        goto 10
     endif
  endif

  !   --- Copy this line into recrd, substituting variable names ---
  if (a(1:1) /= cc1 .AND. noskip) then
     call strncp(recrd,a,1,1,recl)
     if (lex) then
        jr = j0
        ja = 1
        call pvfil1(recl,recl,jr,recrd,cex,nchr,mxchr,ctbl,ja,a, &
             ltrace)
        do   i = ja, recl
           a(i:i) = ' '
        enddo
        if (iecho /= 0) then
           if (mod(iecho,2) /= 0) then
              j = 1
              call skipbl(a,recl,j)
              call skpblb(a,recl,k)
              write(i1mach(4),346) nrd, (a(jr+1:jr+1), jr=j,k)
346           format('#rf', i5, ': ', 256a1)
              !               call cwrite(a,j,k,1)
           endif
           if (iecho/2 == 1) then
              call awrit1('#rf %i: stop encountered', &
                   ' ',80,i1mach(2),nrd)
              call cexit(-1,1)
           endif
           goto 10
        else
           call strncp(recrd,a,1,1,recl)
        endif
     endif
     nr = nr+1
     if (lshow) then
        call skpblb(recrd,recl,k)
        print '(256a1)', (recrd(jr), jr=0,k)
     endif
     return
  endif
  goto 10

99 continue
  if (filstk > 0) then
     if (ltrace > 0) call awrit1('#rf closing include file: ' &
          //'%i lines',' ',80,i1mach(2),nrd)
     close(lunit)
     do  j = 0, nlev
        nlin(j) = nlin(j)-nrd
     enddo
     lunit = unit0(filstk)
     nrd = nrdstk(filstk)
     filstk = filstk-1
     ! ...   And on our merry way, with the old unit
     goto 10
  endif
  if (nkeepv /= -1) call clrsyv(nkeepv)
  nr = -nr
  if (ltrace > 0 .OR. ipr <= 100) return
  call awrit1(cc1//'rdfiln generated %i records.', &
       fname,len(fname),i1mach(2),-nr)
  return

  ! --- Trace entry points ---
  ! ... if-else-endif constructs
81 if (ltrace == 0) goto 10
  call awrit5('%a''%?#n>1# (nesting=%i)#%j#, read following lines: ' &
       //'%?#n#yes#no#%?#n#%2b(blocked by level %i)', &
       aa,80,-i1mach(2),niff,niff,isw(liff(niff)), &
       isw(.not.liff(max(niff,1)-1).and.niff.gt.1),niff-1)
  goto 10
  ! ... while/repeat-end constructs
82 if (ltrace == 0) goto 10
  call awrit5('%a''%?#n>1# (nesting=%i)#%j#, read following lines: ' &
       //'%?#n#yes#no#%?#n#%2b(blocked by level %i)', &
       aa,80,-i1mach(2),nlev,nlev,isw(loop0(nlev)), &
       isw(.not.loop0(max(nlev,1)-1).and.nlev.gt.1),nlev-1)
  goto 10
  ! ... const,var directives
83 if (ltrace > 1) then
     call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
     call numsyv(j)
     do   i = nsyv+1, j
        call watsyv(fname,val,i)
        call awrit1('%a '//fname//'%a%?#n<0#,##',aa,80,0,i-j)
     enddo
     if (nsyv == j) call awrit0('%a ... no new vars',aa,80,0)
     call awrit0('%a',aa,80,-i1mach(2))
  endif
  goto 10
  ! ... char directives
85 if (ltrace > 1) then
     call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
     do  i = nchr0+1, nchr
        call awrit1('%a '//ctbl(i,1)//'%a%?#n<0#,##',aa,80,0,i-nchr)
     enddo
     if (nchr0 == nchr) call awrit0('%a ... no new vars',aa,80,0)
     call awrit0('%a',aa,80,-i1mach(2))
  endif
  goto 10
  ! ... cchar read failed
87 if (ltrace > 1) then
     call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', var '//ctbl(nchr,1)// &
          '%a not changed',aa,80,-i1mach(2),nrd)
  endif
  nchr = nchr0
  goto 10
  ! ... vec directives
91 if (ltrace > 1) then
     call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
     call numsvv(j)
     do  i = nsyv+1, j
        call watsvv(fname,i)
        call awrit1('%a '//fname//'%a%?#n<0#,##',aa,80,0,i-j)
     enddo
     if (nsyv == j) call awrit0('%a ... no new vars',aa,80,0)
     call awrit0('%a',aa,80,-i1mach(2))
  endif
  goto 10
  ! --- General error exit ---
999 call fexit(-1,1,' Exit -1: rdfiln: parse failed, line %i',nrd)
end subroutine rdfiln

subroutine rdfile(unit,cch,recrd,mxrecs,a,recl,nr)
  !- Reads entire disk file into recrd, calling rdfiln.
  implicit none
  integer :: unit,mxrecs,recl,nr
  integer :: mxchr,mxlev,lstsiz,ctlen
  parameter (mxchr=40,mxlev=6,lstsiz=200,ctlen=120)
  character*(*) cch, recrd(0:*)*1, a(*)*1
  character ctbl(mxchr,2)*(ctlen)
  logical loop0(0:mxlev)
  integer :: nlin(0:mxlev),list(lstsiz,mxlev), &
       ilist(mxlev),nlist(0:mxlev)
  character vnam(mxlev)*16
  nr = 0
10 continue
  call rdfiln(unit,cch,mxlev,loop0,nlin,list,lstsiz, &
       ilist,nlist,vnam,ctbl,mxchr,a,recrd(nr*recl),recl,nr)
  if (nr <= 0) then
     nr = -nr
     return
  endif
  if (nr == mxrecs) return
  goto 10
end subroutine rdfile

subroutine rdfilx(a,recl,lif,j,nrd,ctbl,mxchr,nchr,liff0)
  !- Determines whether an 'if' or 'ifdef' expression is true
  !  Inputs
  !    a,recl
  !    lif    T if 'if', F if to parse for 'ifdef'
  !    j      starting position in a
  !    nrd    line number (no longer used)
  !    ctbl,mxchr,nchr character table, leading dimension and number
  !    liff0  F if in a false branch of a conditional-read, else T
  !  Outputs
  !    liff0  Set to F if no T expression found; see Remarks.
  !    j      last character searched
  !  Remarks
  !r  Syntax for 'ifdef': ifdef expr [ op expr ... ]
  !r  Here op is the logical logical OR '|' or the logical AND '&'.
  !r  The entire collection expr op expr op expr ... can be partitioned
  !r  into the following:
  !r    superexpr | superexpr | superexpr ...
  !r  where
  !r    superexpr has the form expr [& expr ...]
  !r
  !r  A superexpr is true if and only if each of its separate expressions
  !r  is true.  The entire
  !r    superexpr | superexpr | superexpr ...
  !r  is true if any one of the the superexpr is true, whether or not
  !r  even they are valid expressions.
  !r  A single 'expr' is true if any of the following are true:
  !r    1.  it is a numerical expression that evaluates to true
  !r    2.  It is an invalid numerical expression, but expr is
  !r        the name of a character variable
  !r    3.  It is of the form    var=='string'
  !r        where var is the name of a character variable whose contents
  !r        are 'string'
  ! ----------------------------------------------------------------
  implicit none
  logical :: lif,liff0
  integer :: recl,j,nrd,nchr,mxchr,ctlen
  parameter (ctlen=120)
  character*(*) ctbl(mxchr,2)*(ctlen),aa*(ctlen)
  character*(*):: a
  ! Local variables
  logical :: a2bin,ltmp,liffa
  integer :: j0,i0,i1,k
  liffa = .true.
10 continue
  call skipbl(a,recl,j)
  j0 = j
  ! ... Attempt to read arithmetic expression
  ltmp = a2bin(a,liff0,0,0,' ',j,recl)
  if (ltmp) j = j-1
  ! ... If invalid expression, liff0 is assumed false
  liff0 =  liff0 .and. ltmp
  ! ... Case invalid-numerical-expr; check for character expression
  if ( .NOT. ltmp) then
     j = j0
     call chrps2(a,' =',2,recl,j,k)
     aa = a(j0+1:j)
     call tokmat(aa,ctbl,nchr,ctlen,' ',i0,i1,.false.)
     ! ...   Check for a character expression
     if (i0 >= 0) then
        !     ... if a ' ' was encountered before an '='
        if (k == 1) then
           liff0 = .true.
           !     ... if an '=' was encountered before a ' '
        else if (k == 2) then
           if (a(j+1:j+3) == '==''') then
              i1 = j+4
              call chrpos(a,'''',recl,i1)
              aa = a(j+4:i1)
              liff0 = ctbl(i0+1,2) .eq. aa
           endif
        endif
        call skp2bl(a,recl,j)
     else
        call skp2bl(a,recl,j)
     endif
  endif
  ! ... Now liff0 is result of current expression.
  !     Nothing more to do for 'if' statement
  if (lif) return
  ! ... logical liff0 with a prior 'and'
  ltmp = liff0
  liff0 = liff0 .and. liffa
  ! ... Case 'and' following this expression
  call skipbl(a,recl,j)
  j = j+1
  if (j >= recl) return
  if (a(j-1:j+1) == ' & ') then
     liffa = liffa .and. liff0
     goto 10
     ! ... Case 'or' following this expression
  elseif (a(j-1:j+1) == ' | ') then
     !   ... The end of possible expr & expr & ...
     if (liff0) return
     liffa = .true.
     goto 10
  endif

  ! ... Put j at last nw char, equivalent to pos of a2bin above
  call skpblb(a,j-1,j)
  j = j+1
end subroutine rdfilx

subroutine pvfil1(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,ja,a,t)
  !- Kernel of rdfiln that parses string with substitutions
  ! ----------------------------------------------------------------
  !i Inputs
  !i   reclr   length of recrd
  !i   recla   length of a
  !i   ja,jr   starting positions of a,recrd
  !i   recrd   string that is to be subsitituted
  !i   ctbl,nchr,mxchr table of character variables, number and dim.
  !i   t       trace, print out intermediate substitutions if t>=3
  !o Outputs
  !o   ja,jr   as final positions of a,recrd
  !o   a       substituted string assembled into a.
  !r Remarks
  !r  Note recrd(0..), a(1..)!
  ! ----------------------------------------------------------------
  implicit none
  ! Passed Parameters
  integer :: reclr,recla,nchr,ja,jr,t,mxchr
  integer :: ctlen
  parameter (ctlen=120)
  character(1) :: recrd(0:reclr-1),a(recla),ctbl(mxchr,2)*(ctlen)
  ! Local variables
  logical :: pvfil2
  character(512) :: aa
  integer :: ia,ja0,jaa
  character:: cex*2

  do  ia = ja, recla
     a(ia) = ' '
  enddo
  ja0 = ja
  if ( .NOT. pvfil2(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,ja,a)) &
       return
10 continue
  aa = ' '
  do   ia = ja0, ja
     aa(ia:ia) = a(ia)
  enddo
  if (t > 2) then
     call skpblb(aa,ja,ja)
     print 333, '#rf subst: "', aa(1:ja+1),'"'
333  format(a,a,a)
  endif
  jaa = ja0-1
  ja = ja0
  if (pvfil2(recla,recla,jaa,aa,cex,nchr,mxchr,ctbl,ja,a)) goto 10
end subroutine pvfil1

logical function pvfil2(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl, &
     ja,a)
  !- Kernel of rdfiln that parses string with substitutions
  !  and returns T if substitution was inside a nested block.
  ! ----------------------------------------------------------------
  !i Inputs
  !i   reclr   length of recrd
  !i   recla   length of a
  !i   ja,jr   starting positions of a,recrd
  !i   recrd   string that is to be subsitituted
  !i   cex     characters demarcating string substitution, eg '{}'
  !i   ctbl,nchr,mxchr table of character variables, number and dim.
  !o Outputs
  !o   ja,jr   as final positions of a,recrd
  !o   a       substituted string assembled into a.
  !r Remarks
  !r  Note recrd(0..), a(1..)!
  ! ----------------------------------------------------------------
  implicit none
  integer :: mxchr,reclr,recla,nchr,ja,jr
  integer,parameter::ctlen=120
  character(1) :: recrd(0:reclr-1),a(recla),ctbl(mxchr,2)*(ctlen)
  ! Local variables
  logical :: a2bin,a2bina,lqchr
  character(1) :: qchr,match
  character:: cex*2,aa*512,aa2*(ctlen),aastr*(ctlen),aa3*(ctlen),aasub*(ctlen)

  ! #if CRAY
  !      parameter (qchr='\')
  ! #else
  parameter (qchr='\\')
  ! #endif
  integer :: i,j,k,m,ival,kkk,i1,i2,ip1,ip2, &
       ip,ii(2),ix(2),istr1,istr2,a2vec,iaa21,iaa22,iaa2,iaa3
  equivalence (i1,ii(1)), (i2,ii(2))
  double precision :: res

  pvfil2 = .false.
  if (cex(1:1) == ' ') return
10 continue
  if (jr >= reclr) return
  !   --- Parse for next expression; quote \cex literally ---
  k = -1
  k = k+1
  call chrpos(recrd(jr),cex,reclr-jr,k)
  !       lqchr is T if quote character literally
  lqchr = k.ge.1 .and. recrd(jr+k-1) .eq. qchr
  !   --- Copy string up to substitution ---
  k = min(k,recla-ja,reclr-jr-1)
  if (k < 0) return
  do  i = 0, k
     a(ja+i) = recrd(jr+i)
  enddo
  jr = jr+k+1
  ja = ja+k+1
  if (lqchr) then
     ja = ja-1
     a(ja-1) = a(ja)
     goto 10
  endif
  !   ... Exit if no substitution
  if (jr >= reclr .OR. ja > recla) return
  !   --- Handle 'nesting' of cex by treating outer cex as a literal ---
  k = 0
  call chrps2(recrd(jr),cex,2,reclr-jr-1,k,kkk)
  !   ... can't find corresponding close of cex
  if (kkk == 0) then
     call pvfil3('missing "',cex(2:2),'" in line:',recrd,1,reclr)
     !   ... cex is nested.  This pass, quote outer cex literally
  elseif (kkk == 1) then
     do  i = 0, k-1
        a(ja+i) = recrd(jr+i)
     enddo
     jr = jr+k
     ja = ja+k
     !         print *, 'testing',(a(kkk), kkk=1,ja)
     pvfil2 = .true.
     goto 10
  endif

  !   --- Parse and replace expression ---
  !       a(ja) is 1st char after \cex
  !   ... Case expr?string1:string2
  if (recrd(jr) == '?') then
     match = recrd(jr+1)
     !         Position of cex(2:2)
     k = k+jr
     !         Start of expression, save for error message
     m = jr
     !         Parse the expression
     jr = jr+2
     if ( .NOT. a2bin(recrd,i,2,0,match,jr,reclr)) &
          call pvfil3('could ','not ','parse:',recrd,m,reclr)
     !         First char after expr
     kkk = jr
     !         Position of matching char
     call chrpos(recrd,match,reclr,jr)
     if (jr >= reclr) call pvfil3( &
          'missing matching "',match,'" in line: ',recrd,m,reclr)
     !         First string
     if (i /= 0) then
        m = kkk
        !         Second string
     else
        m = jr+1
        match = cex(2:2)
     endif
     !         Copy one string
     call strcop(a(ja-1),recrd(m),reclr-m,match,i)
     !         Update ja,jr before resuming parsing
     ja = ja+i-2
     a(ja) = ' '
     jr = k+1
     goto 10
  endif

  !   ... Start of expression should be a variable name; put into aa
  aa = ' '
  call strcop(aa,recrd(jr),reclr-jr,cex(2:2),i)
  match = aa(i:i)
  aa(i:i) = ' '
  call wordg(aa,10,'A-Za-z_0-9',1,i1,i2)
  !       call skipbl(aa,ctlen,i2)
  if (i2+1 < i) then
     match = aa(i2+1:i2+1)
     aa(i2+1:i2+1) = ' '
  endif

  !   ... Case string aa matched an entry in ctbl
  !       i=index to ctbl if start of aa holds member of ctbl
  call tokmat(aa,ctbl,nchr,ctlen,' ',i,m,.false.)
  if (i >= 0) then
     aa2 = ctbl(i+1,2)
     !         Simple character variable substution
     if (match == cex(2:2)) then
        i1 = 1
        i2 = ctlen
        goto 20
        !         If not subexpression, subst does not involve a char variable
     elseif ( match /= '(') then
        i = -1
        goto 20
     endif
     ip = i2+1
     call skipbl(aa,ctlen,ip)
     ip = ip+1
     !         Case aa = ctbl(i+1,1)(/s1/s2/[,n1,n2])
     if (aa(ip:ip+2) == ':e)') then
        call strip(aa2,i1,i2)
        ja = ja-2
        call bin2a(' ',0,0,i2,2,0,recla,a,ja)
        ja = ja+1
        jr = jr+ip+3
        goto 10
     elseif (aa(ip:ip) == '/') then
        istr1 = ip+1
        istr2 = istr1
        call cpstr(aa,ctlen,0,'/',istr2,ip1,aastr)
        if (istr2 <= istr1 .OR. istr2 >= ctlen) goto 999
        istr2 = istr2+1
        call cpstr(aa,ctlen,0,'/',istr2,ip2,aasub)
        if (istr2 <= istr1 .OR. istr2 >= ctlen) goto 999
        m = istr2+1
        if (aa(m:m) == ',') then
           k = a2vec(aa,len(aa),m,2,',)',2,2,2,ix,ii)
           if (i1 < 1) call pvfil3( &
                'illegal',' parameter',' in',recrd,1,reclr)
           if (k /= 2) goto 999
        else
           i1 = 1
           i2 = 1
        endif
        m = m+1
        ip1 = ip1-1
        ip2 = ip2-1
        ip = ctlen
        call skpblb(aa2,ctlen,ip)
        ip = ip+1
        iaa21 = 1
        iaa22 = ip-ip1+1
        iaa3 = 0
        aa3 = ' '
        do  25  j = 1, i2
           do  24  iaa2 = iaa21, iaa22
              iaa3 = iaa3+1
              aa3(iaa3:iaa3) = aa2(iaa2:iaa2)
              if (aa2(iaa2:iaa2+ip1-1) == aastr(1:ip1)) then
                 if (j < i1) then
                    iaa21 = iaa2+1
                    goto 25
                 else
                    aa3(iaa3:iaa3+ip2-1) = aasub(1:ip2)
                    iaa21 = iaa2+ip1
                    iaa3 = iaa3+ip2-1
                    if (j == i2) aa3(iaa3+1:) = aa2(iaa21:)
                    goto 25
                 endif
              endif
24         enddo
           call pvfil3('could not find substr "',aastr(1:ip1),'" in', &
                recrd,1,reclr)
25      enddo
        aa2 = aa3
        i1 = 1
        i2 = ctlen
        call skpblb(aa2,ctlen,i2)
        i2 = i2+1
        !         Case aa = ctbl(i+1,1)('char-list' [,n])
     elseif (aa(ip:ip) == '''' .OR. (aa(ip:ip) == '"')) then
        i1 = ip+1
        i2 = ip
        !           call cpstr(aa,ctlen,1,aa(ip:ip),ip,ip1,aasub)
        call eostr(aa,ctlen,21,aa(ip:ip),i2)
        if (i2 >= ctlen) goto 999
        i2 = i2-1
        kkk = i2+2
        if (aa(kkk:kkk) == ')') then
           k = 1
        elseif (aa(kkk:kkk) == ',') then
           if ( .NOT. a2bin(aa,k,2,0,')',kkk,len(aa))) goto 999
        else
           goto 999
        endif
        call wordg(ctbl(i+1,2),100,aa(i1:i2),k,i1,i2)
        ja = ja-2
        call bin2a(' ',0,0,i2+1,2,0,recla,a,ja)
        ja = ja+1
        jr = jr+kkk+1
        goto 10
        !         Assume aa = ctbl(i+1,1)(n1,n2)
     else
        ip1 = i2+1
        ip = ip1
        k = a2vec(aa,len(aa),ip,2,',)',2,2,2,ix,ii)
        if (k /= 2) goto 999
        !           Add to m so that jr is incremented past (..)
        m = m + ip-ip1+1
     endif
  endif
  !   ... End of test for character variables.  If string substitution,
  !       i>0 and aa2(i1:i2) = (possibly modified) substring of ctbl
  !       and m = number of characters to advance in recrd
20 continue

  !   ... Substitute a character variable
  if (i >= 0) then
     jr = jr+m
     m = i2
     call skpblb(aa2,ctlen,i2)
     m = min(m,i2+1)
     a(ja-1) = ' '
     !         print *, 'before 33',(a(kkk), kkk=1,ja-2+0)
     do  k = i1, m
        a(ja+k-i1-1) = aa2(k:k)
     enddo
     !         print *, 'after 33',k1,k2,(a(kkk), kkk=1,ja+m-i1-1)
     ja = ja+m-i1

     !   ... Numerical expression
  elseif (a2bina(recrd,res,4,0,cex(2:2),jr,reclr)) then
     ja = ja-2
     if (dabs(res) > 1) then
        call bin2a('g9',0,0,res,4,0,recla,a,ja)
     else
        call bin2a('g9:10',0,0,res,4,0,recla,a,ja)
     endif
     ja = ja+1

  else
     !   ... if a name for a vector, replace with ascii rep of entire vector
     call sizsvv(aa,ival,0,k)
     if (ival > 0) then
        ja = ja-2
        do  35  m = 1, k
           call getsvv(' ',ival,1,m,m,res)
           if (dabs(res) > 1) then
              call bin2a('g9',0,0,res,4,0,recla,a,ja)
           else
              call bin2a('g9:10',0,0,res,4,0,recla,a,ja)
           endif
           ja = ja+1
           a(ja) = ' '
35      enddo
        call chrpos(recrd,cex(2:2),reclr,jr)
        jr = jr+1
     else
        print *, 'rdfile: bad expression in line'
        print *, (recrd(j), j=0, jr-1), '  ...  ', &
             (recrd(j), j=jr, reclr-1)
        call cexit(-1,1)
     endif
  endif
  goto 10

999 continue
  call pvfil3('could',' not',' parse',recrd,1,reclr)
end function pvfil2

subroutine pvfil3(s1,s2,s3,recrd,jr,reclr)
  !- Error exit
  implicit none
  integer :: reclr,i,ii,jr
  character(1) :: recrd(0:1)
  character*(*) s1,s2,s3
  character aa*100
  aa = 'rdfiln '// s1 // s2 // s3 // ' ... '
  call skpblb(aa,len(aa),ii)
  ii = ii+3
  do  i = ii, min(reclr-jr,len(aa))
     aa(i:i) = recrd(i+jr-ii)
  enddo
  !      call setpr(20)
  call fexit(-1,009,aa,0d0)
end subroutine pvfil3

logical function rdstrn(unit,a,len,lopt)
  !- Read one line from logical unit
  integer :: unit,len
  logical :: lopt
  character(1) :: a(len)
  rdstrn = .true.
  read(unit,10,end=20,err=20) a
  if (lopt) print 10, a
  return
10 format(2048a1)
20 rdstrn = .false.
  return
end function rdstrn

subroutine macset(args,err)
  !- Routines for macro definitions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   args   :ascii string defining macro, of form:
  !i          :macro_name(arg1,arg2,..) expand_string
  !i   err
  !o Outputs
  !l Local variables
  !l   nmaca  : nmaca(i) = number of argments for macro i; limited to nmxarg
  !l   macarag:macarag(j,i) is the variable names for argument j to macro i
  !l          :First argument is macro name
  !r Remarks
  !r
  !u Updates
  !u   08 Jan 07  Increase the max number of macros
  !u   19 Dec 02  First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nmac,nmacx,nmxarg,err
  character*(*) args,strout

  ! ... Local parameters
  parameter (nmacx=10,nmxarg=5)
  !     logical savvar(nmacx)
  !     double precision locval(nmacx)
  integer :: nmaca(nmacx),nchar(nmxarg,2)
  character(16) :: macarg(nmxarg,nmacx),argsub(nmacx)
  character(256) :: macexp(nmacx),strn,strn2
  integer :: j,j1,j2,i1,i2,iarg,imac,i,ls,errin
  save nmac,nmaca,macarg,macexp
  data nmac /0/

  ! ... too many macros ... exit err=-1
  err = -1
  if (nmac >= nmacx) goto 99

  ! ... No macro defined ... exit err=-2
  call word(args,1,j1,j2)
  err = -2
  if (j1 > j2) goto 99
  ! ... macro name doesn't start with a letter ... exit err=-3
  err = -3
  call wordg(args(j1:),10,'A-Za-z_0-9',1,i1,i2)
  if (i2 < i1) goto 99
  ! ... Assign macro name
  macarg(1,nmac+1) = args(j1+i1-1:j1+i2-1)
  ! ... First letter after macro name isn't '(' ... exit err=-4
  err = -4
  i1 = i2+j1
  call nword(args,1,i1,i2)
  if (args(i1:i1) /= '(') goto 99

  ! --- Assign macro arguments ---
  i1 = i1+1
  iarg = 1
10 continue

  ! ... Error if macro has too many arguments
  iarg = iarg+1
  err = -5
  if (iarg > nmxarg) goto 99

  ! ... Error if no terminator to macro
  call nwordg(args,1,'),',1,i1,i2)
  err = -6
  if (i2 >= len(args)) goto 99

  macarg(iarg,nmac+1) = args(i1:i2)
  i1 = i2+2
  if (args(i2+1:i2+1) /= ')') goto 10

  ! ... get macro string ... error if none
  err = -7
  call nword(args,1,i1,i2)
  if (i2 < i1) goto 99
  macexp(nmac+1) = args(i1:i2)

  ! ... Macro arguments complete ... increment number of macros by 1
  nmac = nmac+1
  nmaca(nmac) = iarg
  err = 0
  goto 99

  ! --- macro evaluation ---
  entry macevl(args,strout,err)
  ! Input err:   0 -> just return index to macro in err if it exists
  ! Input err: <>0 -> return index to macro, and strout = expansion

  errin = err

  ! ... No macro defined ... exit err=-2
  call word(args,1,j1,j2)
  err = -2
  if (j1 > j2) goto 99

  ! ... No match to existing macro names ... exit err=-1
  err = -1
  call wordg(args(j1:),10,'A-Za-z_0-9',1,i1,i2)
  if (i2 < i1) goto 99
  do  imac = 1, nmac
     if (macarg(1,imac) == args(j1+i1-1:j1+i2-1)) goto 15
  enddo
  goto 99

  ! ... macro found
15 continue
  if (errin == 0) then
     err = imac
     return
  endif

  ! ... First letter after macro name isn't '(' ... exit err=-4
  err = -4
  i1 = i2+j1
  call nword(args,1,i1,i2)
  if (args(i1:i1) /= '(') goto 99

  ! --- Assign macro arguments ---
  i1 = i1+1
  iarg = 1
20 continue

  ! ... Error if macro has too many arguments
  iarg = iarg+1
  err = -5
  if (iarg > nmaca(imac)) goto 99

  ! ... Error if no terminator to macro
  call nwordg(args,1,'),',1,i1,i2)
  err = -6
  if (i2 >= len(args)) goto 99

  ! ... For this variable, evaluate argument and load var table with res
  argsub(iarg) = args(i1:i2)
  i1 = i2+2
  if (args(i2+1:i2+1) /= ')') goto 20

  ! ... count number of characters in each macro string
  do  i = 2, nmaca(imac)
     call word(macarg(i,imac),1,j1,j2)
     nchar(i,1) = j2
     call word(argsub(i),1,j1,j2)
     nchar(i,2) = j2
  enddo
  !     i = which string holds current value of macro
  !     j = pointer to current index in string
  i = 1
  j = 1
  strn = macexp(imac)
  ls = len(strn)
  call skpblb(strn,len(strn),i)
  ls = i+1
  ! --- macro substitution ---
30 continue
  do  i = 2, nmaca(imac)
     if (strn(j:j+nchar(i,1)-1) == macarg(i,imac)) then
        !         Also next character cann be part of a word
        call wordg(strn(j+nchar(i,1)-1:),10,'A-Za-z_0-9',1,i1,i2)
        if (i2 <= i1) then
           strn2(1:ls) = strn
           strn(j:j+nchar(i,2)-1) = argsub(i)
           ls = ls + nchar(i,2)-nchar(i,1)
           strn(j+nchar(i,2):ls) = strn2(j+nchar(i,1):)
           if (nchar(i,2)-nchar(i,1) < 0) then
              strn(ls+1:ls+1) = ' '
           endif
           j = j + nchar(i,2)-1
           goto 35
        endif
     endif
  enddo
35 continue
  !     Move past current word
  call nwordg(strn,10,'A-Za-z_0-9',1,j,j2)
  if (j2 > j) then
     j = j2+1
  else
     j = j+1
  endif
  if (j <= ls) goto 30
  err = imac
  strout = strn(1:ls)
  ! --- Error exit ---
99 continue
end subroutine macset

subroutine findctrlstart(nfilin)
  ! if we find 'ctrlstart', locate reading at the next line of ctrlstart.
  !     this is useful if you like to use script, GWinput, ctrl in a file.
  implicit none
  integer:: nfilin
  character(len=9):: strn, ccc
  do
     read(nfilin,"(a)",err=1010,end=1010) strn
     if(strn == 'ctrlstart') return
  enddo
1010 continue
  rewind(nfilin)
end subroutine findctrlstart

subroutine parchv(recrd,recl,mxchr,cex,sep,opts,nchr,ctbl,j)
  !- Parses a string for one or more character variable declarations
  ! ----------------------------------------------------------------
  !i Inputs
  !i   recrd(0:*): string recrd is parsed from j to recl-1
  !i   recl :  record length: characters limited to recrd(0:recl-1)
  !i   mxchr:  maximum number of entries in character table
  !i   cex  :  characters demarcating string substitution, eg '{}'
  !i   sep  :  characters demarcating symbol and string-value; see Remarks
  !i   opts:   0: var declaration of existing variable is supressed
  !i           1: var declaration of existing variable supersedes
  !i        10's digit:
  !i           1  variable's value is substituted for shell environment
  !i              variable of the same name
  !i   nchr :  number of variables defined so far
  !i   ctbl :  table of character variables (col 1) and values (col 2)
  !i   j:      parsing starts at char j (origin at 0)
  !o Outputs
  !o   nchr :  is updated
  !o   j:      last character parsed
  !r Remarks
  !r   Declarations are typically of the form symbol=string-value
  !r   The separator '=' need not be '=', but is defined by 'sep'
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: recl,nchr,j,mxchr,opts,ctlen
  parameter (ctlen=120)
  character(1) :: recrd(0:*),ctbl(mxchr,2)*(ctlen),sep*(*)
  ! Local parameters
  integer :: i,k,j0,jr,opt1,opt2
  character ctmp*1,strn*(ctlen),cex*2

  opt1 = mod(opts,10)
  opt2 = mod(opts/10,10)

  ! --- Parse for next declaration ---
2 call skipbl(recrd,recl,j)
  if (j >= recl .OR. nchr >= mxchr) return
  i = j
  call chrps2(recrd,sep,len(sep),recl,i,k)
  if (i >= recl) return
  !      print *, (recrd(j0),j0=0,j)
  !      print *, (recrd(j0),j0=0,i)
  ctmp = sep(k:k)
  i = i-j
  if (i > ctlen) call fexit(-1,9,'parchv: name too long',0)
  call strcop(strn,recrd(j),i,ctmp,k)
  strn(i+1:i+1) = ' '
  call tokmat(strn,ctbl,nchr,ctlen,' ',i,k,.false.)
  ! ... Replace existing string definition
  if (i >= 0 .AND. opt1 == 1) then
     i = i+1
     ! ... Prior declaration of variable; ignore this declaration
  else if (i >= 0 .AND. opt1 == 0) then
     j = j+k
     call eostr(recrd,recl,11,' ',j)
     !        call skipbl(recrd,recl,j)
     !        call skp2bl(recrd,recl,j)
     goto 2
     ! ... Create new string definition
  else
     i = nchr+1
     ctbl(i,1) = ' '
     call strcop(ctbl(i,1),strn,ctlen,' ',k)
     nchr = i
  endif

  ! --- Copy variable contents into the table ---
  ctbl(i,2) = ' '
  j = j+k
  call skipbl(recrd,recl,j)
  if (j >= recl) return
  j0 = 1
  ctmp = recrd(j)
  if (ctmp /= '"') ctmp = ' '
  if (ctmp == '"') j=j+1
  jr = j
  call chrpos(recrd,ctmp,recl,jr)

  call pvfil1(jr,ctlen,j,recrd,cex,nchr,mxchr,ctbl,j0,ctbl(i,2),0)
  if (ctmp == '"') j=j+1
  if (opt2 == 1) call get_environment_variable(ctbl(i,2),ctbl(i,2))
  goto 2
end subroutine parchv
