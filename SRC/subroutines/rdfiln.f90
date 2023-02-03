module m_rdfiln! preprocessor. ctrl.* is conveter to ctrl_preprocessd.*
  ! This routine depends on a2bin.f90 (ascii to num) and symvar.f90 (console inputs).
  ! File to File, no side effect. 
  ! Not maintain this routine. We should rewrite a substution by python.
  public M_rdfiln_init
  private
contains
  subroutine m_rdfiln_init() !Read recln,recrd,nrecs,mxrecs for m_lmfinit_init
    !     -vfoobar replaced simplified contents of ctrl into recrd
    use m_lgunit,only:stdo
    use m_ext,only: sname
    use m_lmfinit,only: recln 
    integer:: nfilin,master=0,ierr,i,procid=0,mpipid,ncp
    character(8) :: alabl
    character:: strn*1000
    logical:: fileexist,lshowp,lshow,cmdopt0,mlog
    integer,parameter:: mxrecs=10000
    character*(recln):: recrd(mxrecs)
    integer:: nrecs    
    include "mpif.h"
    procid = mpipid(1)
    inquire(file='ctrl.'//trim(sname),exist=fileexist)
    if( .NOT. fileexist) call rx("No ctrl file found! ctrl."//trim(sname))
    open(newunit=nfilin,file='ctrl.'//trim(sname))
    call findctrlstart(nfilin) ! if a tag 'ctrlstart' in ctrl, ctrl is read from the tag.
    alabl = '#{}% ct '
    call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln,nrecs) !read ctrl into recrd
    close(nfilin)
    open(newunit=ncp,file='ctrl_preprocessed.'//trim(sname))
    write(ncp,"(i10, ' #---- preprocessed ctrl file -------')")nrecs
    do i = 1, nrecs
       write(ncp,"(a)")trim(recrd(i))
    enddo
    close(ncp)
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
  logical :: pvfil2,lpv2
  character(512) :: aa
  integer :: ia,ja0,jaa
  character:: cex*2
  do  ia = ja, recla
     a(ia) = ' '
  enddo
  ja0 = ja
  lpv2=pvfil2(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,ja,a)
  if ( .NOT.lpv2 ) return
  do 
     aa = ' '
     do ia = ja0, ja
        aa(ia:ia) = a(ia)
     enddo
     if (t > 2) then
        call skpblb(aa,ja,ja)
        print 333, '#rf subst: "', aa(1:ja+1),'"'
333     format(a,a,a)
     endif
     jaa = ja0-1
     ja = ja0
     lpv2=pvfil2(recla,recla,jaa,aa,cex,nchr,mxchr,ctbl,ja,a)
     if (.not.lpv2) exit 
  enddo
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

subroutine cpstr(strn,lstr,opt,delim,is,io,sout)
  !- Copy a string, excluding delimiters
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   strn    string
  !i   lstr    string length
  !i   opt     1s digit
  !i           1: skip initial blanks
  !i         100s digit
  !i           1: compress repeated blanks into a single blank
  !i        1000s digit
  !i           1: purge delimiter from sout (last char copied unles eos)
  !i   delim   string delimiters marking end of string (any character
  !i           in delim counts as a delimiter
  !i   is      first character to copy (string starts at character 1)
  !o Outputs
  !o   sout    strn copied to sout; see Remarks
  !o   is      smaller of index in strn of last char copied and lstr+1
  !o   io      position of terminator, or 1+ last character copied if
  !o           limit of string reached beforehand
  !r Remarks
  !r   delimiters between " " or ' ' are excluded.  pairs " " and ' '
  !r   are excised in the output string sout
  !u Updates
  !u   01 Nov 01 delimiter may contain more than one character
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: lstr,is,io,opt
  character *(*) strn, delim*(*), sout*(*)
  character ch*60
  integer :: i2,opt0,opt2,opt3,it,ia,nch,ldelim
  sout = ' '
  opt0 = mod(opt,10)
  opt2 = mod(opt/100,10)
  opt3 = mod(opt/1000,10)
  ldelim = len(delim)
  ch = delim // '"'' '
  nch = ldelim + 2
  if (opt2 /= 0) nch = nch+1
  is = is-1
  ia = 1
  ! ... Skip past leading blanks
  if (opt0 /= 0) then
     call skipbl(strn,lstr,is)
  endif
  ! ... Find i2 : points to past last char of the argument
  i2 = is
12 is = i2
  call chrps2(strn,ch,nch,lstr,i2,it)
  if (it <= ldelim .AND. i2 < lstr) i2 = i2+1
  if (i2 > is) sout(ia:ia+i2-is-1) = strn(is+1:i2)
  if (it <= 0 .AND. i2 == lstr) i2 = i2+1
  ia = ia+i2-is
  ! ... A blank encountered ... skip blanks and continue
  if (it == ldelim+3 .AND. i2 < lstr) then
     call skipbl(strn,lstr,i2)
     sout(ia:ia) = ' '
     ia = ia+1
     goto 12
     ! ... A quote encountered ... find match and continue
  elseif (it > ldelim .AND. i2 < lstr) then
     i2 = i2+1
     is = i2
     call chrpos(strn,ch(it:it),lstr,i2)
     call strncp(sout,strn,ia,is+1,i2-is)
     ia = ia+i2-is
     i2 = i2+1
     goto 12
  endif
  io = ia-1
  if (is >= lstr) io = io+1
  !      if (is .ge. lstr) print *, 'hi',io
  if (opt3 == 1) then
     if (i2 <= lstr) sout(io:io) = ' '
  endif
  is = min(i2,lstr+1)
end subroutine cpstr

subroutine eostr(strn,lstr,opt,delim,is)
  !- Mark end of a string
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   strn,is string, and first character (1st character starts at 1)
  !i   delim   string delimiter marking end of string
  !i   opt     1s  digit 1: exclude delimiters between " " or ' '
  !i           10s digit 1: skip initial blanks
  !i                     2: use strn(is) as delimiter;  delim is not used
  !o Outputs
  !i   is      smaller of position in strn of delimiter and lstr+1
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: lstr,is,opt
  character *(*) strn, delim*1
  character ch*3
  integer :: i2,opt0,opt1,it,k
  data ch /' "'''/
  opt0 = mod(opt,10)
  opt1 = mod(opt/10,10)
  ch(1:1) = delim
  if (opt1 == 2) then
     ch(1:1) = strn(is:is)
  else
     is = is-1
  endif
  ! ... Skip past leading blanks
  if (opt1 == 1) then
     call skipbl(strn,lstr,is)
  endif
  ! ... Find i2 : points to past last char of the argument
  k = 1
  if (opt0 /= 0) k = 3
  i2 = is
12 is = i2
  call chrps2(strn,ch,k,lstr,i2,it)
  if (it == 1 .AND. i2 < lstr) i2 = i2+1
  ! ... A quote encountered ... find match and continue
  if (it > 1 .AND. i2 < lstr) then
     i2 = i2+1
     is = i2
     call chrpos(strn,ch(it:it),lstr,i2)
     i2 = i2+1
     goto 12
  endif
  is = min(i2,lstr+1)
end subroutine eostr
subroutine strip(str,i1,i2) !- Returns indices to first and last nonblank characters in a string
  implicit none
  integer :: i1,i2
  character*(*) str
  integer :: i
  i1 = 0
  do  i = 1, len(str)
     if(str(i:i) /= ' ') then
        i1 = i
        goto 2
     endif
  enddo
  i1 = 1
  i2 = 0
  return
2 continue
  i2 = len(str) + 1
  do  i = len(str), 1, -1
     if(str(i:i) /= ' ') then
        i2 = i
        exit
     endif
  enddo
end subroutine strip
subroutine addsvv(nam,nelt,ival)
  !- Add a symbolic vector to list
  ! ----------------------------------------------------------------
  !i Inputs
  !i   nam:  name of variable
  !i   nelt: number of elements of the vector
  !o Outputs
  !o   ival  index to which variable is declared or accessed
  !r Remarks
  !r   addsvv  adds a symbolic name and value to the internal table,
  !r           and allocates memory for the vector.
  !r   lodsvv  sets a range of elements of a vector associated with
  !r           a name or an index, depending on iopt.
  !r           iopt=0: index associated with name
  !r           iopt=1: name associated with index
  !r   getsvv  gets a range of elements of a vector associated with
  !r           a name or an index, depending on iopt.
  !r   sizsvv  returns the length of a vector associated with
  !r           a name or an index, depending on iopt.
  !r   numsvv  returns the number of variables now declared
  !r   watsvv  returns name associated with index
  !r
  !r   Compiler must have either F90 or POINTER capability
  !u Updates
  !u   18 Jan 06 works with F90 compiler
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character*(*) nam
  double precision :: vec(1)
  integer :: ival,first,last,nelt,nvar,ifi,iopt
  ! Local parameters
  integer :: mxnam,namlen
  parameter (mxnam=24,namlen=16)
  character*(namlen) symnam(mxnam), tmpnam
  integer :: size(mxnam)
  integer :: nnam,i,io,i1,i2,i2x
  double precision :: x1,xn
  type row
     real(8), allocatable :: p(:)
  end type row
  type(row) :: symptr(mxnam)
  save symptr
  save symnam, size, nnam
  data nnam /0/
  ! --- Start of addsvv ---
  nnam = nnam+1
  if (nnam > mxnam) call rx('addsvv: too many names')
  symnam(nnam) = nam
  ival = nnam
  call locase(symnam(nnam))
  allocate (symptr(nnam)%p(1:nelt))
  symptr(nnam)%p=0d0 
  size(nnam) = nelt
  return

  ! --- lodsvv, getsvv ---
  entry lodsvv(nam,ival,iopt,i1,i2,vec)

  io = -1
  goto  10

  entry getsvv(nam,ival,iopt,i1,i2,vec)

  io = 1
  goto  10

  entry sizsvv(nam,ival,iopt,i1)

  io = -2
  goto  10

  ! --- lodsvv, getsvv ---
  entry numsvv(nvar)
  nvar = nnam
  return

  ! --- watsvv ---
  entry watsvv(nam,ival)
  nam = ' '
  if (ival <= nnam) nam = symnam(ival)
  return

  ! --- Print out table ---
  entry shosvv(first,last,ifi)
  write(ifi,332)
332 format('  Vec       Name            Size   Val[1..n]')
  do  60  i = max(first,1), min(last,nnam)
     call dpscop(symptr(i)%p,x1,1,1,1,1d0)
     call dpscop(symptr(i)%p,xn,1,size(i),1,1d0)
     write(ifi,333) i, symnam(i), size(i), x1, xn
60 enddo
333 format(i4, 4x, a20, i4, 2g14.5)
  return

  ! --- Find an index associated with a name ---
10 continue
  ! ... If iopt=0, find the index associated with this name
  if (iopt == 0) then
     ival = 0
     tmpnam = nam
     call locase(tmpnam)
     do  16  i = 1, nnam
        if (tmpnam /= symnam(i)) goto 16
        ival = i
        goto 20
16   enddo
  endif
  ! --- Set/Retrieve portions of an array[index], depending on io ---
20 continue
  if (io == 0) return
  if (io == -2) then
     i1 = size(ival)
     return
  endif
  ! ... Return unless ival in range
  if (ival <= 0 .OR. ival > nnam) return
  i2x = min(i2,size(ival))
  if (i2x < i1) return
  if (io == -1) call dpscop(vec,symptr(ival)%p,i2x-i1+1,1,i1,1d0)
  if (io ==  1) call dpscop(symptr(ival)%p,vec,i2x-i1+1,i1,1,1d0)
  return
end subroutine addsvv

subroutine parsvv(recrd,recl,indx,mxelt,i1,ip)
  !- Parses a string for one or more elements of a vector variable
  !     implicit none
  ! Passed parameters
  integer :: recl,ip,mxelt,indx,i1
  character recrd*(100)
  ! Local parameters
  double precision :: res
  integer :: nelt,i,k,ix,a2vec
  nelt = 0
  do  33  i = 1, 999
     call skipbl(recrd,recl,ip)
     if (ip >= recl .OR. nelt >= mxelt) goto 99
     k = a2vec(recrd,recl,ip,4,' ',1,1,1,ix,res)
     if (k == -1) return
     call lodsvv(' ',indx,1,i1+nelt,i1+nelt,res)
     nelt = nelt+k
33 enddo
99 continue
end subroutine parsvv

integer function awrite(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  !- Formatted output, with ascii conversion of binary numbers
  !i ifi: <>0, local output string written to abs(ifi);
  !i      <=0, sout copied to local output string initially;
  !i           local output string copied back to sout on exit
  !i       >0, sout unaltered on exit
  !i mxln: abs(mxln) = maximum number of characters to copy
  !i mxln < 0: suppress trailing blanks when writing to logical unit.
  !o sout:     output string (see ifi, above)
  !r Characters are copied from fmt to the output string sout, which is
  !r   then (optionally) written to logical unit ifi.  Pointer ip keeps
  !r   track of the current position for writing to sout.  Copy
  !r   is literal except when a control char % is encountered.
  !r   % characters do one of several functions:
  !r   %l writes to sout an ascii representation of logical argument a_j
  !r      (NB: j is 1 for first conversion, 2 for second, etc).
  !r   %i writes an integer argument a_j
  !r   %d, %e, %g, %G, %D and %F write ascii rep of double precision a_j
  !r     'd' writes in decimal notation
  !r     'e' writes in exponential notation
  !r     'g' and 'G' take the minimum size of 'd' and 'e'; 'g' is to
  !r                 specify relative precision, 'G' absolute precision.
  !r     All of the above generate 'pretty' ascii representations
  !r     (see pretty.f)
  !r     'D' mimics the fortran 'f' format and is intended for output in
  !r     fixed columns.
  !r     'F' puts the number in a specified space, using whatever form
  !r     produces the most decimal places of precision.
  !r   %% quotes a "%" literally
  !r   %a shifts ip past last nonblank character
  !r   %f shifts ip forward
  !r   %b shifts ip backward
  !r   %p sets   ip to a fixed value
  !r   %t is obsolete
  !r   %x blanks the output string
  !r   %z can suppress leading the leading zero in a decimal fraction
  !r   %W shifts ip forward until a whitespace is encountered
  !r   %w shifts ip forward until a non-whitespace is encountered
  !r   %c closes up whitespace around ip
  !r   %o opens  up whitespace around ip
  !r   %u if numerical argument = NULLI output 'null' instead of number
  !r      Optional argument n1:
  !r      0 turn off null option, for this and future calls
  !r      1 (or default) set null option, for this and future calls
  !r     >1 set null option for this call only
  !r     <0 turn off null option for this call only
  !r   %? conditionally parses one of two strings (see below)
  !r   %j increments argument jumps over call arguments
  !r   %N is turned into a newline, calling nlchar to get newline
  !r Most control characters have optional arguments.
  !r For d,e,g,G,D,F,l,i the general syntax is:
  !r   %[n1][:n2][,n3][;n4][#n5]x, with x=d,e,g,G or F
  !r Here n1..n5 are integer expressions:
  !r   n1 number of values to convert (a_j is regarded as a vector)
  !r   n2 number of blank spaces preceding first character
  !r      n2<0 => subtract one space if argument is negative
  !r   n3 minimum number of digits to display (after '.' for 'd'
  !r      and 'G' and total number for 'e' and 'g')
  !r   n4 round to n4 decimal places (absolute for 'd' and 'G',
  !r      and relative for 'e' and 'g')
  !r   n5 if conversion uses less than n5 characters, append trailing
  !r      blanks to fill (used for lining data in columns)
  !r For D the meanings of n2..n4 are different:
  !r   n2 is not used
  !r   n3 number of digits after decimal
  !r   n4 field width
  !r For F:
  !r   n2 is not used
  !r   n3 is not used
  !r   n4 is the field width
  !r For l and i:
  !r   n3 is the field width
  !r For z, j, p, a, f, o, b, and and the general syntax is:
  !r   %[n1]x, with x=z, p, a, f, b, u
  !r   n1 repeats (f, b)
  !r   n1 1=>suppresses leading 0, 0=>ensures its presence (z)
  !r   For u, see above
  !r NB: there is an option to substitute for any of n1..n4 one of the
  !r arguments a_j.  This is done by using the character 'n' instead
  !r of some integer expression (eg %n:n,5d).  awrite uses the
  !r next argument a_j is used for n, and increments j.  Thus, %n:n,5d
  !r consumes the next three a_j, the first describing the number
  !r of elements to write, the second the number of spaces between
  !r arguments.
  !r For ? the general syntax is
  !r   %?QexprQstr1Qstr2Q
  !r   str1 is parsed if "expr" evaluates to nonzero, otherwise str2 is.
  !r   Q is some character, eg ';'.  It should NOT be some character
  !r   that may be confused as part of "expr", like '?' or '/'.
  !r   The next argument argument a_j is temporarily set to an integer
  !r   value and temporarily named `n', which may be used in 'expr'.
  !r   Also the current number of characters in the string is temporarily
  !r   assigned to `p'.  Finally, as a special case for a expression
  !r   involving strings, the following is permitted:
  !r     %c==X
  !r   where X is some character.  This expression evaluates to nonzero
  !r   if the character output string at the current position is equal
  !r   to X; otherwise it evaluates to zero.
  !r   Example:
  !r     call awrit2('three plus one is %?;n==1;%i;four;, no?',mxlen,
  !r                 s,mxlen,-i1mach(2),m,4)
  !r     prints out "three plus one is 4, no?" if m equals 1; otherwise
  !r     prints out "three plus one is four, no?"
  !u Updates
  !u   13 Oct 07 Modified lnull, for permanent option
  !u   02 Aug 07 Added %u: outputs null string when arg=NULLI (bin2av)
  !u   27 Feb 02 Added %c==X type of conditional expression
  !u    8 May 01 addition of n5 modifier described above
  ! ----------------------------------------------------------------
  implicit none
  integer :: ifi,mxln
  character*(*) fmt,sout
  double precision :: a1(*),a2(*),a3(*),a4(*),a5(*),a6(*),a7(*),a8(*)
  integer :: i,lfmt,ip,jp,cast,ia,iff,i2,iterm,ndec,iv(5),nblk,nx,j, &
       mxl,ls,icond,ivawrt,lens,nsyv,ires,fw
  equivalence (iv(1),i2),(iv(2),nblk),(iv(3),nx),(iv(4),ndec), (iv(5),fw)
  logical :: a2bin,ltmp,lnull,lnulls
  double precision :: holdn,holdp,xx
  parameter (lens=1024)
  character*(lens) s,ss,fchr*29,fm*20,cc*1,ccond*1
  save lnulls
  data fchr /' :;,#irdegGltapfbzDwWcoxu?jNF'/
  data lnulls /.false./
  logical:: l_exec
  ! ... ia,iff,ip: indices to current argument, pos in fmt, pos in s
  mxl = iabs(mxln)
  s = ' '
  if (ifi <= 0) s = sout 
  ip = 0  ! ... ip holds the current position in the output string
  iff = 0 ! ... iff holds the current position in the format string
  ia = 1 ! ... index to current argument in the argument list
  icond = 0 ! ... icond nonzero when in the middle of a conditional expression
  ccond = ' ' ! ... ccond is the terminating character for conditional expression
  call numsyv(nsyv) ! ... hold on to n,p in vars table; we need them as local variables
  call getsyv('p',holdp,j)
  call getsyv('n',holdn,j)
  lfmt = len(fmt)
  ls = len(s)
  lnull = lnulls
 ! --- Parse next character in fmt ---
19 ia = ia-1
20 iff = iff+1
  !  ...  End of fmt
  if (iff > lfmt) goto 10
  !  ...  Character terminating conditional string
  if (icond > 0 .AND. fmt(iff:iff) == ccond) then
     if (icond == 2) then
        call chrpos(fmt,ccond,lfmt,iff)
        iff = iff+1
     endif
     icond = 0
     goto 20
  endif
  !  ...  Any non-% character
  ! ino      if (fmt(iff:iff) .ne. '%' .or. fmt(iff:iff+1) .eq. '%%') then
  l_exec=.false.
  if (fmt(iff:iff) /= '%') l_exec= .TRUE. 
  if (iff+1<=len(fmt)) then
     if (fmt(iff:iff+1) == '%%')  l_exec= .TRUE. 
  endif
  if (l_exec) then
     !         print *, 'parsing non-%:', iff, fmt(1:iff)
     ip = ip+1
     if (ip <= min(ls,mxl)) s(ip:ip) = fmt(iff:iff)
     if (iff+1<=len(fmt)) then
        if (fmt(iff:iff+1) == '%%') iff = iff+1
     endif
     goto 20
     !   --- Parse % ---
  else
     !     ... Default values for %command  !  print *, 'now parsing %:', iff, fmt(iff:min(lfmt,iff+10))
     nblk = 0
     fw = 0
     nx = 99
     ndec = 0
     i2 = 1
     ia = ia+1
     iff = iff+1
     !     ... iterm flags whether cc is ':;,#', to use later
     j = 0
     call chrps2(fmt(iff:iff),fchr,len(fchr),0,j,iterm)
     !     ... Re-entry if to parse another argument
25   if (iterm >= 2 .AND. iterm <= 5) iff = iff+1
     cc = fmt(iff:iff)
     ires = 1 ! ires is integer argument; set to 1 for default value
     if (cc == '(') then ! Check for an integer expression () preceding command:
        j = 1
        do  35  i = iff+1, lfmt
           if (fmt(i:i) == '(') j=j+1
           if (fmt(i:i) == ')') j=j-1
           if (j == 0) then
              xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
              call lodsyv('n',1,xx,j)
              call lodsyv('p',1,dble(ip),j)
              j = iff-1
              ltmp = .not. a2bin(fmt,ires,2,0,' ',j,i-1)
              call rxx(ltmp,'awrite: failed to parse () in format')
              iff = i+1
              goto 36
           endif
35      enddo
        call rx('awrite: missing matching () in format')
36      continue
     else if (cc == 'n') then ! ... Prior integer argument if next char 'n':
        ires = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
        ia = ia+1
        iff = iff+1
        !     ... Prior integer argument if char an integer:
     else if (cc >= '0' .AND. cc <= '9' .OR. cc == '-') then
        do  22  i = iff, lfmt-1
           j = i+1
           if (fmt(j:j) >= '0' .AND. fmt(j:j) <= '9') goto 22
           j = iff-1
           !             call pshpr(130)
           ltmp = .not. a2bin(fmt,ires,2,0,fmt(i+1:i+1),j,i)
           call rxx(ltmp,'awrite: failed to parse format')
           iff = i+1
           goto 23
22      enddo
23      continue
     endif
     !          print 335, ires,iff,fmt(1:iff)
     ! 335     format('*now ires=',i2,' parsed to ',i3,': ',a)
     !     ... If this was an argument to one of ':;,#'
     if (iterm >= 2 .AND. iterm <= 5) then
        iv(iterm) = ires
        !     ... Otherwise ires is an argument to command
     else
        iv(1) = ires
     endif
     cc = fmt(iff:iff)!     ... Next character is the terminator
     j = 0
     call chrps2(cc,fchr,len(fchr),0,j,iterm)
     !     ... If an argument, run through parse again
     if (iterm >= 2 .AND. iterm <= 5) goto 25
     !     ... Otherwise a command:
     cast = 99
     cc = fmt(iff:iff)
     if (cc == 'l') cast=0
     if (cc == 'i') cast=2
     if (cc == 'r') cast=3
     if (cc == 'd' .OR. cc == 'e' .OR. cc == 'D' .OR. &
          cc == 'g' .OR. cc == 'G' .OR. cc == 'F') &
          cast=4
     if (cc == 't') then
        call rx('awrite: use p, not t')
     endif
     if (cc == 'z') then
        call bin2a0(i2)
        goto 19
     elseif (cc == 'j') then
        ia = ia+i2
        goto 19
     elseif (cc == 'a') then
        call skpblb(s,ls,ip)
        ip = ip+i2
        goto 19
     elseif (cc == 'p') then
        ip = i2
        goto 19
     elseif (cc == 'f') then
        ip = ip+i2
        goto 19
     elseif (cc == 'b') then
        ip = ip-i2
        goto 19
     elseif (cc == 'N') then
        call nlchar(1,s(ip+1:ip+1))
        ip = ip+1
        goto 19
        ! ---     Entry point for conditional expression ---
     elseif (cc == '?') then
        if (icond /= 0) call rx('awrite encountered nested "%?"')
        icond = 1
        call lodsyv('p',1,dble(ip),j)
        xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
        call lodsyv('n',1,xx,j)
        ia = ia+1
        iff = iff+1
        !     ...   ccond is character terminating conditional string
        ccond = fmt(iff:iff)
        !     ...   If next char is '%', expression is of the string type:
        if (fmt(iff+1:iff+1) == '%') then
           if (fmt(iff+1:iff+4) == '%c==') then
              ltmp = fmt(iff+5:iff+5) .eq. s(ip:ip)
              iff = iff+6
           else
              call rxs('awrite: failed to parse : ',fmt(iff:))
           endif
           !     ...   Parse expression
        elseif ( .NOT. a2bin(fmt,ltmp,0,0,ccond,iff,lfmt)) then
           call rx('awrite: failed to parse conditional expr')
        endif
        !     ...   Use first string, or skip to second string
        if (ltmp) then
           icond = 2
        else
           icond = 3
           call chrpos(fmt,ccond,lfmt,iff)
           iff = iff+1
        endif
        goto 19
        !     ... clear string
     elseif (cc == 'x') then
        s = ' '
        goto 19
        !     ... toggle on 'null' option
     elseif (cc == 'u') then
        if (i2 < 0) then
           lnull = .false.
        elseif (i2 == 0) then
           lnulls = .false.
           lnull = .false.
        elseif (i2 == 1) then
           lnulls = .true.
           lnull = .true.
        else
           lnull = .true.
        endif
        goto 19
        ! ...     pad whitespace around ip
     elseif (cc == 'o') then
        if (ip == 0) ip = 1
        ss = s(ip:ls)
        s(ip:ip+i2-1) = ' '
        s(ip+i2:ls) = ss
        ip = ip+i2-1
        goto 19
        ! ...     close up whitespace around ip
     elseif (cc == 'c') then
        do  13  j = ip, 1, -1
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 14
           ip = j
13      enddo
14      continue
        jp = ip-1
        do  15  j = jp+1, mxl
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 16
           jp = j
15      enddo
16      continue
        if (jp-ip+1 > 0) then
           ss = s(jp+1:ls)
           s(ip:ls) = ss
        endif
        goto 19
        ! ...     skip to next nw
     elseif (cc == 'w') then
        do  17  j = ip, mxl
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 19
           ip = j
17      enddo
        ! ...     skip to next whitespace
     elseif (cc == 'W') then
        do  18  j = ip, mxl
           if (s(j:j) == ' ' .OR. s(j:j) == '        ') goto 19
           ip = j
18      enddo
     endif
  endif
  if (cast == 99) call rx('awrite: unknown control: ' // cc)
  ! ---   Generate format for bin2a ---
  fm = ' '
  if (cast == 4) then
     fm = cc
     if (cc == 'G') fm = 'g'
     j = 1
     if (nx /= 99) call bin2a(' ',0,0,nx,2,0,20,fm,j)
     if (cc == 'G') call bin2a(':20',0,0,0,1,0,20,fm,j)
  endif
  ! ---   Convert binary numbers ---
  i2 = i2-1
  if (ia == 1) call bin2av(fm,fw,nblk,ndec,a1,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 2) call bin2av(fm,fw,nblk,ndec,a2,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 3) call bin2av(fm,fw,nblk,ndec,a3,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 4) call bin2av(fm,fw,nblk,ndec,a4,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 5) call bin2av(fm,fw,nblk,ndec,a5,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 6) call bin2av(fm,fw,nblk,ndec,a6,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 7) call bin2av(fm,fw,nblk,ndec,a7,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 8) call bin2av(fm,fw,nblk,ndec,a8,cast,0,i2,' ',mxl,lnull,s,ip)
  goto 20
10 continue
  ip = min(ip,mxl)
  if (mxln < 0) then
     call skpblb(s,ip,ip)
     ip = ip+1
  endif
  ia = iabs(ifi)
  if (ifi /= 0 .AND. ip > 0) write(ia,333) s(1:ip)
333 format(a)
  if (ifi <= 0 .AND. ip > 0) sout = s(1:ip)
  awrite = ip
  call lodsyv('p',1,holdp,j) !! --- Restore or undo symbolic variables p,n ---
  call lodsyv('n',1,holdn,j)
  call clrsyv(nsyv)
  return
end function awrite

subroutine bin2av(fmt,w,nblk,ndec,res,cast,i1,i2,sep,mxln,lnull, &
     outs,ip)
  !- Write out a vector of of numbers using bin2a
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   fmt   :format passed to bin2a
  !i    w    :unused if zero.  If >0,
  !i         :w = minimum spacing between successive numbers
  !i   nblk  :number of blanks preceding each value, or if nblk < 0,
  !i         :|nblk| spaces are prepended for positive numbers
  !i         :|nblk|-1 spaces are prepended for negative numbers
  !i   ndec  :retain a mininimum ndec digits after decimal (see bin2a)
  !i   res   :vector of binaries to convert to ascii string
  !i   cast  :0=logical, 1=char, 2=int, 3=real, 4=double
  !i   i1    :convert numbers res(i1..i2)
  !i   i2    :convert numbers res(i1..i2)
  !i   sep   :separator between numbers
  !i   mxln  :maximum allowed value of ip
  !i   ip    :string position pointer
  !i   lnull :if T, numbers equal to NULLI are turned into NULL
  !o Outputs
  !o   outs  :string containing ascii rep'sn of binary numbers
  !i   ip    :string position pointer updated to end of string
  !r Remarks
  !u Updates
  !u   01 Aug 07 new lnull
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) fmt,outs,sep*1
  integer :: nblk,cast,i1,i2,ip,mxln,ndec,w
  double precision :: res(0:1)
  logical :: lnull
  ! ... Local parameters
  logical :: lneg,llnull
  integer :: nblk2,ivalxx,i,k,ip0,NULLI
  real :: rval
  double precision :: dvalxx
  parameter (NULLI=-99999)

  if (mxln <= 0) return
  do  i = i1, i2
     nblk2 = nblk
     if (nblk < 0) then
        nblk2 = -nblk
        lneg = .false.
        if (cast == 2) lneg = ivalxx(res,i+1) < 0
        if (cast == 3) lneg = rval(res,i+1) < 0
        if (cast == 4) lneg = dvalxx(res,i+1) < 0
        if (lneg) nblk2 = nblk2-1
     endif
     !       Set flag llnul if lnull is ON and argument matches NULLI
     llnull = .false.
     if (lnull) then
        if (cast == 2) llnull = ivalxx(res,i+1) == NULLI
        if (cast == 3) llnull = rval(res,i+1) == NULLI
        if (cast == 4) llnull = dvalxx(res,i+1) == dble(NULLI)
        if (llnull) then
           call skpblb(fmt,len(fmt),ip0)
           fmt(2+ip0:) = ':n'
        endif
     endif
     ip0 = ip
     call bin2a(fmt,nblk2,ndec,res,cast,i,mxln,outs,ip)
     if (llnull) then
        !         If fixed width, leave position of null as is
        if (fmt(1:1) == 'D' .OR. fmt(1:1) == 'F' .OR. &
             (cast == 2 .AND. ndec > 0)) then
           !         Skip if not sufficient space for leading blanks + null
        else if (ip-3 <= 1+ip0+iabs(nblk)) then
           !         Otherwise rewrite null starting at 1+ip0+iabs(nblk)
        else
           outs(1+ip0+iabs(nblk):ip) = 'NULL'
           ip = 4+ip0+iabs(nblk)
        endif
        !         print *, outs(1:ip)
     endif
     if (sep /= ' ' .AND. i < i2) then
        ip = ip+1
        outs(ip:ip) = sep
     endif
     if (w /= 0) then
        do  k = ip+1, ip0+w
           outs(k:k) = sep
           ip = ip+1
        enddo
     endif
  enddo
end subroutine bin2av
real(8) function dvalxx(array,index)
  integer :: index
  !- Returns the double precision value of ARRAY(INDEX)
  double precision :: array(index)
  dvalxx = array(index)
END function dvalxx
real function rval(array,index)
  integer :: index
  !- Returns the real value of ARRAY(INDEX)
  real :: array(index)
  rval = array(index)
end function rval

subroutine awrit8(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  !- Subroutine versions of integer function awrite
  !     implicit none
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  character*(*) sout,fmt
  integer :: ifi,mxln,ip,jp,awrite
  save ip
  entry awrit7(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7)
  entry awrit6(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6)
  entry awrit5(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5)
  entry awrit4(fmt,sout,mxln,ifi,a1,a2,a3,a4)
  entry awrit3(fmt,sout,mxln,ifi,a1,a2,a3)
  entry awrit2(fmt,sout,mxln,ifi,a1,a2)
  entry awrit1(fmt,sout,mxln,ifi,a1)
  entry awrit0(fmt,sout,mxln,ifi)
  ip = awrite(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  return
  entry awrip(jp)
  jp = ip
end subroutine awrit8

! subroutine vwrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8,cast,ires,res)
!   !- Writes either integer or double into ires or res, depending on cast
!   ! ----------------------------------------------------------------------
!   !i Inputs
!   !i   ia    :indicates which of arrays a1..a8 to extract element from
!   !i   n     :which entry in array a_ia
!   !i   a1..a8:element is extracted from one of these arrays
!   !i   cast  :array cast
!   !o Outputs
!   !o   ires  :if cast is integer, result poked into ires
!   !o   res   :if cast is double, result poked into res
!   !u Updates
!   ! ----------------------------------------------------------------------
!   !     implicit none
!   integer :: ia,n,cast,ivawrt,ires
!   double precision :: dvawrt,res
!   double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
!   if (cast == 2) then
!      ires = ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
!   elseif (cast == 4) then
!      res = dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
!   else
!      call rxi('vwrt: cannot handle cast',cast)
!   endif
! end subroutine vwrt

integer function ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
  !     implicit none
  integer :: ia,n,ivalxx
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  ivawrt=99999
  if (ia == 1) ivawrt = ivalxx(a1,n)
  if (ia == 2) ivawrt = ivalxx(a2,n)
  if (ia == 3) ivawrt = ivalxx(a3,n)
  if (ia == 4) ivawrt = ivalxx(a4,n)
  if (ia == 5) ivawrt = ivalxx(a5,n)
  if (ia == 6) ivawrt = ivalxx(a6,n)
  if (ia == 7) ivawrt = ivalxx(a7,n)
  if (ia == 8) ivawrt = ivalxx(a8,n)
end function ivawrt
integer function ivalxx(array,index)
  !- Returns the integer value of ARRAY(INDEX)
  integer :: index
  integer :: array(index)
  ivalxx = array(index)
end function ivalxx

! real(8) function dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
!   !     implicit none
!   integer :: ia,n
!   double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
!   dvawrt=1d99
!   if (ia == 1) dvawrt = a1(n)
!   if (ia == 2) dvawrt = a2(n)
!   if (ia == 3) dvawrt = a3(n)
!   if (ia == 4) dvawrt = a4(n)
!   if (ia == 5) dvawrt = a5(n)
!   if (ia == 6) dvawrt = a6(n)
!   if (ia == 7) dvawrt = a7(n)
!   if (ia == 8) dvawrt = a8(n)
! END function dvawrt

subroutine bin2a(fmt,nblk,ndec,res,cast,count,mxlen,outstr,ip)
  !- Converts number to ascii format, stripping leading blanks, trailing 0
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   fmt: cast=1:   holds the string to be appended to outstr)
  !i        cast=0,2: not used
  !i        cast=3,4: syntax X[#][:sw] where X is one of
  !i                  'd' to write in decimal representation
  !i                  'e' to write in exponential format
  !i                  'g' to use the smaller of 'd' and 'e'
  !i                  '(n.m)' fixed format, mimics fortran fmt (fn.m)
  !i                  'D' also mimics fortran fmt (Fn.m)
  !i                      D# => supplies n=#; arg ndec supplies m
  !i                  'F' fixed field, picking between d and e that
  !i                      F# => # is field width
  !i                      generates the most significant digits
  !i                  See Remarks for further description
  !i   nblk:  strip leading blanks, leaving a maximum of nblk
  !i   ndec:  (cast=3,4 only) retain a mininimum ndec digits after decimal
  !i          point, i.e. do not suppress trailing zeros to ndec.
  !i          ndec=0 does nothing.  ndec should not exceed precsn.
  !i          (cast=2 only): ndec specifies a field width
  !i   res:   binary value to be converted into ascii string
  !i   cast:  cast of res: 0=logical, 1=char, 2=int, 3=real, 4=double
  !i   count: res(count) is to be converted.  NB: count=0 for first entry
  !i   mxlen: maximum length of outstr
  ! o Inputs/Outputs
  ! o  ip:    on input, starting position in outstr for write
  ! o         NB: ip=0 points to first character in string
  ! o  ip:    on output, position of final character written to outstr
  !o  Outputs
  !o   outstr:binary res(count) written in ascii form to outstr(ip:..)
  !r Remarks
  !r  *The string representation of floating point numbers is generated
  !r   by a "prettified" modification of the fortran write statement
  !r   (pretty.f), which includes suppression of trailing zeros and the
  !r   option to include or suppress the leading zero in decimal
  !r   fractions less than 1.  Floating-point formats include:
  !r     'd[n][:sw]' for decimal representation,
  !r     'e[n][:sw]' for exponential representation,
  !r     'g[n][:sw]' uses the minimum length of 'd' and 'e'
  !r     'D[n][:sw]' simulates the standard fortran format fn.m
  !r                 Here n follows D, ndec the role of m.  Or:
  !r     'Fn'        fixed field, picking between d and e that generates
  !r                 the most significant digits
  !r      (n.m)      also simulates the standard fortran format.
  !r
  !r  *Optional modifier 'n' is a number specifying how many decimals of
  !r   precision (n=6 if not specified). By default, n means:
  !r      for 'd' format, the absolute precision: i.e.
  !r        number of digits after the decimal point
  !r     for 'e' format, the relative precision , i.e.
  !r        number of digits printed
  !r     for 'D' format, it is the field width n in fortran format fn.m
  !r  *Optional modifier sw is a compound of the 1's and 10's digits.
  !r       1's digit of sw can overwrite the default meaning of 'n' above.
  !r                 sw=0 => n corresponds to absolute precision
  !r                 sw=1 => n corresponds to relative precision
  !r       10's digit nonzero suppresses leading blanks.
  !r  *Entry bin2a0 allows the user to set the default of sw.
  !r  *Examples:
  !r     call bin2a('d2',1,3,1.234951d0,...)    => 1.23
  !r     call bin2a('d4',1,4,1.234951d0,...)    => 1.2350
  !r     call bin2a('d3:11',1,0,1.2349501d-6,4) => .00000123
  !r     call bin2a('e2',1,3,1.2499d7,...)      => 1.2e7
  !r     call bin2a('e5',1,5,1.2349510d7,...)   => 1.2350e7
  !r     call bin2a('e5:0',1,4,1.2349501d5,...) => 1.234950100e5
  !r     call bin2a('g',1,0,1.23d-5,...)        => 1.23e-5
  !r     call bin2a('g3:10',1,3,1.24996d-5,...) => .000
  !r     call bin2a('g4:10',1,4,1.24996d-5,...) => 1e-5
  !r     call bin2a('f4:10',1,4,1.24996d-5,...) => 1e-5
  !u Updates
  !u   02 Aug 07 Added :n outputs null string when res=NULLI
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed Parameters

  !c kino's correctio for ifort was
  !c     character(mxlen):: outstr ! ?---> !character(*) can not check size of outstr.
  !c However, because of a bug in grortran4.3.4, this is not allowed. Thus I now use character(*).
  !c Sep2010
  character(*):: outstr
  !c
  character(*):: fmt
  double precision :: res(0:*)
  integer :: nblk,cast,count,ip,mxlen,ndec,is
  ! Local Variables
  logical :: lD,lS,lF,lnull,llnull,parstr
  integer :: i,j,k,iprint,lsmx,n1,n2,np,precsn,fw, &
       ix(4),iv(4),a2vec,p,isw,isw0,getdig,m,ndig,ndige
  parameter (lsmx=80)
  character(20) :: lfmt, strn*(lsmx), strn2*(lsmx), ss*(lsmx)
  character:: fm*20
  real :: rval
  double precision :: xx
  integer :: NULLI
  parameter (NULLI=-99999)
  save isw0
  data isw0 /0/

  !     write(*,"('enter bin2a: cast,fmt=',i4,1x,a$)") cast,fmt

  ! --- Convert binary to ascii representation (log, int, char) ---
  lnull = .false.
  llnull = .false.
  goto (10,11,12,20,20), cast+1
  call rx('bin2a: bad cast')
10 continue
  call bin2al('(L8)',res,count,strn2)
  goto 15
11 continue
  strn2 = fmt
  goto 15
12 continue
  call bin2ai('(I16)',res,count,strn2,lnull)
  goto 15
  ! --- copy strn2 to strn with appropriate number of blanks ---
15 continue
  i = 0
  call skipbl(strn2,lsmx,i)
  strn = ' '
  !     If a field width specified, overwrite spillover with '*'
  if (ndec /= 0) then
     call skpblb(strn2,lsmx,j)
     j = j-ndec+1
     if (j > i) then
        strn2(j+1:j+ndec) = '****************'
     endif
     i  = j
  endif
  strn(1+nblk:lsmx) = strn2(i+1:lsmx)
  call skpblb(strn,lsmx,n1)
  n1 = n1+1
  if (lnull .AND. fmt /= ' ') then
     i = 0
     if (parstr(fmt,':n',len(fmt)-1,2,'n',i,j)) then
        llnull = .true.
     endif
  endif
  goto 50

  ! --- Entry for setting up or determinining defaults ---
  entry bin2a0(is)
  if (is >= 0) isw0 = is
  if (is < 0) is = isw0
  return

  ! --- Binary->ascii representation, floating-point ---
20 continue
  if (cast == 3) xx = rval(res,count+1)
  if (cast == 4) xx = res(count)
  lnull = xx .eq. dble(NULLI)

  ! ... Determine appropriate format
  lfmt = fmt
  i = 0
  call skipbl(fmt,len(fmt),i)
  if (i >= len(fmt)) then
     lfmt = 'g'
  else
     lfmt = fmt(i+1:len(fmt))
  endif
  i = 0
  if (parstr(lfmt,':n',len(lfmt)-1,2,'n',i,j)) then
     lfmt(i+1:) = ' '
     llnull = .true.
  endif
  ! --- Do the conversion, floating point ---
  if (lfmt(1:1) == '(') then
     write(ss,lfmt) xx
     call pretty(ss,nblk,ndec,20,isw0,strn,n1)
  else
     strn  = ' '
     strn2 = ' '
     lD = .false.
     lF = .false.
     j = 0
     !   ... i=1 => 'd'  i=2 =>  'e'  i=3 => 'g'
     call chrps2(lfmt,'degDF',5,len(lfmt),j,i)
     if (i <= 0) call rx('bin2a: bad format: '//lfmt)
     if (i == 5) then
        i = 3
        lF = .true.
     elseif (i == 4) then
        i = 1
        lD = .true.
     endif
     !   ... Get precsn (or field width for D or F), in iv(1), sw in iv(2)
     j = j+1
     np = a2vec(lfmt,len(lfmt),j,2,': ',2,2,2,ix,iv)
     isw = 1 + isw0
     if (i == 1) isw = 0 + isw0
     !   ... Simulated fortran format: precsn dictated by ndec
     if (lF) then
        if (np <= 0) call rx('bin2a: bad format: '//lfmt)
        fw = iv(1)
     elseif (lD) then
        precsn = ndec
        fw = -1
        if (np >= 1) fw = iv(1)
        !   ... if precsn explicit, use it
     elseif (np >= 1) then
        precsn = iv(1)
        !   ... This is the default
     else
        precsn = 6
     endif
     if (np >= 2) isw = iv(2)
     if (isw >= 20) isw = mod(isw,10) + isw0
     !  21   continue
     !   ... p is the exponent
     p = 0
     if (xx /= 0) then
        p = int(dlog10(dabs(xx)))
        if (dabs(xx) < 1) p = p-1
     endif
     !   ... fortran 'f' format
     if (i == 1 .OR. i == 3) then
        !     ... Estimate total width of format statement for fortran write
        if (lF) then
           !       ... m is the space consumed by a '-' sign
           m = (1-int(dsign(1d0,xx)))/2
           !       ... precsn = # rhs dec = field width - '-' - '.' - (p+1)
           precsn = fw - m - 1 - max(p+1,1)
           !       ... Only works on some compilers
           !            if (mod(isw,10) .ne. 0)
           !     .      precsn = fw - m - 1 - max(p+1,0)
           !       ... ndig = how many nonzero decimals printed
           ndig = max(precsn+p+1,0)
           !       ... Exclude 'f' if it doesn't fit
           if (precsn < 0) then
              ndig = -1
              !       ... Exclude 'e' if it does, and number isn't small
           else if (p > -2) then
              i = 1
           endif
           !       ... Determine how many digits we get from 'e' format
           if (i /= 1) then
              write(ss,'(1pe20.0)') xx
              !             print *, ss
              !         ... We want at least 1 more digit than f format
              call pretty(ss,0,max(ndig+1,1),max(ndig+1,1),1,strn,j)
              !         ... Tack on trailing 'e0' if pretty discarded it
              k = 0
              call chrpos(strn,'e',j,k)
              if (k >= j) j = j+2
              !         ... How many decimals for 'e' format
              ndige = max(ndig+1,1) + fw - j
              !         ... If pretty suppresses '.', add it back if ndige>1
              if (ndige > 1) then
                 k = 0
                 call chrpos(strn,'.',j,k)
                 if (k >= j) ndige=ndige-1
              endif
              !             print *, strn
           else
              ndige = ndig-1
           endif
           !       ... Generate string for F format here.
           if (ndig < 0 .AND. ndige < 0) then
              strn = ' '
              strn(nblk+1:nblk+fw) = '********************************'
              n1 = fw+nblk
              goto 50
           else if (ndig >= ndige) then
              i = 1
           else
              i = 2
              precsn = ndige
              goto 35
           endif
        elseif ( .NOT. lD .OR. (lD .AND. fw == -1)) then
           fw = max(p+3,5) + precsn
           if (getdig(isw,0,10) == 1) &
                fw = max(p+3,3) + max(precsn-p-1,0)
           fw = max(fw,10)
           if (fw > min(lsmx-2,99)) then
              strn = ' '
              strn(nblk+1:nblk+1) = '*'
              n1 = 1+nblk
              goto 35
           endif
        endif
        j = fw
        !     ... Insert leading blanks
        !         if (lF) then
        if (lF .OR. lD) then
           j = j+nblk
        endif
        if (j >= 10) write(fm,'(''(f'',i2,''.'')') j
        if (j < 10) write(fm,'(''( f'',i1,''.'')') j
        k = j
        !     ... Number of decimals for fortran write
        j = precsn
        if ( .NOT. (lD .OR. lF)) then
           if (getdig(isw,0,10) == 1) j = precsn-p-1
           j = max(j,0)
        endif
        !         decimals can't exceed field width - 1
        j = max(min(k-1,j),0)
        if (j >= 10) write(fm(6:8),'(i2,'')'')') j
        if (j < 10) write(fm(6:7),'(i1,'')'')') j
        write(ss,fm) xx
        if (lD .OR. lF) then
           if (nblk <= 0) then
           elseif (ss(1:nblk) /= ' ') then
              ss(1:k) = '*****************************************************'
              ss(1:nblk) = ' '
           endif
           strn = ss
           call skpblb(strn,lsmx,n1)
           n1 = n1+1

           k = 0
           call chrps2(strn,'-.0123456789',12,n1,k,j)
           j = j-1
           lS = j .eq. 0
           if (lS) call chrps2(strn,'.0123456789',11,n1,k,j)
           !     ...   Case fraction should have a leading '0'
           if (j == 1 .AND. getdig(isw,1,10) == 0) then
              if (lS .AND. k > 1) strn(k-1:k) = '-0'
              if ( .NOT. lS .AND. k > 0) strn(k:k) = '0'
              !     ...   Case fraction should have no leading '0'
           elseif (j == 2 .AND. getdig(isw,1,10) /= 0) then
              if (lS) strn(k:k+1) = ' -'
              if ( .NOT. lS)strn(k+1:k+1) = ' '
           endif
        else
           !           print *, 'before pretty ...', fm, ss
           call pretty(ss,nblk,ndec,precsn,isw,strn,n1)
        endif
35      continue
     endif
     !    .. fortran 'e' format
     if (i == 2 .OR. i == 3) then
        j = p + precsn
        if (getdig(isw,0,10) == 1) j = precsn-1
        if (j > 22) then
           strn2 = ' '
           strn2(nblk+1:nblk+1) = '*'
           n2 = 1+nblk
           goto 45
        endif
        j = min(max(j,0),99)
        if (j >= 10) write(fm,'(''(1pe30.'',i2,'')'')') j
        if (j < 10) write(fm,'(''(1pe30.'',i1,'')'')') j
        write(ss,fm) xx
        !         print *, 'before pretty ...', fm, ss
        j = ndec
        if (lF) j = precsn
        call pretty(ss,nblk,j,precsn,isw,strn2,n2)
        !     ... Tack on trailing 'e0' if pretty discarded it
        j = 0
        call chrpos(strn2,'e',n2,j)
        if (j >= n2 .AND. i == 2) then
           strn2(n2+1:n2+2) = 'e0'
           n2 = n2+2
        endif
        !     ... Sometimes the '.' is suppressed; make fw right
        if (lF .AND. n2 < fw+nblk) n2 = fw+nblk
45      continue
     endif
     if (i == 2 .OR. i == 3 .AND. &
          (n2 < n1 .OR. strn(nblk+1:nblk+1) == '*')) then
        strn = strn2
        n1 = n2
     endif
  endif
  ! --- Copy to outstr ---
50 continue
  n1 = max(n1,0)
  n2 = max(min(n1,mxlen-ip),0)
  !     Handle null number: replace ascii string with 'NULL'
  if (lnull .AND. llnull) then
     strn(1:n2) = ' '
     i = max(n2-3,1+nblk)
     strn(i:n2) = 'NULL'
  endif
  if (n2 > 0) outstr(ip+1:ip+n2) = strn(1:n2)
  ip = ip+n2
  if (ip == mxlen .AND. n2 < n1) outstr(ip:ip) = '|'
  if (iprint() > 120) print '(1x,a,a)', 'bin2a:',outstr(1:ip)
end subroutine bin2a

subroutine bin2al(fmt,res,count,strn)
  character*(*) fmt, strn
  integer :: count
  logical res(0:*)
  write(strn,fmt) res(count)
end subroutine bin2al

subroutine bin2ai(fmt,res,count,strn,lnull)
  !- Conversion of integer to ascii string
  ! ----------------------------------------------------------------------
  !i Inputs
  !l         :
  !i   fmt   : fortran format
  !i   res   : res(count) is converted
  !i   count : index to res: res(count) is converted
  !o Outputs
  !o   strn  : ascii representation of integer
  !i   lnull : true if res(count)=NULLI, otherwise false
  ! ----------------------------------------------------------------------
  !     implicit none
  character*(*) fmt, strn
  integer :: count
  integer :: res(0:count)
  logical :: lnull
  integer :: NULLI
  parameter (NULLI=-99999)
  write(strn,fmt) res(count)
  lnull = res(count) .eq. NULLI
end subroutine bin2ai

subroutine nlchar(ich,ps)
  integer :: ich
  character(*) ps
  if (ich==1) then
     ps(1:1) =  char(10)
  endif
end subroutine nlchar

integer function getdig(n,i,base)
  !- Extracts one digit from an integer
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n,i,base
  !o Outputs
  !o   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
  ! ----------------------------------------------------------------
  implicit none
  integer :: n,i,base
  getdig = mod(n/base**i,base)
end function getdig

subroutine pretty(sin,nblk,ndec,precsn,sw,sout,nout)
  !- Prettifies an ascii representation of a floating-point number
  !i  nblk: number of leading blanks in output string
  !i  ndec: minimum number of decimals to display (see sw)
  !i  precsn: truncate to this many decimals of precision (see sw)
  !i  sw: 1's digit: 0 for ndec and precsn absolute (number of places
  !i                   to the right of '.').
  !i                 1 for ndec and precsn number of significant digits.
  !i     10's digit: 1 to suppress possible leading zero
  !o  sout: reformatted ascii representation of number
  !o  nout: length of sout
  !r  ndec and precsn must be >= 0
  !r  If ndec>precsn, ndec truncated to precsn (precsn>0)

  !     implicit none
  !      character*(*) sin,sout*40
  character*(*) sin,sout
  integer :: ndec,precsn,nout,sw,nblk
  integer :: ls,is,iex,ixp,i,j,k,m,i0,i1,ip,nd,im,id,lsin,getdig
  integer :: ndig,idec
  logical :: a2bin,ltmp,isdig
  character ss*40, sexp*5, ch*1, sss*40, fmt*8
  double precision :: xx,dround


  ! --- Copy input string to local ---
  !     print *, 'sin=',sin
  ss = ' '
  ss(2:39) = sin
  sout = ' '
  nout = 1+nblk
  i0 = nout
  lsin = len(sin)

  ! --- Setup, normalize variations of input string ---
  ! ... i1=pos of first char, ip = pos of decimal point
1 continue
  i1 = 0
  ls = len(ss)
  call skipbl(ss,ls,i1)
  call skpblb(ss,ls,is)
  i1 = i1+1
  ! ... Remove leading '+'
  if (ss(i1:i1) == '+') i1 = i1+1
  if (ss(i1:i1) == '-') then
     sout(nout:nout) = '-'
     i1 = i1+1
     nout = nout+1
     i0 = nout
  endif
  ! ... Check for exceptions
  if (ss(i1:i1) == '*') then
     sout(nout:nout) = '*'
     return
  endif
  ! ... Check for 'NaN or Infinity'
  if (ss(i1:ls) == 'nan' .OR. ss(i1:ls) == 'NaN' .OR. &
       ss(i1:ls) == 'NAN' .OR. ss(i1:ls) == 'QNAN') then
     sout(nout:nout+2) = 'NaN'
     nout = nout+2
     return
  endif
  if (ss(i1:ls) == 'inf' .OR. ss(i1:ls) == 'Inf' .OR. &
       ss(i1:ls) == 'INF' .OR. ss(i1:ls) == 'Infinity') then
     sout(nout:nout+2) = 'Inf'
     nout = nout+2
     return
  endif

  is = is+1
  ! ... Location of decimal point (ip); prepend 0 if missing
  ip = 0
  call chrpos(ss,'.',is,ip)
  ip = ip+1
  if (ip > is) then
     is = is+1
     ss(is:is) = '.'
  endif
  if (ip == i1) then
     i1 = i1-1
     ss(i1:i1) = '0'
  endif

  ! --- Find end of mantissa; get exponent (binary in ixp) ---
  iex = 0
  im = ip-1
  ixp = 0
  call chrps2(ss,'dDeE+-',6,is,im,j)
  if (j /= 0) then
     iex = im+2
     if (j >= 5) iex = iex-1
     if (ss(iex:iex) == '+' .OR. ss(iex:iex) == '-') iex = iex+1
     do  10  j = iex, is
        k = j
        if (ss(j:j) /= '0') goto 12
10   enddo
12   continue
     if (ss(iex-1:iex-1) == '-') then
        k = k-1
        ss(k:k) = '-'
     endif
     iex = k
     k = iex-1
     if ( .NOT. a2bin(ss,ixp,2,0,' ',k,-1)) &
          call rx('pretty: parse error')
     !        print *, iex,is, ':', ss(1:iex-1), '|', ss(iex:is), '|', ixp
  endif

  ! --- Truncate string at precsn decimal places, possibly rounding ---
  if (precsn > 0) then
     k = ip+precsn+ixp
     m = i1
     ! ...   If relative precision, k offset from first significant digit
     if (getdig(sw,0,10) == 1) then
        do  14  i = i1, im
           m = i
           if (ss(i:i) >= '1' .AND. ss(i:i) <= '9') goto 15
14      enddo
15      continue
        k = m-1+precsn
        if (m <= ip .AND. k >= ip) k = k+1
     endif
     k = max(min(k,im),m)
     !        print *, k, ' ##|', ss(i1:k), '|'
     ! ...   Round upwards, append exponent and start over
     if (ss(k+1:k+1) >= '5' .AND. ss(k+1:k+1) <= '9') then
        im = min(k,im)
        sss = ss(i1:max(im+1,ip))
        read(sss,'(e30.20)') xx
        fmt = '(f30.'
        if (im-ip >= 10) write(fmt(6:8),'(i2,'')'')') im-ip
        if (im-ip < 10) write(fmt(6:7),'(i1,'')'')') max(im-ip,0)
        ! ...     Next line makes some recalcitrant machines round properly
        if (im > ip) xx = xx + 4d0/10d0**(1+im-ip)
        if (im > ip) write(sss,fmt) xx
        if (im <= ip) write(sss,fmt) dround(xx,k+1-i1)
        !          print *, sss, im-ip
        ! ...     Append exponent
        if (j > 0) then
           im = 0
           call chrps2(sin,'0123456789',10,lsin,im,j)
           call chrps2(sin,'dDeE+-',6,lsin,im,j)
           call skpblb(sss,40,k)
           call strcop(sss(k+2:40),sin(im+1:lsin),lsin-im,' ',j)
           !            print *, sss, im-ip
        endif
        ! ...     Start over
        ss = sss
        goto 1
     endif
     if (k < ip) then
        do   j = k+1, ip-1
           ss(j:j) = '0'
        enddo
     endif
     !        print *, precsn, '(precision):', ss(1:im), '|', k,im
     im = max(min(k,im),ip)
  endif

  ! --- Patch 0.mmmEnnn so that leading digit is nonzero ---
  if (iex /= 0 .AND. ixp /= 0) then
     if (ss(i1:i1+1) == '0.') then
        ixp = ixp-1
        i1 = i1+1
        ss(i1:i1) = ss(i1+1:i1+1)
        ss(i1+1:i1+1) = '.'
        ip = i1+1
        !          print *, 'patch', ss(i1:is)
     endif
     sexp = ' '
     write(sexp,'(i5)') ixp
     !        print *, 'exponent', sexp
  endif

  ! --- Copy mantissa to sout, counting decimals ---
  idec = 0
  ndig = 0
  ltmp = getdig(sw,1,10) .ne. 0
  do  20  j = i1, im
     isdig = ss(j:j) .ge. '0' .and. ss(j:j) .le. '9'
     if (j == i1 .AND. ltmp .AND. ss(j:j) == '0') goto 20
     if (j == ip .OR. isdig) then
        sout(nout:nout) = ss(j:j)
        nout = nout+1
        if (isdig .AND. j > i1) ndig = ndig+1
        if (isdig .AND. j > ip) idec = idec+1
     else
        goto 22
     endif
20 enddo
22 continue
  nout = nout-1
  idec = idec-ixp
  !      print *, 'string before truncation:', sout(1:nout),
  !     .  '...', idec, ' decimals', ndig, ' digits; ndec=',ndec

  ! --- Append or strip trailing zeros, preserving ndec decimals ---
  ip = 0
  call chrpos(sout,'.',nout,ip)
  ip = ip+1
  if (ip > nout) then
     ss = sin
     i1 = 0
     ls = len(ss)
     call skipbl(ss,ls,i1)
     call skpblb(ss,ls,is)
     call rx('pretty:  missing ''.'' in string'//ss(i1:is+1))
  endif
  ! ... Find i1 which points to 1st significant digit
  i1 = 0
  call chrps2(sout,'123456789',9,nout,i1,j)
  i1 = i1+1
  ! ... Case number is zero
  if (i1 > nout) then
     i1 = 0
     call chrpos(sout,'0',nout,i1)
     i1 = i1+1
     if (i1 > nout) call rx('pretty:  missing digit')
  endif
  id = 0
  nd = ndec
  if (precsn > 0) nd = min(nd,precsn)
  if (ndec > 0) then
     id = ip+nd+ixp
     if (getdig(sw,0,10) == 1) id = i1+nd-1
     if (getdig(sw,0,10) == 1 .AND. i1 < ip) id = id+1
  endif
  ! ... Append
  if (nout < id) then
     do   k = nout+1, id
        sout(k:k) = '0'
     enddo
     nout = id
  endif
  ! ... Strip
  if (0<id .AND. id<=len(sout))then
     if (sout(id:id) == '.') id = id-1
  endif
  do  41  k = nout, i0, -1
     j = k
     !        print *, '!!!',sout(1:j),' ',j,id
     ch = sout(j:j)
     if (j <= id) goto 42
     if (ch == '0' .OR. ch == '.') sout(j:j) = ' '
     if (ch /= '0' .OR. ch == '.') goto 42
41 enddo
42 continue
  ! ... Fix up a completely missing mantissa
  nout = j
  if (nout == i0 .AND. sout(i0:i0) == ' ') then
     sout(i0:i0) = '0'
     return
  endif
  if (sout(nout:nout) == ' ') nout = nout-1

  ! --- Append exponent if there is one ---
  if (ixp /= 0) then
     j = 0
     call skipbl(sexp,5,j)
     !        print *, j,':',sexp(j+1:5)
     sout(nout+1:nout+7-j) = 'e' // sexp(j+1:5)
     nout = nout+6-j
  endif

  !      print *, 'final string:', sout(1:nout)
  !      print *, '--------------'

end subroutine pretty

double precision function dround(x,n)
  !- Rounds double precision x after n digits
  !     implicit none
  integer :: n,is,i
  double precision :: x,s
! #if QUAD
!   double precision :: xnint
! #endif
! #if QUAD
!   s = qlog10(dble(abs(x)))
! #else
  dround=1d99
  s = dlog10(dabs(x))
!#endif
  is = s
  if (is > 0) then
     is = -int(is) + n-1
  else
     is = int(-is) + n
  endif
  s = 1
  do    i = 1, iabs(is)
     s = 10*s
  enddo
!#if QUAD
!  if (is > 0) dround = xnint(s*dble(x))/s
!  if (is < 0) dround = xnint(dble(x)/s)*s
!#else
  if (is > 0) dround = dnint(s*x)/s
  if (is < 0) dround = dnint(x/s)*s
!#endif
END function dround

logical function a2bina(instr,res,cast,count,term,j,jmaxi)
  !- Convert ASCII to logical, integer, real, double, with possible assignment
  ! ----------------------------------------------------------------
  !i Inputs
  !i   instr: contains string to be converted
  !i   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
  !i   j:     offset to first character in string to read
  !i   term:  character that terminates the expression
  !i   jmaxi: j is not to exceed jmax
  !i          Also, jmaxi should not exceed (length-of-instr)-1
  !o Outputs
  !o   count'th element of array res is set to converted value
  !o   j is put past last character read
  !r Remarks
  !r   a2bina is identical to the standard expression parser, a2bin, with
  !r   the addition that a valid expression may include assignment
  !r   operators.  The following adds 2 to variable x, and returns value 2
  !r       x+=2
  !r  Sequences of expressions separated by commas, are permitted, viz:
  !r       y=6,x=2,y-=3,y/x
  !r  assigns y to 6, x to 2, subtracts 3 from y, and returns 1.5
  !r
  !r  Allowed assignment operators  'op'  for syntax "var op expr"  are
  !r   var op expr  variable assignment
  !r     =  assign expr to variable
  !r    *=  multiply expr into variable
  !r    /=  divide expr into variable
  !r    +=  add expr into variable
  !r    -=  subtract expr from variable
  !u Updates
  !u   27 Sep 04 Accomodate macro calls
  ! --------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character(1) :: instr(0:*),term
  integer :: cast,count,j,jmaxi
  double precision :: res(0:1)
  ! ... Local parameters
  integer :: tj,jmax,nextop,i0,i,k
  character(1) :: aops(9), aterm(2), ctmp, a(80), lterm
  logical :: lassgn,a2bin,nxtexp,namchr,lmac
  double precision :: dum,dum2,v2dbl
  ! ... Allowed characters in a name
  namchr(ctmp) =ctmp .ge. 'a' .and. ctmp .le. 'z' .or. ctmp .eq. '_' &
       .or. ctmp .ge. 'A' .and. ctmp .le. 'Z' &
       .or. ctmp .ge. '0' .and. ctmp .le. '9'
  data aops/' ','=','*','/','+','-','^',' ',','/
  data aterm/' ',','/

  aops(1) = term
  aterm(1) = term
  a2bina = .false.
  tj = j
  jmax = jmaxi
  if (jmax < 0) jmax = 999999

  ! ... re-entry point for multiple expressions
10 continue
  lmac = .false.
  call skipbl(instr,jmax,tj)

  ! ... Check for assignment statement.  1st char must be allowed name char
  if (tj > jmax) return
  j = tj
  lassgn = &
       namchr(instr(j)) .and. (instr(j) .lt. '0' .or. instr(j) .gt. '9')
  !     Look for first non-word character
  if (lassgn) then
     do  12  tj = j, jmax
        a(tj-j+1:tj-j+2) = instr(tj) // ' '
        if ( .NOT. namchr(instr(tj))) goto 14
12   enddo
     !     Reached last character; no non-word character found
     lassgn = .false.
14   continue

     ! ... If 1st word is a macro, statement is not an assignment
     if (lassgn) then
        k = 0
        call macevl(a,' ',k)
        if (k > 0) lassgn = .FALSE. 
        if (k > 0) lmac = .TRUE. 
     endif

     ! ... Next check for assignment statement. Inspect first non-word char
     if (lassgn) then
        call skipbl(instr,jmax,tj)
        call chrps2(instr,aops,9,tj,tj,nextop)
        !     Set lassgn=.true. if some assignment operator was found
        lassgn = .false.
        if (tj < jmax .AND. nextop > 1 .AND. nextop < 8) then
           i0 = tj
           lassgn = nextop .gt. 2 .and. instr(tj+1) .eq. '='
           if (lassgn) i0 = tj+1
           lassgn = lassgn .or. nextop .eq. 2 .and. instr(tj+1) .ne. '='
           if (lassgn) then
              a = ' '
              call strncp(a,instr(j),1,1,tj-j)
              j = i0+1
           endif
        endif
     endif
  endif
  !     End of checking for assignment statement. If true, lassgn=.true.
  !     j is updated to first character in expression

  ! ... See if there is an expression following this one
  nxtexp = .false.
  lterm = term
  tj = j
  if ( .NOT. lmac) then
     call chrps2(instr,aterm,2,jmax,tj,i)
     if (i == 2)  then
        lterm = aterm(2)
        nxtexp = .true.
     endif
  endif
  a2bina = a2bin(instr,res,cast,count,lterm,j,jmaxi)
  if ( .NOT. a2bina) return

  if (lassgn) then
     dum = v2dbl(res,res,res,res,cast,count+1)
     dum2 = 0
     call getsyv(a,dum2,i)
     if (i == 0) call addsyv(a,0d0,i)
     if (nextop-1 == 2) dum = dum2*dum
     if (nextop-1 == 3) dum = dum2/dum
     if (nextop-1 == 4) dum = dum2+dum
     if (nextop-1 == 5) dum = dum2-dum
     if (nextop-1 == 6) dum = dum2**dum
     call chsyv(a,dum,i)
     !       call shosyv(0,0,0,6)
  endif

  tj = j
  if (nxtexp) goto 10
end function a2bina
real(8) function v2dbl(resL,resI,resR,resd,cast,n)
  !- Converts a number of different casts into a double
  !i   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
  !     implicit none
  integer :: cast,n
  logical :: resL(n)
  integer :: resI(n)
  real :: resR(n)
  double precision :: resD(n)

  if (cast == 0) then
     v2dbl = 0
     if (resL(n)) v2dbl = 1
  elseif (cast == 2) then
     v2dbl = resI(n)
  elseif (cast == 3) then
     v2dbl = dble(resR(n))
  elseif (cast == 4) then
     v2dbl = resD(n)
  else
     call rx('v2dbl: bad cast')
  endif
end function v2dbl


  subroutine mkilst(strn,nlist,list)
    !- Resolve list (ascii string) into a vector of integers
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   strn  :string holding list of integers
    !o Outputs
    !o   nlist :number of integers in list
    !o         :nlist<0 => mkilst failed to parse list
    !o   list  :list of integers
    !r Remarks
    !r   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
    !r   low:high:step
    !r   low, high, and step are integer expressions specifying the sequence
    !r     low, low+step, low+2*step, ... high.
    !r   If :step is missing, the step size defaults to 1.  If also :high
    !r   is missing,  the sequence reduces to a single integer. Thus,
    !r     '5+1'       becomes a single number, 6.
    !r     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
    !r     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
    !r   Sequences may be strung together separated by commas, eg
    !r     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
    !u Updates
    !u   02 Feb 01 strn is now a character string
    ! ----------------------------------------------------------------------
    !     implicit none
    integer :: list(*),nlist
    character*(*) strn
    integer :: it(512),iv(512),a2vec,ip,i,j,k
    ip = 0
    nlist = -1
    call skipbl(strn,len(strn),ip)
    k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
    if (k < 1) return
    if (k >= 99) call rx('mkilst: increase size of iv')
    it(k+1) = 0
    iv(k+1) = iv(k)
    ! ... loop over all iv
    nlist = 0
    i = 0
14  continue
    i = i+1
    ! ... Case iv => a single number
    if (it(i) /= 2) then
       nlist = nlist+1
       list(nlist) = iv(i)
       ! ... Case iv => n1:n2:n3
    elseif (it(i+1) == 2) then
       do  j = iv(i), iv(i+1), iv(i+2)
          nlist = nlist+1
          list(nlist) = j
       enddo
       i = i+2
       ! ... Case iv => n1:n2
    else
       do    j = iv(i), iv(i+1)
          nlist = nlist+1
          list(nlist) = j
       enddo
       i = i+1
    endif
    if (i < k) goto 14
  end subroutine mkilst
! #if TEST
! subroutine fmain
!   implicit none
!   character(20) :: strn
!   integer :: nlist,list(20),i

!   strn = '                 2,1'
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)
!   strn = '                2,1 '
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)
!   strn = '             22:33:3'
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)

! end subroutine fmain
! #endif

