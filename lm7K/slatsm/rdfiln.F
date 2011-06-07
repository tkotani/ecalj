      subroutine rdfiln(unit,cch,mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
C- Read one line of a file, with variables substitution and loop parsing
C ----------------------------------------------------------------
Ci Inputs
Ci   unit: file logical unit number
Ci   recl: maximum record length allowed
Ci   nr:   number of records already read.
Ci         Must set nr=0 for first call to initialize variables.
Ci   cch: comment/command character
Ci   cch(1:1) comment char, eg '#': line ignored if first char=cch(1:1)
Ci   cch(2:3) expression chars, eg '{}'.  See Remarks.
Ci   cch(4:4) directive char, eg '%', signaling rdfiln directives.
Ci   cch(5:5) if 'k', save all declared variables
Ci   cch(6,6) if 'c', try to parse for command-line char variable
Ci            declarations (-c=nam,strn nam,strn ...) when nr=0
Ci   cch(7,7) if 't', replace tabs with spaces
Ci   cch(8,8) if 'p', write contents to stdout
Ci   loop0,nlin,list,ilist,nlist,vnam: arrays to keep track of loop
Ci         nesting levels
Ci   mxlev:   maximum depth allowed for nesting.
Ci   ctbl,mxchr: table of character variables.
Ci   a: character work array of length recl.
Co Outputs
Co   recrd: next record input
Co   nr is incremented by one.  Returns -nr if EOF encounterd.
Cr Remarks
Cr  *Rdfiln reads an ascii string from logical file 'unit' into character
Cr   array 'recrd', with the option of parsing algebraic expressions in
Cr   the line, as described below.  Rdfiln also supplies some facility
Cr   for branching control to skip over reading of certain lines or
Cr   repeatedly reading a block of lines.  The caller calls rdfiln until
Cr   all the lines desired are read in, or until EOF is reached.
Cr
Cr   -------------- Conceptual overview --------------
Cr   For concretness, let's assume the following
Cr     cch(1) = '#' comment character
Cr     cch(2) = '{' character marking start of expression substitution
Cr     cch(3) = '}' character marking end of expression substitution
Cr     cch(4) = '%' rfiln directive character
Cr
Cr  *Expression substitution
Cr   rdfiln parses anything inside `{...}' and substitutes the contents for
Cr   `something else'.  For example, the contents of `{...}' may by an
Cr   algebraic expression.  In that case, rdfiln evaluates the expression
Cr   numerically, and turns it back into a number, in ascii form. Thus
Cr   a string read from 'unit' as
Cr      talk {4/2} me
Cr   becomes
Cr      talk 2 me
Cr   rdfiln evaluates '4/2' as a floating-point expression using a2bin.f,
Cr   turns the result back into an ascii string, using bin2a.f, and
Cr   substitutes the resultant substring for the original.
Cr   Caution: if the substituted string is longer than recl, it is truncated.
Cr
Cr  *Variables
Cr  *Rdfiln permits three kinds of variables, floating point scalar,
Cr   floating-point vector, and character.  It uses a2bin to evaluate
Cr   floating-point expressions and bin2a to recast the result as a string.
Cr   The scalar symbols table is maintained in the standard variables table,
Cr   File symvar.f contains the source code maintaining this table.; file
Cr   symvec.f contains source code maintaining the vector variables table.
Cr   (NB: symvec allocates space for these vectors through the malloc utility;
Cr   your compiler must have pointer capability to use this table.)  The table
Cr   of character variables is maintained in the character array ctbl, which is
Cr   passed as an argument to rdfiln.
Cr
Cr   NB: when EOF is reached, rdfiln eliminates any variables declared within
Cr   the input file using `%' constructs described below (see % save directive
Cr   below for exceptions to this rule).
Cr
Cr  *Comments
Cr   Lines beginning with the comment character are purged.
Cr
Cr  *Rdfiln directives
Cr   Lines beginning with '% directive', where directive is one of:
Cr     const cconst cvar udef var vec char char0 cchar getenv
Cr     if ifdef ifndef iffile else elseif elseifd endif include includo
Cr     while repeat end
Cr     echo show stop exit save trace vfind
Cr   are interpreted by rdfiln not as part of the input, but as a directive to
Cr   do something, such as assign a value to a variable; to conditionally skip
Cr   over a block of lines; to repeatedly read a block of lines using a 'loop'
Cr   construct; and some miscellaneous commands.  Each of these is described
Cr   in the 'Rdfiln directives' section below.
Cr
Cr   -------------- Expression substitution --------------
Cr    A string in `{}', e.g. `{strn}' may contain one of the following.
Cr    These are listed in order of precedence in which rdfiln tries to
Cr    parse the contents of `{...}' .
Cr    NB: The {} can be nested.
Cr
Cr      1. the name of a character variable, say `myvar'
Cr         In that case, rdfiln replaces string `{myvar}' with contents of
Cr         `myvar'.
Cr
Cr      2. a character variable name, say `myvar', followed by a qualifier (...)
Cr         which can be one of the following:
Cr
Cr        *(integer1,integer2) --- returns a substring of `myvar'
Cr         {myvar(n1,n2)} is replaced by the (n1,n2) substring of `myvar'.
Cr
Cr        *('char-list',n) --- marks a position in contents of `myvar'
Cr         {myvar('char-list',n)} is replaced by integer, which is the
Cr         index to the n'th occurence of one element in 'char-list'.
Cr         n is optional.
Cr         Example: If myvar="foo bar", {myvar('abc',2)} evaluates to 6.
Cr
Cr        *(:e) --- returns an integer marking last nonblank character
Cr         Example: If myvar='foo bar', {myvar(:e)} evaluates to 7.
Cr
Cr        *(/'str1'/'str2'/,n1,n2) --- string substitution
Cr         {myvar(/'str1'/'str2'/,n1,n2)} substitutes str2 for str1
Cr         It does it for the n1'th to n2'th occurence.
Cr         n1 and n2 are optional, as are the quotation marks.
Cr         Example: If myvar="foo boor", {myvar(/'oo'/a/,2,2)} = "foo bar"
Cr
Cr      3. The name of a vector variable, say `myvec'
Cr         rdfiln replaces '{myvec}' with a sequence of numbers each separated
Cr         by one space, which are the contents of myvec'  Thus
Cr           % vec myvec[5] 5 4 3 2 1
Cr           {myvec}
Cr         becomes
Cr           5 4 3 2 1
Cr         (The first line declares `myvec' to be vector of length 5 and
Cr         initializes its contents; see description of % vec below)
Cr         Alternatively you can substitute a single element.  Thus
Cr           {myvec(2)}
Cr         is transformed into
Cr           4
Cr
Cr      4. a string consisting an algebraic expression of scalar numbers and
Cr         previously declared variables.  (variables are declared and set with
Cr         '%' directives; see 'Rdfiln directives' section below.)  rdfiln
Cr         parses the expression, turns the result into a string, and
Cr         substitutes the string in place of {expression}.  This is a special
Cr         case of the following:
Cr
Cr      5. A variable assignment, or a sequence of assignments separated by
Cr         commas.  This syntax returns the value of the (last) expression,
Cr         while assigning variables to evaluated expressions.
Cr         NB: the last expression need not have an assignment operator
Cr         {x=3}               ->  is replaced by '3'
Cr         {x=3,y=4}           ->  is replaced by '4'
Cr         {x=3,y=4,x*=y}      ->  is replaced by '4'
Cr         {x=3,y=4,x*=y,x*2}  ->  is replaced by '24'
Cr
Cr         The general syntax is:  {var assignment-op expr [, ... ]} .
Cr         The following are the allowed operators:
Cr         assignment-op         function
Cr           '='            simple assignment
Cr           '*='           replace 'var' by var*expr
Cr           '/='           replace 'var' by var/expr
Cr           '+='           replace 'var' by var+expr
Cr           '-='           replace 'var' by var-expr
Cr           '^-'           replace 'var' by var^expr
Cr         NB: it is permissible to omit the 'var assignment-op' pair;
Cr         may be be useful for the final expression, as in the last example.
Cr
Cr      6. A C-like syntax of the form '{?~expr~strn1~strn2}'
Cr         If expr evaluates to nonzero, the {...} is replaced by strn1
Cr         If expr evaluates to zero, the {...} is replaced by strn2
Cr         NB:  the '~' above can be any character
Cr
Cr   To summarize, emphasizing the order of precedence: rdfiln first looks to
Cr   see if the contents of `{...}' is the name of a character variable, or a
Cr   name followed by qualifiation (...).  If so, `{...}' is replaced by the
Cr   (possibly qualified) value of the variable.  If not, rdfiln sees whether
Cr   the contents of `{...}' is the name of a vector variable.  If so, rdfiln
Cr   substitutes `{...}' into a character representation of the vector as
Cr   described in step 3 above.  If this fails, rdfiln parses `{...}' as a
Cr   scalar expression, or a sequence of expressions, and `{..}' is replaced by
Cr   a character representation of the result (4 and 5 above).
Cr
Cr   Example:  suppose that the variables table looks like:
Cr     Var       Name                 Val
Cr      1        t                   1.0000
Cr      2        f                  0.00000
Cr      3        pi                  3.1416
Cr      4        a                   2.0000
Cr   ...
Cr     Vec       Name            Size   Val[1..n]
Cr      1        firstnums          5    1.0000        5.0000
Cr      2        nextnums           5    6.0000        10.000
Cr   ...
Cr       char symbol                     value
Cr      1 c                               half
Cr      2 a                               whole
Cr      3 blank
Cr
Cr   NB: The scalar variables table always begins with predefined variables
Cr   t=1,f=0 and pi.  It is STRONGLY ADVISED that you never alter any of
Cr   these variables.
Cr
Cr   You can print out the current tables of variables with the 'show'
Cr   command; see below.  (Because the vector variables can have arbitrary
Cr   length, 'show' prints only the size of the vector and the first and
Cr   last entries.  As described in more detail below, you can create such
Cr   a variables table with the following directives:
Cr
Cr   % const a=2
Cr   % char c half a whole blank " "
Cr   % vec firstnums[5] 1 2 3 4 5
Cr   % vec nextnums[5] 6 7 8 9 10
Cr
Cr   Then rdfiln substitutes for the line
Cr    {c} of the {a} {pi} is {pi/2}
Cr   yields the following:
Cr    half of the whole 3.1415926536 is 1.5707963268
Cr
Cr   whereas the line
Cr    one quarter is {1/(nextnums(4)-5)}
Cr   becomes
Cr    one quarter is .25
Cr
Cr   The following illustrates substitution of character substrings:
Cr   % char c half a whole
Cr    To {c(1,3)}ve a cave is to make a {a(2,5)}!
Cr   becomes
Cr    To halve a cave is to make a hole!
Cr
Cr   The following line illustrates substitution of vector name
Cr    {firstnums}, I caught a hare alive, {nextnums} ...
Cr   becomes
Cr    1 2 3 4 5, I caught a hare alive, 6 7 8 9 10 ...
Cr
Cr  *Nesting of {...}.  If the contents of {...} contain an
Cr   inner block of {}, the inner block is subtituted first, as
Cr   the following illustrates.  The following line
Cr      % const xx{1{2+{3+4}1}} = 2
Cr   undergoes substitution in three passes
Cr      % const xx{1{2+71}} = 2
Cr      % const xx{173} = 2
Cr      % const xx173 = 2
Cr
Cr   This line combines nesting and '{?~expr~strn1~strn2}' syntax:
Cr      MODE={?~k~B~C}3
Cr   evaluates to, if k is nonzero
Cr      MODE=B3
Cr   or, if k is zero:
Cr      MODE=C3
Cr
Cr   -------------- Rdfiln directives --------------
Cr  This section describes the syntax for each of the directives rdfiln
Cr  understands.
Cr
Cr  *'const', 'cconst' and 'var' load or alter the variables table.
Cr   A variable 'myvar' is declared eg,  % const  myvar = expr.
Cr   'expr' may be multiplied into, divided into, added into,
Cr   subtracted from or exponentiated into an already-declared variable
Cr   using one of the following C-like syntax:
Cr     myvar*=expr  myvar/=expr  myvar+=expr  myvar-=expr  myvar^=expr
Cr
Cr   'const' and 'var' are equivalent except that, for a variable
Cr   already declared, 'const' ignores a re-declaration of the
Cr   variable (nam=val), thus preserving its original value, while
Cr   'var' alters its value; see example below.
Cr
Cr   'cconst' is a conditional 'const': the first argument following
Cr   'cconst' is an expression; declarations following the expression
Cr   are parsed only if the expression evaluates to true.
Cr
Cr   'cvar' is a conditional 'var': the first argument following
Cr   'cvar' is an expression; declarations following the expression
Cr   are parsed only if the expression evaluates to true.
Cr
Cr   .... Example: for the input file
Cr     % const a = 2 b=3 c=4 d=5
Cr     a={a} b={b} c={c} d={d}
Cr     % const a=3
Cr     % var d=-1
Cr     % const b*=2 c+=3
Cr     a={a} b={b} c={c} d={d}
Cr     % cconst b==6  b+=3 c-=3
Cr     a={a} b={b} c={c} d={d}
Cr     % cconst b==6  b+=3 c-=3
Cr     a={a} b={b} c={c} d={d}
Cr
Cr   generates four lines:
Cr     a=2 b=3 c=4 d=5
Cr     a=2 b=6 c=7 d=-1
Cr     a=2 b=9 c=4 d=-1
Cr     a=2 b=9 c=4 d=-1
Cr   'a' is unchanged from its initial declaration while 'd' changes.
Cr    The two 'cconst' show that 'b' and 'c' are altered in the first
Cr    instance, since then 'b==6' is true, while are unchanged in
Cr    the second instance, since this time 'b==6' is no longer true.
Cr
Cr  *'char' and 'cchar' load or alter the character table. Directive
Cr   % char  c half     a whole      blank
Cr   loads the character table as follows:
Cr       char symbol                     value
Cr      1 c                               half
Cr      2 a                               whole
Cr      3 blank
Cr   The last value may be a blank string. 'cchar' has the syntax
Cr   % cchar nam  expr1 str1 expr2 str2 ...
Cr   expr1 expr2 etc are algebraic expressions and 'nam' takes the
Cr   value 'str1' if expr1 evaluates to true (ie nearest integer is
Cr   nonzero), the value 'str2' if expr2 evaluates to true, etc.
Cr   Re-declaration of any previously defined variable has the effect
Cr   of changing the contents of the variable
Cr  *'char0' is the same as 'char', except re-declaration of existing
Cr   variables is ignored.
Cr  *'getenv' is the same as 'char', except the string char holds
Cr   is used as a name for an environment variable, and its value
Cr   is replaced by the value of the enviroment variable.  Thus
Cr  % getenv myhome HOME
Cr   puts the string of your home directory into variable 'myhome.'
Cr
Cr  *'vec' loads or alters elements in the table of vector variables.
Cr   % vec v[n]                  creates a vector variable of length n
Cr   % vec v[n] n1 n2 n3 ...     ditto, first elements are also set
Cr   NB: Once 'v' is already declared, elements of v may be set with
Cr   the following syntax, which sets all elements bewtween i1..i2
Cr   % vec v(i) n                or
Cr   % vec v(i1:i2)  n1 n2 ... nn
Cr   There must be exactly i2-i1+1 elements n1 ... nn.
Cr   Also, if 'v' is already declared, it is an error to re-declare it.
Cr
Cr  *'vfind' finds the entry in an array that matches a specified value:
Cr   % vfind v(i1:i2)  name match-value
Cr   parses v(i) for i=i1..i2 and sets variable 'name' to i when it
Cr   finds v(i)=match-value.  If no match, 'name' is set to zero.
Cr   .... Example, the lines
Cr   % vec  a[3] 101 2002 30003
Cr   % vfind a(1:3) i 2002    <---- will set  i to 2
Cr   % vfind a(1:3) i 10      <---- will set  i to 0
Cr
Cr  *'save' preserves some variables for future use
Cr   % save              preserves all variables defined to this point
Cr   % save name [name2 ...]                saves only variables named
Cr   NB: only scalar variables may be saved.
Cr
Cr  *'udef' deletes a variable and its definition.
Cr   Only scalar and character variables may be deleted
Cr   rdfiln aborts with error if no variable exists to 'undefine'
Cr  *'udef -f' is equivalent to 'udef' except that 'udef -f' does
Cr   nothing if an attempt is made to udefine a nonexistent variable
Cr
Cr  *'trace' when turned on, chatters about how rdfiln parses the input
Cr           Invoking 'trace' with no argument toggles whether it is
Cr           on or off.
Cr           'trace 0' turns the tracing off (the default)
Cr           'trace 1' turns the tracing to lowest level:
Cr                     all directives having to do with execution flow
Cr                     (if-else-endif, repeat/while-end)
Cr           'trace 2' prints some information about most directives.
Cr
Cr  *'echo' echos the current line to stdout (i1mach(2)).
Cr  *'stop expr msg' aborts with 'msg' if 'expr' evaluates to true
Cr  *'show' prints out various things:
Cr   % show lines       (echos each line generated to the screen until:
Cr   % show stop         is encountered)
Cr   % show vars        (prints out the state of the variables table)
Cr
Cr   Expressions are evaluated for both echo and stop before printout.
Cr
Cr  *'if expr', 'elseif expr', 'else' and 'endif' are conditional read
Cr   blocks.  Lines between these directives are read or not,
Cr   depending on the value of the expression following 'if.'  For
Cr   .... Example, the lines
Cr   % if Quartz
Cr    is clear
Cr   % elseif Ag
Cr    is bright
Cr   % else
Cr    neither is right
Cr   % endif
Cr   generate one line ' is clear', if Quartz is true, ' is bright' if
Cr   Ag is false but Quartz is true, ' neither is right' otherwise.
Cr
Cr  *ifdef is similar to if, but has a more general idea of what
Cr   constitutes an expression.  First, 'if' requires a valid
Cr   expression, while 'ifdef' treats an invalid expression (eg one
Cr   containing an undefined variable) as a valid expression evaluating
Cr   to false.  The syntax of ifdef allows several expressions :
Cr   ifdef expr1 [| expr2 | expr3 ...]
Cr   and if any of expr1, expr2, ... evaluate to true, the result
Cr   is true, whether or not preceding expressions are valid.
Cr   The spaces here are syntatically significant here, since
Cr   expr1|expr2 is only true if both expr1 and expr2 are valid
Cr   expressions, while expr1 | expr2 may be true if either is valid.
Cr   'ifdef'  allows a limited use of character variables in
Cr   expressions. Either of the following are permissible expressions:
Cr     char-variable            (T if char-variable exists, otherwise F)
Cr     char-variable=='string'  (T ifchar-variable equals "string")
Cr   .... Example: ifdef  x1==2 | atom=='Mg'
Cr     is true if scalar 'x1' is 2, or if character variable
Cr     "atom" is equal to "Mg".
Cr   Also, the 'expr' above can be groups of subexpressions of any type
Cr   just mentioned, separated by ' & '.
Cr   .... Example: ifdef  x1==2 & atom=='Mg' | x1===1
Cr     is true if scalar 'x1' is 1, or 'x1' is 2, and character
Cr     variable "atom" is equal to "Mg".
Cr  *'elseifd' is to 'elseif' as 'ifdef' is to 'if'.
Cr  'if' and/or 'ifdef' constructs may be nested to a depth of mxlev.
Cr  *'ifndef' expr ... is equivalent syntatically to
Cr   'ifdef'  expr ... followed immediately by 'else'
Cr
Cr  *iffile file-name  is another conditional read block.  Instead of
Cr   the condition being set by an expression, it it set by whether
Cr   file 'file-name' exists.
Cr
Cr  *'while' and 'end', or 'repeat' and 'end' are looping constructs,
Cr   as illustrated below.  The 'while' construct has the syntax
Cr     % while [assignment assignment ...] expression
Cr      lines here are repeatly read in while expression is true
Cr     % end
Cr   here the (optional) assignments following expression have the
Cr   same syntax and meaning as the 'const' construct.  That is,
Cr   they look like 'nam=expr' or 'nam op= expr'.  As in the 'const'
Cr   case, 'nam=expr' only has effect when the variable 'nam' has
Cr   is not in the variables table.  This is made evident in the
Cr   example below.  'repeat' has the syntax
Cr     % repeat varnam list
Cr       lines here are reread, for each value 'varnam' takes in list;
Cr       list can be an integer, eg '7' or a more complex integer list,
Cr       eg '1:3,6,2' -- see mkilst.f for the syntax of an integer list.
Cr     % end
Cr   .... Example:  note in the 'while' construct the assignment
Cr        'db=-1' is operative only the first time, while 'db+=2'
Cr        is changes db in every loop.
Cr   % const nm=-3 nn=4
Cr   % while db=-1 db+=2 db<=3
Cr   % repeat k= 2,7
Cr   this is k={k} and db={db}
Cr   {db+k+nn+nm} is db + k + nn+nm, where nn+nm={nn+nm}
Cr   % end (loop over k)
Cr   % end (loop over db)
Cr   .... is expanded into
Cr   this is k=2 and db=1
Cr   4 is db + k + nn+nm, where nn+nm=1
Cr   this is k=7 and db=1
Cr   9 is db + k + nn+nm, where nn+nm=1
Cr   this is k=2 and db=3
Cr   6 is db + k + nn+nm, where nn+nm=1
Cr   this is k=7 and db=3
Cr   11 is db + k + nn+nm, where nn+nm=1
Cr
Cr  *include file-name causes rdfiln to open file 'file-name', and
Cr   input is read from the new file until EOF, after which lines are
Cr   read from the calling file.  %include may be nested to a depth
Cr   of 10.  NB:  repeat-end and if-endif constructs MUST reside in the
Cr   same file.  'includo' is identical to 'include', except that the
Cr   rdfiln aborts if the file does not exist.
Cr   Sandwiching include directives inside constructs is permissible.
Cu Updates
Cu   04 Aug 07 Bug fix when trace used together with macros
Cu   21 Dec 05 cch(8) => write contents to stdout
Cu   12 Aug 04 Comment character redefined: characters following comment
Cu             character are replaced by blanks
Cu    3 Aug 04 Changed call to nargc with call to nargf
Cu   19 Dec 02 Added macro declarations
Cu   27 Aug 01 Added option to purge tabs from input lines
Cu   13 Aug 97 {} can be nested
Cu   17 Sep 97 rdfiln allows ifdef expr & expr ...
Cu    6 Oct 97 added % vfind
Cu   24 Nov 98 added {?~expr~string1~string2}
Cu             cchar var ... doesn't alter ctbl unless an expr is true
Cu             udef char-var works
Cu             'stop' with no arguments stops
Cu             first pass at substrings in char variables
Cu   22 Oct 99 rdfiln always loads arg[0] as char variable 'progname'
Cu             added 'exit'
Cu   16 Feb 00 implemented 'iffile' and 'udef -f',
Cu             expression substitution for cchar, and
Cu             {var assignment-op expression [, ...}, {char-var(:e)}
C ----------------------------------------------------------------
      implicit none
      integer unit,recl,nr,mxchr,mxlev,lstsiz,
     .nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),nlist(0:mxlev)
      logical rdstrn,loop0(0:mxlev)
      integer ctlen
      parameter (ctlen=120)
      character*(*) cch, recrd(0:*)*1, a(*)*1, ctbl(mxchr,2)*(ctlen)
      character vnam(mxlev)*16
C Local variables
      integer i,j,nlev,ja,jr,k,lvnam,ipr,nrd,nsyv,niff,iecho,j0,i1mach,
     .nchr,lunit,nkeepv,nelt,ival,i1,i2,nargf,ndir,ltrace,isw,nchr0,
     .mxiff,fopnx
      parameter (ndir=31,mxiff=15)
      character cc1*1,cc2*1,cex*2,dir(ndir)*7,sdir(4)*5,tab*1
C     liff  is the conditional read truth table for each nesting level
C     liff0 is like liff, but is evaluated for each nesting level
C     independently of the others.  It is needed to re-evaluate liff
C     as the nesting level changes.
      logical liff(0:mxiff),liff0(0:mxiff)
      logical lcc,lex,a2bin,noskip,notab
      logical lshow,ltmp,ltmp2,lsequ,namchr
      integer unit0(10),filstk,nrdstk(10)
      double precision val,xx
      character fname*120,ctmp*1,aa*512,asbst*512,strn*512
      character*2 aops(6)
      integer LIF,LFDEF,LLSEID,LLSEIF,LLSE,LNDIF,LEPEAT,LND,LONST,LAR,
     .LHAR,LCHAR,LHOW,LCHO,LTOP,LNCLUD,LAVE,LEC,LENV,LHAR0,
     .LCONST,LVAR,LDEF,LNCLUO,LHILE,LRACE,LFNDEF,LFIND,LXIT,LFILE,
     .LACRO
      parameter (LIF=0,LFDEF=1,LLSEID=2,LLSEIF=3,LLSE=4,LNDIF=5,
     .LEPEAT=6,LND=7,LONST=8,LAR=9,LHAR=10,LCHAR=11,LHOW=12,LCHO=13,
     .LTOP=14,LNCLUD=15,LAVE=16,LEC=17,LENV=18,LHAR0=19,LCONST=20,
     .LVAR=21,LDEF=22,LNCLUO=23,LHILE=24,LRACE=25,LFNDEF=26,LFIND=27,
     .LXIT=28,LFILE=29,LACRO=30)
c #if CSYTLETABCHR
c      parameter (tab='\t')
c #else
      parameter (tab='        ')
c #endif
      save nchr,lshow,cc1,cc2,cex,lcc,lex,notab,nkeepv,liff,liff0,
     .nlev,nrd,niff,lunit,unit0,filstk,nrdstk,ltrace
C ... Allowed characters in a name
      namchr(ctmp) =ctmp .ge. 'a' .and. ctmp .le. 'z' .or. ctmp .eq. '_'
     ..or. ctmp .ge. 'A' .and. ctmp .le. 'Z'
     ..or. ctmp .ge. '0' .and. ctmp .le. '9'
      data dir /'if','ifdef','elseifd','elseif','else','endif','repeat',
     .'end','const','var','char','cchar','show','echo','stop','include',
     .'save','vec','getenv','char0','cconst','cvar','udef','includo',
     .'while','trace','ifndef','vfind','exit','iffile','macro'/
      data sdir /'stop','all','lines','vars'/ ltrace /0/
      data aops/'= ','*=','/=','+=','-=','^='/



      call getpr(ipr)
      if (recl .gt. 512) call rx('rdfiln: increase length of aa')
      aa = ' '
C     10 Oct 02 patch to avoid bug in Intel ifc compiler
C     do i = 1, recl
C       a(i:i) = ' '
C     enddo

C --- if nr=0, initialization for rdfiln ---
      if (nr .eq. 0) then
C   ... The following guarantees basic constants loaded:
        lshow = ipr .ge. 110
        lunit = unit
        nchr = 0
        fname = '-cprogname='
        call getarf(0,fname(12:))
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
C   ... loop0 is true unless within a % repeat or while with null list
        loop0(nlev) = .true.
C   ... liff: truth table for %if blocks; niff is the nesting depth
        liff(0)  = .true.
        liff0(0)  = .true.
C   ... noskip is true unless %if evaluates to F or loop0 is F
        noskip = .true.
        niff = 0
        if (len(cch) .ge. 3) then
          cex = cch(2:3)
          lex = cex .ne. ' '
        endif
        if (len(cch) .ge. 4) then
          lcc = .true.
          cc2 = cch(4:4)
          if (cc2 .eq. ' ') lcc = .false.
        endif
        if (len(cch) .ge. 5) then
          if (cch(5:5) .eq. 'k') nkeepv = -1
        endif
C   ... Look for char variable definitions
        if (len(cch) .ge. 6 .and. cch(6:6) .eq. 'c') then
          do  5  j = 1, nargf()-1
            call getarf(j,fname)
            if (lsequ(fname,'-c',2,' ',k)) then
              k = len(fname)
              jr = 0
              i  = 0
              call chrpos(fname,' ',k,jr)
              call chrpos(fname,'=',jr,i)
              if (i .lt. jr) then
                i = 2
                call parchv(fname,k,mxchr,cex,'=',0,nchr,ctbl,i)
              endif
            endif
    5     continue
        endif
        if (len(cch) .ge. 7) then
          if (cch(7:7) .eq. 't') notab = .true.
        endif
        if (len(cch) .ge. 8) then
          if (cch(8:8) .eq. 'p') lshow = .true.
        endif
      endif

c      call shosyv(0,0,0,6)
      goto 10

C --- Entry point for new value of loop variable ---
   12 continue
      if (nlev .le. 0) call rx('bug in rdfiln ... aborting')
      ilist(nlev) = ilist(nlev)+1
      if (ilist(nlev) .le. nlist(nlev) .and. liff(niff)) then
        val = list(ilist(nlev),nlev)
        call lodsyv(vnam(nlev),1,val,i)
      endif
c      call shosyv(0,0,0,6)
C --- Entry point for new line read ---
   10 continue
      noskip = loop0(nlev) .and. liff(niff)
      if (.not. rdstrn(lunit,a,recl,.false.)) goto 99
C       12 Aug 04 Blank out remining line following comment char
      if (cc1 .ne. ' ') then
        ltmp = .false.
        do  i = 2, recl
          if (a(i) .eq. cc1) ltmp = .true.
          if (ltmp) a(i) = ' '
        enddo
      endif
      if (notab) then
        do  i = 1, recl
          if (a(i) .eq. tab) a(i) = ' '
        enddo
      endif
C       call strncp(aa,a,1,1,recl)
C       Keep track of line number
      iecho = 0
      nrd = nrd+1
      do  14  j = 0, nlev
   14 nlin(j) = nlin(j)+1
c       print *, 'nlev=',nlev,' nlin=',(nlin(j),j=1,nlev),
c     .    '  a=',(a(j),j=1,recl)
      j0 = 0

C   --- Directives to rdfiln ---
      if (lcc .and. a(1) .eq. cc2) then
        j = 1
        call skipbl(a,recl,j)
        k = j-1
        call tokmat(a(j+1),dir,ndir,7,' ',i,j,.false.)
        if (ltrace .gt. 0 .and. i .ge. 0) then
          aa = ' '
          call awrit1('#rf %i: '''//dir(i+1),aa,80,0,nrd)
        endif

        if (i .eq. LIF .or. i .eq. LFDEF .or. i .eq. LFNDEF .or.
     .  i .eq. LFILE) then
          niff = niff+1
          if (niff .gt. mxiff) call rx('increase niff in rdfiln')
          j = j+k
          liff0(niff) = .false.
          liff(niff) = .false.
          if (.not. noskip .and. ltrace .gt. 2) goto 81
          if (.not. noskip) goto 10
C       ... make substitutions for eg % ... {}...
          ja = 1
          asbst = ' '
          jr = 0
          call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja,
     .    asbst,ltrace)
C       ... Determine whether line results in T or F
          if (i .eq. LFILE) then
            i1 = j+1
            call nword(asbst,1,i1,i2)
            liff0(niff) = .false.
            if (i2 .ge. i1) then
              i = fopnx(asbst(i1:i2),172,-1,-1)
              liff0(niff) = i .eq. 1
            endif
          else
            call rdfilx(asbst,ja,i.eq.LIF,j,nrd,ctbl,mxchr,nchr,
     .      liff0(niff))
            if (i .eq. LFNDEF) liff0(niff) = .not. liff0(niff)
          endif
          liff(niff) = liff(niff-1) .and. liff0(niff)
          goto 81
        elseif (i.eq.LLSEID .or. i.eq.LLSEIF .or. i.eq.LLSE) then
          if (niff .le. 0) call fexit(-1,1,' Exit -1 rdfiln:'//
     .    ' else encountered with matching if, line %i',nrd)
          if ((i .eq. LLSEID .or. i .eq. LLSEIF) .and.
     .    .not.liff0(niff)) then
            j = j+k
C         ... make substitutions for eg % ... {}...
            ja = 1
            asbst = ' '
            jr = 0
            call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja,
     .      asbst,ltrace)
            call rdfilx(asbst,ja,i.eq.LLSEIF,j,nrd,ctbl,mxchr,nchr,
     .      liff0(niff))
            liff(niff) = liff(niff-1) .and. liff0(niff)
          else
            liff(niff) = liff(niff-1) .and. .not. liff0(niff)
          endif
          goto 81
        elseif (i .eq. LNDIF) then
          niff = niff-1
          if (niff .lt. 0) call rx
     .    ('rdfiln:  endif  encountered before  if')
          goto 81
        elseif (i .eq. LHILE .and. liff(niff)) then
          nlev = nlev+1
          if (nlev .gt. mxlev) call fexit(-1,1,' Exit -1 RDFILN:  '//
     .    'loop directives nested deeper than %i',mxlev)
C       ... Make %end backspace to this %while line
          nlin(nlev) = 1
C       ... Flag the matching %end that it matches a %while construct
          nlist(nlev) = -1
          ilist(nlev) = 0
C       ... skip while block if in a nested null-list
          loop0(nlev) = .false.
          if (.not. loop0(nlev-1)) goto 82
          j = j+k
C       ... make variable substitutions
          ja = 1
          asbst = ' '
          jr = 0
          call strncp(asbst,a,1,1,recl)
          call pvfil1(recl,recl,jr,asbst,cex,nchr,mxchr,ctbl,ja,
     .    a,ltrace)
C       ... Come back here until next string not an assignment
   13     call skipbl(a,recl,j)
          if (a(j+1) .ne. cch(2:2)) then
            j0 = j+1
C         ... Return here until char not one belonging to 'name'
   23       j0 = j0+1
            if (namchr(a(j0))) goto 23
C         ... An assignment?  If so, returns i>0, or else '='
            call tokmat(a(j0),aops,6,2,' ',i,k,.false.)
            if (i .gt. 0 .or. a(j0).eq.'=' .and. a(j0+1).ne.'=') then
              call parsyv(a,recl,1,0,j)
              call strncp(aa,a,1,1,recl)
              goto 13
            endif
          endif
C       ... Done with assignments, next follows expression
          ltmp = .false.
          call rdfilx(a,recl,.false.,j,nrd,ctbl,mxchr,nchr,ltmp)
          loop0(nlev) = ltmp .and. noskip
          goto 82
        elseif (i .eq. LEPEAT .and. liff(niff)) then
          nlev = nlev+1
          if (nlev .gt. mxlev) call fexit(-1,1,' Exit -1 RDFILN:  '//
     .    'loop directives nested deeper than %i',mxlev)
          ilist(nlev) = 0
          nlin(nlev) = 0
          nlist(nlev) = 0
          loop0(nlev) = .false.
          if (.not. loop0(nlev-1)) goto 82
          j = j+k
C       ... make variable substitutions
          ja = 1
          asbst = ' '
          jr = 0
          call strncp(asbst,a,1,1,recl)
          call pvfil1(recl,recl,jr,asbst,cex,nchr,mxchr,ctbl,ja,
     .    a,ltrace)
C           Pick up variable name
          call skipbl(a,recl,j)
          vnam(nlev) = ' '
          call strcop(vnam(nlev),a(j+1),16,'=',lvnam)
          vnam(nlev)(lvnam:lvnam) = ' '
          call strncp(strn,a,1,1,recl)
          call mkilst(strn(j+1+lvnam:),nlist(nlev),list(1,nlev))
          if (nlist(nlev) .gt. lstsiz)
     .    call fexit2(-1,1,'Exit -1 rdfiln %i: list exceeded'//
     .    ' maximum number of values allowed (%i)',nrd,lstsiz)
          if (nlist(nlev).lt.0 .and. ipr .ge. 40) call shosyv(0,0,0,6)
          if (nlist(nlev).lt.0) call rx('rdfiln: bad or null list')
          loop0(nlev) = nlist(nlev) .gt. 0 .and. loop0(nlev-1)
          lvnam = lvnam-1
          if (ltrace .gt. 0) then
            write (fname,348) nrd, vnam(nlev)(1:lvnam)
  348       format(' rdfiln', i5, ': ', a,' -> ')
            call awrit2(cc1//'%a%n:1i',fname,len(fname),
     .      -i1mach(4),nlist(nlev),list(1,nlev))
          endif
          if (ltrace .gt. 0)
     .    call awrit4('%a'' over %i values; read following lines: '
     .    //'%?#n#yes#no#%?#n#%2b(blocked by level %i)',aa,80,
     .    -i1mach(2),nlist(nlev),isw(loop0(nlev)),
     .    isw(.not.loop0(nlev-1).and.nlev.gt.1),nlev-1)
          goto 12
        elseif (i .eq. LND .and. liff(niff)) then
*           print *, 'nlev=',nlev,'  ','nlist=',(nlist(i),i=1,nlev)
*           print *, 'vnam now ', vnam(nlev)
          if (nlev .le. 0) call fexit(-1,1,' rdfiln %i (stop):'//
     .    ' end encountered before while or repeat',nrd)
C       ... if incomplete a repeat:end sequence or a while-end
          if (ilist(nlev) .lt. nlist(nlev) .or.
     .    nlist(nlev) .eq. -1 .and. loop0(nlev)) then
*             print *, 'nlev=',nlev,' backspacing by',nlin(nlev)
            k = nlin(nlev)
            do  11  i = 1, nlin(nlev)
   11       backspace(lunit)
            nrd = nrd-nlin(nlev)
            call rxx(nrd.lt.0,'rdfiln: attempt to loop across files')
            do  15  i = 0, nlev
   15       nlin(i) = nlin(i)-nlin(nlev)
C         ... Case while-end
            if (nlist(nlev) .eq. -1) then
              nlev = nlev-1
              goto 82
C         ... Case repeat-end
            else
              if (ltrace .gt. 0) call awrit2(
     .        '%a'', repeat loop over '''//vnam(nlev)//
     .        '%a'' with val=%i; reread %i lines',
     .        aa,80,-i1mach(2),list(ilist(nlev),nlev),k)
              goto 12
            endif
          endif
C       ... Case loop has finished
*            nlist(nlev) = 0
*            loop0(nlev) = .true.
          nlev = nlev-1
          goto 82
        elseif ((i .eq. LCONST .or. i .eq. LVAR) .and. noskip) then
          call numsyv(nsyv)
          j = j+k
          ltmp = noskip
          call rdfilx(a,recl,.false.,j,nrd,ctbl,mxchr,nchr,ltmp)
          if (ltmp) then
C         ... make substitutions for eg % const a=nam{xx}
            ja = 1
            aa = ' '
            call pvfil1(recl,len(aa),j,a,cex,nchr,mxchr,ctbl,ja,aa,
     .      ltrace)
            j = 0
            call parsyv(aa,ja+1,999,i-LCONST,j)
          endif
c           call shosyv(0,0,0,6)
          if (noskip) goto 83
          goto 10
        elseif ((i .eq. LONST .or. i .eq. LAR) .and. noskip) then
          call numsyv(nsyv)
          j = j+k
C       ... make substitutions for eg % const a=nam{xx}
          ja = 1
          aa = ' '
          call pvfil1(recl,len(aa),j,a,cex,nchr,mxchr,ctbl,ja,aa,
     .    ltrace)
          j = 0
          call parsyv(aa,ja+1,999,i-LONST,j)
C            call shosyv(0,0,0,6)
C            pause
          if (noskip) goto 83
          goto 10
        elseif ((i .eq. LEC .or. i .eq. LFIND) .and. noskip) then
          call numsvv(nsyv)
          j = j+k
C       ... make substitutions for eg % vec a{xx}...
          ja = 1
          aa = ' '
          jr = 0
          call pvfil1(recl,len(aa),jr,a,cex,nchr,mxchr,ctbl,ja,aa,
     .    ltrace)
          call strncp(a,aa,1,1,min(ja,recl))
          call skipbl(a,recl,j)
          j0 = j
          call chrps2(a,'([',2,recl,j,k)
          if (k .lt. 1) goto 999
          aa = ' '
          call strcop(aa,a(j0+1),j-j0,' ',i1)
C       ... does vector variable aa already exist?
          call getsvv(aa,ival,0,1,1,val)
C       ... case not already declared:
          if (ival .eq. 0) then
            if (k .ne. 2) call fexit(-1,1,' Exit -1: rdfiln: attempt'
     .      //' to use undeclared vector, '//'line %i',nrd)
            j = j+1
            if (.not. a2bin(a,nelt,2,0,']',j,recl)) goto 999
            call addsvv(aa,nelt,ival)
            call parsvv(a,recl,ival,nelt,1,j)
C       ... case already declared:
          else
            if (k .ne. 1) then
              call awrit1('#rf (warning): attempt'
     .        //' to redeclare variable '//aa//'%a, line %i',
     .        ' ',80,i1mach(2),nrd)
              goto 91
              goto 10
            endif
C       ...   Find indices for vector subblock
            j = j+1
            j0 = j
            ltmp2 = a2bin(a,i1,2,0,':',j0,recl)
            if (ltmp2) j = j0
            if (.not. a2bin(a,nelt,2,0,')',j,recl)) goto 999
            if (.not. ltmp2) i1 = nelt
            if (i .eq. LFIND) then
              call strncp(aa,a,1,1,recl)
              call skipbl(a,recl,j)
              fname = ' '
              call strcop(fname,a(j+1),ctlen,' ',k)
              j = j+k
              if (.not. a2bin(a,val,4,0,' ',j,recl)) goto 999
              k = 0
              do  45  i = i1, nelt
                call getsvv(' ',ival,1,i,i,xx)
                if (xx .eq. val) then
                  k = i
                  goto 46
                endif
   45         continue
   46         continue
              call lodsyv(fname,1,dble(k),k)
              goto 91
            else
              call parsvv(a,recl,ival,nelt-i1+1,i1,j)
            endif
          endif
          if (noskip) goto 91
          goto 10
        elseif ((i .eq. LHAR .or. i .eq. LHAR0 .or. i .eq. LENV)
     .  .and. noskip) then
          nchr0 = nchr
          j = j+k
          k = 1
          if (i .eq. LHAR0) k = 0
          if (i .eq. LENV) k = k+10

C       ... make substitutions for eg % ... {}...
          ja = 1
          asbst = ' '
          jr = 0
          call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja,
     .    asbst,ltrace)
          call parchv(asbst,ja,mxchr,cex,'= ',k,nchr,ctbl,j)
          if (noskip) goto 85
          goto 10
        elseif (i .eq. LCHAR .and. noskip) then
          nchr0 = nchr
          j0 = k+j
          i1 = i

          ja = 1
          asbst = ' '
          jr = 0
          call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja,
     .    asbst,ltrace)

          call skipbl(asbst,recl,j0)
          if (j0 .ge. recl) goto 10
          call tokmat(asbst(j0+1:),ctbl,nchr,ctlen,' ',i1,k,.false.)
C ... Replace existing string definition
          if (i1 .ge. 0) then
            i1 = i1+1
C ... Create new string definition
          else
            i1 = nchr+1
            ctbl(i1,1) = ' '
            call strcop(ctbl(i1,1),asbst(j0+1:),ctlen,' ',k)
            nchr = i1
            if (nchr .gt. mxchr)
     .      call rx('rdfiln: too many char variables')
          endif
          j0 = j0+k
C ...       Loop until expression evaluates to true
   73     call skipbl(asbst,recl,j0)
          ltmp2 = .true.
          ltmp = noskip
          j = j0
          call rdfilx(asbst,recl,.true.,j0,nrd,ctbl,mxchr,nchr,ltmp)
c           print *, aa(1:j0)
C           Expression evaluated to false; skip over associated expr.
          if (.not. ltmp) then
C             If at end, give up
            if (j0 .ge. recl) goto 87
C             This ensures movement past end of expression
            call chrpos(asbst,' ',recl,j0)
            call skipbl(asbst,recl,j0)
C             (chrpos works except doesn't allow spaces in string)
C             call chrpos(asbst,' ',recl,j0)
C             If at end, give up
            if (j0 .ge. recl) goto 87
            ctmp = asbst(j0+1:j0+1)
            if (ctmp .ne. '"') ctmp = ' '
            if (ctmp .eq. '"') j0=j0+1
            call chrpos(asbst,ctmp,recl,j0)
            if (ctmp .eq. '"') j0=j0+1
c             print *, (asbst(m), m=1,j0), '|'
C             If at end, give up
            if (j0 .ge. recl) goto 87
            goto 73
          endif
          call skipbl(asbst,recl,j0)
C           Expression evaluates to true, but no string: set to blank
          if (j0 .eq. recl) then
            ctbl(i1,2) = ' '
            if (noskip) goto 85
            goto 10
          endif
c           print *, (asbst(m), m=1,j0), '|'
C           call strcop(ctbl(i1,2),asbst(j0+1),ctlen,' ',k)
          k = 1
          ctmp = asbst(j0+1:j0+1)
          if (ctmp .ne. '"') ctmp = ' '
          if (ctmp .eq. '"') j0=j0+1
          ja = j0
          call chrpos(asbst,ctmp,recl,ja)
          aa = ' '
          call pvfil1(ja,ja,j0,asbst,cex,nchr,mxchr,ctbl,k,aa,ltrace)
          ctbl(i1,2) = aa
          if (ctmp .eq. '"') ja=ja+1
          if (noskip) goto 85
          goto 10
        elseif (i .eq. LHOW .and. noskip) then
          j = j+k
          call skipbl(a,recl,j)
          i = 2
          if (j .ge. recl) goto 32
          call tokmat(a(j+1),sdir,4,5,' ',i,j,.false.)
          goto (31,32,32,34), i+1
          call fexit(-1,1,' Exit -1: rdfiln failed to parse show, '//
     .    'line %i',nrd)
C      ...  'show stop'
   31     continue
          lshow = ipr .gt. 110
          goto 10
C      ...  'show all' or 'show lines'
   32     continue
          lshow = .true.
          if (i .eq. 2) goto 10
C      ...  'show vars'
   34     continue
          call awrit1('#rf %i:',' ',80,i1mach(2),nrd)
          call shosyv(0,0,0,6)
          call numsvv(ival)
          if (ival .gt. 0) then
            print *, '----'
            call shosvv(1,ival,i1mach(2))
          endif
          if (nchr .gt. 0) print *, '---- character variables:'
          do  35  k = 1, nchr
            print '(i4,2x,a,2x,a)', k, ctbl(k,1)(1:20),ctbl(k,2)
ccccccccccccccc
c           call awrit1('%,4i  '//ctbl(k,1)(1:20)//ctbl(k,2)//'%a',
c     .          ' ',ctlen+25,i1mach(2),k)
ccccccccccccccc
   35     continue
          goto 10
        elseif (i .eq. LRACE) then
          j0 = j+k
          ltmp = noskip
          if (.not. ltmp) goto 10
C       ... No expression, toggle trace and exit
          call skipbl(a,recl,j0)
          if (j0 .ge. recl) then
            if (ltrace .le. 0) then
              ltrace = 1
            else
              ltrace = 0
            endif
          else
            j = j0
            call rdfilx(a,recl,.true.,j0,nrd,ctbl,mxchr,nchr,ltmp)
            ltrace = 0
            if (ltmp) ltrace = 1
C         ... If a valid expression, set to numerical value
            if (ltmp) then
              ltmp2 = a2bin(a,k,2,0,' ',j,recl)
              if (ltmp2) ltrace = k
            endif
          endif
          call awrit3(' rdfiln %i: trace %?#n>0#set to %i'//
     .    '#turned off#',' ',80,i1mach(2),nrd,ltrace,ltrace)
          goto 10
        elseif (i .eq. LCHO .and. noskip) then
          iecho = 1
          j0 = j+k
        elseif ((i .eq. LXIT .or. i .eq. LTOP) .and. noskip) then
          j0 = j+k
          call skipbl(a,recl,j0)
          if (j0 .ge. recl) then
            if (i .eq. LXIT) goto 99
            call awrit1('#rf %i: stop encountered',
     .      ' ',80,i1mach(2),nrd)
            call cexit(-1,1)
          endif
          ltmp = a2bin(a,ltmp2,0,0,' ',j0,recl)
          j0 = j0-1
          iecho = 4
          if (ltmp .and. ltmp2) iecho = 3
        elseif ((i .eq. LNCLUD .or. i .eq. LNCLUO) .and. noskip) then
          filstk = filstk+1
          nrdstk(filstk) = nrd
          unit0(filstk) = lunit
          lunit = 99-filstk
          nrd = 0
          j = j+k
          call skipbl(a,recl,j)
          ja = 1
          fname = ' '
          call pvfil1(recl,len(fname),j,a,cex,nchr,mxchr,ctbl,ja,
     .    fname,ltrace)
          call skpblb(fname,len(fname),k)
          if (ltrace .gt. 0)
     .    call awrit1('%x#rf %i: '''//dir(i+1)//'%a'' opening'//
     .    ' file '''//fname//'%a''',aa,120,i1mach(2),nrdstk(filstk))
C       ... Should we use fopna?
          if (i .eq. LNCLUD) then
            open(lunit,FILE=fname(1:k+1),STATUS='UNKNOWN',ERR=98)
          elseif (i .eq. LNCLUO) then
            open(lunit,FILE=fname(1:k+1),STATUS='OLD',ERR=98)
          endif
C       ... Read new line with different logical unit
          goto 10
C       ... Error if cannot open file
   98     call rx('rdfiln: error opening include file "'//
     .    fname(1:k+1)//'"')
        elseif (i .eq. LDEF .and. noskip) then
C       ... For each name, swap with top and remove top variable
          j0 = k+j

          ja = 1
          asbst = ' '
          jr = 0
          call pvfil1(recl,len(asbst),jr,a,cex,nchr,mxchr,ctbl,ja,
     .    asbst,ltrace)

          call skipbl(asbst,recl,j0)
C           set ltmp = .false. unless udef -f , then ltmp = .true.
          ltmp = .false.
          if (j0 .lt. recl-3) then
            if (asbst(j0+1:j0+3) .eq. '-f ') then
              j0 = j0+2
              ltmp = .true.
            endif
          endif
   43     continue
          call skipbl(asbst,recl,j0)
          if (j0 .ge. recl) goto 10
          fname = ' '
          k = 0
          call strcop(fname,asbst(j0+1:),ctlen,' ',k)
          j0 = j0+k
C       ... First try to udef a char variable
          call tokmat(fname,ctbl,nchr,ctlen,' ',i1,k,.false.)
          i1 = i1+1
          if (i1 .gt. 0) then
            do  44  k = i1+1, nchr
              ctbl(k-1,1) = ctbl(k,1)
   44       ctbl(k-1,2) = ctbl(k,2)
            nchr = nchr-1
            goto 43
          endif
C       ... Otherwise, try to udef a scalar variable
          call getsyv(fname,val,k)
          if (k .eq. 0 .and. ltmp) goto 43
          if (k .eq. 0) call fexit(-1,1,' Exit -1 rdfiln line %i:'//
     .    '  no variable "'//fname//'%a" to undefine',nrd)
          call numsyv(nsyv)
          call togsyv(k,nsyv)
          nsyv = nsyv-1
          call clrsyv(nsyv)
          goto 43

        elseif (i .eq. LAVE .and. noskip) then
C       ... If keeping all anyway, forget this command
          if (nkeepv .lt. 0) goto 10
          j0 = k+j
C       ... if no arguments, save everything up to this point
          call skipbl(a,recl,j0)
          if (j0 .ge. recl) then
            call numsyv(nkeepv)
            goto 10
          endif
C       ... For each name, swap with name at nkeepv
   42     continue
          call skipbl(a,recl,j0)
          if (j0 .ge. recl) goto 10
          fname = ' '
          k = 0
          call strcop(fname,a(j0+1),ctlen,' ',k)
          j0 = j0+k
          call getsyv(fname,val,k)
          if (k .eq. 0) call fexit(-1,1,' Exit -1 rdfiln line %i:'//
     .    '  no variable "'//fname//'%a" to save',nrd)
          if (k .le. nkeepv) goto 42
          nkeepv = nkeepv+1
          call togsyv(k,nkeepv)
          goto 42
        elseif (i .eq. LACRO .and. noskip) then
          call strncp(aa,a,1,1,recl)
          j = j+k
          call macset(aa(j+1:),k)
          if (k .lt. 0)
     .    call rxs('rdfiln failed to parse macro',aa(j+1:))
          goto 10
        endif
      endif

C   --- Copy this line into recrd, substituting variable names ---
      if (a(1) .ne. cc1 .and. noskip) then
        call strncp(recrd,a,1,1,recl)
        if (lex) then
          jr = j0
          ja = 1
          call pvfil1(recl,recl,jr,recrd,cex,nchr,mxchr,ctbl,ja,a,
     .    ltrace)
          do  47  i = ja, recl
   47     a(i) = ' '
          if (iecho .ne. 0) then
            if (mod(iecho,2) .ne. 0) then
              j = 1
              call skipbl(a,recl,j)
              call skpblb(a,recl,k)
              write(i1mach(4),346) nrd, (a(jr+1), jr=j,k)
  346         format('#rf', i5, ': ', 256a1)
C               call cwrite(a,j,k,1)
            endif
            if (iecho/2 .eq. 1) then
              call awrit1('#rf %i: stop encountered',
     .        ' ',80,i1mach(2),nrd)
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
      if (filstk .gt. 0) then
        if (ltrace .gt. 0) call awrit1('#rf closing include file: '
     .  //'%i lines',' ',80,i1mach(2),nrd)
        close(lunit)
        do  54  j = 0, nlev
   54   nlin(j) = nlin(j)-nrd
        lunit = unit0(filstk)
        nrd = nrdstk(filstk)
        filstk = filstk-1
C ...   And on our merry way, with the old unit
        goto 10
      endif
      if (nkeepv .ne. -1) call clrsyv(nkeepv)
      nr = -nr
      if (ltrace .gt. 0 .or. ipr .le. 100) return
      call awrit1(cc1//'rdfiln generated %i records.',
     .fname,len(fname),i1mach(2),-nr)
      return

C --- Trace entry points ---
C ... if-else-endif constructs
   81 if (ltrace .eq. 0) goto 10
      call awrit5('%a''%?#n>1# (nesting=%i)#%j#, read following lines: '
     .//'%?#n#yes#no#%?#n#%2b(blocked by level %i)',
     .aa,80,-i1mach(2),niff,niff,isw(liff(niff)),
     .isw(.not.liff(max(niff,1)-1).and.niff.gt.1),niff-1)
      goto 10
C ... while/repeat-end constructs
   82 if (ltrace .eq. 0) goto 10
      call awrit5('%a''%?#n>1# (nesting=%i)#%j#, read following lines: '
     .//'%?#n#yes#no#%?#n#%2b(blocked by level %i)',
     .aa,80,-i1mach(2),nlev,nlev,isw(loop0(nlev)),
     .isw(.not.loop0(max(nlev,1)-1).and.nlev.gt.1),nlev-1)
      goto 10
C ... const,var directives
   83 if (ltrace .gt. 1) then
        call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
        call numsyv(j)
        do  84  i = nsyv+1, j
          call watsyv(fname,val,i)
   84   call awrit1('%a '//fname//'%a%?#n<0#,##',aa,80,0,i-j)
        if (nsyv .eq. j) call awrit0('%a ... no new vars',aa,80,0)
        call awrit0('%a',aa,80,-i1mach(2))
      endif
      goto 10
C ... char directives
   85 if (ltrace .gt. 1) then
        call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
        do  86  i = nchr0+1, nchr
   86   call awrit1('%a '//ctbl(i,1)//'%a%?#n<0#,##',aa,80,0,i-nchr)
        if (nchr0 .eq. nchr) call awrit0('%a ... no new vars',aa,80,0)
        call awrit0('%a',aa,80,-i1mach(2))
      endif
      goto 10
C ... cchar read failed
   87 if (ltrace .gt. 1) then
        call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', var '//ctbl(nchr,1)//
     .  '%a not changed',aa,80,-i1mach(2),nrd)
      endif
      nchr = nchr0
      goto 10
C ... vec directives
   91 if (ltrace .gt. 1) then
        call awrit1('%x#rf %i: '''//dir(i+1)//'%a'', added',aa,80,0,nrd)
        call numsvv(j)
        do  92  i = nsyv+1, j
          call watsvv(fname,i)
   92   call awrit1('%a '//fname//'%a%?#n<0#,##',aa,80,0,i-j)
        if (nsyv .eq. j) call awrit0('%a ... no new vars',aa,80,0)
        call awrit0('%a',aa,80,-i1mach(2))
      endif
      goto 10

C --- General error exit ---
  999 call fexit(-1,1,' Exit -1: rdfiln: parse failed, line %i',nrd)
      end
      subroutine rdfile(unit,cch,recrd,mxrecs,a,recl,nr)
C- Reads entire disk file into recrd, calling rdfiln.
      implicit none
      integer unit,mxrecs,recl,nr
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=40,mxlev=6,lstsiz=200,ctlen=120)
      character*(*) cch, recrd(0:*)*1, a(*)*1
      character ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),
     .ilist(mxlev),nlist(0:mxlev)
      character vnam(mxlev)*16

      nr = 0
   10 call rdfiln(unit,cch,mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr,a,recrd(nr*recl),recl,nr)
      if (nr .le. 0) then
        nr = -nr
        return
      endif
      if (nr .eq. mxrecs) return
      goto 10

      end
      subroutine rdfilx(a,recl,lif,j,nrd,ctbl,mxchr,nchr,liff0)
C- Determines whether an 'if' or 'ifdef' expression is true
C  Inputs
C    a,recl
C    lif    T if 'if', F if to parse for 'ifdef'
C    j      starting position in a
C    nrd    line number (no longer used)
C    ctbl,mxchr,nchr character table, leading dimension and number
C    liff0  F if in a false branch of a conditional-read, else T
C  Outputs
C    liff0  Set to F if no T expression found; see Remarks.
C    j      last character searched
C  Remarks
Cr  Syntax for 'ifdef': ifdef expr [ op expr ... ]
Cr  Here op is the logical logical OR '|' or the logical AND '&'.
Cr  The entire collection expr op expr op expr ... can be partitioned
Cr  into the following:
Cr    superexpr | superexpr | superexpr ...
Cr  where
Cr    superexpr has the form expr [& expr ...]
Cr
Cr  A superexpr is true if and only if each of its separate expressions
Cr  is true.  The entire
Cr    superexpr | superexpr | superexpr ...
Cr  is true if any one of the the superexpr is true, whether or not
Cr  even they are valid expressions.
Cr  A single 'expr' is true if any of the following are true:
Cr    1.  it is a numerical expression that evaluates to true
Cr    2.  It is an invalid numerical expression, but expr is
Cr        the name of a character variable
Cr    3.  It is of the form    var=='string'
Cr        where var is the name of a character variable whose contents
Cr        are 'string'
C ----------------------------------------------------------------
      implicit none
      logical lif,liff0
      integer recl,j,nrd,nchr,mxchr,ctlen
      parameter (ctlen=120)
      character*(*) a,ctbl(mxchr,2)*(ctlen),aa*(ctlen)
C Local variables
      logical a2bin,ltmp,liffa
      integer j0,i0,i1,k

      liffa = .true.

   10 continue
      call skipbl(a,recl,j)
      j0 = j
C ... Attempt to read arithmetic expression
      ltmp = a2bin(a,liff0,0,0,' ',j,recl)
      if (ltmp) j = j-1
C ... If invalid expression, liff0 is assumed false
      liff0 =  liff0 .and. ltmp
C ... Case invalid-numerical-expr; check for character expression
      if (.not. ltmp) then
        j = j0
        call chrps2(a,' =',2,recl,j,k)
        aa = a(j0+1:j)
        call tokmat(aa,ctbl,nchr,ctlen,' ',i0,i1,.false.)
C ...   Check for a character expression
        if (i0 .ge. 0) then
C     ... if a ' ' was encountered before an '='
          if (k .eq. 1) then
            liff0 = .true.
C     ... if an '=' was encountered before a ' '
          else if (k .eq. 2) then
            if (a(j+1:j+3) .eq. '==''') then
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
C ... Now liff0 is result of current expression.
C     Nothing more to do for 'if' statement
      if (lif) return
C ... logical liff0 with a prior 'and'
      ltmp = liff0
      liff0 = liff0 .and. liffa
C ... Case 'and' following this expression
      call skipbl(a,recl,j)
      j = j+1
      if (j .ge. recl) return
      if (a(j-1:j+1) .eq. ' & ') then
        liffa = liffa .and. liff0
        goto 10
C ... Case 'or' following this expression
      elseif (a(j-1:j+1) .eq. ' | ') then
C   ... The end of possible expr & expr & ...
        if (liff0) return
        liffa = .true.
        goto 10
      endif

C ... Put j at last nw char, equivalent to pos of a2bin above
      call skpblb(a,j-1,j)
      j = j+1
      end
      subroutine pvfil1(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,ja,a,t)
C- Kernel of rdfiln that parses string with substitutions
C ----------------------------------------------------------------
Ci Inputs
Ci   reclr   length of recrd
Ci   recla   length of a
Ci   ja,jr   starting positions of a,recrd
Ci   recrd   string that is to be subsitituted
Ci   ctbl,nchr,mxchr table of character variables, number and dim.
Ci   t       trace, print out intermediate substitutions if t>=3
Co Outputs
Co   ja,jr   as final positions of a,recrd
Co   a       substituted string assembled into a.
Cr Remarks
Cr  Note recrd(0..), a(1..)!
C ----------------------------------------------------------------
      implicit none
C Passed Parameters
      integer reclr,recla,nchr,ja,jr,t,mxchr
      integer ctlen
      parameter (ctlen=120)
      character*1 recrd(0:reclr-1),a(recla),cex*2,ctbl(mxchr,2)*(ctlen)
C Local variables
      logical pvfil2
      character*512 aa
      integer ia,ja0,jaa

      do  5  ia = ja, recla
    5 a(ia) = ' '
      ja0 = ja
      if (.not. pvfil2(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,ja,a))
     .return
   10 continue
      aa = ' '
      do  11  ia = ja0, ja
   11 aa(ia:ia) = a(ia)
      if (t .gt. 2) then
        call skpblb(aa,ja,ja)
        print 333, '#rf subst: "', aa(1:ja+1),'"'
  333   format(a,a,a)
      endif
      jaa = ja0-1
      ja = ja0

      if (pvfil2(recla,recla,jaa,aa,cex,nchr,mxchr,ctbl,ja,a)) goto 10

      end
      logical function pvfil2(reclr,recla,jr,recrd,cex,nchr,mxchr,ctbl,
     .ja,a)
C- Kernel of rdfiln that parses string with substitutions
C  and returns T if substitution was inside a nested block.
C ----------------------------------------------------------------
Ci Inputs
Ci   reclr   length of recrd
Ci   recla   length of a
Ci   ja,jr   starting positions of a,recrd
Ci   recrd   string that is to be subsitituted
Ci   cex     characters demarcating string substitution, eg '{}'
Ci   ctbl,nchr,mxchr table of character variables, number and dim.
Co Outputs
Co   ja,jr   as final positions of a,recrd
Co   a       substituted string assembled into a.
Cr Remarks
Cr  Note recrd(0..), a(1..)!
C ----------------------------------------------------------------
      implicit none
      integer mxchr,reclr,recla,nchr,ja,jr,ctlen
      parameter (ctlen=120)
      character*1 recrd(0:reclr-1),a(recla),cex*2,ctbl(mxchr,2)*(ctlen),
     .aa*512,aa2*(ctlen),aa3*(ctlen),aastr*(ctlen),aasub*(ctlen)
C Local variables
      logical a2bin,a2bina,lqchr
      character*1 qchr,match
c #if CRAY
c      parameter (qchr='\')
c #else
      parameter (qchr='\\')
c #endif
      integer i,j,k,m,ival,kkk,i1,i2,ip1,ip2,
     .ip,ii(2),ix(2),istr1,istr2,a2vec,iaa21,iaa22,iaa2,iaa3
      equivalence (i1,ii(1)), (i2,ii(2))
      double precision res

      pvfil2 = .false.
      if (cex .eq. ' ') return
   10 continue
      if (jr .ge. reclr) return
C   --- Parse for next expression; quote \cex literally ---
      k = -1
      k = k+1
      call chrpos(recrd(jr),cex,reclr-jr,k)
C       lqchr is T if quote character literally
      lqchr = k.ge.1 .and. recrd(jr+k-1) .eq. qchr
C   --- Copy string up to substitution ---
      k = min(k,recla-ja,reclr-jr-1)
      if (k .lt. 0) return
      do  15  i = 0, k
   15 a(ja+i) = recrd(jr+i)
      jr = jr+k+1
      ja = ja+k+1
      if (lqchr) then
        ja = ja-1
        a(ja-1) = a(ja)
        goto 10
      endif
C   ... Exit if no substitution
      if (jr .ge. reclr .or. ja .gt. recla) return
C   --- Handle 'nesting' of cex by treating outer cex as a literal ---
      k = 0
      call chrps2(recrd(jr),cex,2,reclr-jr-1,k,kkk)
C   ... can't find corresponding close of cex
      if (kkk .eq. 0) then
        call pvfil3('missing "',cex(2:2),'" in line:',recrd,1,reclr)
C   ... cex is nested.  This pass, quote outer cex literally
      elseif (kkk .eq. 1) then
        do  16  i = 0, k-1
   16   a(ja+i) = recrd(jr+i)
        jr = jr+k
        ja = ja+k
C         print *, 'testing',(a(kkk), kkk=1,ja)
        pvfil2 = .true.
        goto 10
      endif

C   --- Parse and replace expression ---
C       a(ja) is 1st char after \cex
C   ... Case expr?string1:string2
      if (recrd(jr) .eq. '?') then
        match = recrd(jr+1)
C         Position of cex(2:2)
        k = k+jr
C         Start of expression, save for error message
        m = jr
C         Parse the expression
        jr = jr+2
        if (.not. a2bin(recrd,i,2,0,match,jr,reclr))
     .  call pvfil3('could ','not ','parse:',recrd,m,reclr)
C         First char after expr
        kkk = jr
C         Position of matching char
        call chrpos(recrd,match,reclr,jr)
        if (jr .ge. reclr) call pvfil3(
     .  'missing matching "',match,'" in line: ',recrd,m,reclr)
C         First string
        if (i .ne. 0) then
          m = kkk
C         Second string
        else
          m = jr+1
          match = cex(2:2)
        endif
C         Copy one string
        call strcop(a(ja-1),recrd(m),reclr-m,match,i)
C         Update ja,jr before resuming parsing
        ja = ja+i-2
        a(ja) = ' '
        jr = k+1
        goto 10
      endif

C   ... Start of expression should be a variable name; put into aa
      aa = ' '
      call strcop(aa,recrd(jr),reclr-jr,cex(2:2),i)
      match = aa(i:i)
      aa(i:i) = ' '
      call wordg(aa,10,'A-Za-z_0-9',1,i1,i2)
C       call skipbl(aa,ctlen,i2)
      if (i2+1 .lt. i) then
        match = aa(i2+1:i2+1)
        aa(i2+1:i2+1) = ' '
      endif

C   ... Case string aa matched an entry in ctbl
C       i=index to ctbl if start of aa holds member of ctbl
      call tokmat(aa,ctbl,nchr,ctlen,' ',i,m,.false.)
      if (i .ge. 0) then
        aa2 = ctbl(i+1,2)
C         Simple character variable substution
        if (match .eq. cex(2:2)) then
          i1 = 1
          i2 = ctlen
          goto 20
C         If not subexpression, subst does not involve a char variable
        elseif ( match .ne. '(') then
          i = -1
          goto 20
        endif
        ip = i2+1
        call skipbl(aa,ctlen,ip)
        ip = ip+1
C         Case aa = ctbl(i+1,1)(/s1/s2/[,n1,n2])
        if (aa(ip:ip+2) .eq. ':e)') then
          call strip(aa2,i1,i2)
          ja = ja-2
          call bin2a(' ',0,0,i2,2,0,recla,a,ja)
          ja = ja+1
          jr = jr+ip+3
          goto 10
        elseif (aa(ip:ip) .eq. '/') then
          istr1 = ip+1
          istr2 = istr1
          call cpstr(aa,ctlen,0,'/',istr2,ip1,aastr)
          if (istr2 .le. istr1 .or. istr2 .ge. ctlen) goto 999
          istr2 = istr2+1
          call cpstr(aa,ctlen,0,'/',istr2,ip2,aasub)
          if (istr2 .le. istr1 .or. istr2 .ge. ctlen) goto 999
          m = istr2+1
          if (aa(m:m) .eq. ',') then
            k = a2vec(aa,len(aa),m,2,',)',2,2,2,ix,ii)
            if (i1 .lt. 1) call pvfil3(
     .      'illegal',' parameter',' in',recrd,1,reclr)
            if (k .ne. 2) goto 999
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
              if (aa2(iaa2:iaa2+ip1-1) .eq. aastr(1:ip1)) then
                if (j .lt. i1) then
                  iaa21 = iaa2+1
                  goto 25
                else
                  aa3(iaa3:iaa3+ip2-1) = aasub(1:ip2)
                  iaa21 = iaa2+ip1
                  iaa3 = iaa3+ip2-1
                  if (j .eq. i2) aa3(iaa3+1:) = aa2(iaa21:)
                  goto 25
                endif
              endif
   24       continue
            call pvfil3('could not find substr "',aastr(1:ip1),'" in',
     .      recrd,1,reclr)
   25     continue
          aa2 = aa3
          i1 = 1
          i2 = ctlen
          call skpblb(aa2,ctlen,i2)
          i2 = i2+1
C         Case aa = ctbl(i+1,1)('char-list' [,n])
        elseif (aa(ip:ip) .eq. '''' .or. (aa(ip:ip) .eq. '"')) then
          i1 = ip+1
          i2 = ip
C           call cpstr(aa,ctlen,1,aa(ip:ip),ip,ip1,aasub)
          call eostr(aa,ctlen,21,aa(ip:ip),i2)
          if (i2 .ge. ctlen) goto 999
          i2 = i2-1
          kkk = i2+2
          if (aa(kkk:kkk) .eq. ')') then
            k = 1
          elseif (aa(kkk:kkk) .eq. ',') then
            if (.not. a2bin(aa,k,2,0,')',kkk,len(aa))) goto 999
          else
            goto 999
          endif
          call wordg(ctbl(i+1,2),100,aa(i1:i2),k,i1,i2)
          ja = ja-2
          call bin2a(' ',0,0,i2+1,2,0,recla,a,ja)
          ja = ja+1
          jr = jr+kkk+1
          goto 10
C         Assume aa = ctbl(i+1,1)(n1,n2)
        else
          ip1 = i2+1
          ip = ip1
          k = a2vec(aa,len(aa),ip,2,',)',2,2,2,ix,ii)
          if (k .ne. 2) goto 999
C           Add to m so that jr is incremented past (..)
          m = m + ip-ip1+1
        endif
      endif
C   ... End of test for character variables.  If string substitution,
C       i>0 and aa2(i1:i2) = (possibly modified) substring of ctbl
C       and m = number of characters to advance in recrd
   20 continue

C   ... Substitute a character variable
      if (i .ge. 0) then
        jr = jr+m
        m = i2
        call skpblb(aa2,ctlen,i2)
        m = min(m,i2+1)
        a(ja-1) = ' '
         print *, 'before 33',(a(kkk), kkk=1,ja-2+0),"###"
        do  33  k = i1, m
   33   a(ja+k-i1-1) = aa2(k:k)
         print *, 'after 33',(a(kkk), kkk=1,ja+m-i1-1),"###"
        ja = ja+m-i1

C   ... Numerical expression
      elseif (a2bina(recrd,res,4,0,cex(2:2),jr,reclr)) then
        ja = ja-2
        if (dabs(res) .gt. 1) then
          call bin2a('g9',0,0,res,4,0,recla,a,ja)
        else
          call bin2a('g9:10',0,0,res,4,0,recla,a,ja)
        endif
        ja = ja+1

      else
C   ... if a name for a vector, replace with ascii rep of entire vector
        call sizsvv(aa,ival,0,k)
        if (ival .gt. 0) then
          ja = ja-2
          do  35  m = 1, k
            call getsvv(' ',ival,1,m,m,res)
            if (dabs(res) .gt. 1) then
              call bin2a('g9',0,0,res,4,0,recla,a,ja)
            else
              call bin2a('g9:10',0,0,res,4,0,recla,a,ja)
            endif
            ja = ja+1
            a(ja) = ' '
   35     continue
          call chrpos(recrd,cex(2:2),reclr,jr)
          jr = jr+1
        else
          print *, 'rdfile: bad expression in line'
          print *, (recrd(j), j=0, jr-1), '  ...  ',
     .    (recrd(j), j=jr, reclr-1)
          call cexit(-1,1)
        endif
      endif
      goto 10

  999 continue
      call pvfil3('could',' not',' parse',recrd,1,reclr)

      end
      subroutine pvfil3(s1,s2,s3,recrd,jr,reclr)
C- Error exit
      implicit none
      integer reclr,i,ii,jr
      character*1 recrd(0:1)
      character*(*) s1,s2,s3
      character aa*100
      aa = 'rdfiln '// s1 // s2 // s3 // ' ... '
      call skpblb(aa,len(aa),ii)
      ii = ii+3
      do  18  i = ii, min(reclr-jr,len(aa))
   18 aa(i:i) = recrd(i+jr-ii)
      call setpr(20)
      call fexit(-1,009,aa,0d0)

      end
      logical function rdstrn(unit,a,len,lopt)
C- Read one line from logical unit
      integer unit,len
      logical lopt
      character*1 a(len)

      rdstrn = .true.
      read(unit,10,end=20,err=20) a
      if (lopt) print 10, a
      return
   10 format(2048a1)
   20 rdstrn = .false.
      return
      end

