      integer function fgtcat(ifi,cat,descr,swt,lrewnd,mxrec,recln,
     .instr)
C- Read a category from a file
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit
Ci   cat   :Name of category to seek
Ci   descr :string describing the category's purpose (for input mode 0)
Ci   lrewnd:true, rewind file before searching; false, do not
Ci   mxrec :maximum number of records permitted to read
Cio Inputs/Outputs
Cio  swt  : a vector of switches governing how gtcat and rdtok behave.
Cio       : Its internal structure is documented in in stswt, below.
Cio       : User is advised to set these parameters via routine stswt.
Cio       : and query the state of using quswt.  Then there is no need
Cio       : know the internal structure of swt.
Cio
Cio       : As described in detail in stswt below, quantities that can
Cio       : be set and probed are the following:
Cio       : swt()  function (see stswt)
Cio       :  1     start-of category
Cio       :  2     end-of category
Cio       :  3     governs behavior when token or some data is missing,
Cio       :        and what to do if category exceeds mxrec records;
Cio       :        see stswt for description.
Cio       :  5     read and printout mode
Cio       :        0=> do not attempt to read category but show
Cio       :            what category would have been sought
Cio       :        1=> attempt to find category and read contents
Cio       :        2=> print category name
Cio       :        3=> combination of 1 and 2
Ci        :  6     file logical unit for data printout
Ci        :  7     not used by fgtcat
Co Outputs
Co   instr :category read into string `instr'
Co   fgtcat: -1 no category was found
Co         :  0 no attempt was made to read this category
Co         : >mxrec mxrec lines were read; end-of-category not found
Co         : >0 number of records read
Cr Remarks
Cr   this is a file version of subroutine gtcat
Cu Updates
Cu   30 Apr 02 ANSI version that conforms to ANSI standard
Cu   01 Nov 01 First created
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,swt(7),mxrec,recln
      logical lrewnd
      character*(*) cat,descr
#if ANSI | FIXEDLENGTH
      integer recl0
      parameter (recl0=128)
      character*(recl0) instr(mxrec)
#else
      character*(recln) instr(mxrec)
#endif
C ... Local parameters
      integer j2,nrec,ifipr,stderr,lgunit,quswt,lio,lopt,lmxrec
c     The following are specific to rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character ctbl(mxchr,2)*(ctlen),a*(recl),aa*(recl),
     .vnam(mxlev)*16,rdarg*7
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .nlist(0:mxlev)

#if ANSI | FIXEDLENGTH
      if (recln .ne. recl0) call
     .rxi('fgtcat: recompile with recl0=',recln)
#endif
      lio   = quswt('io',swt)
      lopt  = quswt('token,opt',swt)
      lmxrec = quswt('cat,mxrec',swt)
      ifipr = quswt('lunit',swt)
      stderr = lgunit(1)

C --- Show what category would have been sought ---
      if (lio .eq. 0) then
        aa = '%N category  ' // cat //
     .  '%a %?#(n==0)#(optional)#(required)#'
        call awrit1(aa,' ',72,ifipr,lopt)

        if (descr .ne. ' ') then
          aa = '   ' // descr
          call awrit0(aa,' ',-len(aa),ifipr)
        endif
        fgtcat = 0
        return

      elseif (lio .eq. 2) then
        goto 999

      elseif (lio .ne. 1 .and. lio .ne. 3) then
        call rxi('fgtcat : not implemented, mode',lio)

      endif

C --- Find category and read it in ---
      call isanrg(recln,1,recl,'fgtcat:','recln',.true.)
      if (mxrec .le. 0) call rxi('fgtcat: illegal mxrec',mxrec)
      call stswt('start',1,swt)
      if (lrewnd) rewind ifi
      nr = 0
      nrec = 0
      fgtcat = -1
      j2 = len(cat)
      rdarg = '#{}% ct'
   10 continue
      call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr,a,aa,recl,nr)

C ... End-of-file
      if (nr .le. 0) then
        call stswt('stop',nrec*recln,swt)
        if (nrec .gt. 0) goto 999
        if (lopt .eq. 1) then
          write(stderr,'(1x)')
          call rxs2('fgtcat: failed to find category `',cat,'''')
        endif
        return
      endif

C ... End-of-category
      if (nrec .gt. 0 .and. aa(1:1) .ne. ' ') then
        backspace ifi
C       swt(2) = nrec*recln
        call stswt('stop',nrec*recln,swt)
        goto 999
      endif

C ... Append current line to the category
      if (nrec .gt. 0) then
        nrec = nrec+1
        fgtcat = nrec
C       Maximum number of records exceeded
        if (nrec .gt. mxrec) then
          if (lmxrec .ne. 0) call rxi(
     .    'fgtcat : category larger than allowed maximum no lines='
     .    ,mxrec)
C         swt(2) = (nrec-1)*recln
          call stswt('stop',(nrec-1)*recln,swt)
          goto 999
        endif
        instr(nrec) = aa

C ... See whether category was found
      elseif (aa(1:j2) .eq. cat) then
        nrec = 1
        fgtcat = nrec
        instr(nrec) = aa
      endif

C ... Read a new line
      goto 10

C --- Cleanup: check for printout ---
  999 continue

C ... Print category name
      if (lio .ge. 2) then
        write(ifipr,'(a)') cat
      endif

      end

      subroutine stswt(instr,n,swt)
C- Sets flow parameters for fgcat and rdtok
C ----------------------------------------------------------------------
Ci Inputs
Ci   instr: string specifying what parameter in swt to set
Ci        : instr consists of 1 or more substrings separated by spaces.
Ci        : A substring consists of 1 or more words separated by commas.
Ci        : The first word refers to which parameter in swt to set;
Ci        : subsequent words are modifiers specifying how swt is set.
Ci
Ci        : Example:   instr = 'token,reqd,varlen  io,read'
Ci        : has 2 substrings, each setting one element in swt.
Ci        :
Ci        : The following table list the parameters and modifiers.
Ci        : Normally the caller need not know the internal structure of
Ci        : swt but can set states via stswt, and probe them using quswt
Ci
Ci        : NB: the same modifiers apply equally to categories
Ci        :     and tokens, except where noted.
Ci        :    parm  swt()  modifier  function, and assignment to sw()
Ci             ------------------------------------------------------
Ci        :    start   1     ---     set start-of-string index (where
Ci        :                          rdtok begins parsing). Assign index
Ci        :                          to contents of argument n.
Ci        :                          NB: fgtcat sets this parameter when
Ci        :                          a new category is found
Ci        :                          sets swt(1) = contents of n
Ci
Ci        :    stop    2     ---     set end-of-string index (where
Ci        :                          rdtok stop parsing).  Assign index
Ci        :                          to contents of argument n.
Ci        :                          NB: fgtcat sets this parameter when
Ci        :                          a new category is found
Ci        :                          sets swt(2) = contents of n
Ci        :
Ci        :    lunit   6     ---     file logical unit for output
Ci        :                          sets swt(6) = contents of n
Ci        :
Ci        :    io      5    read     parse token and contents
Ci        :                          sets swt(5)=1
Ci        :    io           write    show token and contents
Ci        :                          sets swt(5)=2
Ci        :    io           rw       both read and write
Ci        :                          sets swt(5)=3
Ci        :    io           help     'help' mode: show what token would
Ci        :                          be read and its properties
Ci        :                          sets swt(5)=0
Ci
Ci        :    token   3             has 3 independent parts.
Ci                                  *First part (bits 0,1)
Ci        :    token        opt      token is optional
Ci        :                          sets swt(3),bits(0+1) = 0
Ci        :    token        o        synonymous with 'opt'
Ci        :    token        reqd     token is required
Ci        :                          sets swt(3),bits(0+1) = 1
Ci        :    token        r        synonymous with 'reqd'
Ci        :    token        alt      an alternate token may be
Ci        :                          use in place of this one
Ci        :                          sets swt(3),bits(0+1) = 2
Ci        :    token        ignr     token is ignored
Ci        :                          sets swt(3),bits(0+1) = 3
Ci        :    token        i        synonymous with 'ignr'
Ci
Ci                                  *Second part (bit 2):
Ci        :    token        fixlen   Number of arguments following token
Ci        :                          is fixed.
Ci        :                          sets swt(3),bit(2) = 1
Ci        :    token        varlen   Number of arguments following token
Ci        :                          may be be fewer than
Ci        :                          the total number sought
Ci        :                          sets swt(3),bit(2) = 0
Ci        :    token        mxrec    Synonymous with `varlen'
Ci                                  *Third part (bit 4):
Ci        :    token        mult     tells rdtok to limit endpoint of
Ci        :                          search to smaller of stop (=swt(2))
Ci        :                          and second occurence of same token
Ci        :                          (a flag for multiple occurrences
Ci        :                           of this token in the category)
Ci        :                          sets swt(3),bit(4) = 1
Ci        :    token        nomult   turns off 'mult'
Ci        :                          sets swt(3),bit(4) = 0
Ci
Ci        :    ... The modifiers may be used by both fgtcat and rdtok
Ci        :        though some are not used by only one or the other
Ci        :        and following modifiers are interpreted differently.
Ci        :        Also, you can substitute the string `cat' for `token'
Ci        :    cat          mxrec    About if record length of
Ci        :                          category exceeds mxrec
Ci        :
Ci      n : integer parameter used for some of the specifications
Ci        : above; see description for `instr'
Co Outputs
Co   swt  : set of rdtok switches; some members may be set here
Co        : according to specifications in `instr'.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) instr
      integer n,swt(7)
C ... Local parameters
      integer j10,j1,j2,j3,i0,i2,i3,i4,k

      j1 = 1
C ... Parse for new word
   10 continue
      call nword(instr,1,j1,j2)
      if (j2 .lt. j1) return
      j10 = j1
      if (instr(j1:j1+1) .eq. 'io') then
   20   continue
        call nwordg(instr,1,', ',2,j1,j3)
        if (j3 .ge. j1) then
          if (instr(j1:j3) .eq. 'read') then
            swt(5) = 1
          elseif (instr(j1:j3) .eq. 'write') then
            swt(5) = 2
          elseif (instr(j1:j3) .eq. 'rw') then
            swt(5) = 3
          elseif (instr(j1:j3) .eq. 'help') then
            swt(5) = 0
          else
            goto 99
          endif
          if (j3 .lt. len(instr)) then
            if (instr(j3+1:j3+1) .eq. ',') then
              j1 = j3+1
              goto 20
            endif
          endif
        endif

      elseif (instr(j1:j1+4) .eq. 'token' .or.
     .instr(j1:j1+2) .eq. 'cat') then
        k = 3
C        if (instr(j1:j1+2) .eq. 'cat') k = 4
        i0 = mod(swt(k),4)
        i2 = mod(swt(k)/4,2)
        i3 = mod(swt(k)/8,2)
        i4 = mod(swt(k)/16,2)
   30   continue
        call nwordg(instr,1,', ',2,j1,j3)
        if (j3 .ge. j1) then
          if (instr(j1:j3) .eq. 'opt') then
            i0 = 0
          elseif (instr(j1:j3) .eq. 'o') then
            i0 = 0
          elseif (instr(j1:j3) .eq. 'reqd') then
            i0 = 1
          elseif (instr(j1:j3) .eq. 'r') then
            i0 = 1
          elseif (instr(j1:j3) .eq. 'alt') then
            i0 = 2
          elseif (instr(j1:j3) .eq. 'a') then
            i0 = 2
          elseif (instr(j1:j3) .eq. 'ignr') then
            i0 = 3
          elseif (instr(j1:j3) .eq. 'i') then
            i0 = 3
          elseif (instr(j1:j3) .eq. 'varlen') then
            i2 = 0
          elseif (instr(j1:j3) .eq. 'fixlen') then
            i2 = 1
          elseif (instr(j1:j3) .eq. 'mxrec') then
            i2 = 1
          elseif (instr(j1:j3) .eq. 'nomult') then
            i4 = 0
          elseif (instr(j1:j3) .eq. 'mult') then
            i4 = 1
          else
            goto 99
          endif
          if (j3 .lt. len(instr)) then
            if (instr(j3+1:j3+1) .eq. ',') then
              j1 = j3+1
              goto 30
            endif
          endif
        endif
        swt(k) = swt(k) - mod(swt(k),32) + i4*16 + i2*4 + i0

      elseif (instr(j1:j1+4) .eq. 'lunit') then
        swt(6) = n

      elseif (instr(j1:j1+3) .eq. 'init') then
        call iinit(swt,7)

      elseif (instr(j1:j1+4) .eq. 'start') then
        swt(1) = n

      elseif (instr(j1:j1+3) .eq. 'stop') then
        swt(2) = n

      else
        goto 99
      endif
      j1 = j2+1
      goto 10

   99 continue
      call rxs('stswt : bad mode: ',instr(j10:j2))
      end

      integer function quswt(instr,swt)
C- Reads flow parameters for fgcat and rttok
C ----------------------------------------------------------------------
Ci Inputs
Ci   instr: string specifying what to read
Ci        : instr consists of a word optionally followed by
Ci        : a modifier.  The two are separated by a comma.
Ci        : The first word refers to which parameter is referenced
Ci        : the modifier specifying what information is sought.
Ci
Ci        : Example:   instr = 'token,reqd'
Ci        : returns 0 if the token is required
Ci        : returns 1 if the token is not required
Ci        :
Ci        : See stswt above for a description of all parameters and
Ci        : their modifiers and corresondence between switches and
Ci        : structure of swt.  However, the purpose of this routine
Ci        : is to enable the user to set and probe states of swt
Ci        : without knowing how swt is internally organized.
Ci
Ci        :    parm  swt()  modifier  returns
Ci             ------------------------------------------------------
Ci        :    start   1     ---     start-of-string index where
Ci        :                          rdtok begins parsing.
Ci
Ci        :    stop    2     ---     end-of-string index where
Ci        :                          rdtok stop parsing.
Ci
Ci        :    io      5     ---     0 if help mode
Ci        :                          1 if 'read' mode
Ci        :                          2 if 'write' mode
Ci        :                          3 if 'read+write' mode
Ci        :    io           read     1 if io parameter is
Ci                                   set to 'read' or 'rw' mode,
Ci                                   0 otherwise
Ci        :    io           r        synonymous with 'read'
Ci        :    io           write    1 if io parameter is
Ci                                   set to 'write' or 'rw' mode,
Ci                                   0 otherwise
Ci        :    io           w        synonymous with 'write'
Ci        :    io           help     1 if io parameter is
Ci                                   set to 'help' mode
Ci
Ci        :    lunit   6     ---     file logical unit for output
Ci        :
Ci        :    token   3    opt      0 if token is optional
Ci        :                          1 if token is required
Ci        :                          2 alternate token may be subst.
Ci        :                          3 if token is ignored
Ci        :    token        reqd     1 if token is required, else 0
Ci        :    token        r        synonymous with 'reqd'
Ci        :    token        alt      1 if an alternate token may be
Ci        :                          use in place of this one
Ci        :    token        a        synonymous with 'alt'
Ci        :    token        ignr     1 if token is ignored
Ci        :    token        i        synonymous with 'ignr'
Ci        :    token        fixlen   1 if number of arguments following
Ci        :                          is fixed.
Ci        :                          0 if number of arguments following
Ci        :                          token may be be fewer than
Ci        :                          the total number sought
Ci        :    token        mult     1 'multiple token mode'
Ci
Ci        :    ... The modifiers may be used by both fgtcat and rdtok
Ci        :        though some are not used by only one or the other
Ci        :        and following modifiers are interpreted differently.
Ci        :        Also, you can substitute the string `cat' for `token'
Ci        :    cat          mxrec    maximum number of records allowed
Ci        :                          in category
Ci        :
Ci      n : integer parameter used for some of the specifications
Ci        : above; see description for `instr'
Co Outputs
Co   swt  : set of rdtok switches; some members may be set here
Co        : according to specifications in `instr'.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) instr
      integer swt(7)
C ... Local parameters
      integer j10,j1,j2,j3,i0,i2,i3,i4,k

      quswt = 0
      j1 = 1
      call nword(instr,1,j1,j2)
      if (j2 .lt. j1) return
      j10 = j1
      if (instr(j1:j1+1) .eq. 'io') then
        quswt = swt(5)
        call nwordg(instr,1,', ',2,j1,j3)
        if (j3 .ge. j1) then
          quswt = 0
          if (instr(j1:j3).eq.'read' .or. instr(j1:j3).eq.'r') then
            if (swt(5) .eq. 1 .or. swt(5) .eq. 3) quswt = 1
          elseif (instr(j1:j3).eq.'write' .or. instr(j1:j3).eq.'w') then
            if (swt(5) .eq. 2 .or. swt(5) .eq. 3) quswt = 1
          elseif (instr(j1:j3).eq.'help') then
            if (swt(5) .eq. 0) quswt = 1
          else
            goto 99
          endif
        endif

      elseif (instr(j1:j1+4) .eq. 'token' .or.
     .instr(j1:j1+2) .eq. 'cat') then
        k = 3
C       if (instr(j1:j1+2) .eq. 'cat') k = 4
        i0 = mod(swt(k),4)
        i2 = mod(swt(k)/4,2)
        i3 = mod(swt(k)/8,2)
        i4 = mod(swt(k)/16,2)
        call nwordg(instr,1,', ',2,j1,j3)
        if (j3 .ge. j1) then
          if (.false.) then
            continue
          elseif (instr(j1:j3).eq.'opt' .or. instr(j1:j3).eq.'o') then
            quswt = i0
          elseif (instr(j1:j3).eq.'reqd' .or. instr(j1:j3).eq.'r') then
            if (i0 .eq. 1) quswt = 1
          elseif (instr(j1:j3).eq.'alt' .or. instr(j1:j3).eq.'a') then
            if (i0 .eq. 2) quswt = 1
          elseif (instr(j1:j3).eq.'ignr' .or. instr(j1:j3).eq.'i') then
            if (i0 .eq. 3) quswt = 1
          elseif (instr(j1:j3) .eq. 'fixlen') then
            if (i2 .eq. 1) quswt = 1
          elseif (instr(j1:j3) .eq. 'mxrec') then
            if (i2 .eq. 1) quswt = 1
          elseif (instr(j1:j3) .eq. 'mult') then
            if (i4 .eq. 1) quswt = 1
          else
            goto 99
          endif
        endif

      elseif (instr(j1:j1+4) .eq. 'lunit') then
        quswt = swt(6)

      elseif (instr(j1:j1+4) .eq. 'start') then
        quswt = swt(1)

      elseif (instr(j1:j1+3) .eq. 'stop') then
        quswt = swt(2)

      else
        goto 99
      endif
      return

   99 continue
      call rxs('quswt : bad mode: ',instr(j10:j2))
      end

