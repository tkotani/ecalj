      subroutine query(instr,cast,var)
C- interactive flow control
C ----------------------------------------------------------------
Ci Inputs
Ci   strng: prompt string
Ci   cast:  <0, if nothing to change, otherwise
Ci          cast is  0=,logical, 2=int, 3=real 4=double
Co Outputs
Co   var:   query will change if requested
Cr Remarks
Cr   At the prompt, user enters either nothing, or one of
Cr     'Snnn', where nnn is number (or T or F for logical variable);
Cr     'Vnnn', where nnn is the new verbosity;
Cr     'W' to toggle printing of work array;
Cr     'I' to turn off interactive prompt;
Cr     'A' to abort program execution
Cu Updates
Cu   07 Jul 04 suppress interactive mode in MPI
Cu   19 May 04 Added getqu: returns value of interative mode
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) instr
      integer cast,var
      logical lset
C Local parameters
      integer lprmpt
      parameter (lprmpt=120)
      character*(lprmpt) rdstr, outstr
      logical lquery,lsequ,a2bin,lact
      integer i,ivbs,iprint,j,i1mach,l,fxst,fopna,ifi,awrite
      double precision xx
      save lquery
      character*(1) rdstr2(2)
      equivalence (rdstr,rdstr2)
      external lsequ,a2bin,iprint

#if MPI | MPIK
      return
#endif 
      ifi = 29
      lact = fxst('iact') .ne. 0
      goto 12
C --- Start of interactive loop (read from iact first pass only) ---
   10 continue
      if (lact) call dfclos(ifi)
      lact = .false.
   12 continue
      if (lact) then
        ifi = 29
        ifi = fopna('iact',ifi,1)
        rewind ifi
      endif
      if (.not. lquery .and. .not. lact) return

      outstr = ' QUERY:'
      outstr(9:lprmpt) = instr
      call skpblb(outstr,lprmpt,l)
      if (cast .ge. 0) then
        l = awrite('%a (def=',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
        call bin2a(' ',0,0,var,cast,0,lprmpt,outstr,l)
        call awrit0('%a)',outstr,lprmpt,0)
      endif
      l = awrite('%a?',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
      call cwrite(outstr,0,l,0)
  334 format(a60)
C --- Read string, exit if empty ---
      rdstr = ' '
      if (lact) then
        print *, ' '
        read(ifi,334,end=10,err=10) rdstr
      else
        read(*,334) rdstr
      endif
      j = 0
      call skipbl(rdstr,lprmpt,j)
      if (j .ge. lprmpt) return
      j = j+1
      i = 0
C --- Handle query ---
   14 continue
      if (rdstr2(j) .eq. '?') then
        outstr = ' (A)bort  (I)active  (V)erb  (C)pu  (T)iming  (W)ork'
        if (cast .ge. 0) call awrit0('%a  (S)et value',outstr,lprmpt,0)
        print '(a)',  outstr
        goto 10
      elseif (rdstr2(j) .eq. 'S' .or. rdstr2(j) .eq. 's') then
        if (lact .and. rdstr2(j+1) .eq. '"') then
          l = 0
          call chrpos(rdstr2(j+2),'"',len(rdstr)-j-1,l)
          if (lsequ(rdstr2(j+2),instr,l,'"',l)) then
            j =j+1+l
          else
            lact = .false.
            call fclose(ifi)
            goto 10
          endif
        endif
        if (cast .lt. 0) print *, 'no variable to set--continuing ...'
        if (cast .lt. 0) goto 10
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1))
     .  print *, 'conversion error'
        goto 10
      elseif (rdstr2(j) .eq. 'V' .or. rdstr2(j) .eq. 'v') then
        if (a2bin(rdstr2(j+1),ivbs,2,0,' ',i,-1))  then
          call setpr(ivbs)
        else
          print *, 'conversion error'
        endif
        goto 10
      elseif (rdstr2(j) .eq. 'W' .or. rdstr2(j) .eq. 'w') then
        call wkprnt(2)
        call wkchk('Called from query')
        call wkinfo
        goto 10
      elseif (rdstr2(j) .eq. 'I' .or. rdstr2(j) .eq. 'i') then
        lquery = .not. lquery
      elseif (rdstr2(j) .eq. 'A' .or. rdstr2(j) .eq. 'a') then
        call rx0(outstr)
      elseif (rdstr2(j) .eq. 'Q' .or. rdstr2(j) .eq. 'q') then
        if (lact) call dfclos(ifi)
        call rx0(outstr)
      elseif (rdstr2(j) .eq. 'T' .or. rdstr2(j) .eq. 't') then
        call tc('tog')
        goto 10
      elseif (rdstr2(j) .eq. 'C' .or. rdstr2(j) .eq. 'c') then
        call cpudel(i1mach(2),'called from query:   ',xx)
        goto 10
C ... string not empty, but no keyword recognized.
      elseif (.not. lact) then
        if (cast .lt. 0) return
        j = 0
        outstr = ' '
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1)) then
          call awrit0(' ? input "'//rdstr//'%a" not recognized',
     .    outstr,lprmpt,i1mach(2))
          rdstr2(1) = '?'
          j = 1
          goto 14
        endif
        goto 10
      endif

      if (lact) goto 10
      return

      entry initqu(lset)
      lquery = lset
      return

      entry getqu(lset)
      lset = lquery
#if MPI | MPIK
      lset = .false.
      return
#endif 

      end
C      subroutine fmain
C      double precision x
C      logical l
C
C      call initqu(.true.)
CC      call initqu(.false.)
C      call query('n?',-1,l)
C      call query('l=?',0,l)
C      print *, l
C      call query('i=?',2,i)
C      print *, i
C      x = -2
C      call query('d=?',4,x)
C      print *, x
C      end

