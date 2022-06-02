integer function parg(tok,cast,strn,ip,lstr,sep,itrm,narg,it,res)
  !- Returns vector of binary values from a string
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   tok:  token marking input
  !i  cast:  0=logical, 2=int, 3=real, 4=double
  !i  strn(ip:lstr):string to parse, from (ip=0 for first char)
  !i  lstr:  length of strn
  !i  sep:   string of characters, each of which separates arguments
  !i  itrm:  characters sep(itrm:*) signal the last argument
  !i  narg:  number of values to parse.
  ! o Inputs/Outputs
  !o   ip:   on input, position in strn where to start parsing
  !o         on ouput, position in strn on exit.
  !o Outputs
  !o   res:  Vector of numbers that were converted
  !o   it:   Vector of indices, one for each entry in res, labeling
  !o         which char in 'sep' terminated the expr. for that entry.
  !o parg:   0 if token is not matched in strn.
  !o         n if token match and converted sans error narg numbers
  !o           (for narg=0, returns 1 if token matched)
  !o        -n if error on conversion of argument n
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  integer :: lstr,ip,cast,narg,itrm,it(1)
  character*(*) tok,sep,strn
  double precision :: res(narg)
  ! Local Variables
  logical :: ldum,parstr
  character term*1
  integer :: jp,np,nsep,lentok,a2vec


  nsep = len(sep)
  lentok = len(tok)
  term = tok(lentok:lentok)

  ! --- Find end of string ---
  jp = ip
  it(1) = 0
  if (itrm <= nsep) &
       call chrps2(strn,sep(itrm:nsep),nsep-itrm+1,lstr,jp,it)
  if (it(1) /= 0) then
     np = jp
  else
     np = lstr
  endif

  !     print *, 'np,lstr=',np,lstr,ip,jp

  ! --- Parse for tok within string strn, returning 0 if missing  ---
  !      print *, tok
  !      print *, strn

  if (tok /= ' ') then
     if (narg == 0 .AND. np == lentok) then
        ip = np+1
        parg = 0
        if (strn(1:np) == tok(1:np)) parg = 1
        return
     elseif ( .NOT. parstr(strn,tok,np-lentok,lentok,term,ip,jp)) then
        parg = 0
        ip = np+1
        return
     endif
  else
     jp = ip
  endif

  ! --- Parse for vector of binary values to convert
  if (narg == 0) then
     ip = jp
     parg = 1
     return
  endif

  ip = jp
  parg = a2vec(strn,np,ip,cast,sep,nsep,itrm,narg,it,res)

end function parg
