integer function a2vec(str,lstr,ip,cast,sep,nsep,itrm,nvec,ix,res)
  !- Parses ascii string for a vector of binary values
  ! ----------------------------------------------------------------------
  ! Inputs:
  !   str(ip:lstr-1):string to parse, from (ip=0 for first char)
  !   cast:        0=logical, 2=int, 3=real, 4=double
  !   sep(1:nsep): class of chars that separate arguments
  !   itrm,nvec:   index to sep which indicates last arg: parse terminates
  !                when sep(i>=itrm) encountered or nvec arguments parsed.
  !                nvec<0: sign used to allow nam=expr as an expression
  !                which has the additional effect of loading the nam
  !                with expr into the variables table.
  !                itrm<0:  sign used to skip over white space before
  !                parsing for next token.
  ! Outputs:
  !   res:         binary result
  !   ip:          position in str on exit
  !   ix:          indices to which characters terminated expression
  !   a2vec:       number of values actually converted (- if error)
  ! Remarks
  !u Updates
  !u   26 May 07   When last expr ends at lstr w/out term, return w/out
  !u               error only if a2bin successfully parses expr
  !u   02 Feb 01   return w/out error when last expr ends at lstr w/out term
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  integer :: ip,lstr,cast,nsep,nvec,itrm,ix(1)
  character str*(*), sep(*)*1
  character c_sep
  double precision :: res(1)
  ! Local Variables
  integer :: iv,ip0,iitrm
  logical :: a2bin

  iitrm = iabs(itrm)

  ! --- For each element do ---
  a2vec = 0
  do  10  iv = 1, iabs(nvec)

     !   --- Determine appropriate terminator, possibly exit ---
     !   ... ip0 points to start of next token
12   if (itrm < 0) call skipbl(str,lstr,ip)
     if (ip >= lstr) return
     ip0 = ip
     call chrps2(str,sep,nsep,lstr,ip0,ix(iv))
     if (ip0 >= lstr) ix(iv) = 0
     !       No terminator was found: possibly exit with error
     if (ix(iv) == 0) then
        c_sep=' ' ! any value is OK.
     else
        c_sep=sep(ix(iv))
     endif
     if (ix(iv) == 0) then
        !         Still ok if expr runs to last char in string with no term
        ip0 = ip
        ! i          if (.not.a2bin(str,res,cast,iv-1,sep(ix(iv)),ip0,lstr-1)) then
        if ( .NOT. a2bin(str,res,cast,iv-1,c_sep,ip0,lstr-1)) then
           a2vec = -iv
           return
        endif
     endif
     ! i        if (.not. a2bin(str,res,cast,iv-1,sep(ix(iv)),ip,lstr-1)) then
     if ( .NOT. a2bin(str,res,cast,iv-1,c_sep,ip,lstr-1)) then
        ! ...     Try for nam=expr
        if (nvec < 0) then
           call skipbl(str,lstr,ip)
           ip0 = ip
           call parsyv(str,lstr,1,0,ip)
           if (ip == ip0) then
              a2vec = -iv
              return
           endif
           !           call shosyv(0,0,0,6)
           if (ix(iv) >= iitrm) return
           goto 12
        endif
        ! ...     Error exit
        a2vec = -iv
        return
     endif
     a2vec = iv
     if (ix(iv) >= iitrm) return
10 enddo

end function a2vec

