      integer function a2vec(str,lstr,ip,cast,sep,nsep,itrm,nvec,ix,res)
C- Parses ascii string for a vector of binary values
C ----------------------------------------------------------------------
C Inputs:
C   str(ip:lstr-1):string to parse, from (ip=0 for first char)
C   cast:        0=logical, 2=int, 3=real, 4=double
C   sep(1:nsep): class of chars that separate arguments
C   itrm,nvec:   index to sep which indicates last arg: parse terminates
C                when sep(i>=itrm) encountered or nvec arguments parsed.
C                nvec<0: sign used to allow nam=expr as an expression
C                which has the additional effect of loading the nam
C                with expr into the variables table.
C                itrm<0:  sign used to skip over white space before
C                parsing for next token.
C Outputs:
C   res:         binary result
C   ip:          position in str on exit
C   ix:          indices to which characters terminated expression
C   a2vec:       number of values actually converted (- if error)
C Remarks
Cu Updates
Cu   26 May 07   When last expr ends at lstr w/out term, return w/out
Cu               error only if a2bin successfully parses expr
Cu   02 Feb 01   return w/out error when last expr ends at lstr w/out term
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer ip,lstr,cast,nsep,nvec,itrm,ix(1)
      character str*(*), sep(*)*1
      double precision res(1)
C Local Variables
      integer iv,ip0,iitrm
      logical a2bin

      iitrm = iabs(itrm)

C --- For each element do ---
      a2vec = 0
      do  10  iv = 1, iabs(nvec)

C   --- Determine appropriate terminator, possibly exit ---
C   ... ip0 points to start of next token
   12   if (itrm .lt. 0) call skipbl(str,lstr,ip)
        if (ip .ge. lstr) return
        ip0 = ip
        call chrps2(str,sep,nsep,lstr,ip0,ix(iv))
        if (ip0 .ge. lstr) ix(iv) = 0
C       No terminator was found: possibly exit with error
        if (ix(iv) .eq. 0) then
C         Still ok if expr runs to last char in string with no term
          ip0 = ip
          if (.not.a2bin(str,res,cast,iv-1,sep(ix(iv)),ip0,lstr-1)) then
            a2vec = -iv
            return
          endif
        endif
        if (.not. a2bin(str,res,cast,iv-1,sep(ix(iv)),ip,lstr-1)) then
C ...     Try for nam=expr
          if (nvec .lt. 0) then
            call skipbl(str,lstr,ip)
            ip0 = ip
            call parsyv(str,lstr,1,0,ip)
            if (ip .eq. ip0) then
              a2vec = -iv
              return
            endif
C           call shosyv(0,0,0,6)
            if (ix(iv) .ge. iitrm) return
            goto 12
          endif
C ...     Error exit
          a2vec = -iv
          return
        endif
        a2vec = iv
        if (ix(iv) .ge. iitrm) return
   10 continue

      end

