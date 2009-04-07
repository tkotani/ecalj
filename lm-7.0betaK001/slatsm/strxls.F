      subroutine strxls(slist,nslist,sfind,offs,noffs)
C- Converts an ascii list of parameters to a table of offsets
C ----------------------------------------------------------------
Ci slist   list of entries in strx by name
Ci nslist  number of elements in slist
Ci sfind   list of entries to find in strux, by name
Co Outputs
Co offs    list of parameters to match, indices to offsets in structure
Co noffs   number of entries to find
Co         If ith string is not matched, strxls aborts with noffs = -i
C ----------------------------------------------------------------
C     implicit none
      integer nslist, noffs, offs(1)
      character*(*) sfind, slist(nslist)
C Local variables
      integer i,iw,j1,j2

      iw = 1
      j2 = 0
      noffs = 0
      i = nslist/2
   10 j1 = j2+1
      call nword(sfind,iw,j1,j2)
      if (j1 .gt. j2) return
      call huntst(slist,nslist,sfind(j1:j2),i)
      noffs = noffs+1
      if (i .lt. 1 .or. slist(i) .ne. sfind(j1:j2)) then
        noffs = -noffs
        return
      endif
      offs(noffs) = i
      if (j2 .lt. len(sfind)) goto 10

      end

