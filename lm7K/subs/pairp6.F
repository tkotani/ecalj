      subroutine pairp6(nbas,npadl,npadr,iax,ntab)
C- Find effective matching pairs for padded sites
C     implicit none

      integer nbas,npadl,npadr,niax,ntab(nbas+1)
      parameter (niax=10)
      integer iax(niax,1)
      integer ip,it,jp,ib,npad,jb,jt,ipr,nbasp

      call getpr(ipr)
      nbasp = nbas+npadl+npadr
      npad = nbasp-nbas
      if (ipr .ge. 80) print 332
  332 format(' pairp6: ibp   jbpp   ib   jbp   iax(6)')

C --- For each pair, find matching pairs for padded sites ---
      do  10  ip = nbas+1, nbasp
      do  10  it = ntab(ip)+1, ntab(ip+1)
        jp = iax(2,it)
        if (jp .gt. nbasp) then
C   ...   ip derived from ib by shift as follows; see pgfpad
          if (ip .le. nbas+npadl) then
            ib = ip - nbas
          else
            ib = ip - npad
          endif
C   ...   jb (a doubly padded site) came from jp by the following shift:
          jb = jp - npad
C   ...   Sanity check
          if (ntab(ip+1)-ntab(ip) .ne. ntab(ib+1)-ntab(ib))
     .      call rx('bug in pairp6')
C   ...   Equivalent jp-ip pair is the jb-ib one:
          do  12  jt = ntab(jb)+1, ntab(jb+1)
            if (iax(2,jt) .eq. ib .and.
     .        iax(3,it) .eq. -iax(3,jt) .and.
     .        iax(4,it) .eq. -iax(4,jt) .and.
     .        iax(5,it) .eq. -iax(5,jt))  then
              iax(6,it) = jt
              if (ipr .ge. 80) print 333, iax(1,it),iax(2,it),
     .          ib,jb,iax(6,it)
  333         format(i12,3i6,i7)
              goto 13
            endif
   12     continue
          call rx('bug in pairp6')
   13     continue
        endif

   10 continue

      end

