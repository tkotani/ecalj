      logical function afcmd(catnam,ifi)
C- Adds strings in category catnam, file ifi to command-line 
C     implicit none
      integer ifi
      character catnam*(*)
C Local variables
      integer nr,lc,i
      character*500 s,outs
      logical lsequ
C For rdfiln
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=40,mxlev=6,lstsiz=200,ctlen=120)
      character ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),
     .ilist(mxlev),nlist(0:mxlev)
      character vnam(mxlev)*16

      call word(catnam,1,i,lc)
      afcmd = .false.
      nr = 0
   10 call rdfiln(ifi,'#{}% c',mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr, outs,s,len(s),nr)
C  10 call rdfile(ifi,'#{}% c',s,1,outs,len(s),nr)
c     if (nr .eq. 0) return
      if (nr .le. 0) call rx('afcmd: no category')
      if (.not. lsequ(catnam,s,lc,' ',i)) goto 10
C ... We found start of category.  Now load acmdop until end is found
      afcmd = .true.
   20 if (s(lc+1:) .ne. ' ') call acmdop(s(lc+1:),len(s)-lc,0)
*     call acmdop(s(lc+1:),len(s),1)
C ... For all subsequent lines, start at first
      lc = 0
      call rdfiln(ifi,'#{}% c',mxlev,loop0,nlin,list,lstsiz,
     .ilist,nlist,vnam,ctbl,mxchr, outs,s,len(s),nr)
C ... end of file ... we are done
      if (nr .le. 0) return
C ... end of category ... we are done
      if (s(1:1) .ne. ' ') return
      goto 20
      end

