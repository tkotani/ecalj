      subroutine fmain
      implicit none
      integer nr,recl,iarg,i1mach,ifi,i,j,n,mxchr,mxlev,mxrecs,iprint
      double precision val
      parameter (recl=72,mxchr=20,mxlev=4,mxrecs=500)
      character a*(recl), recrd(mxrecs*recl),first*40,ctbl(mxchr,2)*32
      logical cmdstr,lsequ,a2bin,loop0(0:mxlev)
      integer nlin(0:mxlev),list(100,mxlev),ilist(mxlev),nlist(0:mxlev)
      character vnam(mxlev)*16
      character brack*2, commnt*1, commnd*1, rdparm*4
      
      commnt = '#'
      brack  = '{}'
      commnd = '%'

      iarg = 1
    5 if (.not. cmdstr(iarg,first)) goto 99

      if (lsequ(first,'-v',2,' ',n)) then
        i = 2
        call parsyv(first,len(first),999,0,i)
C        j = 2
C        first = first // ' '
C        call chrpos(first,'=',19,j)
C        n = j+1
C        if (a2bin(first,val,4,0,' ',n,-1)) then
C          call addsyv(first(3:j),val,n)
C        else
C          goto 99
C        endif
C        iarg = iarg+1
        goto 5
      endif
      if (lsequ(first,'-b',3,' ',n)) then
      if (lsequ(first,'-pr',3,' ',n)) then
        i = 3
        if (.not. a2bin(first,j,2,0,' ',i,-1)) goto 99
        call pshprt(j)
        iarg = iarg+1
        goto 5
      endif

      if (first .eq. 'stdin' .or. first .eq. '.') then
        ifi = i1mach(1)
      else
        ifi = 10
        open(ifi, file=first, status='OLD', err=99)
      endif

      rdparm = commnt // brack // commnd
      call rdfile(ifi,rdparm,recrd,mxrecs,a,recl,nr)
      if (iprint() .le. 110) print 20, (recrd(j), j= 1, recl*nr)
   20 format(72a1)

      call cexit(0,1)

   99 continue
      print 333
  333 format('usage: rdfile [-b -vvar=# ...] [-pr#] filename'/
     .       '       -s use [] to demarcate string substitution')
      end
