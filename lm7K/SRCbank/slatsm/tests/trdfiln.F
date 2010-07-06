C to test, run on file 'ctrl.ss'
      subroutine fmain
      implicit none
      integer nr,recl,iarg,i1mach,ifi,i,j,n,mxchr,mxlev,iprint,ctlen
C     double precision val
      parameter (recl=500,mxchr=40,mxlev=6,ctlen=120)
      character a*(recl), recrd*(recl), first*512, ctbl(mxchr,2)*(ctlen)
      logical cmdstr,lsequ,a2bin,loop0(0:mxlev)
      integer nlin(0:mxlev),list(100,mxlev),ilist(mxlev),nlist(0:mxlev)
      character vnam(mxlev)*16
      character brack*2, commnt*1, commnd*1, rdparm*6, keep*1, cvar*1
      integer w(1000)
      common /w/ w
      
      commnt = '#'
      brack  = '{}'
      commnd = '%'
      keep = ' '
      cvar = 'c'

      call pshpr(0)
      call wkinit(1000)
      call poppr
C     call finits(2,0,0,i)

      iarg = 1
    5 if (.not. cmdstr(iarg,first)) goto 99
      if (lsequ(first,'-v',2,' ',n)) then
        i = 2
        call parsyv(first,len(first),999,0,i)
        iarg = iarg+1
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

      if (lsequ(first,'-s ',3,' ',n)) then
        brack = '[]'
        iarg = iarg+1
        goto 5
      endif

      if (lsequ(first,'-c',2,' ',n)) then
        iarg = iarg+1
        goto 5
      endif
        
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

C --- read and print file ---
      rdparm = commnt // brack // commnd // keep // cvar
      nr = 0
   10 call rdfiln(ifi,rdparm,mxlev,loop0,nlin,list,100,ilist,nlist,
     .  vnam,ctbl,mxchr,a,recrd,recl,nr)
      if (nr .le. 0) goto 30
      call awrit0('%a',recrd,recl,-i1mach(2))
      goto 10

C --- cleanup ---
   30 continue
      if (iprint() .ge. 50) then
        print *, '------------------------'
        print *, 'variables still defined:'
        call shosyv(0,0,0,i1mach(2))
      endif
C     call wkchk('testing')
      call cexit(0,1)

   99 continue
      print 333
  333 format('usage: rdfile [-s] [-cnam=strn] [-vvar=# ...] [-pr#]',
     .  ' filename'/
     .  '       -s use brackets [] to demarcate string substitution')
      end
