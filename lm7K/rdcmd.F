c I (takao) think this rdcmd is probably nor for MPI|MPIK 
c This routine is kept only for tests in lm*K pacages.
c
      program main
      implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer ierr, nsize, id 
      character(14):: testbuf
      call mpi_init(ierr)
cstart--- check mpi. You can delete this. 
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      write(*,*) 'mpi: size and id=',nsize,id
cend-----
#endif      
      call fmain()
#if MPI|MPIK
      call mpi_finalize(ierr)
#endif
      end program main

C to test, run on file 'ctrl.ss'
      subroutine fmain()
      implicit none
      integer nr,recl,iarg,i1mach,ifi,i,j,mxchr,mxlev,ctlen
      parameter (recl=500,mxchr=40,mxlev=6,ctlen=120)
      character a*(recl), recrd*(recl), first*512, ctbl(mxchr,2)*(ctlen)
      logical cmdstr,a2bin,loop0(0:mxlev),noerr,noexec
      integer nlin(0:mxlev),list(100,mxlev),ilist(mxlev),nlist(0:mxlev)
      character vnam(mxlev)*16, catnam*8,shorec*(recl),fmt*(mxchr)
      character brack*2, commnt*1, commnd*1, rdparm*7, keep*1, cvar*1
C     incat: -1 => do not look for category to start
C             0 => read from category, but category not read yet
C             1 => read from category, and category has started
C             2 => category has finished
      integer ires,ipr,cn1,cn2,incat,stdo,i1,i2,j1,j2,awrite
      integer w(1000)
      common /w/ w

      commnt = '#'
      brack  = '{}'
      commnd = '%'
      keep = ' '
      cvar = 'c'
      noerr = .false.
      noexec = .false.
      incat = -1
      stdo = i1mach(2)
      fmt = '>>%f'

      call pshpr(0)
      call wkinit(1000)
      call poppr
C     call finits(2,0,0,i)

      iarg = 1
    5 if (.not. cmdstr(iarg,first)) goto 99

      if (first(1:2) .eq. '-v') then
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

      elseif (first(1:4) .eq. '--s ') then
        brack = '[]'
        iarg = iarg+1
        goto 5

      elseif (first(1:5) .eq. '-cat:') then
        catnam = first(6:)
        call word(catnam,1,cn1,cn2)
        incat = 0
        iarg = iarg+1
        goto 5

      elseif (first(1:2) .eq. '-c') then
        iarg = iarg+1
        goto 5

      elseif (first(1:4) .eq. '--pr') then
        i = 4
        if (.not. a2bin(first,j,2,0,' ',i,-1)) goto 99
        call pshprt(j)
        iarg = iarg+1
        goto 5

      elseif (first(1:7) .eq. '--noerr') then
        noerr = .true.
        iarg = iarg+1
        goto 5

      elseif (first(1:3) .eq. '--n') then
        noexec = .true.
        iarg = iarg+1
        goto 5

      elseif (first(1:3) .eq. '-f:') then
        fmt = first(4:)
        iarg = iarg+1
        goto 5

      elseif (first .eq. 'stdin' .or. first .eq. '.') then
        ifi = i1mach(1)

      else
        ifi = 10
        open(ifi, file=first, status='OLD', err=99)
      endif

C --- read and print file ---
      rdparm = commnt // brack // commnd // keep // cvar // 't'
      nr = 0
      call getpr(ipr)
   10 call rdfiln(ifi,rdparm,mxlev,loop0,nlin,list,100,ilist,nlist,
     .  vnam,ctbl,mxchr,a,recrd,recl,nr)
      if (nr .le. 0) goto 30
C     We are using lines within a category
      if (incat .eq. 0) then
        if (recrd(cn1:cn2+1) .eq. catnam(cn1:cn2+1)) incat = 1
        if (incat .eq. 0) then
          if (ipr .ge. 100/1)
     .      call awrit0('%a','>> skipping: '//recrd,recl+13,-stdo)
          goto 10
        else
          recrd(cn1:cn2) = ' '
          if (recrd .eq. ' ') goto 10
        endif
      elseif (incat .eq. 1) then
        if (recrd(1:1) .ne. ' ') goto 30
      endif
      if (ipr .gt. 20 .or. noexec) then
        shorec = ' '
        call strip(fmt,i1,i2)
        j1 = 1
        j2 = awrite(fmt(i1:i2),shorec,recl,0,0,0,0,0,0,0,0,0)
        call strip(recrd,i1,i2)
        write(shorec(1+j2:),'(a)') recrd(i1:min(i2,recl-j2))
        call strip(shorec,i1,i2)
        write(stdo,'(a)') shorec(1:i2)
        call ftflsh(stdo)
      endif
C        if (ipr .gt. 20) then
C          call awrit0('%a','>>'//recrd,recl+2,-stdo)
C          call ftflsh(stdo)
C        endif
      if (.not. noexec) then
        ires = 0
        if (recrd .ne. ' ') call fsystm(recrd,ires)
        if (ires .ne. 0 .and. noerr) then
          if (ipr .gt. 30) print 345, commnt, ires
  345     format(a,' rdcmd aborting because of nonzero exit status',i4)
          call cexit(-1,1)
        endif
      endif
      goto 10

C --- cleanup ---
   30 continue
      if (ipr .ge. 50) then
        print *, '------------------------'
        print *, 'rdcmd: variables still defined:'
        call shosyv(0,0,0,stdo)
      endif
C     call wkchk('testing')
      call cexit(0,1)

   99 continue
      print 333
  333 format('usage: rdcmd [--pr#] [--noerr] [--n] [--s]',
     .  ' [-f:string] [-cat:name]'/
     .  t13,' [-cnam=strn] [-vvar=# ...] filename'/
     .  t7,'Switches:'/
     .  t7,'--noerr',t17,'Abort when nozero exit status encountered'/
     .  t7,'--n',t17,'Print commands that would be executed, ',
     .               'but do not execute them'/
     .  t7,'--s',t17,'Use brackets [] to demarcate string substitution'/
     .  t7,'-cat:name',t17,'Execute commands within category "name"'/
     .  t7,'-f:string',t17,'Set format string for command printout')
      end
