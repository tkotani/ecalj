      subroutine tcinit(ltrace,lprof,levelinit)
C- Setup for timer and workspace routines
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltrace : >0 turn on running account of cpu usage
Ci            The value of ltrace is the depth to which quantities shown
Ci            <0 also turn on running account of workspace usage
Ci   lprof  : >0 print out summary of cpu and workspace usage
Ci            The value of lprof is the depth to which quantities shown
Cr Remarks
Cr   This routine logs cpu and workspace usage.
Cu Updates
Cu   30 Jun 06 ltrace <0 => running workspace printout
Cu   15 Feb 02 tcprt now writes to a logical unit, rather than stdout
C ----------------------------------------------------------------------
      implicit none
      integer ltrace,lprof,levelinit,iwnow
      integer i,ii,init,it,jcb,jlev,job,jt0,l,level,levx,lok,lpr,ltop,
     .ltr,mfree,mmax,mnow,nn,nx
      parameter (nx = 200)
      character*20 name(nx),nam0
      character*(*) str1,str2
      character*3  tt(0:20)
      integer ncall(nx),lev(nx),it1(nx),icb(nx),icall(-2:30) 
                       !bugfux this was icall(0:30) sep2010
      integer iwk0(nx),iwk1(nx),iwk2(nx),iwk3(nx)
      integer lvl,lv,last,j1,j2,iwmax,ilev,stdo,lgunit
      integer ifi,iprint
C$    integer mp_my_threadnum
      real ttot(nx),tavg,tm
      double precision cpusec
      save
      data nn/0/, ii/0/, init/0/, ltr/0/, lpr/0/ !level/0/, 

      level=levelinit

      stdo = lgunit(1)
      if (ltrace .ne. 0) ltr = ltrace
      if (lprof .gt. 0) lpr = lprof
      ltop = max0(iabs(ltr),lpr)
      jt0 = 1000*cpusec()
C     jt1 = 1000*cpusec()
      init = 1
      return

      entry tcn(str1)
C$    if(mp_my_threadnum().ne.0) return
      levx = level
      level = level+1
      if (level .gt. ltop) return
      nam0 = str1
      job = 1
      goto 95

      entry tclev(str1,ilev)
      ilev = level-1
      str1 = ' '
      if (ilev .ge. 0) str1 = name(icall(ilev))
      return

      entry tcx(str2)
C$    if(mp_my_threadnum().ne.0) return
      levx = level-1
      level = level-1
      if (level+1 .gt. ltop) return
      nam0 = str2
      job = 2

  95  continue
      if (init .eq. 0) jt0 = 1000*cpusec()
      init = 1
c      call wquery(mfree,mnow,mmax) !takao no w(x) anymore.
      mfree=0
      mnow=0
      mmax=0

c ... look for entry with same call path
      ii = 0
      do  10  i = 1, nn, 1
        if (nam0.eq.name(i) .and. levx.eq.lev(i)) then
          jcb = icb(i)
          lok = 1
          do  12  jlev = levx-1, 0, -1
            if (jcb .ne. icall(jlev)) lok = 0
            jcb = icb(jcb)
   12     continue
          if (lok .eq. 1) then
            ii = i
            goto 91
          endif
        endif
   10 continue

c ... add a new entry
      if (job .eq. 2) call rxs('tcx: missing call to tcn: ',str2)
      nn = nn+1
      if (nn .gt. nx) call rx('tcn: overflow')
      ncall(nn) = 0
      iwk0(nn)=mnow*0.001+0.5
      iwk1(nn)=-10.0
      iwk2(nn)=mmax*0.001+0.5
      iwk3(nn)=-10.0
      ttot(nn) = 0
      name(nn) = nam0
      lev(nn) = level-1
c      print *,'level-2=',level-2
      icb(nn) = icall(level-2)
      ii = nn

c ... do real operations here
  91  continue
      if (job .eq. 1) then
        it = 1000*cpusec()
        it1(ii) = it
        icall(level-1) = ii
        if (ltr .gt. levx .and. iprint() .gt. 0) then
          write(stdo,100) level,0.001*(it-jt0),str1
        elseif (iabs(ltr) .gt. levx .and. iprint() .gt. 0) then
          write(stdo,100) level,0.001*(it-jt0),str1,mnow,mmax
        endif
  100   format(' >> level:',i2,'  CPUsec=',f10.2,'  enter ',a:t43,
     .  '  mnow =',i11,'  mmax=',i11)
C       jt1 = it
        return
      else
        it = 1000*cpusec()
        tm = (it-it1(ii))/1000d0
        ncall(ii) = ncall(ii)+1
        ttot(ii) = ttot(ii)+tm
        if (iwk1(ii).lt.0.0) iwk1(ii)=mnow*0.001+0.5
        if (iwk3(ii).lt.0.0) iwk3(ii)=mmax*0.001+0.5
        if (iabs(ltr).gt.levx .and. iprint() .gt. 0) then
          if (ltr .lt. 0) then
            write(stdo,101) 0.001d0*(it-jt0),str2,tm,mnow,mmax
          elseif (iwk0(ii).ne.iwk1(ii)) then
            write(stdo,102) 0.001*(it-jt0),str2,tm,iwk1(ii)-iwk0(ii)
  102       format(' >>',f10.2,'   exit  ',a,t35,f8.2,
     .      1x,i7,'K left on stack')
          else
            write(stdo,101) 0.001d0*(it-jt0),str2,tm
  101       format(' >>',f10.2,'   exit  ',a,t35,f8.2:
     .      '  mnow =',i11,'  mmax=',i11)
          endif
        endif
        return
      endif

      entry tcprt(ifi)
      if (lpr .eq. 0) return

      if (iprint() .gt. 0) write (ifi,800) lpr
  800 format(/
     .'  ==== xxxxxxxxx ====     calls      == cpu time ===   depth',i2/
     .'  entry   xxxx  xxxx                per call  total  (depth is by TIM= in ctrl.*.)')

      it = 1000*cpusec()
      tm = (it-jt0)/1000d0
      
c      call wquery(mfree,mnow,mmax) !takao no common/w anymore
      mfree=0
      mnow=0
      mmax=0
      iwmax = 0.001*mmax+0.5
      iwnow = 0.001*mnow+0.5
      if (iprint() .gt. 0)
c     .write (ifi,700) 0,0,iwmax,1,tm,tm,'main'
     .write (ifi,700) 0,0,0,1,tm,tm,'main'
      lev(nn+1)=-1

      lvl=-1
      do  20  i = 1, nn
        lv = lev(i)
        if (lv+1 .le. lpr) then

          if (lv .gt. lvl) then
            tt(lv) = '|--'
            if (lv.gt.0 .and. tt(lv-1).eq.'|--') tt(lv-1) = '|  '
            if (lv.gt.0 .and. tt(lv-1).eq.'`--') tt(lv-1) = '   '
          endif
          if (lv .lt. lvl) tt(lv) = '|--'
          last = 1
          do  22  ii = i+1,nn+1
            if (lev(ii) .eq. lv) last = 0
            if (lev(ii) .lt. lv) goto 99
   22     continue
   99     continue
          if (last .eq. 1) tt(lv) = '`--'

          tavg = ttot(i)/max(ncall(i),1)
          call strip(name(i),j1,j2)
          if (iprint() .gt. 0)
     .    write (ifi,700) iwk0(i),iwk2(i),iwk3(i),ncall(i),
     .    tavg,ttot(i),(tt(l),l=0,lev(i)),name(i)(j1:j2)
  700     format(3i7,i9,2f11.2,3x,30a)
          lvl = lv
        endif
   20 continue

      end
      subroutine tc(string)
C- Routine for tracing program flow
C     implicit none
      integer i1mach,ltr,ltc,ncall
      character*(*) string
      double precision x
      save ltr,ncall
      data ltr/0/ ncall/0/
      ncall = ncall+1
      if (string .eq. 'on' .or. string .eq. 'ON') then
        ltr = 1
        call cpudel(i1mach(2),'set tc ...',x)
      elseif (string .eq. 'off' .or. string .eq. 'OFF') then
        ltr = 0
      elseif (string .eq. 'tog' .or. string .eq. 'TOG') then
        ltr = 1-ltr
      else
        if (ltr .eq. 1) call cpudel(i1mach(2),string,x)
      endif
      return

      entry tcget(ltc)
      ltc = ltr
      end




