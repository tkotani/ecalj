subroutine tcinit(ltrace,lprof,levelinit)
  use m_lgunit,only:stdo
  !- Setup for timer and workspace routines
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ltrace : >0 turn on running account of cpu usage
  !i            The value of ltrace is the depth to which quantities shown
  !i            <0 also turn on running account of workspace usage
  !i   lprof  : >0 print out summary of cpu and workspace usage
  !i            The value of lprof is the depth to which quantities shown
  !r Remarks
  !r   This routine logs cpu and workspace usage.
  !u Updates
  !u   30 Jun 06 ltrace <0 => running workspace printout
  !u   15 Feb 02 tcprt now writes to a logical unit, rather than stdout
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ltrace,lprof,levelinit,iwnow
  integer :: i,ii,init,it,jcb,jlev,job,jt0,l,level,levx,lok,lpr,ltop, &
       ltr,mfree,mmax,mnow,nn,nx
  parameter (nx = 200)
  character(20) :: name(nx),nam0
  character*(*) str1,str2
  character(3) ::  tt(0:20)
  integer :: ncall(nx),lev(nx),it1(nx),icb(nx),icall(-2:30)
  integer :: iwk0(nx),iwk1(nx),iwk2(nx),iwk3(nx)
  integer :: lvl,lv,last,j1,j2,iwmax,ilev
  integer :: ifi,iprint
  real :: ttot(nx),tavg,tm
  double precision :: cpusec
  save
  data nn/0/, ii/0/, init/0/, ltr/0/, lpr/0/ !level/0/,
  level=levelinit
  !      stdo = lgunit(1)
  if (ltrace /= 0) ltr = ltrace
  if (lprof > 0) lpr = lprof
  ltop = max0(iabs(ltr),lpr)
  jt0 = 1000*cpusec()
  init = 1
  return
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccc
  entry tcn(str1)
  levx = level
  level = level+1
  if (level > ltop) return
  nam0 = str1
  job = 1
  goto 95
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccc
  entry tclev(str1,ilev)
  ilev = level-1
  str1 = ' '
  if (ilev >= 0) str1 = name(icall(ilev))
  return
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccc
  entry tcx(str2)
  levx = level-1
  level = level-1
  if (level+1 > ltop) return
  nam0 = str2
  job = 2
95 continue
  if (init == 0) jt0 = 1000*cpusec()
  init = 1
  mfree=0
  mnow=0
  mmax=0
  ! ... look for entry with same call path
  ii = 0
  do  10  i = 1, nn, 1
     if (nam0 == name(i) .AND. levx == lev(i)) then
        jcb = icb(i)
        lok = 1
        do  12  jlev = levx-1, 0, -1
           if (jcb /= icall(jlev)) lok = 0
           jcb = icb(jcb)
12      enddo
        if (lok == 1) then
           ii = i
           goto 91
        endif
     endif
10 enddo
  ! ... add a new entry
  if (job == 2) call rxs('tcx: missing call to tcn: ',str2)
  nn = nn+1
  if (nn > nx) call rx('tcn: overflow')
  ncall(nn) = 0
  iwk0(nn)=mnow*0.001+0.5
  iwk1(nn)=-10.0
  iwk2(nn)=mmax*0.001+0.5
  iwk3(nn)=-10.0
  ttot(nn) = 0
  name(nn) = nam0
  lev(nn) = level-1
  !      print *,'level-2=',level-2
  icb(nn) = icall(level-2)
  ii = nn
  ! ... do real operations here
91 continue
  if (job == 1) then
     it = 1000*cpusec()
     it1(ii) = it
     icall(level-1) = ii
     if (ltr > levx .AND. iprint() > 0) then
        write(stdo,100) level,0.001*(it-jt0),str1
     elseif (iabs(ltr) > levx .AND. iprint() > 0) then
        write(stdo,100) level,0.001*(it-jt0),str1,mnow,mmax
     endif
100  format(' >> level:',i2,'  CPUsec=',f10.2,'  enter ',a:t43, &
          '  mnow =',i11,'  mmax=',i11)
     !       jt1 = it
     return
  else
     it = 1000*cpusec()
     tm = (it-it1(ii))/1000d0
     ncall(ii) = ncall(ii)+1
     ttot(ii) = ttot(ii)+tm
     if (iwk1(ii) < 0.0) iwk1(ii)=mnow*0.001+0.5
     if (iwk3(ii) < 0.0) iwk3(ii)=mmax*0.001+0.5
     if (iabs(ltr) > levx .AND. iprint() > 0) then
        if (ltr < 0) then
           write(stdo,101) 0.001d0*(it-jt0),str2,tm,mnow,mmax
        elseif (iwk0(ii) /= iwk1(ii)) then
           write(stdo,102) 0.001*(it-jt0),str2,tm,iwk1(ii)-iwk0(ii)
102        format(' >>',f10.2,'   exit  ',a,t35,f8.2, &
                1x,i7,'K left on stack')
        else
           write(stdo,101) 0.001d0*(it-jt0),str2,tm
101        format(' >>',f10.2,'   exit  ',a,t35,f8.2:'  mnow =',i11,'  mmax=',i11)
        endif
     endif
     return
  endif
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccc
  entry tcprt(ifi)
  if (lpr == 0) return
  if (iprint() > 0) write (ifi,800) lpr
800 format(/'  ==== procid=0 ====     calls      == cpu time ===   depth',i2/ &
       '  entry   xxxx  xxxx                per call  total  (depth is by TIM= in ctrl.*.)')
  it = 1000*cpusec()
  tm = (it-jt0)/1000d0
  mfree=0
  mnow=0
  mmax=0
  iwmax = 0.001*mmax+0.5
  iwnow = 0.001*mnow+0.5
  if (iprint() > 0)write (ifi,700) 0,0,0,1,tm,tm,'main'
  lev(nn+1)=-1
  lvl=-1
  do  20  i = 1, nn
     lv = lev(i)
     if (lv+1 <= lpr) then
        if (lv > lvl) then
           tt(lv) = '|--'
           if (lv > 0 .AND. tt(lv-1) == '|--') tt(lv-1) = '|  '
           if (lv > 0 .AND. tt(lv-1) == '`--') tt(lv-1) = '   '
        endif
        if (lv < lvl) tt(lv) = '|--'
        last = 1
        do  22  ii = i+1,nn+1
           if (lev(ii) == lv) last = 0
           if (lev(ii) < lv) goto 99
22      enddo
99      continue
        if (last == 1) tt(lv) = '`--'
        tavg = ttot(i)/max(ncall(i),1)
        call strip(name(i),j1,j2)
        if (iprint() > 0) write (ifi,700) iwk0(i),iwk2(i),iwk3(i),ncall(i), &
             tavg,ttot(i),(tt(l),l=0,lev(i)),name(i)(j1:j2)
700     format(3i7,i9,2f11.2,3x,30a)
        lvl = lv
     endif
20 enddo
end subroutine tcinit
! sssssssssssssssss
subroutine tc(string)
  !- Routine for tracing program flow
  !     implicit none
  integer :: i1mach,ltr,ltc,ncall
  character*(*) string
  double precision :: x
  save ltr,ncall
  data ltr/0/ ncall/0/
  ncall = ncall+1
  if (string == 'on' .OR. string == 'ON') then
     ltr = 1
     call cpudel(i1mach(2),'set tc ...',x)
  elseif (string == 'off' .OR. string == 'OFF') then
     ltr = 0
  elseif (string == 'tog' .OR. string == 'TOG') then
     ltr = 1-ltr
  else
     if (ltr == 1) call cpudel(i1mach(2),string,x)
  endif
  return
  ! sssssssssssssssss
  entry tcget(ltc)
  ltc = ltr
end subroutine tc
