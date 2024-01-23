subroutine cputid(ifilein) ! ifilein = file number ==> screen (id=6)
  use m_lgunit,only:stdo
  implicit none
  real(8) :: cpuetime, etw(2),cpulast, etime
  real(8) :: cpu0=-1d0
  save cpu0
  real(8):: cpusec,cpumax,cpumin
  integer:: ierr,rank ,ifile,ifilein
  logical,save::firsttime=.true.
  integer,save::i1
  integer:: i2,irate,imax
  real(8)::diff
  ifile=merge(ifilein,stdo,ifilein/=0)
  if (firsttime) then
     call system_clock(i1)
     firsttime=.false.
     cpusec=0.0d0
     cpumin=0.0d0
  else
     call system_clock(i2,irate,imax)
     diff=i2-i1
     if (diff<0) diff=imax-i1+i2
     diff=diff/dble(irate)
     cpusec=diff
     cpumin=cpusec/60.d0
  endif
  write(ifile,"('  CPU: ',f12.4,' secs = ',f7.1,' mins')")cpusec,cpumin
end subroutine cputid

subroutine cpudel(unit,strn,delt)  !- incremental cup time, in seconds
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   unit>=0: printout of time to unit; else no printing
  !i   unit<-1: delt not calculated
  !o Outputs:
  !o   delt: incremental cpu time since last call, seconds
  !r Remarks
  !r   Uses cpusec
  ! ----------------------------------------------------------------------
  character*(*) strn
  integer :: unit
  double precision :: delt
  character(1) :: timeu, outs(80)
  double precision :: cpusec,told,tnew
  save told
  data told /0d0/
  if (unit < -1) return
  tnew = cpusec()
  delt = tnew - told
  told = tnew
  timeu = 's'
  if (tnew > 60) then
     timeu = 'm'
     tnew = tnew/60
     if (tnew > 60) then
        timeu = 'h'
        tnew = tnew/60
     endif
  endif
  if (unit >= 0 .AND. tnew > 0d0) write(unit,333) strn,delt,tnew,timeu
333 format(' cpudel',a25,'  time(s):',g10.3,'  total:',f8.3,a1)
end subroutine cpudel

real(8) function cpusec()
  implicit none
  logical,save::firsttime=.true.
  integer,save::i1
  integer:: i2,irate,imax
  real(8)::diff
  if (firsttime) then
     call system_clock(i1)
     firsttime=.false.
     cpusec=0.0d0
  else
     call system_clock(i2,irate,imax)
     diff=i2-i1
     if (diff<0) diff=imax-i1+i2
     diff=diff/dble(irate)
     cpusec=diff
  endif
end function cpusec
