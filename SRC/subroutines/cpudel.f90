subroutine cpudel(unit,strn,delt)
  !- incremental cup time, in seconds
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
  if (unit >= 0 .AND. tnew > 0d0) &
       write(unit,333) strn,delt,tnew,timeu
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
     ! cpumin=0.0d0
  else
     call system_clock(i2,irate,imax)
     diff=i2-i1
     if (diff<0) diff=imax-i1+i2
     diff=diff/dble(irate)
     cpusec=diff
     ! cpumin=cpusec/60.d0
  endif
end function cpusec
