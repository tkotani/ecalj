subroutine iosave(flg,vstrn,vars,ifi,nvario)
  use m_ftox
  !- Write line to save file
  ! ----------------------------------------------------------------
  !i Inputs
  !i   flg: 1-character label
  !i        optional second character can be a digit specifying
  !i                 precision to save output variables.
  !i                 If character is missing or not between
  !i                 '0' and '9', defaults to '6'
  !i        optional third character can be a digit specifying
  !i                 precision to save vstrn variables
  !i                 If character is missing or not between
  !i                 '0' and '9', defaults to '6'
  !i        optional fourth character, set to 'f' or 'F'
  !i                 suppresses rewinding file.
  !i   vstrn: list of variable names, separated by commas, eg
  !i          "time,mmom,etot"
  !i   vars:  list of variables corresponding to string
  !i   ifi:   file unit (<0 means output)
  !i          ifi>0 not implemented
  !i   nvario: number of variables to i/o
  !o Outputs
  !o   variables and ehk,ehf written to -ifi
  !r Remarks
  !u Updates
  !u   17 Aug 01 Added (optional) fourth character to flg
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character*(*) vstrn
  character*(*) flg
  double precision :: vars(1)
  integer :: ifi,nvario
  ! Local parameters
  logical :: leof
  integer :: nvar,i,ib
  double precision :: val
  integer :: lstr,k,k0,j,lgunit
  parameter (lstr=1024)
  character outstr*(lstr), nam*16, p1*1, p2*1
  integer::p1x,p2x
  call rxx(ifi.gt.0,'iosave:  attempt to read from save file')
  ! ... Print fractions w/out leading 0, to save space
  ib = -1
  !      call bin2a0(ib)
  !      call bin2a0(10)
  call numsyv(nvar)
  outstr = flg(1:1)
  p1 = '6'
  p2 = '6'
  ! --- Output variables up to nvario ---
  do  10  i = 4, min(nvario,nvar)
     nam = ' '
     call watsyv(nam,val,i)
     outstr=trim(outstr)//' '//trim(nam)//'='//ftom(val)
10 enddo
  ! --- Output all variables in string ---
  k0 = 0
  call skipbl(vstrn,len(vstrn),k0)
  i = 0
20 k = k0
  i = i+1
  call chrps2(vstrn,', ',2,len(vstrn),k,j)
  nam = vstrn(k0+1:k)
  outstr=trim(outstr)//' '//trim(nam)//'='//ftof(vars(i),6)
  k0 = k+1
  if (j == 1 .AND. k < len(vstrn)) goto 20
  if (leof) call poseof(-ifi)
  write(-ifi,"(a)") trim(outstr)
  !      call bin2a0(ib)
end subroutine iosave

