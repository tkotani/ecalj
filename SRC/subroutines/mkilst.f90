subroutine mkilst(strn,nlist,list)
  !- Resolve list (ascii string) into a vector of integers
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   strn  :string holding list of integers
  !o Outputs
  !o   nlist :number of integers in list
  !o         :nlist<0 => mkilst failed to parse list
  !o   list  :list of integers
  !r Remarks
  !r   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
  !r   low:high:step
  !r   low, high, and step are integer expressions specifying the sequence
  !r     low, low+step, low+2*step, ... high.
  !r   If :step is missing, the step size defaults to 1.  If also :high
  !r   is missing,  the sequence reduces to a single integer. Thus,
  !r     '5+1'       becomes a single number, 6.
  !r     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
  !r     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
  !r   Sequences may be strung together separated by commas, eg
  !r     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
  !u Updates
  !u   02 Feb 01 strn is now a character string
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: list(*),nlist
  character*(*) strn
  integer :: it(512),iv(512),a2vec,ip,i,j,k
  ip = 0
  nlist = -1
  call skipbl(strn,len(strn),ip)
  k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
  if (k < 1) return
  if (k >= 99) call rx('mkilst: increase size of iv')
  it(k+1) = 0
  iv(k+1) = iv(k)
  ! ... loop over all iv
  nlist = 0
  i = 0
14 continue
  i = i+1
  ! ... Case iv => a single number
  if (it(i) /= 2) then
     nlist = nlist+1
     list(nlist) = iv(i)
     ! ... Case iv => n1:n2:n3
  elseif (it(i+1) == 2) then
     do  j = iv(i), iv(i+1), iv(i+2)
        nlist = nlist+1
        list(nlist) = j
     enddo
     i = i+2
     ! ... Case iv => n1:n2
  else
     do    j = iv(i), iv(i+1)
        nlist = nlist+1
        list(nlist) = j
     enddo
     i = i+1
  endif
  if (i < k) goto 14
end subroutine mkilst

! #if TEST
! subroutine fmain
!   implicit none
!   character(20) :: strn
!   integer :: nlist,list(20),i

!   strn = '                 2,1'
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)
!   strn = '                2,1 '
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)
!   strn = '             22:33:3'
!   call mkilst(strn,nlist,list)
!   print *, nlist, (list(i), i=1,nlist)

! end subroutine fmain
! #endif

