! These routines are taken from Ferdi's rw.f
!-----------------------------------------------------------------
subroutine readx(ifil,n)
  integer::ifil,n,i,j
  character(72) :: rchar
  do 10 i = 1,n
     read(ifil,5000)rchar
     j       = 0
     !        call rmvbl (rchar,72,j)
     !        rchar      = rchar(j+1:72)
     rchar=trim(adjustl(rchar))

     if(rchar(1:3) == '***')return
     if(rchar(1:3) == '###')return

10 enddo
  ! top2rx 2013.08.09 kino      stop 'readx: cannot find the string(rw.f)'
  call rx( 'readx: cannot find the string(rw.f)')
5000 format(a72)
  !     return
end subroutine readx
!-------------------------------------------------------------------
subroutine rwdd (ifi, &
     ldim,n, &
     a)

  ! 92.02.07
  ! direct access read (ifi>0) or write (ifi<0)

  ! ldim = leading dimension of a

  implicit real*8  (a-h,o-z)
  integer:: ldim,n,ifi,j,i
  dimension a(ldim,n)

  ! top2rx 2013.08.09 kino      if (ifi .eq. 0) stop 'rwdd: ifi .eq. 0'
  if (ifi == 0) call rx( 'rwdd: ifi == 0')

  ! read
  if (ifi > 0) then
     do       j = 1,n
        read (ifi,rec=j) (a(i,j),i=1,ldim)
     end do
  endif

  ! write
  if (ifi < 0) then
     do       j = 1,n
        write (-ifi,rec=j) (a(i,j),i=1,ldim)
     end do
  endif

  return
end subroutine rwdd
!-------------------------------------------------------------------
subroutine rwdd1 (ifi,irec, &
     ldim, &
     a)

  ! 92.02.07
  ! direct access read (ifi>0) or write (ifi<0) for record irec

  ! irec = record number
  ! ldim = leading dimension of a

  implicit real*8  (a-h,o-z)
  integer:: ifi,irec,i,ldim
  real*8 :: a(ldim)

  ! top2rx 2013.08.09 kino      if (ifi .eq. 0) stop 'rwdd1: ifi .eq. 0'
  if (ifi == 0) call rx( 'rwdd1: ifi == 0')

  ! read
  if (ifi > 0) then
     read (ifi,rec=irec) (a(i),i=1,ldim)
  endif

  ! write
  if (ifi < 0) then
     write (-ifi,rec=irec) (a(i),i=1,ldim)
  endif

  return
end subroutine rwdd1

!--------------------------------------------------------------------
subroutine wkpnt (ifkp,qbz,wbz,nqbz)

  ! 92.04.21
  ! write k-points

  implicit real*8 (a-h,o-z)
  dimension qbz(3,nqbz),wbz(nqbz)
  integer:: ifkp,nqbz,k,i
  ! top2rx 2013.08.09 kino      if (ifkp .lt. 0) stop 'wkpnt: unit file < 0'
  if (ifkp < 0) call rx( 'wkpnt: unit file < 0')
  write (ifkp,*) ' label  k-vector  weight '
  do       k = 1,nqbz
     write (ifkp,6000)k,(qbz(i,k),i=1,3),wbz(k)
  end do

6000 format (1x,i4,4f12.6)
  return
end subroutine wkpnt
