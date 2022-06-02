program eout2
  !-------------------------------------------------------------
  ! Gather datas from files and write final output.
  !------------------------------------------------------------
  implicit none
  integer(4) :: ifx,ierr,n,i
  real(8) :: rydberg,eks, rhoexc, rhoexc2, exc_lda,exc_lda2, &
       alat,etot,ex, exxcc, exxcv, exxvv, ecorr(1000), ex_lda, ec_lda &
       ,ecut(1000),ecuts(1000),etotx

  call getexx2('TEEXXcc',exxcc,ierr)
  call getexx2('TEEXXcv',exxcv,ierr)
  call getexx2('TEEXXvv',exxvv,ierr)

  !      ifx = 101
  open(newunit=ifx,file='TEECORR2')
  n=0
  do
     call getexx3(ifx,ecorr(n+1),ecut(n+1),ecuts(n+1),ierr)
     if(ierr==1) goto 1111
     n=n+1
     !        write(6,*) ' n ecorr=',n,ecorr(n),ecut(n),ecuts(n)
  enddo
1111 continue
  close(ifx)
  !       print *,' n=',n
  !      stop 'xxxxxxxxxxxxxxxxx test end xxxxxxxxxxxxxx'
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      ifx = ifile_handle()
  open(newunit=ifx,file='LATTC',form='unformatted')
  read(ifx)alat
  close(ifx)

  open(newunit=ifx,file='ETOTLDA')
  read(ifx,*) eks
  read(ifx,*) exc_lda
  close(ifx)

  open(newunit=ifx,file='RoVxcLDA')
  read(ifx,*) exc_lda2
  read(ifx,*) ex_lda
  read(ifx,*) ec_lda
  if(exc_lda /= exc_lda2) then
     write(6,*) exc_lda, exc_lda2
     call rx( ' eout:rhoexc in ETOTLDA and RoVxcLDA are not consistent')
  endif

  open(newunit=ifx,file='ETOTeV.dat',access='append')

  Ex = exxcc + 2d0*exxcv + exxvv
  ! ccccccccccccccccccccc
  !      print *, 'ex =',ex
  ! ccccccccccccccccccc
  do i=1,n
     Etot = eks*rydberg() - exc_lda*rydberg() + ex + ecorr(i)
     Etotx= eks*rydberg() - exc_lda*rydberg() + ex
     if(ecut(i) == 1d10) then
        write(ifx,"(f12.6, 2f18.6,f16.6,f15.6, f18.6,f16.6,f15.6, &
             ' | alat   Etot Etotx Ex Ec  EtotLDA ExLDA EcLDA')") &
             alat, Etot, Etotx, ex, ecorr(i), &
             Eks*rydberg(), ex_lda*rydberg(),ec_lda*rydberg()
     else
        write(ifx,"(f12.6, 2f18.6,f16.6,f15.6, f18.6,f16.6,f15.6,5x,2f8.4, &
             ' | alat   Etot Etotx Ex Ec  EtotLDA ExLDA EcLDA ecut ecuts')") &
             alat, Etot, Etotx, ex, ecorr(i), &
             Eks*rydberg(), ex_lda*rydberg(),ec_lda*rydberg(),ecut(i),ecuts(i)
     endif
  enddo
  close(ifx)
  ! top2rx 2013.08.09 kino      stop ' OK!: eout. Add data to ETOTeV.dat'
  call rx0( ' OK!: eout. Add data to ETOTeV.dat')
END PROGRAM eout2

!------------------------------------
subroutine getexx2(tefil,exx,ierr)
  character*(*) tefil
  integer(4):: ifx,ierr
  real(8):: exx
  open(newunit=ifx, file =tefil)
  call readxx3(ifx,ierr)
  read (ifx,*) exx
  close(ifx)
end subroutine getexx2
!------------------------------------
subroutine getexx3(ifx,exx,ecut,ecuts,ierr)
  !     character*(*) tefil
  implicit none
  integer(4):: ifx,ierr
  real(8):: exx,a2,a3,ecut,ecuts
  call readxx3(ifx,ierr)
  if(ierr==1) return
  read (ifx,*) exx,a2,a3,ecut,ecuts
  !      close(ifx)
end subroutine getexx3

! This is taken from Ferdi's rw.f
!----------------------------------------
subroutine readxx3(ifil,ierr)
  implicit none
  character(72) :: rchar
  integer(4):: n=1000,i,j,ierr,ifil
  ierr=0
  do 10 i = 1,n
     read(ifil,'(a)',end=1011,err=1011) rchar
     !      print *, rchar
     j       = 0
     call rmvbl (rchar,72,j)
     rchar      = rchar(j+1:72)
     if(rchar(1:3) == '***') return
10 enddo
1011 continue
  ierr=1
  !     stop 'readx: cannot find the string'
end subroutine readxx3
subroutine rmvbl(t,nt,i)
  !- Parses string T(I) for blanks
  integer :: nt,i
  character(1) :: t(0:nt)
99 if (t(i) /= ' ') return
  i = i + 1
  if (i >= nt) return
  goto 99
end subroutine rmvbl
!-------------------------------------------------------------------
