subroutine reindx (noccv,nunoccv,nindxv, &
     noccc,nunoccc,nindxc, &
     nl,nn,nnv,nnc,nclass, &
     nocc,nunocc,nindx)
  ! taken from index.f
  ! 92.03.14
  ! incorporates core states into nocc,nunocc,nindx

  implicit real*8 (a-h,o-z)
  implicit integer(i-n)
  dimension noccv(0:nl-1,nnv,nclass),nunoccv(0:nl-1,nnv,nclass), &
       noccc(0:nl-1,nnc,nclass),nunoccc(0:nl-1,nnc,nclass), &
       nindxv(0:nl-1,nclass),nindxc(0:nl-1,nclass), &
       nocc(0:nl-1,nn,nclass),nunocc(0:nl-1,nn,nclass), &
       nindx(0:nl-1,nclass)

  do      ic = 1,nclass
     do       l = 0,nl-1

        ncore      = nindxc(l,ic)
        nval       = nindxv(l,ic)
        nindx(l,ic)= ncore + nval
        ! top2rx 2013.08.09 kino          if (ncore+nval .gt. nn) stop 'reindx: ncore+nval > nn'
        if (ncore+nval > nn) call rx( 'reindx: ncore+nval > nn')

        do       n = 1,ncore
           nocc(l,n,ic)   = noccc(l,n,ic)
           nunocc(l,n,ic) = nunoccc(l,n,ic)
        end do

        do       n = 1,nval
           nocc(l,ncore+n,ic)   = noccv(l,n,ic)
           nunocc(l,ncore+n,ic) = nunoccv(l,n,ic)
        end do

     end do
  end do

  return
end subroutine reindx
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! taken from basn.f
integer function nodnum(f,n)
  integer::i,n
  real(8):: f(n)
  nodnum=0
  do i=2,n-1
     if(f(i)*f(i+1)<0) nodnum=nodnum+1
  enddo
END function nodnum
!-------------------------------------------------------------------
subroutine phivc (phiv,phivd,phic, &
     nindxv,nindxc, &
     nrx,nl,nn,nnv,nnc,nclass, &
     phi,phidot)

  ! 92.03.15
  ! combines valence and core phi

  ! phiv,phivd = valence phi and phidot
  ! phic = core phi
  ! nindxv,nindxc = n index for valence and core
  ! nrx = max. no. radial points
  ! nl,nn = max. l,n
  ! nnv = max. n for valence phi
  ! nnc = max. n for core

  ! phi,phidot = phi and phidot including core
  ! phidot core  = 0

  implicit real*8 (a-h,o-z)
  implicit integer(i-n)
  dimension phiv (nrx,0:nl-1,nnv,nclass), &
       phivd(nrx,0:nl-1,nnv,nclass), &
       phic (nrx,0:nl-1,nnc,nclass), &
       nindxv(0:nl-1,nclass), &
       nindxc(0:nl-1,nclass), &
       phi(nrx,0:nl-1,nn,nclass), &
       phidot(nrx,0:nl-1,nn,nclass)

  do      ic = 1,nclass

     ! core
     do       l = 0,nl-1
        do       n = 1,nindxc(l,ic)
           do       i = 1,nrx
              phi(i,l,n,ic) = phic(i,l,n,ic)
              phidot(i,l,n,ic)= 0.d0
           end do
        end do
     end do

     ! valence
     do       l = 0,nl-1
        ncore      = nindxc(l,ic)
        nval       = nindxv(l,ic)
        ! top2rx 2013.08.09 kino          if (ncore+nval .gt. nn) stop 'phivc: ncore+nval > nn'
        if (ncore+nval > nn) call rx( 'phivc: ncore+nval > nn')
        do       n = 1,nval
           do       i = 1,nrx
              phi(i,l,ncore+n,ic) = phiv(i,l,n,ic)
              phidot(i,l,ncore+n,ic) = phivd(i,l,n,ic)
           end do
        end do
     end do

  end do

  return
end subroutine phivc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! aken from rw.f
subroutine rwphia (ifil,nindx,nl,nn, &
     z,nclass, &
     nrx, &
     a,r,phi,phidot,nrofi)

  ! read ( ifil > 0 ) and write ( ifil < 0 ) phi and phidot
  ! and radial mesh for all classes

  ! ifil  = file unit where phi and phidot are stored
  ! nindx  = no. orbitals per l channel
  ! nl     = no. l
  ! nn     = max. no. orbitals per l channel
  ! z      = atomic number, for precaution
  ! nrofi  = no. radial points, for precaution
  ! nclass = no. class
  ! nrx    = max. no. radial points

  ! r(i)   = b(exp(i-1)a -1)
  ! phi    = rR

  implicit real*8(a-h,o-z)
  implicit integer(i-n)
  dimension nindx(0:nl-1,nclass), &
       z(nclass), &
       nrofi(nclass)
  dimension a(nclass), &
       r(nrx,nclass), &
       phi(nrx,0:nl-1,nn,nclass), &
       phidot(nrx,0:nl-1,nn,nclass)

  if(ifil > 0)then
     do      ic = 1,nclass
        read(ifil)ict,zt,nrt,a(ic),b
        ! top2rx 2013.08.09 kino          if(ict .ne. ic)stop 'rwphi: wrong class'
        if(ict /= ic)call rx( 'rwphi: wrong class')
        ! top2rx 2013.08.09 kino          if(zt .ne. z(ic)) stop 'rwphi: wrong atom'
        if(zt /= z(ic)) call rx( 'rwphi: wrong atom')
        !     if(nrt.ne. nrofi(ic)) stop 'rwphi: wrong radial mesh'
        nrofi(ic)  = nrt
        ! top2rx 2013.08.09 kino          if(nrt .gt. nrx) stop 'rwphi: too many radial mesh pts'
        if(nrt > nrx) call rx( 'rwphi: too many radial mesh pts')
        read(ifil)( r(l,ic),l=1,nrt )
        do       l = 0,nl-1
           do       n = 1,nindx(l,ic)
              read(ifil)ict,lt,nt
              ! top2rx 2013.08.09 kino              if(lt .ne. l)stop 'rwphi: wrong l'
              if(lt /= l)call rx( 'rwphi: wrong l')
              ! top2rx 2013.08.09 kino              if(nt .ne. n)stop 'rwphi: wrong n'
              if(nt /= n)call rx( 'rwphi: wrong n')
              ! top2rx 2013.08.09 kino              if(n .gt. nn)stop 'rwphi: wrong nn'
              if(n > nn)call rx( 'rwphi: wrong nn')
              read(ifil)( phi(i,l,n,ic),    i=1,nrt )
              read(ifil)( phidot(i,l,n,ic), i=1,nrt )
           end do
        end do
     end do
     rewind ifil
  end if

  ! write
  if(ifil < 0)then
     do      ic = 1,nclass
        write(-ifil)ic,z(ic),nrofi(ic),a(ic)
        write(-ifil)( r(l,ic),l=1,nrofi(ic) )
        do       l = 0,nl-1
           do       n = 1,nindx(l,ic)
              write(-ifil)ic,l,n
              write(-ifil)( phi(i,l,n,ic),    i=1,nrofi(ic) )
              write(-ifil)( phidot(i,l,n,ic), i=1,nrofi(ic) )
           end do
        end do
     end do
  end if

  return
end subroutine rwphia
!-----------------------------------------------------------------------
subroutine rwphic (ifil,nindx,nl,nn, &
     z,nclass, &
     nrx, &
     a,r,phi,nrofi)

  ! 92.03.15 from rwphia
  ! read ( ifil > 0 ) and write ( ifil < 0 ) phi and phidot
  ! and radial mesh for all classes

  ! ifil  = file unit where phi and phidot are stored
  ! nindx  = no. orbitals per l channel
  ! nl     = no. l
  ! nn     = max. no. orbitals per l channel
  ! z      = atomic number, for precaution
  ! nrofi  = no. radial points, for precaution
  ! nclass = no. class
  ! nrx    = max. no. radial points

  ! r(i)   = b(exp(i-1)a -1)
  ! phi    = rR

  implicit real*8(a-h,o-z)
  implicit integer(i-n)
  dimension nindx(0:nl-1,nclass), &
       z(nclass), &
       nrofi(nclass)
  dimension a(nclass), &
       r(nrx,nclass), &
       phi(nrx,0:nl-1,nn,nclass)

  if(ifil > 0)then
     do      ic = 1,nclass
        read(ifil)ict,zt,nrt,a(ic),b
        ! top2rx 2013.08.09 kino          if(ict .ne. ic)stop 'rwphic: wrong class'
        if(ict /= ic)call rx( 'rwphic: wrong class')
        ! top2rx 2013.08.09 kino          if(zt .ne. z(ic)) stop 'rwphic: wrong atom'
        if(zt /= z(ic)) call rx( 'rwphic: wrong atom')
        ! top2rx 2013.08.09 kino          if(nrt.ne. nrofi(ic)) stop 'rwphic: wrong radial mesh'
        if(nrt /= nrofi(ic)) call rx( 'rwphic: wrong radial mesh')
        nrofi(ic)  = nrt
        ! top2rx 2013.08.09 kino          if(nrt .gt. nrx) stop 'rwphic: too many radial mesh pts'
        if(nrt > nrx) call rx( 'rwphic: too many radial mesh pts')
        read(ifil)( r(l,ic),l=1,nrt )
        do       l = 0,nl-1
           do       n = 1,nindx(l,ic)
              read(ifil)ict,lt,nt
              ! top2rx 2013.08.09 kino              if(lt .ne. l)stop 'rwphic: wrong l'
              if(lt /= l)call rx( 'rwphic: wrong l')
              ! top2rx 2013.08.09 kino              if(nt .ne. n)stop 'rwphic: wrong n'
              if(nt /= n)call rx( 'rwphic: wrong n')
              ! top2rx 2013.08.09 kino              if(n .gt. nn)stop 'rwphic: wrong nn'
              if(n > nn)call rx( 'rwphic: wrong nn')
              read(ifil)( phi(i,l,n,ic),    i=1,nrt )
              !     read(ifil)( phidot(i,l,n,ic), i=1,nrt )
           end do
        end do
     end do
     rewind ifil
  end if

  ! write
  if(ifil < 0)then
     do      ic = 1,nclass
        write(-ifil)ic,z(ic),nrofi(ic),a(ic)
        write(-ifil)( r(l,ic),l=1,nrofi(ic) )
        do       l = 0,nl-1
           do       n = 1,nindx(l,ic)
              write(-ifil)ic,l,n
              write(-ifil)( phi(i,l,n,ic),    i=1,nrofi(ic) )
              !     write(-ifil)( phidot(i,l,n,ic), i=1,nrofi(ic) )
           end do
        end do
     end do
  end if

  return
end subroutine rwphic
