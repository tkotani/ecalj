      subroutine reindx (noccv,nunoccv,nindxv,
     i                   noccc,nunoccc,nindxc,
     d                   nl,nn,nnv,nnc,nclass,
     o                   nocc,nunocc,nindx)
c taken from index.f
c 92.03.14
c incorporates core states into nocc,nunocc,nindx

      implicit real*8 (a-h,o-z)
      dimension noccv(0:nl-1,nnv,nclass),nunoccv(0:nl-1,nnv,nclass),
     i          noccc(0:nl-1,nnc,nclass),nunoccc(0:nl-1,nnc,nclass),
     i          nindxv(0:nl-1,nclass),nindxc(0:nl-1,nclass),
     o          nocc(0:nl-1,nn,nclass),nunocc(0:nl-1,nn,nclass),
     o          nindx(0:nl-1,nclass)

      do      ic = 1,nclass
        do       l = 0,nl-1

          ncore      = nindxc(l,ic)
          nval       = nindxv(l,ic)
          nindx(l,ic)= ncore + nval
Cstop2rx 2013.08.09 kino          if (ncore+nval .gt. nn) stop 'reindx: ncore+nval > nn'
          if (ncore+nval .gt. nn) call rx( 'reindx: ncore+nval > nn')

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
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c taken from basn.f
      integer(4) function nodnum(f,n)
      real(8):: f(n)
      nodnum=0
      do i=2,n-1
        if(f(i)*f(i+1)<0) nodnum=nodnum+1
      enddo
      end
c-------------------------------------------------------------------
      subroutine phivc (phiv,phivd,phic,
     i                  nindxv,nindxc,
     d                  nrx,nl,nn,nnv,nnc,nclass,
     o                  phi,phidot)

c 92.03.15
c combines valence and core phi

c phiv,phivd = valence phi and phidot
c phic = core phi
c nindxv,nindxc = n index for valence and core
c nrx = max. no. radial points
c nl,nn = max. l,n
c nnv = max. n for valence phi
c nnc = max. n for core

c phi,phidot = phi and phidot including core
c phidot core  = 0

      implicit real*8 (a-h,o-z)
      dimension phiv (nrx,0:nl-1,nnv,nclass),
     i          phivd(nrx,0:nl-1,nnv,nclass),
     i          phic (nrx,0:nl-1,nnc,nclass),
     i          nindxv(0:nl-1,nclass),
     i          nindxc(0:nl-1,nclass),
     o          phi(nrx,0:nl-1,nn,nclass),
     o          phidot(nrx,0:nl-1,nn,nclass)

      do      ic = 1,nclass

c core
        do       l = 0,nl-1
          do       n = 1,nindxc(l,ic)
            do       i = 1,nrx
              phi(i,l,n,ic) = phic(i,l,n,ic)
              phidot(i,l,n,ic)= 0.d0
            end do
          end do
        end do

c valence
        do       l = 0,nl-1
          ncore      = nindxc(l,ic)
          nval       = nindxv(l,ic)
Cstop2rx 2013.08.09 kino          if (ncore+nval .gt. nn) stop 'phivc: ncore+nval > nn'
          if (ncore+nval .gt. nn) call rx( 'phivc: ncore+nval > nn')
          do       n = 1,nval
            do       i = 1,nrx
              phi(i,l,ncore+n,ic) = phiv(i,l,n,ic)
              phidot(i,l,ncore+n,ic) = phivd(i,l,n,ic)
            end do
          end do
        end do

      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctaken from rw.f
      subroutine rwphia (ifil,nindx,nl,nn,
     i                   z,nclass,
     i                   nrx,
     o                   a,r,phi,phidot,nrofi)

c read ( ifil > 0 ) and write ( ifil < 0 ) phi and phidot
c and radial mesh for all classes

c ifil  = file unit where phi and phidot are stored
c nindx  = no. orbitals per l channel
c nl     = no. l
c nn     = max. no. orbitals per l channel
c z      = atomic number, for precaution
c nrofi  = no. radial points, for precaution
c nclass = no. class
c nrx    = max. no. radial points

c r(i)   = b(exp(i-1)a -1)
c phi    = rR

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass),
     i          z(nclass),
     i          nrofi(nclass)
      dimension a(nclass),
     o          r(nrx,nclass),
     o          phi(nrx,0:nl-1,nn,nclass),
     o          phidot(nrx,0:nl-1,nn,nclass)

      if(ifil .gt. 0)then
        do      ic = 1,nclass
          read(ifil)ict,zt,nrt,a(ic),b
Cstop2rx 2013.08.09 kino          if(ict .ne. ic)stop 'rwphi: wrong class'
          if(ict .ne. ic)call rx( 'rwphi: wrong class')
Cstop2rx 2013.08.09 kino          if(zt .ne. z(ic)) stop 'rwphi: wrong atom'
          if(zt .ne. z(ic)) call rx( 'rwphi: wrong atom')
c     if(nrt.ne. nrofi(ic)) stop 'rwphi: wrong radial mesh'
          nrofi(ic)  = nrt
Cstop2rx 2013.08.09 kino          if(nrt .gt. nrx) stop 'rwphi: too many radial mesh pts'
          if(nrt .gt. nrx) call rx( 'rwphi: too many radial mesh pts')
          read(ifil)( r(l,ic),l=1,nrt )
          do       l = 0,nl-1
            do       n = 1,nindx(l,ic)
              read(ifil)ict,lt,nt
Cstop2rx 2013.08.09 kino              if(lt .ne. l)stop 'rwphi: wrong l'
              if(lt .ne. l)call rx( 'rwphi: wrong l')
Cstop2rx 2013.08.09 kino              if(nt .ne. n)stop 'rwphi: wrong n'
              if(nt .ne. n)call rx( 'rwphi: wrong n')
Cstop2rx 2013.08.09 kino              if(n .gt. nn)stop 'rwphi: wrong nn'
              if(n .gt. nn)call rx( 'rwphi: wrong nn')
              read(ifil)( phi(i,l,n,ic),    i=1,nrt )
              read(ifil)( phidot(i,l,n,ic), i=1,nrt )
            end do
          end do
        end do
        rewind ifil
      end if

c write
      if(ifil .lt. 0)then
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
      end
c-----------------------------------------------------------------------
      subroutine rwphic (ifil,nindx,nl,nn,
     i                   z,nclass,
     i                   nrx,
     o                   a,r,phi,nrofi)

c 92.03.15 from rwphia
c read ( ifil > 0 ) and write ( ifil < 0 ) phi and phidot
c and radial mesh for all classes

c ifil  = file unit where phi and phidot are stored
c nindx  = no. orbitals per l channel
c nl     = no. l
c nn     = max. no. orbitals per l channel
c z      = atomic number, for precaution
c nrofi  = no. radial points, for precaution
c nclass = no. class
c nrx    = max. no. radial points

c r(i)   = b(exp(i-1)a -1)
c phi    = rR

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass),
     i          z(nclass),
     i          nrofi(nclass)
      dimension a(nclass),
     o          r(nrx,nclass),
     o          phi(nrx,0:nl-1,nn,nclass)

      if(ifil .gt. 0)then
        do      ic = 1,nclass
          read(ifil)ict,zt,nrt,a(ic),b
Cstop2rx 2013.08.09 kino          if(ict .ne. ic)stop 'rwphic: wrong class'
          if(ict .ne. ic)call rx( 'rwphic: wrong class')
Cstop2rx 2013.08.09 kino          if(zt .ne. z(ic)) stop 'rwphic: wrong atom'
          if(zt .ne. z(ic)) call rx( 'rwphic: wrong atom')
Cstop2rx 2013.08.09 kino          if(nrt.ne. nrofi(ic)) stop 'rwphic: wrong radial mesh'
          if(nrt.ne. nrofi(ic)) call rx( 'rwphic: wrong radial mesh')
          nrofi(ic)  = nrt
Cstop2rx 2013.08.09 kino          if(nrt .gt. nrx) stop 'rwphic: too many radial mesh pts'
          if(nrt .gt. nrx) call rx( 'rwphic: too many radial mesh pts')
          read(ifil)( r(l,ic),l=1,nrt )
          do       l = 0,nl-1
            do       n = 1,nindx(l,ic)
              read(ifil)ict,lt,nt
Cstop2rx 2013.08.09 kino              if(lt .ne. l)stop 'rwphic: wrong l'
              if(lt .ne. l)call rx( 'rwphic: wrong l')
Cstop2rx 2013.08.09 kino              if(nt .ne. n)stop 'rwphic: wrong n'
              if(nt .ne. n)call rx( 'rwphic: wrong n')
Cstop2rx 2013.08.09 kino              if(n .gt. nn)stop 'rwphic: wrong nn'
              if(n .gt. nn)call rx( 'rwphic: wrong nn')
              read(ifil)( phi(i,l,n,ic),    i=1,nrt )
c     read(ifil)( phidot(i,l,n,ic), i=1,nrt )
            end do
          end do
        end do
        rewind ifil
      end if

c write
      if(ifil .lt. 0)then
        do      ic = 1,nclass
          write(-ifil)ic,z(ic),nrofi(ic),a(ic)
          write(-ifil)( r(l,ic),l=1,nrofi(ic) )
          do       l = 0,nl-1
            do       n = 1,nindx(l,ic)
              write(-ifil)ic,l,n
              write(-ifil)( phi(i,l,n,ic),    i=1,nrofi(ic) )
c     write(-ifil)( phidot(i,l,n,ic), i=1,nrofi(ic) )
            end do
          end do
        end do
      end if

      return
      end
