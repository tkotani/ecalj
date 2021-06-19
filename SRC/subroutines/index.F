      integer function nocc (ek,ef,locc,
     d                       n)
c ek   = eigenvalues
c ef   = fermi level
c locc = .true.  ==> occupied
c      = .false. ==> unoccupied
c n    = no. states
      implicit none
      logical, intent(in) :: locc
      integer, intent(in) :: n
      real(8), intent(in) :: ek(n), ef
      integer :: i
      nocc       = 0
      if (locc) then
        do       i = 1,n
          if (ek(i) .le. ef) nocc = nocc + 1
        end do
      else
        do       i = 1,n
          if (ek(i) .gt. ef) nocc = nocc + 1
        end do
      endif
      return
      end function nocc

      subroutine maxdim (nocc,nunocc,nindx,nl,nn,nclass,
     o nprodx,nlnx,nlnmx,nlnax)
c largest number of product basis, (l,n) and (l,n,m)
      implicit real*8(a-h,o-z)
      integer:: nocc(nl*nn,nclass),
     i          nunocc(nl*nn,nclass),
     i     nindx(0:nl-1,nclass),ic,nclass,nn,nl,nln,nlnm,nlna,
     i     nprodx,nlnmx,nlnx,nlnax,nprod,nallow,noflnm,nalwln,nofln
      nprodx     = 0
      nlnx       = 0
      nlnmx      = 0
      nlnax      = 0
      do      ic = 1,nclass
        nprod      = nallow (nocc(1,ic),nunocc(1,ic),
     .                     nindx(0,ic),nl,nn)
        nln        = nofln(nindx(0,ic),nl)
        nlnm       = noflnm(nindx(0,ic),nl)
        nlna       = nalwln (nocc(1,ic),nunocc(1,ic),
     .                     nindx(0,ic),nl,nn)
        if(nprod .gt. nprodx) nprodx = nprod
        if(nln   .gt. nlnx)   nlnx   = nln
        if(nlnm  .gt. nlnmx)  nlnmx  = nlnm
        if(nlna  .gt. nlnax)  nlnax  = nlna
      end do
      return
      end

      integer function noccx1 (ekt,nk,nt,ef)
c 92.02.27
c find the highest number of occupied states
c ekt   = eigenvalues for all k-points and states
c nk    = number of k-points in the 1st BZ
c nt    = number of states
c ef    = fermi level
c noccx1= highest number of occupied states
      implicit real*8(a-h,o-z)
      dimension ekt(nt,nk)
      logical:: allocc
      integer:: noccx,k,it,nk,nt
!! Find the highest number of occupied states ---
!! mar2016. (Seungwoo)
!!   Previous version works,But this is necessary when all bands occuied (for model calculaitons). -->I had bug -->see 
!! aug2017 below (okumura test)
!!   bugfix (probably affects little only for metallic case). For single band model, it stops since previous version givs noccx1=0.
!!
      noccx = 0
      allocc=.true.
      do   k  = 1,nk !k index
        do it = 1,nt !band index
          if(ekt(it,k) > ef) then
             allocc=.false.
            if(it-1> noccx) noccx = it-1 !bugfix aug2017takao => previous verion check noccx only at it=nt because of bug.
            exit
          endif
        enddo
      enddo
      if(allocc) noccx=nt   !completely occupied case. aug2017takao
      noccx1 = noccx
      return
      end

      subroutine nolnma  (nindx,nl,nclass,
     o nlnm )
c 92.jan.07
c number of l,n,m for all classes
      implicit real*8(a-h,o-z)
      integer:: ic,nclass,noflnm,l,nl,
     &     nindx(0:nl-1,nclass), nlnm(nclass)
      do     ic = 1,nclass
        noflnm    = 0
        do 1    l = 0,nl-1
          noflnm    = noflnm + nindx(l,ic)*(2*l+1)
    1   continue
        nlnm(ic)  = noflnm
      end do
      return
      end

      integer function maxocc2 (nspin,ef, nband, qbz,nqbz)
c maximum no. occupied states
      use m_readeigen, only: readeval
      implicit none
      integer(4):: nspin,nqbz,nband,noccx,is,iq,noccxt,noccx1
      real(8) :: qbz(3,nqbz),ef
      real(8),allocatable :: ekt(:,:)
      allocate( ekt(nband, nqbz ) )
      noccx      = 0
      do  is = 1,nspin
        do iq = 1,nqbz
          ekt(:,iq)= readeval(qbz(:,iq),is)
        enddo
        noccxt  = noccx1 (ekt,nqbz,nband,ef)
        if (noccxt .gt. noccx) noccx = noccxt
      enddo
      maxocc2 = noccx
      deallocate( ekt )
      end

c$$$c------------------------------------------------------------------
c$$$
c$$$
c$$$      subroutine idxlmto (nindxv,iclass,
c$$$     d                    nl,nnv,nlmto,natom,nclass,
c$$$     o il,in,im,ilmto)
c$$$
c$$$c 92.03.14
c$$$c 92.10.13 from idxlmto (see above) modified for all atoms
c$$$c indexing of LMTO basis functions
c$$$c follows that in TB-LMTO program
c$$$c
c$$$c l    0     1                 2
c$$$c n    1  2  1        2        1               2
c$$$c m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
c$$$c ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
c$$$c
c$$$c il(ind),in(ind),im(ind) = l,n,m
c$$$c ilmto(n,lm,atom) = index of lmto
c$$$c lm = l*l + l + m + 1
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindxv(0:nl-1,nclass),iclass(natom),
c$$$     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto)
c$$$
c$$$      ind       = 0
c$$$      do 1   ia = 1,natom
c$$$        ic        = iclass(ia)
c$$$      do 1    l = 0,nl-1
c$$$        l2        = l*l
c$$$      do 1    n = 1,nindxv(l,ic)
c$$$      do 1    m = -l,l
c$$$        ind       = ind + 1
c$$$        lm        = l2 + l + m + 1
c$$$        il(ind)   = l
c$$$        in(ind)   = n
c$$$        im(ind)   = m
c$$$        ilmto(n,lm,ia) = ind
c$$$    1 continue
c$$$Cstop2rx 2013.08.09 kino      if (ind .ne. nlmto) stop 'idxlmto: wrong no. l,n,m'
c$$$      if (ind .ne. nlmto) call rx( 'idxlmto: wrong no. l,n,m')
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine idxlmt1 (nindxv,iclass,
c$$$     d                    nl,nnv,nlmto,natom,nclass,
c$$$     o il,in,im,ilmto,ialmto)
c$$$
c$$$c 93.08.?? from idxlmto (see above) modified to include atom index
c$$$c indexing of LMTO basis functions
c$$$c follows that in TB-LMTO program
c$$$c
c$$$c l    0     1                 2
c$$$c n    1  2  1        2        1               2
c$$$c m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
c$$$c ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
c$$$c
c$$$c il(ind),in(ind),im(ind) = l,n,m
c$$$c ilmto(n,lm,atom) = index of lmto
c$$$c ialmto(ind) = ia
c$$$c lm = l*l + l + m + 1
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindxv(0:nl-1,nclass),iclass(natom),
c$$$     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto),
c$$$     o          ialmto(nlmto)
c$$$
c$$$      ind       = 0
c$$$      do 1   ia = 1,natom
c$$$        ic        = iclass(ia)
c$$$      do 1    l = 0,nl-1
c$$$        l2        = l*l
c$$$      do 1    n = 1,nindxv(l,ic)
c$$$      do 1    m = -l,l
c$$$        ind       = ind + 1
c$$$        lm        = l2 + l + m + 1
c$$$        il(ind)   = l
c$$$        in(ind)   = n
c$$$        im(ind)   = m
c$$$        ialmto(ind) = ia
c$$$        ilmto(n,lm,ia) = ind
c$$$    1 continue
c$$$Cstop2rx 2013.08.09 kino      if (ind .ne. nlmto) stop 'idxlmto: wrong no. l,n,m'
c$$$      if (ind .ne. nlmto) call rx( 'idxlmto: wrong no. l,n,m')
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine idxlna (nindx,
c$$$     d                   nl,nn,nclass,nlnx,
c$$$     o                   ilr,inr,ilnr,nlnr)
c$$$
c$$$c 92.jan.13
c$$$c index of l,n for radial phi(l,n) for all classes
c$$$
c$$$c nlnx = max. number of l,n  from maxdim.f
c$$$
c$$$c ilr(index) = l
c$$$c inr(index) = n
c$$$c ilnr(l,n)  = index
c$$$c nlnr       = no. l,n = max(index)
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindx(0:nl-1,nclass)
c$$$      dimension ilr(nlnx,nclass),inr(nlnx,nclass),
c$$$     o          ilnr(0:nl-1,nn,nclass),nlnr(nclass)
c$$$
c$$$      do      ic = 1,nclass
c$$$        i          = 0
c$$$        do       l = 0,nl-1
c$$$          do       n = 1,nindx(l,ic)
c$$$            i          = i + 1
c$$$            ilr(i,ic)   = l
c$$$            inr(i,ic)   = n
c$$$            ilnr(l,n,ic)= i
c$$$          end do
c$$$        end do
c$$$Cstop2rx 2013.08.09 kino        if(i .gt. nlnx)stop 'idxlna: nln .gt. nlnx'
c$$$        if(i .gt. nlnx)call rx( 'idxlna: nln .gt. nlnx')
c$$$        nlnr(ic)    = i
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine idxlnp (nocc,nunocc,nindx,
c$$$     i                   nl,nn,nclass,nlnax,
c$$$     o                   il1a,in1a,il2a,in2a,ilna,nlna)
c$$$
c$$$c gives an index for the allowed product phi
c$$$
c$$$c nocc(l,n) = 0,1 ==> unoccupied, occupied
c$$$c nlna      = number of allowed product phi,
c$$$c             from integer function nalwln
c$$$
c$$$c in1(lbas) = n1 etc.
c$$$c ilna(l1,n1,l2,n2) = index of allowed product phi
c$$$c ilna(l1,n1,l2,n2) = 0 ==> not allowed
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
c$$$     i          nindx(nl,nclass)
c$$$      dimension il1a(nlnax,nclass),in1a(nlnax,nclass),
c$$$     o          il2a(nlnax,nclass),in2a(nlnax,nclass),
c$$$     o          ilna(nl*nn*nl*nn,nclass),nlna(nclass)
c$$$
c$$$      do      ic = 1,nclass
c$$$        nlna(ic)   = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)
c$$$
c$$$        call indxlna (nocc(1,ic),nunocc(1,ic),nindx(1,ic),
c$$$     d              nl,nn,nlna(ic),
c$$$     o              il1a(1,ic),in1a(1,ic),il2a(1,ic),in2a(1,ic),
c$$$     o              ilna(1,ic) )
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c---------------------------------------------------------------
c$$$c---------------------------------------------------------------
c$$$      subroutine indxbas(nocc,nunocc,nindx,nl,nn,
c$$$     i                   nallow,
c$$$     o                   il1,in1,im1,
c$$$     o                   il2,in2,im2,indxb)
c$$$
c$$$c gives an index for the allowed product basis
c$$$
c$$$c nocc(l,n) = 0,1 ==> unoccupied, occupied
c$$$c nallow    = number of allowed product basis,
c$$$c             from integer function nallow,
c$$$c             must be equal to final lbas
c$$$
c$$$c in1(lbas) = n1 etc.
c$$$c indxb(n1,lm1,n2,lm2) = index of allowed product basis
c$$$c indxb(n1,lm1,n2,lm2) = 0 ==> not allowed
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
c$$$     i          nindx(0:nl-1)
c$$$      dimension indxb(nn,nl*nl,nn,nl*nl),
c$$$     o          in1(nallow),il1(nallow),im1(nallow),
c$$$     o          in2(nallow),il2(nallow),im2(nallow)
c$$$
c$$$c$$$c write to  file PRODUCT if it exists
c$$$c$$$      ifprod     = ifile('PRODUCT')
c$$$c$$$      if (ifprod .gt. 0)
c$$$c$$$     .write (ifprod,*)'lbas,n1,l1,im1,n2,l2,im2'
c$$$
c$$$c initialise indxb
c$$$c      call iinit   (indxb,nn*nn*nl**4)
c$$$      indxb = 0
c$$$      lbas       = 0
c$$$      do 10   l1 = 0,nl-1
c$$$        ll1        = l1*l1
c$$$      do 10   n1 = 1,nindx(l1)
c$$$      do 10   m1 = 1,2*l1+1
c$$$      do 10   l2 = 0,nl-1
c$$$        ll2        = l2*l2
c$$$      do 10   n2 = 1,nindx(l2)
c$$$      do 10   m2 = 1,2*l2+1
c$$$c     indxb(n1,ll1+m1,n2,ll2+m2) = 0
c$$$        if (nocc(l1,n1) .eq. 0)   goto 10
c$$$        if (nunocc(l2,n2) .eq. 0) goto 10
c$$$        if (l1 .ne. l2 .or. n1 .ne. n2) then
c$$$          if (indxb(n2,ll2+m2,n1,ll1+m1) .ne. 0) goto 10
c$$$        endif
c$$$
c$$$c temporary
c$$$        if (l1 .eq. l2 .and. n1 .eq. n2) then
c$$$          if (indxb(n2,ll2+m2,n1,ll1+m1) .ne. 0) goto 10
c$$$        endif
c$$$
c$$$        lbas       = lbas + 1
c$$$        indxb(n1,ll1+m1,n2,ll2+m2) = lbas
c$$$        in1(lbas)  = n1
c$$$        il1(lbas)  = l1
c$$$        im1(lbas)  = m1 - l1 - 1
c$$$        in2(lbas)  = n2
c$$$        il2(lbas)  = l2
c$$$        im2(lbas)  = m2 - l2 - 1
c$$$c$$$        if (ifprod .gt. 0)
c$$$c$$$     .write (ifprod,*)lbas,n1,l1,im1(lbas),n2,l2,im2(lbas)
c$$$   10 continue
c$$$Cstop2rx 2013.08.09 kino      if(lbas .ne. nallow)stop 'indxbas: lbas .ne. nallow'
c$$$      if(lbas .ne. nallow)call rx( 'indxbas: lbas .ne. nallow')
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine indxk(n1,n2,n3,qbas,
c$$$     i                 qbz,nqbz,
c$$$     o                 ginv,indk)
c$$$
c$$$c generates indexing of k-points inside 1st BZ
c$$$c given k inside BZ indk(ia,ib,ic) = k-point inside 1st BZ
c$$$c where ia=a*n1, ib=b*n2, ic=c*n3 and
c$$$c k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors
c$$$
c$$$c n1,n2,n3 = number of divisions along G1,G2,G3
c$$$c qbz      = k-points cartesian coordinates in the 1st BZ
c$$$c nqbz     = number of k-points
c$$$c qbas     = reciprocal base vectors
c$$$
c$$$c indk     = s. above
c$$$c ginv     = inverse of base reciprocal vector
c$$$c notes:   all vectors are in units of 2*pi/a
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension qbz(3,nqbz),qbas(3,3)
c$$$      dimension ginv(3,3),
c$$$     o          indk(-n1:n1,-n2:n2,-n3:n3)
c$$$      dimension work(3)
c$$$
c$$$c write to unit ifkp
c$$$c     ifkp       = ifile('KPNT')
c$$$c     if(ifkp .gt. 0)then
c$$$c     write(ifkp,*)' k-points in the 1st BZ '
c$$$c     write(ifkp,*)' nG1,nG2,nG3, qx,qy,qz, index '
c$$$c     endif
c$$$
c$$$c initialise indk
c$$$c      call iinit(indk,(2*n1+1)*(2*n2+1)*(2*n3+1))
c$$$      indk=0
c$$$c invert qbas
c$$$      call dinv33(qbas,0,ginv,det)
c$$$
c$$$      do       k = 1,nqbz
c$$$        call dmv(3,3,ginv,3,qbz(1,k),work)
c$$$        ia         = idnint(n1*work(1))
c$$$        ib         = idnint(n2*work(2))
c$$$        ic         = idnint(n3*work(3))
c$$$        indk(ia,ib,ic) = k
c$$$c     if (ifkp .gt. 0)
c$$$c    .write(ifkp,6000)ia,ib,ic,qbz(1,k),qbz(2,k),qbz(3,k),
c$$$c    .indk(ia,ib,ic)
c$$$c6000 format(3i4,3f10.4,i5)
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine indxlm(nl,
c$$$     o ila,ima)
c$$$
c$$$c index of l,m
c$$$c i = (l+1)**2 + m - l
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension ila(nl*nl),ima(nl*nl)
c$$$
c$$$      i          = 0
c$$$      do       l = 0,nl-1
c$$$        do       m = -l,l
c$$$          i          = i + 1
c$$$          ila(i)     = l
c$$$          ima(i)     = m
c$$$        end do
c$$$      end do
c$$$Cstop2rx 2013.08.09 kino      if(i .ne. nl*nl)stop 'indxlm: wrong no. of l,m'
c$$$      if(i .ne. nl*nl)call rx( 'indxlm: wrong no. of l,m')
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine indxln(nindx,nl,
c$$$     i                  nn,nofln,
c$$$     o il,in,iln)
c$$$
c$$$c index of l,n
c$$$c nofln = number of l,n  from integer function nofln
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindx(0:nl-1)
c$$$      dimension il(nofln),in(nofln),iln(0:nl-1,nn)
c$$$
c$$$      i          = 0
c$$$      do       l = 0,nl-1
c$$$        do       n = 1,nindx(l)
c$$$          i          = i + 1
c$$$          il(i)      = l
c$$$          in(i)      = n
c$$$          iln(l,n)   = i
c$$$        end do
c$$$      end do
c$$$Cstop2rx 2013.08.09 kino      if(i .ne. nofln)stop 'indxln: wrong no. of ln'
c$$$      if(i .ne. nofln)call rx( 'indxln: wrong no. of ln')
c$$$
c$$$      return
c$$$      end
c$$$c---------------------------------------------------------------
c$$$      subroutine indxlna(nocc,nunocc,nindx,nl,nn,
c$$$     i                   nlna,
c$$$     o                   il1a,in1a,il2a,in2a,ilna)
c$$$
c$$$c gives an index for the allowed product phi
c$$$
c$$$c nocc(l,n) = 0,1 ==> unoccupied, occupied
c$$$c nlna      = number of allowed product phi,
c$$$c             from integer function nalwln
c$$$
c$$$c in1(lbas) = n1 etc.
c$$$c ilna(l1,n1,l2,n2) = index of allowed product phi
c$$$c ilna(l1,n1,l2,n2) = 0 ==> not allowed
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
c$$$     i          nindx(0:nl-1)
c$$$      dimension il1a(nlna),in1a(nlna),
c$$$     o          il2a(nlna),in2a(nlna),ilna(0:nl-1,nn,0:nl-1,nn)
c$$$
c$$$c      call iinit   (ilna,nl*nl*nn*nn)
c$$$      ilna=0
c$$$      lbas       = 0
c$$$      do 10   l1 = 0,nl-1
c$$$      do 10   n1 = 1,nindx(l1)
c$$$      do 10   l2 = 0,nl-1
c$$$      do 10   n2 = 1,nindx(l2)
c$$$c     ilna(l1,n1,l2,n2) = 0
c$$$        if(nocc(l1,n1) .eq. 0)goto 10
c$$$        if(nunocc(l2,n2) .eq. 0)goto 10
c$$$        if((l1.ne.l2 .or. n1.ne.n2) .and. ilna(l2,n2,l1,n1) .ne. 0)goto 10
c$$$        lbas       = lbas + 1
c$$$        il1a(lbas) = l1
c$$$        in1a(lbas) = n1
c$$$        il2a(lbas) = l2
c$$$        in2a(lbas) = n2
c$$$        ilna(l1,n1,l2,n2) = lbas
c$$$   10 continue
c$$$Cstop2rx 2013.08.09 kino      if(lbas .ne. nlna)stop 'indxlna: lbas .ne. nlna'
c$$$      if(lbas .ne. nlna)call rx( 'indxlna: lbas .ne. nlna')
c$$$
c$$$      return
c$$$      end
c$$$c--------------------------------------------------------------------
c$$$      subroutine indxlnm(nindx,nl,
c$$$     i                   nn,noflnm,
c$$$     o il,in,im,ilnm)
c$$$
c$$$c indexing of LMTO basis functions for a given class,
c$$$c follows that in TB-LMTO program
c$$$c il,in,im = l,n,m
c$$$c ilnm(n,lm) = index of n,l,m
c$$$c lm = l*l + l + m + 1
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindx(0:nl-1)
c$$$      dimension ilnm(nn,nl*nl),in(noflnm),il(noflnm),im(noflnm)
c$$$
c$$$      ind       = 0
c$$$      do 1    l = 0,nl-1
c$$$        l2        = l*l
c$$$      do 1    n = 1,nindx(l)
c$$$      do 1    m = 1,2*l+1
c$$$        ind       = ind + 1
c$$$        lm        = l2 + m
c$$$        il(ind)   = l
c$$$        in(ind)   = n
c$$$        im(ind)   = m - l - 1
c$$$        ilnm(n,lm) = ind
c$$$    1 continue
c$$$Cstop2rx 2013.08.09 kino      if (ind .ne. noflnm) stop 'indxlnm: wrong no. l,n,m'
c$$$      if (ind .ne. noflnm) call rx( 'indxlnm: wrong no. l,n,m')
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine idxlnma(nindx,nl,nclass,
c$$$     i                   nn,nlnmx,
c$$$     o il,in,im,ilnm)
c$$$
c$$$c 92.jan.07
c$$$c indexing of LMTO basis functions for all classes,
c$$$c follows that in TB-LMTO program
c$$$c il,in,im = l,n,m
c$$$c ilnm(n,lm) = index of n,l,m
c$$$c lm = l*l + l + m + 1
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindx(0:nl-1,nclass)
c$$$      dimension ilnm(nn,nl*nl,nclass),
c$$$     o          in(nlnmx,nclass),
c$$$     o          il(nlnmx,nclass),
c$$$     o          im(nlnmx,nclass)
c$$$
c$$$      do     ic = 1,nclass
c$$$        ind       = 0
c$$$        do 1    l = 0,nl-1
c$$$          l2        = l*l
c$$$        do 1    n = 1,nindx(l,ic)
c$$$        do 1    m = 1,2*l+1
c$$$          ind       = ind + 1
c$$$          lm        = l2 + m
c$$$          il(ind,ic)= l
c$$$          in(ind,ic)= n
c$$$          im(ind,ic)= m - l - 1
c$$$          ilnm(n,lm,ic) = ind
c$$$    1   continue
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$c-------------------------------------------------------------------
c$$$      subroutine indxmdm (mdim,nclass,
c$$$     i                    iclass,natom,
c$$$     o imdim )
c$$$
c$$$c given an atom, imdim points to the starting address of the
c$$$c total optimal product basis
c$$$c   1 2 3 ...  m1, m1+1, ...  m2, m2+1, ...  m3, m3+1 ...
c$$$c   ----atom1----  ----atom2----  ----atom2----  ----atom4----
c$$$c imdim(1) = 1
c$$$c imdim(2) = m1 + 1  etc.
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension mdim(nclass),iclass(natom)
c$$$      dimension imdim(natom)
c$$$
c$$$      imdim(1)   = 1
c$$$      do       i = 2,natom
c$$$cmod 93.06.07
c$$$c     ic         = iclass(i)
c$$$c     imdim(i)   = imdim(i-1) + mdim(i-1)
c$$$        ic         = iclass(i-1)
c$$$        imdim(i)   = imdim(i-1) + mdim(ic)
c$$$c-----
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      integer function indxpp (l1,n1,l2,n2,
c$$$     i                         iln,nl,nn,nln)
c$$$
c$$$c index for phi(l1,n1) phi(l2,n2), ln1 >= ln2 (lower triangle)
c$$$c indxpp = (l1,n1 ; l2,n2)
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension iln(0:nl-1,nn)
c$$$
c$$$      ln1        = iln(l1,n1)
c$$$      ln2        = iln(l2,n2)
c$$$      irow       = max0(ln1,ln2)
c$$$      icol       = ln1 + ln2 - irow
c$$$      indxpp     = lowtri(nln,irow,icol)
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine indxql (lmx,
c$$$     i                   il,in,nln,
c$$$     i                   nql,
c$$$     o iql,iql1,iqn1,iql2,iqn2,iqlln )
c$$$
c$$$c indexing for qls
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension il(nln),in(nln)
c$$$      dimension iql(nql),iql1(nql),iqn1(nql),iql2(nql),iqn2(nql),
c$$$     o          iqlln(0:lmx,nln,nln)
c$$$
c$$$      do       l = 0,lmx
c$$$        do     ln1 = 1,nln
c$$$          do     ln2 = 1,nln
c$$$            iqlln(l,ln1,ln2) = 0
c$$$          end do
c$$$        end do
c$$$      end do
c$$$
c$$$      i          = 0
c$$$      do       l = 0,lmx
c$$$        do     ln1 = 1,nln
c$$$          l1         = il(ln1)
c$$$          n1         = in(ln1)
c$$$          do     ln2 = ln1,nln
c$$$            l2         = il(ln2)
c$$$            n2         = in(ln2)
c$$$            if(l .gt. l1+l2) goto 1111
c$$$            i          = i + 1
c$$$            iql(i)     = l
c$$$            iql1(i)    = l1
c$$$            iqn1(i)    = n1
c$$$            iql2(i)    = l2
c$$$            iqn2(i)    = n2
c$$$            iqlln(l,ln1,ln2) = i
c$$$ 1111       continue
c$$$          end do
c$$$        end do
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c--------------------------------------------------------------------
c$$$      subroutine indxrk(n1,n2,n3,qp,nq,g,ng,qbas,
c$$$     o                  ginv,indxr,indxk)
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$c generates indexing of symmetry operators and k-point
c$$$c inside the irreducible BZ (IBZ).
c$$$c given k inside BZ indxr(ia,ib,ic) = symmetry operator
c$$$c                   indxk(ia,ib,ic) = k-point inside IBZ
c$$$c where ia=a*n1, ib=b*n2, ic=c*n3 and
c$$$c k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors
c$$$c input:
c$$$c n1,n2,n3 = number of divisions along G1,G2,G3
c$$$c qp       = k-points cartesian coordinates in the IBZ
c$$$c nq       = number of irreducible k-points
c$$$c g        = symmetry operators
c$$$c ng       = number of symmetry operators
c$$$c qbas     = reciprocal base vectors
c$$$c output:
c$$$c indxr, indxk
c$$$c remarks:
c$$$c all vectors are in units of 2*pi/a
c$$$
c$$$      dimension qp(3,nq),g(3,3,ng),qbas(3,3),ginv(3,3),
c$$$     .          indxr(-n1:n1,-n2:n2,-n3:n3),work(3*3),
c$$$     .          indxk(-n1:n1,-n2:n2,-n3:n3)
c$$$
c$$$c initialise indxr and indxk
c$$$c      call iinit(indxr,(2*n1+1)*(2*n2+1)*(2*n3+1))
c$$$c      call iinit(indxk,(2*n1+1)*(2*n2+1)*(2*n3+1))
c$$$      indxr=0
c$$$      indxk=0
c$$$c invert qbas
c$$$      call dinv33(qbas,0,ginv,det)
c$$$
c$$$c rotate each k-point and label it with the symmetry operator
c$$$      do       i = 1,ng
c$$$        do       k = 1,nq
c$$$          call dmv(3,3,ginv,3,qp(1,k),work)
c$$$          ia         = idnint(n1*work(1))
c$$$          ib         = idnint(n2*work(2))
c$$$          ic         = idnint(n3*work(3))
c$$$          indxr(ia,ib,ic) = i
c$$$          indxk(ia,ib,ic) = k
c$$$        end do
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c--------------------------------------------------------------------
c$$$c------------------------------------------------------------------
c$$$      integer function maxocc (ifev,ispin,ef,
c$$$     d                         nlmto,nqbz)
c$$$c kino removed common w
c$$$c 92.02.28
c$$$c finds the max. no. occupied states from eigenvalues at all k
c$$$c NOTE: assumes that the eigenvalues are stored in direct access file
c$$$
c$$$c ifev(2) = direct access unit file for eigenvalues
c$$$c ispin   = 1, paramagnetic
c$$$c           2, ferromagnetic
c$$$c ef      = fermi level
c$$$c nlmto   = no. states
c$$$c nqbz    = no. k-points
c$$$
c$$$      implicit real*8 (a-h,o-z)
c$$$      integer   ifev(2)
c$$$
c$$$c integer work array
c$$$c      integer w
c$$$c      common /w/ w(1)
c$$$      real*8,allocatable :: w_iekt(:)
c$$$
c$$$c maximum no. occupied states
c$$$c      call defdr   (iekt,nqbz*nlmto)
c$$$      allocate( w_iekt( nqbz*nlmto ) )
c$$$      noccx      = 0
c$$$      do      is = 1,ispin
c$$$c> read eigenvalues
c$$$        call rwdd    (ifev(is),
c$$$     d              nlmto,nqbz, 
c$$$     o              w_iekt )
c$$$        noccxt     = noccx1 (w_iekt,nqbz,nlmto,ef)
c$$$        if (noccxt .gt. noccx) noccx = noccxt
c$$$      end do
c$$$      maxocc     = noccx
c$$$
c$$$c      call rlse    (iekt)
c$$$      deallocate( w_iekt )
c$$$      return
c$$$      end
c$$$c------------------------------------------------------------------
c$$$c------------------------------------------------------------------
c$$$
c$$$!> counts the no. occupied/unoccupied states
c$$$
c$$$c-------------------------------------------------------------------
c$$$      subroutine ncorex (nindxc,
c$$$     d                   nl,nclass,
c$$$     o                   ncore,ncx)
c$$$
c$$$c 92.03.18
c$$$c counts the no. core states/class
c$$$
c$$$c nindxc(l,cl) = no. core states
c$$$
c$$$c ncore(cl) = no. core states
c$$$c ncx       = max. no. core states
c$$$
c$$$      implicit real*8 (a-h,o-z)
c$$$      dimension nindxc(0:nl-1,nclass)
c$$$      dimension ncore(nclass)
c$$$
c$$$      do      ic = 1,nclass
c$$$        ncore(ic)  = ivsum (nindxc(0,ic),nl)
c$$$      end do
c$$$      ncx        = maxim (ncore,nclass)
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine noccmax(ekt,nk,nt,ef,
c$$$     o                   noccx)
c$$$
c$$$c find the highest number of occupied states
c$$$
c$$$c ekt   = eigenvalues for all k-points and states
c$$$c nk    = number of k-points in the 1st BZ
c$$$c nt    = number of states
c$$$c ef    = fermi level
c$$$
c$$$c noccx = highest number of occupied states
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension ekt(nk,nt)
c$$$
c$$$c find the highest number of occupied states
c$$$      noccx      = 0
c$$$      do       k = 1,nk
c$$$        do      it = 1,nt
c$$$          if(ekt(k,it) .gt. ef)goto 1111
c$$$        end do
c$$$ 1111   if(it .gt. noccx) noccx = it
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c------------------------------------------------------------------
c$$$c------------------------------------------------------------------
c$$$      integer function nofblch(mdim,iclass,nclass,natom)
c$$$
c$$$c number of bloch optimal product basis functions
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension mdim(nclass),iclass(natom)
c$$$
c$$$      nofblch    = 0
c$$$      do      ia = 1,natom
c$$$        ic         = iclass(ia)
c$$$        nofblch    = nofblch + mdim(ic)
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$c--------------------------------------------------------------------
c$$$c-------------------------------------------------------------------
c$$$      integer function noflnmx(nindx,nl,nclass)
c$$$
c$$$c maximum number of l,n,m
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nindx(0:nl-1,nclass)
c$$$
c$$$      noflnmx   = 0
c$$$      do     ic = 1,nclass
c$$$        nlnm      = noflnm(nindx(0,ic),nl)
c$$$        if(nlnm .gt. noflnmx) noflnmx = nlnm
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
c$$$      subroutine nofpppp(nofln,
c$$$     o npppp,npppd,nppdd)
c$$$
c$$$c number of
c$$$c <phi(l1,n1) phi(l2,n2) phi(l3,n3) phi(l4,n4)>
c$$$c (l1,n1) =< (l2,n2) =< (l3,n3) =< (l4,n4)
c$$$
c$$$c <phi(l1,n1) phi(l2,n2) phi(l3,n3) phidot(l4,n4)>
c$$$c (l1,n1) =< (l2,n2) =< (l3,n3)
c$$$
c$$$c <phi(l1,n1) phi(l2,n2) phidot(l3,n3) phidot(l4,n4)>
c$$$c (l1,n1) =< (l2,n2)  ;  (l3,n3) =< (l4,n4)
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$      npppp      = 0
c$$$      do 10   i1 = 1, nofln
c$$$      do 10   i2 = i1,nofln
c$$$      do 10   i3 = i2,nofln
c$$$      do 10   i4 = i3,nofln
c$$$        npppp      = npppp + 1
c$$$   10 continue
c$$$
c$$$      npppd      = 0
c$$$      do 20   i1 = 1, nofln
c$$$      do 20   i2 = i1,nofln
c$$$      do 20   i3 = i2,nofln
c$$$      do 20   i4 = 1 ,nofln
c$$$        npppd      = npppd + 1
c$$$   20 continue
c$$$
c$$$      nppdd      = 0
c$$$      do 30   i1 = 1, nofln
c$$$      do 30   i2 = i1,nofln
c$$$      do 30   i3 = 1 ,nofln
c$$$      do 30   i4 = i3,nofln
c$$$        nppdd      = nppdd + 1
c$$$   30 continue
c$$$
c$$$      return
c$$$      end
c$$$c-----------------------------------------------------------------
c$$$      integer function nofyyyy(nl)
c$$$
c$$$c number of ylm(i1)*ylm(i2)*ylm(i3)*ylm(i4)
c$$$c i1 .le. i2 .le. i3 .le. i4
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$      nofyyyy    = 0
c$$$      ll         = nl*nl
c$$$      do 10   i1 = 1,ll
c$$$      do 10   i2 = i1,ll
c$$$      do 10   i3 = i2,ll
c$$$      do 10   i4 = i3,ll
c$$$        nofyyyy    = nofyyyy + 1
c$$$   10 continue
c$$$
c$$$      return
c$$$      end
c$$$c------------------------------------------------------------------
c$$$      subroutine nolna (nocc,nunocc,nindx,
c$$$     d                  nl,nn,nclass,
c$$$     o                  nlna)
c$$$
c$$$c gives the number of allowed product radial phi
c$$$c nocc(l,n)   = 0,1 ==> unoccupied, occupied
c$$$c nunocc(l,n) = 1,0 ==> unoccupied,occupied
c$$$c nalwln    = number of allowed phi(l1,n1) phi(l2,n2)
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
c$$$     i          nindx(nl,nclass)
c$$$      dimension nlna(nclass)
c$$$
c$$$      do     ic = 1,nclass
c$$$        nlna(ic)  = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$c-------------------------------------------------------------------
