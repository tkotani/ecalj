      integer function maxocc2 (nspin,ef, nband, qbz,nqbz)
      use m_readeigen
      implicit none
      integer(4):: nspin,nqbz,nband,noccx,is,iq,noccxt,noccx1
      real(8) :: qbz(3,nqbz),ef
      real(8),allocatable :: ekt(:,:)
c maximum no. occupied states
      allocate( ekt(nband, nqbz ) )
      noccx      = 0
      do  is = 1,nspin
        do iq = 1,nqbz
          call readeval(qbz(:,iq),is, ekt(1,iq))
        enddo
        noccxt  = noccx1 (ekt,nqbz,nband,ef)
        if (noccxt .gt. noccx) noccx = noccxt
      enddo
      maxocc2 = noccx
      deallocate( ekt )
      end
c------------------------------------------------------------------


      subroutine idxlmto (nindxv,iclass,
     d                    nl,nnv,nlmto,natom,nclass,
     o il,in,im,ilmto)

c 92.03.14
c 92.10.13 from idxlmto (see above) modified for all atoms
c indexing of LMTO basis functions 
c follows that in TB-LMTO program
c
c l    0     1                 2 
c n    1  2  1        2        1               2
c m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
c ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
c
c il(ind),in(ind),im(ind) = l,n,m
c ilmto(n,lm,atom) = index of lmto
c lm = l*l + l + m + 1

      implicit real*8(a-h,o-z)
      dimension nindxv(0:nl-1,nclass),iclass(natom),
     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto)

      ind       = 0
      do 1   ia = 1,natom
      ic        = iclass(ia)
      do 1    l = 0,nl-1
      l2        = l*l
      do 1    n = 1,nindxv(l,ic)
      do 1    m = -l,l
      ind       = ind + 1
      lm        = l2 + l + m + 1
      il(ind)   = l
      in(ind)   = n 
      im(ind)   = m 
    1 ilmto(n,lm,ia) = ind
      if (ind .ne. nlmto) stop 'idxlmto: wrong no. l,n,m'

      return
      end
c-------------------------------------------------------------------
      subroutine idxlmt1 (nindxv,iclass,
     d                    nl,nnv,nlmto,natom,nclass,
     o il,in,im,ilmto,ialmto)

c 93.08.?? from idxlmto (see above) modified to include atom index
c indexing of LMTO basis functions 
c follows that in TB-LMTO program
c
c l    0     1                 2 
c n    1  2  1        2        1               2
c m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
c ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
c
c il(ind),in(ind),im(ind) = l,n,m
c ilmto(n,lm,atom) = index of lmto
c ialmto(ind) = ia
c lm = l*l + l + m + 1

      implicit real*8(a-h,o-z)
      dimension nindxv(0:nl-1,nclass),iclass(natom),
     o          ilmto(nnv,nl*nl,natom),in(nlmto),il(nlmto),im(nlmto),
     o          ialmto(nlmto)

      ind       = 0
      do 1   ia = 1,natom
      ic        = iclass(ia)
      do 1    l = 0,nl-1
      l2        = l*l
      do 1    n = 1,nindxv(l,ic)
      do 1    m = -l,l
      ind       = ind + 1
      lm        = l2 + l + m + 1
      il(ind)   = l
      in(ind)   = n 
      im(ind)   = m 
      ialmto(ind) = ia
    1 ilmto(n,lm,ia) = ind
      if (ind .ne. nlmto) stop 'idxlmto: wrong no. l,n,m'

      return
      end
c-------------------------------------------------------------------
      subroutine idxlna (nindx,
     d                   nl,nn,nclass,nlnx,
     o                   ilr,inr,ilnr,nlnr)

c 92.jan.13
c index of l,n for radial phi(l,n) for all classes

c nlnx = max. number of l,n  from maxdim.f

c ilr(index) = l
c inr(index) = n
c ilnr(l,n)  = index
c nlnr       = no. l,n = max(index)

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass)    
      dimension ilr(nlnx,nclass),inr(nlnx,nclass),
     o          ilnr(0:nl-1,nn,nclass),nlnr(nclass)

      do      ic = 1,nclass
      i          = 0
      do       l = 0,nl-1
      do       n = 1,nindx(l,ic)
      i          = i + 1
      ilr(i,ic)   = l
      inr(i,ic)   = n
      ilnr(l,n,ic)= i
      end do
      end do
      if(i .gt. nlnx)stop 'idxlna: nln .gt. nlnx'
      nlnr(ic)    = i
      end do
      
      return
      end
c-------------------------------------------------------------------
      subroutine idxlnp (nocc,nunocc,nindx,
     i                   nl,nn,nclass,nlnax,
     o                   il1a,in1a,il2a,in2a,ilna,nlna)

c gives an index for the allowed product phi

c nocc(l,n) = 0,1 ==> unoccupied, occupied
c nlna      = number of allowed product phi,
c             from integer function nalwln

c in1(lbas) = n1 etc.
c ilna(l1,n1,l2,n2) = index of allowed product phi
c ilna(l1,n1,l2,n2) = 0 ==> not allowed
  
      implicit real*8(a-h,o-z)
      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
     i          nindx(nl,nclass)
      dimension il1a(nlnax,nclass),in1a(nlnax,nclass),
     o          il2a(nlnax,nclass),in2a(nlnax,nclass),
     o          ilna(nl*nn*nl*nn,nclass),nlna(nclass)

      do      ic = 1,nclass
      nlna(ic)   = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)

      call indxlna (nocc(1,ic),nunocc(1,ic),nindx(1,ic),
     d              nl,nn,nlna(ic),
     o              il1a(1,ic),in1a(1,ic),il2a(1,ic),in2a(1,ic),
     o              ilna(1,ic) )
      end do

      return
      end
c---------------------------------------------------------------
      subroutine incor  (ncwf,nindxc,iclass,
     d                   nl,nnc,nclass,natom,
     o                   icore,ncore,nctot)

c 92.03.18
c sorts out allowed core states and count the number of core states

c ncwf(l,n,cl) = 1 ==> allowed, 0 ==> not allowed
c nindxc(l,cl)  = no. core states/l,class
c nl,nnc = max. no. l,n

c icore(i,cl) = index for allowed core states
c ncore(cl)   = no. allowed core states
c nctot       = total no. allowed core states

      implicit real*8 (a-h,o-z)
      dimension ncwf(0:nl-1,nnc,nclass),nindxc(0:nl-1,nclass),
     i          iclass(natom)
      dimension icore(nl*nl*nnc,nclass),ncore(nclass)

      ncx        = nl*nl*nnc
      do      ic = 1,nclass
      i          = 0
      j          = 0
      do       l = 0,nl-1
      do       n = 1,nindxc(l,ic)
      do       m = -l,l
      j          = j + 1
      if (ncwf(l,n,ic) .eq. 1) then
      i          = i + 1
      if (i .gt. ncx) stop 'incore: wrong ncx'
      icore(i,ic)= j
      endif
      end do
      end do
      end do
      ncore(ic)  = i
      end do

c total no. allowed core states
      nctot      = 0
      do       i = 1,natom
      ic         = iclass(i)
      nctot      = nctot + ncore(ic)
      end do

      return
      end
c---------------------------------------------------------------
      subroutine indxbas(nocc,nunocc,nindx,nl,nn,
     i                   nallow,
     o                   il1,in1,im1,
     o                   il2,in2,im2,indxb)

c gives an index for the allowed product basis

c nocc(l,n) = 0,1 ==> unoccupied, occupied
c nallow    = number of allowed product basis,
c             from integer function nallow, 
c             must be equal to final lbas

c in1(lbas) = n1 etc.
c indxb(n1,lm1,n2,lm2) = index of allowed product basis
c indxb(n1,lm1,n2,lm2) = 0 ==> not allowed
  
      implicit real*8(a-h,o-z)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension indxb(nn,nl*nl,nn,nl*nl),
     o          in1(nallow),il1(nallow),im1(nallow),
     o          in2(nallow),il2(nallow),im2(nallow)

c write to  file PRODUCT if it exists
      ifprod     = ifile('PRODUCT')
      if (ifprod .gt. 0)
     .write (ifprod,*)'lbas,n1,l1,im1,n2,l2,im2'

c initialise indxb
      call iinit   (indxb,nn*nn*nl**4)

      lbas       = 0
      do 10   l1 = 0,nl-1
      ll1        = l1*l1
      do 10   n1 = 1,nindx(l1)
      do 10   m1 = 1,2*l1+1
      do 10   l2 = 0,nl-1
      ll2        = l2*l2
      do 10   n2 = 1,nindx(l2)
      do 10   m2 = 1,2*l2+1
c     indxb(n1,ll1+m1,n2,ll2+m2) = 0
      if (nocc(l1,n1) .eq. 0)   goto 10
      if (nunocc(l2,n2) .eq. 0) goto 10
      if (l1 .ne. l2 .or. n1 .ne. n2) then
      if (indxb(n2,ll2+m2,n1,ll1+m1) .ne. 0) goto 10
      endif

c temporary
      if (l1 .eq. l2 .and. n1 .eq. n2) then
      if (indxb(n2,ll2+m2,n1,ll1+m1) .ne. 0) goto 10
      endif

      lbas       = lbas + 1
      indxb(n1,ll1+m1,n2,ll2+m2) = lbas
      in1(lbas)  = n1
      il1(lbas)  = l1
      im1(lbas)  = m1 - l1 - 1
      in2(lbas)  = n2
      il2(lbas)  = l2
      im2(lbas)  = m2 - l2 - 1
      if (ifprod .gt. 0)
     .write (ifprod,*)lbas,n1,l1,im1(lbas),n2,l2,im2(lbas)
   10 continue
      if(lbas .ne. nallow)stop 'indxbas: lbas .ne. nallow'

      return
      end
c-------------------------------------------------------------------
      subroutine indxk(n1,n2,n3,qbas,
     i                 qbz,nqbz,
     o                 ginv,indk)

c generates indexing of k-points inside 1st BZ
c given k inside BZ indk(ia,ib,ic) = k-point inside 1st BZ
c where ia=a*n1, ib=b*n2, ic=c*n3 and
c k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors

c n1,n2,n3 = number of divisions along G1,G2,G3
c qbz      = k-points cartesian coordinates in the 1st BZ
c nqbz     = number of k-points
c qbas     = reciprocal base vectors

c indk     = s. above
c ginv     = inverse of base reciprocal vector
c notes:   all vectors are in units of 2*pi/a

      implicit real*8(a-h,o-z)
      dimension qbz(3,nqbz),qbas(3,3)
      dimension ginv(3,3),
     o          indk(-n1:n1,-n2:n2,-n3:n3)
      dimension work(3)

c write to unit ifkp
c     ifkp       = ifile('KPNT')
c     if(ifkp .gt. 0)then
c     write(ifkp,*)' k-points in the 1st BZ '
c     write(ifkp,*)' nG1,nG2,nG3, qx,qy,qz, index '
c     endif

c initialise indk
      call iinit(indk,(2*n1+1)*(2*n2+1)*(2*n3+1))

c invert qbas
      call dinv33(qbas,0,ginv,det)

      do       k = 1,nqbz
      call dmv(3,3,ginv,3,qbz(1,k),work)
      ia         = idnint(n1*work(1))
      ib         = idnint(n2*work(2))
      ic         = idnint(n3*work(3))
      indk(ia,ib,ic) = k
c     if (ifkp .gt. 0) 
c    .write(ifkp,6000)ia,ib,ic,qbz(1,k),qbz(2,k),qbz(3,k),
c    .indk(ia,ib,ic)
c6000 format(3i4,3f10.4,i5)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine indxlm(nl,
     o ila,ima)

c index of l,m
c i = (l+1)**2 + m - l

      implicit real*8(a-h,o-z)
      dimension ila(nl*nl),ima(nl*nl)

      i          = 0
      do       l = 0,nl-1
      do       m = -l,l
      i          = i + 1
      ila(i)     = l
      ima(i)     = m
      end do
      end do
      if(i .ne. nl*nl)stop 'indxlm: wrong no. of l,m'
      
      return
      end
c-------------------------------------------------------------------
      subroutine indxln(nindx,nl,
     i                  nn,nofln,
     o il,in,iln)

c index of l,n
c nofln = number of l,n  from integer function nofln

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)    
      dimension il(nofln),in(nofln),iln(0:nl-1,nn)

      i          = 0
      do       l = 0,nl-1
      do       n = 1,nindx(l)
      i          = i + 1
      il(i)      = l
      in(i)      = n
      iln(l,n)   = i
      end do
      end do
      if(i .ne. nofln)stop 'indxln: wrong no. of ln'
      
      return
      end
c---------------------------------------------------------------
      subroutine indxlna(nocc,nunocc,nindx,nl,nn,
     i                   nlna,
     o                   il1a,in1a,il2a,in2a,ilna)

c gives an index for the allowed product phi

c nocc(l,n) = 0,1 ==> unoccupied, occupied
c nlna      = number of allowed product phi,
c             from integer function nalwln

c in1(lbas) = n1 etc.
c ilna(l1,n1,l2,n2) = index of allowed product phi
c ilna(l1,n1,l2,n2) = 0 ==> not allowed
  
      implicit real*8(a-h,o-z)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension il1a(nlna),in1a(nlna),
     o          il2a(nlna),in2a(nlna),ilna(0:nl-1,nn,0:nl-1,nn)

      call iinit   (ilna,nl*nl*nn*nn)
      lbas       = 0
      do 10   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
      do 10   l2 = 0,nl-1
      do 10   n2 = 1,nindx(l2)
c     ilna(l1,n1,l2,n2) = 0
      if(nocc(l1,n1) .eq. 0)goto 10
      if(nunocc(l2,n2) .eq. 0)goto 10
      if((l1.ne.l2 .or. n1.ne.n2) .and. ilna(l2,n2,l1,n1) .ne. 0)goto 10
      lbas       = lbas + 1
      il1a(lbas) = l1
      in1a(lbas) = n1
      il2a(lbas) = l2
      in2a(lbas) = n2
      ilna(l1,n1,l2,n2) = lbas
   10 continue
      if(lbas .ne. nlna)stop 'indxlna: lbas .ne. nlna'

      return
      end
c--------------------------------------------------------------------
      subroutine indxlnm(nindx,nl,
     i                   nn,noflnm,
     o il,in,im,ilnm)

c indexing of LMTO basis functions for a given class, 
c follows that in TB-LMTO program
c il,in,im = l,n,m
c ilnm(n,lm) = index of n,l,m
c lm = l*l + l + m + 1

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)
      dimension ilnm(nn,nl*nl),in(noflnm),il(noflnm),im(noflnm)

      ind       = 0
      do 1    l = 0,nl-1
      l2        = l*l
      do 1    n = 1,nindx(l)
      do 1    m = 1,2*l+1
      ind       = ind + 1
      lm        = l2 + m
      il(ind)   = l
      in(ind)   = n
      im(ind)   = m - l - 1
    1 ilnm(n,lm) = ind
      if (ind .ne. noflnm) stop 'indxlnm: wrong no. l,n,m'

      return
      end
c-------------------------------------------------------------------
      subroutine idxlnma(nindx,nl,nclass,
     i                   nn,nlnmx,
     o il,in,im,ilnm)

c 92.jan.07
c indexing of LMTO basis functions for all classes, 
c follows that in TB-LMTO program
c il,in,im = l,n,m
c ilnm(n,lm) = index of n,l,m
c lm = l*l + l + m + 1

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass)
      dimension ilnm(nn,nl*nl,nclass),
     o          in(nlnmx,nclass),
     o          il(nlnmx,nclass),
     o          im(nlnmx,nclass)

      do     ic = 1,nclass
      ind       = 0
      do 1    l = 0,nl-1
      l2        = l*l
      do 1    n = 1,nindx(l,ic)
      do 1    m = 1,2*l+1
      ind       = ind + 1
      lm        = l2 + m
      il(ind,ic)= l
      in(ind,ic)= n
      im(ind,ic)= m - l - 1
    1 ilnm(n,lm,ic) = ind
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine idxlnmc(nindxv,nindxc,
     d                   nl,nn,nnv,nnc,nlnmx,nlnmxv,nlnmxc,nclass,
     o                   il,in,im,ilnm,
     o                   ilv,inv,imv,ilnmv,
     o                   ilc,inc,imc,ilnmc)

c 92.jan.07
c 92.03.17 include core states
c indexing of core states and LMTO basis functions for all classes, 
c follows that in TB-LMTO program
c il,in,im = l,n,m
c ilnm(n,lm) = index of n,l,m
c lm = l*l + l + m + 1
c NOTE: the indexing starts with core first and then valence on top
c       of core (not the same as index generated from nindx)

      implicit real*8(a-h,o-z)
      dimension nindxv(0:nl-1,nclass),nindxc(0:nl-1,nclass)
      dimension ilnm(nn,nl*nl,nclass),
     o          ilnmv(nnv,nl*nl,nclass),
     o          ilnmc(nnc,nl*nl,nclass),
     o          in(nlnmx,nclass),il(nlnmx,nclass),im(nlnmx,nclass),
     o         inv(nlnmxv,nclass),ilv(nlnmxv,nclass),imv(nlnmxv,nclass),
     o         inc(nlnmxc,nclass),ilc(nlnmxc,nclass),imc(nlnmxc,nclass)

      do     ic = 1,nclass
      ind       = 0

c core
      do      l = 0,nl-1
      l2        = l*l
      do      n = 1,nindxc(l,ic)
      do      m = 1,2*l+1
      ind       = ind + 1
      if (ind .gt. nlnmx) stop 'idxlnmc: ind > nlnmx'
      lm        = l2 + m
      il(ind,ic)= l
      in(ind,ic)= n
      im(ind,ic)= m - l - 1
      ilnm(n,lm,ic) = ind
      ilc(ind,ic)= l
      inc(ind,ic)= n
      imc(ind,ic)= m - l - 1
      ilnmc(n,lm,ic)= ind
      end do
      end do
      end do

c valence
      indv      = 0
      do      l = 0,nl-1
      l2        = l*l
      ncore     = nindxc(l,ic)
      do      n = 1,nindxv(l,ic)
      if (ncore+n .gt. nn) stop 'idxlnmc: ncore+n > nn'
      do      m = 1,2*l+1
      ind       = ind + 1
      indv      = indv + 1
      if (ind .gt. nlnmx) stop 'idxlnmc: ind > nlnmx'
      lm        = l2 + m
      il(ind,ic)= l
      in(ind,ic)= ncore + n
      im(ind,ic)= m - l - 1
      ilnm(ncore+n,lm,ic) = ind
      ilv(indv,ic)= l
      inv(indv,ic)= n
      imv(indv,ic)= m - l - 1
      ilnmv(n,lm,ic) = indv
      end do
      end do
      end do

      end do

      return
      end
c-------------------------------------------------------------------
      subroutine indxmdm (mdim,nclass,
     i                    iclass,natom,
     o imdim )

c given an atom, imdim points to the starting address of the
c total optimal product basis
c   1 2 3 ...  m1, m1+1, ...  m2, m2+1, ...  m3, m3+1 ...
c   ----atom1----  ----atom2----  ----atom2----  ----atom4----
c imdim(1) = 1
c imdim(2) = m1 + 1  etc.

      implicit real*8(a-h,o-z)
      dimension mdim(nclass),iclass(natom)
      dimension imdim(natom)

      imdim(1)   = 1
      do       i = 2,natom
cmod 93.06.07
c     ic         = iclass(i)
c     imdim(i)   = imdim(i-1) + mdim(i-1)
      ic         = iclass(i-1)
      imdim(i)   = imdim(i-1) + mdim(ic)
c-----
      end do

      return
      end
c-------------------------------------------------------------------
      integer function indxpp (l1,n1,l2,n2,
     i                         iln,nl,nn,nln)

c index for phi(l1,n1) phi(l2,n2), ln1 >= ln2 (lower triangle)
c indxpp = (l1,n1 ; l2,n2)

      implicit real*8(a-h,o-z)
      dimension iln(0:nl-1,nn)

      ln1        = iln(l1,n1)
      ln2        = iln(l2,n2)
      irow       = max0(ln1,ln2)
      icol       = ln1 + ln2 - irow
      indxpp     = lowtri(nln,irow,icol)

      return
      end
c-------------------------------------------------------------------
      subroutine indxql (lmx,
     i                   il,in,nln,
     i                   nql,
     o iql,iql1,iqn1,iql2,iqn2,iqlln )

c indexing for qls

      implicit real*8(a-h,o-z)
      dimension il(nln),in(nln)
      dimension iql(nql),iql1(nql),iqn1(nql),iql2(nql),iqn2(nql),
     o          iqlln(0:lmx,nln,nln)

      do       l = 0,lmx
      do     ln1 = 1,nln
      do     ln2 = 1,nln
      iqlln(l,ln1,ln2) = 0
      end do
      end do
      end do

      i          = 0
      do       l = 0,lmx
      do     ln1 = 1,nln
      l1         = il(ln1)
      n1         = in(ln1)
      do     ln2 = ln1,nln
      l2         = il(ln2)
      n2         = in(ln2)
      if(l .gt. l1+l2) goto 1111
      i          = i + 1
      iql(i)     = l
      iql1(i)    = l1
      iqn1(i)    = n1
      iql2(i)    = l2
      iqn2(i)    = n2
      iqlln(l,ln1,ln2) = i
 1111 continue
      end do
      end do
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine indxrk(n1,n2,n3,qp,nq,g,ng,qbas,
     o                  ginv,indxr,indxk)
      implicit real*8(a-h,o-z)

c generates indexing of symmetry operators and k-point
c inside the irreducible BZ (IBZ).
c given k inside BZ indxr(ia,ib,ic) = symmetry operator
c                   indxk(ia,ib,ic) = k-point inside IBZ
c where ia=a*n1, ib=b*n2, ic=c*n3 and
c k = a*G1 + b*G2 + c*G3, G1,G2,G3 are the reciprocal base vectors
c input:
c n1,n2,n3 = number of divisions along G1,G2,G3
c qp       = k-points cartesian coordinates in the IBZ
c nq       = number of irreducible k-points
c g        = symmetry operators
c ng       = number of symmetry operators
c qbas     = reciprocal base vectors
c output:
c indxr, indxk
c remarks:
c all vectors are in units of 2*pi/a

      dimension qp(3,nq),g(3,3,ng),qbas(3,3),ginv(3,3),
     .          indxr(-n1:n1,-n2:n2,-n3:n3),work(3*3),
     .          indxk(-n1:n1,-n2:n2,-n3:n3)

c initialise indxr and indxk
      call iinit(indxr,(2*n1+1)*(2*n2+1)*(2*n3+1))
      call iinit(indxk,(2*n1+1)*(2*n2+1)*(2*n3+1))

c invert qbas
      call dinv33(qbas,0,ginv,det)

c rotate each k-point and label it with the symmetry operator
      do       i = 1,ng
      do       k = 1,nq
      call dmv(3,3,ginv,3,qp(1,k),work)
      ia         = idnint(n1*work(1))
      ib         = idnint(n2*work(2))
      ic         = idnint(n3*work(3))
      indxr(ia,ib,ic) = i
      indxk(ia,ib,ic) = k
      end do
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine maxdim (nocc,nunocc,nindx,nl,nn,nclass,
     o nprodx,nlnx,nlnmx,nlnax)

c largest number of product basis, (l,n) and (l,n,m)

      implicit real*8(a-h,o-z)
      dimension nocc(nl*nn,nclass),
     i          nunocc(nl*nn,nclass),
     i          nindx(0:nl-1,nclass)

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
c------------------------------------------------------------------
      integer function maxocc (ifev,ispin,ef,
     d                         nlmto,nqbz)
c kino removed common w
c 92.02.28
c finds the max. no. occupied states from eigenvalues at all k
c NOTE: assumes that the eigenvalues are stored in direct access file

c ifev(2) = direct access unit file for eigenvalues
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level 
c nlmto   = no. states
c nqbz    = no. k-points

      implicit real*8 (a-h,o-z)
      integer   ifev(2)

c integer work array
c      integer w
c      common /w/ w(1)
      real*8,allocatable :: w_iekt(:)

c maximum no. occupied states
c      call defdr   (iekt,nqbz*nlmto) 
      allocate( w_iekt( nqbz*nlmto ) )
      noccx      = 0
      do      is = 1,ispin
c> read eigenvalues
      call rwdd    (ifev(is),
     d              nlmto,nqbz, 
     o              w_iekt )
      noccxt     = noccx1 (w_iekt,nqbz,nlmto,ef)
      if (noccxt .gt. noccx) noccx = noccxt
      end do
      maxocc     = noccx

c      call rlse    (iekt)
      deallocate( w_iekt )
      return
      end
c------------------------------------------------------------------
      integer function nallow (nocc,nunocc,nindx,nl,nn)

c gives the number of allowed product basis
c nocc(n,l) = 0,1 ==> unoccupied, occupied
c nallow    = number of allowed product basis
  
      implicit real*8(a-h,o-z)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)

      if(nl-1 .gt. lmax) stop 'nallow: increase lmax'
      if(nn .gt. nnx) stop 'nallow: increase nnx'
      call iinit   (icheck,nnx*nnx*(lmax+1)*(lmax+1))
      do      l1 = 0,nl-1
      do      n1 = 1,nindx(l1)
      do      l2 = 0,nl-1
      do      n2 = 1,nindx(l2)
      icheck(l1,n1,l2,n2) = nocc(l1,n1)*nunocc(l2,n2)
      if (l1 .ne. l2 .or. n1 .ne. n2) then
      if (icheck(l1,n1,l2,n2)*icheck(l2,n2,l1,n1) .ne. 0)
     . icheck(l1,n1,l2,n2) = 0
      endif
      end do
      end do
      end do
      end do

      nallow     = 0
      do 10   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
      do 10   m1 = 1,2*l1+1

      do 10   l2 = 0,nl-1
      do 10   n2 = 1,nindx(l2)
      do 10   m2 = 1,2*l2+1

c     if (nocc(l1,n1) .eq. 0)goto 10
c     if (nunocc(l2,n2) .eq. 0)goto 10
      if (icheck(l1,n1,l2,n2) .eq. 0) goto 10

c temporary
      if (l1 .eq. l2 .and. n1.eq.n2 .and. m1.lt.m2)goto 10

      nallow     = nallow + 1
   10 continue

      return
      end
c------------------------------------------------------------------
      integer function nalwln (nocc,nunocc,nindx,nl,nn)

c gives the number of allowed product radial phi
c nocc(l,n)   = 0,1 ==> unoccupied, occupied
c nunocc(l,n) = 1,0 ==> unoccupied,occupied
c nalwln    = number of allowed phi(l1,n1) phi(l2,n2)
  
      implicit real*8(a-h,o-z)
      parameter (lmax=6,nnx=10)
      dimension nocc(0:nl-1,nn),nunocc(0:nl-1,nn),
     i          nindx(0:nl-1)
      dimension icheck(0:lmax,nnx,0:lmax,nnx)

      if (nl-1 .gt. lmax) stop 'nalwln: increase lmax'
      if (nn .gt. nnx) stop 'nalwln: increase nnx'
      call iinit   (icheck,nnx*nnx*(lmax+1)*(lmax+1))

      nalwln     = 0
      do 10   l1 = 0,nl-1
      do 10   n1 = 1,nindx(l1)
      if(nocc(l1,n1) .eq. 0)goto 10
      do 20   l2 = 0,nl-1
      do 20   n2 = 1,nindx(l2)
      if(nunocc(l2,n2) .eq. 0)goto 20
      if((l1.ne.l2 .or. n1.ne.n2) .and. icheck(l2,n2,l1,n1).ne.0)
     . goto 20
      nalwln     = nalwln + 1
      icheck(l1,n1,l2,n2) = nalwln
   20 continue
   10 continue

      return
      end
c-------------------------------------------------------------------
      integer function nocc (ek,ef,locc,
     d                       n)

c 92.02.18
c counts the no. occupied/unoccupied states

c ek   = eigenvalues
c ef   = fermi level
c locc = .true.  ==> occupied
c      = .false. ==> unoccupied
c n    = no. states

      implicit real*8 (a-h,o-z)
      logical   locc
      dimension ek(n)

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
      end
c-------------------------------------------------------------------
      subroutine ncorex (nindxc,
     d                   nl,nclass,
     o                   ncore,ncx)

c 92.03.18
c counts the no. core states/class

c nindxc(l,cl) = no. core states

c ncore(cl) = no. core states
c ncx       = max. no. core states

      implicit real*8 (a-h,o-z)
      dimension nindxc(0:nl-1,nclass)
      dimension ncore(nclass)

      do      ic = 1,nclass
      ncore(ic)  = ivsum (nindxc(0,ic),nl)
      end do
      ncx        = maxim (ncore,nclass)

      return
      end 
c-------------------------------------------------------------------
      subroutine noccmax(ekt,nk,nt,ef,
     o                   noccx)

c find the highest number of occupied states

c ekt   = eigenvalues for all k-points and states
c nk    = number of k-points in the 1st BZ
c nt    = number of states
c ef    = fermi level

c noccx = highest number of occupied states

      implicit real*8(a-h,o-z)
      dimension ekt(nk,nt)

c find the highest number of occupied states
      noccx      = 0
      do       k = 1,nk
      do      it = 1,nt
      if(ekt(k,it) .gt. ef)goto 1111
      end do
 1111 if(it .gt. noccx) noccx = it
      end do

      return
      end
c------------------------------------------------------------------
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

c find the highest number of occupied states
      noccx      = 0
      do       k = 1,nk
      do      it = 1,nt
      if(ekt(it,k) .gt. ef)goto 1111
      end do
 1111 if(it .gt. noccx) noccx = it
      end do
ctm
      if (noccx.gt.nt) noccx = nt

      noccx1     = noccx
      return
      end
c------------------------------------------------------------------
      integer function nofblch(mdim,iclass,nclass,natom)

c number of bloch optimal product basis functions

      implicit real*8(a-h,o-z)
      dimension mdim(nclass),iclass(natom)

      nofblch    = 0
      do      ia = 1,natom
      ic         = iclass(ia)
      nofblch    = nofblch + mdim(ic)
      end do

      return
      end 
c-------------------------------------------------------------------
      integer function noflmto(nindx,iclass,nl,nclass,natom)

c total number of LMTO basis functions

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass),iclass(natom)

      noflmto   = 0
      do 1    i = 1,natom 
      ic        = iclass(i)
      do 1    l = 0,nl-1
      noflmto   = noflmto + (2*l+1)*nindx(l,ic)
    1 continue

      return
      end
c--------------------------------------------------------------------
      integer function nofln(nindx,nl)

c count the number of l,n

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)    

      nofln      = 0
      do       l = 0,nl-1
      nofln      = nofln + nindx(l) 
      end do
 
      return
      end
c------------------------------------------------------------------
      integer function noflnm(nindx,nl)

c number of l,n,m

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1)

      noflnm    = 0
      do 1    l = 0,nl-1
      noflnm    = noflnm + nindx(l)*(2*l+1)
    1 continue

      return
      end
c-------------------------------------------------------------------
      integer function noflnmx(nindx,nl,nclass)

c maximum number of l,n,m

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass)

      noflnmx   = 0
      do     ic = 1,nclass
      nlnm      = noflnm(nindx(0,ic),nl)
      if(nlnm .gt. noflnmx) noflnmx = nlnm
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine nofpppp(nofln,
     o npppp,npppd,nppdd)

c number of
c <phi(l1,n1) phi(l2,n2) phi(l3,n3) phi(l4,n4)>
c (l1,n1) =< (l2,n2) =< (l3,n3) =< (l4,n4)

c <phi(l1,n1) phi(l2,n2) phi(l3,n3) phidot(l4,n4)>
c (l1,n1) =< (l2,n2) =< (l3,n3)

c <phi(l1,n1) phi(l2,n2) phidot(l3,n3) phidot(l4,n4)>
c (l1,n1) =< (l2,n2)  ;  (l3,n3) =< (l4,n4)

      implicit real*8(a-h,o-z)
      
      npppp      = 0
      do 10   i1 = 1, nofln
      do 10   i2 = i1,nofln
      do 10   i3 = i2,nofln
      do 10   i4 = i3,nofln
      npppp      = npppp + 1
   10 continue 

      npppd      = 0
      do 20   i1 = 1, nofln
      do 20   i2 = i1,nofln
      do 20   i3 = i2,nofln
      do 20   i4 = 1 ,nofln
      npppd      = npppd + 1
   20 continue 

      nppdd      = 0
      do 30   i1 = 1, nofln
      do 30   i2 = i1,nofln
      do 30   i3 = 1 ,nofln
      do 30   i4 = i3,nofln
      nppdd      = nppdd + 1
   30 continue 

      return
      end
c-----------------------------------------------------------------
      integer function nofyyyy(nl)

c number of ylm(i1)*ylm(i2)*ylm(i3)*ylm(i4)
c i1 .le. i2 .le. i3 .le. i4

      implicit real*8(a-h,o-z)

      nofyyyy    = 0
      ll         = nl*nl
      do 10   i1 = 1,ll
      do 10   i2 = i1,ll
      do 10   i3 = i2,ll
      do 10   i4 = i3,ll
      nofyyyy    = nofyyyy + 1
   10 continue
    
      return
      end            
c------------------------------------------------------------------
      subroutine nolna (nocc,nunocc,nindx,
     d                  nl,nn,nclass,
     o                  nlna)

c gives the number of allowed product radial phi
c nocc(l,n)   = 0,1 ==> unoccupied, occupied
c nunocc(l,n) = 1,0 ==> unoccupied,occupied
c nalwln    = number of allowed phi(l1,n1) phi(l2,n2)
  
      implicit real*8(a-h,o-z)
      dimension nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
     i          nindx(nl,nclass)
      dimension nlna(nclass)

      do     ic = 1,nclass
      nlna(ic)  = nalwln (nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn)
      end do
      
      return
      end
c-------------------------------------------------------------------
      subroutine nolnma  (nindx,nl,nclass,
     o nlnm )

c 92.jan.07
c number of l,n,m for all classes

      implicit real*8(a-h,o-z)
      dimension nindx(0:nl-1,nclass)
      dimension nlnm(nclass)

      do     ic = 1,nclass
      noflnm    = 0
      do 1    l = 0,nl-1
      noflnm    = noflnm + nindx(l,ic)*(2*l+1)
    1 continue
      nlnm(ic)  = noflnm
      end do

      return
      end
