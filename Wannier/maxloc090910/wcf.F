      subroutine wcf(nbloch, vcoul, zx0, ixc, zw) 
C- x0 (ixc==0), W-v (ixc==1) or epsi (ixc==2) calculation ------------------
Ci ixc=-1 is inverse of zx0 to zw May2005
ci Vcou 
ci Vcou 
ci x0
co W-v  
      implicit none
      integer(4) :: i,nbloch,ixc,ipl1,nev,nmx,lwork,info,itype
      complex(8), dimension(nbloch,nbloch) :: vcoul,zw,zx0
      real(8), allocatable :: eb(:),rwork(:), 
     &  reps(:,:),ceps(:,:), repsi(:,:),cepsi(:,:)  !, rx0,cx0, rvx0,cvx0 
      complex(8),allocatable :: zvx0(:,:),oo(:,:),hh(:,:),work(:)

c------------------------------------
c calculate v x0
c      call mmulc   (dreal(vcoul),dimag(vcoul),nbloch,     ! cou
c     i              dreal(zx0), dimag(zx0), nbloch,     ! x0 
c     i              nbloch,nbloch,nbloch,nbloch,
c     o              rvx0,cvx0 )           ! cou \times x0
c      allocate( zvx0(nbloch,nbloch),
c     &          reps(nbloch,nbloch), ceps(nbloch,nbloch),
c     &          repsi(nbloch,nbloch),cepsi(nbloch,nbloch))
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if(ixc==2) then
c        zw=zx0
c        return
c      endif
c      if(.false.) then
c        LWORK = nbloch*nbloch
c          allocate(hh(nbloch,nbloch),oo(nbloch,nbloch)
c     &     ,eb(nbloch),WORK(nbloch*nbloch)
c     &     ,RWORK(max(1, 3*nbloch-2)) )
c          hh  = zx0 !dcmplx(reps,ceps)
c          oo  = 0d0
c          do ipl1=1,nbloch
c            oo(ipl1,ipl1)=1d0
c          enddo
c        ITYPE =1
c cxml of lapack. See http://www.compaq.com/math/documentation/cxml/zhegv.3lapack.html
c        call ZHEGV( ITYPE, "N", "U", nbloch, hh, nbloch, oo, nbloch, eb, 
c     &    WORK, LWORK, RWORK, INFO )
c          write(6,*) 'info=',info
c          do ipl1=1,5
c            write(6,'(i4,d24.16)')ipl1, eb(ipl1)
c          enddo
c          do ipl1= nbloch-4,nbloch
c            write(6,'(i4,d24.16)')ipl1, eb(ipl1)
c          enddo
c        stop ' eigen test end'
c      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( ixc==-1) then
        allocate( zvx0(nbloch,nbloch) )
        zvx0 = zx0
        call matcinv(nbloch,zvx0)     !  inverse
        zw = zvx0
        deallocate(zvx0)
      elseif(ixc==0) then
        zw = zx0
      else
        allocate( zvx0(nbloch,nbloch) )
        call matm( -vcoul,zx0, zvx0,nbloch,nbloch,nbloch) !  zvx0 <---  -v X0
        do i=1, nbloch
          zvx0(i,i) = 1d0 + zvx0(i,i) !  zvx0 <--- e = 1 - v X0
        enddo
        call matcinv(nbloch,zvx0)     !  zvx0 <--- 1/e

c-- calculate W - v
        if(ixc==1) then
          call matm(zvx0,vcoul,zw,nbloch,nbloch,nbloch) ! zw = 1/e v
          zw = zw - vcoul 
        else
          zw = zvx0                                     ! zw =1/e
        endif
        deallocate( zvx0)
      endif

c calculate e = 1 - v X0
c      reps = - dreal(zvx0)
c      ceps = - dimag(zvx0)
c      do i=1, nbloch
c        reps(i,i) = 1d0   + reps(i,i)
c      enddo
c
c invert e
c      call minvc90  (reps,ceps,
c     d               nbloch,nbloch,
c     o               repsi,cepsi )
cccccccccccccccccccccccccccccccccccccccccccccccccc
c  normalization test
c      repsi = 0d0
c      cepsi = 0d0
c      do i=1, nbloch
c        repsi(i,i) = 1d0  
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      if(ixc==1) then
c-- calculate W - v
c        call matm(dcmplx(repsi,cepsi),vcoul,zw,nbloch,nbloch,nbloch) 
c        zw = zw - vcoul 
c      else
c        zw = dcmplx(repsi,cepsi)
c      endif
c      deallocate( zvx0, reps, ceps, repsi,cepsi)
      end

c-------------------------------------------------------------------
      subroutine wcf_maxloc(nbloch, vcoul, zx0, zw_full, ixc, zw) 
C- x0 (ixc==0), W-v (ixc==1) or epsi (ixc==2) calculation ------------------
ci Vcou 
ci x0
ci Wfull - v
co U  
      implicit none
      integer(4) :: i,nbloch,ixc,ipl1,nev,nmx,lwork,info,itype
      complex(8), dimension(nbloch,nbloch) :: vcoul,zw,zx0,zw_full
      real(8), allocatable :: eb(:),rwork(:), 
     &  reps(:,:),ceps(:,:), repsi(:,:),cepsi(:,:)  !, rx0,cx0, rvx0,cvx0 
      complex(8),allocatable :: zvx0(:,:),oo(:,:),hh(:,:),work(:,:)

c------------------------------------
      if(ixc==0) then
        zw = zx0
      else
        allocate( zvx0(nbloch,nbloch),work(nbloch,nbloch) )
        work = vcoul + zw_full ! work = W <--- v + Wc
        call matm( work,zx0, zvx0,nbloch,nbloch,nbloch) !  zvx0 <---  W X0
        do i=1, nbloch
          zvx0(i,i) = 1d0 + zvx0(i,i) !  zvx0 <--- z = 1 + W X0
        enddo
        call matcinv(nbloch,zvx0)     !  zvx0 <--- 1/z

c-- calculate W
        if(ixc==1) then
          call matm(zvx0,work,zw,nbloch,nbloch,nbloch) ! zw = 1/z W
        else
          zw = zvx0                                    ! zw =1/z
        endif
        deallocate( zvx0,work)
      endif

      end

c-------------------------------------------------------------------
c      subroutine minvc90(r,c,  !takao omitted work area from minvc90
c     d                  ldim,n,
c     o ri,ci )
c
c 91.11.29
c invert a complex matrix
c r,c   = real and imaginary parts of the matrix
c ldim  = leading dimension of matrix
c n     = dimension of the problem
c work,ipvt,w1,w2 = work arrays 
c ri,ci = real and imaginary parts of the inverse matrix
c
c      implicit double precision (a-h,o-z)
c      dimension r(ldim,n),c(ldim,n),
c     w          work(n),ipvt(n),w1(n,n),w2(n,n)
c      dimension ri(n,n),ci(n,n)
c
c invert real part
c      call minv    (r,ldim,n,
c     w              work,ipvt,
c     o              w1)
c
c real part of inverse
c      call mmul    (c,ldim,w1,n,n,n,n,n,
c     o              w2 )
c      call mmul    (w2,n,c,ldim,n,n,n,n,
c     o              ri )
c      call madd    (ri,n,r,ldim,n,n,n,
c     o              ci )
c      call minv    (ci,n,n,
c     w              work,ipvt,
c     o              ri )
c
c imaginary part of inverse
c      call mmul    (ri,n,w2,n,n,n,n,n,
c     o              ci )
c      call cv      (-1.d0,ci,n*n,
c     o              ci )
c
c      return
c      end

c------------
      subroutine chknbnb( n1b,n2b,   nbnb,
     &                    n1b0,n2b0, nbnb0,
     &                    nqbz,nbnbx)
      implicit none 
      integer(4) :: nband, niwtnwt, nqbz, nbnbx ,iq,ibib
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),nbnb(nqbz)
      integer(4) :: n1b0(nbnbx,nqbz),n2b0(nbnbx,nqbz),nbnb0(nqbz)
      do iq   = 1,nqbz
        if( nbnb(iq)/=nbnb0(iq) ) then
           print *, 'chknbnb: nbnb(iq)/=',nbnb(iq),nbnb0(iq)
           stop 'chknbnb: nbnb(iq)/=nbnb0(iq)'
        endif
      do ibib = 1,nbnb(iq)
        if( n1b(ibib,iq)/=n1b0(ibib,iq) .or.
     &      n2b(ibib,iq)/=n2b0(ibib,iq) ) then
           print *, 'chknbnb: n1b',ibib,iq,n1b(ibib,iq),n1b0(ibib,iq)
           print *, 'chknbnb: n2b',ibib,iq,n2b(ibib,iq),n2b0(ibib,iq)
           stop 'chknbnb: '
        endif
      enddo
      enddo
      end
c---------------


c-------------------------------
c
c
c      if(iread) then
c        ifi = iopen ('ANFcond',1,0,0)
c        read ( ifi, *)                !skip lines
c        read ( ifi, *)
c        read ( ifi, *)
c        write( 6, "(' anti ferro condition from ANFcond')")
c        do ia1 = 1, natom
c          read ( ifi, *,      err=999) ix, iaf(ia1),(trans(i,ia1),i=1,3)
c          write( 6, "(2i3,2x,3f10.5)") ix, iaf(ia1),(trans(i,ia1),i=1,3)
c          if(ix.ne.ia1) goto 999
c        enddo
c        ifi = iclose('ANFcond')
c        iread=.false.
c      endif
c----------------------------------------------------------------------

      subroutine anfx0k(natom,nclass,mdim,iclass, bas,nbloch,ngci,
     i  q,ngvecc,qlat,    ! for q+G
     i  anfvec,iaf,  ! these are antiferro informations.
     i  zxq)  ! i/o 
C- antiferro part is added to x0k
c We assume that the crystal has a magnetic symmetry described by (translataion + spin flip).
c The translation is specified by a vector, 
c   AFvector = anfvec(1:3)*alat, which is the true real vector in Cartesian coodinate.
c   The anfvec(1:3) is given in ANFcond file.
c
c Each mixed basis is mapped to the other mixed basis.
c E.g. the product basis B({\bf r}-{\bf a}) is mapped to 
c  B({\bf r}-{\bf a}-{\bf A}) = B({\bf r}-{\bf a}'-{\bf T}_0),
c  by the translation specified by AFvec={\bf A}.
c  Here {\bf T}_0 is some crystal tralslation vector.
c  In this code you see,  
c      bas(1:3,ia1)+ anfvec  = bas(1:3,iaf(ia1)) + transaf(1:3,ia1)  
c  ==   {\bf a}    +{\bf A}  = {\bf a}'          + {\bf T}_0
c
c ---- The corresponding atoms should have the same product basis.
c
      implicit none
      integer natom, nbloch, nclass,ngci,ngb
      integer mdim(nclass), iaf(natom), iof(natom),iclass(natom),
     &  ix,ia1,ia2,ic1,ic2,ifi,i,j, iaf1,iaf2,igp,im1,im2
      integer(4):: iam(nbloch), ngvecc(1:3,ngci),imf(nbloch+ngci)
      real(8) :: qt(natom), rf,cf,q(3),qlat(3,3),bas(3,natom),
     &         anfvec(3),transaf(3,natom),qg(3)
      complex(8) :: zxq (nbloch+ngci,nbloch+ngci),fac(nbloch+ngci)
      complex(8):: imagtwopi ,imag=(0d0,1d0)
      complex(8),allocatable :: zxqw(:,:)
      real(8) :: pi=3.1415926535897932
      integer(4):: verbose
cr
cr    True_q(1:3)     = 2*pi/alat * q(1:3)
cr  True G is given by 
cr    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvecc(1:3,igp)) ,igp=1,ngp
c------------------------
      imagtwopi = 2d0*(0d0,3.1415926535897932)
c
      do ia1 = 1, natom
        transaf(1:3,ia1)= bas(1:3,ia1)+ anfvec - bas(1:3,iaf(ia1))
        if(verbose()>=200)
     &  write(6,"(' ia1 transaf=',i3, 3f13.5)") ia1,transaf(1:3,ia1)
      enddo
c
      iof(1) = 0
      do ia1 = 1, natom-1
        iof(ia1+1)= iof(ia1) + mdim( iclass(ia1) )
        iam(iof(ia1)+1:iof(ia1+1)) = ia1
      enddo
      iam(iof(natom)+1:nbloch) = natom
      write(6,*) (ia1, mdim( iclass(ia1) ),ia1 = 1, natom)
      if( nbloch .ne. iof(natom) +mdim(iclass(natom)) )
     &  stop ' anfx0k: nbloch.ne....'
c phase shifts 
      do ia1 = 1, natom
        qt(ia1) = 2d0*pi*sum(q*transaf(1:3,ia1))
        if(verbose()>=200)
     &          write( 6, "(i3,2x,f10.5)") iaf(ia1),  qt(ia1)
      enddo
c
        if(verbose()>=200)
     &        write(6,*) ' anfx0k:  nbloch=',nbloch
c
      do im1 = 1,nbloch
        ia1  = iam(im1)
        i    = im1 - iof(ia1)
        imf (im1) = i + iof(iaf(ia1))
        fac (im1) = exp( imag*qt(ia1) )
cccccccccccccccccccccccc
c        write( 669, "(5i4,2d13.5)") ia1,iaf(ia1), i, im1, imf(im1)
c     &  ,fac(im1)
cccccccccccccccccccccccc
      enddo

      do igp = 1,ngci
        im1 = nbloch+igp
        imf (im1) = im1
        qg = q + matmul(qlat, ngvecc(1:3,igp))
        fac(im1) = exp( imagtwopi*sum(qg*anfvec(1:3)))
      enddo
c
c      write(6,*) ' anfx0k: 1'
c
      ngb=nbloch + ngci
      allocate(zxqw(ngb,ngb))
      zxqw = zxq
c
c      write(6,*) ' anfx0k: 2'
c
      do im1 = 1,ngb
      do im2 = 1,ngb
        zxq(im1, im2)= zxq (im1,im2) 
     &               + zxqw(imf(im1),imf(im2))*fac(im1)*dconjg(fac(im2))
      enddo
      enddo
      deallocate(zxqw)
      write(6,*) ' anfx0k: end '
      return
      end

