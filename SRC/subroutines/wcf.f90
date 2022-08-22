subroutine wcf(nbloch, vcoul, zx0, ixc, zw)
  !- x0 (ixc==0), W-v (ixc==1) or epsi (ixc==2) calculation.------------------
  !i ixc=-1 is inverse of zx0 to zw May2005
  !i Vcou
  !i x0
  !o W-v
  implicit none
  integer(4) :: i,nbloch,ixc,ipl1,nev,nmx,lwork,info,itype
  complex(8), dimension(nbloch,nbloch) :: vcoul,zw,zx0
  real(8), allocatable :: eb(:),rwork(:), &
       reps(:,:),ceps(:,:), repsi(:,:),cepsi(:,:)  !, rx0,cx0, rvx0,cvx0
  complex(8),allocatable :: zvx0(:,:),oo(:,:),hh(:,:),work(:)

  !------------------------------------
  ! calculate v x0
  !      call mmulc   (dreal(vcoul),dimag(vcoul),nbloch,     ! cou
  !     i              dreal(zx0), dimag(zx0), nbloch,     ! x0
  !     i              nbloch,nbloch,nbloch,nbloch,
  !     o              rvx0,cvx0 )           ! cou \times x0
  !      allocate( zvx0(nbloch,nbloch),
  !     &          reps(nbloch,nbloch), ceps(nbloch,nbloch),
  !     &          repsi(nbloch,nbloch),cepsi(nbloch,nbloch))


  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      if(ixc==2) then
  !        zw=zx0
  !        return
  !      endif
  !      if(.false.) then
  !        LWORK = nbloch*nbloch
  !          allocate(hh(nbloch,nbloch),oo(nbloch,nbloch)
  !     &     ,eb(nbloch),WORK(nbloch*nbloch)
  !     &     ,RWORK(max(1, 3*nbloch-2)) )
  !          hh  = zx0 !dcmplx(reps,ceps)
  !          oo  = 0d0
  !          do ipl1=1,nbloch
  !            oo(ipl1,ipl1)=1d0
  !          enddo
  !        ITYPE =1
  ! cxml of lapack. See http://www.compaq.com/math/documentation/cxml/zhegv.3lapack.html
  !        call ZHEGV( ITYPE, "N", "U", nbloch, hh, nbloch, oo, nbloch, eb,
  !     &    WORK, LWORK, RWORK, INFO )
  !          write(6,*) 'info=',info
  !          do ipl1=1,5
  !            write(6,'(i4,d24.16)')ipl1, eb(ipl1)
  !          enddo
  !          do ipl1= nbloch-4,nbloch
  !            write(6,'(i4,d24.16)')ipl1, eb(ipl1)
  !          enddo
  !        stop ' eigen test end'
  !      endif
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
     !-- calculate W - v
     if(ixc==1) then
        call matm(zvx0,vcoul,zw,nbloch,nbloch,nbloch) ! zw = 1/e v
        zw = zw - vcoul
     else
        zw = zvx0                                     ! zw =1/e
     endif
     deallocate( zvx0)
  endif

  ! calculate e = 1 - v X0
  !      reps = - dreal(zvx0)
  !      ceps = - dimag(zvx0)
  !      do i=1, nbloch
  !        reps(i,i) = 1d0   + reps(i,i)
  !      enddo

  ! invert e
  !      call minvc90  (reps,ceps,
  !     d               nbloch,nbloch,
  !     o               repsi,cepsi )
  ! cccccccccccccccccccccccccccccccccccccccccccccccc
  !  normalization test
  !      repsi = 0d0
  !      cepsi = 0d0
  !      do i=1, nbloch
  !        repsi(i,i) = 1d0
  !      enddo
  ! cccccccccccccccccccccccccccccccccccccccccccccccc

  !      if(ixc==1) then
  !-- calculate W - v
  !        call matm(dcmplx(repsi,cepsi),vcoul,zw,nbloch,nbloch,nbloch)
  !        zw = zw - vcoul
  !      else
  !        zw = dcmplx(repsi,cepsi)
  !      endif
  !      deallocate( zvx0, reps, ceps, repsi,cepsi)
end subroutine wcf

!-------------------------------------------------------------------
!      subroutine minvc90(r,c,  !takao omitted work area from minvc90
!     d                  ldim,n,
!     o ri,ci )

! 91.11.29
! invert a complex matrix
! r,c   = real and imaginary parts of the matrix
! ldim  = leading dimension of matrix
! n     = dimension of the problem
! work,ipvt,w1,w2 = work arrays
! ri,ci = real and imaginary parts of the inverse matrix

!      implicit double precision (a-h,o-z)
!      dimension r(ldim,n),c(ldim,n),
!     w          work(n),ipvt(n),w1(n,n),w2(n,n)
!      dimension ri(n,n),ci(n,n)

! invert real part
!      call minv    (r,ldim,n,
!     w              work,ipvt,
!     o              w1)

! real part of inverse
!      call mmul    (c,ldim,w1,n,n,n,n,n,
!     o              w2 )
!      call mmul    (w2,n,c,ldim,n,n,n,n,
!     o              ri )
!      call madd    (ri,n,r,ldim,n,n,n,
!     o              ci )
!      call minv    (ci,n,n,
!     w              work,ipvt,
!     o              ri )

! imaginary part of inverse
!      call mmul    (ri,n,w2,n,n,n,n,n,
!     o              ci )
!      call cv      (-1.d0,ci,n*n,
!     o              ci )

!      return
!      end

!------------
subroutine chknbnb( n1b,n2b,nbnb,n1b0,n2b0,nbnb0,nqbz,nbnbx)
  implicit none
  integer(4) :: nband, niwtnwt, nqbz, nbnbx ,iq,ibib
  integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),nbnb(nqbz)
  integer(4) :: n1b0(nbnbx,nqbz),n2b0(nbnbx,nqbz),nbnb0(nqbz)
  do iq   = 1,nqbz
     if( nbnb(iq)/=nbnb0(iq) ) then
        print *, 'chknbnb: nbnb(iq)/=',nbnb(iq),nbnb0(iq)
        call rx( 'chknbnb: nbnb(iq)/=nbnb0(iq)')
     endif
     do ibib = 1,nbnb(iq)
        if( n1b(ibib,iq)/=n1b0(ibib,iq) .OR. &
             n2b(ibib,iq)/=n2b0(ibib,iq) ) then
           print *, 'chknbnb: n1b',ibib,iq,n1b(ibib,iq),n1b0(ibib,iq)
           print *, 'chknbnb: n2b',ibib,iq,n2b(ibib,iq),n2b0(ibib,iq)
           call rx( 'chknbnb: ')
        endif
     enddo
  enddo
end subroutine chknbnb
!---------------

subroutine anfx0k(natom,nclass,mdim,iclass,bas,nbloch,ngci,q,ngvecc,qlat,anfvec,iaf, zxq) 
  !- antiferro part is added to x0k
  ! We assume that the crystal has a magnetic symmetry described by (translataion + spin flip).
  ! The translation is specified by a vector,
  !   AFvector = anfvec(1:3)*alat, which is the true real vector in Cartesian coodinate.
  ! Each mixed basis is mapped to the other mixed basis.
  ! E.g. the product basis B({\bf r}-{\bf a}) is mapped to
  !  B({\bf r}-{\bf a}-{\bf A}) = B({\bf r}-{\bf a}'-{\bf T}_0),
  !  by the translation specified by AFvec={\bf A}.
  !  Here {\bf T}_0 is some crystal tralslation vector.
  !  In this code you see,
  !      bas(1:3,ia1)+ anfvec  = bas(1:3,iaf(ia1)) + transaf(1:3,ia1)
  !  ==   {\bf a}    +{\bf A}  = {\bf a}'          + {\bf T}_0
  ! ---- The corresponding atoms should have the same product basis.
  implicit none
  integer :: natom, nbloch, nclass,ngci,ngb
  integer :: mdim(nclass), iaf(natom), iof(natom),iclass(natom), &
       ix,ia1,ia2,ic1,ic2,ifi,i,j, iaf1,iaf2,igp,im1,im2
  integer(4):: iam(nbloch), ngvecc(1:3,ngci),imf(nbloch+ngci)
  real(8) :: qt(natom), rf,cf,q(3),qlat(3,3),bas(3,natom),anfvec(3),transaf(3,natom),qg(3)
  complex(8) :: zxq (nbloch+ngci,nbloch+ngci),fac(nbloch+ngci)
  complex(8):: imagtwopi ,imag=(0d0,1d0)
  complex(8),allocatable :: zxqw(:,:)
  real(8) :: pi=3.1415926535897932
  integer(4):: verbose
  !r    True_q(1:3)     = 2*pi/alat * q(1:3)
  !r  True G is given by
  !r    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvecc(1:3,igp)) ,igp=1,ngp
  !------------------------
  imagtwopi = 2d0*(0d0,3.1415926535897932)
  do ia1 = 1, natom
     transaf(1:3,ia1)= bas(1:3,ia1)+ anfvec - bas(1:3,iaf(ia1))
     if(verbose()>=200)write(6,"(' ia1 transaf=',i3, 3f13.5)") ia1,transaf(1:3,ia1)
  enddo
  iof(1) = 0
  do ia1 = 1, natom-1
     iof(ia1+1)= iof(ia1) + mdim( iclass(ia1) )
     iam(iof(ia1)+1:iof(ia1+1)) = ia1
  enddo
  iam(iof(natom)+1:nbloch) = natom
  write(6,*) (ia1, mdim( iclass(ia1) ),ia1 = 1, natom)
  if( nbloch /= iof(natom) +mdim(iclass(natom)) ) call rx( ' anfx0k: nbloch /= ...')
  ! phase shifts
  do ia1 = 1, natom
     qt(ia1) = 2d0*pi*sum(q*transaf(1:3,ia1))
     if(verbose()>=200) write( 6, "(i3,2x,f10.5)") iaf(ia1),  qt(ia1)
  enddo
  if(verbose()>=200) write(6,*) ' anfx0k:  nbloch=',nbloch
  do im1 = 1,nbloch
     ia1  = iam(im1)
     i    = im1 - iof(ia1)
     imf (im1) = i + iof(iaf(ia1))
     fac (im1) = exp( imag*qt(ia1) )
  enddo
  do igp = 1,ngci
     im1 = nbloch+igp
     imf (im1) = im1
     qg = q + matmul(qlat, ngvecc(1:3,igp))
     fac(im1) = exp( imagtwopi*sum(qg*anfvec(1:3)))
  enddo
  ngb=nbloch + ngci
  allocate(zxqw(ngb,ngb))
  zxqw = zxq
  do im1 = 1,ngb
     do im2 = 1,ngb
        zxq(im1, im2)= zxq (im1,im2) + zxqw(imf(im1),imf(im2))*fac(im1)*dconjg(fac(im2))
     enddo
  enddo
  deallocate(zxqw)
  write(6,*) ' anfx0k: end '
  return
end subroutine anfx0k

