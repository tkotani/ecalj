module m_rpaq
  contains
subroutine calcpv2 (zw,rv,cv, nbloch, rpv,cpv,rw,cw )
  !- calculate v^{1/2}Pv^{1/2} for given P and v --------
  !i zw : P, the polarization function. <\tilde{M}_i | P | \tilde{M}_j >
  !i rv+ i cv: v, the Coulomb matrix    <M_i |v |M_j>
  !i nbloch  = dimension of vc,P
  !o rw  + i cw: (1-Pv)(1-vP)
  !o rpv +i cpv:  Pv
  !--- History ---------
  ! tkotani modifyed from Miayke's calcpv Nov2004
  ! m 01.07.24
  ! pv is (1-Pv)(1-vP), not v^{1/2}Pv^{1/2}

  ! m 01.07.01
  ! 00.07.03
  !---------------------------------------------------------
  implicit double precision (a-h,o-z)
  implicit integer(i-n)
  dimension   rv(nbloch,nbloch),cv(nbloch,nbloch), &
       rpv(nbloch,nbloch),cpv(nbloch,nbloch), &
       rw(nbloch,nbloch),cw(nbloch,nbloch)
  complex(8) zw(nbloch,nbloch)
  real(8),allocatable:: rp(:,:),cp(:,:), rw1(:,:),cw1(:,:),rw2(:,:),cw2(:,:)
  logical :: debug = .false.
  integer(4):: verbose
  !-----------------------------------------
  if(verbose()>100) debug= .TRUE. 
  nbloch2    = nbloch*nbloch
  pi         = 4.d0*datan(1.d0)
  allocate ( rp(nbloch,nbloch), cp(nbloch,nbloch), &
       rw1(nbloch,nbloch),cw1(nbloch,nbloch), &
       rw2(nbloch,nbloch),cw2(nbloch,nbloch))
  !... Polarization function is decompsed into rp + i cp
  rp = dreal(zw(1:nbloch,1:nbloch))
  cp = dimag(zw(1:nbloch,1:nbloch))
  !... symmetrize the matrix
  call hermat  (rv,cv,nbloch)
  call hermat  (rp,cp,nbloch)
  !... calculate Pv
  if (debug) write(*,*)'debug> mmulc vP in'
  call mmulc   (rp,cp,nbloch, &
       rv,cv,nbloch, nbloch,nbloch,nbloch,nbloch, &
       rpv,cpv )
  rw1       = - rpv
  cw1       = - cpv
  do i      = 1,nbloch
     rw1(i,i)  = 1d0 + rw1(i,i)
  end do
  !... calculate vP
  if (debug) write(*,*)'debug> mmulc vP in'
  call mmulc   (rv,cv,nbloch, &
       rp,cp,nbloch, &
       nbloch,nbloch,nbloch,nbloch, &
       rw2,cw2 )
  rw2       = - rw2
  cw2       = - cw2
  do i      = 1,nbloch
     rw2(i,i)  = 1d0 + rw2(i,i)
  end do

  !... calculate (1-Pv)(1-vP)
  if (debug) write(*,*)'debug> mmulc vP in'
  call mmulc   (rw1,cw1,nbloch, &
       rw2,cw2,nbloch, &
       nbloch,nbloch,nbloch,nbloch, &
       rw,cw )

  !      if (debug) write(*,*)'debug> mmulc out'
  !      if (debug) then
  !        write(*,*)'*** P'
  !        call chhermat (rp,cp,nbloch)
  !        write(*,*)'*** v'
  !        call chhermat (rv,cv,nbloch)
  !        write(*,*)'*** Pv'
  !        call chhermat (rpv,cpv,nbloch)
  !        pause 'calcpv: check hermite'
  !      endif

  deallocate (rp,cp,rw1,cw1,rw2,cw2)
  return
end subroutine calcpv2
!---------------------------------------------------------------------
subroutine diagno0(nbloch,wpv,ovlp,evec,eval)
  implicit none
  integer(4):: nbloch,nmx,nev
  !      call diagno(nbloch,wpv,ovlp,wdiag,iwdiag,evec,eval)
  complex(8),allocatable:: ovlpc(:,:),wpvc(:,:),evecc(:,:)
  real(8)::emx
  real(8):: wpv(nbloch,nbloch,2),ovlp(nbloch,nbloch,2) &
       ,evec(nbloch,nbloch,2),eval(nbloch)
  allocate( ovlpc(nbloch,nbloch),wpvc(nbloch,nbloch), &
       evecc(nbloch,nbloch))
  ovlpc= dcmplx(ovlp(:,:,1), ovlp(:,:,2))
  wpvc = dcmplx( wpv(:,:,1), wpv (:,:,2))
  nev  = nbloch
  nmx  = nbloch
  !      print *,' goto diagcv-----------'
  call diagcv(ovlpc,wpvc, evecc, nbloch, eval, nmx, 1d99, nev)
  evec(:,:,1) =dreal(evecc)
  evec(:,:,2) =dimag(evecc)
  deallocate(wpvc,ovlpc,evecc)
end subroutine diagno0

!-------------------------------------------
subroutine getwk(ip,wibz,wqt,nqbz, nqibz,nstibz, nq0i, &
     wk)
  implicit none
  integer(4):: ip,nqibz,nq0i
  integer(4):: nqbz,nstibz(nqibz)
  real(8):: wibz(nqibz),wqt(nq0i),wk
  real(8)   :: wgtq0p
  !$$$      integer(4):: bzcase
  !- weight for k-point sampling
  !$$$      if(bzcase()==2)then
  !$$$        if(ip <=nqibz) then
  !$$$          wk = wibz(ip)*.5d0
  !$$$          if(nstibz(ip)/=0) wk = wibz(ip)*.5d0 *(1d0-wgtq0p()/nstibz(ip))
  !$$$        elseif(ip>nqibz) then
  !$$$          wk = wqt(ip-nqibz)* 1d0/dble(nqbz) * wgtq0p()
  !$$$        endif
  !$$$      else
  If(ip <= nqibz) then
     wk = wibz(ip)*.5d0 ! 0.5 for the normalization of wibz
  else
     wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
  endif
  if(abs(wibz(1)-2d0/dble(nqbz)) >1d-10) then
     print *,' wibz(1) nqbz=',wibz(1),nqbz
     print *,' sum wibz=',sum(wibz)
     print *, ' ecorq2: Bug Stop! this may be a bug?  wibz(1) /= 1/dble(2*nqbz)'
     ! top2rx 2013.08.09 kino          stop     ' ecorq2: Bug Stop! this may be a bug?  wibz(1) /= 1/dble(2*nqbz)'
     call rx( ' ecorq2: Bug Stop! this may be a bug?  wibz(1) /= 1/dble(2*nqbz)')
  endif
  !$$$      endif
  if(abs(sum(wibz)-2d0)>1d-10) then
     print *,' sum(wibz)=',sum(wibz)
     print *,' ecorq2: Bug Stop! this may be a bug? abs(sum(wibz)-2d0) /=0 '
     ! top2rx 2013.08.09 kino        stop    ' ecorq2: Bug Stop! this may be a bug? abs(sum(wibz)-2d0) /=0 '
     call rx( ' ecorq2: Bug Stop! this may be a bug? abs(sum(wibz)-2d0) /=0 ')
  endif
end subroutine getwk

!----------------------------------------------------------------
subroutine  ecorq2 (zv,zw, nbloch,  iq,iw,ieceig, &
     erpaqw, trpvqw, trlogqw)
  !- Contribution to Ec from given (q(ip),iw(ix)), Tr(log(1-Pv) +Pv) -----------
  ! Takao modified from ecorq2. Work arrays are embedded.
  ! from rpaq.f, July 24, 2001
  ! diagonalize (1-Pv)(1-vP) instead of Pv

  ! from GW/lw/hecor.F, July 05, 2001
  ! TM

  ! 00.07.03
  ! calculates the correlated part of the total energy
  ! Erpa  = 1/4pi Int[iw=-inf:inf] Tr{log(1-Pv)+Pv}
  !       = 1/2pi Int[iw=   0:inf] Tr{log(1-Pv)+Pv}
  ! Erpa is negative

  ! nbloch  = ngb total number of Bloch basis functions
  implicit real*8(a-h,o-z)
  implicit integer(i-n)
  real(8)::  ecqw
  real(8),allocatable :: rpv(:,:),cpv(:,:), &
       rw(:,:),cw(:,:),rv(:,:),cv(:,:), &
       wpv(:,:,:),ovlp(:,:,:), evec(:,:,:),eval(:),eval2(:)
  complex(8) zw(nbloch,nbloch), zv(nbloch,nbloch)
  data      nsngl,ndble/4,8/
  logical :: debug = .false.
  integer(4)::verbose
  !-----------------------------------------
  pi = 4d0*datan(1d0)
  if(verbose()>100) debug= .TRUE. 
  ngb= nbloch
  allocate( &
       rv(ngb,ngb),  cv(ngb,ngb), &
       rpv(ngb,ngb), cpv(ngb,ngb), &
       rw(ngb,ngb), cw(ngb,ngb), &
       wpv(ngb,ngb,2), ovlp(ngb,ngb,2), &
       !     &          wdiag(11*ngb), iwdiag(ngb),
       evec(ngb,ngb,2), eval(ngb), eval2(ngb) )
  rv = dreal(zv)
  cv = dimag(zv)

  !--- ix loop is now out of this routine.
  !'      do      ix = 1,niw
  call calcpv2 (zw,rv,cv, &
       nbloch, &
       rpv,cpv,rw,cw )
  if(debug) write(*,*)' calcpv2_out sumcheck =', sum(abs(rpv)),sum(abs(cpv))

  !--- diagonalize Pv
  do j = 1,nbloch
     do i = 1,nbloch
        wpv(i,j,1) = .5d0*(rpv(i,j)+rpv(j,i))
        wpv(i,j,2) = .5d0*(cpv(i,j)-cpv(j,i))
     enddo
  enddo

  ! cccccccccccccccccccccccccccccccccccccccccccccc
  ! eigenvlaue test1
  !        wpv(:,:,1) = dreal(zv)
  !        wpv(:,:,2) = dimag(zv)
  ! eigenvalue test2
  !        wpv(:,:,1) = dreal(zw)
  !        wpv(:,:,2) = dimag(zw)
  ! cccccccccccccccccccccccccccccccccccccccccccccc

  ovlp = 0d0
  do i = 1,nbloch
     ovlp(i,i,1)= 1d0
  enddo
  evec   = 0d0
  eval   = 0d0
  if (debug) write(*,*)'debug> diagno Pv in'
  call diagno0(nbloch,wpv,ovlp,evec,eval)

  !      wdiag  = 0.d0
  !      iwdiag = 0
  ! ccccccccccccccccccc
  !      call diagno(nbloch,wpv,ovlp,wdiag,iwdiag,evec,eval)
  !      do i = 1,nbloch
  !         write(6,"(i5,256d10.2)")i,(wpv(i,j,1),j=1,10)
  !         write(6,"(i5,256d10.2)")i,(wpv(i,j,2),j=1,10)
  !      enddo
  ! ccccccccccccccccccccc
  !      do i=1,nbloch
  !      do j=1,nbloch
  !         wpv(i,j,1:2)=0d0
  !         ovlp(i,j,1:2)=0d0
  !         if(i==j) wpv(i,j,1)=id0
  !         if(i==j) ovlp(i,j,1)=1d0
  !      enddo
  !      enddo
  !      evec=0d0
  !      eval=0d0
  !      stop 'xxxxxxxxxxxxxxxxxxxx test end xxxxxxx'
  ! ccccccccccccccccccc

  !--- diagonalize (1-Pv)(1-vP)
  do       j = 1,nbloch
     do       i = 1,nbloch
        wpv(i,j,1) = rw(i,j)
        wpv(i,j,2) = cw(i,j)
        ovlp(i,j,1)= 0d0
        ovlp(i,j,2)= 0d0
     enddo
  enddo
  do       i = 1,nbloch
     ovlp(i,i,1)= 1.d0
  enddo
  !      wdiag  = 0.d0 !      iwdiag = 0
  evec   = 0d0
  eval2  = 0d0
  if (debug) write(*,*)'debug> diagno (1-Pv)(1-vP) in'
  call diagno0(nbloch,wpv,ovlp,evec,eval2)

  trlogqw = 0.5d0*sum(dlog(eval2(1:nbloch)))
  trpvqw  = sum (eval(1:nbloch))
  erpaqw  = sum (.5d0*dlog(eval2(1:nbloch)) + eval(1:nbloch) )

  ! check write
  if(ieceig>0) then
     close(ieceig)
     open (ieceig, file='rpa_eigen.chk',access='append')
     write(ieceig,*)
     write(ieceig, &
          "('--- iq iw: ',2i6,' Eigen of  Pv+vP (1-Pv)(1-vP) ')") iq,iw
     do i = 1,nbloch
        write(ieceig,"(i5,d17.6,d17.6)")i,eval(i),eval2(i)
     enddo
  endif
  deallocate(rpv,cpv,rw,cw, wpv,ovlp,evec,eval,eval2)
  return
end subroutine ecorq2




!---------------------------------------------------------------------
subroutine sqrtmat (rmat,cmat, &
     wmat,ovlp,wdiag,iwdiag, &
     evec,eval,rw1,cw1, &
     ldim, &
     rmat2,cmat2)

  implicit double precision (a-h,o-z)
  implicit integer(i-n)

  dimension   rmat(ldim,ldim),cmat(ldim,ldim), &
       wmat(ldim,ldim,2),ovlp(ldim,ldim,2), &
       wdiag(11*ldim),iwdiag(ldim), &
       evec(ldim,ldim,2),eval(ldim), &
       rw1(ldim,ldim),cw1(ldim,ldim), &
       rmat2(ldim,ldim),cmat2(ldim,ldim)

  data tol /1.d-8/

  logical :: debug = .false.

  !-----------------------------------------
  do       j  = 1,ldim
     do       i  = 1,ldim
        wmat(i,j,1) = rmat(i,j)
        wmat(i,j,2) = cmat(i,j)
        ovlp(i,j,1) = 0.d0
        ovlp(i,j,2) = 0.d0
     enddo
  enddo
  do        i = 1,ldim
     ovlp(i,i,1) = 1.d0
  enddo

  wdiag       = 0.d0
  iwdiag      = 0
  evec        = 0.d0
  eval        = 0.d0

  ! cccccccccccccccccccccccccccccccc
  !      call diagno(ldim,wmat,ovlp,wdiag,iwdiag,evec,eval)
  call diagno0(ldim,wmat,ovlp,evec,eval)
  ! ccccccccccccccccccccccccccc

  ! v^{1/2}
  do        i = 1,ldim
     eval(i)     = dsqrt(eval(i))
  enddo

  ! U{-1}
  do       j  = 1,ldim
     do       i  = 1,ldim
        wmat(i,j,1) =  evec(j,i,1)
        wmat(i,j,2) = -evec(j,i,2)
     enddo
  enddo

  !      if (debug) then
  ! debug> check if evec is unitary
  !      do       i = 1,ldim
  !      do       j = 1,ldim
  !      rtmp       = 0.d0
  !      ctmp       = 0.d0
  !      do       k = 1,ldim
  !      rtmp       = rtmp
  !     .     + evec(i,k,1)*wmat(k,j,1) - evec(i,k,2)*wmat(k,j,2)
  !      ctmp       = ctmp
  !     .     + evec(i,k,1)*wmat(k,j,2) + evec(i,k,2)*wmat(k,j,1)
  !      enddo
  !      if ((dabs(rtmp).gt.tol) .or. (dabs(ctmp).gt.tol))
  !     .write(*,*)i,j,rtmp,ctmp
  !      enddo
  !      enddo
  !      pause 'sqrtmat: unitary'
  !      endif

  ! v^{1/2} U^{-1}
  do       j  = 1,ldim
     do       i  = 1,ldim
        rw1(i,j)    =  eval(i) * wmat(i,j,1)
        cw1(i,j)    =  eval(i) * wmat(i,j,2)
     enddo
  enddo

  ! U v^{1/2} U^{-1}
  do        j = 1,ldim
     do        i = 1,ldim
        rtmp        = 0.d0
        ctmp        = 0.d0
        do        k = 1,ldim
           rtmp        = rtmp &
                + wmat(k,i,1)*rw1(k,j) + wmat(k,i,2)*cw1(k,j)
           ctmp        = ctmp &
                + wmat(k,i,1)*cw1(k,j) - wmat(k,i,2)*rw1(k,j)
        end do
        rmat2(i,j)  = rtmp
        cmat2(i,j)  = ctmp
     enddo
  enddo

  !      if (debug) then
  ! ebug> check mat2*mat2 = mat
  !      do        i = 1,ldim
  !      write(*,*)'i=',i
  !      do        j = 1,ldim
  !      wmat(i,j,1) = 0.d0
  !      wmat(i,j,2) = 0.d0
  !      do        k = 1,ldim
  !      wmat(i,j,1) = wmat(i,j,1)
  !     .             + rmat2(i,k)*rmat2(k,j) - cmat2(i,k)*cmat2(k,j)
  !      wmat(i,j,2) = wmat(i,j,2)
  !     .             + rmat2(i,k)*cmat2(k,j) + cmat2(i,k)*rmat2(k,j)
  !      enddo
  !      rdiff       = wmat(i,j,1) - rmat(i,j)
  !      cdiff       = wmat(i,j,2) - cmat(i,j)
  !      if ((dabs(rdiff).gt.tol) .or. (dabs(cdiff).gt.tol))
  !     .write(*,*)i,j,rdiff,cdiff
  !      enddo
  !      enddo
  !      pause 'check mat2 out'
  !      endif

  return
end subroutine sqrtmat
!---------------------------------------------------------------------
subroutine hermat  (rmat,cmat,ldim)
  implicit double precision (a-h,o-z)
  implicit integer(i-n)
  real(8)   rmat(ldim,ldim),cmat(ldim,ldim)

  do      i = 1,ldim
     do      j = i,ldim
        w1        = 0.5d0 * (rmat(i,j) + rmat(j,i))
        w2        = 0.5d0 * (cmat(i,j) - cmat(j,i))
        rmat(i,j) =  w1
        cmat(i,j) =  w2
        rmat(j,i) =  w1
        cmat(j,i) = -w2
     end do
  end do

  return
end subroutine hermat
!$$$c---------------------------------------------------------------------
!$$$      subroutine chhermat (rmat,cmat,ldim)
!$$$
!$$$      implicit double precision (a-h,o-z)
!$$$
!$$$      dimension   rmat(ldim,ldim),cmat(ldim,ldim)
!$$$      data tol /1.d-8/
!$$$
!$$$      do i = 1,ldim
!$$$        do j = 1,ldim
!$$$          w1   = rmat(i,j) - rmat(j,i)
!$$$          w2   = cmat(i,j) + cmat(j,i)
!$$$          if ((dabs(w1).gt.tol).or.(dabs(w2).gt.tol)) then
!$$$            write(*,*)i,j,w1,w2
!$$$            pause 'chhermat: non-hermite!'
!$$$          endif
!$$$        end do
!$$$      end do
!$$$
!$$$      return
!$$$      end
!$$$c---------------------------------------------------------------------
subroutine mmulc (ra,ca,lda, &
     rb,cb,ldb, &
     nrow,nmul,ncol,ldc, &
     rc,cc)

  ! 91.11.29
  ! multiply two complex matrices a b = c

  ! ra,ca = real and imaginary parts of a
  ! rb,cb =                             b
  ! lda,ldb,ldc = leading dimensions of a,b and c
  ! nrow,ncol   = no. rows and coulmns of c
  ! nmul  = no. contractions

  ! rc,cc = real and imaginary parts of c

  implicit double precision (a-h,o-z)
  implicit integer(i-n)

  dimension ra(lda,1),ca(lda,1), &
       rb(ldb,1),cb(ldb,1)
  dimension rc(ldc,1),cc(ldc,1)

  ! top2rx 2013.08.09 kino      if(nrow .gt. lda) stop 'mmulc: lda too small'
  if(nrow > lda) call rx( 'mmulc: lda too small')
  ! top2rx 2013.08.09 kino      if(nmul .gt. ldb) stop 'mmulc: ldb too small'
  if(nmul > ldb) call rx( 'mmulc: ldb too small')
  ! top2rx 2013.08.09 kino      if(nmul .gt. ldc) stop 'mmulc: ldc too small'
  if(nmul > ldc) call rx( 'mmulc: ldc too small')
  !     do      ir = 1,nrow
  !     do      ic = 1,ncol
  !     rsum       = 0.d0
  !     csum       = 0.d0
  !     do       i = 1,nmul
  !     rsum       = rsum + ra(ir,i)*rb(i,ic)
  !    .                  - ca(ir,i)*cb(i,ic)
  !     csum       = csum + ra(ir,i)*cb(i,ic)
  !    .                  + ca(ir,i)*rb(i,ic)
  !     end do
  !     rc(ir,ic)  = rsum
  !     cc(ir,ic)  = csum
  !     end do
  !     end do

  do      ic = 1,ncol

     do      ir = 1,nrow
        rc(ir,ic)  = 0.d0
        cc(ir,ic)  = 0.d0
     enddo

     do       i = 1,nmul
        rbic       = rb(i,ic)
        cbic       = cb(i,ic)
        do      ir = 1,nrow
           rc(ir,ic)  = rc(ir,ic) + rbic*ra(ir,i) - cbic*ca(ir,i)
           cc(ir,ic)  = cc(ir,ic) + cbic*ra(ir,i) + rbic*ca(ir,i)
        end do
     end do

  end do

  return
end subroutine mmulc

end module m_rpaq
