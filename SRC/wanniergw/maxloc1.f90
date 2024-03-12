subroutine init_unkg(is,qbz,ginv,ef,lein, &
     iko_ix,iko_fx,iko_i,iko_f, &
     iki_i,iki_f, &
     nwf,nband,nqbz, &
     amnk,cnk)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  logical :: lein
  !      complex(8),allocatable :: psig(:,:,:)
  complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz), &
       cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
  real(8) :: qbz(3,nqbz),ginv(3,3), &
       ovlp(iko_ix:iko_fx)
  integer(4) :: iko_i(nqbz),iko_f(nqbz), &
       iki_i(nqbz),iki_f(nqbz)
  ! initialize
  cnk = (0d0,0d0)
  amnk = (0d0,0d0)

  ! read psig(it,iwf,iqbz) = < psi(it,iqbz) | g(iwf) >
  print *,'goto readpsig'
  call readpsig(is,iko_ix,iko_fx, &
       nqbz,nwf, &
       !     o              psig)
       amnk)

  call amnk2unk(amnk, &
       iko_ix,iko_fx,iko_i,iko_f, &
       nwf,nqbz, &
       cnk)

  ! inner energy window
  if (lein) &
       call init_iew(iko_ix,iko_fx,iko_i,iko_f, &
       iki_i,iki_f, &
       nwf,nband,nqbz, &
       cnk)

  return
end subroutine init_unkg
!-----------------------------------------------------------------------
subroutine pick_nwf(ovlp,iti,itf,nwf,isort)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  real(8) :: ovlp(iti:itf),otmp(iti:itf)
  integer(4) :: isort(nwf),istate(iti:itf)

  ! initial
  do it = iti,itf
     istate(it) = it
     otmp(it) = ovlp(it)
  enddo

  ! sorting
  do it1 = iti,itf-1
     do it2 = it1+1,itf
        if (ovlp(it1) < ovlp(it2)) then
           tmp = ovlp(it2)
           ovlp(it2) = ovlp(it1)
           ovlp(it1) = tmp
           itmp = istate(it2)
           istate(it2) = istate(it1)
           istate(it1) = itmp
        endif
     enddo
  enddo

  ! sort check
  do i1 = iti,itf-1
     it1 = istate(i1)
     it2 = istate(i1+1)
     if (otmp(it1) < otmp(it2)) stop 'pick_nwf: sort error'
  enddo

  ! pick largest nwf states
  do it = 1,nwf
     itmp = iti - 1 + it
     isort(it) = istate(itmp)
  enddo

  return
end subroutine pick_nwf
!-----------------------------------------------------------------------
subroutine read_cnq0(ifhoev,is,qwf0,qbz,ginv,ef, &
     itq, &
     nwf,nband,nqbz, &
     cnq0)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:)
  complex(8) :: cnq0(nband,nwf)
  real(8),allocatable :: eval(:)
  real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
  real(8) :: rydberg
  integer(4) :: itq(nwf)

  iq0 = iqindx(qwf0,ginv,qbz,nqbz)
  cnq0 = (0d0,0d0)

  ! open
  if (is == 1) then
     ifhoev = iopen('HOEV.UP',0,0,0)
  elseif (is == 2) then
     ifhoev = iopen('HOEV.DN',0,0,0)
  else
     stop 'read_cnq0: iopen error'
  endif

  ! read
  read(ifhoev)ndimh,nqtot
  if (ndimh /= nband) stop 'read_cnq0: nband error'
  if (nqtot < nqbz) stop 'read_cnq0: nqbz error'

  allocate(hks(nband,nband),oks(nband,nband), &
       cks(nband,nband),eval(nband))

  do iq = 1,nqbz
     read(ifhoev)iq2,q(1:3)
     read(ifhoev)hks(1:nband,1:nband)
     read(ifhoev)oks(1:nband,1:nband)
     read(ifhoev)cks(1:nband,1:nband)
     read(ifhoev)eval(1:nband)

     iq3 = iqindx(q,ginv,qbz,nqbz)
     if (iq3 /= iq) stop 'read_cnq0: iqindx error'
     if (iq3 == iq0) then
        do it = 1,nwf
           cnq0(:,it) = cks(:,itq(it))
           !               ev = (eval(itq(it)) - ef) * rydberg()
           !               write(*,*)'iwf,nwf,ev',ev
        enddo
        goto 99
     endif
  enddo
  stop 'read_cnq0: cannot find q0'
99 continue

  deallocate(hks,oks,cks,eval)

  ! close
  if (is == 1) then
     ifi = iclose('HOEV.UP')
  else
     ifi = iclose('HOEV.DN')
  endif

  return
end subroutine read_cnq0
!-----------------------------------------------------------------------
subroutine get_amnk(ifhoev,is,qwf0,qbz,ginv, &
     cnq0, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nwf,nband,nqbz, &
     amnk)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:), &
       wmat(:,:)
  complex(8) :: cnq0(nband,nwf), &
       amnk(iko_ix:iko_fx,nwf,nqbz)
  real(8),allocatable :: eval(:)
  real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
  integer(4) :: iko_i(nqbz),iko_f(nqbz)


  ! open
  if (is == 1) then
     ifhoev = iopen('HOEV.UP',0,0,0)
  elseif (is == 2) then
     ifhoev = iopen('HOEV.DN',0,0,0)
  else
     stop 'get_amnk: iopen error'
  endif

  ! read
  read(ifhoev)ndimh,nqtot
  if (ndimh /= nband) stop 'get_amnk: nlmto error'
  if (nqtot < nqbz) stop 'get_amnk: nqbz error'

  allocate(hks(nband,nband),oks(nband,nband), &
       cks(nband,nband),eval(nband), &
       wmat(nband,nwf))

  ! initialize
  amnk = (0d0,0d0)

  do iq = 1,nqbz
     read(ifhoev)iq2,q(1:3)
     read(ifhoev)hks(1:nband,1:nband)
     read(ifhoev)oks(1:nband,1:nband)
     read(ifhoev)cks(1:nband,1:nband)
     read(ifhoev)eval(1:nband)

     iq3 = iqindx(q, ginv,qbz,nqbz)
     if (iq3 /= iq) stop 'get_amnk: iqindx error'


     ! wmat = cnq0 * oks
     wmat = (0d0,0d0)
     do in = 1,nwf
        do ij = 1,nband
           do ii = 1,nband
              wmat(ij,in) = wmat(ij,in) &
                   + cnq0(ii,in) * oks(ij,ii)
           enddo
        enddo
     enddo

     ! amnk = cks^{*} * wmat
     do im = iko_i(iq),iko_f(iq)
        do in = 1,nwf
           do ij = 1,nband
              amnk(im,in,iq) = amnk(im,in,iq) + &
                   dconjg(cks(ij,im)) * wmat(ij,in)
           enddo
        enddo
     enddo
  enddo
  deallocate(hks,oks,cks,eval,wmat)
  ! close
  if (is == 1) then
     ifi = iclose('HOEV.UP')
  else
     ifi = iclose('HOEV.DN')
  endif
999 format(i5,3f16.8)
  return
end subroutine get_amnk
!-----------------------------------------------------------------------
subroutine amnk2unk(amnk,iko_ix,iko_fx,iko_i,iko_f,nwf,nqbz, cnk0)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  parameter (eps = 1d-3)
  intent(in)::     amnk,iko_ix,iko_fx,iko_i,iko_f,nwf,nqbz
  intent(out)::                                              cnk0
  complex(8),allocatable :: aa(:,:),cc(:,:),zz(:,:),vv(:,:)
  real(8),allocatable :: dd(:)
  complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz), cnk0(iko_ix:iko_fx,nwf,nqbz),ctmp
  integer :: iko_i(nqbz),iko_f(nqbz)
  !! singular value decomposition, 061003. see around Eq.23 of Ref.[2]
  nks = iko_fx - iko_ix + 1
  allocate (aa(nks,nks),zz(nks,nks),vv(nwf,nwf),dd(nwf))
  cnk0 = 0d0
  do iq = 1,nqbz
     aa(1:nks,1:nwf) = amnk(iko_ix:iko_fx,1:nwf,iq)
     call zgesvdmn(nks,nwf,aa,dd,zz,vv)
     do ij = iko_i(iq),iko_f(iq)
        jj = ij - iko_ix + 1
        do ii = 1,nwf
           do ik = 1,nwf
              cnk0(ij,ii,iq) = cnk0(ij,ii,iq) + zz(jj,ik)*vv(ik,ii)
           enddo
        enddo
     enddo
  enddo 
  deallocate (aa,zz,vv,dd)
  return
end subroutine amnk2unk
!-----------------------------------------------------------------------
subroutine init_iew(iko_ix,iko_fx,iko_i,iko_f, iki_i,iki_f, nwf,nband,nqbz, cnk) !Modify cnk to include all the inner space explictly. See Souza Eq.27 around. 
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8),allocatable :: cnk2(:,:),vnk(:,:),mat(:,:),evec(:,:)
  complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
  real(8),allocatable :: eval(:)
  integer(4),allocatable :: it(:)
  integer(4) :: iko_i(nqbz),iko_f(nqbz), iki_i(nqbz),iki_f(nqbz)
  allocate(cnk2(iko_ix:iko_fx,nwf), vnk(iko_ix:iko_fx,nwf))
  nox = iko_fx - iko_ix + 1
!  print *, 'initnnnnnnnnox=',nox,iko_ix, iko_fx,nwf
  do iq = 1,nqbz
 !    print *,'initnnnnnnnnn iq=',iq,iki_i(iq),iki_f(iq)
     ! no innner energy window
     if (iki_i(iq) == 0) then
        if (iki_f(iq) /= -1) stop 'init_iew: iki_f error'
     else        ! inner energy window
        nout = iko_f(iq) - iko_i(iq) + 1
        nin  = iki_f(iq) - iki_i(iq) + 1
        if (nin < 1) stop 'init_iew: nin error'
        cnk2 = 0d0
        do il = 1, nin ! Nin(k) states in the inner energy window
           cnk2(iki_i(iq)-1+il,il) = 1d0
        enddo
        ! Nwf - Nin(k) states out of Nout states         ! |vnk>
        vnk(:,:) = cnk(:,:,iq)
        vnk(iki_i(iq):iki_f(iq),:) = 0d0
!        print *,'initnnnn222 ',nox,nwf,iko_ix
        allocate(mat(nox,nox),evec(nox,nox),eval(nox))
        do ii = 1,nox !nox=nouter-ninner
           do ij = 1,nox
              mat(ii,ij) = sum(vnk(ii+iko_ix-1,1:nwf)*dconjg(vnk(ij+iko_ix-1,1:nwf))) !QPQ matrix of Eq.27 in Souza paper.
           enddo
        enddo
        forall (ii=1:nox) mat(ii,ii)=mat(ii,ii)+1d-16 !stabilize calculaiton 2023-6-6
        call chk_hm(mat,nox)
        call diag_hm(mat,nox,eval,evec)
 !       do i=1,nox
 !          write(6,*)'nnnnnnneval=',i,eval(i)
 !       enddo
        do il = 1, nwf-nin
           do ij = iko_ix, iko_fx
              cnk2(ij,il+nin) = evec(ij-iko_ix+1,nox-il+1)
           enddo
        enddo
        deallocate(mat,evec,eval)
        call chk_on(cnk2(iko_ix:iko_fx,:),nox,nwf)
        ! new cnk(:,:,q)
        cnk(:,:,iq) = cnk2(:,:)
        ! end of if (iki_i == 0)
     endif     ! end of iq-loop
     !xxxx debug if(iq==6) stop 'xxxxxxxxxxxx'
  enddo
  deallocate(cnk2,vnk)
end subroutine init_iew
!-----------------------------------------------------------------------
subroutine   getupu(isc, &
     uumat,cnk, &
     lein,alpha_in,iq,ikbidx, &
     iko_ix,iko_fx, &
     iko_i,iko_f, &
     iki_i,iki_f, & !inner window range
     ikbo_i,ikbo_f, &
     ikbi_i,ikbi_f, &
     nwf,nbb,nqbz, &
     upu)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8),allocatable :: wmat(:,:),wmat2(:,:)
  complex(8) :: upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb), &
       uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb), &
       cnk(iko_ix:iko_fx,nwf,nqbz), &
       ctmp
  integer :: ikbidx(nbb), &
       ikbo_i(nbb),ikbo_f(nbb), &
       ikbi_i(nbb),ikbi_f(nbb)
  logical :: lein
  if (isc == 1) then
     alpha = 1d0
     upu(:,:,:) = 0d0
  else
     alpha = alpha_in
     upu(:,:,:) = upu(:,:,:) * (1d0-alpha)
  endif
!  nin = iki_f - iki_i + 1 !starting index of outer-inner (ex. When iki_i=1, nin= iki_f+1)
!  if (nin >= nwf) stop 'getupu: Nin >= Nwf'
  allocate(wmat(iko_ix:iko_fx,iko_ix:iko_fx), wmat2(iko_ix:iko_fx,iko_ix:iko_fx))
  do ibb = 1,nbb
     iqb = ikbidx(ibb)
     i1= ikbo_i(ibb)
     i2= ikbo_f(ibb)
     do concurrent(inp=i1:i2, imp=i1:i2) !wmat = cnk * cnk^{*} is projector to 'wannier space'.
        wmat(inp,imp)= sum(dconjg(cnk(inp,1:nwf,iqb))*cnk(imp,1:nwf,iqb)) !miyake BUG-> range of sum was nin+1,nwf before 2023-6-8
     enddo
     do concurrent(in=iko_i:iko_f, imp=i1:i2)
        wmat2(imp,in)= sum(dconjg(uumat(in,i1:i2,ibb)) * wmat(i1:i2,imp))
     enddo
     do concurrent(im=iko_i:iko_f, in=iko_i:iko_f)! uumat* wmat * uumat
        upu(im,in,ibb) = upu(im,in,ibb) +  sum(uumat(im,i1:i2,ibb)*wmat2(i1:i2,in)) * alpha  ! upu = upu + uumat * wmat2 * alpha
     enddo
  enddo
  deallocate(wmat,wmat2)
end subroutine getupu
!-----------------------------------------------------------------------
subroutine dimz(lein,iko_i,iko_f,iki_i,iki_f, ndz,nin)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  logical :: lein
  nout = iko_f - iko_i + 1
  nin  = iki_f - iki_i + 1
  ndz  = nout - nin
  ! eck
  if (iki_i == 0) then
     if (iki_f /= -1) stop 'dimz: iki_f error'
  else
     if (iko_i > iki_i) stop 'dimz: ik_i error'
     if (iko_f < iki_f) stop 'dimz: ik_f error'
     if (iki_i > iki_f) stop 'dimz: iki error'
  endif
  return
end subroutine dimz
!-----------------------------------------------------------------------
subroutine getzmn(upu,wbb,lein, &
     iko_ix,iko_fx, &
     iko_i,iko_f, &
     iki_i,iki_f, &
     nwf,nbb,nqbz,ndz, &
     zmn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8) :: upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb),    zmn(ndz,ndz)
  real(8) :: wbb(nbb)
  integer(4) :: it(ndz)
  logical :: lein
  zmn = (0d0,0d0)
  ! no inner energy window
  if (iki_i == 0) then
!     no = iko_f - iko_i + 1
!     if (no /= ndz) stop 'getzmn: ndz error'
     do ibb = 1, nbb
        zmn(1:ndz,1:ndz) = zmn(1:ndz,1:ndz) +  wbb(ibb) * upu(iko_i:iko_f,iko_i:iko_f,ibb)
     enddo
  else ! inner energy window
     j = 0
     do i = iko_i,iki_i-1
        j = j + 1
        it(j) = i
     enddo
     do i = iki_f+1,iko_f
        j = j + 1
        it(j) = i
     enddo
     if (j /= ndz) stop 'getzmn: ndz error'
     do im = 1,ndz
        do in = 1,ndz
           do ibb = 1, nbb
              zmn(im,in) = zmn(im,in) + wbb(ibb) * upu(it(im),it(in),ibb)
           enddo
        enddo
     enddo
  endif
end subroutine getzmn
!-----------------------------------------------------------------------
subroutine chk_hm(zmat,ndim)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps = 1d-4)
  complex(8):: zmat(ndim,ndim)

  do i = 1,ndim
     do j = 1,ndim
        dr = dabs(dreal(zmat(i,j)) - dreal(zmat(j,i)))
        di = dabs(dimag(zmat(i,j)) + dimag(zmat(j,i)))
        dc = abs(zmat(i,j))
        if (dr > eps) stop 'chk_hm: real part error'
        if (di > eps) stop 'chk_hm: imag part error'
        !         if (dr/dc .gt. eps) stop 'chk_hm: real part error'
        !         if (di/dc .gt. eps) stop 'chk_hm: imag part error'
        !         if (dr/dc .gt. eps) then
        !             write(*,*)zmat(i,j),i,j
        !             write(*,*)zmat(j,i),i,j
        !             stop 'chk_hm: real part error'
        !         endif
        !         if (di/dc .gt. eps) then
        !             write(*,*)zmat(i,j),i,j
        !             write(*,*)zmat(j,i),i,j
        !             stop 'chk_hm: imag part error'
        !         endif
     enddo
  enddo

  return
end subroutine chk_hm
!-----------------------------------------------------------------------
subroutine chk_um(zmat,ndim)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps = 1d-4)
  complex(8):: zmat(ndim,ndim),cij

  do i = 1,ndim
     do j = 1,ndim
        cij = (0d0,0d0)
        do k = 1,ndim
           cij = cij + dconjg(zmat(k,i))*zmat(k,j)
        enddo
        if (i == j) cij = cij - 1d0
        rij = abs(cij)
        if (rij > eps) stop 'chk_um: error'
     enddo
  enddo

  return
end subroutine chk_um
!-----------------------------------------------------------------------
subroutine chk_on(zmat,n1,n2)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-4)
  complex(8):: zmat(n1,n2),cij
  do i = 1,n2
     do j = 1,n2
        cij = sum(dconjg(zmat(:,i))*zmat(:,j))
        if (i == j) cij = cij - 1d0
        rij = abs(cij)
        if (rij > eps) stop 'chk_on: error'
     enddo
  enddo
end subroutine chk_on
!-----------------------------------------------------------------------
subroutine diag_hm(zmat,ndim,eval,evecc)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
  complex(8):: zmat(ndim,ndim),evecc(ndim,ndim)
  real(8):: eval(ndim)

  allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))

  nev  = ndim
  nmx  = ndim

  zmat2 = zmat

  ovlpc = (0d0,0d0)
  do i=1,ndim
     ovlpc(i,i) = (1d0,0d0)
  enddo

  evecc = (0d0,0d0)
  eval = 0d0

  call diagcv(ovlpc,zmat2, evecc, ndim, eval, nmx, 1d99, nev)

  deallocate(zmat2,ovlpc)

  return
end subroutine diag_hm
!-----------------------------------------------------------------------
subroutine chk_eval(wbb,evz,nbb,ndz)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps = 1d-4)
  real(8):: wbb(nbb),evz(ndz)


  ! check order
  do i = 1,ndz-1
     e1 = evz(i)
     e2 = evz(i+1)
     if (e1 > e2) stop 'chk_eval: order is wrong'
  enddo

  ! check
  ws = sum(wbb)
  if (evz(ndz) > ws) then
     write(*,*) 'chk_eval: eval is too large'
     write(*,*)'sum(wbb) =',ws
     do i = ndz,1,-1
        write(*,*)'i,evz(i)=',i,evz(i)
     enddo
     stop
  endif

  return
end subroutine chk_eval
!-----------------------------------------------------------------------
subroutine new_cnk(cnk,evecc,iq, &
     iko_ix,iko_fx, &
     iko_i,iko_f, &
     iki_i,iki_f, &
     nwf,ndz, &
     cnk2)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8) :: cnk(iko_ix:iko_fx,nwf), cnk2(iko_ix:iko_fx,nwf), evecc(ndz,ndz)
  integer(4) :: it(ndz)
  cnk2(:,:) = 0d0
  if (iki_i == 0) then! no inner energy window
     nout = iko_f - iko_i + 1
     if (nout /= ndz) stop 'new_cnk: ndz error'
     if (nwf > ndz) stop 'new_cnk: nwf error'
     il2 = ndz + 1
     do il = 1,nwf
        il2 = il2 - 1
        cnk2(iko_i:iko_f,il) = evecc(1:ndz,il2)
     enddo
  else ! with inner energy window
     nout = iko_f - iko_i + 1
     nin  = iki_f - iki_i + 1
     if (nout-nin /= ndz) stop 'new_cnk: ndz error'
     if (nwf-nin > ndz) stop 'new_cnk: nwf error'
     cnk2(:,1:nin) = cnk(:,1:nin) !inner window
     j = 0
     do i = iko_i,iki_i-1
        j = j + 1
        it(j) = i
     enddo
     do i = iki_f+1,iko_f
        j = j + 1
        it(j) = i
     enddo
     il2 = ndz + 1
     do il = nin+1,nwf ! pick nwf-nin states with the largest eigenvalues
        il2 = il2 - 1
        do in = 1,ndz
           cnk2(it(in),il) = evecc(in,il2)
        enddo
     enddo
  endif
end subroutine new_cnk
!-----------------------------------------------------------------------
subroutine  get_omgik(wbb,evz, &
     iko_i,iko_f, &
     iki_i,iki_f, &
     nbb,nwf,ndz, &
     omgik)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-4)
  real(8):: wbb(nbb),evz(ndz)
  nin = iki_f - iki_i + 1
  nn  = nwf - nin
  if (ndz < nn) stop 'get_omgik: ndz error'
  esum = 0d0
  do i = 1,nn
     j = ndz + 1 - i
     esum = esum + evz(j)
  enddo
  omgik = dble(nn)*sum(wbb) - esum
end subroutine get_omgik
!-----------------------------------------------------------------------
subroutine zgesvdmn(ngb1,ngb2,zzz, SS,UU,VT)
  implicit none
  integer(4)::lwork,info,ngb1,ngb2,i
  complex(8):: zzz(ngb1,ngb2),UU(ngb1,ngb1),VT(ngb2,ngb2)
  real(8):: ss(ngb2)
  real(8),allocatable:: rwork(:)
  complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
  lwork=4*ngb1
  allocate(work(LWORK),rwork(5*ngb1))
  call zgesvd('A','A',ngb1,ngb2,zzz,ngb1,SS,UU,ngb1,VT,ngb2,work,lwork,rwork,info)
  deallocate(work,rwork)
end subroutine zgesvdmn
!$$$!--------------------------------------------------------------------------
!$$$      subroutine zgesvdnn(ngb,zzz, SS,UU,VT)
!$$$C--- SVD zzz= matmul(UU,matmul(SS,VT)) ------------
!$$$c$$$C--- SVD of chipm test !SVD procedure is not required to calculate <m|chi|m>
!$$$c$$$      lwork=4*ngb
!$$$c$$$      allocate(work(LWORK),rwork(5*ngb),zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb),VTT(ngb,ngb),ss0(ngb))
!$$$c$$$      zzz= matmul(transpose(conjg(ovlpi)), matmul(zxq(1:ngb,1:ngb,iw), ovlpi) )
!$$$c$$$      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS0,UU,ngb,VT,ngb,work,lwork,rwork,info)
!$$$c$$$      write(6,*)' -------SVD: Oinv*chipm*Oinv ----------------'
!$$$c$$$      aaax = 0d0
!$$$c$$$      do i=1,ngb
!$$$c$$$        addx= sum(svec(1:nbloch)*uu(1:nbloch,i)) *ss0(i)* sum(VT(i,1:nbloch)*svec(1:nbloch))
!$$$c$$$        write(6,"(' SVD OcO: eig_k <m|chi|m>_k=',i4,2x, d13.5,2x,2d14.6)")i,SS0(i),addx
!$$$c$$$        if(i<25) aaax= aaax+ addx
!$$$c$$$      enddo
!$$$c$$$      aaax= mmnorm**2/aaax
!$$$c$$$      deallocate(work,rwork,zzz,uu,vt,vtt)
!$$$c$$$      deallocate(ovlpi)
!$$$      implicit none
!$$$      integer(4)::lwork,info,ngb,i
!$$$      complex(8):: zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb)
!$$$      real(8):: ss(ngb)
!$$$      real(8),allocatable:: rwork(:)
!$$$      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
!$$$      lwork=4*ngb
!$$$      allocate(zw0bk(ngb,ngb))
!$$$      allocate(work(LWORK),rwork(5*ngb)) !,VTT(ngb,ngb))
!$$$      zw0bk = zzz
!$$$!      write(6,*)' zgesvdnn: singular value decomp '
!$$$      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS,UU,ngb,VT,ngb,work,lwork,rwork,info)
!$$$!      do i=1,ngb
!$$$!         write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!$$$!         vtt(i,:)=ss(i)*vt(i,:)
!$$$!      enddo
!$$$!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!$$$!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!$$$!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
!$$$!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!$$$!      deallocate(vtt)
!$$$      end
!$$$!--------------------------------------------------------------------------
