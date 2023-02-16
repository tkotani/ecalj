subroutine diag_unk(is,qbz, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nband,nwf,nqbz, &
     cnk, &
     eunk)
  use m_readeigen,only:readeval
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-3)
  complex(8),allocatable :: evecc(:,:),ham(:,:),cnk2(:,:)
  real(8),allocatable :: eval(:),eks(:)
  complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz), cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
  real(8) :: qbz(3,nqbz),q(3),eunk(nwf,nqbz)
  integer :: iko_i(nqbz),iko_f(nqbz)
  allocate (ham(nwf,nwf),eks(nband), evecc(nwf,nwf),eval(nwf), cnk2(iko_ix:iko_fx,nwf))
  do iq = 1,nqbz
     ! read eigenvalues
     q(:) = qbz(:,iq)
     eks= readeval (q,is)

     ! construct H
     ham = (0d0,0d0)
     do ii = 1,nwf
        do ij = 1,nwf
           do ik = iko_i(iq),iko_f(iq)
              ham(ii,ij) = ham(ii,ij) + &
                   dconjg(cnk(ik,ii,iq))*cnk(ik,ij,iq)*eks(ik)
           enddo
        enddo
     enddo

     ! diagonalize H
     call chk_hm(ham,nwf)
     call diag_hm(ham,nwf,eval,evecc)


     ! cnk(n,l) = S[m] evecc(m,l)*c(n,m)
     cnk2 = (0d0,0d0)
     do il = 1,nwf
        do im = 1,nwf
           do in = iko_i(iq),iko_f(iq)
              cnk2(in,il) = cnk2(in,il) + &
                   evecc(im,il) * cnk(in,im,iq)
           enddo
        enddo
     enddo

     cnk(:,:,iq) = cnk2
     eunk(:,iq) = eval(:)

     ! end of iq-loop
  enddo

  deallocate (ham,eks,evecc,eval,cnk2)

  return
end subroutine diag_unk
!-----------------------------------------------------------------------
subroutine chk_eunk(is,qbz,eunk,ef, &
     nqbz,nband,nwf)

  use m_readeigen,only:readeval
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  real(8) :: qbz(3,nqbz),eunk(nwf,nqbz)
  real(8),allocatable :: eks(:)

  allocate(eks(nband))

  do iq = 1,nqbz
     eks= readeval (qbz(:,iq),is)

     write(*,*)'Diag energy',iq
     do iband = 1,nwf
        eev = (eunk(iband,iq)-ef)*rydberg()
        write(*,*)iband,eunk(iband,iq),eev
     enddo

     write(*,*)'KS energy  ',iq
     do iband = 1,nband
        eev = (eks(iband)-ef)*rydberg()
        write(*,*)iband,eks(iband),eev
     enddo
  enddo

  deallocate(eks)

  return
end subroutine chk_eunk
!-----------------------------------------------------------------------
subroutine chk_cnk(cnk, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nband,nwf,nqbz)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps = 1d-4)
  complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
  integer :: iko_i(nqbz),iko_f(nqbz)

  do iq = 1,nqbz

     do im = 1,nwf
        do in = 1,nwf
           ctmp = (0d0,0d0)
           do ii = iko_i(iq),iko_f(iq)
              ctmp = ctmp + dconjg(cnk(ii,im,iq))*cnk(ii,in,iq)
           enddo
           if (in == im) ctmp = ctmp - 1d0
           a = dabs(dreal(ctmp))
           b = dabs(dimag(ctmp))
           if (a > eps) then
              write(*,*)'chk_cnk: real error,iq,im,in, Re'
              write(*,*)iq,im,in,a
           endif
           if (b > eps) then
              write(*,*)'chk_cnk: imag error,iq,im,in, Im'
              write(*,*)iq,im,in,b
           endif
        enddo
     enddo

  enddo

  stop 'chk_cnk: end'

  return
end subroutine chk_cnk
!-----------------------------------------------------------------------
subroutine init_mmn(cnk,uumat,ikbidx, &
     iko_ix,iko_fx,iko_i,iko_f,ikbo_i,ikbo_f, &
     nwf,nqbz,nbb, &
     mmn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8),allocatable :: wmat(:,:)
  complex(8):: cnk(iko_ix:iko_fx,nwf,nqbz), &
       uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz), &
       mmn(nwf,nwf,nbb,nqbz)
  integer(4):: ikbidx(nbb,nqbz), &
       iko_i(nqbz),iko_f(nqbz), &
       ikbo_i(nbb,nqbz),ikbo_f(nbb,nqbz)

  allocate(wmat(iko_ix:iko_fx,nwf))

  mmn = (0d0,0d0)
  do iq = 1,nqbz
     do ibb = 1,nbb
        iqb = ikbidx(ibb,iq)

        ! wmat = cnk * uumat
        wmat = (0d0,0d0)
        do imp = iko_i(iq),iko_f(iq)
           do in = 1,nwf
              do inp = ikbo_i(ibb,iq),ikbo_f(ibb,iq)
                 wmat(imp,in) = wmat(imp,in) &
                      + cnk(inp,in,iqb) * uumat(imp,inp,ibb,iq)
              enddo
           enddo
        enddo

        ! mmn = cnk^{*} * wmat
        do im = 1,nwf
           do in = 1,nwf
              do imp = iko_i(iq),iko_f(iq)
                 mmn(im,in,ibb,iq) = mmn(im,in,ibb,iq) &
                      + dconjg(cnk(imp,im,iq)) * wmat(imp,in)
              enddo
           enddo
        enddo

     enddo
  enddo

  deallocate(wmat)

  return
end subroutine init_mmn
!-----------------------------------------------------------------------
subroutine init_Umnk(amnk,cnk, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nwf,nqbz, &
     umnk)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-3)
  !      complex(8),allocatable :: amn(:,:),
  !     &                          evecc(:,:),smat(:,:),wmat(:,:)
  !      real (8),allocatable :: eval(:)
  complex(8),allocatable :: amn(:,:),cc(:,:),zz(:,:),vv(:,:)
  real(8),allocatable :: dd(:)
  complex(8) :: umnk(nwf,nwf,nqbz), &
       amnk(iko_ix:iko_fx,nwf,nqbz), &
       cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
  integer :: iko_i(nqbz),iko_f(nqbz)
  allocate (amn(nwf,nwf),zz(nwf,nwf),vv(nwf,nwf),dd(nwf))
  umnk = (0d0,0d0)
  do iq = 1,nqbz
     ! construct A
     amn = (0d0,0d0)
     do in = 1,nwf
        do im = 1,nwf
           do imp = iko_i(iq),iko_f(iq)
              amn(im,in) = amn(im,in) + &
                   dconjg(cnk(imp,im,iq))*amnk(imp,in,iq)
           enddo
        enddo
     enddo
     ! singular value decomposition: amn = zz*dd*vv
     call zgesvdmn(nwf,nwf,amn,dd,zz,vv)
     ! U(m,n) = ( A S^{-1/2} )_mn = (zz*vv)_mn
     do im = 1,nwf
        do in = 1,nwf
           do ik = 1,nwf
              umnk(im,in,iq) = umnk(im,in,iq) + zz(im,ik)*vv(ik,in)
           enddo
        enddo
     enddo

  enddo ! iq

  deallocate(amn,zz,vv,dd)

  return
end subroutine init_Umnk
!-----------------------------------------------------------------------
subroutine updt_mmn(umnk,mmn0,ikbidx, &
     nwf,nqbz,nbb, &
     mmn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps = 1d-3)
  complex(8),allocatable :: wmat(:,:),wmat2(:,:)
  complex(8) :: umnk(nwf,nwf,nqbz), &
       mmn(nwf,nwf,nbb,nqbz), &
       mmn0(nwf,nwf,nbb,nqbz)
  integer :: ikbidx(nbb,nqbz)

  allocate (wmat(nwf,nwf),wmat2(nwf,nwf))
  wmat = (0d0,0d0)

  do iq = 1,nqbz
     do ibb = 1,nbb
        iqb = ikbidx(ibb,iq)

        ! wmat = M * U
        wmat = (0d0,0d0)
        do ii = 1,nwf
           do in = 1,nwf
              do ij = 1,nwf
                 wmat(ii,in) = wmat(ii,in) + &
                      mmn0(ii,ij,ibb,iq)*umnk(ij,in,iqb)
              enddo
           enddo
        enddo

        ! wmat2 = U^{-1} * wmat
        wmat2 = (0d0,0d0)
        do im = 1,nwf
           do in = 1,nwf
              do ii = 1,nwf
                 wmat2(im,in) = wmat2(im,in) + &
                      dconjg(umnk(ii,im,iq)) * wmat(ii,in)
              enddo
           enddo
        enddo

        ! M = wmat2
        mmn(:,:,ibb,iq) = wmat2(:,:)

        ! end of ibb-loop
     enddo
     ! end of iq-loop
  enddo

  deallocate (wmat,wmat2)

  return
end subroutine updt_mmn
!-----------------------------------------------------------------------
subroutine get_rn(mmn,bb,wbb,wbz, &
     nwf,nqbz,nbb, &
     rn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: mmn(nwf,nwf,nbb,nqbz)
  real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf)

  rn = 0d0

  do in = 1,nwf
     do ibb = 1,nbb
        rtmp = 0d0
        do iq = 1,nqbz
           rtmp = rtmp &
                + dimag(log(mmn(in,in,ibb,iq)))*wbz(iq)
           ! ccccccccccccccccccccccccccc
           !            if(abs(log(mmn(in,in,ibb,iq)))>0.9) then
           !               write(6,"('mmmnnn in ibb iq log(mmn)',3i5,2f13.4)")
           !     &         in,ibb,iq,log(mmn(in,in,ibb,iq))
           !            endif
           ! ccccccccccccccccccccccccccc
        enddo
        rn(:,in) = rn(:,in) + wbb(ibb)*bb(:,ibb)*rtmp
        !            write(*,*)in,ibb,rtmp
     enddo
  enddo

  rn(:,:) = - rn(:,:)

  return
end subroutine get_rn
!--------------------------------------------------------r---------------
subroutine get_rnm(mmn,bb,wbb,wbz, &
     nwf,nqbz,nbb, &
     rnm)! complex. See Eq.(4) and Eq.(22) in [1].
  implicit none
  integer:: in,im,ibb,iq,nbb,nqbz,nwf
  complex(8) :: mmn(nwf,nwf,nbb,nqbz)
  real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz)
  complex(8):: img=(0d0,1d0),rnm(3,nwf,nwf),rtmp
  rnm = 0d0
  do in = 1,nwf
     do im = 1,nwf
        if(im==in) cycle
        do ibb = 1,nbb
           rtmp = 0d0
           do iq = 1,nqbz
              rtmp = rtmp + img*mmn(in,im,ibb,iq)*wbz(iq)
           enddo
           rnm(:,in,im) = rnm(:,in,im) + wbb(ibb)*bb(:,ibb)*rtmp
        enddo
     enddo
  enddo
end subroutine get_rnm
!-----------------------------------------------------------------------
subroutine getrmn(mmn, &
     nwf, &
     rmn)

  ! [1] eq.(45)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  complex(8) :: mmn(nwf,nwf),rmn(nwf,nwf)

  rmn = (0d0,0d0)

  do im = 1,nwf
     do in = 1,nwf
        rmn(im,in) = mmn(im,in) * dconjg(mmn(in,in))
     enddo
  enddo

  return
end subroutine getrmn
!-----------------------------------------------------------------------
subroutine getamn(bmn, &
     nwf, &
     amn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: bmn(nwf,nwf),amn(nwf,nwf)

  amn = (0d0,0d0)

  do im = 1,nwf
     do in = 1,nwf
        amn(im,in) = (bmn(im,in) - dconjg(bmn(in,im))) * 0.5d0
     enddo
  enddo

  return
end subroutine getamn
!-----------------------------------------------------------------------
subroutine getsmn(bmn, &
     nwf, &
     smn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: bmn(nwf,nwf),smn(nwf,nwf),ci

  smn = (0d0,0d0)
  ci = (0d0,1d0)

  do im = 1,nwf
     do in = 1,nwf
        smn(im,in) = -(bmn(im,in) + dconjg(bmn(in,im))) * ci * 0.5d0
     enddo
  enddo

  return
end subroutine getsmn
!-----------------------------------------------------------------------
subroutine gettmn(rn,mmn,bb, &
     nwf, &
     qn,tmn)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  complex(8),allocatable :: rtmn(:,:)
  complex(8) :: mmn(nwf,nwf),tmn(nwf,nwf)
  !      real(8),allocatable :: qn(:)
  real(8) :: qn(nwf)
  real(8) :: rn(3,nwf),bb(3)

  allocate (rtmn(nwf,nwf))

  ! qn ([1] eq.47)
  do in = 1,nwf
     qn(in) = dimag(log(mmn(in,in))) + sum(bb(:)*rn(:,in))
  enddo

  ! R~mn ([1] eq.48)
  do im = 1,nwf
     do in = 1,nwf
        rtmn(im,in) = mmn(im,in) / mmn(in,in)
     enddo
  enddo

  ! Tmn = R~mn * qn
  do im = 1,nwf
     do in = 1,nwf
        tmn(im,in) = rtmn(im,in) * qn(in)
     enddo
  enddo

  deallocate(rtmn)

  return
end subroutine gettmn
!-----------------------------------------------------------------------
subroutine updt_uk(wmn, &
     nwf, &
     umn)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8),allocatable :: iwmn(:,:),ewmn(:,:),evecc(:,:), &
       zmn(:,:),evalc(:)
  complex(8) :: wmn(nwf,nwf),umn(nwf,nwf),ci
  real(8),allocatable :: eval(:)

  allocate (iwmn(nwf,nwf),ewmn(nwf,nwf), &
       evecc(nwf,nwf),zmn(nwf,nwf),evalc(nwf), &
       eval(nwf))

  ci = (0d0,1d0)

  ! diagonalize W
  ! notice! W is anti-Hermitian, not Hermitian
  iwmn(:,:) = ci * wmn(:,:)

  call chk_hm(iwmn,nwf)
  call diag_hm(iwmn,nwf,eval,evecc)

  ! eW = exp(W) = exp(-i * iW)
  do ii = 1, nwf
     !         write(*,*)'ev',ii,eval(ii)
     evalc(ii) = exp(-ci*eval(ii))
  enddo

  ewmn = (0d0,0d0)
  do im = 1,nwf
     do in = 1,nwf
        do il = 1,nwf
           ewmn(im,in) = ewmn(im,in) + &
                evecc(im,il) * evalc(il) * dconjg(evecc(in,il))
        enddo
     enddo
  enddo

  call chk_um(ewmn,nwf)

  ! U = U * exp(W)
  zmn = (0d0,0d0)
  do im = 1,nwf
     do in = 1,nwf
        do il = 1,nwf
           zmn(im,in) = zmn(im,in) + umn(im,il) * ewmn(il,in)
        enddo
     enddo
  enddo
  umn = zmn

  deallocate(iwmn,ewmn,evecc,zmn,evalc,eval)

  return
end subroutine updt_uk
!-----------------------------------------------------------------------
subroutine getOmg(mmn,rn,bb,wbb,wbz, &
     nwf,nqbz,nbb, &
     omgi,omgd,omgod,omgdod,omgidod, &
     alat,iprint)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: mmn(nwf,nwf,nbb,nqbz)
  real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf),sss,r2n,alat,pi=4d0*atan(1d0)
  logical :: iprint
  ! omgi ([1] eq.34)
  omgi = 0d0
  do iq = 1,nqbz
     do ibb = 1,nbb
        a2mmn = 0d0
        do im = 1,nwf
           do in = 1,nwf
              a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
        enddo
        omgi = omgi + wbz(iq)*wbb(ibb)*(dble(nwf)-a2mmn)
     enddo
  enddo

  ! omgod ([1] eq.35)
  omgod = 0d0
  do iq = 1,nqbz
     do ibb = 1,nbb
        a2mmn = 0d0
        do im = 1,nwf
           do in = 1,nwf
              if (im /= in) &
                   a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
        enddo
        omgod = omgod + wbz(iq)*wbb(ibb)*a2mmn
     enddo
  enddo

  ! omgd ([1] eq.36)
  omgd = 0d0
  do ibb = 1,nbb
     rtmp = 0d0
     do in = 1,nwf
        brn = sum(bb(:,ibb)*rn(:,in))
        !            write(6,*)'brn=',brn,sum(abs(bb(:,ibb))),sum(abs(rn(:,in)))
        do iq = 1,nqbz
           rtmp = rtmp + &
                wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
        enddo
     enddo
     !         sss=0d0
     !         do in=1,nwf
     !         sss=sss+sum(abs(imag(mmn(in,in,ibb,1:nqbz))))
     !         enddo
     !         write(6,*)'sss=',sss
     omgd = omgd + wbb(ibb)*rtmp
  enddo
  !     !
  if(iprint) then
     do in = 1,nwf
        r2n=0d0
        do ibb = 1,nbb
           brn = sum(bb(:,ibb)*rn(:,in))
           do iq = 1,nqbz
              r2n = r2n + &
                   wbz(iq)*wbb(ibb)* ( (1d0-abs(mmn(in,in,ibb,iq))**2) + (dimag(log(mmn(in,in,ibb,iq))))**2 )
           enddo
        enddo
        write(6,"('dddd spread: iwf r2n(bohr^2) rn(bohr)=',i3,f13.6,3x,3(f9.4))") &
             in,(alat/(2d0*pi))**2*(r2n-sum(rn(:,in)**2)), alat/(2d0*pi)*rn(:,in)
     enddo
  endif

  ! sum up
  omgdod = omgd + omgod
  omgidod = omgi + omgdod

  return
end subroutine getOmg
!-----------------------------------------------------------------------
subroutine writeOmg(is,mmn,rn,bb,wbb,wbz,tpia, &
     nwf,nqbz,nbb)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: mmn(nwf,nwf,nbb,nqbz)
  real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf),omg(nwf)


  ! omgi ([1] eq.34)
  omgi = 0d0
  do iq = 1,nqbz
     do ibb = 1,nbb
        a2mmn = 0d0
        do im = 1,nwf
           do in = 1,nwf
              a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
        enddo
        omgi = omgi + wbz(iq)*wbb(ibb)*(dble(nwf)-a2mmn)
     enddo
  enddo
  print *,'ooooooooooomgi=', omgi
  omg = 0d0
  do im = 1,nwf
     do iq = 1,nqbz
        do ibb = 1,nbb
           a2mmn = 0d0
           do in = 1,nwf
              a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
           omg(im) = omg(im) + wbz(iq)*wbb(ibb)*(1d0-a2mmn)
        enddo
     enddo
  enddo

  ! omgod ([1] eq.35)
  omgod = 0d0
  do iq = 1,nqbz
     do ibb = 1,nbb
        a2mmn = 0d0
        do im = 1,nwf
           do in = 1,nwf
              if (im /= in) &
                   a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
        enddo
        omgod = omgod + wbz(iq)*wbb(ibb)*a2mmn
     enddo
  enddo

  do im = 1,nwf
     do iq = 1,nqbz
        do ibb = 1,nbb
           a2mmn = 0d0
           do in = 1,nwf
              if (im /= in) &
                   a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
           enddo
           omg(im) = omg(im) + wbz(iq)*wbb(ibb)*a2mmn
        enddo
     enddo
  enddo

  ! omgd ([1] eq.36)
  omgd = 0d0
  do ibb = 1,nbb
     rtmp = 0d0
     do in = 1,nwf
        brn = sum(bb(:,ibb)*rn(:,in))
        do iq = 1,nqbz
           rtmp = rtmp + &
                wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
        enddo
     enddo
     omgd = omgd + wbb(ibb)*rtmp
  enddo

  do in = 1,nwf
     do ibb = 1,nbb
        rtmp = 0d0
        brn = sum(bb(:,ibb)*rn(:,in))
        do iq = 1,nqbz
           rtmp = rtmp + &
                wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
        enddo
        omg(in) = omg(in) + wbb(ibb)*rtmp
     enddo
  enddo

  ! multiply (a/2pi)^2
  omgd = omgd / tpia**2
  omgod = omgod / tpia**2
  omgi = omgi / tpia**2
!  omgdod = omgdod / tpia**2
  omg = omg / tpia**2

  ! sum up
  omgdod = omgd + omgod
  omgidod = omgi + omgdod

  ! output
  if (is == 1) then
     ifbnd = iopen('spread.up',1,-1,0)
  else
     ifbnd = iopen('spread.dn',1,-1,0)
  endif

  write(ifbnd,*)'Total   omega(a.u.)     omega(ang.^2)'
  write(ifbnd,"(5x,2f16.8)")omgidod,omgidod*0.5292d0**2
  write(ifbnd,*)'band    omega(a.u.)     omega(ang.^2)'
  do in = 1,nwf
     write(ifbnd,"(i5,2f16.8)")in,omg(in),omg(in)*0.5292d0**2
  enddo
  omgall = sum(omg(:))
  !      write(ifbnd,"(2x,a3,2f16.8)")'all',omgall,omgall*0.5292d0**2

  write(ifbnd,*)'(a/2pi)^2 corrected'
  if (is == 1) then
     isx = iclose('spread.up')
  else
     isx = iclose('spread.dn')
  endif

  return
end subroutine writeOmg
!-----------------------------------------------------------------------
subroutine writermn(is,mmn,bb,wbb,qbz,qbz0,wbz,rt, &
     nwf,nqbz,nbb,n1,n2,n3)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr &
       ,ctmp1,ctmp2 &
       ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz) &
       ,c1,c2,c3
  real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),qbz0(3,nqbz),wbz(nqbz), &
       rt(3,nqbz),dmn,rtmp(3)

  pi = 4d0*atan(1d0)
  ci = (0d0,1d0)
  czero = (0d0,0d0)

  ! Berry connecrion (in the Wannier gauge)
  amn = czero
  do im = 1,nwf
     do in = 1,nwf
        dmn = 0d0
        if (im == in) dmn = 1d0
        do iq = 1,nqbz
           do ibb = 1,nbb
              !            ctmp1 = (mmn(im,in,ibb,iq) - dmn)*wbb(ibb)
              r1 = dimag(0.5d0*(mmn(im,in,ibb,iq)+mmn(in,im,ibb,iq))-dmn)
              r2 = dreal(0.5d0*(mmn(im,in,ibb,iq)-mmn(in,im,ibb,iq)))
              ctmp1 = dcmplx(r2,r1)
              ctmp1 = log(ctmp1+1d0)*wbb(ibb)
              do ix = 1,3 ! x,y,z
                 amn(ix,im,in,iq) = amn(ix,im,in,iq) &
                      + ci*imag(ctmp1) * bb(ix,ibb)
              enddo ! ix
           enddo ! ibb
        enddo ! iq
     enddo ! in
  enddo ! im
  amn = amn * ci

  ! <0m | r | Rn>
  rmn = czero
  do ir = 1,nqbz
     do iq = 1,nqbz
        !         rk = sum(rt(:,ir)*qbz(:,iq))
        rk = sum(rt(:,ir)*qbz0(:,iq))
        cikr = -ci * 2d0 * pi * rk
        ceikr = exp(cikr) * wbz(iq)
        do in = 1,nwf
           do im = 1,nwf
              do ix = 1,3
                 rmn(ix,im,in,ir) = rmn(ix,im,in,ir) + &
                      ceikr * amn(ix,im,in,iq)
              enddo ! ix
           enddo ! im
        enddo ! in
     enddo ! iq
  enddo ! ir

  ! output
  if (is == 1) then
     ifrmn = iopen('rmn.up',1,-1,0)
  else
     ifrmn = iopen('rmn.dn',1,-1,0)
  endif

  write(ifrmn,*)'*** nwf,nsite'
  write(ifrmn,*)nwf,nqbz
  write(ifrmn,*)'*** rsite'
  write(ifrmn,*)rt
  write(ifrmn,*)'*** rmn'
  write(ifrmn,*)rmn
  write(ifrmn,*)'*** amn'
  write(ifrmn,*)amn

  if (is == 1) then
     isx = iclose('rmn.up')
  else
     isx = iclose('rmn.dn')
  endif

  return
end subroutine writermn
!-----------------------------------------------------------------------
subroutine writemmn(is,mmn,bb,wbb,qbz,wbz,rt, &
     nwf,nqbz,nbb,n1,n2,n3)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr &
       ,ctmp1,ctmp2 &
       ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz) &
       ,c1,c2,c3
  real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),wbz(nqbz), &
       rt(3,nqbz),dmn

  pi = 4d0*atan(1d0)
  ci = (0d0,1d0)
  czero = (0d0,0d0)

  ! output
  if (is == 1) then
     ifrmn = iopen('mmn.up',1,-1,0)
  else
     ifrmn = iopen('mmn.dn',1,-1,0)
  endif

  write(ifrmn,*)'*** nwf,nsite,nb'
  write(ifrmn,*)nwf,nqbz,nbb
  write(ifrmn,*)'*** mmn'
  write(ifrmn,*)mmn
  write(ifrmn,*)'*** bb'
  write(ifrmn,*)bb
  write(ifrmn,*)'*** wbb'
  write(ifrmn,*)wbb
  write(ifrmn,*)'*** qbz'
  write(ifrmn,*)qbz
  write(ifrmn,*)'*** wbz'
  write(ifrmn,*)wbz

  if (is == 1) then
     isx = iclose('mmn.up')
  else
     isx = iclose('mmn.dn')
  endif

  return
end subroutine writemmn
!-----------------------------------------------------------------------
subroutine wmaxloc(ifmlw,ifmlwe, &
     qbz,umnk,cnk,eunk, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nwf,nqbz,nband,nlmto,isp, dnk)
  use m_ftox
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8):: dnk(iko_ix:iko_fx,nwf,nqbz)
  !      complex(8),allocatable :: dnk(:,:,:)
  complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz), &
       umnk(nwf,nwf,nqbz)
  real(8) :: qbz(3,nqbz),eunk(nwf,nqbz)
  integer(4) :: iko_i(nqbz),iko_f(nqbz),ib,isp
  character*4::fname

  !      allocate(dnk(iko_ix:iko_fx,nwf,nqbz))

  write(ifmlw)nqbz,nwf,iko_ix,iko_fx
  !      write(ifmlwe)nqbz,nwf,iko_ix,iko_fx

  !$$$      if (isp.eq.1) fname='pkmU'
  !$$$      if (isp.eq.2) fname='pkmD'
  !$$$      ifi=9034
  !$$$      open(ifi,file=fname,form='formatted',status='unknown')
  !$$$      write(ifi,"('== p_km^alpha in PRB83,121101 ! weight in l-subspace ==')")
  !$$$      write(ifi,"(' (this is called as c^sigma_km in book of 45th IFFK by Ersoy)')")
  !$$$      write(ifi,"(8i8)") nqbz,nwf,iko_ix,iko_fx
  !$$$      write(ifi,"('       |pkm|**2          ib      iq     q(1:3)')")
  dnk = (0d0,0d0)
  do iq = 1,nqbz

     do imp = iko_i(iq),iko_f(iq)
        do in = 1,nwf
           do im = 1,nwf
              dnk(imp,in,iq) = dnk(imp,in,iq) &
                   + umnk(im,in,iq) * cnk(imp,im,iq)
           enddo
        enddo
     enddo

     !         do iwf=1,nwf
     !            write(6,ftox)' wwww   dnk',isp,iq,ftof(dnk(iko_ix:iko_fx,iwf,iq),2)
     !         enddo

     write(ifmlw)iq,qbz(1:3,iq)
     write(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iq)

     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !! jan2015 new section. Write out weight in l-subsapce.
     !         do iwf=1,nwf
     !            do iko=iko_ix,iko_fx
     !               write(*,"('ib iwf iq=',3i4,f10.5)")iko,iwf,iq,abs(dnk(iko,iwf,iq))**2
     !            enddo
     !            write(*,"('iwf iq sumweight=',2i4,f10.5)")
     !     &        iwf,iq, sum(abs(dnk(iko_ix:iko_fx,iwf,iq))**2)
     !         enddo
     !         do ib = iko_ix,iko_fx
     !            write(*,"('ib iq winsubspace=',2i4,f10.5)")
     !     &      ib,iq, sum(abs(dnk(ib,1:nwf,iq))**2)
     !     enddo

     !         write(*,"(' iq sumcheck winsub=',i4,f10.5)")
     !     &      iq, sum(abs(dnk(iko_ix:iko_fx,1:nwf,iq))**2)

     !$$$         do ib = iko_ix,iko_fx
     !$$$            write(ifi,"(f19.15, 2i8, 3f13.6 )")
     !$$$     &      sum(abs(dnk(ib,1:nwf,iq))**2),ib, iq, qbz(1:3,iq)
     !$$$         enddo

     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !         write(ifmlwe)iq,qbz(1:3,iq)
     !         write(ifmlwe)eunk(1:nwf,iq)

  enddo
  !$$$      close(ifi)

  !      call chk_cnkweight(qbz,iko_ix,iko_fx,dnk,
  !     &     nqbz,nwf,nband,nlmto)

  !      deallocate(dnk)

  return
end subroutine wmaxloc
!-----------------------------------------------------------------------
subroutine chk_dnk(is,eunk,qbz, &
     umnk,cnk, &
     iko_ix,iko_fx,iko_i,iko_f, &
     nband,nwf,nqbz)

  use m_readeigen,only:readeval
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  parameter (eps=1d-3)
  complex(8),allocatable :: dnk(:,:)
  complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz), &
       umnk(nwf,nwf,nqbz), &
       ctmp,ham(nwf,nwf),evecc(nwf,nwf)
  real(8) :: eunk(nwf,nqbz),qbz(3,nqbz),eks(nband),eval(nwf)
  integer(4) :: iko_i(nqbz),iko_f(nqbz)

  allocate(dnk(iko_ix:iko_fx,nwf))

  do iq = 1,nqbz
     eks= readeval (qbz(:,iq),is)
     dnk = (0d0,0d0)
     do imp = iko_i(iq),iko_f(iq)
        do in = 1,nwf
           do im = 1,nwf
              dnk(imp,in) = dnk(imp,in) &
                   + umnk(im,in,iq) * cnk(imp,im,iq)
           enddo
        enddo
     enddo


     ham = 0d0
     do ii=1,nwf
        do ij=1,nwf
           ctmp = 0d0
           do ik = iko_i(iq),iko_f(iq)
              ctmp = ctmp + eks(ik) * &
                   dconjg(dnk(ik,ii))*dnk(ik,ij)
           enddo
           ham(ii,ij) = ctmp
           write(97,"(3i5,2f12.6)")iq,ii,ij,dreal(ctmp),dimag(ctmp)
        enddo
     enddo
     call diag_hm(ham,nwf,eval,evecc)
     do ii = 1,nwf
        write(98,"(2i5,2f12.6)")iq,ii,eval(ii),eunk(ii,iq)
     enddo

  enddo

  deallocate(dnk)

  return
end subroutine chk_dnk
!-----------------------------------------------------------------------
