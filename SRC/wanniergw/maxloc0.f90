!$$$      subroutine wan_input(leout,lein,lbin,ieo_swt,iei_swt,
!$$$     &    eomin,eomax,itout_i,itout_f,nbbelow,nbabove,
!$$$     &    eimin,eimax,itin_i,itin_f,
!$$$     &    nsc1,nsc2,conv1,conv2,alpha1,alpha2,rcut)
!$$$      use m_keyvalue,only: getkeyvalue
!$$$
!$$$      implicit none
!$$$      real(8) :: eomin,eomax,eimin,eimax,conv1,conv2,alpha1,alpha2,rcut
!$$$      real(8) :: alph1,alp2,rc1
!$$$      integer(4) :: ieo_swt,iei_swt,itout_i,itout_f,itin_i,itin_f,
!$$$     &    nbbelow,nbabove,nsc1,nsc2
!$$$      logical :: leout,lein,lbin
!$$$!
!$$$      ieo_swt = 0
!$$$      eomin   = 0d0
!$$$      eomax   = 0d0
!$$$      itout_i = 0
!$$$      itout_f = 0
!$$$      iei_swt = 0
!$$$      eimin   = 0d0
!$$$      eimax   = 0d0
!$$$      itin_i  = 0
!$$$      itin_f  = 0
!$$$      nbbelow = 0
!$$$      nbabove = 0
!$$$      call getkeyvalue("GWinput","wan_out_ewin",leout,default=.true.)
!$$$      call getkeyvalue("GWinput","wan_in_ewin",lein,default=.false.)
!$$$      call getkeyvalue("GWinput","wan_in_bwin",lbin,default=.false.)
!$$$      if (leout) then
!$$$        call getkeyvalue("GWinput","wan_out_emin",eomin,default=999d0 )
!$$$        call getkeyvalue("GWinput","wan_out_emax",eomax,default=-999d0 )
!$$$        if (eomin.gt.eomax) stop 'hmaxloc: eomin > eomax'
!$$$        ieo_swt = 1
!$$$      else
!$$$        call getkeyvalue("GWinput","wan_out_bmin",itout_i,default=999 )
!$$$        call getkeyvalue("GWinput","wan_out_bmax",itout_f,default=-999 )
!$$$        if (itout_i.gt.itout_f) stop 'hmaxloc: itout_i > itout_f'
!$$$      endif
!$$$      if (lein) then
!$$$        call getkeyvalue("GWinput","wan_in_emin",eimin,default=999d0 )
!$$$        call getkeyvalue("GWinput","wan_in_emax",eimax,default=-999d0 )
!$$$        if (eimin.gt.eimax) stop 'hmaxloc: eimin > eimax'
!$$$        iei_swt = 1
!$$$      endif
!$$$      if (lbin) then
!$$$        call getkeyvalue("GWinput","wan_in_bmin",itin_i,default=999 )
!$$$        call getkeyvalue("GWinput","wan_in_bmax",itin_f,default=-999 )
!$$$        if (itin_i.gt.itin_f) stop 'hmaxloc: itin_i > itin_f'
!$$$        iei_swt = 2
!$$$      endif
!$$$      call getkeyvalue("GWinput","wan_maxit_1st",nsc1,default=100)
!$$$      call getkeyvalue("GWinput","wan_conv_1st",conv1,default=1d-5)
!$$$      call getkeyvalue("GWinput","wan_mix_1st",alpha1,default=0.1d0)
!$$$      call getkeyvalue("GWinput","wan_maxit_2nd",nsc2,default=100)
!$$$      call getkeyvalue("GWinput","wan_conv_2nd",conv2,default=1d-5)
!$$$      call getkeyvalue("GWinput","wan_mix_2nd",alpha2,default=0.1d0)
!$$$      call getkeyvalue("GWinput","wan_tb_cut",rcut,default=1.01d0)
!$$$      call getkeyvalue("GWinput","wan_nb_below",nbbelow,default=0)
!$$$      call getkeyvalue("GWinput","wan_nb_above",nbabove,default=0)
!$$$!
!$$$      return
!$$$      end
!-----------------------------------------------------------------------
subroutine getrt(qbz,qbas,plat,n1,n2,n3,nqbz, &
     rt,rt8,qbz0)
  implicit integer (i-n)

  implicit real*8(a-h,o-z)

  parameter (eps=1d-4)
  real(8) :: qbz(3,nqbz),qbas(3,3),plat(3,3),r1(3),r2(3), &
       rt(3,nqbz),rt8(3,8,nqbz),qbz0(3,nqbz)

  ! rt
  irt = 0
  do i1 = 0,n1-1
     a1 = dble(i1)
     do i2 = 0,n2-1
        a2 = dble(i2)
        do i3 = 0,n3-1
           a3 = dble(i3)

           irt = irt + 1
           rt(1,irt) = a1
           rt(2,irt) = a2
           rt(3,irt) = a3
        enddo
     enddo
  enddo

  if (irt /= nqbz) stop 'getrt: nqbz error'

  ! rt8
  irt = 0
  do i1 = 0,n1-1
     j1 = i1-n1
     a1 = dble(i1)
     b1 = dble(j1)
     if (abs(j1) < abs(i1)) then
        a1 = b1
     elseif (abs(i1) < abs(j1)) then
        b1 = a1
     elseif (i1+j1 /= 0) then
        stop 'getrt: i1 error'
     endif
     do i2 = 0,n2-1
        j2 = i2-n2
        a2 = dble(i2)
        b2 = dble(j2)
        if (abs(j2) < abs(i2)) then
           a2 = b2
        elseif (abs(i2) < abs(j2)) then
           b2 = a2
        elseif (i2+j2 /= 0) then
           stop 'getrt: i2 error'
        endif
        do i3 = 0,n3-1
           j3 = i3-n3
           a3 = dble(i3)
           b3 = dble(j3)
           if (abs(j3) < abs(i3)) then
              a3 = b3
           elseif (abs(i3) < abs(j3)) then
              b3 = a3
           elseif (i3+j3 /= 0) then
              stop 'getrt: i3 error'
           endif

           irt = irt + 1

           rt8(1,1,irt) = a1
           rt8(1,2,irt) = b1
           rt8(1,3,irt) = a1
           rt8(1,4,irt) = b1
           rt8(1,5,irt) = a1
           rt8(1,6,irt) = b1
           rt8(1,7,irt) = a1
           rt8(1,8,irt) = b1

           rt8(2,1,irt) = a2
           rt8(2,2,irt) = a2
           rt8(2,3,irt) = b2
           rt8(2,4,irt) = b2
           rt8(2,5,irt) = a2
           rt8(2,6,irt) = a2
           rt8(2,7,irt) = b2
           rt8(2,8,irt) = b2

           rt8(3,1,irt) = a3
           rt8(3,2,irt) = a3
           rt8(3,3,irt) = a3
           rt8(3,4,irt) = a3
           rt8(3,5,irt) = b3
           rt8(3,6,irt) = b3
           rt8(3,7,irt) = b3
           rt8(3,8,irt) = b3
        enddo
     enddo
  enddo

  if (irt /= nqbz) stop 'getrt: nqbz error'

  ! check qbas(i) * plat(j) = delta(i,j)
  do ii=1,3
     do ij=1,3
        rtmp = sum(qbas(:,ii)*plat(:,ij))
        if (ii == ij) rtmp = rtmp - 1d0
        if (abs(rtmp) > eps) stop 'getrt: qbas*plat error'
     enddo
  enddo

  ! calc qbz0
  do iq = 1,nqbz
     call q2q0(qbz(:,iq),plat,qbz0(:,iq))
  enddo

  return
end subroutine getrt
!-----------------------------------------------------------------------
subroutine q2q0(q,plat,q0)
  implicit integer (i-n)

  implicit real*8(a-h,o-z)
  parameter (eps=1d-4)
  real(8) :: q(3),q0(3),plat(3,3)

  do ii = 1,3
     qp = sum(q(:)*plat(:,ii))
     qp = qp - dble(int(qp)) + 1d0 + eps
     qp = qp - dble(int(qp)) - eps
     if (qp > 1d0 .OR. qp < -eps) stop 'q2q0: qp error'
     q0(ii) = qp
  enddo

  return
end subroutine q2q0
!-----------------------------------------------------------------------
subroutine  getbb(plat,alat,n1,n2,n3, nbb,wbb,wbbsum,bb)
  ! determine bb matrix
  ! finite difference method for k-space grids
  ! Ref.
  ! Appendix B, Marzari and Vanderbilt, PRB56, 12847 (1997)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-4)
  real(8),allocatable :: glat(:,:)
  real(8) :: plat(3,3),qlat(3,3),wbb(12),bb(3,12), &
       vtmp(3,7),avtmp(7)
  logical :: lcub,lort,labc
  ! for inverting matrix
  integer(4):: ndim,nmx,nev
  complex(8),allocatable :: cmat(:,:)
  real(8),   allocatable :: rmat(:,:),work(:)
  ! begin
  call minv33tp(plat,qlat) ! inverse and tranpose
  ! eck write
  !      write(*,*)'n1,n2,n3'
  !      write(*,*)n1,n2,n3
  !      write(*,*)'plat'
  !      write(*,990)plat(:,1)
  !      write(*,990)plat(:,2)
  !      write(*,990)plat(:,3)
  !      write(*,*)'qlat'
  !      write(*,990)qlat(:,1)
  !      write(*,990)qlat(:,2)
  !      write(*,990)qlat(:,3)
  ! check if 1D or not
  if (n1 == 1 .AND. n2 == 1) then
     write(*,*)'getbb: 1D'
     nbb = 2
     bb(1:3,1) = qlat(1:3,1) / n3
     bb(1:3,2) = -bb(1:3,1)
     db = dsqrt(sum(qlat(1:3,3)**2)) / dble(n3)
     wbb(:) = 1d0 / (2d0*db*db)

     wbbsum = 0d0
     do ibb = 1,nbb
        wbbsum = wbbsum + wbb(ibb)
     enddo

     ! check if S[i=1,2] wbb(i)*bb(3,i)*bb(3,i) = 1
     allocate (work(1))
     work = 0d0
     do i=1,nbb
        work(1) = work(1) + bb(3,i)*bb(3,i)*wbb(i)
     enddo

     ! eck write
     tmp = dabs(work(1)-1d0)
     if (tmp > eps) stop 'getbb: wbb is wrong'

     deallocate (work)

     return
  endif

  ! check orthogonal, a=b=c or not
  dq1 = dsqrt(sum(qlat(1:3,1)**2))
  dq2 = dsqrt(sum(qlat(1:3,2)**2))
  dq3 = dsqrt(sum(qlat(1:3,3)**2))
  d12 = dabs(dq1 - dq2)
  d13 = dabs(dq1 - dq3)

  q1q2 = sum(qlat(1:3,1)*qlat(1:3,2)) / (dq1*dq2)
  q1q3 = sum(qlat(1:3,1)*qlat(1:3,3)) / (dq1*dq3)
  q2q3 = sum(qlat(1:3,2)*qlat(1:3,3)) / (dq2*dq3)

  lort = .false.
  if (dabs(q1q2) < eps .AND. &
       dabs(q1q3) < eps .AND. &
       dabs(q2q3) < eps) then
     lort = .true.
  endif

  labc = .false.
  if (d12 < eps .AND. d13 < eps .AND. &
       n1 == n2 .AND. n1 == n3) then
     labc = .true.
  endif

  ! if a=b=c
  if (labc) then

     ! simple cubic
     if (lort) then

        write(*,*)'getbb: simple cubic'
        nbb = 6
        bb(1:3,1) = qlat(1:3,1) / n1
        bb(1:3,2) = qlat(1:3,2) / n2
        bb(1:3,3) = qlat(1:3,3) / n3
        do i = 1,3
           bb(1:3,i+3) = -bb(1:3,i)
        enddo
        wbb(:) = 3d0 / (nbb*dq1*dq1) * n1*n1

        goto 900
     endif

     ! fcc (bcc in k-space)
     tmp12 = dabs(q1q2 + 1d0/3d0)
     tmp13 = dabs(q1q3 + 1d0/3d0)
     tmp23 = dabs(q2q3 + 1d0/3d0)
     if (tmp12 < eps .AND. &
          tmp13 < eps .AND. &
          tmp23 < eps) then
        write(*,*)'getbb: fcc'
        nbb = 8
        bb(1:3,1) = qlat(1:3,1) / n1
        bb(1:3,2) = qlat(1:3,2) / n2
        bb(1:3,3) = qlat(1:3,3) / n3
        bb(1:3,4) = - bb(1:3,1) - bb(1:3,2) - bb(1:3,3)
        do i = 1,4
           bb(1:3,i+4) = -bb(1:3,i)
        enddo
        wbb(:) = 3d0 / (nbb*dq1*dq1) * n1*n1 !3/(Z*|b|*|b|) See Marzari (B1) around

        goto 900
     endif

     ! bcc (fcc in k-space)
     tmp12 = dabs(q1q2 - 0.5d0)
     tmp13 = dabs(q1q3 - 0.5d0)
     tmp23 = dabs(q2q3 - 0.5d0)
     if (tmp12 < eps .AND. &
          tmp13 < eps .AND. &
          tmp23 < eps) then
        write(*,*)'getbb: bcc'
        nbb = 12
        bb(1:3,1) = qlat(1:3,1) / n1
        bb(1:3,2) = qlat(1:3,2) / n2
        bb(1:3,3) = qlat(1:3,3) / n3
        bb(1:3,4) = bb(1:3,1) - bb(1:3,2)
        bb(1:3,5) = bb(1:3,2) - bb(1:3,3)
        bb(1:3,6) = bb(1:3,3) - bb(1:3,1)
        do i = 1,6
           bb(1:3,i+6) = -bb(1:3,i)
        enddo
        wbb(:) = 3d0 / (nbb*dq1*dq1) * n1*n1

        goto 900
     endif

     ! end of if (labc)
  endif

  ! hexagonal
  if ((dq1 == dq2) .AND. (abs(abs(q1q2)-0.5d0) < eps) .AND. &
       (abs(q1q3) < eps) .AND. (abs(q2q3) < eps)) then
     write(*,*)'getbb: hexagonal'
     nbb = 8
     bb(1:3,1) = qlat(1:3,1) / n1
     bb(1:3,2) = qlat(1:3,2) / n2
     if (q1q2 > 0.0d0) then
        bb(1:3,3) = bb(1:3,1) - bb(1:3,2)
     else
        bb(1:3,3) = bb(1:3,1) + bb(1:3,2)
     endif
     bb(1:3,4) = qlat(1:3,3) / n3

     bb(1:3,5) = -bb(1:3,1)
     bb(1:3,6) = -bb(1:3,2)
     bb(1:3,7) = -bb(1:3,3)
     bb(1:3,8) = -bb(1:3,4)
     wbb(:) = 1d0 / (3d0*dq1*dq1) * n1*n1
     wbb(4) = 1d0 / (2d0*dq3*dq3) * n3*n3
     wbb(8) = 1d0 / (2d0*dq3*dq3) * n3*n3

     goto 900
  endif

  ! orthrombic, but not cubic
  if (lort) then

     write(*,*)'getbb: orthorhombic'
     nbb = 6
     bb(1:3,1) = qlat(1:3,1) / n1
     bb(1:3,2) = qlat(1:3,2) / n2
     bb(1:3,3) = qlat(1:3,3) / n3
     wbb(1) = 0.5d0 / (dq1*dq1) * n1*n1
     wbb(2) = 0.5d0 / (dq2*dq2) * n2*n2
     wbb(3) = 0.5d0 / (dq3*dq3) * n3*n3

     do i = 1,3
        bb(1:3,i+3) = -bb(1:3,i)
        wbb(   i+3) = wbb(i)
     enddo

     goto 900

     ! non-orthorhombic cell
  else

     write(*,*)'getbb: non-orthorhombic'
     nbb = 12

     ! m 041218
     !         bb(1:3,1) = qlat(1:3,1) / n1
     !         bb(1:3,2) = qlat(1:3,2) / n2
     !         bb(1:3,3) = qlat(1:3,3) / n3
     nshell = 2
     ntmp = (2*nshell + 1)**3
     allocate (glat(3,ntmp))
     call getglat(n1,n2,n3,nshell,ntmp,qlat,glat)
     !         write(*,*) 'getglat done'
     call sortvec(ntmp,glat)
     !         write(*,*) 'sortvec done'
     call get3vec(ntmp,glat,bb(1:3,1:3))
     deallocate (glat)

     !         write(*,*) 'get3vec done'
     !         write(*,500)bb(1:3,1)
     !         write(*,500)bb(1:3,2)
     !         write(*,500)bb(1:3,3)
     ! 500     format(3f12.6)
     !         stop

     vtmp(1:3,1) = bb(1:3,1) + bb(1:3,2)
     vtmp(1:3,2) = bb(1:3,1) + bb(1:3,3)
     vtmp(1:3,3) = bb(1:3,2) + bb(1:3,3)
     vtmp(1:3,4) = bb(1:3,1) - bb(1:3,2)
     vtmp(1:3,5) = bb(1:3,1) - bb(1:3,3)
     vtmp(1:3,6) = bb(1:3,2) - bb(1:3,3)
     do i = 1,6
        avtmp(i) = dsqrt(sum(vtmp(:,i)**2))
     enddo
     if (avtmp(1) < avtmp(4)) then
        bb(1:3,4) = vtmp(1:3,1)
     else
        bb(1:3,4) = vtmp(1:3,4)
     endif
     if (avtmp(2) < avtmp(5)) then
        bb(1:3,5) = vtmp(1:3,2)
     else
        bb(1:3,5) = vtmp(1:3,5)
     endif
     if (avtmp(3) < avtmp(6)) then
        bb(1:3,6) = vtmp(1:3,3)
     else
        bb(1:3,6) = vtmp(1:3,6)
     endif
     !         do i = 1,6
     !            do j = i+1,6
     !               if (avtmp(i).gt.avtmp(j)) then
     !                  avtmp(7) = avtmp(j)
     !                  avtmp(j) = avtmp(i)
     !                  avtmp(i) = avtmp(7)
     !                  vtmp(:,7) = vtmp(:,j)
     !                  vtmp(:,j) = vtmp(:,i)
     !                  vtmp(:,i) = vtmp(:,7)
     !               endif
     !            enddo
     !         enddo
     !         do i = 1,5
     !            if (avtmp(i) .gt. avtmp(j)) stop 'getbb: sorting error'
     !         enddo
     !         bb(1:3,4) = vtmp(1:3,1)
     !         bb(1:3,5) = vtmp(1:3,2)
     !         bb(1:3,6)=  vtmp(1:3,3)

     do i = 1,6
        bb(1:3,i+6) = -bb(1:3,i)
     enddo

     ndim = nbb / 2
     allocate(cmat(ndim,ndim),rmat(ndim,ndim),work(ndim))

     rmat = 0d0
     do i=1,6
        rmat(1,i) = bb(1,i)*bb(1,i)   ! (j,k) = (x,x)
        rmat(2,i) = bb(1,i)*bb(2,i)   ! (j,k) = (x,y)
        rmat(3,i) = bb(1,i)*bb(3,i)   ! (j,k) = (x,z)
        rmat(4,i) = bb(2,i)*bb(2,i)   ! (j,k) = (y,y)
        rmat(5,i) = bb(2,i)*bb(3,i)   ! (j,k) = (y,z)
        rmat(6,i) = bb(3,i)*bb(3,i)   ! (j,k) = (z,z)
     enddo
     cmat(1:ndim,1:ndim) = dcmplx(rmat(1:ndim,1:ndim),0d0)

     call matcinv(ndim,cmat(1:ndim,1:ndim))
     rmat(:,:) = dreal(cmat(:,:))

     work = 0d0
     work(1) = 0.5d0
     work(4) = 0.5d0
     work(6) = 0.5d0
     do i = 1,6
        wbb(i) = 0d0
        do j =1,6
           wbb(i)   = wbb(i) + cmat(i,j)*work(j)
        enddo
        wbb(i+6) = wbb(i)
     enddo
     deallocate(cmat,rmat,work)

     goto 900
  endif ! lort

  write(*,*)'getbb: something is wrong.'
  stop

900 continue

  wbbsum = 0d0
  do ibb = 1,nbb
     wbbsum = wbbsum + wbb(ibb)
  enddo

  ! check if S[i] wbb(i)*bb(j,i)*bb(k,i) = delta(j,k)
  allocate (work(6))
  work = 0d0
  do i=1,nbb
     work(1) = work(1) + bb(1,i)*bb(1,i)*wbb(i) ! (j,k) = (x,x)
     work(2) = work(2) + bb(1,i)*bb(2,i)*wbb(i) ! (j,k) = (x,y)
     work(3) = work(3) + bb(1,i)*bb(3,i)*wbb(i) ! (j,k) = (x,z)
     work(4) = work(4) + bb(2,i)*bb(2,i)*wbb(i) ! (j,k) = (y,y)
     work(5) = work(5) + bb(2,i)*bb(3,i)*wbb(i) ! (j,k) = (y,z)
     work(6) = work(6) + bb(3,i)*bb(3,i)*wbb(i) ! (j,k) = (z,z)
  enddo

  tmp = dabs(work(1)-1d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'
  tmp = dabs(work(2)-0d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'
  tmp = dabs(work(3)-0d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'
  tmp = dabs(work(4)-1d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'
  tmp = dabs(work(5)-0d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'
  tmp = dabs(work(6)-1d0)
  if (tmp > eps) stop 'getbb: wbb is wrong'

  deallocate (work)

990 format(3f12.6)

  return
end subroutine getbb
!-----------------------------------------------------------------------
subroutine  getglat(n1,n2,n3,nshell,ng,qlat, &
     glat)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: qlat(3,3),glat(3,ng)

  nv = 0
  do i1 = -nshell,nshell
     do i2 = -nshell,nshell
        do i3 = -nshell,nshell
           nv = nv + 1
           glat(1:3,nv) = dble(i1)*qlat(1:3,1)/dble(n1) &
                + dble(i2)*qlat(1:3,2)/dble(n2) &
                + dble(i3)*qlat(1:3,3)/dble(n3)
        enddo
     enddo
  enddo

  if (nv /= ng) stop "getglat: wrong ng"

  return
end subroutine getglat
!-----------------------------------------------------------------------
subroutine  sortvec(ndat,vec)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: vec(3,ndat),vtmp(3,ndat),dist(ndat)
  integer(4) :: idat(ndat)

  vtmp = vec
  do i = 1,ndat
     dist(i) = dsqrt(sum(vtmp(:,i)**2))
     idat(i) = i
  enddo

  do j = 2,ndat
     d = dist(j)
     do i = j-1,1,-1
        if (dist(i)<=d) goto 999
        dist(i+1) = dist(i)
        idat(i+1) = idat(i)
     enddo
     i = 0
999  dist(i+1) = d
     idat(i+1) = j
  enddo

  do i = 1,ndat
     vec(1:3,i) = vtmp(1:3,idat(i))
  enddo

  do i = 1,ndat-1
     d1 = dsqrt(sum(vec(:,i)**2))
     d2 = dsqrt(sum(vec(:,i+1)**2))
     if (d1 > d2) stop 'sortvec: sorting error!'
  enddo

  return
end subroutine sortvec
!-----------------------------------------------------------------------
subroutine  sortvec2(ndat,vec,dist,idat)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: vec(3,ndat),vtmp(3,ndat),dist(ndat)
  integer :: idat(ndat)

  vtmp = vec
  do i = 1,ndat
     dist(i) = dsqrt(sum(vtmp(:,i)**2))
     idat(i) = i
  enddo

  do j = 2,ndat
     d = dist(j)
     do i = j-1,1,-1
        if (dist(i)<=d) goto 999
        dist(i+1) = dist(i)
        idat(i+1) = idat(i)
     enddo
     i = 0
999  dist(i+1) = d
     idat(i+1) = j
  enddo

  do i = 1,ndat
     vec(1:3,i) = vtmp(1:3,idat(i))
  enddo

  do i = 1,ndat-1
     d1 = dsqrt(sum(vec(:,i)**2))
     d2 = dsqrt(sum(vec(:,i+1)**2))
     if (d1 > d2) stop 'sortvec: sorting error!'
     if (abs(d1-dist(i)) > 1.d-4) &
          stop 'sortvec: sorting error in d!'
  enddo

  return
end subroutine sortvec2
!-----------------------------------------------------------------------
subroutine  get3vec(nv,vec,vv)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps1 = 0.001d0, eps2 = 0.25d0, eps3 = 0.7d0)
  real(8) :: vec(3,nv),vv(3,3),vtmp(3),vtmp2(3)

  ! vv1
  do i1 = 1,nv-2
     vtmp(1:3) = vec(1:3,i1)
     c1 = sum(vtmp(:)*vtmp(:))
     if (c1 > eps1) goto 777
  enddo
776 continue
  stop 'get3vec: cannot find vv1'
777 continue
  vv(1:3,1) = vtmp(1:3)

  ! vv2
  ini = i1+1
  do i2 = ini,nv-1
     vtmp(1:3) = vec(1:3,i2)
     c1 = sum(vtmp(:)*vv(:,1))
     c2 = dsqrt(sum(vtmp(:)*vtmp(:)))
     c3 = dsqrt(sum(vv(:,1)*vv(:,1)))
     c4 = 1d0 - dabs(c1/(c2*c3))
     if (c4 > eps2) goto 888
  enddo
887 continue
  stop 'get3vec: cannot find vv2'
888 continue
  vv(1:3,2) = vtmp(1:3)

  ! vv3
  vtmp2(1) = vv(2,1)*vv(3,2) - vv(3,1)*vv(2,2)
  vtmp2(2) = vv(3,1)*vv(1,2) - vv(1,1)*vv(3,2)
  vtmp2(3) = vv(1,1)*vv(2,2) - vv(2,1)*vv(1,2)

  ini = i2 + 1
  do i3 = ini,nv
     vtmp(1:3) = vec(1:3,i3)
     c1 = sum(vtmp(:)*vtmp2(:))
     c2 = dsqrt(sum(vtmp(:)*vtmp(:)))
     c3 = dsqrt(sum(vtmp2(:)*vtmp2(:)))
     c4 = dabs(c1/(c2*c3))
     if (c4 > eps3) goto 999
  enddo
998 continue
  stop 'get3vec: cannot find vv3'
999 continue
  vv(1:3,3) = vtmp(1:3)


  return
end subroutine get3vec
!-----------------------------------------------------------------------
subroutine  writebb2(ifbb,wbb,bb, ikbidx,ku,kbu, nqbz,nbb)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-4)
  integer :: ikbidx(nbb,nqbz) !iko_ixs(2),iko_fxs(2),noxs(2)
  real(8) :: wbb(nbb),bb(3,nbb),bb2(3),ku(3,nqbz),kbu(3,nbb,nqbz)
  open(newunit=ifbb,file='BBVEC')
  write(ifbb,*)'nbb,nqbz'
  write(ifbb,*)nbb,nqbz
  do i = 1,nbb
     write(ifbb,*)bb(1:3,i),wbb(i)
  enddo
  do iq = 1,nqbz
     write(ifbb,999)iq,ku(:,iq)
     do ib = 1,nbb
        itmp = ikbidx(ib,iq)
        if (itmp > nqbz .OR. itmp < 1) stop 'writebb: ikbidx error!'
        write(ifbb,998)iq,ib,ikbidx(ib,iq),kbu(:,ib,iq)
     enddo
  enddo
  close(ifbb)
998 format(3i5,3f24.16)
999 format(i5,3f24.16)
end subroutine writebb2
!------------------------------------------------------------
subroutine  writebb(ifbb,wbb,bb, &
     ikbidx,ku,kbu, &
     iko_ixs,iko_fxs,noxs, &
     nspin,nqbz,nbb)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-4)
  integer :: iopen, &
       ikbidx(nbb,nqbz), &
       iko_ixs(2),iko_fxs(2),noxs(2)
  real(8) :: wbb(nbb),bb(3,nbb),bb2(3), &
       ku(3,nqbz),kbu(3,nbb,nqbz)
  ifbb = iopen('BBVEC',1,-1,0)
  write(ifbb,*)'nbb,nqbz'
  write(ifbb,*)nbb,nqbz
  do i = 1,nbb
     write(ifbb,*)bb(1:3,i),wbb(i)
  enddo
  do iq = 1,nqbz
     write(ifbb,999)iq,ku(:,iq)
     do ib = 1,nbb
        itmp = ikbidx(ib,iq)
        if (itmp > nqbz .OR. itmp < 1) &
             stop 'writebb: ikbidx error!'
        write(ifbb,998)iq,ib,ikbidx(ib,iq),kbu(:,ib,iq)
     enddo
  enddo
  write(ifbb,*)'nspin'
  write(ifbb,*)nspin
  do is = 1,nspin
     write(ifbb,*)iko_ixs(is),iko_fxs(is),noxs(is)
  enddo

998 format(3i5,3f24.16)
999 format(i5,3f24.16)

  return
end subroutine writebb
!-----------------------------------------------------------------------
subroutine kbbindx(qbz,ginv,bb, &
     nqbz,nbb, &
     ikbidx,ku,kbu)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  real(8) :: qbz(3,nqbz),ginv(3,3),bb(3,nbb), &
       ku(3,nqbz),kbu(3,nbb,nqbz), &
       q(3),qu(3),g(3)
  integer :: ikbidx(nbb,nqbz)

  do ik = 1,nqbz
     q(:) = qbz(:,ik)
     call iqindx2(q,ginv,qbz,nqbz, iq, qu) ! q-qu= G vectors.
     ku(:,ik) = qu(:)
     do ib = 1,nbb
        q(:) = qbz(:,ik) + bb(:,ib)
        call iqindx2(q,ginv,qbz,nqbz, iq, qu) ! q-qu= G vectors.
        ikbidx(ib,ik) = iq
        kbu(:,ib,ik) = qu(:)
     enddo
  enddo

  return
end subroutine kbbindx
!-----------------------------------------------------------------------
subroutine ewindow(is,ieo_swt,iei_swt,itout_i,itout_f,itin_i,itin_f, &
     eomin,eomax,eimin,eimax,ef,qbz,ikbidx, &
     nbbelow,nbabove, &
     nqbz,nbb,nband,nwf,nspin, &
     iko_i,iko_f,iki_i,iki_f, &
     ikbo_i,ikbo_f,ikbi_i,ikbi_f, &
     iko_ix,iko_fx,nox, &
     leout,lein)

  use m_readeigen,only:readeval
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  real(8),allocatable :: eqt(:)
  real(8) :: qbz(3,nqbz),q(3)
  integer:: ikbidx(nbb,nqbz)
  integer:: iko_i(nqbz),iko_f(nqbz), &
       iki_i(nqbz),iki_f(nqbz), &
       ikbo_i(nbb,nqbz),ikbo_f(nbb,nqbz), &
       ikbi_i(nbb,nqbz),ikbi_f(nbb,nqbz)
  logical :: leout,lein
  real(8):: enwfmax,eeee
  ! energy window switches
  leout = .false.
  if (ieo_swt == 1) then
     leout = .true.
     write(*,*)'Outer energy window, on'
     write(*,*)'eomin, eomax =',eomin,eomax
  else
     write(*,*)'Outer band window'
     write(*,*)'itout_i, itout_f =',itout_i,itout_f
  endif

  lein  = .false.
  !      if (.not.leout) iei_swt = 0
  if (iei_swt == 1) then
     lein = .true.
     write(*,*)'Inner energy window, on'
     write(*,*)'eimin, eimax =',eimin,eimax
  elseif (iei_swt == 2) then
     lein = .true.
     write(*,*)'Inner band window, on'
     write(*,*)'itin_i, itin_f =',itin_i,itin_f
  else
     write(*,*)'Inner energy window, off'
  endif

  ! outer energy window ---------------------------------------
  ! no energy window
  if ( .NOT. leout) then
     iko_i(:) = itout_i
     iko_f(:) = itout_f
     ikbo_i(:,:) = itout_i
     ikbo_f(:,:) = itout_f

     iko_ix = itout_i
     iko_fx = itout_f
     if (iko_ix > iko_fx) stop 'ewindow: error!'
     nox = iko_fx - iko_ix + 1
  else ! leout
     allocate (eqt(nband))
     ! k
     do iq = 1,nqbz
        q = qbz(:,iq)
        eqt=readeval(q,is)
        call erange(1,eomin,eomax,ef,eqt,nband,nwf,iti,itf)
        iko_i(iq) = iti
        iko_f(iq) = itf
     enddo
     ! k + b
     do iq = 1,nqbz
        do ib = 1,nbb
           ikb = ikbidx(ib,iq)
           q = qbz(:,ikb)
           eqt= readeval(q,is)
           call erange(1,eomin,eomax,ef,eqt,nband,nwf,iti,itf)
           ikbo_i(ib,iq) = iti
           ikbo_f(ib,iq) = itf
        enddo
     enddo

     deallocate (eqt)
  endif ! leout

  ! search iko_ix and iko_fx
  iko_ix = nband
  iko_fx = 1
  do iq = 1,nqbz
     if (iko_i(iq) < iko_ix) iko_ix = iko_i(iq)
     if (iko_f(iq) > iko_fx) iko_fx = iko_f(iq)
  enddo
  ! m, 080110
  iko_ix = iko_ix - nbbelow
  if (iko_ix < 1) iko_ix = 1
  iko_fx = iko_fx + nbabove
  if (iko_fx > nband) iko_fx = nband
  if (iko_ix > iko_fx) stop 'ewindow: error!'
  nox = iko_fx - iko_ix + 1


  ! inner energy window ---------------------------------------
  ! no energy window
  if (iei_swt == 0) then
     iki_i(:) = 0
     iki_f(:) = -1
     ikbi_i(:,:) = 0
     ikbi_f(:,:) = -1
  elseif (iei_swt == 2) then
     iki_i(:) = itin_i
     iki_f(:) = itin_f
     ikbi_i(:,:) = itin_i
     ikbi_f(:,:) = itin_f
  elseif (iei_swt == 1) then
     allocate (eqt(nband))
     ! k
     do iq = 1,nqbz
        q = qbz(:,iq)
        eqt=readeval(q,is)
        call erange(2,eimin,eimax,ef,eqt,nband,nwf,iti,itf)
        iki_i(iq) = iti
        iki_f(iq) = itf
     enddo
     ! k + b
     do iq = 1,nqbz
        do ib = 1,nbb
           ikb = ikbidx(ib,iq)
           q = qbz(:,ikb)
           eqt= readeval(q,is)
           call erange(2,eimin,eimax,ef,eqt,nband,nwf,iti,itf)
           ikbi_i(ib,iq) = iti
           ikbi_f(ib,iq) = itf
        enddo
     enddo

     deallocate (eqt)
  else ! iei_swt
     stop 'ewindow: wrong switch for innner window!'
  endif ! iei_swt

  return
end subroutine ewindow
!-----------------------------------------------------------------------
subroutine erange(iflg,emin,emax,ef,eqt,nband,nwf, &
     iti,itf)
  ! iflg = 1: for outer energy window
  ! iflg = 2: for inner energy window
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: eqt(nband),eev(nband)

  eev(1:nband) = (eqt(1:nband) - ef) * rydberg()

  iti = nband
  itf = 1

  do it = nband,1,-1
     if (eev(it) >= emin) iti = it
  enddo

  do it = 1,nband
     if (eev(it) <= emax) itf = it
  enddo
  !!
  if (iflg == 1) then
     !         if (iti .gt. itf)
     !     &       stop 'erange: outer energy window too narrow'
     !         if (eev(iti) .gt. emax)
     !     &       stop 'erange: outer enrgy window too low'
     !         if (eev(itf) .lt. emin)
     !     &       stop 'erange: outer enrgy window too high'
     neo = itf - iti + 1
     if (neo < nwf) then
        print *,'neo,itf,iti=',neo,itf,iti,nwf
        stop 'energy range: outer energy window error'
     endif
  elseif (iflg == 2) then
     nei = itf - iti + 1
     if (nei > nwf) then
        stop 'energy range: inner energy window too wide'
     else

        evi = eev(iti)
        evf = eev(itf)
        if (iti > itf .OR. evi > emax .OR. evf < emin) then
           iti = 0
           itf = -1
        endif
     endif
  else
     stop 'erange: iflg error'
  endif

  return
end subroutine erange
!-----------------------------------------------------------------------
subroutine chk_ewindow(ifbb,ispin,nspin,nqbz,nbb,iko_ix,iko_fx)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)

  logical :: lbb,luu
  real(8) :: dummy(4)
  integer :: idummy(4),iti(nspin),itf(nspin),nt(nspin)


  inquire(file='BBVEC',exist=lbb)
  if (ispin == 1) then
     inquire(file='UUU',exist=luu)
  else
     inquire(file='UUD',exist=luu)
  endif

  if ( .NOT. lbb) stop 'chk_ewindow: Cannot find BBVEC'
  if ( .NOT. luu) stop 'chk_ewindow: Cannot find UUU/UUD'

  ifbb = iopen('BBVEC',1,0,0)
  read(ifbb,*)
  read(ifbb,*)nbb2,nqbz2
  if (nqbz /= nqbz2) stop 'chk_ewindow: nqbz is wrong!'
  if (nbb /= nbb2) stop 'chk_ewindow: nbb is wrong!'
  do i = 1,nbb
     read(ifbb,*)dummy(1:4)
  enddo
  do iq = 1,nqbz
     read(ifbb,*)idummy(1),dummy(1:3)
     do ib = 1,nbb
        read(ifbb,*)idummy(1:3),dummy(1:3)
     enddo
  enddo
  read(ifbb,*)
  read(ifbb,*)nspin2
  if (nspin2 /= nspin) stop 'chk_ewindow: nspin is wrong!'
  do is = 1,nspin
     read(ifbb,*)iti(is),itf(is),nt(is)
  enddo

  ifbb = iclose('BBVEC')

  iko_ix2 = iti(ispin)
  iko_fx2 = itf(ispin)

  if (iko_ix2 < iko_ix) then
     iko_ix = iko_ix2
  elseif (iko_ix2 > iko_ix) then
     stop 'chk_ewindow: iko_ix error'
  endif

  if (iko_fx2 > iko_fx) then
     iko_fx = iko_fx2
  elseif (iko_fx2 < iko_fx) then
     stop 'chk_ewindow: iko_fx error'
  endif

  return
end subroutine chk_ewindow
!-----------------------------------------------------------------------
subroutine readuu(is,iti,itf,ikbidx, &
     nqbz,nbb, &
     uum)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  integer :: ikbidx(nbb,nqbz)
  complex(8) :: uum(iti:itf,iti:itf,nbb,nqbz)
  complex(8),allocatable::uumread(:,:)
  if (is == 1) then
     ifuu      = iopen('UUU',0,0,0)
  else
     ifuu      = iopen('UUD',0,0,0)
  endif

  read(ifuu)
  read(ifuu)nqbz2,nbb2,iti2,itf2
  !      if (nqbz2 .ne. nqbz) stop 'readuu: nqbz error'
  !      if (nbb2 .ne. nbb) stop 'readuu: nbb error'
  !      if (iti2 .ne. iti) stop 'readuu: iti error'
  !      if (itf2 .ne. itf) stop 'readuu: itf error'
  if (nqbz2 /= nqbz)call rx('readuu: nqbz2/= nqbz')
  if (nbb2  /= nbb) call rx('readuu: nbb2 /= nbb')
  if (iti2 > iti)   call rx('readuu: iti error')
  if (itf2 < itf)   call rx('readuu: itf error')
  allocate(uumread(iti2:itf2,iti2:itf2))

  do iqbz = 1,nqbz
     do ibb = 1,nbb
        read(ifuu)iflg
        if (iflg == -10) then
           read(ifuu) iqbz2,ibb2,ikbidx2
           read(ifuu) ((uumread(j1,j2),j1=iti2,itf2),j2=iti2,itf2)
           !     &         ((uum(j1,j2,ibb,iqbz),j1=iti,itf),j2=iti,itf)
           uum(iti:itf,iti:itf,ibb,iqbz)=uumread(iti:itf,iti:itf)
           if (iqbz2 /= iqbz) stop 'readuu: iqbz error'
           if (ibb2 /= ibb) stop 'readuu: ibb error'
           if (ikbidx2 /= ikbidx(ibb,iqbz)) &
                stop 'readuu: ikbidx error'
        elseif(iflg == -20) then
           read(ifuu)iqbz2,ibb2,iqtmp,ibbtmp
           if (iqbz2 /= iqbz) stop 'readuu: iqbz error'
           if (ibb2 /= ibb) stop 'readuu: ibb error'
           if (iqtmp >= iqbz) stop 'readuu: iqtmp error'
           if (ibbtmp <= 0) stop 'readuu: ibbtmp error'
           if (ibbtmp > nbb) stop 'readuu: ibbtmp error'

           do j1 = iti,itf
              do j2 = iti,itf
                 uum(j1,j2,ibb,iqbz) = &
                      dconjg(uum(j2,j1,ibbtmp,iqtmp))
              enddo
           enddo
        else
           stop 'readuu: iflg error'
        endif
     enddo
  enddo

  if (is == 1) then
     ifu = iclose('UUU')
  else
     ifu = iclose('UUD')
  endif

  return
end subroutine readuu
!-----------------------------------------------------------------------
subroutine chkuu(is,iti,itf,ikbidx,uum, &
     nqbz,nbb)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-6)
  integer(4) :: ikbidx(nbb,nqbz)
  complex(8) :: uum(iti:itf,iti:itf,nbb,nqbz)

  do iqbz = 1,nqbz
     do ibb = 1,nbb

        iqb = ikbidx(ibb,iqbz)

        if (iqb < iqbz) then
           iqtmp = iqb
           do ibbtmp = 1,nbb
              iqbtmp = ikbidx(ibbtmp,iqtmp)
              if (iqbtmp == iqbz) goto 900
           enddo
           stop 'chkuu: iqbtmp error'
900        continue
           do j1 = iti,itf
              do j2 = iti,itf
                 ! real part
                 rtmp = dabs(dreal(uum(j1,j2,ibb,iqbz)) &
                      - dreal(uum(j2,j1,ibbtmp,iqtmp)))
                 ! imag part
                 ctmp = dabs(dimag(uum(j1,j2,ibb,iqbz)) &
                      + dimag(uum(j2,j1,ibbtmp,iqtmp)))
                 if (rtmp > eps) stop 'chkuu: real part error'
                 if (ctmp > eps) stop 'chkuu: imag part error'
              enddo
           enddo

        endif

     enddo
  enddo

  return
end subroutine chkuu
!-----------------------------------------------------------------------
subroutine readuu0(is,iti,itf, &
     nqbz, &
     uum0)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8) :: uum0(iti:itf,iti:itf,1,nqbz)

  if (is == 1) then
     ifuu      = iopen('UU0U',0,0,0)
  else
     ifuu      = iopen('UU0D',0,0,0)
  endif

  read(ifuu)
  read(ifuu)nqbz2,nbb2,iti2,itf2
  if (nqbz2 /= nqbz) stop 'readuu: nqbz error'
  if (iti2 /= iti) stop 'readuu: iti error'
  if (itf2 /= itf) stop 'readuu: itf error'

  do iqbz = 1,nqbz
     do ibb = 1,1

        read(ifuu)iflg
        if (iflg /= -10) stop 'readuu0: iflg error'
        read(ifuu) iqbz2
        read(ifuu) &
             ((uum0(j1,j2,ibb,iqbz),j1=iti,itf),j2=iti,itf)
        if (iqbz2 /= iqbz) stop 'readuu0: iqbz error'

     enddo
  enddo

  if (is == 1) then
     ifu = iclose('UU0U')
  else
     ifu = iclose('UU0D')
  endif

  return
end subroutine readuu0
!-----------------------------------------------------------------------
subroutine readpsig(is,iti,itf, &
     nqbz,nwf, &
     psig)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  complex(8) :: psig(iti:itf,nwf,nqbz)
  complex(8),allocatable:: psigread(:,:,:)
  print *,'readpsig:'
  if (is == 1) then
     ifpsig      = iopen('PSIGU',0,0,0)
  else
     ifpsig      = iopen('PSIGD',0,0,0)
  endif
  read(ifpsig)nqbz2,iti2,itf2,nwf2
  if (nqbz2 /= nqbz)call rx('readpsig: nqbz error')
  if (nwf2  /= nwf) call rx('readpsig: nwf  error')
  if (iti2 > iti)   call rx('readpsig: iti range error')
  if (itf2 < itf)   call rx('readpsig: itf range error')
  allocate(psigread(iti2:itf2,nwf2,nqbz2))
  do iqbz = 1,nqbz
     ! ccccccccccccc
     !         print *,'readpsig:',iti2,itf2,nwf,iqbz
     ! ccccccccccccc
     read(ifpsig)iqbz2
     if(iqbz2 /= iqbz) call rx('readpsig: iqbz error')
     read(ifpsig)((psigread(j1,j2,iqbz),j1=iti2,itf2),j2=1,nwf)
  enddo
  psig(iti:itf,:,:)=psigread(iti:itf,:,:)
  if (is == 1) then
     ifu = iclose('PSIGU')
  else
     ifu = iclose('PSIGD')
  endif

  return
end subroutine readpsig
!-----------------------------------------------------------------------
subroutine wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer :: n1,n2,n3,nrws
  real(8) :: alat,plat(3,3)
  integer :: irws(n1*n2*n3*8)
  real(8) :: rws(3,n1*n2*n3*8),drws(n1*n2*n3*8)

  integer(4):: ii0(3,8),isort(8), &
       iwork1(n1*n2*n3*8),iwork2(n1*n2*n3*8)
  real(8):: rr(3,8),dd(8)
  parameter (tol=1.d-6)


  nrws = 0
  do i1=0,n1-1
     do i2=0,n2-1
        do i3=0,n3-1

           n = 0
           do j1=0,1
              do j2=0,1
                 do j3=0,1
                    n = n+1
                    ii0(1,n) = i1 - j1*n1
                    ii0(2,n) = i2 - j2*n2
                    ii0(3,n) = i3 - j3*n3
                 enddo ! j3
              enddo ! j2
           enddo ! j1

           do n=1,8
              rr(1:3,n) =  ( plat(1:3,1)*dble(ii0(1,n)) &
                   +   plat(1:3,2)*dble(ii0(2,n)) &
                   +   plat(1:3,3)*dble(ii0(3,n)) )
           enddo
           call sortvec2(8,rr,dd,isort)
           ndegen = 1
           do n=2,8
              if ((dd(n)-dd(1)) <= tol) ndegen = n
           enddo

           do n=1,ndegen
              nrws = nrws + 1
              rws(1:3,nrws) = rr(1:3,n)
              drws(nrws) = dd(n)
              irws(nrws) = ndegen
           enddo

        enddo ! i3
     enddo ! i2
  enddo ! i1

  call sortvec2(nrws,rws,drws,iwork1)
  iwork2(1:nrws) = irws(1:nrws)
  do n=1,nrws
     irws(n) = iwork2(iwork1(n))
  enddo

  return
end subroutine wigner_seitz
!-----------------------------------------------------------------------
subroutine super_cell(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4):: n1,n2,n3,nrws
  real(8):: alat,plat(3,3)
  integer(4) :: irws(n1*n2*n3)
  real(8) :: rws(3,n1*n2*n3),drws(n1*n2*n3)
  real(8):: rr(3),dd
  parameter (tol=1.d-6)

  nrws = 0
  do i1=0,n1-1
     do i2=0,n2-1
        do i3=0,n3-1
           rr(1:3) =  ( plat(1:3,1)*dble(i1) &
                + plat(1:3,2)*dble(i2) &
                + plat(1:3,3)*dble(i3) )
           nrws = nrws + 1
           rws(1:3,nrws) = rr(1:3)
           drws(nrws) = dsqrt(sum(rr(:)**2))
           irws(nrws) = 1
        enddo ! i3
     enddo ! i2
  enddo ! i1

  if (nrws /= n1*n2*n3) stop 'super_cell: nrws error!'

  return
end subroutine super_cell
!-----------------------------------------------------------------------
