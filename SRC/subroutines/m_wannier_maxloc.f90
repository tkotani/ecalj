!> utilities for maxloc
module m_maxloc0
contains
  subroutine getrt(qbz,qbas,plat,n1,n2,n3,nqbz, rt,rt8,qbz0)
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
999   dist(i+1) = d
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
999   dist(i+1) = d
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
  subroutine  writebb(ifbb,wbb,bb, ikbidx,ku,kbu, iko_ixs,iko_fxs,noxs, nspin,nqbz,nbb)
    use m_ftox
    implicit integer (i-n)
    implicit real*8(a-h,o-z)
    parameter (eps = 1d-4)
    integer :: ikbidx(nbb,nqbz), iko_ixs(2),iko_fxs(2),noxs(2)
    real(8) :: wbb(nbb),bb(3,nbb),bb2(3),ku(3,nqbz),kbu(3,nbb,nqbz)
    open(newunit=ifbb,file='BBVEC')
    write(ifbb,ftox)'nbb,nqbz'
    write(ifbb,ftox)nbb,nqbz
    do i = 1,nbb
      write(ifbb,ftox)bb(1:3,i),wbb(i)
    enddo
    do iq = 1,nqbz
      write(ifbb,ftox)iq,ku(:,iq)
      do ib = 1,nbb
        itmp = ikbidx(ib,iq)
        if (itmp > nqbz .OR. itmp < 1) call rx('writebb: ikbidx error!')
        write(ifbb,ftox)iq,ib,ikbidx(ib,iq),kbu(:,ib,iq)
      enddo
    enddo
    write(ifbb,ftox)'nspx'
    write(ifbb,ftox)nspin
    do is = 1,nspin
      write(ifbb,ftox)iko_ixs(is),iko_fxs(is),noxs(is)
    enddo
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
       nqbz,nbb,nband,nwf, &
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
  ! subroutine chk_ewindow(ifbb,ispin,nspin,nqbz,nbb,iko_ix,iko_fx)
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   logical :: lbb,luu
  !   real(8) :: dummy(4)
  !   integer :: idummy(4),iti(nspin),itf(nspin),nt(nspin)


  !   inquire(file='BBVEC',exist=lbb)
  !   if (ispin == 1) then
  !      inquire(file='UUU',exist=luu)
  !   else
  !      inquire(file='UUD',exist=luu)
  !   endif

  !   if ( .NOT. lbb) stop 'chk_ewindow: Cannot find BBVEC'
  !   if ( .NOT. luu) stop 'chk_ewindow: Cannot find UUU/UUD'

  ! !  ifbb = iopen('BBVEC',1,0,0)
  !   open(newunit=ifbb,file='BBVEC')
  !   read(ifbb,*)
  !   read(ifbb,*)nbb2,nqbz2
  !   if (nqbz /= nqbz2) stop 'chk_ewindow: nqbz is wrong!'
  !   if (nbb /= nbb2) stop 'chk_ewindow: nbb is wrong!'
  !   do i = 1,nbb
  !      read(ifbb,*)dummy(1:4)
  !   enddo
  !   do iq = 1,nqbz
  !      read(ifbb,*)idummy(1),dummy(1:3)
  !      do ib = 1,nbb
  !         read(ifbb,*)idummy(1:3),dummy(1:3)
  !      enddo
  !   enddo
  !   read(ifbb,*)
  !   read(ifbb,*)nspin2
  !   if (nspin2 /= nspin) stop 'chk_ewindow: nspin is wrong!'
  !   do is = 1,nspin
  !      read(ifbb,*)iti(is),itf(is),nt(is)
  !   enddo

  !   !ifbb = iclose('BBVEC')
  !   close(ifbb)
  !   iko_ix2 = iti(ispin)
  !   iko_fx2 = itf(ispin)

  !   if (iko_ix2 < iko_ix) then
  !      iko_ix = iko_ix2
  !   elseif (iko_ix2 > iko_ix) then
  !      stop 'chk_ewindow: iko_ix error'
  !   endif

  !   if (iko_fx2 > iko_fx) then
  !      iko_fx = iko_fx2
  !   elseif (iko_fx2 < iko_fx) then
  !      stop 'chk_ewindow: iko_fx error'
  !   endif

  !   return
  ! end subroutine chk_ewindow
  !-----------------------------------------------------------------------
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
900       continue
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
  ! subroutine readuu0(is,iti,itf, &
  !      nqbz, &
  !      uum0)
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)
  !   complex(8) :: uum0(iti:itf,iti:itf,1,nqbz)

  !   if (is == 1) then
  ! !     ifuu      = iopen('UU0U',0,0,0)
  !     open(newunit=ifuu,file='UU0U',form='unformatted')
  !   else
  ! !     ifuu      = iopen('UU0D',0,0,0)
  !     open(newunit=ifuu,file='UU0D',form='unformatted')
  !   endif

  !   read(ifuu)
  !   read(ifuu)nqbz2,nbb2,iti2,itf2
  !   if (nqbz2 /= nqbz) stop 'readuu: nqbz error'
  !   if (iti2 /= iti) stop 'readuu: iti error'
  !   if (itf2 /= itf) stop 'readuu: itf error'

  !   do iqbz = 1,nqbz
  !      do ibb = 1,1

  !         read(ifuu)iflg
  !         if (iflg /= -10) stop 'readuu0: iflg error'
  !         read(ifuu) iqbz2
  !         read(ifuu) &
  !              ((uum0(j1,j2,ibb,iqbz),j1=iti,itf),j2=iti,itf)
  !         if (iqbz2 /= iqbz) stop 'readuu0: iqbz error'

  !      enddo
  !   enddo
  !   close(ifu)
  !   ! if (is == 1) then
  !   !    close(ifu)! = iclose('UU0U')
  !   ! else
  !   !    close(ifu)! = iclose('UU0D')
  !   ! endif
  !   return
  ! end subroutine readuu0
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
      !    ifpsig      = iopen('PSIGU',0,0,0)
      open(newunit=ifpsig,file='PSIGU',form='unformatted')
    else
      !     ifpsig      = iopen('PSIGD',0,0,0)
      open(newunit=ifpsig,file='PSIGD',form='unformatted')
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
    close(ifu)
    !if (is == 1) then
    !   ! = iclose('PSIGU')
    !else
    !   ifu = iclose('PSIGD')
    !endif

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
  end subroutine wigner_seitz
  !-----------------------------------------------------------------------
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
  ! subroutine pick_nwf(ovlp,iti,itf,nwf,isort)

  !   implicit real*8(a-h,o-z)
  !   implicit integer (i-n)
  !   real(8) :: ovlp(iti:itf),otmp(iti:itf)
  !   integer(4) :: isort(nwf),istate(iti:itf)

  !   ! initial
  !   do it = iti,itf
  !      istate(it) = it
  !      otmp(it) = ovlp(it)
  !   enddo

  !   ! sorting
  !   do it1 = iti,itf-1
  !      do it2 = it1+1,itf
  !         if (ovlp(it1) < ovlp(it2)) then
  !            tmp = ovlp(it2)
  !            ovlp(it2) = ovlp(it1)
  !            ovlp(it1) = tmp
  !            itmp = istate(it2)
  !            istate(it2) = istate(it1)
  !            istate(it1) = itmp
  !         endif
  !      enddo
  !   enddo

  !   ! sort check
  !   do i1 = iti,itf-1
  !      it1 = istate(i1)
  !      it2 = istate(i1+1)
  !      if (otmp(it1) < otmp(it2)) stop 'pick_nwf: sort error'
  !   enddo

  !   ! pick largest nwf states
  !   do it = 1,nwf
  !      itmp = iti - 1 + it
  !      isort(it) = istate(itmp)
  !   enddo

  !   return
  ! end subroutine pick_nwf
  !-----------------------------------------------------------------------
  ! subroutine read_cnq0(ifhoev,is,qwf0,qbz,ginv,ef, &
  !      itq, &
  !      nwf,nband,nqbz, &
  !      cnq0)

  !   implicit real*8(a-h,o-z)
  !   implicit integer (i-n)
  !   complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:)
  !   complex(8) :: cnq0(nband,nwf)
  !   real(8),allocatable :: eval(:)
  !   real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
  !   real(8) :: rydberg
  !   integer(4) :: itq(nwf)

  !   iq0 = iqindx(qwf0,ginv,qbz,nqbz)
  !   cnq0 = (0d0,0d0)

  !   ! open
  !   if (is == 1) then
  ! !     ifhoev = iopen('HOEV.UP',0,0,0)
  !      open(newunit=ifhoev,file='HOEV.UP',form='unformatted')
  !   elseif (is == 2) then
  ! !     ifhoev = iopen('HOEV.DN',0,0,0)
  !      open(newunit=ifhoev,file='HOEV.DN',form='unformatted')
  !   else
  !      stop 'read_cnq0: open error'
  !   endif

  !   ! read
  !   read(ifhoev)ndimh,nqtot
  !   if (ndimh /= nband) stop 'read_cnq0: nband error'
  !   if (nqtot < nqbz) stop 'read_cnq0: nqbz error'

  !   allocate(hks(nband,nband),oks(nband,nband), &
  !        cks(nband,nband),eval(nband))

  !   do iq = 1,nqbz
  !      read(ifhoev)iq2,q(1:3)
  !      read(ifhoev)hks(1:nband,1:nband)
  !      read(ifhoev)oks(1:nband,1:nband)
  !      read(ifhoev)cks(1:nband,1:nband)
  !      read(ifhoev)eval(1:nband)

  !      iq3 = iqindx(q,ginv,qbz,nqbz)
  !      if (iq3 /= iq) stop 'read_cnq0: iqindx error'
  !      if (iq3 == iq0) then
  !         do it = 1,nwf
  !            cnq0(:,it) = cks(:,itq(it))
  !            !               ev = (eval(itq(it)) - ef) * rydberg()
  !            !               write(*,*)'iwf,nwf,ev',ev
  !         enddo
  !         goto 99
  !      endif
  !   enddo
  !   stop 'read_cnq0: cannot find q0'
  ! 99 continue

  !   deallocate(hks,oks,cks,eval)

  !   close(ifi)
  !   ! if (is == 1) then
  !   !    ifi = iclose('HOEV.UP')
  !   ! else
  !   !    ifi = iclose('HOEV.DN')
  !   ! endif

  !   return
  ! end subroutine read_cnq0
  !-----------------------------------------------------------------------
  ! subroutine get_amnk(ifhoev,is,qwf0,qbz,ginv, &
  !      cnq0, &
  !      iko_ix,iko_fx,iko_i,iko_f, &
  !      nwf,nband,nqbz, &
  !      amnk)
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:), &
  !        wmat(:,:)
  !   complex(8) :: cnq0(nband,nwf), &
  !        amnk(iko_ix:iko_fx,nwf,nqbz)
  !   real(8),allocatable :: eval(:)
  !   real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
  !   integer(4) :: iko_i(nqbz),iko_f(nqbz)


  !   ! open
  !   if (is == 1) then
  ! !     ifhoev = iopen('HOEV.UP',0,0,0)
  !      open(newunit=ifhoev,file='HOEV.UP',form='unformatted')
  !   elseif (is == 2) then
  ! !     ifhoev = iopen('HOEV.DN',0,0,0)
  !      open(newunit=ifhoev,file='HOEV.DN',form='unformatted')
  !   else
  !      stop 'get_amnk: open error'
  !   endif

  !   ! read
  !   read(ifhoev)ndimh,nqtot
  !   if (ndimh /= nband) stop 'get_amnk: nlmto error'
  !   if (nqtot < nqbz) stop 'get_amnk: nqbz error'

  !   allocate(hks(nband,nband),oks(nband,nband), &
  !        cks(nband,nband),eval(nband), &
  !        wmat(nband,nwf))

  !   ! initialize
  !   amnk = (0d0,0d0)

  !   do iq = 1,nqbz
  !      read(ifhoev)iq2,q(1:3)
  !      read(ifhoev)hks(1:nband,1:nband)
  !      read(ifhoev)oks(1:nband,1:nband)
  !      read(ifhoev)cks(1:nband,1:nband)
  !      read(ifhoev)eval(1:nband)

  !      iq3 = iqindx(q, ginv,qbz,nqbz)
  !      if (iq3 /= iq) stop 'get_amnk: iqindx error'


  !      ! wmat = cnq0 * oks
  !      wmat = (0d0,0d0)
  !      do in = 1,nwf
  !         do ij = 1,nband
  !            do ii = 1,nband
  !               wmat(ij,in) = wmat(ij,in) &
  !                    + cnq0(ii,in) * oks(ij,ii)
  !            enddo
  !         enddo
  !      enddo

  !      ! amnk = cks^{*} * wmat
  !      do im = iko_i(iq),iko_f(iq)
  !         do in = 1,nwf
  !            do ij = 1,nband
  !               amnk(im,in,iq) = amnk(im,in,iq) + &
  !                    dconjg(cks(ij,im)) * wmat(ij,in)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   deallocate(hks,oks,cks,eval,wmat)
  !   close(ifi)
  ! !  if (is == 1) then
  ! !     ifi = iclose('HOEV.UP')
  ! !  else
  ! !     ifi = iclose('HOEV.DN')
  ! !  endif
  ! 999 format(i5,3f16.8)
  !   return
  ! end subroutine get_amnk
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
  subroutine init_iew(iko_ix,iko_fx,iko_i,iko_f, iki_i,iki_f, nwf,nband,nqbz, cnk)
    !Modify cnk to include all the inner space explictly. See Souza Eq.27 around.
    !2025-2-25, tk tests for NiO_magnon, where we faied for diagonalization when we set innerwindow -5eV to 1eV. -10eV works.
    ! It seems iki_f(iq)=2 is not accepted. I don't resolve this yet.
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
      q(:) = qbz(:,iq)
      eks= readeval (q,is)! read eigenvalues
      ham = 0d0
      do ii = 1,nwf
        do ij = 1,nwf
          do ik = iko_i(iq),iko_f(iq)
            ham(ii,ij) = ham(ii,ij) + dconjg(cnk(ik,ii,iq))*cnk(ik,ij,iq)*eks(ik)
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
            cnk2(in,il) = cnk2(in,il) + evecc(im,il) * cnk(in,im,iq)
          enddo
        enddo
      enddo
      cnk(:,:,iq) = cnk2
      eunk(:,iq) = eval(:)
    enddo
    deallocate (ham,eks,evecc,eval,cnk2)
    return
  end subroutine diag_unk
  !-----------------------------------------------------------------------
  ! subroutine chk_eunk(is,qbz,eunk,ef, &
  !      nqbz,nband,nwf)

  !   use m_readeigen,only:readeval
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   real(8) :: qbz(3,nqbz),eunk(nwf,nqbz)
  !   real(8),allocatable :: eks(:)

  !   allocate(eks(nband))

  !   do iq = 1,nqbz
  !      eks= readeval (qbz(:,iq),is)

  !      write(*,*)'Diag energy',iq
  !      do iband = 1,nwf
  !         eev = (eunk(iband,iq)-ef)*rydberg()
  !         write(*,*)iband,eunk(iband,iq),eev
  !      enddo

  !      write(*,*)'KS energy  ',iq
  !      do iband = 1,nband
  !         eev = (eks(iband)-ef)*rydberg()
  !         write(*,*)iband,eks(iband),eev
  !      enddo
  !   enddo

  !   deallocate(eks)

  !   return
  ! end subroutine chk_eunk
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
  ! subroutine get_rnm(mmn,bb,wbb,wbz, &
  !      nwf,nqbz,nbb, &
  !      rnm)! complex. See Eq.(4) and Eq.(22) in [1].
  !   implicit none
  !   integer:: in,im,ibb,iq,nbb,nqbz,nwf
  !   complex(8) :: mmn(nwf,nwf,nbb,nqbz)
  !   real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz)
  !   complex(8):: img=(0d0,1d0),rnm(3,nwf,nwf),rtmp
  !   rnm = 0d0
  !   do in = 1,nwf
  !      do im = 1,nwf
  !         if(im==in) cycle
  !         do ibb = 1,nbb
  !            rtmp = 0d0
  !            do iq = 1,nqbz
  !               rtmp = rtmp + img*mmn(in,im,ibb,iq)*wbz(iq)
  !            enddo
  !            rnm(:,in,im) = rnm(:,in,im) + wbb(ibb)*bb(:,ibb)*rtmp
  !         enddo
  !      enddo
  !   enddo
  ! end subroutine get_rnm
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
      !     ifbnd = iopen('spread.up',1,-1,0)
      open(newunit=ifbnd,file='spread.up')
    else
      !     ifbnd = iopen('spread.dn',1,-1,0)
      open(newunit=ifbnd,file='spread.dn')
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
    close(ifbnd)
    !  if (is == 1) then
    !     isx = iclose('spread.up')
    !  else
    !     isx = iclose('spread.dn')
    !  endif

    return
  end subroutine writeOmg
  !-----------------------------------------------------------------------
  ! subroutine writermn(is,mmn,bb,wbb,qbz,qbz0,wbz,rt, &
  !      nwf,nqbz,nbb,n1,n2,n3)
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr &
  !        ,ctmp1,ctmp2 &
  !        ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz) &
  !        ,c1,c2,c3
  !   real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),qbz0(3,nqbz),wbz(nqbz), &
  !        rt(3,nqbz),dmn,rtmp(3)

  !   pi = 4d0*atan(1d0)
  !   ci = (0d0,1d0)
  !   czero = (0d0,0d0)

  !   ! Berry connecrion (in the Wannier gauge)
  !   amn = czero
  !   do im = 1,nwf
  !      do in = 1,nwf
  !         dmn = 0d0
  !         if (im == in) dmn = 1d0
  !         do iq = 1,nqbz
  !            do ibb = 1,nbb
  !               !            ctmp1 = (mmn(im,in,ibb,iq) - dmn)*wbb(ibb)
  !               r1 = dimag(0.5d0*(mmn(im,in,ibb,iq)+mmn(in,im,ibb,iq))-dmn)
  !               r2 = dreal(0.5d0*(mmn(im,in,ibb,iq)-mmn(in,im,ibb,iq)))
  !               ctmp1 = dcmplx(r2,r1)
  !               ctmp1 = log(ctmp1+1d0)*wbb(ibb)
  !               do ix = 1,3 ! x,y,z
  !                  amn(ix,im,in,iq) = amn(ix,im,in,iq) &
  !                       + ci*imag(ctmp1) * bb(ix,ibb)
  !               enddo ! ix
  !            enddo ! ibb
  !         enddo ! iq
  !      enddo ! in
  !   enddo ! im
  !   amn = amn * ci

  !   ! <0m | r | Rn>
  !   rmn = czero
  !   do ir = 1,nqbz
  !      do iq = 1,nqbz
  !         !         rk = sum(rt(:,ir)*qbz(:,iq))
  !         rk = sum(rt(:,ir)*qbz0(:,iq))
  !         cikr = -ci * 2d0 * pi * rk
  !         ceikr = exp(cikr) * wbz(iq)
  !         do in = 1,nwf
  !            do im = 1,nwf
  !               do ix = 1,3
  !                  rmn(ix,im,in,ir) = rmn(ix,im,in,ir) + &
  !                       ceikr * amn(ix,im,in,iq)
  !               enddo ! ix
  !            enddo ! im
  !         enddo ! in
  !      enddo ! iq
  !   enddo ! ir

  !   ! output
  !   if (is == 1) then
  ! !     ifrmn = iopen('rmn.up',1,-1,0)
  !      open(newunit=ifrmn,file='rmn.up')
  !   else
  ! !     ifrmn = iopen('rmn.dn',1,-1,0)
  !      open(newunit=ifrmn,file='rmn.dn')
  !   endif

  !   write(ifrmn,*)'*** nwf,nsite'
  !   write(ifrmn,*)nwf,nqbz
  !   write(ifrmn,*)'*** rsite'
  !   write(ifrmn,*)rt
  !   write(ifrmn,*)'*** rmn'
  !   write(ifrmn,*)rmn
  !   write(ifrmn,*)'*** amn'
  !   write(ifrmn,*)amn
  !   close(ifrmn)
  ! !  if (is == 1) then
  ! !     isx = iclose('rmn.up')
  ! !  else
  ! !     isx = iclose('rmn.dn')
  ! !  endif

  !   return
  ! end subroutine writermn
  ! !-----------------------------------------------------------------------
  ! subroutine writemmn(is,mmn,bb,wbb,qbz,wbz,rt, &
  !      nwf,nqbz,nbb,n1,n2,n3)
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr &
  !        ,ctmp1,ctmp2 &
  !        ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz) &
  !        ,c1,c2,c3
  !   real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),wbz(nqbz), &
  !        rt(3,nqbz),dmn

  !   pi = 4d0*atan(1d0)
  !   ci = (0d0,1d0)
  !   czero = (0d0,0d0)

  !   ! output
  !   if (is == 1) then
  ! !     ifrmn = iopen('mmn.up',1,-1,0)
  !      open(newunit=ifrmn,file='mmn.up')
  !   else
  ! !     ifrmn = iopen('mmn.dn',1,-1,0)
  !      open(newunit=ifrmn,file='mmn.up')
  !   endif

  !   write(ifrmn,*)'*** nwf,nsite,nb'
  !   write(ifrmn,*)nwf,nqbz,nbb
  !   write(ifrmn,*)'*** mmn'
  !   write(ifrmn,*)mmn
  !   write(ifrmn,*)'*** bb'
  !   write(ifrmn,*)bb
  !   write(ifrmn,*)'*** wbb'
  !   write(ifrmn,*)wbb
  !   write(ifrmn,*)'*** qbz'
  !   write(ifrmn,*)qbz
  !   write(ifrmn,*)'*** wbz'
  !   write(ifrmn,*)wbz
  !   close(ifrmn)
  ! !  if (is == 1) then
  ! !     isx = iclose('mmn.up')
  ! !  else
  ! !     isx = iclose('mmn.dn')
  ! !  endif

  !   return
  ! end subroutine writemmn
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
  ! subroutine chk_dnk(is,eunk,qbz, &
  !      umnk,cnk, &
  !      iko_ix,iko_fx,iko_i,iko_f, &
  !      nband,nwf,nqbz)

  !   use m_readeigen,only:readeval
  !   implicit integer (i-n)
  !   implicit real*8(a-h,o-z)

  !   parameter (eps=1d-3)
  !   complex(8),allocatable :: dnk(:,:)
  !   complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz), &
  !        umnk(nwf,nwf,nqbz), &
  !        ctmp,ham(nwf,nwf),evecc(nwf,nwf)
  !   real(8) :: eunk(nwf,nqbz),qbz(3,nqbz),eks(nband),eval(nwf)
  !   integer(4) :: iko_i(nqbz),iko_f(nqbz)

  !   allocate(dnk(iko_ix:iko_fx,nwf))

  !   do iq = 1,nqbz
  !      eks= readeval (qbz(:,iq),is)
  !      dnk = (0d0,0d0)
  !      do imp = iko_i(iq),iko_f(iq)
  !         do in = 1,nwf
  !            do im = 1,nwf
  !               dnk(imp,in) = dnk(imp,in) &
  !                    + umnk(im,in,iq) * cnk(imp,im,iq)
  !            enddo
  !         enddo
  !      enddo


  !      ham = 0d0
  !      do ii=1,nwf
  !         do ij=1,nwf
  !            ctmp = 0d0
  !            do ik = iko_i(iq),iko_f(iq)
  !               ctmp = ctmp + eks(ik) * &
  !                    dconjg(dnk(ik,ii))*dnk(ik,ij)
  !            enddo
  !            ham(ii,ij) = ctmp
  !            write(97,"(3i5,2f12.6)")iq,ii,ij,dreal(ctmp),dimag(ctmp)
  !         enddo
  !      enddo
  !      call diag_hm(ham,nwf,eval,evecc)
  !      do ii = 1,nwf
  !         write(98,"(2i5,2f12.6)")iq,ii,eval(ii),eunk(ii,iq)
  !      enddo

  !   enddo

  !   deallocate(dnk)

  !   return
  ! end subroutine chk_dnk
  ! !-----------------------------------------------------------------------
  subroutine rot_hmnk(umnk,eunk, &
       nwf,nqbz, &
       hrotk)

    ! see Ref.[2] eq.24
    implicit integer (i-n)
    implicit real*8(a-h,o-z)

    complex(8) :: umnk(nwf,nwf,nqbz),hrotk(nwf,nwf,nqbz),ctmp
    real(8) :: eunk(nwf,nqbz)

    hrotk = (0d0,0d0)

    do iq = 1,nqbz
      do im = 1,nwf
        do in = 1,nwf
          do ii = 1,nwf
            hrotk(im,in,iq) = hrotk(im,in,iq) + &
                 dconjg(umnk(ii,im,iq))*eunk(ii,iq)*umnk(ii,in,iq)
          enddo
        enddo
      enddo
    enddo

    return
  end subroutine rot_hmnk
  !-----------------------------------------------------------------------
  subroutine get_hrotr_ws(hrotk,qbz,wbz, &
       rws,irws,drws, &
       nwf,nqbz,nrws, &
       hrotr)

    ! see Ref.[2] eq.25
    implicit integer (i-n)
    implicit real*8(a-h,o-z)

    complex(8) :: hrotk(nwf,nwf,nqbz),hrotr(nwf,nwf,nrws), &
         ci,cikr,ceikr,ctmp
    real(8) :: qbz(3,nqbz),wbz(nqbz),q(3),r(3) &
         ,rws(3,nrws),drws(nrws)
    integer(4) :: irws(nrws)

    pi = 4d0* atan(1d0)
    ci = (0d0,1d0)

    hrotr = (0d0,0d0)

    ir = 0
    do ir = 1,nrws
      do iq = 1,nqbz
        rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
        ceikr = exp(-ci*rk)
        do im = 1,nwf
          do in = 1,nwf
            hrotr(im,in,ir) = hrotr(im,in,ir) + &
                 ceikr * hrotk(im,in,iq) / dble(nqbz)
          enddo ! in
        enddo ! im
      enddo ! iq
    enddo ! ir

    return
  end subroutine get_hrotr_ws
  !-----------------------------------------------------------------------
  subroutine get_hrotkp_ws(hrotr,rws,drws,irws,q, &
       nwf,nqbz,nrws, &
       hrotkp)

    ! see Ref.[2] eq.26
    implicit integer (i-n)
    implicit real*8(a-h,o-z)

    complex(8) :: hrotr(nwf,nwf,nrws),hrotkp(nwf,nwf), &
         ci,cikr,ceikr,ctmp
    real(8) :: q(3),rws(3,nrws),drws(nrws)
    integer(4) :: irws(nrws)

    pi = 4d0* atan(1d0)
    ci = (0d0,1d0)

    hrotkp = (0d0,0d0)

    do ir = 1,nrws
      rk = sum(rws(:,ir)*q(:))
      cikr = ci * 2d0 * pi * rk
      ceikr = exp(cikr) / dble(irws(ir))
      do im = 1,nwf
        do in = 1,nwf
          hrotkp(im,in) = hrotkp(im,in) + &
               ceikr * hrotr(im,in,ir)
        enddo
      enddo
    enddo

    return
  end subroutine get_hrotkp_ws
  !-----------------------------------------------------------------------
  subroutine get_hrotkp_tb_ws(rcut,plat,alat, &
       hrotr,rws,drws,irws,q,  ibasiwf,bas,natom, &
       nwf,nqbz,nrws, &
       hrotkp)

    ! truncate long-range part of hrotr
    ! from get_hrotkp_ws
    ! see Ref.[2] eq.26
    implicit integer (i-n)
    implicit real*8(a-h,o-z)

    parameter (delta = 1d-3)
    complex(8) :: hrotr(nwf,nwf,nrws),hrotkp(nwf,nwf), &
         ci,cikr,ceikr,ctmp
    real(8) :: q(3),rws(3,nrws),drws(nrws),plat(3,3)
    integer(4) :: irws(nrws)

    integer:: ibasiwf(nwf),natom
    real(8):: bas(3,natom)
    !!
    ! m
    !      rc = rcut + delta
    rc = rcut

    pi = 4d0* atan(1d0)
    ci = (0d0,1d0)

    hrotkp = (0d0,0d0)

    do ir = 1,nrws
      ceikr = (0d0,0d0)
      !         rtmp = alat*dsqrt(sum(rws(:,ir)**2))
      !         if (rtmp.le.rc) then
      rk = sum(rws(:,ir)*q(:))
      cikr = ci * 2d0 * pi * rk
      ceikr = exp(cikr) / dble(irws(ir))
      !         endif
      do im = 1,nwf
        imp= ibasiwf(im)
        do in = 1,nwf
          inp= ibasiwf(in)
          !           rtmp = alat*dsqrt( sum( (rws(:,ir)+bas(:,imp)-bas(:,inp))**2 ) )
          !! Rn-0m
          rtmp = alat*dsqrt( sum( (rws(:,ir)+bas(:,inp)-bas(:,imp))**2 ) )
          if(rtmp<rc) then
            hrotkp(im,in) = hrotkp(im,in) + &
                 ceikr * hrotr(im,in,ir)
          endif
        enddo
      enddo
    enddo

    return
  end subroutine get_hrotkp_tb_ws
  !-----------------------------------------------------------------------
  subroutine wmaxloc_diag(ifmlw,ifmlwe, &
       iq,q,umnk,cnk,eunk,evecc,eval, &
       iko_ix,iko_fx,iko_i,iko_f, &
       nwf,nqbz)
    implicit integer (i-n)
    implicit real*8(a-h,o-z)

    complex(8),allocatable :: dnk(:,:)
    complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz), &
         umnk(nwf,nwf,nqbz),evecc(nwf,nwf)
    real(8) :: q(3),eunk(nwf,nqbz),eval(nwf)
    integer(4) :: iko_i(nqbz),iko_f(nqbz)


    write(ifmlw)iq,q(1:3)
    write(ifmlw)evecc(1:nwf,1:nwf)

    write(ifmlwe)iq,q(1:3)
    write(ifmlwe)eval(1:nwf)

    return
  end subroutine wmaxloc_diag
  !-----------------------------------------------------------------------
  ! subroutine chk_diag(q,ham,nwf,eval,evec)
  !   implicit real*8(a-h,o-z)
  !   implicit integer (i-n)
  !   complex(8) :: ham(nwf,nwf),evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf)
  !   real(8) :: eval(nwf),q(3)


  !   write(7100,*)'***'
  !   write(7100,"(3f12.6)")q

  !   ! zm1(i1,iwf2)   = S[i2] ham(i1,i2) * evec(i2,iwf2)
  !   zm1 = 0d0
  !   do iwf2 = 1,nwf
  !      do i1   = 1,nwf
  !         do i2 = 1,nwf
  !            zm1(i1,iwf2) = zm1(i1,iwf2) + ham(i1,i2) * evec(i2,iwf2)
  !         enddo
  !      enddo
  !   enddo

  !   ! zm2(iwf1,iwf2) = S[i1] conjg(evec(i1,iwf1)) * zm1(i1,iwf2)
  !   zm2 = 0d0
  !   do iwf1 = 1,nwf
  !      do iwf2 = 1,nwf
  !         do i1 = 1,nwf
  !            zm2(iwf1,iwf2) = zm2(iwf1,iwf2) &
  !                 + dconjg(evec(i1,iwf1)) * zm1(i1,iwf2)
  !         enddo
  !      enddo
  !   enddo

  !   ! output
  !   do iwf1 = 1,nwf
  !      do iwf2 = 1,nwf
  !         tmp = 0d0
  !         if (iwf1 == iwf2) tmp = eval(iwf1)
  !         write(7100,"(2i5,3f12.6)")iwf1,iwf2, &
  !              dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2)),tmp
  !      enddo
  !   enddo

  !   return
  ! end subroutine chk_diag
  !-----------------------------------------------------------------------
  ! subroutine chk_umnk(q,ham,nwf,eval,umn)
  !   implicit real*8(a-h,o-z)
  !   implicit integer (i-n)
  !   complex(8) :: ham(nwf,nwf),evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf), &
  !        umn(nwf,nwf)
  !   real(8) :: eval(nwf),q(3)


  !   do i1 = 1,nwf
  !      do i2 = 1,nwf
  !         evec(i1,i2) = dconjg(umn(i2,i1))
  !      enddo
  !   enddo

  !   write(7000,*)'***'
  !   write(7000,"(3f12.6)")q

  !   ! zm1(i1,iwf2)   = S[i2] ham(i1,i2) * evec(i2,iwf2)
  !   zm1 = 0d0
  !   do iwf2 = 1,nwf
  !      do i1   = 1,nwf
  !         do i2 = 1,nwf
  !            zm1(i1,iwf2) = zm1(i1,iwf2) + ham(i1,i2) * evec(i2,iwf2)
  !         enddo
  !      enddo
  !   enddo

  !   ! zm2(iwf1,iwf2) = S[i1] conjg(evec(i1,iwf1)) * zm1(i1,iwf2)
  !   zm2 = 0d0
  !   do iwf1 = 1,nwf
  !      do iwf2 = 1,nwf
  !         do i1 = 1,nwf
  !            zm2(iwf1,iwf2) = zm2(iwf1,iwf2) &
  !                 + dconjg(evec(i1,iwf1)) * zm1(i1,iwf2)
  !         enddo
  !      enddo
  !   enddo

  !   ! output
  !   do iwf1 = 1,nwf
  !      do iwf2 = 1,nwf
  !         tmp = 0d0
  !         if (iwf1 == iwf2) tmp = eval(iwf1)
  !         write(7000,"(2i5,3f12.6)")iwf1,iwf2, &
  !              dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2)),tmp
  !      enddo
  !   enddo

  !   return
  ! end subroutine chk_umnk
  !-----------------------------------------------------------------------
  ! subroutine cmp_umn_evec(q,umn,evec,eval,nwf)
  !   implicit real*8(a-h,o-z)
  !   implicit integer (i-n)
  !   complex(8) :: evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf), &
  !        umn(nwf,nwf)
  !   real(8) :: eval(nwf),q(3)


  !   do i1 = 1,nwf
  !      do i2 = 1,nwf
  !         zm1(i1,i2) = dconjg(umn(i2,i1))
  !      enddo
  !   enddo

  !   zm2 = zm1 - evec

  !   ! output
  !   write(7000,*)'***'
  !   write(7000,"(3f12.6)")q
  !   do iwf = 1,nwf
  !      write(7000,*)iwf,eval(iwf)
  !   enddo
  !   do iwf1 = 1,nwf
  !      do iwf2 = 1,nwf
  !         write(7000,"(2i5,3f12.6)")iwf1,iwf2, &
  !              dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2))
  !      enddo
  !   enddo

  !   return
  ! end subroutine cmp_umn_evec
  !-----------------------------------------------------------------------
  subroutine writeham(ifi,is,ef,alat,plat,pos,qbz,wbz,rws,irws,hrotk, &
       nspin,natom,nwf,nqbz,nrws)
    implicit none
    integer(4) :: nspin,nwf,natom,nqbz,nrws,is,irws(nrws)
    double precision :: ef,alat,plat(3,3),pos(3,natom),qbz(3,nqbz), &
         wbz(nqbz),rws(3,nrws)
    complex(8) :: hrotk(nwf,nwf,nqbz)

    integer(4) :: ifi

    if (is == 1) then
      !    ifi = iopen('HMLWF',1,-1,0)
      open(newunit=ifi,file='HMLWF')
      write(ifi,*)nspin,natom,nwf,nqbz,nrws,ef
      write(ifi,*)alat
      write(ifi,*)plat
      write(ifi,*)pos
      write(ifi,*)qbz
      write(ifi,*)wbz
      write(ifi,*)rws
      write(ifi,*)irws
    endif

    write(ifi,*)hrotk
    close(ifi)
    !  if (is == nspin) ifi = iclose('HMLWF')

    return
  end subroutine writeham
  !-----------------------------------------------------------------------

  !  write *rws and hrotr to ifh

  subroutine write_hrotr(ifh, hrotr, &
       rws,irws,drws, &
       nwf,nrws )
    implicit none
    complex(8),intent(in) :: hrotr(nwf,nwf,nrws)
    integer,intent(in):: ifh, nwf,nrws
    real(8),intent(in) :: rws(3,nrws),drws(nrws)
    integer,intent(in) :: irws(nrws)

    integer:: i,j,k
    integer:: ir,im,in
    real(8):: rtmp,heps2

    write(ifh,*) 'nwf ',nwf
    write(ifh,*) 'nrws ',nrws

    write(ifh,*) '<rws>'
    do i=1,nrws
      write(ifh,*) i, rws(:,i), drws(i) , irws(i)
    enddo
    write(ifh,*)'</rws>'

    write(ifh,*) '<hrotr>'
    do i=1,nwf
      do j=1,nwf
        write(ifh,*)  i,j, 'i,j , the next line is hrotr(i,j,:)'
        write(ifh,'(10E20.10)')   hrotr(i,j,:)
      enddo
    enddo
    write(ifh,*) '</hrotr>'

    write(ifh,*) '<hrotr.abs>'
    do i=1,nwf
      do j=1,nwf
        write(ifh,*)  i,j, 'i,j , the next line is hrotr(i,j,:)'
        write(ifh,'(10E20.10)')  ( abs(hrotr(i,j,k)),k=1,nrws)
      enddo
    enddo
    write(ifh,*) '</hrotr.abs>'


  end subroutine write_hrotr


  subroutine read_hrotr(filename,nwf,nrws, &
       hrotr)
    use m_keyvalue,only: getkeyvalue

    implicit none
    character(*),intent(in):: filename
    integer,intent(in):: nwf,nrws
    complex(8):: hrotr(nwf,nwf,nrws)

    integer:: nwf_, nrws_
    character(10):: thisfunc='read_hrotr'
    integer:: ierror, ifh
    integer:: i,j,i_,j_
    character(120):: str

    write(*,*) 'reading ',filename
    call getkeyvalue(filename,'nwf',nwf_)
    call getkeyvalue(filename,'nrws',nrws_)

    ierror=0
    if ( nwf /= nwf_) then
      write(*,*) thisfunc,': data inconsistent nwf=', nwf, ' nwf(file)=',nwf_
      ierror=ierror+1
    endif

    if ( nrws /= nrws_ ) then
      write(*,*) thisfunc,': data inconsistent nrws=', nrws, ' nrws(file)=',nrws_
      ierror=ierror+1
    endif

    if (ierror /= 0) then
      goto 999
    endif

    call getkeyvalue(filename,'<hrotr>',unit=ifh,status=ierror,errstop='on')
    write(*,*) 'ifh,ierror=',ifh,ierror
    if (ierror == 0) then
      write(*,*) thisfunc,': failed to read <hrotr>'
      goto 999
    endif

    do i=1,nwf
      do j=1,nwf
        read(ifh,'(a120)')  str
        write(*,*) 'str=',str(:len_trim(str))
        read(str,*)i_,j_
        write(*,*) '1)',i_,j_
        read(ifh,'(10E20.10)')   hrotr(i,j,1:nrws)
        !          read(ifh,*)   hrotr(i,j,1:nrws)
        write(*,*) '2)',i_,j_
      enddo
    enddo

    close(ifh)

    return
999 write(*,*) 'abnormal exit'
    stop 'in read_hrotr'

  end subroutine read_hrotr


  subroutine make_hrotrcut( hrotr, &
       rws,irws,drws, &
       rcut,heps, &
       nwf,nrws, &
       hrotrcut )
    implicit none
    complex(8):: hrotr(nwf,nwf,nrws)
    real(8):: rws(3,nrws),drws(nrws)
    integer:: irws(nrws)
    real(8):: rcut, heps
    integer:: nwf,nrws
    complex(8),intent(out):: hrotrcut(nwf,nwf,nrws)
    integer:: ir,im,in
    real(8):: heps2,rtmp
    heps2=heps*heps
    hrotrcut=hrotr
    do ir = 1,nrws
      rtmp = dsqrt(sum(rws(:,ir)**2))  ! unit of alat
      write(*,"('cut:',i5,2f10.3)") ir,rtmp,rcut
      if (rtmp>rcut) then
        hrotrcut(:,:,ir)=0.0d0
      endif
      do im = 1,nwf
        do in = 1,nwf
          if ( im == in ) continue
          write(*,"('  ',2i5,'(',d13.5,',',d13.5,')',d13.5,d13.5)") &
               im,in,hrotr(im,in,ir) ,dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2, heps2
          if ( dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2 < heps2 ) then
            hrotrcut(im,in,ir)=0.0d0
          endif

        enddo
      enddo
    enddo
  end subroutine make_hrotrcut





  subroutine cart_to_fract (cart, fract_coord, qlat)
    implicit none
    real(8) :: qlat(3,3)
    real(8) :: cart(3), fract_coord(3)
    integer :: i
    do i=1,3
      fract_coord(i) = sum(cart(1:3)*qlat(1:3,i))
    enddo

    return

  end subroutine cart_to_fract


  subroutine write_hopping_output(is, hrotr, &
       rws,irws,alat,plat,qlat,pos,natom, &
       ibasiwf, nwf,nrws, spid, m_indx, l_indx, &
       nphix, iphi, ldim2)
    implicit none
    integer:: natom, is
    real(8) :: alat,plat(3,3),pos(3,natom),qlat(3,3)
    real(8),allocatable :: cart_coord(:,:), fract_coord(:,:)
    integer:: i,j,k, ldim2
    complex(8),intent(in) :: hrotr(nwf,nwf,nrws)
    integer,intent(in):: nwf,nrws
    real(8),intent(in) :: rws(3,nrws)
    integer,intent(in) :: irws(nrws)

    integer(4) :: ibasiwf(nwf), nphix, iphi(nphix,nwf)
    integer :: quantum_l, quantum_sym
    character(4) :: quantum_n, spin
    character(9), dimension(0:3,-5:5) :: orbital_sym
    character(8) :: spid(natom)
    integer(4) ::  m_indx(ldim2), l_indx(ldim2)
    integer:: ir, iwf, iwf1, iwf2, ifh1,ifh

    open(newunit=ifh, file='Hopping.dat')
    if(is==1)open(newunit=ifh1,file='Hopping.up')
    if(is==2)open(newunit=ifh1,file='Hopping.dn')

    write(ifh,"('Name : Need to be modified')")
    write(ifh1,"('Name : Need to be modified')")
    allocate (cart_coord(3,3))
    do i=1,3
      cart_coord(1:3,i) = alat*plat(1:3,i)*0.529177249
      write(ifh,"(3F16.9)") cart_coord(1:3,i)
      write(ifh1,"(3F16.9)") cart_coord(1:3,i)
    enddo
    write(ifh,"(2I6,I12)") nwf,nrws,nwf*nwf*nrws
    write(ifh,"(1A,I5)") "spin",is
    write(ifh1,"(2I6,I12)") nwf,nrws,nwf*nwf*nrws
    write(ifh1,"(1A,I5)") "spin",is
    deallocate (cart_coord)

    write(ifh,"(9(1A,7X))") &
         "leg","atom","n","l","sym","spin","x","y","z"
    write(ifh1,"(9(1A,7X))") &
         "leg","atom","n","l","sym","  spin","x  ","y  ","z  "


    orbital_sym(0,0)="s"
    orbital_sym(1,-1)="y"
    orbital_sym(1,0)="z"
    orbital_sym(1,1)="x"
    orbital_sym(2,-2)="xy"
    orbital_sym(2,-1)="yz"
    orbital_sym(2,0)="z2"
    orbital_sym(2,1)="xz"
    orbital_sym(2,2)="x2y2"
    orbital_sym(3,-3)="y(3x2-y2)" !see wikipedia!
    orbital_sym(3,-2)="xyz"
    orbital_sym(3,-1)="y(5z2-1)"
    orbital_sym(3,0)="z(5z2-1)"
    orbital_sym(3,1)="x(5z2-1)"
    orbital_sym(3,2)="z(x2-y2)"
    orbital_sym(3,3)="x(x2-3y2)"

    quantum_n = "--"

    allocate (fract_coord(3,nwf))
    do iwf=1,nwf
      if (is == 1) spin = "up"
      if (is == 2) spin = "down"

      call cart_to_fract(pos(1:3,ibasiwf(iwf)), fract_coord(1:3,iwf), &
           qlat)

      write(ifh,"(I3,1A16,A6,I6,2A12,3F10.5)") &
           iwf, spid(ibasiwf(iwf)), quantum_n, l_indx(iphi(1,iwf)), &
           orbital_sym(l_indx(iphi(1,iwf)),m_indx(iphi(1,iwf))), &
           spin, fract_coord(1:3,iwf)

      write(ifh1,"(I3,1A16,A6,I6,2A12,3F10.5)") &
           iwf, spid(ibasiwf(iwf)), quantum_n, l_indx(iphi(1,iwf)), &
           orbital_sym(l_indx(iphi(1,iwf)),m_indx(iphi(1,iwf))), &
           spin, fract_coord(1:3,iwf)

    end do
    deallocate (fract_coord)

    allocate (fract_coord(3,nrws))
    do ir=1,nrws
      call cart_to_fract(rws(:,ir), fract_coord(1:3,ir), &
           qlat)

      do iwf1=1,nwf
        do iwf2=1,nwf
          write(ifh,"(1x, 3i6, 3f12.6, 2x,2i4,3f12.6)") &
               nint(fract_coord(:,ir)), rws(:,ir), iwf1, iwf2, &
               hrotr(iwf1,iwf2,ir)*13.605698066
          write(ifh1,"(1x, 3i6, 3f12.6, 2x,2i4,3f12.6)") &
               nint(fract_coord(:,ir)), rws(:,ir), iwf1, iwf2, &
               hrotr(iwf1,iwf2,ir)*13.605698066

        enddo
      enddo
      !        do iwf1=1,nwf
      !          do iwf2=1,nwf
      !            write(ifh,"(1x, 6f12.6, 2i4,3f12.6)")
      !     &             rws(:,ir), fract_coord(:,ir), iwf1, iwf2,  hrotr(iwf1,iwf2,ir)
      !          enddo
      !        enddo
    enddo
    deallocate (fract_coord)
    close(ifh1)
    close(ifh)
  end subroutine write_hopping_output
end module m_maxloc0


subroutine readuu(is,iti,itf,ikbidx, &
       nqbz,nbb, &
       uum)
    implicit integer (i-n)
    implicit real*8(a-h,o-z)
    integer :: ikbidx(nbb,nqbz)
    complex(8) :: uum(iti:itf,iti:itf,nbb,nqbz)
    complex(8),allocatable::uumread(:,:)
    if (is == 1) then
      !    ifuu      = iopen('UUU',0,0,0)
      open(newunit=ifuu,file='UUU',form='unformatted')
    else
      !     ifuu      = iopen('UUD',0,0,0)
      open(newunit=ifuu,file='UUD',form='unformatted')
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
    if (itf2 < itf)   call rxii('readuu: itf error',itf2,itf)
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
    close(ifu)
    !if (is == 1) then
    !   close(ifu)! = iclose('UUU')
    !else
    !   close(ifu)! = iclose('UUD')
    !endif

    return
end subroutine readuu
