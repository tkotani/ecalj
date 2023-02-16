subroutine get_nwf_MPI(nwf)
  use rsmpi !RS
  implicit none
  integer(4):: nwf,nqbz,iko_ix,iko_fx,ifi
  integer(4):: iopen,iclose

  if (Is_IO_Root_RSMPI()) then
     ifi  = iopen('MLWU',0,0,0)
     read(ifi)nqbz,nwf,iko_ix,iko_fx
     ifi = iclose('MLWU')
  endif
  call MPI_Bcast(nwf,1,MPI_INTEGER,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(nwf)",ierror_rsmpi)
end subroutine get_nwf_MPI
!===================================================================
subroutine choose_wanband_MPI(iwbnd,nwf,nqbze,nspin)
  use rsmpi !RS
  implicit none
  integer(4):: nwf,nqbze,nspin,iwbnd(nwf,nqbze,nspin)
  integer(4):: iopen,iclose
  ! local
  integer(4):: is,ifmlw,ifuu,nq0i,nwf2,nqbz,nqbz2,iko_ix,iko_fx, &
       iko_ix2,iko_fx2,iqbz,iqbz2,iq0i,iq0i2,j1,j2, &
       ikp,ib,iwf,nbnd
  integer(4),allocatable:: isort(:)
  real(8):: q(3)
  real(8),allocatable:: wbnd(:)
  complex(8),allocatable:: dnk(:,:,:,:),uum(:,:,:,:,:),cbwf(:,:)

  if (Is_IO_Root_RSMPI()) then
     ! --- Readin MLWU/D and UUq0U/D
     do is = 1,nspin

        ! fileopen
        if (is == 1) then
           ifmlw  = iopen('MLWU',0,0,0)
           ifuu   = iopen('UUq0U',0,0,0)
        else ! is
           ifmlw  = iopen('MLWD',0,0,0)
           ifuu   = iopen('UUq0D',0,0,0)
        endif ! is

        ! nqbz mesh-points
        read(ifmlw)nqbz,nwf2,iko_ix,iko_fx
        if (nwf2 /= nwf) &
             call RSMPI_Stop("choose_wanband: nwf error")
        !         if (nqbze.ne.nqbz)
        !     >      call RSMPI_Stop("choose_wanband: nqbz error")
        if (is == 1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nspin))
        do iqbz = 1,nqbz
           read(ifmlw)iqbz2,q(1:3)
           if (iqbz2 /= iqbz) &
                call RSMPI_Stop('choose_wanband: iqbz error')
           read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo ! iqbz

        ! shifted mesh points
        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is == 1) &
             allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nspin))
        if (nqbz2 /= nqbz) &
             call RSMPI_Stop("choose_wanband: nqbz2 error")
        if (iko_ix2 /= iko_ix) &
             call RSMPI_Stop("choose_wanband: iko_ix2 error")
        if (iko_fx2 /= iko_fx) &
             call RSMPI_Stop("choose_wanband: iko_fx2 error")
        do iqbz = 1,nqbz
           do iq0i =1,nq0i
              read(ifuu)
              read(ifuu)iqbz2,iq0i2
              if (iqbz2 /= iqbz) &
                   call RSMPI_Stop('choose_wanband: iqbz error')
              if (iq0i2 /= iq0i) &
                   call RSMPI_Stop('choose_wanband: iq0i error')
              read(ifuu) &
                   ((uum(j1,j2,iqbz,iq0i,is), &
                   j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
           enddo ! iq0i
        enddo ! iqbz

        ! fileclose
        if (is == 1) then
           ifmlw  = iclose('MLWU')
           ifuu   = iclose('UUq0U')
        else ! is
           ifmlw  = iclose('MLWD')
           ifuu   = iclose('UUq0D')
        endif ! is

     enddo ! is

     allocate(cbwf(iko_ix:iko_fx,nwf),wbnd(iko_ix:iko_fx), &
          isort(iko_fx-iko_ix+1))
     do ikp = 1,nqbze
        iqbz = mod(ikp,nqbz)
        if (iqbz == 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is = 1,nspin
           cbwf = 0d0
           if (iq0i == 0) then
              cbwf(:,:) = dnk(:,:,iqbz,is)
           else ! iq0i
              !   <psi(k+q0,n) | psi(k+q0,m)^W>
              ! = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
              !      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^W>
              ! ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^W>

              ! psi^W : bloch fn. in the Wannier gauge
              do ib = iko_ix,iko_fx
                 do iwf= 1,nwf
                    cbwf(ib,iwf) = &
                         sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is)) &
                         *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
                 enddo ! iwf
              enddo ! ib
           endif ! iq0i

           ! choose bands
           wbnd = 0.0d0
           nbnd = iko_fx - iko_ix + 1
           do ib = iko_ix,iko_fx
              do iwf = 1,nwf
                 wbnd(ib) = wbnd(ib) + abs(cbwf(ib,iwf))**2
              enddo ! iwf
           enddo ! ib
           call sortbnd(wbnd(iko_ix:iko_fx),iko_fx-iko_ix+1, &
                isort)
           do iwf = 1,nwf
              iwbnd(iwf,ikp,is) = isort(iwf) + iko_ix - 1
           enddo ! iwf
        enddo ! is
     enddo ! ikp
     deallocate(dnk,uum,cbwf,wbnd,isort)

  endif ! Is_IO_Root_RSMPI

  call MPI_Bcast(iwbnd,nwf*nqbze*nspin,MPI_INTEGER,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(iwbnd)",ierror_rsmpi)

end subroutine choose_wanband_MPI
!===================================================================
subroutine sortbnd(rin,n,isort)
  implicit none
  integer(4):: n,isort(n)
  integer(4):: i,j,itmp
  real(8):: rin(n),r(n),rtmp

  do j = 1,n
     isort(j) = j
  enddo ! j

  r(:) = rin(:)
  do j = 2,n
     rtmp = r(j)
     do i = j-1,1,-1
        if (r(i) >= rtmp) goto 999
        r(i+1) = r(i)
        isort(i+1) = isort(i)
     enddo ! i
     i = 0
999  continue
     r(i+1) = rtmp
     isort(i+1) = j
  enddo ! j

  do j = 1,n-1
     if (rin(isort(j)) < rin(isort(j+1))) &
          stop "hx0fp0: sortbnd error"
  enddo

end subroutine sortbnd
!===================================================================
subroutine sortr(a,n,isort)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  real(8) :: eps=1.0d-6
  real(8) :: a(n),b(n)
  integer:: isort(n)

  b = a
  do i = 1,n
     isort(i) = i
  enddo

  do j = 2,n
     c = b(j)
     do i = j-1,1,-1
        if (b(i)<=c) goto 999
        b(i+1) = b(i)
        isort(i+1) = isort(i)
     enddo
     i = 0
999  b(i+1) = c
     isort(i+1) = j
  enddo

  do i = 1,n-1
     if (b(i) > b(i+1)) stop 'sortr: sorting error!'
  enddo
  do i = 1,n
     if (abs(b(i)-a(isort(i))) > eps) stop 'sortr: sorting error!'
  enddo

  return
end subroutine sortr
!-----------------------------------------------------------------------
subroutine q2q0xx(q,plat,q0)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps=1d-4)
  real(8) :: q(3),q0(3),plat(3,3)

  do ii = 1,3
     q0(ii) = sum(q(:)*plat(:,ii))
  enddo

  return
end subroutine q2q0xx
!-----------------------------------------------------------------------
subroutine q02q(q0,qbas,q)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: q(3),q0(3),qbas(3,3)

  q(:) = 0d0
  do ii = 1,3
     q(:) = q(:) + qbas(:,ii)*q0(ii)
  enddo

  return
end subroutine q02q
!-----------------------------------------------------------------------
subroutine q02q0g0(qin,qout,ng)
  real(8) :: qin(3), qout(3)
  integer :: ng(3)
  ng = nint(qin)
  qout=qin-ng
end subroutine q02q0g0
!-----------------------------------------------------------------------
subroutine read_syml(qbandx,nqbandx,nqband)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4) :: nqbandx,nqband,ifsyml,nline,nlinemax
  integer(4),allocatable :: nqq(:)
  real(8) :: qbandx(3,nqbandx),qqx(3)
  real(8),allocatable :: qq1(:,:),qq2(:,:)

  nlinemax = 50
  allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
  ifsyml = ifile_handle() !3001
  open(ifsyml,file='SYML')
  nline = 0
  do
     nline = nline + 1
     read(ifsyml,*,err=601,end=601) &
          nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
  enddo
601 continue
  close(ifsyml)
  nline = nline - 1
  write(6,"(/' Symmetry lines:'/' points',12x,'start',22x,'end')")
  do is=1,nline
     write(6,"(i6,2x,3f8.4,2x,3f8.4)") &
          nqq(is),(qq1(i,is),i=1,3),(qq2(i,is),i=1,3)
  enddo
  nqnumx = sum(nqq(1:nline))
  iqq = 0
  do is = 1,nline
     nk = nqq(is)
     do iq=1,nk
        xx = 0d0
        if(nk>1) xx=(iq-1d0)/(nk-1d0)
        qqx = xx*qq2(1:3,is)+(1d0-xx)*qq1(1:3,is)
        !          if(iqq>1 ) then
        !            if(abs(sum(qqx-qq_rsband(:,iqq)))<1d-10) cycle
        !          endif
        iqq = iqq + 1
        qbandx(1:3,iqq) = qqx
        write (6,"('  q=',3f7.3)") qbandx(1:3,iqq)
     enddo
  enddo
  nqband = iqq
  write (6,"(' Total number of q-points:',i5/)") nqband
  deallocate(nqq,qq1,qq2)
  return
end subroutine read_syml
!-----------------------------------------------------------------------
subroutine findk(q1,qbz,plat,qbas,nqbz,iq2)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  real(8) :: q1(3),qbz(3,nqbz),plat(3,3),qbas(3,3), &
       dqx(3),dqx0(3),dq0(3),dq(3),ddq(nqbz)
  integer(4) :: ndg(3),isort(nqbz)

  do i = 1,nqbz
     dqx = q1 - qbz(:,i)
     call q2q0xx(dqx,plat,dqx0)
     call q02q0g0(dqx0,dq0,ndg)
     call q02q(dq0,qbas,dq)
     ddq(i) = sum(dq*dq)
  enddo
  call sortr(ddq,nqbz,isort)
  iq2 = isort(1)

  return
end subroutine findk
!-----------------------------------------------------------------------
subroutine findk4(qband,qbz,plat,qbas,nqband,nqbz,wqk4,iqk4)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  real(8) :: qband(3,nqband),qbz(3,nqbz),plat(3,3),qbas(3,3),q1(3), &
       wqk4(4,nqband),dqx(3),dqx0(3),dq0(3),dq(3),ddq(nqbz), &
       amat(4,4),amat2(4,4),bvec(4),vol(3,3),qtmp(3,4),ipiv(4)
  integer(4) :: iqk4(4,nqband),ndg(3),isort(nqbz)
  logical :: lddq

  eps = 1.0d-4
  nn = min(nqbz,1000)
  iqk4 = -1
  wqk4 = 0.0d0
  do iq = 1,nqband
     q1 = qband(:,iq)
     nerr = 0
     ! sorting
     do i = 1,nqbz
        dqx = q1 - qbz(:,i)
        call q2q0xx(dqx,plat,dqx0)
        call q02q0g0(dqx0,dq0,ndg)
        call q02q(dq0,qbas,dq)
        ddq(i) = sum(dq*dq)
     enddo
     call sortr(ddq,nqbz,isort)
     iqk4(1:4,iq) = isort(1:4)

     lddq = .false.
     if (ddq(isort(1)) < eps) lddq = .TRUE. 
     if (lddq) then
        wqk4(1,iq) = 1.0d0
        ! main part
     else ! ddq
200     continue
        do ii=1,4
           call q2qg(qbz(:,iqk4(ii,iq)),q1,plat,qbas,qtmp(:,ii))
           do jj=1,3
              amat(ii,jj)=qtmp(jj,ii)
           enddo ! jj
           amat(ii,4) = 1.d0
        enddo ! ii
        !-----------------------------------------------------------------------
        !     Check if the quadruplet of points form a tetrahedron with non-zero
        !     volume

        do ii=1,3
           do jj=1,3
              vol(ii,jj) = amat(ii+1,jj) - amat(1,jj)
           enddo
        enddo
        avol = vol(1,1)*vol(2,2)*vol(3,3) + vol(1,2)*vol(2,3)*vol(3,1) + &
             vol(1,3)*vol(2,1)*vol(3,2) - vol(1,1)*vol(2,3)*vol(3,2) - &
             vol(1,2)*vol(2,1)*vol(3,3) - vol(1,3)*vol(2,2)*vol(3,1)
        if (abs(avol) < eps*eps*eps) then
300        continue
           nerr = nerr + 1
           if (nerr == (nn-3)**3) then
              write(6,*) ' warning: nn parameter in findk4 too small!'
              write(6,"(3f10.5)")q1
              write(6,"(i5,6f10.5)")iqk4(1,iq),qbz(:,iqk4(1,iq)),qtmp(:,1)
              write(6,"(i5,6f10.5)")iqk4(2,iq),qbz(:,iqk4(2,iq)),qtmp(:,2)
              write(6,"(i5,6f10.5)")iqk4(3,iq),qbz(:,iqk4(3,iq)),qtmp(:,3)
              write(6,"(i5,6f10.5)")iqk4(4,iq),qbz(:,iqk4(4,iq)),qtmp(:,4)
              write(6,"('amat',4f12.6)")amat(1:4,1)
              write(6,"('amat',4f12.6)")amat(1:4,2)
              write(6,"('amat',4f12.6)")amat(1:4,3)
              write(6,"('amat',4f12.6)")amat(1:4,4)
              write(6,600) avol,q1,nn,nerr
600           format(4f9.4,2i10)
              write(6,*) (isort(jj),jj=1,nn)

              !     impossible to form a tetrahedron with non-zero volume from the
              !     q-points listed in isort; at this point, the last quadruplet
              !     of q-points has the 1st, (nn-2)-th, (nn-1)-th and nn-th nearest
              !     q-points

              !            goto 500
              write(6,*) 'ERROR: Cannot form tetrahedron'
              write(6,*) '  The program stops...'
              stop
           endif ! nerr

           !     change the 2nd, 3rd and/or 4th q-points in the quadruplet and try
           !     to find a tetrahedron again

           ii = int( nerr/(nn-3)**2 )
           jj = int( (nerr - ii*(nn-3)**2)/(nn-3) )
           kk = nerr - ii*(nn-3)**2 - jj*(nn-3)
           if (ii > jj .OR. jj > kk) goto 300
           iqk4(2,iq) = isort(2 + ii)
           iqk4(3,iq) = isort(3 + jj)
           iqk4(4,iq) = isort(4 + kk)
           goto 200
        endif ! abs(avol)
        !-----------------------------------------------------------------------
        do jj=1,4
           amat2=amat
           bvec=0d0
           bvec(jj)=1d0

           !  call lapack subroutine

           call DGESV(4,1,amat2,4,ipiv,bvec,4,info)
           if (info /= 0) then
              write(6,*) 'failure to determine wqk4. info =',info, &
                   q1(:),iqk4(1:4,iq),nerr
              write(6,*) 'ERROR: interpolation'
              write(6,*) '  The program stops...'
              stop
              !            goto 500
           endif
           asum=0.d0
           do ii=1,3
              asum = asum + bvec(ii)*q1(ii)
           enddo
           wqk4(jj,iq) = asum + bvec(4)
        enddo     ! jj=1,4
500     continue

     endif ! ddq
  enddo ! iq

  ! eck
  do iq = 1,nqband
     wsum = -1.0d0
     do i4 = 1,4
        if (iqk4(i4,iq) < 1 .OR. iqk4(i4,iq) > nqbz) &
             stop 'findk4: iqk4 error'
        wsum = wsum + wqk4(i4,iq)
     enddo ! i4
     if (abs(wsum) > eps) stop 'findk4: wqk4 error'
  enddo ! iq

  return
end subroutine findk4
!-----------------------------------------------------------------------
subroutine q2qg(q,q0,plat,qbas,qg)
  implicit none

  double precision :: q(3),q0(3),qg(3),diff(3),p,eps, &
       plat(3,3),qbas(3,3),qq(3),qq0(3),qqg(3)
  integer :: i

  eps=1d-4

  call q2q0xx(q,plat,qq)
  call q2q0xx(q0,plat,qq0)
  do i=1,3
     qqg(i) = qq(i) - idnint(qq(i)-qq0(i))
  enddo

  diff(:) = qqg(:)-qq0(:)
  do i=1,3
     if (abs(diff(i)) > 0.5+eps) stop 'q2qg error'
  enddo

  call q02q(qqg,qbas,qg)

  return
end subroutine q2qg
!-----------------------------------------------------------------------
subroutine readuuq0(is,iko_ix,iko_fx,nqbz,nq0i,uuq0)
  use rsmpi !RS
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  complex(8) :: uuq0(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i)

  if (is == 1) then
     ifuu      = iopen('UUq0U',0,0,0)
  else
     ifuu      = iopen('UUq0D',0,0,0)
  endif

  read(ifuu)
  read(ifuu)nqbz2,nq0i_2,iko_ix2,iko_fx2

  if (nqbz2 /= nqbz) &
       call RSMPI_Stop("readuuq0: nqbz2 error")
  if (nq0i_2 /= nq0i) &
       call RSMPI_Stop("readuuq0: nq0i_2 error")
  if (iko_ix2 /= iko_ix) &
       call RSMPI_Stop("readuuq0: iko_ix2 error")
  if (iko_fx2 /= iko_fx) &
       call RSMPI_Stop("readuuq0: iko_fx2 error")
  do iqbz = 1,nqbz
     do iq0i =1,nq0i
        read(ifuu)
        read(ifuu)iqbz2,iq0i2
        if (iqbz2 /= iqbz) &
             call RSMPI_Stop('readuuq0: iqbz error')
        if (iq0i2 /= iq0i) &
             call RSMPI_Stop('readuuq0: iq0i error')
        read(ifuu) &
             ((uuq0(j1,j2,iqbz,iq0i), &
             j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
     enddo ! iq0i
  enddo ! iqbz

  if (is == 1) then
     ifu = iclose('UUq0U')
  else
     ifu = iclose('UUq0D')
  endif

  return
end subroutine readuuq0
!-----------------------------------------------------------------------
