subroutine getgv2(alat,plat,qlat,q, QpGcut,job,    ng, ngvec)
  !-  Set up a list of recip vectors within cutoff |Q+G| < QpGcut a.u.
  ! takao
  ! job==1 -> return ng (number of G ) only
  ! job==2 -> return ng and ngvec
  ! True G is given by
  !    G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
  implicit none
  integer :: ng
  real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qpg(3),enor(3) &
       ,pi,alat,tpiba,QpGmax,QpGmax2,QpGcut,Qenor
  integer :: &
       n1max,n1min,n2max,n2min,n3max,n3min, i1,i2,i3,ig,job,imx
  integer ::  ngvec(3,ng)
  pi=4d0*datan(1d0)
  tpiba=2*pi/alat
  QpGmax   = QpGcut/tpiba  ! QpGcut in a.u.= tpiba*Qcut
  QpGmax2  = QpGmax**2
  call eprod(qlat(1:3,2),qlat(1:3,3),enor)
  Qenor = sum(qlat(1:3,1)*enor)
  n1max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
  n1min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  call eprod(qlat(1:3,3),qlat(1:3,1),enor)
  Qenor = sum(qlat(1:3,2)*enor)
  n2max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
  n2min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  call eprod(qlat(1:3,1),qlat(1:3,2),enor)
  Qenor = sum(qlat(1:3,3)*enor)
  n3max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
  n3min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  ! get ngvec within the limit.
  ig=0
  imx=-9999
  do i1 = n1min, n1max
     do i2 = n2min, n2max
        do i3 = n3min, n3max
           qpg(1:3)= q(1:3) + &
                qlat(1:3,1)*i1 +qlat(1:3,2)*i2 +qlat(1:3,3)*i3
           if( sum(qpg(1:3)**2) < QpGmax2) then
              ig = ig+1
              if(job==2) ngvec(1:3,ig) = (/i1,i2,i3/)
              if(job==1) imx=max(imx,abs(i1),abs(i2),abs(i3))
              !    if(job==2) write(1116,'(f8.4,3i3)') tpiba*sqrt(sum(qpg(1:3)**2)),ngvec(1:3,ig) ! check write
           endif
        enddo
     enddo
  enddo
  ng = ig
  if(job==1) ngvec(1,1)=imx !mar2012takao
end subroutine getgv2
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine eprod(a,b, c)
  ! c gives normalized normal vector for a and b.
  real(8) :: a(3),b(3),c(3),cnorm !anorm,bnorm,
  c(1)= a(2)*b(3)-a(3)*b(2)
  c(2)= a(3)*b(1)-a(1)*b(3)
  c(3)= a(1)*b(2)-a(2)*b(1)
  cnorm = sqrt(sum(c(1:3)**2))
  c = c/cnorm
end subroutine eprod
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine sortea(ea,ieaord,n,isig)
  ! mini-sort routine.
  implicit real*8(a-h,o-z)
  real(8)::        ea(n)
  ! ino delete integer(4) def.      integer(4):: ieaord(n)
  integer:: ieaord(n),n,isig,itmp,i,ix
  ! sorting of ea
  isig = 1
  do i = 1,n
     ieaord(i) = i
  enddo
  do ix= 2,n
     do i=ix,2,-1
        if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
           itmp = ieaord(i-1)
           ieaord(i-1) = ieaord(i)
           ieaord(i) = itmp
           isig= -isig
           cycle
        endif
        exit
     enddo
  enddo
end subroutine sortea

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine pairtakao_unuused(job,pos,nbas,plat,nk,npair, npairmx,range, nlat,qwgt)
  ! Obtain pair table of atomic sites.
  !    number of pairs are larger than nk. range(ib1,ib2) is calculated.
  ! job=1
  !i   pos,nbas,plat,nk
  !o   npairmx,range
  ! job=2
  !i   pos,nbas,plot,nk,npairmx,range.
  !o   npair(rewritten), nlat,qwgt.
  !r   qwgt is the weight for integration. Total number of qwgt is nk =nk1*nk2*nk3.

  implicit none
  integer::nbas,job,nk,npair(nbas,nbas),npairmx,ibas,ib1,ib2,ni,isig,nnn,iwend,nend,iwinit
  real(8):: pos(3,nbas),plat(3,3), xx,pi,qlat(3,3),q(3),rrr(3),rx
  integer :: nlat(3,npairmx,nbas,nbas)
  real(8)::     qwgt(npairmx,nbas,nbas)
  real(8),allocatable:: rr(:,:,:),rrx(:)
  integer,allocatable:: iord(:),nlatx(:,:)
  real(8):: range(nbas,nbas),rangex
  !      logical:: exiton
  pi=4d0*datan(1d0)
  call dinv33(plat,1,qlat,xx)
  print *,'job=',job,nk,npairmx
  if(job==1) then !determine range,and npairmx
     npair=0
     do ib1=1,nbas
        do ib2=1,nbas
           !            exiton=.false.
           rangex=0d0
           do
              rangex= rangex+.2d0 !We will improve this to give some better guess.
              call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), rangex,job, &
                   npair(ib1,ib2), nlat(:,:,ib1,ib2) )
              if(npair(ib1,ib2)>=nk) then
                 rangex= rangex+.1d0 !We will improve this to give some better guess.
                 call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), rangex,job, &
                      npair(ib1,ib2), nlat(:,:,ib1,ib2) )
                 range(ib1,ib2)=rangex ! rangex + 0.01d0 is to make degeneracy safe.
                 print *,' ib1 ib2 rangex nk npair',ib1,ib2,rangex,nk,npair(ib1,ib2)
                 exit
              endif
           enddo
        enddo
     enddo
     npairmx = maxval(npair)
     npair=-99999 !not the output
     return
  elseif(job==2) then
     allocate(rr(npairmx,nbas,nbas))
     do ib1=1,nbas
        do ib2=1,nbas
           call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), range(ib1,ib2),job, &
                npair(ib1,ib2), nlat(:,:,ib1,ib2) )
           nnn = npair(ib1,ib2)
           do ni=1,nnn
              rr(ni,ib1,ib2)= sqrt(sum((pos(:,ib1)-pos(:,ib2)+ matmul(plat, nlat(:,ni,ib1,ib2)))**2))
           enddo
           ! sort
           allocate(iord(nnn),nlatx(3,nnn),rrx(nnn))
           call sortea(rr(:,ib1,ib2),iord,nnn,isig)
           nlatx = nlat(:,1:nnn,ib1,ib2)
           rrx   = rr(1:nnn,ib1,ib2)
           do ni = 1,nnn
              nlat(:,ni,ib1,ib2) = nlatx(:,iord(ni))
              rr(ni,ib1,ib2)     = rrx(iord(ni))
           enddo
           deallocate(iord,nlatx,rrx)
           ! weight
           rx = rr(nk,ib1,ib2)
           iwinit=-99
           iwend=iwinit-1
           qwgt(:,ib1,ib2)=0d0
           do ni=1,nnn
              if(abs(rr(ni,ib1,ib2)-rx)< 1d-6 .AND. iwinit==-99) iwinit=ni
              if(abs(rr(ni,ib1,ib2)-rx)< 1d-6) iwend=ni
           enddo
           qwgt(1:iwinit-1,ib1,ib2)=1d0
           qwgt(iwinit:iwend,ib1,ib2)=dble(nk-iwinit+1)/(iwend-iwinit+1)
           npair(ib1,ib2)=iwend
           qwgt(:,ib1,ib2)=qwgt(:,ib1,ib2)/nk
           if(abs(sum(qwgt(:,ib1,ib2))-1d0)>1d-8) then
              stop 'paritakao: abs(sum(qwgt(:,ib1,ib2))-nk)/=0'
           endif
           ! check write
           write(6,"(a,2i4,2x,i6)") ' ib1 ib2 npair=',ib1,ib2,npair(ib1,ib2)
           do ni = 1,npair(ib1,ib2)
              write(6,"(i6,3x,3i3,2f8.3)") ni,nlat(1:3,ni,ib1,ib2),rr(ni,ib1,ib2),qwgt(ni,ib1,ib2)
           enddo
           !          do ni = npair(ib1,ib2)+1,nnn
           !            write(6,"(i6,3x,3i3,2f8.3)") ni,nlat(1:3,ni,ib1,ib2),rr(ni,ib1,ib2)
           !          enddo
        enddo
     enddo
     deallocate(rr)
  endif
end subroutine pairtakao_unuused

