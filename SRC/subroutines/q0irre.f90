subroutine auxfunlm(q, alp, alat, qlat, ngc, ngvect,lx,  auxfun,skipg0)
  !!  Give auxfun(L) = \sum_G  exp(-alp*(q+G)**2)*Y_L(q+G)
  implicit none
  real(8) :: alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2 &
       ,pi=3.1415926535897932D0,alp
  integer :: ig, ngc, ngvect(3,ngc),lx
  real(8):: r2s, auxfun((lx+1)**2),alpqg2
  real(8),allocatable:: cy(:),yl(:)
  logical:: skipg0
  allocate(cy((lx+1)**2),yl((lx+1)**2))
  call sylmnc(cy,lx)
  tpiba  = 2d0*pi/alat
  auxfun = 0d0
  do ig = 1,ngc
     if(skipg0 .AND. sum(ngvect(1:3,ig)**2)==0) cycle
     qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
     qg2     = sum(qg(1:3)**2)
     alpqg2= alp* qg2
     call sylm(qg/sqrt(qg2),yl,lx,r2s) !spherical factor Y( q+G )
     auxfun = auxfun + exp(-alpqg2)/qg2 *cy(:)*yl(:) !cy*yl =Y_L(qg/|qg|)
  enddo
end subroutine auxfunlm

subroutine auxfunqg(q, alp, alat, qlat, ngc, ngvect, fout)
  !!  Give auxfun(L) = \sum_G  exp(-alp*(q+G)**2)*Y_L(q+G)
  implicit none
  real(8) :: alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2 &
       ,pi=3.1415926535897932D0,alp,qg2smallest
  integer :: ig, ngc, ngvect(3,ngc)
  real(8):: fout,alpqg2
  tpiba  = 2d0*pi/alat
  qg2smallest=1d10
  do ig = 1,ngc
     qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
     qg2     = sum(qg(1:3)**2)
     if(qg2<qg2smallest) then
        qg2smallest=qg2
     endif
  enddo
  alpqg2= alp* qg2smallest
  fout = exp(-alpqg2)/qg2smallest  !y00 =1/sqrt(4d0*pi)
end subroutine auxfunqg


!! -----------------------------------------------------------------------------
subroutine getwklm(alat,vol,plat,qlat,alp,qbz,nnn,ngc,ngcmx,ngvect,lx,n1q,n2q,n3q, &
     wklm)!,wqfac)
  !! == spherical integration weight Klm with reference auxially functions ==
  !! Output
  !!   wklm: this means K_lm defined around Eq.35 in Copmuter Physics Comm. 176(2007)1-13.
  implicit none
  integer :: iq,i,nnn,lm,lx
  integer :: ngc(nnn), ngcmx, ngvect(3,ngcmx,nnn)
  real(8)    :: alp, alat,qbz(3,nnn),vol,volinv, &
       qlat(3,3) ,plat(3,3),wtrue00, &
       pi=3.1415926535897932D0, wklm((lx+1)**2)
  real(8),allocatable:: funa(:,:),wsumau(:),funac(:)
  logical ::skipg0
  integer:: iii,n1q,n2q,n3q,ndiv,iq1,iq2,iq3,lxx
  !! ----------------

  !$$$      real(8):: fcenter,fmean,fout(1),qmic(3,3),qx(3),fsum,wg,ftot,ftotc
  !$$$      real(8):: wqfac(nnn)
  !$$$!! discrete sum problem. Central value v.s. Mean value for each cell.
  !$$$      ndiv=5
  !$$$      print *,' getwklm n12 n2q n3q 2*ndiv',n1q,n2q,n3q,2*ndiv
  !$$$      skipg0=.false.
  !$$$      lxx=0
  !$$$      ftotc=0d0
  !$$$      ftot=0d0
  !$$$      do iq = 2,nnn ! omit q=0 point iq=1
  !$$$         fsum=0d0
  !$$$         qmic(:,1)=qlat(:,1)/(n1q*ndiv*2)
  !$$$         qmic(:,2)=qlat(:,2)/(n2q*ndiv*2)
  !$$$         qmic(:,3)=qlat(:,3)/(n3q*ndiv*2)
  !$$$         do iq1=-ndiv,ndiv
  !$$$         do iq2=-ndiv,ndiv
  !$$$         do iq3=-ndiv,ndiv
  !$$$c            print *,'iq1 ',iq1,iq2,iq3
  !$$$            qx=matmul(qmic,(/iq1,iq2,iq3/)) + qbz(1:3,iq)
  !$$$c            call auxfunlm(qx, alp,alat, qlat,
  !$$$c     &             ngc(iq),ngvect(1:3,1:ngc(iq),iq),lxx,fout,skipg0)
  !$$$            call auxfunqg(qx, alp,alat, qlat,
  !$$$     &           ngc(iq),ngvect(1:3,1:ngc(iq),iq),fout(1))
  !$$$            wg=1d0
  !$$$            if(mod(iq1+ndiv,2*ndiv)==0) wg=0.5d0
  !$$$            if(mod(iq2+ndiv,2*ndiv)==0) wg=wg*0.5d0
  !$$$            if(mod(iq3+ndiv,2*ndiv)==0) wg=wg*0.5d0
  !$$$            fsum = fsum + wg*fout(1)
  !$$$            if(iq1==0 .and.iq2==0.and.iq3==0 ) fcenter=fout(1)
  !$$$         enddo
  !$$$         enddo
  !$$$         enddo
  !$$$         fmean=fsum/(8d0*ndiv**3)
  !$$$         ftot = ftot+fmean
  !$$$         ftotc= ftotc+fcenter
  !$$$         wqfac(iq) =fmean/fcenter
  !$$$         write(*,"('iq q fcenter fmean ratio=',i4,3f7.3,2x,3f8.4)")
  !$$$     &   iq,qbz(1:3,iq),fcenter,fmean,fmean/fcenter
  !$$$      enddo
  !$$$      write(*,"('ftot ftotc ftot/ftotc: ',2f13.4)") ftot/(nnn-1),ftotc/(nnn-1),ftot/ftotc


  allocate(funa((lx+1)**2,nnn),wsumau((lx+1)**2),funac((lx+1)**2))
  !! ==== true integal of auxially function ====
  !! integral \int d^3k/vol * exp(-alp*k**2) *Y_0
  volinv  = (2*pi)**3/vol
  !     wtrue00 = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alp**3) /sqrt(4d0*pi) ! 1d0/sqrt(4pi) is Y00.
  wtrue00 = 4d0*pi/volinv *sqrt(pi)/2d0/sqrt(alp)/sqrt(4d0*pi)!1/sqrt(4pi)=Y00. !bugfix 19nov2012. alp

  !! ==== discrete sum (except q=0) of auxially functions ====
  ! constant part
  skipg0=.true.
  iq=1
  call auxfunlm((/0d0,0d0,0d0/), alp,alat, qlat, &
       ngc(iq),ngvect(1:3,1:ngc(iq),iq),lx,funac(:),skipg0)
  !      do lm=1,(lx+1)**2
  !        if(abs(funac(lm))>1d-6) write(6,"('  getwklm: lm funac =',i3,f18.9)") lm,funac(lm)
  !      enddo

  do iq = 2,nnn ! omit q=0 point iq=1
     skipg0=.false.
     call auxfunlm(qbz(1:3,iq), alp,alat, qlat, &
          ngc(iq),ngvect(1:3,1:ngc(iq),iq),lx,funa(:,iq),skipg0)
     !        funa(:,iq)= funa(:,iq) - funac(:) !mar2016takao we now consider
  enddo

  do lm=1,(lx+1)**2
     wsumau(lm) = sum(funa(lm,2:nnn))/dble(nnn)
     !c        wsumau(lm) = sum(wqfac(2:nnn)*funa(lm,2:nnn))/dble(nnn)
     !         print *
     !         do iii=2,nnn
     !            print *,' zzz lm iii funa=',lm,iii,qbz(1:3,iii),funa(lm,iii)
     !         enddo
     !         if(abs(wsumau(lm))>1d-6) write(6,"('  wsum fnua=',i3,8f10.5)") lm,wsumau(lm),funac(lm)
  enddo

  !! wklm(lm) = f_L in Eq.(28) T.Kotani,JPSJ83,094711,(2014). mar2016takao:
  !! We rewrite code equivalently, following Eq.(28).
  !! To keep positive definiteness of integral, we only use wklm(1),
  !! that is, f_L (L=lm=(0,0)) only.
  lm=1
  wklm(lm) = wtrue00- wsumau(lm) - funac(1)/dble(nnn)
  !      wklm(2:(lx+1)**2) = 0d0    - wsumau(2:(lx+1)**2) - funac(2:(lx+1)**2)/dble(nnn)
  wklm(2:(lx+1)**2) = 0d0 ! 14march12016
  write(*,"('  lm=1 Klm=wtrue00 - wsumau',3f9.4)") wklm(lm),wtrue00-funac(1),wsumau(lm)
  deallocate(funa)
end subroutine getwklm

!========== version 1 ==============
subroutine diele_invariant(q0x,nq0x,symops,ngrp,  epinv,q0i,nq0i,wq0i)
  !! == invariant dielectric tensor given by symmetrization ==
  !! In addition, we have generate it.
  !! Output::
  !!      epinv(3,nq0i): inequivalent tensor.
  !!      q0i(3,nq0i)  : irreducible q-point for tensor.
  !!      nq0i: number of inequivalent k point.
  !!
  implicit none
  integer :: nq0x,ngrp,nq0i,ik,i,j,ig,ix,ixx,jk,ixxin
  real(8) :: q0x(1:3,nq0x),q0i(1:3,nq0x),symops(3,3,ngrp),sym(3,3),q0xoi(3,nq0x)
  real(8) :: epinv(3,3,nq0x),emat(3,3),qr(3),fac,qnorm,wq0i(nq0x)
  real(8),allocatable:: epinv_(:,:,:)

  integer,parameter:: nxxx=9
  real(8):: zzz(nxxx,nxxx),UU(nxxx,nxxx),VT(nxxx,nxxx)
  real(8):: ss(nxxx),sij
  real(8):: tolq!=1d-6

  !! Generate invariant tensor for each q0x
  write(*,*) ' diele_invariant: nq0x=',nq0x
  allocate(epinv_(3,3,nq0x))
  write(6,"(a,i5)")'  === epinv_: all invariant tensor generaged from q0i(3,1:nq0x) ngrp= ===',ngrp
  epinv_= 0d0
  do ik=1,nq0x
     do ig = 1,ngrp
        sym = symops(:,:,ig)
        qr=matmul(sym,q0x(:,ik))
        do i=1,3
           do j=1,3
              epinv_(i,j,ik) = epinv_(i,j,ik) + qr(i)*qr(j)
           enddo
        enddo
        !          write(*,"(i20,3f19.14)")ig,epinv_(1,:,ik)
        !          write(*,"(i20,3f19.14)")ig,epinv_(2,:,ik)
        !          write(*,"(i20,3f19.14)")ig,epinv_(3,:,ik)
        !          write(*,*)
     enddo
     epinv_(:,:,ik) = epinv_(:,:,ik)/sqrt(sum(epinv_(:,:,ik)**2))
     !        write(*,"(3f19.14)")epinv_(1,:,ik)
     !        write(*,"(3f19.14)")epinv_(2,:,ik)
     !        write(*,"(3f19.14)")epinv_(3,:,ik)
     !        write(*,*)
  enddo

  !! obtain independent epsinv_
  write(6,"(a)")'  === epinv: invariant tensor ==='
  epinv=0d0
  ixx=0
  do ik=1,nq0x
     ixxin=ixx
     call gsorth(9,ixx,epinv,epinv_(:,:,ik))
     if(ixx==ixxin+1) then
        q0i(:,ixx)=q0x(:,ik)
     endif
  enddo
  !      do i=1,ixx
  !        write(6,"('epinv=',i3,9f9.4)")i,epinv(:,:,i)
  !      enddo
  nq0i = ixx
  wq0i= 1d0 !dummy

  !$$$!! check agreement.
  !$$$      ixx=0
  !$$$      do ik = 1,nq0x
  !$$$        do ix = 1,ixx
  !$$$          fac= epinv_(1,1,ik)/epinv(1,1,ix)
  !$$$          if( sum(abs(epinv_(:,:,ik)-fac*epinv(:,:,ix))) < 1d-6*sum(abs(epinv_(:,:,ik))) ) then
  !$$$            goto 990
  !$$$          endif
  !$$$        enddo
  !$$$        ixx=ixx+1
  !$$$        q0i(:,ixx)=q0x(:,ik)
  !$$$        qnorm=sum(q0x(:,ik)* matmul(epinv_(:,:,ik),q0x(:,ik)))/sum(q0x(:,ik)**2)
  !$$$c        print *,' qnorm=',qnorm
  !$$$c        print *,' epinv=',epinv_(:,:,ik)
  !$$$        epinv(:,:,ix)=epinv_(:,:,ik)/qnorm
  !$$$  990   continue
  !$$$        wq0i(ix)=wq0i(ix)+1d0/nq0x
  !$$$      enddo
  !$$$      nq0i = ixx

  !! clean zero
  do ik = 1,nq0i
     do i=1,3
        if(abs(q0i(i,ik))<tolq()) then
           q0i(i,ik)=0d0
        endif
     enddo
     do i=1,3
        do j=1,3
           if(abs(epinv(i,j,ik))<1d-8) then
              epinv(i,j,ik)=0d0
           endif
        enddo
     enddo
  enddo
  !!
  do ik=1,nq0i
     write(*,"(2x,3d22.14)")epinv(1,:,ik)
     write(*,"(2x,3d22.14)")epinv(2,:,ik)
     write(*,"(2x,3d22.14)")epinv(3,:,ik)
     write(*,"('  ---')")
  enddo
  do ik=1,nq0i
     do jk=1,nq0i
        sij=sum(epinv(:,:,ik)*epinv(:,:,jk))
        !          print *,ik,jk,sij
        if(ik==jk) sij=sij-1d0
        if(abs(sij)>1d-7) then
           write(*,"(2i3,3d22.14)")ik,jk,sij
           call rx( 'epsinv_invariant: epsinv are not normalized')
        endif
     enddo
  enddo
end subroutine diele_invariant

subroutine q0irre(qibz,nqibz,q0,wt0,nx06,symops,ngrp, q0i,nq0i,wt,plat,ltrans,ifix,irr,nqbz)
  ! inequivalent points
  !! input:  qibz,q0,
  !! output: q0i,nq0i
  ! xxxxxxx NOTE:input q0i(1:3,1:nq0i) is kept. Output nq0i is larger than input nq0i.
  implicit none
  integer :: ixx,ix,i,ngrp,ig,nq0i,nx06,ifix,irr(nx06),nqirr,nqibz,ib,nqbz,ierr,retval
  real(8) :: q0(1:3,nx06),q0i(1:3,nx06),symops(3,3,ngrp),sym(3,3), &
       qt(3) ,wt0(nx06), wt(nx06),qx(3),qibz(3,nqibz)
  logical:: ltrans
  real(8):: plat(3,3)
  real(8):: platt(3,3)
  real(8):: tolq
  logical:: cmdopt0
  irr=0
  if(ltrans) then
     platt=transpose(plat)
  endif
  wt(1:nx06)=0d0
  ixx = 0 !nq0i !=0 apr2013
  do i = 1,nx06
     qt = q0(:,i)
     do ib = 1,nqibz
        if(sum(abs(qibz(:,ib)-qt))<tolq() .OR. (cmdopt0('--allqbz') .AND. i<=nqbz)) then
           ixx = ixx+1
           q0i(:,ixx) = qt
           irr(i)=1  !this is irreducible
           goto 980
        endif
     enddo
980  continue
  enddo

  !! for genMLWFdipole 2022
  if(cmdopt0('--allqbz')) then
     nq0i=ixx
     return
  endif
  !!
  do i = 1,nx06
     qt = q0(:,i)
     ! equivalence check
     do ix = 1,ixx
        do ig = 1,ngrp
           sym = symops(:,:,ig)
           if(ltrans) then
              call rangedq(matmul(platt,(qt-matmul(sym,q0i(:,ix)))), qx)
              if(sum(abs(qx))<tolq()) then
                 ! rite(6,"(a,3i3,3f10.5,' --> ',3f10.5)") ' q0=sym(ig)*q0i(:,ix) q0 ig q0i= ',i, ig, ix, q0i(:,ix), qt
                 ! rite(6,"(a,2i3,3f10.5)") ' q0(i)=sym(ig)*q0i(ix): ix ig q0= ',ix,ig, qt
                 ! rite(ifix) qt,q0i(:,ix),ix,ig !qtarget, q0i(ix),ig
                 wt(ix) = wt(ix)+wt0(i)
                 goto 990
              endif
           else
              if(sum(abs(q0i(:,ix)-matmul(sym,qt)))<tolq()) then !2000 Nov.
                 wt(ix) = wt(ix)+wt0(i)
                 goto 990
              endif
           endif
        enddo
     enddo
     ixx = ixx+1
     q0i(:,ixx) = qt
     irr(i)=1  !this is irreducible
     wt(ixx) = wt(ixx)+wt0(i)
990  continue
  enddo
  nq0i = ixx
end subroutine q0irre

!----------------------
real(8) function aufcc(q)
  !- auxiliarry function for fcc lattice. from Gygi PRB34 4405
  implicit none
  real(8)    :: cosx,cosy,cosz,q(3), pi=3.1415926535897932D0
  cosx = cos(pi*q(1))
  cosy = cos(pi*q(2))
  cosz = cos(pi*q(3))
  aufcc = 1d0/(3d0 - cosx*cosy -cosy*cosz - cosz*cosx)
END function aufcc

!------------------------------------------------------------
subroutine setq0_2(icase,alat,vol,plat,qlat,alpv,qbz,nstbz,nnn,ngc,ngcmx,ngvect, &
     nq0x,nx0,xn,n1q,n2q,n3q, &
     q0x,wt0,nq00i)
  !- search Q points q0x instead of q=0
  !r  wtt=1d0/(n1q*n2q*n3q) is the weight for each q-points.
  !r  We have to determine Q so that wtrue  = wsumau + wtt*aufcc(Q) .
  ! -------------------------------------------------------------------
  implicit none
  integer :: iq,i,nnn,nx0,icase,iqini,nq0x,q0pc
  real(8)    :: q0x(3,nq0x),qax(3,nq0x),qbx(3,nq0x),qcx(3,nq0x),alpv(3), &
       alat,qbz(3,nnn),auxfun,auxfun6x,qa0, &
       qlat(3,3) ,funa(nnn),qini,plat(3,3), &
       snorm, qmm,vol,volinv,auxf0,auxfun6xnx, &
       wtt,wsumau,wtrue,qa,qb,fa,fb,fc,fder,qc,qx,wt0(nq0x), &
       q2oq1,wtq1,wtq2,www,wgtq0p,qout(3), &
       pi=3.1415926535897932D0, escale,qmin
  integer :: ngc(nnn), ngcmx, ngvect(3,ngcmx,nnn),nstbz(nnn) !,ixtest
  real(8) :: xn    !ratio of weight is xn:1-xn for nx0=2 case
  integer :: q0pchoice,iix,nq,nq00i,verbose,iq1,iq2,iq3,nxx,ni,ne !,auxfunq0p
  integer:: n1q,n2q,n3q
  !      real(8)::alp2
  !---------------------------------
  !      alp = alpv(1)

  !----------------------------
  !  q0chice= 10n1n2 mode 10th intger--->n1
  !                        1th intger--->n2
  integer::nn1,nn2,i1,i2,ix,nsc,nx
  real(8):: dq1(3),dq2(3),ez(3,nq0x),wres,eee,aaa,qmax,vol1,vol2,w1,w2,qm(3,3)

  !      logical:: mode1d=.false.
  !----------------------------
  q0x=0d0
  iix= q0pchoice()
  print *, ' setq0_2: q0pchoice() icase=',iix,icase
  if(q0pchoice()<0 .AND. icase==2) call rx('this branch is not maintained now')
  ! if(q0pchoice()<0 .AND. icase==2) then !see wgtq0p in switch!
  !    nx = abs(q0pchoice())
  !    if( .FALSE. ) then
  !       !------------------------
  !       ! testing case
  !       !        nxx=(nnn)**(1/3d0)+1d-10 !this is a test case for n1=n2=n3
  !       qm(:,1) = qlat(:,1)/dble(2*nx)/n1q
  !       qm(:,2) = qlat(:,2)/dble(2*nx)/n2q
  !       qm(:,3) = qlat(:,3)/dble(2*nx)/n3q
  !       print *,'n1 n2 n3 2*nx qm=',n1q,n2q,n3q,nnn,2*nx,qm
  !       ni= -nx
  !       ne=  nx-1
  !       !        if(mod(nx,2)==0) then
  !       !          ni= -nx/2
  !       !          ne= nx/2-1
  !       !        else
  !       !           stop ' not implimented xxxxxxxxxxxxxxxxxxxxxxx'
  !       !          ni= -(nx+1)/2
  !       !          ne=  (nx+1)/2-1
  !       !        endif
  !       ix=0
  !       do iq1 = ni,ne
  !          do iq2 = ni,ne
  !             do iq3 = ni,ne
  !                ix=ix+1
  !                q0x(:,ix) = qm(:,1)*(iq1+.5) &
  !                     + qm(:,2)*(iq2+.5) &
  !                     + qm(:,3)*(iq3+.5)
  !             enddo
  !          enddo
  !       enddo
  !       nq00i = ix
  !       wt0(1:ix) = 1d0/ix
  !       return
  !       !-----------------------------
  !    endif
  !    ! ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !    !--- scaling box case :
  !    !        goto 1120 ! Not used now
  !    ! wgtq0p in swiches.f should be given correctly --->it should be consistent with wt0 given here
  !    eee= escale()
  !    nsc=0
  !    qmin=1d10
  !    do iq = 1,nnn ! omit q=0 point
  !       if(nstbz(iq)/=0) then
  !          nsc=nsc+1
  !          call shorbz(qbz(1:3,iq),qout,qlat,plat)
  !          if(sum(qout**2)>qmax) qmax=sum(qout**2)
  !       endif
  !    enddo
  !    print *,' nearest iq qmin=',nsc,qmin
  !    if(nsc/=8 ) call rx( ' setq0_2: err.we assumce nsc==8')
  !    call wgtscale(eee,w1,w2)
  !    ix=0
  !    do iq = 1,nnn
  !       call shorbz(qbz(1:3,iq),qout,qlat,plat)
  !       !          if(nstbz(iq)/=0) then
  !       ! only shortest
  !       if(nstbz(iq)/=0 .AND. sum(qout**2)>qmax-1d-4 ) then
  !          ! c            write(6,"(i3,2x,3f16.7,2x,f8.3)") iq,qbz(1:3,iq),sum(qbz(1:3,iq)**2)
  !          write(6,"(i3,2x,3f16.7,2x,f8.3)") iq,qout,sum(qout**2)
  !          do i=1,nx
  !             ix=ix+1
  !             q0x(:,ix) = qout * eee**(nx+1-i)
  !             if(i==1) then
  !                !                wt0(ix) = 1d0
  !                wt0(ix) = (eee**(nx+1-i))**3
  !                !2                vol2    = (eee**(nx+1-i))**3*(1d0/eee**3 - 1d0 )
  !                !2                wt0(ix) = (eee**(nx+1-i))**3 + vol2*w1
  !             else
  !                !                wt0(ix) = 1d0
  !                wt0(ix) = (eee**(nx+1-i))**3 *( 1d0 - eee**3  )
  !                !2                vol2    = (eee**(nx+1-i))**3*(1d0/eee**3 - 1d0 )
  !                !2                vol1    = (eee**(nx+1-i))**3*( 1d0 - eee**3  )
  !                !2                wt0(ix) = vol2*w1 + vol1*w2
  !             endif
  !             !              write(6,"(i3,f16.7,2x,3f16.7,2x,f9.4)") ix,wt0(ix),qout,sum(qout**2)
  !          enddo
  !       endif
  !    enddo
  !    wt0(1:ix)=wt0(1:ix)/sum(wt0(1:ix))
  !    nq00i = ix
  !    return
  !    ! 1120   continue
  !    ! ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! endif

  !      write(6,*) ' newq0p=',iix
  q0pc=q0pchoice()
  if(q0pc==0) then
     ! search q0 ----------------------
     snorm = sqrt( sum ((qlat(:,1)**2)) )
     q0x(:,1)=  qlat(:,1)/snorm
     q0x(:,2)= -qlat(:,1)/snorm
     snorm = sqrt( sum ((qlat(:,2)**2)) )
     q0x(:,3)=  qlat(:,2)/snorm
     q0x(:,4)= -qlat(:,2)/snorm
     snorm = sqrt( sum ((qlat(:,3)**2)) )
     q0x(:,5)=  qlat(:,3)/snorm
     q0x(:,6)= -qlat(:,3)/snorm
     nq=6
  elseif(q0pc==1) then  ! .OR. q0pc==5 .OR. q0pc==6) then
     q0x(:,1)=  (/1d0,0d0,0d0/)
     q0x(:,2)= -(/1d0,0d0,0d0/)
     q0x(:,3)=  (/0d0,1d0,0d0/)
     q0x(:,4)= -(/0d0,1d0,0d0/)
     q0x(:,5)=  (/0d0,0d0,1d0/)
     q0x(:,6)= -(/0d0,0d0,1d0/)
     nq=6
  elseif(q0pc==2) then
     snorm = sqrt( sum ((qlat(:,1)**2)) )
     q0x(:,1)=  qlat(:,1)/snorm
     q0x(:,2)= -qlat(:,1)/snorm
     snorm = sqrt( sum ((qlat(:,2)**2)) )
     q0x(:,3)=  qlat(:,2)/snorm
     q0x(:,4)= -qlat(:,2)/snorm
     nq=4
  elseif(q0pc==3) then
     snorm = sqrt( sum ((qlat(:,1)**2)) )
     q0x(:,1)=  qlat(:,3)/snorm
     q0x(:,2)= -qlat(:,3)/snorm
     nq=2
     !      elseif(q0pchoice()==4) then
     !        snorm = sqrt( sum ((qlat(:,3)**2)) )
     !        q0x(:,1)=   qlat(:,3)/snorm
     !        q0x(:,2)=  -qlat(:,3)/snorm
     !        nq=2
     !------------------------------------------

     !      elseif(q0pchoice()==5) then
     !        q0x(:,1)=  (/1d0,0d0,0d0/)
     !        q0x(:,2)= -(/1d0,0d0,0d0/)
     !        nq=2
     !      elseif(q0pchoice()==6) then
     !        q0x(:,1)=  (/0d0,1d0,0d0/)
     !        q0x(:,2)= -(/0d0,1d0,0d0/)
     !        q0x(:,3)=  (/0d0,0d0,1d0/)
     !        q0x(:,4)= -(/0d0,0d0,1d0/)
     !        nq=4

     !      elseif(q0pchoice()/1000==1) then
     !        mode1d=.true.
     !        if(nx0/=1) stop 'setq0_2: Wrong! q0pchoice()/1000==1 and nx0/=1'
     !        nn1= (q0pchoice()-1000)/10
     !        nn2= mod(q0pchoice()-1000,10)
     !        print *,' mode1d nn1 nn2=',nn1,nn2
     !        dq1= qlat(:,1)/2/nn1
     !        dq2= qlat(:,2)/2/nn2
     !        do i1=1,2*nn1
     !        do i2=1,2*nn2
     !          ix = ix+1
     !          q0x(:,ix) = dq1*(i1 - nn1 -.5d0) + dq2 *(i2 - nn2 -.5d0)
     !          write(6,"(' mode1d ix q0x=',i3,3f9.4)")ix,q0x(:,ix)
     !        enddo
     !        enddo
     !        nq=ix
     !        do ix=1,nq
     !        ez(1:3,ix)= (/0d0,0d0,1d0/)
     !        enddo
  else
     call rx( "q0irre: wrong q0pchoice! ")
  endif

  !-- check wirte near gamma points of the periodic 1/q^2 function, auxfun.
  !      iq=1
  !      do i = 1,99 ! omit q=0 point
  !        qx  =i*0.01
  !        write(1024,*) qx, auxfun((/qx,qx,qx/), alp,alat, qlat,
  !     &                        ngc(iq),ngvect(1:3,1:ngc(iq),iq))
  !      enddo
  !------
  volinv  = (2*pi)**3/vol
  funa(1) = 0d0
  iqini = 2
  if(icase==2) iqini=1
  print *,' iqini nnn=',iqini,nnn
  do iq = iqini,nnn ! omit q=0 point
     www=1d0
     if(icase==2 .AND. nstbz(iq)/=0) www = 1d0-wgtq0p()/nstbz(iq)
     if(verbose() >200) print *,' iq www=',iq,www
     funa(iq) = auxfun(qbz(1:3,iq), alpv,alat, qlat, &
          ngc(iq),ngvect(1:3,1:ngc(iq),iq)) * www
     !        write(6,"('  fnua=',i3,d13.6,3d12.5)") iq,funa(iq),alpv
  enddo
  wsumau = sum(funa)/dble(nnn)

  ! cccccccccccccccccccccccccccccccccccccccc
  ! TEST
  !      if(auxfunQ0P()==1) then
  !        wtrue  = 4*pi/volinv *( 1d0/alp)
  !      alp2 = -5d0
  !      wtrue  = 4*pi/volinv *(  alp2* 1d0/(2d0*alp) + sqrt(pi/alp)/2d0)
  !      else !original version
  wtrue  = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))
  !        if(ixtest()==1) then
  !          wtrue = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))/3d0
  !        elseif(ixtest()==2) then
  !          wtrue = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))*2d0/3d0
  !        endif
  !      endif
  ! cccccccccccccccccccccccccccccccccccccccc

  if(icase==2) then
     wtt = wgtq0p()/dble(nnn)
  else
     wtt = 1d0/dble(nnn)
  endif

  auxf0  = (wtrue - wsumau)/wtt        ! seach q so that auxf0 = aufcc(q).
  write(6,*) " wsum0 wtrue   =", wsumau, wtrue
  write(6,*) " (wtrue - wsumau)/wtt =", auxf0
  if(icase==2 .AND. wtrue - wsumau<0) then
     write(6,*)' Negative Weight CASE!!! wgtQ0p=',wgtq0p()
  endif
  if(auxf0 <0d0) then
     call rx( 'setq0_2: (wtrue - wsumau)/wtt')
  endif

  qini = sqrt(1/auxf0)*alat/(2*pi)/sqrt(minval(alpv))   ! initial guess for q

  if(nx0==2) then
     if(nq/=6) call rx( ' nx0==2 .AND. nq/=6 not implimented yet')
     q0x(:,7:12) = q0x(:,1:6)
     nq00i=12
     wt0(1:6)  = xn/6d0
     wt0(7:12) = (1d0-xn)/6d0
     q2oq1 = sqrt(xn/(xn-1d0))
  elseif(nx0==1) then
     nq00i=nq
     wt0  =1d0/nq
  else
     call rx( ' setq0: I only support nx0=1 and nx0=2')
  endif

  print *,' setq0_2: q0pchoice nq0x=',q0pchoice(),nq0x
  !      if(mode1d) then
  !        print *, ' 1-dimentional mode '
  !      endif

  ! --- find qa for auxf0 = aufcc(qa).
  write(6,*) "  qini = ", qini
  iq=1
  qa = qini
  qb = qini + 0.001d0

  do i= 1,100
     !        if(mode1d) then
     !          qax= q0x + ez*qa
     !          qbx= q0x + ez*qb
     !        else
     qax= qa*q0x
     qbx= qb*q0x
     !        endif
     if(nx0==1) then
        fa = auxfun6x(qax,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),nq)
        fb = auxfun6x(qbx,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),nq)
     else
        fa = auxfun6xnx(qax,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
        fb = auxfun6xnx(qbx,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
     endif
     if(fa-auxf0 < -max(0.1*abs(auxf0),1d0) ) then
        qa = qa *0.5
        qb = qa *(1d0+ 0.001d0)
        print *,' qa --->qa*.5'
        cycle
     endif
     ! cccccccccccccccccccccccccc
     write(6,*)' qa fa-auxf0 =',qa,fa-auxf0
     write(6,*)' qb fb-auxf0 =',qb,fb-auxf0
     ! ccccccccccccccccccccccccccc
     fder = (fb-fa)/(qb-qa)
     qc = qa + (auxf0-fa)/fder
     qb = qa
     qa = qc
     if(abs(qa-qb) <1d-8) exit
  enddo
  write(6,*) " i qmm and fdiff  = ",i,qa, (fa-fb)



  !  refine qa
  qa0 = qa
  qb= qa0+1d-6
  qa= qa0-1d-6
  do i= 1,100
     ! cccccccccccccccccccccccccc
     !        write(6,"(' repeat i qa fa-auxf0 =',i4,d18.10,d13.5)')") i,qa,fa-auxf0
     !        write(6,"('          qb fb-auxf0 =',i4,d18.10,d13.5)')") i,qb,fb-auxf0
     ! ccccccccccccccccccccccccccc
     if(nx0==1) then
        qax=qa*q0x
        fa = auxfun6x(qax,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),nq)
        qbx=qb*q0x
        fb = auxfun6x(qbx,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),nq)
        qc =0.5d0*(qa+qb)
        qcx=qc*q0x
        fc = auxfun6x(qcx,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),nq)
     else
        fa = auxfun6xnx(qa*q0x,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
        fb = auxfun6xnx(qb*q0x,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
        qc =0.5d0*(qa+qb)
        fc = auxfun6xnx(qc*q0x,alpv,alat,qlat,ngc(iq) &
             ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
     endif

     if( (auxf0-fa)*(auxf0-fb) >0d0 ) then
        write(6,*) i,fa,fb,auxf0
        write(6,*) qa,qb
        call rx( ' setq0: something wrong (auxf0-fa)*(auxf0-fb) >0d0')
     elseif((auxf0-fa)*(auxf0-fc) <0d0 ) then
        qb=qc
     elseif((auxf0-fb)*(auxf0-fc) <0d0 ) then
        qa=qc
     endif
     if(abs(fa-fb) <1d-10) exit
  enddo
  ! ----------------------------------------
  ! seach q so that auxf0 = aufcc(q).
  !      write(6,*) " --- qmm norm = ", qa
  !      write(6,*) " sum check  i=",i
  !      write(6,*)
  !     &  auxfun((/qa,qa,qa/), alp,alat,qlat,ngc(iq),ngvect(:,:,iq))
  !     &, alp,alat,qlat,ngc(iq),ngvect(:,:,iq)), auxf0

  write(6,'(" i fdiff qdiff= ",i4,3d24.16)')i,(fa-fb),(qa-qb)

  if(nx0==1) then
     qax= qa*q0x
     write(6,*) " wsum  wtrue   =", &
          wsumau + auxfun6x(qax, alpv,alat,qlat,ngc(iq) &
          ,ngvect(1:3,1:ngc(iq),iq),nq)*wtt, wtrue

     q0x(:,1:nq)=  qa* q0x(:,1:nq)
     write(6,'(" q0x(: 1) = ",3d24.16)') q0x(1:3,1)
  else
     wtq1  = xn/6d0
     wtq2  = (1d0-xn)/6d0
     q2oq1 = sqrt(xn/(xn-1d0))
     qb    = q2oq1*qa
     write(6,*) " wsum  wtrue   =", &
          wsumau &
          + auxfun6xnx(qa*q0x, alpv,alat,qlat,ngc(iq) &
          ,ngvect(1:3,1:ngc(iq),iq),xn,nq)*wtt &
          , wtrue
     q0x(:,1:6)=   qa* q0x(:,1:6)
     write(6,'(" q0x(: 1) = ",3d24.16)') q0x(1:3,1)
     q0x(:,7:12)=  qb* q0x(:,7:12)
     write(6,'(" q0x(: 7) = ",3d24.16)') q0x(1:3,7)
  endif
  !--------------------
  write (6,*) 'sumcheck q0x',sum(abs(q0x(1:3,1:nq)))


end subroutine setq0_2


!-------------------------------------------------------------
real(8) function auxfun(q, alpv, alat, qlat, ngc, ngvect)
  implicit none
  real(8) :: alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2 &
       ,pi=3.1415926535897932D0,alpv(3),alp2,ccc
  integer :: ig, ngc, ngvect(3,ngc),auxfunq0p,ix !,ixtest
  complex(8) :: formfac_test
  tpiba  = 2*pi/alat
  auxfun = 0d0

  !      ix=ixtest()
  !      if(ix==0) then
  do ig = 1,ngc
     qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
     qg2     = sum(alpv(1:3)*qg(1:3)**2)
     auxfun = auxfun + exp(-qg2)/qg2
  enddo
  !      elseif(ix==1) then
  !        do ig = 1,ngc
  !        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
  !        qg2     = sum(alpv(1:3)*qg(1:3)**2)
  !        ccc = qg(3)**2/sum(qg(1:3)**2)
  !        auxfun = auxfun + ccc*exp(-qg2)/qg2
  !        enddo
  !      elseif(ix==2) then
  !        do ig = 1,ngc
  !        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
  !        qg2     = sum(alpv(1:3)*qg(1:3)**2)
  !        ccc = (qg(1)**2+qg(2)**2)/sum(qg(1:3)**2)
  !        auxfun = auxfun + ccc*exp(-qg2)/qg2
  !        enddo
  !      endif
END function auxfun

! cccccccccccccccccccccc
!      alp2=-5d0
! cccccccccccccccccccccc
! ccc TEST cccccccccccccccc
!      if(auxfunq0p()==1) then
!       auxfun = auxfun + exp(-alp*sqrt(qg2))/qg2
! c       auxfun = auxfun  + alp2* exp(-alp*qg2)/sqrt(qg2) + exp(-alp*qg2)/qg2
! c       auxfun = auxfun + abs(formfac_test(qg,alat))**2/qg2
!      else
!       auxfun = auxfun + exp(-qg2)/qg2
!      endif
! cccccccccccccccccccccccccccccccccccccccc
!      enddo
!      end

real(8) function auxfun6x(q0x, alpv, alat, qlat, ngc, ngvect,nq)
  implicit none
  integer :: i,ngc, ngvect(3,ngc),nq
  real(8) :: q0x(3,nq)
  real(8) :: qa,alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2 &
       ,pi=3.1415926535897932D0,alpv(3),auxfun

  auxfun6x = 0d0
  do i=1,nq
     auxfun6x = auxfun6x + &
          auxfun( (/q0x(1,i), q0x(2,i), q0x(3,i)/) &
          , alpv,alat,qlat,ngc,ngvect)
  enddo
  auxfun6x = auxfun6x/nq
END function auxfun6x

real(8) function auxfun6xnx(q0x, alpv, alat, qlat, ngc, ngvect &
     ,xn,nq)
  implicit none
  integer :: i,ngc, ngvect(3,ngc),nq
  real(8) :: q0x(3,nq),xn,wtq1,wtq2,q2oq1,qb
  real(8) :: qa,alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2 &
       ,pi=3.1415926535897932D0,alpv(3),auxfun
  wtq1  = xn/nq
  wtq2 = (1d0-xn)/nq
  q2oq1 = sqrt(xn/(xn-1d0))
  qb= q2oq1 !*qa

  auxfun6xnx = 0d0
  do i=1,nq
     auxfun6xnx = auxfun6xnx &
          + wtq1*auxfun( (/q0x(1,i), q0x(2,i), q0x(3,i)/) &
          , alpv,alat,qlat,ngc,ngvect) &
          + wtq2*auxfun( (/qb*q0x(1,i), qb*q0x(2,i), qb*q0x(3,i)/) &
          , alpv,alat,qlat,ngc,ngvect)
  enddo
END function auxfun6xnx


!-------------------------------------------------------------
complex(8) function formfac_test(qg,alat)
  complex(8):: imag=(0d0,1d0)
  integer,parameter:: natom=2
  integer::ia
  real(8):: absqg,qg(3),bas(3,natom),z(natom),valn(natom),atomform0,alat
  formfac_test=0d0
  absqg  = sqrt(sum(qg(1:3)**2))
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      natom=2
  bas(:,1)=(/0D0,0d0, 0.2408438061041292D+00/)
  bas(:,2)=(/0D0,0d0,-0.2408438061041292D+00/)
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      alp=1d0
  !      formfac_test =  sqrt(exp(-alp*absqg*absqg))
  !      return
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccc
  do ia=1,natom
     formfac_test= formfac_test + atomform0(absqg)*exp(-imag*sum(qg*alat*bas(:,ia)))
  enddo
END function formfac_test

!--------------------------------------------------------------------------------
complex(8) function formfac(qg,bas,natom,z,valn)
  complex(8):: imag=(0d0,1d0)
  integer:: natom,ia
  real(8):: absqg,qg(3),bas(3,natom),z(natom),valn(natom),atomform0
  formfac=0d0
  do ia=1,natom
     absqg  = sqrt(sum(qg(1:3)**2))
     formfac= formfac + atomform0(absqg)*exp(-imag*sum(qg*bas(:,ia)))
  enddo
END function formfac

subroutine getnvaln(konfig,z,natom,nl,iclass,nclass,valn)
  real(8):: valn(natom),z(natom)
  integer:: ia,ic,l,nl,nclass,natom, konfig(0:nl-1,nclass),iclass(natom)
  valn  = 0d0
  do ia = 1,natom
     ic  = iclass(ia)
     valn(ia)  = z(ic)
     print *,' ia z(ic)=',ia, z(ic)
     do    l = 0,nl-1
        print *,' l (konfig(l,ic)-l-1) 2*(2l+1)=',l,(konfig(l,ic)-l-1),( 2*l +1)*2
        valn(ia)  = valn(ia) - (konfig(l,ic)-l-1) *( 2*l +1)*2
     end do
     print *,' valn(ia)=',ia,valn(ia)
  enddo
end subroutine getnvaln

real(8) function atomform0(aqg) !(aqg,z,valn)
  real(8):: kappa2 =(1d0/3d0)**2,aqg !aqg was missing. Add aqg here on 24Apr2012.
  !      atomform0 = 1d0/(kappa2 +aqg**2)**2/aqg**2
  atomform0 = 1d0/(kappa2 +aqg**2)**2/aqg**2
END function atomform0






!------------------------------------------------------------
subroutine setq0x_notused(alat, qbas,  qbz,n1q,n2q,n3q, &
     q0x )
  !- search Q points q0x instead of q=0
  !r  wtt=1d0/(n1q*n2q*n3q) is the weight for each q-points.
  !r  We have to determine Q so that wtrue  = wsumau + wtt*aufcc(Q) .
  !r  This version is only for fcc --------------------
  ! -------------------------------------------------------------------
  implicit none
  integer :: n1q,n2q,n3q,iq,i
  real(8)    :: q0x(3,6),qx, &
       alat,qbz(3,n1q*n2q*n3q),qlat,aufcc,qa0, &
       qbas(3,3), pi=3.1415926535897932D0 ,funa(n1q*n2q*n3q),qini, &
       snorm, qmm,  wtt,wsumau,wtrue,aufccq0,qa,qb,fa,fb,fc,fder,qc
  print *, ' setq0x:'
  ! cccccccccccccc
  iq=1
  do i = 1,99 ! omit q=0 point
     qx  =i*0.01
     write(1025,*) qx, aufcc((/qx,qx,qx/))
  enddo
  ! ccccccccccccccccc

  funa(1)=0d0
  do iq = 2,n1q*n2q*n3q ! omit q=0 point
     funa(iq)= aufcc( qbz(1:3,iq))
  enddo
  wtt     = 1d0/dble(n1q*n2q*n3q)
  wsumau  = sum(funa)*wtt
  wtrue   = 4.423758d0/pi**2
  aufccq0 = (wtrue - wsumau)/wtt ! seach q so that aufccq0 = aufcc(q).
  qini  = sqrt(4/(aufccq0*alat**2))/sqrt(3d0) !initial guess for q

  print *,"  wsum ", wsumau*(alat/2d0)**2, wtrue*(alat/2d0)**2
  write(6,*) " qini = ", qini


  ! --- find qa for aufccq0 = aufcc(qa).
  qa = qini
  qb = qini + 0.001d0
  do i= 1,100
     fa = aufcc((/qa,qa,qa/))
     fb = aufcc((/qb,qb,qb/))
     fder = (fb-fa)/(qb-qa)
     qc = qa + (aufccq0-fa)/fder
     qb = qa
     qa = qc
     if(abs(qa-qb) <1d-8) exit
  enddo
  !  refine qa
  qa0 = qa
  qb= qa0+1d-6
  qa= qa0-1d-6
  do i= 1,100
     fa = aufcc((/qa,qa,qa/))
     fb = aufcc((/qb,qb,qb/))
     qc =0.5d0*(qa+qb)
     fc = aufcc((/qc,qc,qc/))
     if( (aufccq0-fa)*(aufccq0-fb) >0d0 ) then
        write(6,*) i,fa,fb,aufccq0
        write(6,*) qa,qb
        call rx( ' setq0x: something wrong (aufccq0-fa)*(aufccq0-fb) >0d0')
     elseif((aufccq0-fa)*(aufccq0-fc) <0d0 ) then
        qb=qc
     elseif((aufccq0-fb)*(aufccq0-fc) <0d0 ) then
        qa=qc
     endif
     if(abs(fa-fb) <1d-14) exit
  enddo

  ! ----------------------------------------
  write(6,*) " i qmm and err    = ",i,qa, (qa-qb)
  write(6,*) " wtrue and wsum = ", &
       wtrue, wsumau+ aufcc((/qa,qa,qa/))*wtt ! seach q so that aufccq0 = aufcc(q).
  !      write(6,*) " sum check  i=",i
  write(6,*) aufcc((/qa,qa,qa/)), aufccq0

  qa = sqrt(3d0)*qa
  write(6,*) " --- qmm norm = ", qa
  write(6,*)
  !c      qmm = (6/pi/(n1q*n2q*n3q))**(1d0/3)/sqrt(3d0) !!! ????
  !     print *, ' qmm=',qmm

  ! ccccccccccccccccccccccccccccccccccccccccccccc
  ! --- 6-points near q=0 q0(1:3,6)
  !      qbas(1:3,1)/n1q/2
  !      qbas(1:3,2)/n2q/2
  !      qbas(1:3,3)/n3q/2
  !      vol0  = abs(tripl(qbas,qbas(1,2),qbas(1,3)))/(n1q*n2q*n3q)
  !      qsint = ???  ! \int 1/q~2 within the region of the central microcell.
  !      q0smean =  qsint/vol
  !     sumt    = 1/3d0 * ( 1d0/sum((qbas(1:3,1)/2/n1q)**2) +
  !     &                    1d0/sum((qbas(1:3,2)/2/n2q)**2) +
  !     &                    1d0/sum((qbas(1:3,3)/2/n3q)**2) )
  !     qmm = sqrt(q0smean/sumt)

  qmm = qa
  snorm = sqrt( sum ((qbas(:,1)**2)) )
  q0x(:,1)=  qmm* qbas(:,1)/snorm
  q0x(:,2)= -qmm* qbas(:,1)/snorm
  q0x(:,3)=  qmm* qbas(:,2)/snorm
  q0x(:,4)= -qmm* qbas(:,2)/snorm
  q0x(:,5)=  qmm* qbas(:,3)/snorm
  q0x(:,6)= -qmm* qbas(:,3)/snorm
  write(6,*) " --- q0x = ",q0x(:,1)
  call rx( " setq0x: test end")
end subroutine setq0x_notused

!      integer function ixtest()
!      integer:: q0pchoice
!      ixtest=0
!      if(q0pchoice()==5) ixtest=1
!      if(q0pchoice()==6) ixtest=2
!      end

subroutine zccopy(n,nb,zzz)
  integer::n,nb,i
  real(8)::rrr(n)
  complex(8)::zzz(n)
  do i=1,n
     zzz(i)=rrr(i)
  enddo
end subroutine zccopy


subroutine gsorth(ndim,mx,aset,b)
  !! for gram-schmit diagonalization
  implicit none
  integer::i,mx,ndim
  real(8):: aset(ndim,ndim),b(ndim),bout(ndim)
  bout=b
  do i=1,mx
     bout=bout-aset(:,i)*sum(aset(:,i)*bout)
  enddo
  if(sum(bout**2)<1d-10) return
  mx=mx+1
  !      bout=bout/sqrt(sum(bout**2))
  !      aset(:,mx)=bout
  aset(:,mx)=bout/sqrt(sum(bout**2))
  !      print *,' mx bout=',mx,bout
end subroutine gsorth

