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
  real(8):: tolq

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
!----------------------
subroutine qqirre(qibz,nqibz,symops,ngrp,plat,nqbz, qq,nqq, qqi,nqi,irr)
  implicit none 
  intent(in)::    qibz,nqibz,symops,ngrp,plat,nqbz, qq,nqq
  intent(inout)::                                           qqi,nqi,irr
  integer :: ixx,ix,i,ngrp,ig,nqqi,nqq,irr(nqq),nqirr,nqibz,ib,nqbz,nqi
  real(8) :: qq(1:3,nqq),qqi(1:3,nqq),symops(3,3,ngrp),sym(3,3),qt(3),qz(3),qibz(3,nqibz) !
  real(8):: plat(3,3),platt(3,3),tolq
  platt=transpose(plat)
  ixx=nqi
  do i = 1,nqq
     qt = qq(:,i)
     Equivalencecheck: do ix = 1,ixx
        do ig = 1,ngrp
           sym = symops(:,:,ig)
           call rangedq(matmul(platt,(qt-matmul(sym,qqi(:,ix)))), qz)
           if(sum(abs(qz))<tolq()) then
              goto 990
           endif
        enddo
     enddo Equivalencecheck
     irr(i)=1  !this is irreducible i= for 1,nqq
     ixx = ixx+1
     qqi(:,ixx) = qt
990  continue
  enddo
  nqi = ixx ! nqi is the number of irreducible mesh point.  
end subroutine qqirre

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

