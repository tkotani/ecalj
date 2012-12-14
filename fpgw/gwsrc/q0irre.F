      subroutine auxfunlm(q, alp, alat, qlat, ngc, ngvect,lx,  auxfun,skipg0)
!!  Give auxfun(L) = \sum_G  exp(-alp*(q+G)**2)*Y_L(q+G)
      implicit none
      real(8) :: alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2
     & ,pi=3.1415926535897932D0,alp
      integer(4) :: ig, ngc, ngvect(3,ngc),lx
      real(8):: r2s, auxfun((lx+1)**2),alpqg2
      real(8),allocatable:: cy(:),yl(:)
      logical:: skipg0
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
      tpiba  = 2d0*pi/alat
      auxfun = 0d0
      do ig = 1,ngc
        if(skipg0.and.sum(ngvect(1:3,ig)**2)==0) cycle
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
        qg2     = sum(qg(1:3)**2)
        alpqg2= alp* qg2
        call sylm(qg/sqrt(qg2),yl,lx,r2s) !spherical factor Y( q+G )
        auxfun = auxfun + exp(-alpqg2)/qg2 *cy(:)*yl(:) !cy*yl =Y_L(qg/|qg|)
      enddo
      end
!! -----------------------------------------------------------------------------
      subroutine getwklm(alat,vol,plat,qlat,alp,qbz,nnn,ngc,ngcmx,ngvect,lx,n1q,n2q,n3q,
     o    wklm,wqfac)
!! == spherical integration weight Klm with reference auxially functions ==
!! Output
!!   wklm: K_lm is defined around Eq.35 in Copmuter Physics Comm. 176(2007)1-13.
      implicit none
      integer :: iq,i,nnn,lm,lx
      integer :: ngc(nnn), ngcmx, ngvect(3,ngcmx,nnn)
      real(8)    :: alp, alat,qbz(3,nnn),vol,volinv,
     &      qlat(3,3) ,plat(3,3),wtrue00,
     &  pi=3.1415926535897932D0, wklm((lx+1)**2)
      real(8),allocatable:: funa(:,:),wsumau(:),funac(:)
      logical ::skipg0
      integer:: iii,n1q,n2q,n3q,ndiv,iq1,iq2,iq3,lxx
      real(8):: fcenter,fmean,fout(1),qmic(3,3),qx(3),fsum,wg,ftot,ftotc
      real(8):: wqfac(nnn),tripl

      real(8),allocatable:: q0x(:,:),yl(:),cy(:),yll(:)
      real(8):: qg(3),delta5,delta8,tpiba,kmax(6),qcell1(3),qcell2(3),qcell3(3),r2s
      integer:: nq00i,ll
!! ----------------
!!
      nq00i=216
      allocate(q0x(3,nq00i))
      call sphericaltdesign(q0x,nq00i) !now only for nq0i=216
!!
      qcell1=qlat(:,1)/(n1q*2d0)
      qcell2=qlat(:,2)/(n2q*2d0)
      qcell3=qlat(:,3)/(n3q*2d0)
      print *,'qcell1=',qcell1
      print *,'qcell2=',qcell2
      print *,'qcell3=',qcell3
      tpiba = 2d0*pi/alat
      wklm = 0d0
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      call sylmnc(cy,lx)
      do iq = 1,nq00i
        qg(1:3) = q0x(:,iq)
        call sylm(qg/sum(qg**2),yl,lx,r2s) !spherical factor Y( q+G )
c take distance from six planes. We have to take smallest and positive
        call distancecrossplane(qg, qcell1,qcell2,qcell3,kmax(1))
        call distancecrossplane(qg,-qcell1,qcell2,qcell3,kmax(2))
        call distancecrossplane(qg, qcell2,qcell3,qcell1,kmax(3))
        call distancecrossplane(qg,-qcell2,qcell3,qcell1,kmax(4))
        call distancecrossplane(qg, qcell3,qcell1,qcell2,kmax(5))
        call distancecrossplane(qg,-qcell3,qcell1,qcell2,kmax(6))
        print *,'distance=',iq,minval(kmax)
        wklm(:) = wklm(:) + cy(:)*yl(:)* minval(kmax) !minval(kmax)-->kmax(Omega) !Eq.(34)
c        wklm(:) = wklm(:) + cy(:)*yl(:)* cy(:)*yl(:)
      enddo  
c      wklm(:) = 4*pi/nq00i*wklm(:) 
      print *,'voltot vol*nnn=',vol,8*n1q*n2q*n3q*abs(alat**3*tripl(qcell1,qcell2,qcell3))
      wklm= tpiba*wklm* 4d0*pi/nq00i*vol/((2*pi)**3)
      do i=1,(lx+1)**2
        !if(abs(wklm(i))>1d-10) 
        write(6,'(" getwklm: l lm Wklm=",2i3,f9.4)')ll(i),i,wklm(i) !
      enddo


!! discrete sum problem. Central value v.s. Mean value for each cell. 
      ndiv=10
      print *,' getwklm: n1q n2q n3q 2*ndiv',n1q,n2q,n3q,2*ndiv
      skipg0=.false.
      lxx=0
      ftotc=0d0
      ftot=0d0
      qmic(:,1)=qlat(:,1)/(n1q*ndiv*2)
      qmic(:,2)=qlat(:,2)/(n2q*ndiv*2)
      qmic(:,3)=qlat(:,3)/(n3q*ndiv*2)
      do iq = 2,nnn ! omit q=0 point iq=1
         fsum=0d0
         do iq1=-ndiv,ndiv
         do iq2=-ndiv,ndiv
         do iq3=-ndiv,ndiv
c            print *,'iq1 ',iq1,iq2,iq3
            qx=matmul(qmic,(/iq1,iq2,iq3/)) + qbz(1:3,iq)
            call auxfunlm(qx, alp,alat, qlat, 
     &             ngc(iq),ngvect(1:3,1:ngc(iq),iq),lxx,fout,skipg0)
            wg=1d0
            if(mod(iq1+ndiv,2*ndiv)==0) wg=0.5d0
            if(mod(iq2+ndiv,2*ndiv)==0) wg=wg*0.5d0
            if(mod(iq3+ndiv,2*ndiv)==0) wg=wg*0.5d0
            fsum = fsum + wg*fout(1)
            if(iq1==0 .and.iq2==0.and.iq3==0 ) fcenter=fout(1)
         enddo
         enddo
         enddo
         fmean=fsum/(8d0*ndiv**3)
         ftot = ftot+fmean
         ftotc= ftotc+fcenter
         wqfac(iq) =fmean/fcenter
         write(*,"('iq q fcenter fmean ratio=',i4,3f7.3,2x,3f8.4)")
     &   iq,qbz(1:3,iq),fcenter,fmean,fmean/fcenter
      enddo
      write(*,"('ftot ftotc ftot/ftotc: ',2f13.4)") ftot/(nnn-1),ftotc/(nnn-1),ftot/ftotc

ccccccccccccccccccccccccccccc
      return
ccccccccccccccccccccccccccccc


      allocate(funa((lx+1)**2,nnn),wsumau((lx+1)**2),funac((lx+1)**2))
!! ==== true integal of auxially function ====
!! integral \int d^3k/vol * exp(-alp*k**2) *Y_0
      volinv  = (2*pi)**3/vol
c     wtrue00  = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alp**3) /sqrt(4d0*pi) ! 1d0/sqrt(4pi) is Y00.
      wtrue00 = 4d0*pi/volinv *sqrt(pi)/2d0/sqrt(alp)/sqrt(4d0*pi)!1/sqrt(4pi)=Y00. !bugfix 19nov2012. alp

!! ==== discrete sum (except q=0) of auxially functions ====
c constant part
      skipg0=.true.
      iq=1
      call auxfunlm((/0d0,0d0,0d0/), alp,alat, qlat, 
     &             ngc(iq),ngvect(1:3,1:ngc(iq),iq),lx,funac(:),skipg0)
      do lm=1,(lx+1)**2
        if(abs(funac(lm))>1d-6) print *,' xxx lm funac =',lm,funac(lm)
      enddo

      do iq = 2,nnn ! omit q=0 point iq=1
         skipg0=.false.
        call auxfunlm(qbz(1:3,iq), alp,alat, qlat, 
     &             ngc(iq),ngvect(1:3,1:ngc(iq),iq),lx,funa(:,iq),skipg0)
c        write(6,"('  fnua=',i3,8f10.5)") iq,funa(28:30,iq)
        funa(:,iq)= funa(:,iq) - funac(:)
      enddo

      do lm=1,(lx+1)**2
c        wsumau(lm) = sum(funa(lm,2:nnn))/dble(nnn)
        wsumau(lm) = sum(wqfac(2:nnn)*funa(lm,2:nnn))/dble(nnn)
        if(lm==1.or.lm==7.or.lm==43) then
           do iii=2,nnn
           print *,' zzz lm iii funa=',lm,iii,qbz(1:3,iii),funa(lm,iii)
           enddo
        endif   
        if(abs(wsumau(lm))>1d-6) write(6,"('  wsum fnua=',i3,8f10.5)") lm,wsumau(lm),funac(lm)
      enddo
      lm=1
      wklm(lm)          = wtrue00- wsumau(lm)  -funac(1) 
c      wklm(2:(lx+1)**2) = 0d0    
      wklm(2:(lx+1)**2) = 0d0    - wsumau(2:(lx+1)**2)  -funac(2:(lx+1)**2)
      write(*,"('lm=1 Klm=wtrue00 - wsumau',3f9.4)") wklm(lm),wtrue00-funac(1),wsumau(lm)
      deallocate(funa)

      do i=1,(lx+1)**2
        !if(abs(wklm(i))>1d-10) 
        write(6,'(" getwklm: l lm Wklm=",2i3,f9.4)')ll(i),i,wklm(i) !
      enddo
      stop 'vvvvvvvvvvvvvvvvvvvvvvvvv'
      end

C========== version 1 ==============
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
      real(8) :: q0x(1:3,nq0x),q0i(1:3,nq0i),symops(3,3,ngrp),sym(3,3),q0xoi(3,nq0x)
      real(8) :: epinv(3,3,nq0x),emat(3,3),qr(3),fac,qnorm,wq0i(nq0x)
      real(8),allocatable:: epinv_(:,:,:)

      integer,parameter:: nxxx=9
      real(8):: zzz(nxxx,nxxx),UU(nxxx,nxxx),VT(nxxx,nxxx)
      real(8):: ss(nxxx),sij

!! Generate invariant tensor for each q0x
      write(*,*) ' diele_invariant: nq0x=',nq0x
      allocate(epinv_(3,3,nq0x))
      write(*,*)'  === epinv_: all invariant tensor generaged from q0i(3,1:nq0x) ngrp= ===',ngrp
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
c          write(*,"(i20,3f19.14)")ig,epinv_(1,:,ik)
c          write(*,"(i20,3f19.14)")ig,epinv_(2,:,ik)
c          write(*,"(i20,3f19.14)")ig,epinv_(3,:,ik)
c          write(*,*)
        enddo
        epinv_(:,:,ik) = epinv_(:,:,ik)/sqrt(sum(epinv_(:,:,ik)**2))
        write(*,"(3f19.14)")epinv_(1,:,ik)
        write(*,"(3f19.14)")epinv_(2,:,ik)
        write(*,"(3f19.14)")epinv_(3,:,ik)
        write(*,*)
      enddo

!! obtain independent epsinv_
      write(*,*)'  === epinv: invariant tensor ==='
      epinv=0d0
      ixx=0
      do ik=1,nq0x
        ixxin=ixx
        call gsorth(9,ixx,epinv,epinv_(:,:,ik))
        if(ixx==ixxin+1) then
          q0i(:,ixx)=q0x(:,ik)
        endif
      enddo
c      do i=1,ixx
c        write(6,"('epinv=',i3,9f9.4)")i,epinv(:,:,i)
c      enddo
      nq0i = ixx
      wq0i= 1d0 !dummy

c$$$!! check agreement.
c$$$      ixx=0
c$$$      do ik = 1,nq0x
c$$$        do ix = 1,ixx
c$$$          fac= epinv_(1,1,ik)/epinv(1,1,ix)
c$$$          if( sum(abs(epinv_(:,:,ik)-fac*epinv(:,:,ix))) < 1d-6*sum(abs(epinv_(:,:,ik))) ) then
c$$$            goto 990
c$$$          endif
c$$$        enddo
c$$$        ixx=ixx+1
c$$$        q0i(:,ixx)=q0x(:,ik)
c$$$        qnorm=sum(q0x(:,ik)* matmul(epinv_(:,:,ik),q0x(:,ik)))/sum(q0x(:,ik)**2)
c$$$c        print *,' qnorm=',qnorm
c$$$c        print *,' epinv=',epinv_(:,:,ik)
c$$$        epinv(:,:,ix)=epinv_(:,:,ik)/qnorm
c$$$  990   continue
c$$$        wq0i(ix)=wq0i(ix)+1d0/nq0x
c$$$      enddo
c$$$      nq0i = ixx

!! clean zero
      do ik = 1,nq0i
        do i=1,3
          if(abs(q0i(i,ik))<1d-6) then
            q0i(i,ik)=0d0
          endif
        enddo
        do i=1,3
          do j=1,3
            if(abs(epinv(i,j,ik))<1d-6) then
              epinv(i,j,ik)=0d0
            endif
          enddo
        enddo
      enddo
!!
      do ik=1,nq0i
        write(*,"(3d22.14)")epinv(1,:,ik)
        write(*,"(3d22.14)")epinv(2,:,ik)
        write(*,"(3d22.14)")epinv(3,:,ik)
        write(*,*)
      enddo
      do ik=1,nq0i
        do jk=1,nq0i
          sij=sum(epinv(:,:,ik)*epinv(:,:,jk))
          if(ik==jk) sij=sij-1d0
          if(abs(sij)>1d-10) then
            write(*,"(2i3,3d22.14)")ik,jk,sij
            stop 'epsinv_invariant: epsinv are not normalized'
          endif
        enddo
      enddo
      end


      subroutine q0irre(qibz,nqibz,q0,wt0,nx06,symops,ngrp, q0i,nq0i,wt,plat,ltrans,ifix,irr)
c inequivalent points
      implicit none
      integer(4) :: ixx,ix,i,ngrp,ig,nq0i,nx06,ifix,irr(nx06),nqirr,nqibz,ib
      real(8) :: q0(1:3,nx06),q0i(1:3,nx06),symops(3,3,ngrp),sym(3,3),
     &   qt(3) ,wt0(nx06), wt(nx06),qx(3),qibz(3,nqibz)
      logical:: ltrans
      real(8):: plat(3,3)
      real(8):: platt(3,3)
c      print *,'presense trans=',present(ltrans)
ccccccccc
c        print *, ' --- rotation ---'
c        do ig = 1,ngrp
c         print *, ' ig=',ig
c       do i=1,3
c         write(6,"(3f14.6)") symops(i,1:3,ig)
c        enddo
c       enddo
ccccccccc
      irr=0
      if(ltrans) then
        platt=transpose(plat)
      endif

      wt=0d0
      ixx = 0
      do i = 1,nx06
        qt = q0(:,i)
        do ib = 1,nqibz
          if(sum(abs(qibz(:,ib)-qt))<1d-8) then
c            write(ifix) qt,q0i(:,ix),ix,ig !qtarget, q0i(ix),ig
c            wt(ix) = wt(ix)+wt0(i)
            ixx = ixx+1
            q0i(:,ixx) = qt
            irr(i)=1  !this is irreducible
            goto 980
          endif
        enddo
 980    continue
      enddo
!!
      do i = 1,nx06
        qt = q0(:,i)

c equivalence check
        do ix = 1,ixx
          do ig = 1,ngrp
            sym = symops(:,:,ig)
ccccccccccccccccccccc
c       print *, 'q0irre: i ix sum=',
c     &  i,ix, sum(abs(q0i(:,ix)-matmul(sym,qt)))
cccccccccccccccccccccc
c         if(sum(abs(q0i(:,ix)-matmul(sym,qt)))<1d-10) then
            if(ltrans) then
              call rangedq(matmul(platt,(qt-matmul(sym,q0i(:,ix)))), qx)
              if(sum(abs(qx))<1d-8) then
c             write(6,"(a,3i3,3f10.5,' --> ',3f10.5)") ' q0=sym(ig)*q0i(:,ix) q0 ig q0i= ',i, ig, ix, q0i(:,ix), qt
c             write(6,"(a,2i3,3f10.5)") ' q0(i)=sym(ig)*q0i(ix): ix ig q0= ',ix,ig, qt
c              write(ifix) qt,q0i(:,ix),ix,ig !qtarget, q0i(ix),ig
                wt(ix) = wt(ix)+wt0(i)
                goto 990
              endif
            else
              if(sum(abs(q0i(:,ix)-matmul(sym,qt)))<1d-6) then !2000 Nov.
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
  990   continue
      enddo
      nq0i = ixx
      end

c----------------------
      real(8) function aufcc(q)
c- auxiliarry function for fcc lattice. from Gygi PRB34 4405
      implicit none
      real(8)    :: cosx,cosy,cosz,q(3), pi=3.1415926535897932D0
      cosx = cos(pi*q(1))
      cosy = cos(pi*q(2))
      cosz = cos(pi*q(3))
      aufcc = 1d0/(3d0 - cosx*cosy -cosy*cosz - cosz*cosx)
      end

c------------------------------------------------------------
      subroutine setq0_2(icase,alat,vol,plat,qlat,alpv,qbz,nstbz,nnn,ngc,ngcmx,ngvect,
     i    nq0x,nx0,xn,n1q,n2q,n3q,
     o    q0x,wt0,nq00i)
c- search Q points q0x instead of q=0
cr  wtt=1d0/(n1q*n2q*n3q) is the weight for each q-points.
cr  We have to determine Q so that wtrue  = wsumau + wtt*aufcc(Q) .
c -------------------------------------------------------------------
      implicit none
      integer(4) :: iq,i,nnn,nx0,icase,iqini,nq0x,q0pc
      real(8)    :: q0x(3,nq0x),qax(3,nq0x),qbx(3,nq0x),qcx(3,nq0x),alpv(3), !alp,
     &              alat,qbz(3,nnn),auxfun,auxfun6x,qa0,
     &      qlat(3,3) ,funa(nnn),qini,plat(3,3),
     &      snorm, qmm,vol,volinv,auxf0,auxfun6xnx,
     &  wtt,wsumau,wtrue,qa,qb,fa,fb,fc,fder,qc,qx,wt0(nq0x),
     &  q2oq1,wtq1,wtq2,www,wgtq0p,qout(3),
     &  pi=3.1415926535897932D0, escale,qmin
      integer(4) :: ngc(nnn), ngcmx, ngvect(3,ngcmx,nnn),nstbz(nnn) !,ixtest

      real(8) :: xn    !ratio of weight is xn:1-xn for nx0=2 case

      integer(4) :: q0pchoice,iix,nq,nq00i,verbose,iq1,iq2,iq3,nxx,ni,ne !,auxfunq0p
      integer(4):: n1q,n2q,n3q
c      real(8)::alp2
c---------------------------------
c      alp = alpv(1)

c----------------------------
c  q0chice= 10n1n2 mode 10th intger--->n1
c                        1th intger--->n2
      integer(4)::nn1,nn2,i1,i2,ix,nsc,nx
      real(8):: dq1(3),dq2(3),ez(3,nq0x),wres,eee,aaa,qmax,vol1,vol2,w1,w2,qm(3,3)

c      logical:: mode1d=.false.
c----------------------------
      q0x=0d0
      iix= q0pchoice()
      print *, ' setq0_2: q0pchoice() icase=',iix,icase
c
      if(q0pchoice()<0.and.icase==2) then !see wgtq0p in switch!
        nx = abs(q0pchoice())

        if(.false.) then
c------------------------
c testing case
c        nxx=(nnn)**(1/3d0)+1d-10 !this is a test case for n1=n2=n3
          qm(:,1) = qlat(:,1)/dble(2*nx)/n1q
          qm(:,2) = qlat(:,2)/dble(2*nx)/n2q
          qm(:,3) = qlat(:,3)/dble(2*nx)/n3q
          print *,'n1 n2 n3 2*nx qm=',n1q,n2q,n3q,nnn,2*nx,qm
          ni= -nx
          ne=  nx-1
c        if(mod(nx,2)==0) then
c          ni= -nx/2
c          ne= nx/2-1
c        else
c           stop ' not implimented xxxxxxxxxxxxxxxxxxxxxxx'
c          ni= -(nx+1)/2
c          ne=  (nx+1)/2-1
c        endif
          ix=0
          do iq1 = ni,ne
            do iq2 = ni,ne
              do iq3 = ni,ne
                ix=ix+1
                q0x(:,ix) = qm(:,1)*(iq1+.5)
     &             + qm(:,2)*(iq2+.5)
     &             + qm(:,3)*(iq3+.5)
              enddo
            enddo
          enddo
          nq00i = ix
          wt0(1:ix) = 1d0/ix
          return
c-----------------------------
        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--- scaling box case :
c        goto 1120 ! Not used now
c wgtq0p in swiches.f should be given correctly --->it should be consistent with wt0 given here
        eee= escale()
        nsc=0
        qmin=1d10
        do iq = 1,nnn ! omit q=0 point
          if(nstbz(iq)/=0) then
            nsc=nsc+1
            call shorbz(qbz(1:3,iq),qout,qlat,plat)
            if(sum(qout**2)>qmax) qmax=sum(qout**2)
          endif
        enddo
        print *,' nearest iq qmin=',nsc,qmin
        if(nsc/=8 ) stop ' setq0_2: err.we assumce nsc==8'
        call wgtscale(eee,w1,w2)
        ix=0
        do iq = 1,nnn
          call shorbz(qbz(1:3,iq),qout,qlat,plat)
c          if(nstbz(iq)/=0) then
c only shortest
          if(nstbz(iq)/=0.and. sum(qout**2)>qmax-1d-4 ) then
ccc            write(6,"(i3,2x,3f16.7,2x,f8.3)") iq,qbz(1:3,iq),sum(qbz(1:3,iq)**2)
            write(6,"(i3,2x,3f16.7,2x,f8.3)") iq,qout,sum(qout**2)
            do i=1,nx
              ix=ix+1
              q0x(:,ix) = qout * eee**(nx+1-i)
              if(i==1) then
c                wt0(ix) = 1d0
                wt0(ix) = (eee**(nx+1-i))**3
c2                vol2    = (eee**(nx+1-i))**3*(1d0/eee**3 - 1d0 )
c2                wt0(ix) = (eee**(nx+1-i))**3 + vol2*w1
              else
c                wt0(ix) = 1d0
                wt0(ix) = (eee**(nx+1-i))**3 *( 1d0 - eee**3  )
c2                vol2    = (eee**(nx+1-i))**3*(1d0/eee**3 - 1d0 )
c2                vol1    = (eee**(nx+1-i))**3*( 1d0 - eee**3  )
c2                wt0(ix) = vol2*w1 + vol1*w2
              endif
c              write(6,"(i3,f16.7,2x,3f16.7,2x,f9.4)") ix,wt0(ix),qout,sum(qout**2)
            enddo
          endif
        enddo
        wt0(1:ix)=wt0(1:ix)/sum(wt0(1:ix))
        nq00i = ix
        return
c 1120   continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif

c      write(6,*) ' newq0p=',iix
      q0pc=q0pchoice()
      if(q0pc==0) then
c search q0 ----------------------
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
      elseif(q0pc==1) then  !.or.q0pc==5.or.q0pc==6) then
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
c      elseif(q0pchoice()==4) then
c        snorm = sqrt( sum ((qlat(:,3)**2)) )
c        q0x(:,1)=   qlat(:,3)/snorm
c        q0x(:,2)=  -qlat(:,3)/snorm
c        nq=2
c------------------------------------------

c      elseif(q0pchoice()==5) then
c        q0x(:,1)=  (/1d0,0d0,0d0/)
c        q0x(:,2)= -(/1d0,0d0,0d0/)
c        nq=2
c      elseif(q0pchoice()==6) then
c        q0x(:,1)=  (/0d0,1d0,0d0/)
c        q0x(:,2)= -(/0d0,1d0,0d0/)
c        q0x(:,3)=  (/0d0,0d0,1d0/)
c        q0x(:,4)= -(/0d0,0d0,1d0/)
c        nq=4

c      elseif(q0pchoice()/1000==1) then
c        mode1d=.true.
c        if(nx0/=1) stop 'setq0_2: Wrong! q0pchoice()/1000==1 and nx0/=1'
c        nn1= (q0pchoice()-1000)/10
c        nn2= mod(q0pchoice()-1000,10)
c        print *,' mode1d nn1 nn2=',nn1,nn2
c        dq1= qlat(:,1)/2/nn1
c        dq2= qlat(:,2)/2/nn2
c        do i1=1,2*nn1
c        do i2=1,2*nn2
c          ix = ix+1
c          q0x(:,ix) = dq1*(i1 - nn1 -.5d0) + dq2 *(i2 - nn2 -.5d0)
c          write(6,"(' mode1d ix q0x=',i3,3f9.4)")ix,q0x(:,ix)
c        enddo
c        enddo
c        nq=ix
c        do ix=1,nq
c        ez(1:3,ix)= (/0d0,0d0,1d0/)
c        enddo
      else
        stop "q0irre: wrong q0pchoice! "
      endif

c-- check wirte near gamma points of the periodic 1/q^2 function, auxfun.
c      iq=1
c      do i = 1,99 ! omit q=0 point
c        qx  =i*0.01
c        write(1024,*) qx, auxfun((/qx,qx,qx/), alp,alat, qlat,
c     &                        ngc(iq),ngvect(1:3,1:ngc(iq),iq))
c      enddo
c------
      volinv  = (2*pi)**3/vol
      funa(1) = 0d0
      iqini = 2
      if(icase==2) iqini=1
      print *,' iqini nnn=',iqini,nnn
      do iq = iqini,nnn ! omit q=0 point
        www=1d0
        if(icase==2.and.nstbz(iq)/=0) www = 1d0-wgtq0p()/nstbz(iq)
        if(verbose() >200) print *,' iq www=',iq,www
        funa(iq) = auxfun(qbz(1:3,iq), alpv,alat, qlat, 
     &                   ngc(iq),ngvect(1:3,1:ngc(iq),iq)) * www
c        write(6,"('  fnua=',i3,d13.6,3d12.5)") iq,funa(iq),alpv
      enddo
      wsumau = sum(funa)/dble(nnn)

cccccccccccccccccccccccccccccccccccccccccc
c TEST
c      if(auxfunQ0P()==1) then
c        wtrue  = 4*pi/volinv *( 1d0/alp)
c      alp2 = -5d0
c      wtrue  = 4*pi/volinv *(  alp2* 1d0/(2d0*alp) + sqrt(pi/alp)/2d0)
c      else !original version
      wtrue  = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))
c        if(ixtest()==1) then
c          wtrue = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))/3d0
c        elseif(ixtest()==2) then
c          wtrue = 4*pi/volinv *sqrt(pi)/2d0/sqrt(alpv(1)*alpv(2)*alpv(3))*2d0/3d0
c        endif
c      endif
cccccccccccccccccccccccccccccccccccccccccc

      if(icase==2) then
        wtt = wgtq0p()/dble(nnn)
      else
        wtt = 1d0/dble(nnn)
      endif

      auxf0  = (wtrue - wsumau)/wtt        ! seach q so that auxf0 = aufcc(q).
      write(6,*) " wsum0 wtrue   =", wsumau, wtrue
      write(6,*) " (wtrue - wsumau)/wtt =", auxf0
      if(icase==2.and.wtrue - wsumau<0) then
        write(6,*)' Negative Weight CASE!!! wgtQ0p=',wgtq0p()
      endif
      if(auxf0 <0d0) then
        stop  'setq0_2: (wtrue - wsumau)/wtt'
      endif

      qini = sqrt(1/auxf0)*alat/(2*pi)/sqrt(minval(alpv))   ! initial guess for q

      if(nx0==2) then
        if(nq/=6) stop ' nx0==2.and.nq/=6 not implimented yet'
        q0x(:,7:12) = q0x(:,1:6)
        nq00i=12
        wt0(1:6)  = xn/6d0
        wt0(7:12) = (1d0-xn)/6d0
        q2oq1 = sqrt(xn/(xn-1d0))
      elseif(nx0==1) then
        nq00i=nq
        wt0  =1d0/nq
      else
        stop ' setq0: I only support nx0=1 and nx0=2'
      endif

      print *,' setq0_2: q0pchoice nq0x=',q0pchoice(),nq0x
c      if(mode1d) then
c        print *, ' 1-dimentional mode '
c      endif

c --- find qa for auxf0 = aufcc(qa).
      write(6,*) "  qini = ", qini
      iq=1
      qa = qini
      qb = qini + 0.001d0

      do i= 1,100
c        if(mode1d) then
c          qax= q0x + ez*qa
c          qbx= q0x + ez*qb
c        else
        qax= qa*q0x
        qbx= qb*q0x
c        endif
        if(nx0==1) then
          fa = auxfun6x(qax,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),nq)
          fb = auxfun6x(qbx,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),nq)
        else
          fa = auxfun6xnx(qax,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
          fb = auxfun6xnx(qbx,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
        endif
        if(fa-auxf0 < -max(0.1*abs(auxf0),1d0) ) then
          qa = qa *0.5
          qb = qa *(1d0+ 0.001d0)
          print *,' qa --->qa*.5'
          cycle
        endif
cccccccccccccccccccccccccccc
        write(6,*)' qa fa-auxf0 =',qa,fa-auxf0
        write(6,*)' qb fb-auxf0 =',qb,fb-auxf0
ccccccccccccccccccccccccccccc
        fder = (fb-fa)/(qb-qa)
        qc = qa + (auxf0-fa)/fder
        qb = qa
        qa = qc
        if(abs(qa-qb) <1d-8) exit
      enddo
      write(6,*) " i qmm and fdiff  = ",i,qa, (fa-fb)



c  refine qa
      qa0 = qa
      qb= qa0+1d-6
      qa= qa0-1d-6
      do i= 1,100
cccccccccccccccccccccccccccc
c        write(6,"(' repeat i qa fa-auxf0 =',i4,d18.10,d13.5)')") i,qa,fa-auxf0
c        write(6,"('          qb fb-auxf0 =',i4,d18.10,d13.5)')") i,qb,fb-auxf0
ccccccccccccccccccccccccccccc
        if(nx0==1) then
          qax=qa*q0x
          fa = auxfun6x(qax,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),nq)
          qbx=qb*q0x
          fb = auxfun6x(qbx,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),nq)
          qc =0.5d0*(qa+qb)
          qcx=qc*q0x
          fc = auxfun6x(qcx,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),nq)
        else
          fa = auxfun6xnx(qa*q0x,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
          fb = auxfun6xnx(qb*q0x,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
          qc =0.5d0*(qa+qb)
          fc = auxfun6xnx(qc*q0x,alpv,alat,qlat,ngc(iq)
     &              ,ngvect(1:3,1:ngc(iq),iq),xn,nq)
        endif
c
        if( (auxf0-fa)*(auxf0-fb) >0d0 ) then
          write(6,*) i,fa,fb,auxf0
          write(6,*) qa,qb
          stop ' setq0: something wrong (auxf0-fa)*(auxf0-fb) >0d0'
        elseif((auxf0-fa)*(auxf0-fc) <0d0 ) then
          qb=qc
        elseif((auxf0-fb)*(auxf0-fc) <0d0 ) then
          qa=qc
        endif
        if(abs(fa-fb) <1d-10) exit
      enddo
c ----------------------------------------
! seach q so that auxf0 = aufcc(q).
c      write(6,*) " --- qmm norm = ", qa
c      write(6,*) " sum check  i=",i
c      write(6,*)
c     &  auxfun((/qa,qa,qa/), alp,alat,qlat,ngc(iq),ngvect(:,:,iq))
c     &, alp,alat,qlat,ngc(iq),ngvect(:,:,iq)), auxf0
c
      write(6,'(" i fdiff qdiff= ",i4,3d24.16)')i,(fa-fb),(qa-qb)
c
      if(nx0==1) then
        qax= qa*q0x
        write(6,*) " wsum  wtrue   =",
     &  wsumau + auxfun6x(qax, alpv,alat,qlat,ngc(iq)
     &               ,ngvect(1:3,1:ngc(iq),iq),nq)*wtt, wtrue
c
        q0x(:,1:nq)=  qa* q0x(:,1:nq)
        write(6,'(" q0x(: 1) = ",3d24.16)') q0x(1:3,1)
      else
        wtq1  = xn/6d0
        wtq2  = (1d0-xn)/6d0
        q2oq1 = sqrt(xn/(xn-1d0))
        qb    = q2oq1*qa
        write(6,*) " wsum  wtrue   =",
     &  wsumau 
     &      + auxfun6xnx(qa*q0x, alpv,alat,qlat,ngc(iq)
     &               ,ngvect(1:3,1:ngc(iq),iq),xn,nq)*wtt
     & , wtrue
        q0x(:,1:6)=   qa* q0x(:,1:6)
        write(6,'(" q0x(: 1) = ",3d24.16)') q0x(1:3,1)
        q0x(:,7:12)=  qb* q0x(:,7:12)
        write(6,'(" q0x(: 7) = ",3d24.16)') q0x(1:3,7)
      endif
c--------------------
      write (6,*) 'sumcheck q0x',sum(abs(q0x(1:3,1:nq)))


      end


c-------------------------------------------------------------
      real(8) function auxfun(q, alpv, alat, qlat, ngc, ngvect)
      implicit none
      real(8) :: alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2
     & ,pi=3.1415926535897932D0,alpv(3),alp2,ccc
      integer(4) :: ig, ngc, ngvect(3,ngc),auxfunq0p,ix !,ixtest
      complex(8) :: formfac_test
      tpiba  = 2*pi/alat
      auxfun = 0d0

c      ix=ixtest()
c      if(ix==0) then
      do ig = 1,ngc
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
        qg2     = sum(alpv(1:3)*qg(1:3)**2)
        auxfun = auxfun + exp(-qg2)/qg2
      enddo
c      elseif(ix==1) then
c        do ig = 1,ngc
c        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
c        qg2     = sum(alpv(1:3)*qg(1:3)**2)
c        ccc = qg(3)**2/sum(qg(1:3)**2)
c        auxfun = auxfun + ccc*exp(-qg2)/qg2
c        enddo
c      elseif(ix==2) then
c        do ig = 1,ngc
c        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
c        qg2     = sum(alpv(1:3)*qg(1:3)**2)
c        ccc = (qg(1)**2+qg(2)**2)/sum(qg(1:3)**2)
c        auxfun = auxfun + ccc*exp(-qg2)/qg2
c        enddo
c      endif
      end

cccccccccccccccccccccccc
c      alp2=-5d0
cccccccccccccccccccccccc
ccccc TEST cccccccccccccccc
c      if(auxfunq0p()==1) then
c       auxfun = auxfun + exp(-alp*sqrt(qg2))/qg2
ccc       auxfun = auxfun  + alp2* exp(-alp*qg2)/sqrt(qg2) + exp(-alp*qg2)/qg2
ccc       auxfun = auxfun + abs(formfac_test(qg,alat))**2/qg2
c      else
c       auxfun = auxfun + exp(-qg2)/qg2
c      endif
cccccccccccccccccccccccccccccccccccccccccc
c      enddo
c      end
c
      real(8) function auxfun6x(q0x, alpv, alat, qlat, ngc, ngvect,nq)
      implicit none
      integer(4) :: i,ngc, ngvect(3,ngc),nq
      real(8) :: q0x(3,nq)
      real(8) :: qa,alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2
     & ,pi=3.1415926535897932D0,alpv(3),auxfun
c
      auxfun6x = 0d0
      do i=1,nq
        auxfun6x = auxfun6x +
     &    auxfun( (/q0x(1,i), q0x(2,i), q0x(3,i)/)
     &    , alpv,alat,qlat,ngc,ngvect)
      enddo
      auxfun6x = auxfun6x/nq
      end
c
      real(8) function auxfun6xnx(q0x, alpv, alat, qlat, ngc, ngvect
     &                   ,xn,nq)
      implicit none
      integer(4) :: i,ngc, ngvect(3,ngc),nq
      real(8) :: q0x(3,nq),xn,wtq1,wtq2,q2oq1,qb
      real(8) :: qa,alat,q(3),tpiba,aucfun,qg(3),qlat(3,3),qg2
     & ,pi=3.1415926535897932D0,alpv(3),auxfun
      wtq1  = xn/nq
      wtq2 = (1d0-xn)/nq
      q2oq1 = sqrt(xn/(xn-1d0))
      qb= q2oq1 !*qa
c
      auxfun6xnx = 0d0
      do i=1,nq
        auxfun6xnx = auxfun6xnx 
     &    + wtq1*auxfun( (/q0x(1,i), q0x(2,i), q0x(3,i)/)
     &    , alpv,alat,qlat,ngc,ngvect)
     &    + wtq2*auxfun( (/qb*q0x(1,i), qb*q0x(2,i), qb*q0x(3,i)/)
     &    , alpv,alat,qlat,ngc,ngvect)
      enddo
      end


c-------------------------------------------------------------
      complex(8) function formfac_test(qg,alat)
      complex(8):: imag=(0d0,1d0)
      integer(4),parameter:: natom=2
      integer(4)::ia
      real(8):: absqg,qg(3),bas(3,natom),z(natom),valn(natom),atomform0,alat
      formfac_test=0d0
      absqg  = sqrt(sum(qg(1:3)**2))
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      natom=2
      bas(:,1)=(/0D0,0d0, 0.2408438061041292D+00/)
      bas(:,2)=(/0D0,0d0,-0.2408438061041292D+00/)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      alp=1d0
c      formfac_test =  sqrt(exp(-alp*absqg*absqg))
c      return
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ia=1,natom
        formfac_test= formfac_test + atomform0(absqg)*exp(-imag*sum(qg*alat*bas(:,ia)))
      enddo
      end

c--------------------------------------------------------------------------------
      complex(8) function formfac(qg,bas,natom,z,valn)
      complex(8):: imag=(0d0,1d0)
      integer(4):: natom,ia
      real(8):: absqg,qg(3),bas(3,natom),z(natom),valn(natom),atomform0
      formfac=0d0
      do ia=1,natom
        absqg  = sqrt(sum(qg(1:3)**2))
        formfac= formfac + atomform0(absqg)*exp(-imag*sum(qg*bas(:,ia)))
      enddo
      end
c
      subroutine getnvaln(konfig,z,natom,nl,iclass,nclass,valn)
      real(8):: valn(natom)
      integer(4):: ia,ic,l,nl,nclass,natom, konfig(0:nl-1,nclass),iclass(natom)
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
      end
c
      real(8) function atomform0(aqg) !(aqg,z,valn)
      real(8):: kappa2 =(1d0/3d0)**2,aqg !aqg was missing. Add aqg here on 24Apr2012.
c      atomform0 = 1d0/(kappa2 +aqg**2)**2/aqg**2
      atomform0 = 1d0/(kappa2 +aqg**2)**2/aqg**2
      end






c------------------------------------------------------------
      subroutine setq0x_notused(alat, qbas,  qbz,n1q,n2q,n3q, 
     o    q0x )
c- search Q points q0x instead of q=0
cr  wtt=1d0/(n1q*n2q*n3q) is the weight for each q-points.
cr  We have to determine Q so that wtrue  = wsumau + wtt*aufcc(Q) .
cr  This version is only for fcc --------------------
c -------------------------------------------------------------------
      implicit none
      integer(4) :: n1q,n2q,n3q,iq,i
      real(8)    :: q0x(3,6),qx,
     &              alat,qbz(3,n1q*n2q*n3q),qlat,aufcc,qa0,
     &      qbas(3,3), pi=3.1415926535897932D0 ,funa(n1q*n2q*n3q),qini,
     &      snorm, qmm,  wtt,wsumau,wtrue,aufccq0,qa,qb,fa,fb,fc,fder,qc
      print *, ' setq0x:'
cccccccccccccccc
      iq=1
      do i = 1,99 ! omit q=0 point
        qx  =i*0.01
        write(1025,*) qx, aufcc((/qx,qx,qx/))
      enddo
ccccccccccccccccccc

      funa(1)=0d0
      do iq = 2,n1q*n2q*n3q ! omit q=0 point
        funa(iq)= aufcc( qbz(1:3,iq))
      enddo
      wtt     = 1d0/dble(n1q*n2q*n3q)
      wsumau  = sum(funa)*wtt
      wtrue   = 4.423758d0/pi**2
      aufccq0 = (wtrue - wsumau)/wtt ! seach q so that aufccq0 = aufcc(q).
      qini  = sqrt(4/(aufccq0*alat**2))/sqrt(3d0) !initial guess for q
c
      print *,"  wsum ", wsumau*(alat/2d0)**2, wtrue*(alat/2d0)**2
      write(6,*) " qini = ", qini


c --- find qa for aufccq0 = aufcc(qa).
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
c  refine qa
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
          stop ' setq0x: something wrong (aufccq0-fa)*(aufccq0-fb) >0d0'
        elseif((aufccq0-fa)*(aufccq0-fc) <0d0 ) then
          qb=qc
        elseif((aufccq0-fb)*(aufccq0-fc) <0d0 ) then
          qa=qc
        endif
        if(abs(fa-fb) <1d-14) exit
      enddo

c ----------------------------------------
      write(6,*) " i qmm and err    = ",i,qa, (qa-qb)
      write(6,*) " wtrue and wsum = ",
     &  wtrue, wsumau+ aufcc((/qa,qa,qa/))*wtt ! seach q so that aufccq0 = aufcc(q).
c      write(6,*) " sum check  i=",i
      write(6,*) aufcc((/qa,qa,qa/)), aufccq0
c
      qa = sqrt(3d0)*qa
      write(6,*) " --- qmm norm = ", qa
      write(6,*)
cc      qmm = (6/pi/(n1q*n2q*n3q))**(1d0/3)/sqrt(3d0) !!! ????
c     print *, ' qmm=',qmm

ccccccccccccccccccccccccccccccccccccccccccccccc
c --- 6-points near q=0 q0(1:3,6)
c      qbas(1:3,1)/n1q/2
c      qbas(1:3,2)/n2q/2
c      qbas(1:3,3)/n3q/2
c      vol0  = abs(tripl(qbas,qbas(1,2),qbas(1,3)))/(n1q*n2q*n3q)
c      qsint = ???  ! \int 1/q~2 within the region of the central microcell.
c      q0smean =  qsint/vol
c     sumt    = 1/3d0 * ( 1d0/sum((qbas(1:3,1)/2/n1q)**2) +
c     &                    1d0/sum((qbas(1:3,2)/2/n2q)**2) +
c     &                    1d0/sum((qbas(1:3,3)/2/n3q)**2) )
c     qmm = sqrt(q0smean/sumt)

      qmm = qa
      snorm = sqrt( sum ((qbas(:,1)**2)) )
      q0x(:,1)=  qmm* qbas(:,1)/snorm
      q0x(:,2)= -qmm* qbas(:,1)/snorm
      q0x(:,3)=  qmm* qbas(:,2)/snorm
      q0x(:,4)= -qmm* qbas(:,2)/snorm
      q0x(:,5)=  qmm* qbas(:,3)/snorm
      q0x(:,6)= -qmm* qbas(:,3)/snorm
      write(6,*) " --- q0x = ",q0x(:,1)
      stop  " setq0x: test end"
      end

c      integer(4) function ixtest()
c      integer(4):: q0pchoice
c      ixtest=0
c      if(q0pchoice()==5) ixtest=1
c      if(q0pchoice()==6) ixtest=2
c      end

      subroutine zccopy(n,nb,zzz)
      real(8)::rrr(n)
      complex(8)::zzz(n)
      do i=1,n
        zzz(i)=rrr(i)
      enddo
      end

c
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
c      bout=bout/sqrt(sum(bout**2))
c      aset(:,mx)=bout
      aset(:,mx)=bout/sqrt(sum(bout**2))
c      print *,' mx bout=',mx,bout
      end

      subroutine distancecrossplane(qinr,qmic1,qmic2,qmic3,kmax)
      implicit none
      real(8):: qin(3),qmic1(3),qmic2(3),qmic3(3),ipx,qinr(3),kmax
      real(8):: qlat(3,3),plat(3,3),pp,qq,x,xx
      qin=qinr/(sum(qinr**2))**.5
      qlat(:,1)=qmic1
      qlat(:,2)=qmic2
      qlat(:,3)=qmic3
      call dinv33x (qlat,plat) 
      xx = sum(qin*plat(:,1))
      if(xx<1d-3) then
         kmax=1d10 !
         return
      endif
c      pp = xx*sum(qin*plat(:,2))
c      qq = xx*sum(qin*plat(:,3))
      kmax=sum((qin/xx)**2)**.5
      end

      subroutine sphericaltdesign(q0x,nq00i)
!! spherical t-design
!! http://www2.research.att.com/~njas/sphdesigns/dim3/des.3.216.20.txt
      implicit none
      integer ::i
      integer::nq00i
      real(8):: q0x(3,nq00i)
      character*30000 :: designtxt
      if(nq00i/=216) stop "sphericaltdesin for 3.216.20 set nq00i/=216"
      designtxt="
     & -.2515812993559380E+00
     & .9657024628131560E+00
     & -.6423085809004400E-01
     & .9657024628129730E+00
     & -.6423085809016300E-01
     & -.2515812993566090E+00
     & .2515812993562100E+00
     & -.9657024628130760E+00
     & -.6423085809018400E-01
     & -.6423085809003699E-01
     & -.2515812993564690E+00
     & .9657024628130180E+00
     & -.9657024628129880E+00
     & -.6423085809021201E-01
     & .2515812993565400E+00
     & -.6423085809028301E-01
     & .2515812993562130E+00
     & -.9657024628130680E+00
     & -.9657024628131290E+00
     & .6423085809034999E-01
     & -.2515812993559620E+00
     & .6423085809020899E-01
     & .2515812993563220E+00
     & .9657024628130449E+00
     & .9657024628130900E+00
     & .6423085808991100E-01
     & .2515812993562260E+00
     & .6423085809020000E-01
     & -.2515812993563000E+00
     & -.9657024628130511E+00
     & -.2515812993566000E+00
     & -.9657024628129920E+00
     & .6423085808991900E-01
     & .2515812993565160E+00
     & .9657024628129810E+00
     & .6423085809040199E-01
     & -.7742655338457720E+00
     & .3815151823433970E+00
     & -.5049346975005830E+00
     & .3815151823431970E+00
     & -.5049346975006570E+00
     & -.7742655338458230E+00
     & .7742655338455831E+00
     & -.3815151823433860E+00
     & -.5049346975008830E+00
     & -.5049346975007970E+00
     & -.7742655338456810E+00
     & .3815151823433000E+00
     & -.3815151823431530E+00
     & -.5049346975008050E+00
     & .7742655338457480E+00
     & -.5049346975006220E+00
     & .7742655338458870E+00
     & -.3815151823431140E+00
     & -.3815151823432720E+00
     & .5049346975008830E+00
     & -.7742655338456390E+00
     & .5049346975008080E+00
     & .7742655338456150E+00
     & .3815151823434190E+00
     & .3815151823434900E+00
     & .5049346975006210E+00
     & .7742655338457030E+00
     & .5049346975005800E+00
     & -.7742655338458060E+00
     & -.3815151823433330E+00
     & -.7742655338457189E+00
     & -.3815151823433210E+00
     & .5049346975007230E+00
     & .7742655338458940E+00
     & .3815151823429800E+00
     & .5049346975007110E+00
     & .6218920898658570E+00
     & .4517167996942610E+00
     & -.6396891131137470E+00
     & .4517167996941910E+00
     & -.6396891131139180E+00
     & .6218920898657310E+00
     & -.6218920898656480E+00
     & -.4517167996942250E+00
     & -.6396891131139760E+00
     & -.6396891131139010E+00
     & .6218920898654990E+00
     & .4517167996945350E+00
     & -.4517167996940080E+00
     & -.6396891131137999E+00
     & -.6218920898659860E+00
     & -.6396891131138790E+00
     & -.6218920898656550E+00
     & -.4517167996943510E+00
     & -.4517167996943470E+00
     & .6396891131136750E+00
     & .6218920898658690E+00
     & .6396891131137880E+00
     & -.6218920898659950E+00
     & .4517167996940130E+00
     & .4517167996945870E+00
     & .6396891131139550E+00
     & -.6218920898654060E+00
     & .6396891131140610E+00
     & .6218920898659001E+00
     & -.4517167996937570E+00
     & .6218920898658890E+00
     & -.4517167996942810E+00
     & .6396891131137010E+00
     & -.6218920898658979E+00
     & .4517167996937130E+00
     & .6396891131140940E+00
     & .2818110426750910E+00
     & .8580478476961970E+00
     & -.4293441827838140E+00
     & .8580478476964080E+00
     & -.4293441827836590E+00
     & .2818110426746880E+00
     & -.2818110426751140E+00
     & -.8580478476963060E+00
     & -.4293441827835810E+00
     & -.4293441827833150E+00
     & .2818110426749470E+00
     & .8580478476964950E+00
     & -.8580478476963860E+00
     & -.4293441827833290E+00
     & -.2818110426752570E+00
     & -.4293441827839790E+00
     & -.2818110426747930E+00
     & -.8580478476962130E+00
     & -.8580478476961360E+00
     & .4293441827839480E+00
     & .2818110426750750E+00
     & .4293441827835740E+00
     & -.2818110426750020E+00
     & .8580478476963470E+00
     & .8580478476964301E+00
     & .4293441827834320E+00
     & -.2818110426749640E+00
     & .4293441827834070E+00
     & .2818110426754000E+00
     & -.8580478476963000E+00
     & .2818110426747800E+00
     & -.8580478476965150E+00
     & .4293441827833830E+00
     & -.2818110426751930E+00
     & .8580478476962270E+00
     & .4293441827836880E+00
     & -.6496120041073690E+00
     & -.6153110840694710E+00
     & .4465383678261700E+00
     & -.6153110840695750E+00
     & .4465383678265440E+00
     & -.6496120041070140E+00
     & .6496120041073380E+00
     & .6153110840692740E+00
     & .4465383678264870E+00
     & .4465383678262900E+00
     & -.6496120041072340E+00
     & -.6153110840695259E+00
     & .6153110840696310E+00
     & .4465383678261890E+00
     & .6496120041072050E+00
     & .4465383678263000E+00
     & .6496120041072230E+00
     & .6153110840695309E+00
     & .6153110840693370E+00
     & -.4465383678262700E+00
     & -.6496120041074280E+00
     & -.4465383678262480E+00
     & .6496120041073460E+00
     & -.6153110840694390E+00
     & -.6153110840693730E+00
     & -.4465383678265360E+00
     & .6496120041072110E+00
     & -.4465383678262860E+00
     & -.6496120041073030E+00
     & .6153110840694570E+00
     & -.6496120041071211E+00
     & .6153110840697230E+00
     & -.4465383678261830E+00
     & .6496120041071251E+00
     & -.6153110840695510E+00
     & -.4465383678264150E+00
     & .9933631163195030E+00
     & -.1134687281482460E+00
     & -.1882994605477500E-01
     & -.1134687281480350E+00
     & -.1882994605463900E-01
     & .9933631163195290E+00
     & -.9933631163195230E+00
     & .1134687281482040E+00
     & -.1882994605396400E-01
     & -.1882994605390300E-01
     & .9933631163195540E+00
     & -.1134687281479430E+00
     & .1134687281480660E+00
     & -.1882994605432300E-01
     & -.9933631163195320E+00
     & -.1882994605474300E-01
     & -.9933631163195330E+00
     & .1134687281479860E+00
     & .1134687281482190E+00
     & .1882994605448500E-01
     & .9933631163195110E+00
     & .1882994605434400E-01
     & -.9933631163195100E+00
     & -.1134687281482540E+00
     & -.1134687281481780E+00
     & .1882994605424600E-01
     & -.9933631163195210E+00
     & .1882994605448500E-01
     & .9933631163195030E+00
     & .1134687281482870E+00
     & .9933631163195400E+00
     & .1134687281479850E+00
     & .1882994605438200E-01
     & -.9933631163195310E+00
     & -.1134687281480370E+00
     & .1882994605454200E-01"
      read(designtxt,*)(q0x(1:3,i),i=1,72)

      designtxt="
     & .2463988858915690E+00
     & -.7208015696498040E+00
     & .6478677999575010E+00
     & -.7208015696493920E+00
     & .6478677999578860E+00
     & .2463988858917620E+00
     & -.2463988858916820E+00
     & .7208015696496320E+00
     & .6478677999576490E+00
     & .6478677999574370E+00
     & .2463988858916630E+00
     & -.7208015696498290E+00
     & .7208015696498640E+00
     & .6478677999575770E+00
     & -.2463988858911920E+00
     & .6478677999576580E+00
     & -.2463988858916790E+00
     & .7208015696496251E+00
     & .7208015696496560E+00
     & -.6478677999577340E+00
     & .2463988858913890E+00
     & -.6478677999579040E+00
     & -.2463988858914330E+00
     & -.7208015696494890E+00
     & -.7208015696498650E+00
     & -.6478677999573730E+00
     & -.2463988858917270E+00
     & -.6478677999574740E+00
     & .2463988858911660E+00
     & .7208015696499660E+00
     & .2463988858917940E+00
     & .7208015696495070E+00
     & -.6478677999577450E+00
     & -.2463988858914560E+00
     & -.7208015696496660E+00
     & -.6478677999576971E+00
     & -.7935442048021790E+00
     & -.3876287734012690E+00
     & -.4690751848651830E+00
     & -.3876287734013530E+00
     & -.4690751848647940E+00
     & -.7935442048023680E+00
     & .7935442048021710E+00
     & .3876287734015360E+00
     & -.4690751848649750E+00
     & -.4690751848650970E+00
     & -.7935442048020340E+00
     & -.3876287734016680E+00
     & .3876287734016800E+00
     & -.4690751848649880E+00
     & .7935442048020930E+00
     & -.4690751848651100E+00
     & .7935442048021040E+00
     & .3876287734015120E+00
     & .3876287734014250E+00
     & .4690751848652980E+00
     & -.7935442048020350E+00
     & .4690751848650680E+00
     & .7935442048023370E+00
     & -.3876287734010840E+00
     & -.3876287734014910E+00
     & .4690751848649310E+00
     & .7935442048022190E+00
     & .4690751848647840E+00
     & -.7935442048022960E+00
     & .3876287734015120E+00
     & -.7935442048022650E+00
     & .3876287734012240E+00
     & .4690751848650750E+00
     & .7935442048021850E+00
     & -.3876287734018230E+00
     & .4690751848647150E+00
     & .1649450576530030E+00
     & -.9583769097171539E+00
     & .2330382519605870E+00
     & -.9583769097169350E+00
     & .2330382519611260E+00
     & .1649450576535120E+00
     & -.1649450576532380E+00
     & .9583769097170010E+00
     & .2330382519610480E+00
     & .2330382519606680E+00
     & .1649450576535040E+00
     & -.9583769097170480E+00
     & .9583769097171020E+00
     & .2330382519605140E+00
     & -.1649450576534090E+00
     & .2330382519607420E+00
     & -.1649450576532880E+00
     & .9583769097170670E+00
     & .9583769097170990E+00
     & -.2330382519608270E+00
     & .1649450576529820E+00
     & -.2330382519611220E+00
     & -.1649450576532260E+00
     & -.9583769097169860E+00
     & -.9583769097170930E+00
     & -.2330382519606320E+00
     & -.1649450576532930E+00
     & -.2330382519604340E+00
     & .1649450576532610E+00
     & .9583769097171470E+00
     & .1649450576534940E+00
     & .9583769097169650E+00
     & -.2330382519610150E+00
     & -.1649450576534580E+00
     & -.9583769097170310E+00
     & -.2330382519607690E+00
     & .5604842504669760E+00
     & .8132526494836950E+00
     & -.1564529740408340E+00
     & .8132526494836900E+00
     & -.1564529740414460E+00
     & .5604842504668130E+00
     & -.5604842504672400E+00
     & -.8132526494834310E+00
     & -.1564529740412630E+00
     & -.1564529740410300E+00
     & .5604842504670470E+00
     & .8132526494836090E+00
     & -.8132526494838200E+00
     & -.1564529740407260E+00
     & -.5604842504668260E+00
     & -.1564529740410970E+00
     & -.5604842504667780E+00
     & -.8132526494837810E+00
     & -.8132526494836300E+00
     & .1564529740409670E+00
     & .5604842504670350E+00
     & .1564529740412850E+00
     & -.5604842504670530E+00
     & .8132526494835550E+00
     & .8132526494834810E+00
     & .1564529740411510E+00
     & -.5604842504671990E+00
     & .1564529740408810E+00
     & .5604842504669960E+00
     & -.8132526494836720E+00
     & .5604842504668360E+00
     & -.8132526494837370E+00
     & .1564529740411220E+00
     & -.5604842504667400E+00
     & .8132526494838230E+00
     & .1564529740410180E+00
     & .3666300586513120E+00
     & .9220188325509330E+00
     & -.1243530157042820E+00
     & .9220188325508800E+00
     & -.1243530157047620E+00
     & .3666300586512840E+00
     & -.3666300586517610E+00
     & -.9220188325507080E+00
     & -.1243530157046290E+00
     & -.1243530157043770E+00
     & .3666300586515770E+00
     & .9220188325508150E+00
     & -.9220188325509330E+00
     & -.1243530157042030E+00
     & -.3666300586513410E+00
     & -.1243530157045340E+00
     & -.3666300586511110E+00
     & -.9220188325509791E+00
     & -.9220188325508830E+00
     & .1243530157044780E+00
     & .3666300586513720E+00
     & .1243530157046300E+00
     & -.3666300586515370E+00
     & .9220188325507970E+00
     & .9220188325507450E+00
     & .1243530157044630E+00
     & -.3666300586517230E+00
     & .1243530157042990E+00
     & .3666300586515630E+00
     & -.9220188325508309E+00
     & .3666300586512860E+00
     & -.9220188325509230E+00
     & .1243530157044380E+00
     & -.3666300586512290E+00
     & .9220188325509380E+00
     & .1243530157044920E+00
     & -.8046719536517351E+00
     & -.7083625075572700E-01
     & .5894788143650050E+00
     & -.7083625075605800E-01
     & .5894788143650030E+00
     & -.8046719536517070E+00
     & .8046719536519210E+00
     & .7083625075538300E-01
     & .5894788143647920E+00
     & .5894788143647260E+00
     & -.8046719536519410E+00
     & -.7083625075571399E-01
     & .7083625075593900E-01
     & .5894788143647760E+00
     & .8046719536518840E+00
     & .5894788143650180E+00
     & .8046719536517150E+00
     & .7083625075584600E-01
     & .7083625075560100E-01
     & -.5894788143648110E+00
     & -.8046719536518880E+00
     & -.5894788143647840E+00
     & .8046719536518840E+00
     & -.7083625075587500E-01
     & -.7083625075555100E-01
     & -.5894788143649440E+00
     & .8046719536517950E+00
     & -.5894788143649780E+00
     & -.8046719536517590E+00
     & .7083625075567000E-01
     & -.8046719536518360E+00
     & .7083625075619300E-01
     & -.5894788143648110E+00
     & .8046719536517311E+00
     & -.7083625075576400E-01
     & -.5894788143650060E+00"
      read(designtxt,*)(q0x(1:3,i),i=1+72,72+72)

      designtxt="
     & -.8305971377714630E+00
     & -.4813562216367220E+00
     & .2800081831259090E+00
     & -.4813562216367630E+00
     & .2800081831263240E+00
     & -.8305971377712990E+00
     & .8305971377714670E+00
     & .4813562216366280E+00
     & .2800081831260560E+00
     & .2800081831258640E+00
     & -.8305971377713430E+00
     & -.4813562216369560E+00
     & .4813562216370750E+00
     & .2800081831258990E+00
     & .8305971377712620E+00
     & .2800081831260040E+00
     & .8305971377713510E+00
     & .4813562216368590E+00
     & .4813562216366530E+00
     & -.2800081831258590E+00
     & -.8305971377715200E+00
     & -.2800081831260120E+00
     & .8305971377715200E+00
     & -.4813562216365640E+00
     & -.4813562216367410E+00
     & -.2800081831261120E+00
     & .8305971377713840E+00
     & -.2800081831260530E+00
     & -.8305971377713141E+00
     & .4813562216368940E+00
     & -.8305971377713660E+00
     & .4813562216368400E+00
     & -.2800081831259940E+00
     & .8305971377711940E+00
     & -.4813562216370290E+00
     & -.2800081831261780E+00
     & .6225761054046420E+00
     & .2744190843023600E-01
     & -.7820779594393990E+00
     & .2744190843027600E-01
     & -.7820779594394310E+00
     & .6225761054046010E+00
     & -.6225761054049630E+00
     & -.2744190843004500E-01
     & -.7820779594391500E+00
     & -.7820779594391180E+00
     & .6225761054049880E+00
     & .2744190843039700E-01
     & -.2744190843020100E-01
     & -.7820779594392960E+00
     & -.6225761054047740E+00
     & -.7820779594394081E+00
     & -.6225761054046280E+00
     & -.2744190843028900E-01
     & -.2744190843023800E-01
     & .7820779594392210E+00
     & .6225761054048660E+00
     & .7820779594392629E+00
     & -.6225761054048200E+00
     & .2744190843008300E-01
     & .2744190843041900E-01
     & .7820779594392691E+00
     & -.6225761054047980E+00
     & .7820779594394510E+00
     & .6225761054045910E+00
     & -.2744190842992800E-01
     & .6225761054047880E+00
     & -.2744190843038000E-01
     & .7820779594392780E+00
     & -.6225761054045720E+00
     & .2744190842986800E-01
     & .7820779594394680E+00
     & -.9318695934738700E+00
     & .3187128632820320E+00
     & -.1733238919982290E+00
     & .3187128632819440E+00
     & -.1733238919982580E+00
     & -.9318695934738940E+00
     & .9318695934737440E+00
     & -.3187128632820510E+00
     & -.1733238919988710E+00
     & -.1733238919988410E+00
     & -.9318695934738360E+00
     & .3187128632817990E+00
     & -.3187128632819240E+00
     & -.1733238919986170E+00
     & .9318695934738340E+00
     & -.1733238919982450E+00
     & .9318695934739750E+00
     & -.3187128632817140E+00
     & -.3187128632819970E+00
     & .1733238919985150E+00
     & -.9318695934738280E+00
     & .1733238919985010E+00
     & .9318695934738011E+00
     & .3187128632820840E+00
     & .3187128632820890E+00
     & .1733238919985390E+00
     & .9318695934737929E+00
     & .1733238919984430E+00
     & -.9318695934738240E+00
     & -.3187128632820500E+00
     & -.9318695934738650E+00
     & -.3187128632819280E+00
     & .1733238919984480E+00
     & .9318695934738970E+00
     & .3187128632818020E+00
     & .1733238919985030E+00
     & .8838481768527030E+00
     & .2014238044752130E+00
     & .4221858018276850E+00
     & .2014238044757030E+00
     & .4221858018276610E+00
     & .8838481768526020E+00
     & -.8838481768525340E+00
     & -.2014238044755540E+00
     & .4221858018278750E+00
     & .4221858018279100E+00
     & .8838481768524840E+00
     & .2014238044757010E+00
     & -.2014238044754720E+00
     & .4221858018277440E+00
     & -.8838481768526150E+00
     & .4221858018276230E+00
     & -.8838481768526470E+00
     & -.2014238044755860E+00
     & -.2014238044753970E+00
     & -.4221858018278330E+00
     & .8838481768525900E+00
     & -.4221858018279300E+00
     & -.8838481768525230E+00
     & .2014238044754890E+00
     & .2014238044754790E+00
     & -.4221858018276820E+00
     & -.8838481768526431E+00
     & -.4221858018275140E+00
     & .8838481768527690E+00
     & -.2014238044752800E+00
     & .8838481768524760E+00
     & -.2014238044756140E+00
     & -.4221858018279670E+00
     & -.8838481768527100E+00
     & .2014238044755630E+00
     & -.4221858018275020E+00
     & .2042750399564050E+00
     & .7187705698842261E+00
     & .6645604381236631E+00
     & .7187705698843340E+00
     & .6645604381234740E+00
     & .2042750399566370E+00
     & -.2042750399562600E+00
     & -.7187705698842650E+00
     & .6645604381236641E+00
     & .6645604381238100E+00
     & .2042750399561560E+00
     & .7187705698841600E+00
     & -.7187705698843250E+00
     & .6645604381235790E+00
     & -.2042750399563280E+00
     & .6645604381234920E+00
     & -.2042750399563730E+00
     & -.7187705698843930E+00
     & -.7187705698843611E+00
     & -.6645604381235540E+00
     & .2042750399562800E+00
     & -.6645604381235540E+00
     & -.2042750399566620E+00
     & .7187705698842540E+00
     & .7187705698840900E+00
     & -.6645604381238021E+00
     & -.2042750399564320E+00
     & -.6645604381235050E+00
     & .2042750399566820E+00
     & -.7187705698842930E+00
     & .2042750399563760E+00
     & -.7187705698841650E+00
     & -.6645604381237380E+00
     & -.2042750399563670E+00
     & .7187705698845380E+00
     & -.6645604381233370E+00
     & -.8988479274720690E+00
     & .4377008233682800E+00
     & .2214480756061700E-01
     & .4377008233679230E+00
     & .2214480756096300E-01
     & -.8988479274722340E+00
     & .8988479274721820E+00
     & -.4377008233680650E+00
     & .2214480756027000E-01
     & .2214480756031500E-01
     & -.8988479274722930E+00
     & .4377008233678340E+00
     & -.4377008233677660E+00
     & .2214480756062300E-01
     & .8988479274723190E+00
     & .2214480756055900E-01
     & .8988479274721600E+00
     & -.4377008233680940E+00
     & -.4377008233682550E+00
     & -.2214480756032700E-01
     & -.8988479274720880E+00
     & -.2214480756066100E-01
     & .8988479274721030E+00
     & .4377008233682070E+00
     & .4377008233680300E+00
     & -.2214480756060700E-01
     & .8988479274721910E+00
     & -.2214480756073300E-01
     & -.8988479274721950E+00
     & -.4377008233680150E+00
     & -.8988479274722450E+00
     & -.4377008233679080E+00
     & -.2214480756079600E-01
     & .8988479274723130E+00
     & .4377008233677780E+00
     & -.2214480756063400E-01 "
      read(designtxt,*)(q0x(1:3,i),i=1+2*72,72+72+72)
c      do i=1,216
c      write(*,*) i,q0x(1:3,i)
c      enddo
      end
