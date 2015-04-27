      module m_q0p !this is used in mkqg2 defined in this file.
!! These are outputs when getq0p is called.
      integer,allocatable::irr(:)
      real(8),allocatable:: q0i(:,:),wt(:)
      integer:: nq0i !number of offset Gamma point.

      integer,private:: nq0x,nmm
      contains
      subroutine getq0p(newoffsetG,alat,plat,qlat,n1q,n2q,n3q,alp,alpv,
     & ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect)
!!> Get q0p and write 'EPSwklm' file for offset Gamma method.
!! all are inputs. 
      use keyvalue,only: getkeyvalue
      implicit none
      logical,intent(in):: newoffsetG
      integer,intent(in):: n1q,n2q,n3q,nstbz(*),nqbz,nqibz,ngcxx,ngcx(nqbz),ngrp
      real(8),intent(in):: alat,qlat(3,3),alp,alpv(3),plat(3,3)
     & ,qbz(3,nqbz),qibz(3,nqibz),symops(3,3,ngrp)
      integer,intent(in) ::ngvect(3,ngcxx,nqbz)

      integer::nq00ix,nq0x,nx0,nq00i, xyz2lm(3),nnn
      real(8):: xn,www,wgtq0
      logical:: noq0p,timereversal
      real(8),allocatable:: q0x(:,:),wt0(:)
      real(8):: deltaq,deltaq_scale,delta8,delta5,emat(3,3)
      real(8):: pi=4d0*atan(1d0)
      real(8),allocatable:: wti(:),qi(:,:),epinv(:,:,:),cg(:,:,:),matxxl(:,:,:),
     &     dmlx(:,:),cy(:),yl(:),epinvq0i(:,:),wklm(:) !,norq0x(:) !,wqfac(:)
      integer:: bzcase=1,i,iq0i,ifidmlx,iopen,lmxax,lx,lxklm,j,iclose,llxxx
      real(8):: rrr(3),r2s,qxx(3),voltot,wgtq0p,tripl

      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))

!! number of spherical points.
cccc  nq00ix=12 !spherical points
c         nq00ix=12              !spherical points
c     nq00ix=3   !spherical points
      nq00ix=6                  !spherical points
      nx0 = 1
      if(nx0==2) xn=3d0         ! ratio parameter for Q2 and Q1, 
! only effective for nx0=2 case
      nq0x=nq00ix*nx0
c     nq0x=4*nx0
c     if(q0pchoice()/1000==1) then
c     nn1= (q0pchoice()-1000)/10
c     nn2= mod(q0pchoice()-1000,10)
c     nq0x= 4*nn1*nn2
c     print *,' mkqg: q0pchoice nq0x=',q0pchoice(),nq0x
c     endif
c$$$         if(newaniso) then    !feb2012
      nq0x=nq00ix
c$$$         elseif( q0pchoice()<0) then
c$$$c     nq0x = 8*abs (q0pchoice())
c$$$            nq0x = max( (2*abs(q0pchoice()))**3, 8*abs(q0pchoice()))
c$$$         endif
      www = wgtq0p()
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if(noq0p) then
        nq00i=0
        print *,' TestNoQ0P=.true. '
        nq0i=0
      else
        nmm=1
        if(.not.timereversal()) nmm=2
        allocate( q0x(3,nq0x), wt0(nq0x), irr(nq0x), wt(nq0x), q0i(3,nq0x*nmm)) 
c!,norq0x(nq0x))
        if(newoffsetG) then
          deltaq=deltaq_scale()*alat/(2*pi) !dq is 0.01 a.u.
          if(nq00ix==3) then
            nq00i=3
c     q0x(:,1)=(/-deltaq, deltaq,  deltaq/)
c     q0x(:,2)=(/deltaq, -deltaq,  deltaq/)
c     q0x(:,3)=(/deltaq, deltaq,  -deltaq/)
            q0x(:,1)= qlat(:,1)/n1q/2d0*deltaq_scale()
            q0x(:,2)= qlat(:,2)/n2q/2d0*deltaq_scale()
            q0x(:,3)= qlat(:,3)/n3q/2d0*deltaq_scale()
          elseif(nq00ix==6) then
!! six independent direction is required to calculate full dielectric matrix (symmetric -->six components).
            nq00i=6
            q0x(:,1)= qlat(:,1)/n1q/2d0*deltaq_scale()
            q0x(:,2)= qlat(:,2)/n2q/2d0*deltaq_scale()
            q0x(:,3)= qlat(:,3)/n3q/2d0*deltaq_scale()
c                  norq0x(1)=sqrt(sum(q0x(:,1)**2))
c                  norq0x(2)=sqrt(sum(q0x(:,2)**2))
c                  norq0x(3)=sqrt(sum(q0x(:,3)**2))
c before 21dec2012
c     q0x(:,4)= (q0x(:,1)-q0x(:,2))/2d0
c     q0x(:,5)= (q0x(:,2)-q0x(:,3))/2d0
c     q0x(:,6)= (q0x(:,3)-q0x(:,1))/2d0
c     norq0x(4)=sqrt(sum(q0x(:,4)**2))
c     norq0x(5)=sqrt(sum(q0x(:,5)**2))
c     norq0x(6)=sqrt(sum(q0x(:,6)**2))
c     q0x(:,4)= (q0x(:,1)-q0x(:,2))/norq0x(4)*(norq0x(1)+norq0x(2))/2d0
c     q0x(:,5)= (q0x(:,2)-q0x(:,3))/norq0x(5)*(norq0x(2)+norq0x(3))/2d0
c     q0x(:,6)= (q0x(:,3)-q0x(:,1))/norq0x(6)*(norq0x(3)+norq0x(1))/2d0
!! shorter ones. no normalization. dec2012
            if(sum((q0x(:,1)-q0x(:,2))**2)<sum((q0x(:,1)+q0x(:,2))**2)) then
              q0x(:,4)= (q0x(:,1)-q0x(:,2))/2d0
            else
              q0x(:,4)= (q0x(:,1)+q0x(:,2))/2d0
            endif
            if(sum((q0x(:,2)-q0x(:,3))**2)<sum((q0x(:,2)+q0x(:,3))**2)) then
              q0x(:,5)= (q0x(:,2)-q0x(:,3))/2d0
            else  
              q0x(:,5)= (q0x(:,2)+q0x(:,3))/2d0
            endif
            if(sum((q0x(:,3)-q0x(:,1))**2)<sum((q0x(:,3)+q0x(:,1))**2)) then
              q0x(:,6)= (q0x(:,3)-q0x(:,1))/2d0
            else  
              q0x(:,6)= (q0x(:,3)+q0x(:,1))/2d0
            endif
c     q0x(:,1)=(/-deltaq, deltaq,  deltaq/)
c     q0x(:,2)=(/deltaq, -deltaq,  deltaq/)
c     q0x(:,3)=(/deltaq, deltaq,  -deltaq/)
c     q0x(:,4)=(/deltaq, -deltaq,  -deltaq/)
c     q0x(:,5)=(/-deltaq, deltaq,  -deltaq/)
c     q0x(:,6)=(/-deltaq, -deltaq,  deltaq/)
c     nq00i=6
c     q0x(:,1)=(/deltaq, 0d0,    0d0/)
c     q0x(:,2)=(/0d0, deltaq,    0d0/)
c     q0x(:,3)=(/0d0,    0d0, deltaq/)
c     q0x(:,4)=(/0d0,   deltaq, deltaq/)
c     q0x(:,5)=(/deltaq, 0d0,   deltaq/)
c     q0x(:,6)=(/deltaq,deltaq, 0d0/)

          elseif(nq00ix==12) then
!! spherical design des.3.12.5
!! des.3.12.5
            nq00i=12
            delta8=0.850650808352d0*deltaq
            delta5=0.525731112119d0*deltaq
            q0x(:,1)=(/delta8, 0d0, -delta5/)
            q0x(:,2)=(/delta5, -delta8, 0d0/)
            q0x(:,3)=(/0d0,-delta5, delta8/)
            
            q0x(:,4)=(/delta8, 0d0, delta5/)
            q0x(:,5)=(/-delta5,-delta8,0d0/)
            q0x(:,6)=(/0d0,delta5,-delta8/)
            
            q0x(:,7)=(/-delta8,0d0,-delta5/)
            q0x(:,8)=(/-delta5,delta8,0d0/)
            q0x(:,9)=(/0d0,delta5,delta8/)
            
            q0x(:,10)=(/-delta8,0d0,delta5/)
            q0x(:,11)=(/delta5,delta8,0d0/)
            q0x(:,12)=(/0d0,-delta5,-delta8/)
          else
            call rx( 'mkqg: not implemented nq00i')
          endif
          do i=1,nq00i
            write(*,'(" initial q0x=",i3,3f9.3)')i,q0x(:,i)
          enddo
!! invariante dielectoric tensor.
          allocate(epinv(3,3,nq0x))
          call diele_invariant(q0x,nq0x,symops,ngrp,  epinv,q0i,nq0i, wt)
          print *,' nq0x,nmm nq0i=',nq0x,nmm,nq0i
!! == To convert invariant tensor on YL representation (Y00 and Y2m) ==
          lmxax=1
          allocate( cg((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2) )
          allocate( matxxl(3,3,(2*lmxax+1)**2) )
          call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg)
          xyz2lm( 2)=-1         !y
          xyz2lm( 3)= 0         !z
          xyz2lm( 1)= 1         !x
!! matxxl(i,j,L) = \int d\Omega x_i x_j  Y_L(\Omega), where x_i are nomlized.
          do i=1,3
            do j=1,3
              matxxl(i,j,:) = cg(xyz2lm(i)+3,xyz2lm(j)+3,:)*4d0*pi/3d0
!sqrt(4*pi/3) comes from normalization of Y_l=1.
            enddo
          enddo
!! epinv is expanded as
!!   <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat)
          allocate(dmlx(nq0i,9))
          do iq0i=1,nq0i
            do lx=1,9
              dmlx(iq0i,lx)=sum(epinv(:,:,iq0i)*matxxl(:,:,lx))
            enddo
          enddo
c$$$  !! check xxxxxxxxxxxxxxxxxxxx
c$$$  do lx=5,9
c$$$  do i=2,4
c$$$  do j=2,4
c$$$  write(*,"(' l1 l2 l= cg=',3i3,f9.5)")i-1,j-1,lx-7,cg(i,j,lx)
c$$$  enddo
c$$$  enddo
c$$$  write(*,*)
c$$$  enddo
c$$$  do lx=5,9
c$$$  do i=1,3
c$$$  do j=1,3
c$$$  write(*,"(' matxxl l1 l2 l= cg=',3i3,f9.5)")i,j,lx,matxxl(i,j,lx)
c$$$  enddo
c$$$  enddo
c$$$  write(*,*)
c$$$  enddo
c$$$  do lx=1,1
c$$$  do i=1,3
c$$$  do j=1,3
c$$$  write(*,"(' matxxl l1 l2 l= cg=',3i3,f9.5)")i,j,lx,matxxl(i,j,lx)
c$$$  enddo
c$$$  enddo
c$$$  write(*,*)
c$$$  enddo
c$$$  do iq0i=1,nq0i
c$$$  do lx=1,9
c$$$  write(*,"(' iq0i lx dmlx=',2i3,f9.3)")iq0i,lx,dmlx(iq0i,lx)
c$$$  enddo
c$$$  enddo

!! Test for one r vector as for <ehat|epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat)
!! === generate YL for a test vector rrr (rrr is ehat above).====
          lx=2
          allocate(cy((lx+1)**2),yl((lx+1)**2))
          call sylmnc(cy,lx)
          rrr=(/.5d0,-.1d0,-0.7d0/) !test data
          rrr=rrr/sqrt(sum(rrr**2))
          call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
!! ===== check (rrr*emat*rrr = sum(dmlx* YL)
c     do lm=1,9; write(*,"('r lm=',3f8.3,i4,' ylm=',f8.3)") rrr,lm,cy(lm)*yl(lm) ;   enddo
          write(*,*)
          write(*,"(  ' test: r=',3f10.5)") rrr
!!
          do iq0i=1,nq0i
            write(*,"(' test: ylm   expansion=',i3,f10.5)") 
     &       iq0i,sum(dmlx(iq0i,:)*cy(:)*yl(:))
            emat=epinv(:,:,iq0i)
            write(*,"(' test: epinv expansion=',i3,f10.5)") 
     &       iq0i,sum(rrr*matmul(emat,rrr))
          enddo
          write(*,*)
          allocate( epinvq0i(nq0i,nq0i))
          do i=1,nq0i
          do j=1,nq0i         !epinvq0i= <q0i/|q0i|| epinv(:,:,iq0j)|q0i/|q0i|>
           epinvq0i(i,j)=sum(q0i(:,i)*matmul(epinv(:,:,j),q0i(:,i)))/sum(q0i(:,i)**2)
          enddo
          enddo
          deallocate(cy,yl)
          lxklm=6         !this is used for inversion procedure in hx0fp0.sc.m.f
          nnn=n1q*n2q*n3q
          allocate(wklm((lxklm+1)**2)) !wklm-->Klm in Comp.Phys. Comm 176(1007)1-13
          call getwklm(alat,voltot,plat,qlat,alp,qbz,nqbz,ngcx,ngcxx,ngvect, lxklm,n1q,n2q,n3q,
     o     wklm)                !,wqfac)
!! Takao think anisotropic treatment in Computer Phys. Comm 176(1007)1-13 
!! (in our version with auxially function) can be numerically prorematic.
!! We only keep wklm only up to l=2. (I observed high wklm(lm) components 
!! are largely negative --->it may cause numerical error).
!! 
!! From the begining, we can only excpect "virtual convergence on Nz" for 
!! NxNyNz for Si100 slab model in the paper.
!! (I still not understand why it does not show divergent behevior in the anisotropic case).
!!
c     print *,' set wklm=0 for l>2. But lxklm(for inversion of epsioln)=',lxklm
          do i=1,(lxklm+1)**2
            if(abs(wklm(i))>1d-6 ) write(6,'(" mkqg: l lm Wklm=",2i3,f9.4)')llxxx(i),i,wklm(i)
          enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     ! spherical design des.3.12.5 check. Because of angular momentum synsesize,
!     ! des.3.12.5 gives correct normalization of product up to l=2 (lm<=9)
c$$$  deallocate(cy,yl,q0x)
c$$$  nq00i=12
c$$$  allocate(cy((lxklm+1)**2),yll((lxklm+1)**2,nq00i))
c$$$  tpiba  = 2d0*pi/alat
c$$$  call sylmnc(cy,lxklm)
c$$$  allocate( q0x(3,nq00i) )
c$$$  delta8=0.850650808352d0 !*deltaq
c$$$  delta5=0.525731112119d0 !*deltaq
c$$$  q0x(:,1)=(/delta8, 0d0, -delta5/)
c$$$  q0x(:,2)=(/delta5, -delta8, 0d0/)
c$$$  q0x(:,3)=(/0d0,-delta5, delta8/)
c$$$  q0x(:,4)=(/delta8, 0d0, delta5/)
c$$$  q0x(:,5)=(/-delta5,-delta8,0d0/)
c$$$  q0x(:,6)=(/0d0,delta5,-delta8/)
c$$$  q0x(:,7)=(/-delta8,0d0,-delta5/)
c$$$  q0x(:,8)=(/-delta5,delta8,0d0/)
c$$$  q0x(:,9)=(/0d0,delta5,delta8/)
c$$$  q0x(:,10)=(/-delta8,0d0,delta5/)
c$$$  q0x(:,11)=(/delta5,delta8,0d0/)
c$$$  q0x(:,12)=(/0d0,-delta5,-delta8/)
c$$$  do iq=1,nq00i
c$$$  qg(1:3) = q0x(:,iq)
c$$$  call sylm(qg/sum(qg**2),yll(:,iq),lxklm,r2s) !spherical factor Y( q+G )
c$$$  c        print *,' qg for yll=',iq, qg
c$$$  c        print *,' yll=',cy(1:20)*yll(1:20,iq)
c$$$  enddo

!! normalization check
c     do lm1=1,(3+1)**2
c     do lm2=lm1,(3+1)**2
c     aaa=sum(cy(lm1)*cy(lm2)*yll(lm1,:)*yll(lm2,:))/12d0*4d0*pi
c     if(abs(aaa)>1d-6) write(6,"('ylm*ylm=',2i3,d13.5)")lm1,lm2,aaa
c     enddo
c     enddo
c     do lm1=1,(5+1)**2
c     aaa=sum(cy(lm1)*yll(lm1,:))
c     if(abs(aaa)>1d-6) write(6,"('ylm*ylm=',i3,d13.5)")lm1,aaa
c     enddo
c     stop ' xxxxxxxxx spherical normalization xxxxxxxx'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$  !! test function generation
c$$$  deallocate(cy,yl)
c$$$  allocate(cy((lxklm+1)**2),yl((lxklm+1)**2))
c$$$  tpiba  = 2d0*pi/alat
c$$$  call sylmnc(cy,lxklm)
c$$$  do iq=1,nqbz
c$$$  funa(:,iq)=0d0
c$$$  do ig=1,ngcx(iq)
c$$$  qg(1:3) = tpiba * (qbz(1:3,iq)+ matmul(qlat, ngvect(1:3,ig,iq)))
c$$$  qg2     = sum(qg(1:3)**2)
c$$$  alpqg2= alp* qg2
c$$$  call sylm(qg/sqrt(qg2),yl,lxklm,r2s) !spherical factor Y( q+G )
c$$$  funa(:,iq) = funa(:,iq) + exp(-alpqg2)*(1d0+1/(qg2+5d0))/qg2*cy(:)*yl(:) !cy*yl =Y_L(qg/|qg|)
c$$$  enddo
c$$$  enddo
c$$$  c what is wtrue???
c$$$  ccccccccccccccccccccccccc
c$$$  do lm=1,(lxklm+1)**2
c$$$  wsumau(lm) = sum(funa(lm,2:nqbz))/dble(nqbz)
c$$$  c        write(6,"('  wsum fnua=',i3,8f10.5)") lm,wsumau(lm)
c$$$  if(lm==1) then
c$$$  write(*,"('lm l wklm wtrue wsum wsummesh',2i3,4f12.8)")
c$$$  &     lm,llxxx(lm),wklm(lm), wtrue00,wklm(lm)+wsumau(lm), wsumau(lm)!,wklm(lm)+wsumau(lm)-wtrue00
c$$$  else
c$$$  write(*,"('lm l wklm wtrue wsum wsummesh',2i3,4f12.8)")
c$$$  &     lm,llxxx(lm),wklm(lm), 0d0,  wklm(lm)+wsumau(lm), wsumau(lm)
c$$$  endif
c$$$  enddo
c$$$  stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
c$$$  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ifidmlx = iopen('EPSwklm',0,-1,0)
          write(ifidmlx) nq0i,lxklm
          write(ifidmlx) dmlx, epinv(:,:,1:nq0i),epinvq0i
          write(ifidmlx) wklm
          ifidmlx = iclose('EPSwklm')
        else
          call setq0_2(bzcase, alat, voltot,plat, qlat,alpv, qbz, nstbz, nqbz,
     i     ngcx, ngcxx, ngvect,   nq0x,nx0,xn,n1q,n2q,n3q,
     o     q0x,wt0,nq00i)
!     ! ... inequivalent q0x points ...
          nq0i=0
          call q0irre(qibz,0,q0x,wt0,nq00i,symops,ngrp,   q0i,nq0i,wt,plat,.false.,0,irr)
        endif
        deallocate(irr)
      endif
      write(6,"('i wt q0i=',i3,f16.7,2x,3d23.15)")(i,wt(i),q0i(1:3,i),i=1,nq0i)
      end subroutine getq0p
      end module m_q0p


!! ----------------------------------------------------------------------
      subroutine mkQG2(alat,plat,symops,ngrp,n1q,n2q,n3q,iq0pin,
     &     QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc)
      use m_get_bzdata1,only:  getbzdata1,
     &     nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm,nqibz_r,
     &     qbz,wbz,qibz,wibz,
     &     qbzw, qbasmc, dq_bzcase2, qibz_r,
     &     idtetf, ib1bz, idteti,
     &     irk, nstar, nstbz,
     &     qbzm, qbzwm
      use keyvalue,only: getkeyvalue
      use m_q0p,only: getq0p, irr,q0i,wt,nq0i
!! == Make required q and G in the expantion of GW. ==
!!     |q+G| < QpGcut_psi for eigenfunction psi.
!!     |q+G| < QpGcut_Cou for coulomb interaction
!!
!! OUTPUT
!!     file handle= ifiqg,  which contains q and G points for eigenfunction psi. --> QGpsi
!!     file handle= ifiqgc, which contains q and G points for Coulomb            --> QGcou
!!
!!     QGpsi(ifiqg), QGcou(ifiqgc), Q0P are written.
!!     See the end of console output.
!! ---------------------------------------------------
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,nnn,ngcxx,
     &     ngrp,i,j,iq,iq00,ngp,ngpmx,ngc,ngcmx,nqnum,iq0pin,
     &     nline,nlinemax,ifsyml,iqq,is,nk,ix,nqnumx,i1,ifkpt
      real(8)  :: plat(3,3),qlat(3,3),q(3),dummy,qp(3), 
     &     QpGcut_psi, QpGcut_Cou,QpGcut,alpv(3),q0smean,sumt,alp,
     &     volum,voltot,pi,q0(3),qlat0(3,3), alat,tripl,
     &     symops(3,3,ngrp),xx,qqx(3),alpm
      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:),
     &     ngpn(:),ngcn(:),ngvect(:,:,:),ngcx(:), nqq(:)   
      real(8),allocatable :: q0x(:,:), 
     &     qq(:,:),qq1(:,:),qq2(:,:),qqm(:,:)
      real(8) :: qbas(3,3),vol,ginv(3,3),aaa,det,dq(3),www
      integer(4) :: mxkp,ifiqibz,iqibz,ifigwin,mtet(3),nm1,nm2,nm3
      logical ::tetrai,tetraf,tetra_hsfp0    
      integer(4):: ifbz
      integer(4):: bzcase=1
c     logical :: readgwinput
      integer(4):: nqnumm,ifiqmtet,verbose,q0pchoice,nn1,nn2,ifiqbz,iqbz !,auxfunq0p
      real(8)::aaij,bbij,wgtq0p, qdum(6)
      logical:: qbzreg

      logical :: qreduce ,qreduce0 ,ibzqq
      real(8),allocatable:: qsave(:,:),   qmin(:,:),qmax(:,:),qbzrr(:,:)
      integer:: imx,ifinin,il,nq0i0,ni,nq0i00,imx0
      integer,allocatable :: ndiv(:),ngvecprev(:,:,:),ngveccrev(:,:,:)

      real(8):: ddq(3)
      logical :: offmesh=.false. ,offmeshg=.false.
      logical :: regmesh=.false. ,regmeshg=.false. ,  timereversal

      logical ::anyq                , caca,debug=.true. !,newaniso
      real(8),allocatable:: qany(:,:)
      integer(4):: nany,ifqpnt,ret,imxc,nnn3(3),imx0c,imx11(1,1)
      real(8):: deltaq,delta5,delta8,deltaq_scale!=1d0/3.0**.5d0 

      integer:: nqi,ifix,ig,iq0i,lm
      real(8),allocatable:: wti(:),qi(:,:)
      integer:: ifidmlx,iclose,iopen !,ifiwqfac

      integer:: llxxx,lm1,lm2
      real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:)
      real(8)::volinv,wtrue00,qg(3),alpqg2,qg2,tpiba
      character*7:: q0pf        !nov2012
      integer:: dummyia(1,1),iimx,irradd,nmax
      real(8):: epstol=1d-8,tolq=1d-6,qx(3),qxx(3)
      logical :: newoffsetG !july2014
      real(8),allocatable:: wt0(:)
c------------------------------------------------
      print *,' mkqg2: '
      qreduce0 = qreduce()
      newoffsetG=.true. !newaniso()
      if(iq0pin == 101) then
         iq0pin=1
         newoffsetG=.false. !for old oldset Gamma case
      endif
!!     band case --- iq0pin == 3
!!     nqq(is),qq1(1:3,is),qq2(1:3,is),is =1,nline
      if(iq0pin == 3) then
         qreduce0=.false.
         nlinemax = 50
         allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
         ifsyml = 3001
         open(ifsyml,file='SYML')
         nline = 0
         do
            nline = nline + 1
            read(ifsyml,*,err=601,end=601) 
     &           nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
         enddo
 601     continue
         close(ifsyml)
         nline = nline - 1
         write(6,"(/' Symmetry lines:'/' points',12x,'start',22x,'end')")
         do is=1,nline
            write(6,"(i6,2x,3f8.4,2x,3f8.4)") 
     &           nqq(is),(qq1(i,is),i=1,3),(qq2(i,is),i=1,3)
         enddo
         nqnumx = sum(nqq(1:nline))
         allocate( qq(1:3,nqnumx),irr(nqnumx) )
         iqq = 0
         do is = 1,nline
            nk = nqq(is)
            do iq=1,nk
               xx = 0d0
               if(nk>1) xx=(iq-1d0)/(nk-1d0)
               qqx = xx*qq2(1:3,is)+(1d0-xx)*qq1(1:3,is)
               iqq = iqq + 1
               qq(1:3,iqq) = qqx
               write (6,"('  q=',3f7.3)") qq(1:3,iqq)
            enddo
         enddo
         nqnum = iqq
         write (6,"(' Total number of q-points:',i5/)") nqnum
         call dinv33x (plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
         goto 2001
      endif

!! --- Ordinary case --- iq0pin == 1 or 2
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      call dinv33x (plat,qlat)  !it was dinv33(plat,1,qlat) by Ferdi
      call getkeyvalue("GWinput","delta",aaa)
      if(aaa<0d0) then
         print * ,'READ GWIN_V2 --->: tetrahedron method for x0'
         tetraf=.true.
      else
         print * ,'READ GWIN_V2 --->: not use tetrahedron method for x0'
         tetraf=.false.
      endif
      tetrai = .true.           !used in heftet tetra_hsfp0()

!! ---  See indxk in index.f \in genbz2 \in genallc_v2
      call dinv33(qlat,0,ginv,det)
      write(6,*)'=== plat ==='
      write(6,"(3d23.15)") plat
      write(6,*)'=== qlat ==='
      write(6,"(3d23.15)") qlat
      write(6,*)'=== ginv==='
      write(6,"(3f9.4)") ginv
      do i=1,3
         do j=1,3
            aaij=sum(qlat(:,i)*plat(:,j))
            bbij=sum(qlat(:,i)*ginv(j,:))
            if(verbose()>=40) print *,' i j aaij bbij', i,j,aaij,bbij
            if(i==j) then
               if(abs(aaij-1d0) >1d-10) call rx( 'bug 1 qg4gw')
               if(abs(bbij-1d0) >1d-10) call rx( 'bug 2 qg4gw')
            else
               if(abs(aaij) >1d-10) call rx( 'bug 3 qg4gw')
               if(abs(bbij) >1d-10) call rx( 'bug 4 qg4gw')
            endif
         enddo
      enddo
      mtet=(/1,1,1/)
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
      if(sum(abs(mtet))<3) then
         print *, ' we use regular meshing for tetrahedron scheme '
      endif

!! getbzdata1 allocate dataset in m_get_bzdata1
      print *
      print *,'goto getbzdata1...'
      call cputid (0)
      call getbzdata1(bzcase,plat,qlat,ginv,n1q,n2q,n3q
     &     ,symops,ngrp,tetrai,tetraf,mtet) !all are inputs.
!! Write BZDATA
      ifbz = 6661
      open (ifbz, file='BZDATA')
      write(ifbz,"(10i10)")  nqbz,nqibz, nqbzw, ntetf, nteti,ngrp,nqibz_r
      write(ifbz,"(10i10)")  n1q,n2q,n3q
      print *,' writing BZDATA...'
      call cputid (0)
      call rwbzdata(ifbz,-1,
     &     ngrp,qlat, ginv, qbasmc,
     i     qbz, wbz, nstbz,  nqbz, 
     i     qibz,wibz, nstar,irk,  nqibz,
     i     idtetf, ntetf, qbzw,ib1bz, nqbzw,
     i     idteti, nteti,dq_bzcase2, qibz_r,nqibz_r)
      close(ifbz)
      write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2
!! Write QIBZ
      write(6,*)' qibz are written in QIBZ file...'
      ifiqibz = 6661
      open (ifiqibz, file='QIBZ') !write q-points in IBZ.
      write(ifiqibz,"(i10)") nqibz
      do iqibz = 1,nqibz
         write(ifiqibz,"(3d24.16,3x,d24.16)") qibz(1:3,iqibz),wibz(iqibz)
      enddo
      close(ifiqibz)
!! Write QBZ
      ifiqbz = 6661
      open (ifiqbz, file='QBZ') !write q-points in IBZ.
      write(ifiqbz,"(i10)") nqbz
      do iqbz = 1,nqbz
         write(ifiqbz,"(3d24.16,3x,d24.16)") qbz(1:3,iqbz)
      enddo
      close(ifiqbz)
!!  Write to file KPNTin1BZ
      ifkpt = 335
      open(ifkpt,file='KPTin1BZ.mkqg.chk')
      write(ifkpt,*)"  qbz --> shoten(qbz)"
      do      i1 = 1,nqbz
         call shorbz(qbz(1,i1),qp,qlat,plat)
         write (ifkpt,"(1x,i7,4f10.5,'   ',3f10.5)") 
     &        i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
      end do
      close (ifkpt)
      write(6,*) ' --- TOTAL num of q =',nqbz
      write(6,*)
      write(6,"( '  ngrp  = ',i3)")ngrp
      write(6,'("  qibz=",i6,3f12.5)')(i,qibz(1:3,i),i=1,min(10,nqibz))
      write(6,*)" ... QIBZ is written in QIBZ file ..."

      call getkeyvalue("GWinput","alpha_OffG",alp,default=-1d60)
      alpv(:)=alp
      if(alp==-1d60) then
       call getkeyvalue("GWinput","alpha_OffG_vec",alpv,3,default=(/-1d50,0d0,0d0/))
       if(alpv(1)==-1d50) then
         call rx( ' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput')
       endif
      endif
      print *
      print *,' alpv=',alpv
      print *
      alpm = minval(alpv)
      if(alpm<=0d0) call rx( 'alpha_offG or alpha_offG_vec <=0')
!! === Large if start. ===
      if(iq0pin==1) then   ! --- get q0x (offsetted q=0 point) ----------------------
c  QpGcut = 15d0/alpm  !a.u. !exp( -  alp * QpGcut) !alp * QpGcut = 10
c  QpGcut = sqrt(25d0/alpm) !a.u. !exp( -  alp * QpGcut**2) !alp * QpGcut**2 = 22
c  QpGcut = sqrt(100d0/alp)
c  QpGcut = sqrt(150d0/alp)
c  QpGcut = sqrt(300d0/alp)
         QpGcut = sqrt(25d0/alpm) !a.u. !exp( -alp*QpGcut**2) !alp * QpGcut**2 = 22
         allocate( ngcx(nqbz) )
         ngcx=1
         do iq = 1, nqbz
            q   = qbz(1:3,iq)
            call getgv2(alat,plat,qlat,q, QpGcut, 1, ngcx(iq),  dummyia)
         enddo
         ngcxx = maxval(ngcx)
         allocate( ngvect(3,ngcxx,nqbz) )
         print *,' goto getgv2: ngcxx=',ngcxx
         do iq = 1, nqbz
            q  = qbz(1:3,iq)
            call getgv2( alat,plat,qlat, q, QpGcut, 2, 
     &           ngcx(iq), ngvect(1:3,1:ngcx(iq),iq) )
         enddo
!! all inputs
         call getq0p(newoffsetG,alat,plat,qlat,n1q,n2q,n3q,alp,alpv,
     i    ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect)
         open (1101,file='Q0P')
         www = wgtq0p()
         write(1101,"(2i5,' !nq0i iq0pin' )") nq0i,iq0pin
         write(1101,"(d24.16,3x, 3d24.16)" ) (wt(i),q0i(1:3,i),i=1,nq0i)
         close(1101)

      elseif(iq0pin==2) then
         call getkeyvalue("GWinput","QforEPSIBZ",ibzqq,default=.false.)
         if(ibzqq) then
            write(6,*)'=== Find QforEPSIBZ=on === '
            nq0i= nqibz
            allocate( q0i(3,nq0i) )
            q0i = qibz
         else
            write(6,*)'==== Readin <QforEPS>or<QforEPS> in GWinput === '
            call getkeyvalue("GWinput","<QforEPS>", unit=ifinin,status=nq0i00,errstop='off')
            nq0i00 =max(nq0i00,0)
            if(nq0i00>0) close(ifinin)
            print *,' end of reaing QforEPS nq0i00',nq0i00,ifinin

            call getkeyvalue("GWinput","<QforEPSL>",unit=ifinin,status=nq0i0,errstop='off')
            nq0i0  =max(nq0i0,0)
            print *,' end of reaing QforEPSL nq0i0',nq0i0,ifinin
            if(nq0i0>0) then
               allocate( ndiv(nq0i0) )
               do i=1,nq0i0
                  read(ifinin,*) qdum(1:6), ndiv(i)
               enddo
               nq0i = nq0i00 + sum(ndiv)
               close(ifinin)
            else
               nq0i = nq0i00
            endif
            if(nq0i <=0) call rx( 'There are neither <QforEPS> nor <QforEPS>.')
            allocate( q0i(3,nq0i) )
            print *,' nq0i=',nq0i
            if(nq0i00>0) then
               call getkeyvalue("GWinput","<QforEPS>",unit=ifinin,status=nq0i00)
               do i=1,nq0i00
                  read (ifinin,*) q0i(1:3,i)
                  write (6,"('<QforEPS> ' 3f12.8)") q0i(:,i)
               enddo
               close(ifinin)    !25jan2006
            endif
            if(nq0i0>0) then
               call getkeyvalue("GWinput","<QforEPSL>",unit=ifinin,status=nq0i0)
               allocate( qmin(3,nq0i0), qmax(3,nq0i0) )
               do i=1, nq0i0
               read(ifinin,*)qmin(:,i), qmax(:,i), ndiv(i)
               write(6,"('<QforEPSL>',3f12.8,2x,3f12.8,i5)")qmin(:,i),qmax(:,i),ndiv(i)
               enddo
               close(ifinin)
               ni = nq0i00
               do il=1, nq0i0
                  do i=1, ndiv(il)
                     q0i(:,i+ni)= qmin(:,il)+ (qmax(:,il)-qmin(:,il))/ndiv(il) * i
                  enddo
                  ni= ni + ndiv(il)
               enddo
               deallocate(qmin,qmax,ndiv)
            endif
         endif
         allocate( wt(nq0i) )
         wt = 0d0
         open (1101,file='Q0P')
         write(1101,"(2i5,a)") nq0i,iq0pin, " !nq0i iq0pin ---"//
     &        "This is readin Q0P from GWinput <QforEPS> ---"
         write(1101,"(d24.16,3x, 3d24.16)") (wt(i),q0i(1:3,i),i=1,nq0i)
         close(1101)
      endif
      print *,' end fo writing Q0P'
      call cputid (0)

!! === AnyQ mechanism. ===
      call getkeyvalue("GWinput","AnyQ",anyq,default=.false.)
      if(anyq.and.iq0pin==1) then
         print *,'AnyQ (read <QPNT> section =T'
!!     read q-points and states
         call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
         call readx   (ifqpnt,10)
         call readx   (ifqpnt,100)
         call readx   (ifqpnt,100)
         read (ifqpnt,*) nany
         print *,'  nany=',nany
         allocate(qany(3,nany))
         do ix=1,nany
            read (ifqpnt,*) i, qany(:,ix)
            write(6,'(i3,3f13.6)') ix,qany(:,ix)
         enddo
         nany =ix-1
         write(6,*)" Anyq mode: nany=",nany
         allocate(qsave(3,nq0i+nany))
         qsave(:,    1 :nq0i)     = q0i (:,1:nq0i)
         qsave(:,nq0i+1:nq0i+nany)= qany(:,1:nany)
         nq0i = nq0i+nany
         deallocate(q0i)
         allocate(q0i(3,nq0i))
         q0i=qsave
         deallocate(qsave)
         close(ifqpnt)
      endif
!! qreduce0 will reduce the number of q points in anyway.
      if(.not.timereversal().and.iq0pin==1) then
         write(6,*)" timereversal==off : add -Q0P points"
         do iq=1,nq0i
            q0i(:,iq+nq0i)= -q0i(:,iq)
         enddo
         nq0i=nq0i*2
      endif
!! ---  Get qq; q-point for eigenfunctions'
      regmeshg = qbzreg()       !Gamma mesh based on regular mesh
      offmeshg = .not.qbzreg()  !Gamma mesh based on off-regular mesh

      if(iq0pin==2) then        !this is just for dielectric case
         regmesh = qbzreg()
         offmesh = (.not.qbzreg())
      else
         regmesh = .true.
         offmesh = (bzcase==2) .or. (.not.qbzreg())
      endif
      print *,' regmesh offmeshg=', regmesh,regmeshg !regular, regular+shifted
      print *,' offmesh offmeshg=', offmesh,offmeshg !offregmesh, offregular+shifted
      nqnum = 0
      if(regmesh) nqnum = nqnum + nqbz
      if(offmesh) nqnum = nqnum + nqbz
!! We check wether all q0i \in qbz or not. 
      allocate( qq(1:3,nqnum),irr(nqnum) )
      qq(1:3,1:nqbz) = qbz(1:3,1:nqbz) 
      do iq0i=1,nq0i
         do iq=1,nqbz
            if(sum(abs(q0i(:,iq0i)-qq(:,iq)))<tolq) goto 2112
            call rangedq( matmul(ginv,q0i(:,iq0i)-qq(:,iq)), qx)
            if(sum(abs(qx))< tolq) goto 2112
         enddo
         goto 2111
 2112    continue
         qq(:,iq) = q0i(:,iq0i) !replaced with equivalent q0i.
      enddo
      print *,' --- We find all q0i in qbz. Skip qreduce.'
      goto 2001
 2111 continue
!!
      deallocate(qq,irr)
      nqnum = nqnum + nqbz*nq0i
      nqnum = nqnum + 1         !add Gamma
      nqnum = nqnum + nq0i         !add q0i
      allocate( qq(1:3,nqnum),irr(nqnum) )
!!!   nqnum = ix !=n1q*n2q*n3q*(nq0i+1) ! total q number
C...  Regular mesh (qq are the same in the both bzcase)
      ddq = 0d0
c      if(bzcase()==2) ddq=  dq_bzcase2
      allocate(qbzrr(1:3,nqbz))
      do iq = 1, nqbz
         qbzrr(1:3,iq) = qbz(1:3,iq) + ddq
! qbzzr is regular mesh including gamma.
! qbzzr are identical in the both bzcase.
!--- Sorry about a little confusion; qbz is dependent on bzcase().---
      enddo
      ix = 0
      if(regmesh) then
         qq(1:3,1:nqbz) = qbzrr(1:3,1:nqbz)
         ix = ix+ nqbz
      endif
C...  Off Regular mesh.
      if(offmesh) then
         do iq = 1, nqbz
            ix = ix+1
            qq(1:3,ix) = qbzrr(1:3,iq) - dq_bzcase2
         enddo
      endif
      nnn   =  ix               !n1q*n2q*n3q!      if(offmesh) nnn = 2*n1q*n2q*n3q
      print *,' nnn=',nnn       !This is the number to calcualte Vxc
C...  Shifted mesh
      if(regmeshg) then
         do iq00 = 1, nq0i
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbzrr(1:3,iq) +  q0i(1:3,iq00)
            enddo
         enddo
      endif
      if(offmeshg) then
         do iq00 = 1, nq0i
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbzrr(1:3,iq) - dq_bzcase2 + q0i(1:3,iq00)
            enddo
         enddo
      endif
!! Add offg anyway.
      do iq00 = 1, nq0i
         ix = ix+1
         qq(1:3,ix) = q0i(1:3,iq00)
      enddo
!! Add gamma point for anycase.
      ix=ix+1
      qq(1:3,ix)=0d0
!! Get qqm; q point for eigenvalues. Saved to Qmtet. Not so much used now... Maybe we need check.
      if(sum(abs(mtet))/=3) then
         nqnumm= nqbzm * (nq0i+1)
         allocate( qqm(1:3,nqnumm) )
         ix=0
         do iq00 = 1, 1 + nq0i
            do iq   = 1, nqbzm
               ix = ix+1
               if(iq00==1) then
                  qqm(1:3,ix) = qbzm(1:3,iq)
               else
                  qqm(1:3,ix) = q0i(1:3,iq00-1) + qbzm(1:3,iq)
               endif
            enddo
         enddo
         ifiqmtet=501
         open(ifiqmtet, file='Qmtet')
         write(ifiqmtet,"(i10)") nqnumm
         do iq=1,nqnumm
            write(ifiqmtet,"(3d24.16)") qqm(1:3,iq)
         enddo
         close(ifiqmtet)
         deallocate(qqm)
      endif

!! == Reduce number of q. by modulo of translational symmetry ==
      if( qreduce0 ) then       !June2005 remove duplicated qqm
         print *,'goto qqsave nq0i nqnum',nq0i,nqnum
         call cputid (0)
         nmax=nq0i+nqnum
         allocate(qsave(3,nmax)) !,qsavel(nmax))
         imx=0
         if(iq0pin /=1) then
            do iq=1,nq0i
               call qqsave(q0i(1:3,iq),nmax,ginv,qsave,imx)
            enddo
c         call qqsave2(q0i,nq0i,nmax,ginv,qsave,qsavel,imx)
         endif
         do iq=1,nqnum
            call qqsave(qq(1:3,iq),nmax,ginv,qsave,imx)
         enddo
c         call qqsave2(qq,nqnum,nmax,ginv,qsave,qsavel,imx)
         nqnum = imx
         qq(:,1:imx)=qsave(:,1:imx)
         deallocate(qsave)
      endif
 2001 continue

!! ------------------------------------------
      if(allocated(wt0)) deallocate(wt0)
      allocate(wt0(nqnum+nq0i),qi(3,nqnum+nq0i),wti(nqnum+nq0i))
      wt0=1d0
!! Set irreducible k-point flag. irr=1 for (irredusible point) flag, otherwise =0. 
!! irr(iq)=1 for irreducile qq(:,iq), iq=1,nqnum
      call q0irre(qibz,nqibz,qq,wt0,nqnum,symops,ngrp, qi,nqi,wti,plat,.true.,0,irr)
c      write(6,*)' --- Irreducible qpoint --- nqi=',nqi
c     &     '             qx             qy             qz          count'
c      write(6,"(i3,3f15.7,2x,f12.4)")(i,qi(1:3,i),wti(i),i=1,nqi)
      allocate(ngpn(nqnum), ngcn(nqnum))
      if(debug) write(6,*) ' *** q vector in 1st BZ + q0 shift.   ngp ***'
      imx=0
      imxc=0
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q                  !use q for Q0P points.
         if(iq0pin==1) then
            do i=1,nq0i
               if(sum(abs(q0i(:,i)-qxx))<1d-6) then
                  qxx=(/0d0,0d0,0d0/)
                  exit
               endif
            enddo
         endif
         ngpn(iq)=1
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi,1,ngpn(iq),imx11) !imx11
         imx0=imx11(1,1)
         if(imx0>imx) imx=imx0
         ngcn(iq)=1
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou,1,ngcn(iq),imx11) !imx11 to avoid warning.
         imx0c=imx11(1,1)
         if(imx0c>imxc) imxc=imx0c
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngpn(iq) !,ngcn(iq,iq00)
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngcn(iq) !,ngcn(iq,iq00)
      enddo
      ngpmx = maxval(ngpn)
      ngcmx = maxval(ngcn)
      write(ifiqg ) nqnum,ngpmx,QpGcut_psi,nqbz,nqi,imx,nqibz
      write(ifiqgc) nqnum,ngcmx,QpGcut_cou,nqbz,nqi,imxc
!! :nqi:   The number of irreducible points (including irr. of offset points). irr=1.
!! ::       We calcualte eigenfunction and Vxc for these points.
!! :nqnum: total number of q points.
!! :imx:   to allocate ngvecprev as follows.
      print *,' number of irrecucible points nqi=',nqi
      print *,' imx nqnum=',imx,nqnum
      write(6,*) ' --- Max number of G for psi =',ngpmx
      write(6,*) ' --- Max number of G for Cou =',ngcmx

!! === Get G vector for given q. ===
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) )       !mar2012takao
      allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc) ) !mar2012takao
      ngvecprev=9999
      ngveccrev=9999
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q
         q0pf=''
         do i=1,nq0i
            if(sum(abs(q0i(:,i)-qxx))<1d-6) then
               if(iq0pin==1) then
                  qxx=(/0d0,0d0,0d0/)
               endif
               q0pf=' <--Q0P'
               exit
            endif
         enddo
         ngp = ngpn(iq)
         ngc = ngcn(iq)
         write(6,"(' iq=',i8,' q=',3f7.3,' ngp ngc= ',2i6,' irr.=',i2,a)")
     &        iq, q, ngp, ngc, irr(iq),q0pf
         allocate( ngvecp(3,max(ngp,1)), ngvecc(3,max(ngc,1)) ) !larger than 1. Aug2005
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi, 2, ngp,  ngvecp) ! for eigenfunction
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou, 2, ngc,  ngvecc) ! for Coulomb
         write (ifiqg) q, ngp, irr(iq)
!! Reverse table of ngvecp ---> ngvecprev
         do ig = 1,ngp
            nnn3 = ngvecp(1:3, ig) 
            ngvecprev( nnn3(1), nnn3(2),nnn3(3)) = ig
         enddo
         write (ifiqg)  ngvecp,ngvecprev !ngvecprev is added on mar2012takao
         do ig = 1,ngc
            nnn3 = ngvecc(1:3, ig) 
            ngveccrev( nnn3(1), nnn3(2),nnn3(3)) = ig
         enddo
         write (ifiqgc) q, ngc
         write (ifiqgc) ngvecc,ngveccrev
         deallocate(ngvecp,ngvecc)
      enddo
      deallocate(ngpn,ngcn,ngvecprev,ngveccrev)
      if(iq0pin==1) deallocate(ngvect)
      if(debug) print *,'--- end of mkqg ---'
      end

!! #############################################################################
      double precision function tripl(a,b,c)
!! == tripl (determinant of 3x3 matrix) ==
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(0:2),b(0:2),c(0:2)
c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c     .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c     ... g77 needs this rewriting
      tmp = 0.d0
      do i = 0,2
         j = mod(i + 1, 3)
         k = mod(i + 2, 3)
         tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
      enddo
      tripl = tmp
      end

!! #############################################################################
      subroutine genbzdata1_xxxxxxx(plat,n1q,n2q,n3q,symops,ngrp,
     &     tetraf, tetrai,
     o     qbz, wbz,              nqbz, 
     o     qibz,wibz, nstar,irk,  nqibz,
     o     idtetf, ntetf, qbzw,ib1bz, nqbzw,
     o     idteti, nteti )
!     ! == genbzdata1_xxxxxxx: Generate BZ mesh data reuired for GW ==
C     iput
C     plat,qlat,n1q,n2,n3q,symops,ngrp,
C     tetraf: if true, this generate the BZ data for full-1stBZ tetrahedron method
C     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C     tetrai: if true, this generate the BZ data for tetrahedron method in irreducible BZ.
C     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C     output
C     o  qbz, wbz,              nqbz,
C     o  qibz,wibz, nstar,irk,  nqibz,
C     ouptput(optional)
C     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C     o  idteti, nteti
C-------------------------------------------------------------------------
      implicit none
      integer(4):: nqbz,n1q,n2q,n3q,ntetf,nteti,nqbzw,iqbz,ifbz
     &     ,nqibz,iqibz,itet,ngrp
      real(8) :: plat(3,3),qbas(3,3),ginv(3,3),vol,symops(3,3,ngrp)
     &     ,qbasmc(3,3)
      logical tetraf,tetrai
      integer(4),allocatable:: ipq(:),iw1(:)
      real(8):: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:)
      integer(4) :: idtetf(:,:),ib1bz(:),idteti(:,:),irk(:,:)
     &     ,nstar(:)
      end

!! #############################################################################
      integer(4) function llxxx(ilm)
      integer(4),parameter :: lmx=50
      integer(4),save:: lla((lmx+1)**2)
      logical:: init=.true.
Cstop2rx 2013.08.09 kino      if(ilm>(lmx+1)**2) stop 'll: ilm too large'
      if(ilm>(lmx+1)**2) call rx( 'll: ilm too large')
      if(init) then
         do l=0,lmx
            lini= l**2 + 1
            lend=(l+1)**2
            lla(lini:lend)=l
         enddo
      endif
      llxxx = lla(ilm)
ccccccccccccccccccccccc
c     do i=1,(lmx+1)**2
c     write(117,"(' i ll=',2i5)")i,lla(i)
c     enddo
c     stop 'test end'
cccccccccccccccccccccccc
      return
      end
