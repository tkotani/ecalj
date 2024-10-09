!> Q0P (offset Gamma points) and thir weight generator.
module m_q0p   
  !! All of them are outputs when getq0p is called.
  !! Read GWinput
  !! Write Q0P: Q0P contains offset Gamma, or k points specified in GWinput
  !!       EPSwklm: Sperical integration weitht. But currently only l=0 is used to keep positive definitness.
  ! note
  !! TKkotani think anisotropic treatment in Computer Phys. Comm 176(1007)1-13
  !! (in our version with auxially function) can be numerically prorematic.
  !! We now take wklm for l=0 only because I observed high wklm(lm) components
  !! can be negative--->it may cause numerical error.
  !!
  !! From the begining, we can only excpect "virtual convergence on Nz" for
  !! NxNyNz for Si100 slab model in the paper.
  !! (I still not understand why it does not show divergent behevior in the anisotropic case).
  !!
  implicit none
  integer,public,protected:: nq0i=0,nq0iadd=0,nany=0 ! Number of Q0P !,nq0itrue=0
  integer,public,protected,allocatable:: ixyz(:)  ! ixyz(1:nq0i+nq0iadd) q0i for x,y,z directions
  real(8),public,allocatable,protected:: q0i(:,:),wt(:) ! Q0P and its weight.
  integer,public,allocatable,protected:: epslgroup(:) !EPSwklm
  integer,public,protected:: lxklm=0                  !EPSwklm
  real(8),public,allocatable,protected:: epinv(:,:,:),wklm(:), dmlx(:,:), epinvq0i(:,:) !EPSwklm
  public:: Getallq0p
  private
contains
  subroutine getallq0p(iq0pin,alat,plat,qlat,nnn,alp,alpv,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,lnq0iadd) !! All arguments are input.
    use m_keyvalue,only: getkeyvalue
    intent(in)         iq0pin,alat,plat,qlat,nnn,alp,alpv,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,lnq0iadd
    integer:: iq0pin !    logical:: newoffsetG
    integer:: nnn(3),nstbz(*),nqbz,nqibz,ngcxx,ngcx(nqbz),ngrp !n1q,n2q,n3q,
    !      integer::ngvect(3,ngcxx,nqbz)
    real(8):: alat,qlat(3,3),alp,alpv(3),plat(3,3),qbz(3,nqbz),qibz(3,nqibz),symops(3,3,ngrp)
    integer::nq00ix,nx0,nq00i, xyz2lm(3)
    real(8):: xn !,www,wgtq0
    logical:: noq0p,timereversal
    real(8),allocatable:: q0x(:,:),wt0(:)
    real(8):: deltaq,deltaq_scale,delta8,delta5,emat(3,3)
    real(8),allocatable:: wti(:),qi(:,:),cg(:,:,:),matxxl(:,:,:), cy(:),yl(:) !,norq0x(:) !,wqfac(:)
    integer:: bzcase=1,i,iq0i,ifidmlx,lmxax,lx,j
    real(8):: rrr(3),r2s,qxx(3),voltot,tripl
    integer,allocatable::irrx(:),ngvect(:,:,:)

    real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:),alpqg2,qg2,tpiba,wtrue00
    real(8):: qg(3), qdum(6),alpm,QpGcut,q(3)
    integer:: ig,lm,iq
    integer:: nmm !not output
    integer:: ifi0,ifi00,il,ix,ni,ifinin
    integer:: nq0i0,nq0i00,ifqpnt,ret
    integer,allocatable :: ndiv(:)
    real(8),allocatable:: qsave(:,:),   qmin(:,:),qmax(:,:)
    real(8),allocatable:: qany(:,:)
    logical:: anyq,ibzqq,lnq0iadd,unita
    !      integer,allocatable:: epslgroup(:)
    integer:: dummyia(1,1)
    real(8),parameter:: pi=4d0* atan(1d0)
    real(8):: tpioa
    tpioa=2*pi/alat
    alpm = minval(alpv)
    if(alpm<=0d0) call rx( 'alpha_offG or alpha_offG_vec <=0')
    !      if(debug) write(6,"('  alpm=',3f12.6)") alpm
    if(iq0pin==1) then
       !! Determine G vectors for q points set by getgv2
       QpGcut = sqrt(25d0/alpm) !a.u. !exp( -alp*QpGcut**2) !alp * QpGcut**2 = 22
       ngcx=1
       do iq = 1, nqbz
          q   = qbz(1:3,iq)
          call getgv2(alat,plat,qlat,q, QpGcut, 1, ngcx(iq), dummyia)
       enddo
       ngcxx = maxval(ngcx)
       allocate( ngvect(3,ngcxx,nqbz) )
       do iq = 1, nqbz
          q  = qbz(1:3,iq)
          call getgv2(alat,plat,qlat,q, QpGcut, 2, ngcx(iq), ngvect(1:3,1:ngcx(iq),iq) )
       enddo
       !! Normal mode. all inputs
       call getq0p(alat,plat,qlat,nnn,alp,alpv, &
       ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect,lnq0iadd)
       !! Get Q0P from GWinput
    elseif(iq0pin==2) then
       call getkeyvalue("GWinput","QforEPSunita",unita,default=.false.)
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
          allocate(epslgroup(nq0i))
          epslgroup=0
          allocate( q0i(3,nq0i) )
          print *,' nq0i=',nq0i
          if(nq0i00>0) then
             call getkeyvalue("GWinput","<QforEPS>",unit=ifinin,status=nq0i00)
             do i=1,nq0i00
                read (ifinin,*) q0i(1:3,i)
                if(unita) q0i(:,i)=q0i(:,i)/tpioa !2023-5-19fixed
                write (6,"('<QforEPS> ' 3f12.8)") q0i(:,i)
             enddo
             close(ifinin)    !25jan2006
          endif
          if(nq0i0>0) then
             call getkeyvalue("GWinput","<QforEPSL>",unit=ifinin,status=nq0i0)
             allocate( qmin(3,nq0i0), qmax(3,nq0i0) )
             do i=1, nq0i0
                read(ifinin,*)qmin(:,i), qmax(:,i), ndiv(i)
                if(unita) qmin(:,i)=qmin(:,i)/tpioa
                if(unita) qmax(:,i)=qmax(:,i)/tpioa
                write(6,"('<QforEPSL>',3f12.8,2x,3f12.8,i5)")qmin(:,i),qmax(:,i),ndiv(i)
             enddo
             close(ifinin)
             ni = nq0i00
             do il=1, nq0i0
                do i=1, ndiv(il)
                   q0i(:,i+ni)= qmin(:,il)+ (qmax(:,il)-qmin(:,il))/ndiv(il) * i
                enddo
                epslgroup(ni+1:ni+ndiv(il)) = il !!group of QforEPSL
                ni= ni + ndiv(il)
             enddo
             deallocate(qmin,qmax,ndiv)
          endif
       endif
       allocate( wt(nq0i),source=0d0 )
    elseif(iq0pin==3) then
       nq0i=5
       allocate(epslgroup(nq0i))
       epslgroup=0
       allocate( q0i(3,nq0i) )
       print *,' nq0i=',nq0i
       do i=1,nq0i
          q0i(:,i)= [0d0,0d0,0.01d0]*i
       enddo   
       allocate( wt(nq0i),source=0d0 )
    endif
    print *,' end fo writing Q0P'
    call cputid (0)
    !! Timereversal may require q0i. Anyway, qreduce0 will reduce the number of q points by symops.
    if( .NOT. timereversal() .AND. iq0pin==1) then
       write(6,*)" timereversal==off : add -Q0P points"
       do iq=1,nq0i
          q0i(:,iq+nq0i)= -q0i(:,iq)
       enddo
       nq0i=nq0i*2
    endif
    !! === AnyQ mechanism. === q0i is extended. nq0i/=nq0itrue
    call getkeyvalue("GWinput","AnyQ",anyq,default=.false.)
    if(anyq .AND. iq0pin==1) then
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
       !nq0itrue=nq0i !nov2015
       !nq0i = nq0i+nany
       deallocate(q0i)
       allocate(q0i(3,nq0i+nany))
       q0i=qsave
       deallocate(qsave)
       close(ifqpnt)
    else
       !nq0itrue=nq0i !nov2015
    endif
  end subroutine getallq0p
  !! ==================================================================
  subroutine getq0p(alat,plat,qlat,nnn,alp,alpv,ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect,lnq0iadd) !! All arguments are input.
    use m_keyvalue,only: getkeyvalue
    intent(in)::    alat,plat,qlat,nnn,alp,alpv,ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect,lnq0iadd
    !! output
    !!  q0i (offset Gamma point)
    !!  wt  (weight of q0i)
    !!  EPSwklm (file)
    !!
    !! In addition, we write a file EPSwklm, which is key for new offset Gamma method.
    !! deltaq_scale() given by Q0Pchoice in GWinput change the of offset Gamma method.
    logical :: lnq0iadd
    integer :: nnn(3),nstbz(*),nqbz,nqibz,ngcxx,ngcx(nqbz),ngrp
    real(8) :: alat,qlat(3,3),alp,alpv(3),plat(3,3) ,qbz(3,nqbz),qibz(3,nqibz),symops(3,3,ngrp)
    real(8):: deltaq,deltaq_scale,delta8,delta5,emat(3,3), pi=4d0*atan(1d0)
    real(8):: rrr(3),r2s,qxx(3),voltot,tripl
    integer ::ngvect(3,ngcxx,nqbz),nq00ix,nx0,nq00i, xyz2lm(3)!,nnnt
    real(8):: xn !,www,wgtq0
    logical:: noq0p,timereversal
    real(8),allocatable:: q0x(:,:),wt0(:)
    real(8),allocatable:: wti(:),qi(:,:),cg(:,:,:),matxxl(:,:,:), cy(:),yl(:) !,norq0x(:) !,wqfac(:)
    integer:: bzcase=1,i,iq0i,ifidmlx,lmxax,lx,j
    integer,allocatable::irrx(:)
    real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:),alpqg2,qg2,tpiba,wtrue00
    real(8):: qg(3)
    integer:: ig,lm,iq, nq0x,nmm ,iq0,iqx!not output
    real(8):: tolw
    write(6,"(    'getq0p: offset Gamma method by T.Kotani,JPSJ83,094711,(2014)')")
    write(6,"(a)")'        We use only wklm(lm=0) so as to keep positive definite weight(from mar2016)'
    voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
    ! number of spherical points.
    nq00ix=6 !we only use this now 2024-9-17. See previous codes for other experimental cases)
    nx0 = 1
    if(nx0==2) xn=3d0 ! ratio parameter for Q2 and Q1,
    nq0x=nq00ix
    !    call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
    !    if(noq0p) then
    !       nq00i=0
    !       print *,' TestNoQ0P=.true. '
    !       nq0i=0
    !    else
    nmm=1
    if( .NOT. timereversal()) nmm=2
    allocate( q0x(3,nq0x), wt0(nq0x), irrx(nq0x), wt(nq0x), q0i(3,nq0x*nmm))
    deltaq=deltaq_scale()*alat/(2*pi) !dq is 0.01 a.u.
    ! six independent direction is required to calculate full dielectric matrix (symmetric -->six components).
    nq00i=6
    do i=1,3
       q0x(:,i)= qlat(:,i)/nnn(i)/2d0*deltaq_scale()
    enddo
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
    !! invariante dielectoric tensor.
    allocate(epinv(3,3,nq0x))
    call diele_invariant(q0x,nq0x,symops,ngrp,  epinv,q0i,nq0i, wt)
    write(6,"(a,3i5)")'  nq0x,nmm nq0i=',nq0x,nmm,nq0i
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
          ! qrt(4*pi/3) comes from normalization of Y_l=1.
       enddo
    enddo
    !! epinv is expanded as ! <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat)
    allocate(dmlx(nq0i,9))
    do iq0i=1,nq0i
       do lx=1,9
          dmlx(iq0i,lx)=sum(epinv(:,:,iq0i)*matxxl(:,:,lx))
       enddo
    enddo
    !! Test for one r vector as for <ehat|epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat)
    !! === generate YL for a test vector rrr (rrr is ehat above).====
    lx=2
    allocate(cy((lx+1)**2),yl((lx+1)**2))
    call sylmnc(cy,lx)
    rrr=(/.5d0,-.1d0,-0.7d0/) !test data
    rrr=rrr/sqrt(sum(rrr**2))
    call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
    !! ===== check (rrr*emat*rrr = sum(dmlx* YL)
    !     do lm=1,9; write(*,"('r lm=',3f8.3,i4,' ylm=',f8.3)") rrr,lm,cy(lm)*yl(lm) ;   enddo
    write(*,"('  test: sample r=',3f10.5)") rrr
    do iq0i=1,nq0i
       write(*,"('  test: ylm   expansion=',i3,f10.5)") iq0i,sum(dmlx(iq0i,:)*cy(:)*yl(:))
       emat=epinv(:,:,iq0i)
       write(*,"('  test: epinv expansion=',i3,f10.5)") iq0i,sum(rrr*matmul(emat,rrr))
    enddo
    allocate( epinvq0i(nq0i,nq0i))
    do i=1,nq0i
       do j=1,nq0i         !epinvq0i= <q0i/|q0i|| epinv(:,:,iq0j)|q0i/|q0i|>
          epinvq0i(i,j)=sum(q0i(:,i)*matmul(epinv(:,:,j),q0i(:,i)))/sum(q0i(:,i)**2)
       enddo
    enddo
    deallocate(cy,yl)
    !          lxklm=6         !this is used for inversion procedure in hx0fp0.sc.m.f
    !          nnnt=n1q*n2q*n3q
    allocate(wklm((lxklm+1)**2)) !wklm-->Klm in Comp.Phys. Comm 176(1007)1-13
    call getwklm(alat,voltot,plat,qlat,alp,qbz,nqbz,ngcx,ngcxx, ngvect,lxklm,nnn(1),nnn(2),nnn(3),  wklm)
    !    print *,' set wklm=0 for l>2. But lxklm(for inversion of epsioln)=',lxklm
    do i=1,(lxklm+1)**2
       if(abs(wklm(i))>1d-6 ) write(6,'("  l lm Wklm=",2i3,f9.4)')llxxx(i),i,wklm(i)
    enddo
    ! spherical design des.3.12.5 check. Because of angular momentum synsesize, des.3.12.5 gives correct normalization of product up to l=2 (lm<=9)    --->removed 2024-9-17
    !$$$!! test wklm
    !$$$          if(allocated(cy)) deallocate(cy,yl)
    !$$$          allocate(cy((lxklm+1)**2),yl((lxklm+1)**2),funa((lxklm+1)**2,nqbz))
    !$$$          tpiba  = 2d0*pi/alat
    !$$$          call sylmnc(cy,lxklm)
    !$$$          do iq=1,nqbz
    !$$$            funa(:,iq)=0d0
    !$$$            do ig=1,ngcx(iq)
    !$$$              qg(1:3) = tpiba * (qbz(1:3,iq)+ matmul(qlat, ngvect(1:3,ig,iq)))
    !$$$              qg2     = sum(qg(1:3)**2)
    !$$$              alpqg2= alp* qg2
    !$$$              call sylm(qg/sqrt(qg2),yl,lxklm,r2s) !spherical factor Y( q+G )
    !$$$              if(qg2<1d-6) cycle !for iq=1 remove divergent part
    !$$$              funa(:,iq) = funa(:,iq) + exp(-alpqg2)/qg2*cy(:)*yl(:) !cy*yl =Y_L(qg/|qg|)
    !$$$            enddo
    !$$$          enddo
    !$$$          allocate(wsumau((lxklm+1)**2))
    !$$$          do lm=1,(lxklm+1)**2
    !$$$            wsumau(lm) = sum(funa(lm,2:nqbz))/dble(nqbz)
    !$$$c            if(abs(funa(lm,1))>1d-8) write(6,*)' funa1=',lm,funa(lm,1)
    !$$$c        write(6,"('  wsum fnua=',i3,8f10.5)") lm,wsumau(lm)
    !$$$            if(lm==1) then
    !$$$c              write(*,"('lm l wklm wtrue wsum wsummesh',2i3,4f12.8)")
    !$$$c     &         lm,llxxx(lm),wklm(lm), wtrue00,wklm(lm)+wsumau(lm), wsumau(lm) !,wklm(lm)+wsumau(lm)-wtrue00
    !$$$            else
    !$$$              write(*,"('lm l wsumau+wklm+w0const = wtotal ',2i3,15f12.8)")
    !$$$     &         lm,llxxx(lm), wklm(lm), wsumau(lm), funa(lm,1)/dble(nqbz), wklm(lm)+wsumau(lm)+funa(lm,1)/dble(nqbz)
    !$$$            endif
    !$$$          enddo
    !$$$          stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
    deallocate(irrx)
    write(6,"('  i wt(i) q0i(:,i)=',i3,f16.7,2x,3d23.15)")(i,wt(i),q0i(1:3,i),i=1,nq0i)
    !! Add q0
    nq0iadd=0
    allocate(ixyz(nq0i+3)) !nq0i+3 is large enough
    ixyz=0
    if(lnq0iadd) then
       if(nq00ix/=6 ) then
          call rx('mkqg: we assumes q0x(:,i)= qlat(:,i)/nnn(i)/2d0*deltaq_scale() for 1=1,3' )
       endif
       do iq=1,3           !we assume q0x(:,i)= qlat(:,i)/nnn(i)/2d0*deltaq_scale() for i=1,3
          do iq0=1,nq0i
             if(sum(abs(q0x(:,iq)-q0i(:,iq0)))<tolw()) then
                ixyz(iq0)=iq
                goto 1011
             endif
          enddo
          nq0iadd = nq0iadd+1
          ixyz( nq0i+nq0iadd)= iq
          q0i(:,nq0i+nq0iadd)=q0x(:,iq)
1011      continue
       enddo
    endif
  end subroutine getq0p
  integer function llxxx(ilm)
    integer,parameter :: lmx=50
    integer,save:: lla((lmx+1)**2)
    logical:: init=.true.
    integer:: ilm,l,lend,lini
    if(ilm>(lmx+1)**2) call rx( 'll: ilm too large')
    if(init) then
       do l=0,lmx
          lini= l**2 + 1
          lend=(l+1)**2
          lla(lini:lend)=l
       enddo
    endif
    llxxx = lla(ilm)
    return
  end function llxxx
end module m_q0p
