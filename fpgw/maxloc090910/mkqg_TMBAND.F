      subroutine mkQG2_TMBAND(alat,plat,symops,ngrp,n1q,n2q,n3q,iq0pin,
     &                   QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc)
      use m_get_bzdata1
      use keyvalue
ctm, 081030, from mkqg_RSBAND.F
C- Make required q and G in the expantion of GW.
c   |q+G| < QpGcut_psi for eigenfunction psi.
c   |q+G| < QpGcut_Cou for coulomb interaction
c output 
c    file handle= ifiqg,  which contains q and G points for eigenfunction psi.
c    file handle= ifiqgc, which contains q and G points for Coulomb
c    Q0P 
c
co    QIBZ       : q points in the IBZ of iq0pin==1 or 2
c---------------------------------------------------
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,nnn,ngcxx,nx0,
     &  ngrp,nq0i,i,j,iq,iq00,ngp,ngpmx,ngc,ngcmx,nqnum,iq0pin,
     &  nline,nlinemax,ifsyml,iqq,is,nk,ix,nqnumx,i1,ifkpt
      real(8)  :: plat(3,3),qlat(3,3),q(3),dummy,qp(3), 
     &  QpGcut_psi, QpGcut_Cou,QpGcut,alpv(3),q0smean,sumt,xn,alp,
     &  volum,voltot,pi,q0(3),qlat0(3,3), alat,tripl,
     &  symops(3,3,ngrp),xx,qqx(3),alpm
      integer(4),allocatable:: ngvecp(:,:), ngvecc(:,:),
     &  ngpn(:),ngcn(:),ngvect(:,:,:),ngcx(:), nqq(:)
      real(8),allocatable :: q0x(:,:), q0i(:,:),wt(:),wt0(:),
     &  qq(:,:),qq1(:,:),qq2(:,:),qqm(:,:)
      real(8) :: qbas(3,3),vol,ginv(3,3),aaa,det,dq(3),www
      integer(4) :: mxkp,ifiqibz,iqibz,ifigwin,nq00i,mtet(3),nm1,nm2,nm3
      logical ::tetrai,tetraf,tetra_hsfp0    ,noq0p
      integer(4):: ifbz
      integer(4):: bzcase
c      logical :: readgwinput
      integer(4):: nqnumm,ifiqmtet,verbose, nq0x,q0pchoice,nn1,nn2 !,auxfunq0p
      real(8)::aaij,bbij,wgtq0p, qdum(6)
      logical:: qbzreg

      logical :: qreduce ,qreduce0 ,ibzqq
      real(8),allocatable:: qsave(:,:),   qmin(:,:),qmax(:,:),qbzrr(:,:)
      integer:: imx,ifinin,il,nq0i0,ni,nq0i00
      integer(4),allocatable :: ndiv(:)

      real(8):: ddq(3)
      logical :: offmesh=.false. ,offmeshg=.false.
      logical :: regmesh=.false. ,regmeshg=.false. ,  timereversal
      integer(4)::nmm,nqnum0
c------------------------------------------------
C    modified by T. Miyake, 081111
C    Band-structure mode by R. Sakuma 
C    May 2007
      integer(4) :: nqnum_rsband
      real(8),allocatable :: qq_rsband(:,:)
c------------------------------------------------
      qreduce0 = qreduce()
c      qreduce0 = .false.

C RS: read q-points for band-plotting
C      if(iq0pin == 3) then
      if(.true.) then

c band case --- iq0pin == 3
c nqq(is),qq1(1:3,is),qq2(1:3,is),is =1,nline
c        qreduce0=.false.
        nlinemax = 50
        allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
        ifsyml = 3001
        open(ifsyml,file='SYML')
        nline = 0
        do 
        nline = nline + 1
        read(ifsyml,*,err=601,end=601) 
     &    nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
        enddo
  601   continue
        close(ifsyml)
        nline = nline - 1
        write(6,"(/' Symmetry lines:'/' points',12x,'start',22x,'end')")
        do is=1,nline
          write(6,"(i6,2x,3f8.4,2x,3f8.4)") 
     &     nqq(is),(qq1(i,is),i=1,3),(qq2(i,is),i=1,3)
        enddo
        nqnumx = sum(nqq(1:nline))
        allocate( qq_rsband(1:3,nqnumx) )
        iqq = 0
        do is = 1,nline
         nk = nqq(is)
        do iq=1,nk
          xx = 0d0
          if(nk>1) xx=(iq-1d0)/(nk-1d0)
          qqx = xx*qq2(1:3,is)+(1d0-xx)*qq1(1:3,is) 
c          if(iqq>1 ) then
c            if(abs(sum(qqx-qq_rsband(:,iqq)))<1d-10) cycle
c          endif
          iqq = iqq + 1
          qq_rsband(1:3,iqq) = qqx
          write (6,"('  q=',3f7.3)") qq_rsband(1:3,iqq)
        enddo
        enddo
        nqnum_rsband = iqq
        write (6,"(' Total number of q-points:',i5/)") nqnum_rsband
c        call dinv33x (plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
c        goto 2001
      endif
 
C --- Ordinary case --- iq0pin == 1 or 2
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      call dinv33x (plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
c      print *,'plat qlat='
c      print *, plat,qlat
c----------------------

c      if(readgwinput()) then
       call getkeyvalue("GWinput","delta",aaa) 
c      else 
c       ifigwin = 788
c       open (ifigwin, file='GWIN_V2')
c       do i=1,6;read(ifigwin,*) ;enddo
c       read(ifigwin,*) aaa
c      endif

      if(aaa<0d0) then
         print * ,'READ GWIN_V2 --->: tetrahedron method for x0'
         tetraf=.true.
      else
         print * ,'READ GWIN_V2 --->: not use tetrahedron method for x0'
         tetraf=.false.
      endif
      tetrai = .true. !used in heftet tetra_hsfp0()

C ... Get BZ data --- case1
c      call dinv33 (plat,1,qlat,vol)
c      call dinv33 (qlat,1,ginv,vol)

c--- See indxk in index.f \in genbz2 \in genallc_v2
      call dinv33(qlat,0,ginv,det)
      write(6,*)'=== plat ==='
      write(6,"(3f9.4)") plat
      write(6,*)'=== qlat ==='
      write(6,"(3f9.4)") qlat
      write(6,*)'=== ginv==='
      write(6,"(3f9.4)") ginv
      do i=1,3
      do j=1,3
        aaij=sum(qlat(:,i)*plat(:,j))
        bbij=sum(qlat(:,i)*ginv(j,:))
        if(verbose()>=40) print *,' i j aaij bbij', i,j,aaij,bbij
        if(i==j) then
           if(abs(aaij-1d0) >1d-10) stop 'bug 1 qg4gw' 
           if(abs(bbij-1d0) >1d-10) stop 'bug 2 qg4gw' 
        else
          if(abs(aaij) >1d-10) stop 'bug 3 qg4gw' 
          if(abs(bbij) >1d-10) stop 'bug 4 qg4gw' 
        endif
      enddo
      enddo


      mtet=(/1,1,1/)
c      if(readgwinput()) then
        call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
        if(sum(abs(mtet))<3) then
          print *, ' we use regular meshing for tetrahedron scheme '
        endif
c      endif  
c--- getbzdata1 allocate all these data
c Output lists of getbzdata1
c        integer(4):: nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm,nqibz_r
c        real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:),wibz(:)
c     &   ,qbzw(:,:) ,qbasmc(3,3),dq_bzcase2(3),qibz_r(:,:)
c        integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &   ,irk(:,:),nstar(:),nstbz(:)
c        real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
      call getbzdata1(bzcase(),plat,qlat,ginv,n1q,n2q,n3q
     &  ,symops,ngrp,tetrai,tetraf,mtet) !all are inputs.

C ... Write BZDATA
      ifbz = 6661
      open (ifbz, file='BZDATA') 
      write(ifbz,"(10i10)")  nqbz,nqibz, nqbzw, ntetf, nteti,ngrp,nqibz_r
      write(ifbz,"(10i10)")  n1q,n2q,n3q
      call rwbzdata(ifbz,-1,
     &  ngrp,qlat, ginv, qbasmc,
     i  qbz, wbz, nstbz,  nqbz, 
     i  qibz,wibz, nstar,irk,  nqibz,
     i  idtetf, ntetf, qbzw,ib1bz, nqbzw,
     i  idteti, nteti,dq_bzcase2, qibz_r,nqibz_r)
      close(ifbz)

      write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2

C ... Write QIBZ
      ifiqibz = 6661
      open (ifiqibz, file='QIBZ') !write q-points in IBZ.
      write(ifiqibz,"(i10)") nqibz
      do iqibz = 1,nqibz
        write(ifiqibz,"(3d24.16,3x,d24.16)") qibz(1:3,iqibz),wibz(iqibz)
      enddo
      close(ifiqibz)

C ... Write to file KPNTin1BZ
      ifkpt = 335  
      open(ifkpt,file='KPTin1BZ.mkqg.chk') 
      do      i1 = 1,nqbz
        call shorbz(qbz(1,i1),qp,qlat,plat)
        write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5)") 
     &   i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
      end do
      close (ifkpt)
c      deallocate(wbz)
      write(6,*) ' --- TOTAL num of q =',nqbz

c ...
      write(6,*)
      write(6,"( '  ngrp  = ',i3)")ngrp
      write(6,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)
c
      
      call getkeyvalue("GWinput","alpha_OffG",alp,default=-1d60)
      alpv(:)=alp
      if(alp==-1d60) then
        call getkeyvalue("GWinput","alpha_OffG_vec",alpv,3,default=(/-1d50,0d0,0d0/))
        if(alpv(1)==-1d50) then
          print *,' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput' 
          stop    ' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput' 
        endif
      endif
      print *
      print *,' alpv=',alpv
      print *
      alpm= minval(alpv)
      if(alpm<=0d0) stop 'alpha_offG or alpha_offG_vec <=0'
      if(iq0pin==1) then ! --- get q0x (offsetted q=0 point) ----------------------
c        if(auxfunq0p()==1) then
c          QpGcut = 15d0/alpm  !a.u. !exp( -  alp * QpGcut) !alp * QpGcut = 10
c        else
          QpGcut = sqrt(25d0/alpm)  !a.u. !exp( -  alp * QpGcut**2) !alp * QpGcut**2 = 22
c        QpGcut = sqrt(100d0/alp)
c        QpGcut = sqrt(150d0/alp)
c        QpGcut = sqrt(300d0/alp)
c        endif
        allocate( ngcx(nqbz) )
        do iq = 1, nqbz
        q   = qbz(1:3,iq)  
        call getgv2(alat,plat,qlat,q, QpGcut, 1, ngcx(iq),  dummy)
        enddo
        ngcxx = maxval(ngcx)
        allocate( ngvect(3,ngcxx,nqbz) )
        print *,' goto getgv2: ngcxx=',ngcxx
        do iq = 1, nqbz
          q  = qbz(1:3,iq)  
          call getgv2( alat,plat,qlat, q, QpGcut, 2, 
     &          ngcx(iq), ngvect(1:3,1:ngcx(iq),iq) ) 
        enddo
c 
        nx0 = 1 
        if(nx0==2) xn  = 3d0 !ratio parameter for Q2 and Q1, only effective for nx0=2 case 
        nq0x=6*nx0

c        if(q0pchoice()/1000==1) then
c          nn1= (q0pchoice()-1000)/10 
c          nn2= mod(q0pchoice()-1000,10) 
c          nq0x= 4*nn1*nn2
c          print *,' mkqg: q0pchoice nq0x=',q0pchoice(),nq0x
c        endif
        if( q0pchoice()<0) then
c          nq0x = 8*abs (q0pchoice())
          nq0x = max( (2*abs(q0pchoice()))**3, 8*abs(q0pchoice()))
        endif

        www=wgtq0p()

        call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
        if(noq0p) then
           nq00i=0
           print *,' TestNoQ0P=.true. '
           nq0i=0
        else
          nmm=1
          if(.not.timereversal()) nmm=2
          allocate( q0x(3,nq0x), wt0(nq0x), wt(nq0x), q0i(3,nq0x*nmm))
          call setq0_2(bzcase(), alat, voltot,plat, qlat,alpv, qbz, nstbz, nqbz,
     i         ngcx, ngcxx, ngvect,   nq0x,nx0,xn,n1q,n2q,n3q,
     o         q0x,wt0,nq00i)
c ... inequivalent q0x points ----
          call q0irre(q0x,wt0,nq00i,symops,ngrp,   q0i,nq0i,wt)
          write(6,"(i3,f16.7,2x,3f16.7)")(i,wt(i),q0i(1:3,i),i=1,nq0i)
        endif
        print *,' xxxxxxx nq0i=',nq0i

        open (1101,file='Q0P')
        www=wgtq0p()
        if(bzcase()==2) then
          write(1101,"(2i5,' !nq0i iq0pin :wgtq0p()=', d23.15)") nq0i,iq0pin,www
        else
          write(1101,"(2i5,' !nq0i iq0pin' )") nq0i,iq0pin
        endif
        write(1101,"(d24.16,3x, 3d24.16)" )
     &     (wt(i),q0i(1:3,i),i=1,nq0i)
c
c        wtx = 0d0; q0xx(1) = 0.01d0; q0xx(2:3) =0d0
c        write(101,"(d24.16,3x, 3d24.16,a)" ) wtx,q0xx(1:3),
c     &  " ! Zero-weight lines used  for \eps calculation.
c     &   To make it effective, enlarge nq0i and do echo 2|qg4gw!"
        close(1101)
c ... end of get q0x
      elseif(iq0pin==2) then
C RS: not implemented!
         stop 'RS: Error - iq0pin==2 is not implemented yet.'
CAug2005
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
            if(nq0i <=0) stop 'There are neither <QforEPS> nor <QforEPS>. '
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
                  read(ifinin,*) qmin(:,i), qmax(:,i), ndiv(i)
                  write(6,"('<QforEPSL>' 3f12.8,2x,3f12.8,i5)") qmin(:,i), qmax(:,i), ndiv(i)
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
c        read (1101,*) nq0i
c        allocate( wt(nq0i),q0i(3,nq0i) )
c        read (1101,*) ( wt(i),q0i(1:3,i),i=1,nq0i)
c        rewind(1101)
c        write(1101,"(i5,a)") nq0i," ---This is readin Q0P --- "
         write(1101,"(2i5,a)") nq0i,iq0pin
     &        ," !nq0i iq0in ---"//
     &        "This is readin Q0P from GWinput <QforEPS> ---"
         write(1101,"(d24.16,3x, 3d24.16)") (wt(i),q0i(1:3,i),i=1,nq0i)
         close(1101)
      endif ! elseif(iq0pin==2) then
cccccccccccccccccccc
c      print *,' end fo writing QOP'

C... feb2006 enlarge number of q
      if(.not.timereversal()) then
        write(6,*)" timereversal==off : take -Q0P"
        do iq=1,nq0i
          q0i(:,iq+nq0i)= -q0i(:,iq)
        enddo 
        nq0i=nq0i*2
      endif

c--- Get qq; q-point for eigenfunctions'
      regmeshg = qbzreg()        !Gamma mesh based on regular mesh
      offmeshg = .not.qbzreg()   !Gamma mesh based on off-regular mesh
C
C RS: regmeshg==F is not meaningful in band-plot mode..(Maybe,I'm not sure)
      if (.not.regmeshg) then
         stop "RS: qbzreg()==F is not implemented!"
      endif

      if(iq0pin==2) then !this is just for dielectric case
C RS: not implemented!
         stop "RS: iq0pin==2"
         regmesh = qbzreg()
         offmesh = (.not.qbzreg())
      else
        regmesh = .true.
c RS:
c        offmesh = (bzcase()==2) .or. (.not.qbzreg())
        offmesh = (bzcase()==2)
      endif  
      print *,' regmesh offmeshg=', regmesh,regmeshg
      print *,' offmesh offmeshg=', offmesh,offmeshg

      nqnum = 0
      if(regmesh) nqnum = nqnum + nqbz
      if(offmesh) nqnum = nqnum + nqbz
      nqnum = nqnum + nqbz*nq0i
      nqnum = nqnum + 1 !add Gamma
      nqnum = nqnum + nqnum_rsband ! for q-points along symmetry-lines 

      allocate( qq(1:3,nqnum) )

!!!      nqnum = ix !=n1q*n2q*n3q*(nq0i+1) ! total q number
C... Regular mesh (qq are the same in the both bzcase)
      ddq = 0d0
      if(bzcase()==2) ddq=  dq_bzcase2
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
C... Off Regular mesh.  
      if(offmesh) then
        do iq = 1, nqbz
          ix = ix+1
          qq(1:3,ix) = qbzrr(1:3,iq) - dq_bzcase2
        enddo
      endif
      nnn   =  ix !n1q*n2q*n3q!      if(offmesh) nnn = 2*n1q*n2q*n3q
      print *,' nnn=',nnn !This is the number to calcualte Vxc
C... Shifted mesh
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
C Add gamma
      ix=ix+1
      qq(1:3,ix)=0d0

ccccc
c RS: set q      
      do iq=1,nqnum_rsband
        ix = ix + 1
        qq(1:3,ix) = qq_rsband(1:3,iq)
      enddo

C RS:
      print *,' nqbnu_rsband=',nqnum_rsband

C... Get qqm; q point for eigenvalues. Saved to Qmtet
      if(sum(abs(mtet))/=3) then
c RS: not implemented and omitted
        stop 'sum(abs(mtet))/=3'
      endif

 2001 continue
cccccccccccccccccccccccccccccccccccc
c      stop ' --- test end ----'
cccccccccccccccccccccccccccccccccccc


C RS: 

      write(6,*) "RS: qreduce0 = ",qreduce0
c ... get QpG for given q; get number of plane waves, ngpn ngcn.


      if( qreduce0 ) then !June2005 remove duplicated qqm
         write(6,*) "RS: Before qsave() nqnum = ",nqnum
        nqnum0 = nqnum-nqnum_rsband
        allocate(qsave(3,nqnum))
        qsave(:,nqnum0+1:nqnum) = qq(:,nqnum0+1:nqnum)
        imx=0
c RS: qq(1:3,nqnum-nqnum_rsband:nqnum) remains the same
c        do iq=1,nqnum
        do iq=1,nqnum0
          call qqsave(qq(1:3,iq),nqnum0,ginv,qsave,imx)
        enddo
        qq(:,1:imx)=qsave(:,1:imx)
        qq(:,imx+1:imx+nqnum_rsband) = qsave(:,nqnum0+1:nqnum)
        nqnum=imx+nqnum_rsband
        deallocate(qsave)
        write(6,*) "RS: After qsave() nqnum = ",nqnum
      endif
      
      
      allocate(ngpn(nqnum), ngcn(nqnum))
      if(verbose()>150) 
     &  write(6,*) ' *** q vector in 1st BZ + q0 shift.   ngp ***'
      do iq = 1, nqnum
        q = qq(1:3,iq)
        call getgv2(alat,plat,qlat, q, QpGcut_psi,1,ngpn(iq),dummy)
        call getgv2(alat,plat,qlat, q, QpGcut_Cou,1,ngcn(iq),dummy)
        if(verbose()>150) write(6,'(3f12.5,3x,2i4)') q ,ngpn(iq) !,ngcn(iq,iq00)
      enddo
      ngpmx = maxval(ngpn)
      ngcmx = maxval(ngcn)

      write(ifiqg ) nqnum,ngpmx,QpGcut_psi, nnn !takao 2003 Sep.
      write(ifiqgc) nqnum,ngcmx,QpGcut_cou
      write(6,*) ' --- Maxmum number of G for psi =',ngpmx
      write(6,*) ' --- Maxmum number of G for Cou =',ngcmx
c
c ... Get G vector for given q.
      do iq = 1, nqnum
        q = qq(1:3,iq)
        ngp = ngpn(iq)
        ngc = ngcn(iq)
        write (6,"(' iq=',i5,' q=',3f7.3,' ngp ngc = ',2i6)")
     &  iq, q, ngp, ngc
        allocate( ngvecp(3,max(ngp,1)), ngvecc(3,max(ngc,1)) ) !larger than 1. Aug2005
        call getgv2(alat,plat,qlat, q, QpGcut_psi, 2, ngp,  ngvecp)  ! for eigenfunction
        call getgv2(alat,plat,qlat, q, QpGcut_Cou, 2, ngc,  ngvecc)  ! for Coulomb
        write (ifiqg) q, ngp 
        write (ifiqg) ngvecp
        write (ifiqgc) q, ngc
        write (ifiqgc) ngvecc
        deallocate(ngvecp,ngvecc) 
      enddo
      deallocate(ngpn,ngcn)
      if(iq0pin==1) deallocate(ngvect)
c      print *,'end of mkqg'
      end

c --- tripl
      doubleprecision function tripl(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
c     dimension a(3),b(3),c(3)
c ... g77 needs this rewriting for opt > O1
      dimension a(0:2),b(0:2),c(0:2)
c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c    .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c ... g77 needs this rewriting
      tmp = 0.d0
      do i = 0,2
          j = mod(i + 1, 3)
          k = mod(i + 2, 3)
          tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
      enddo
      tripl = tmp
      end

      subroutine genbzdata1_xxxxxxx(plat,n1q,n2q,n3q,symops,ngrp,
     &  tetraf, tetrai,
     o  qbz, wbz,              nqbz, 
     o  qibz,wibz, nstar,irk,  nqibz,
     o  idtetf, ntetf, qbzw,ib1bz, nqbzw,
     o  idteti, nteti )
C- Generate BZ mesh data reuired for GW ---
Ciput
C  plat,qlat,n1q,n2,n3q,symops,ngrp,
C  tetraf: if true, this generate the BZ data for full-1stBZ tetrahedron method
C               o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C  tetrai: if true, this generate the BZ data for tetrahedron method in irreducible BZ.
C               o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C output
C     o  qbz, wbz,              nqbz, 
C     o  qibz,wibz, nstar,irk,  nqibz,
C ouptput(optional)  
C     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
C     o  idteti, nteti 
C-------------------------------------------------------------------------
      implicit none
      integer(4):: nqbz,n1q,n2q,n3q,ntetf,nteti,nqbzw,iqbz,ifbz
     & ,nqibz,iqibz,itet,ngrp
      real(8) :: plat(3,3),qbas(3,3),ginv(3,3),vol,symops(3,3,ngrp)
     &           ,qbasmc(3,3)
      logical tetraf,tetrai
      integer(4),allocatable:: ipq(:),iw1(:)
      real(8):: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:)
      integer(4) :: idtetf(:,:),ib1bz(:),idteti(:,:),irk(:,:)
     &       ,nstar(:)
      end
