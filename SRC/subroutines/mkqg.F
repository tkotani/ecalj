      subroutine mkQG2(iq0pin, gammacellctrl,lnq0iadd,lmagnon)
!! 'call getbzdata1' gives all follwing data 
      use m_get_bzdata1,only: Getbzdata1, nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm,
     &  qbz,wbz,qibz,wibz, qbzw, idtetf, ib1bz, idteti, irk, nstar, nstbz, qbzm, qbzwm
      use m_q0p,only: Getallq0p, q0i,wt,nq0i,nq0itrue, nq0iadd
      use m_keyvalue,only: getkeyvalue
      use m_hamindex0,only: Readhamindex0, symops,ngrp,alat,plat,qlat
      implicit none
      intent(in)::     iq0pin, gammacellctrl,lnq0iadd,lmagnon
!! == Make required q and G in the expansion of GW. ==
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
      integer ::nnn(3),ifiqg,ifiqgc,ngcxx,
     &     i,j,iq,iq00,ngp,ngpmx,ngc,ngcmx,nqnum,iq0pin,
     &     nline,nlinemax,ifsyml,iqq,is,nk,ix,nqnumx,i1,ifkpt
      integer :: mxkp,ifiqibz,iqibz,ifigwin,mtet(3),nm1,nm2,nm3
      integer:: nqnumm,ifiqmtet,verbose,nn1,nn2,ifiqbz,iqbz !,auxfunq0p
      real(8)  :: q(3),dummy,qp(3), 
     &     QpGcut_psi, QpGcut_Cou,QpGcut,alpv(3),q0smean,sumt,alp,
     &     volum,voltot,q0(3),qlat0(3,3),tripl,
     &     xx,qqx(3),alpm
      real(8)::aaij,bbij
      real(8) :: vol,ginv(3,3),dq(3) !,www
      integer,allocatable:: ngvecp(:,:), ngvecc(:,:),
     &     ngpn(:),ngcn(:), nqq(:)   
      real(8),allocatable ::
     &     qq(:,:),qq1(:,:),qq2(:,:),qqm(:,:)
      logical ::tetrai,tetraf,tetra_hsfp0    
      integer :: ifbz
      logical:: qbzreg, qreduce ,qreduce0
      real(8),allocatable:: qsave(:,:)
      integer:: imx,ifinin,il,imx0
      integer,allocatable :: ngvecprev(:,:,:),ngveccrev(:,:,:)
c      integer(4):: bzcase=1
c     logical :: readgwinput
      real(8):: ddq(3)
      logical :: offmesh=.false. ,offmeshg=.false.
      logical :: regmesh=.false. ,regmeshg=.false. ,  timereversal
      logical :: caca,debug=.false. !,newaniso
      integer:: imxc,nnn3(3),imx0c,imx11(1,1)
      real(8):: deltaq,delta5,delta8,deltaq_scale!=1d0/3.0**.5d0 
      integer:: nqi,ifix,ig,iq0i,lm,ifi
      real(8),allocatable:: wti(:),qi(:,:) !,symops(:,:,:)
      integer:: ifidml!,iclose,iopen !,ifiwqfac
      integer:: llxxx,lm1,lm2
      real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:)
      real(8)::volinv,wtrue00,qg(3),alpqg2,qg2,tpiba
      character*99:: q0pf        !nov2012
      integer:: dummyia(1,1),iimx,irradd,nmax
      real(8):: qx(3),qxx(3),tolq
      logical :: newoffsetG !july2014
      real(8),allocatable:: wt0(:)
      integer,allocatable::irr(:)
      real(8):: dq_(3),qlatbz(3,3)
      integer:: gammacellctrl,nnng(3),ifi0,itet !, ngrp
      real(8)::imat33(3,3),unit,QpGx2,aaa
      logical:: lnq0iadd, lmagnon
      logical ::unit2=.false. 
      real(8),parameter:: pi=4d0* atan(1d0)  
c------------------------------------------------
      write(6,"('mkqg2: ')")
! initial set up      
c      open (newunit=ifi, file='LATTC',form='unformatted')
c      read(ifi) alat,plat
c      close(ifi)
c      
c      open (newunit=ifi, file='SYMOPS')
c      read(ifi,*) ngrp
c      allocate(symops(3,3,ngrp))
c      do ig = 1,ngrp
c        read(ifi,*)
c        do i=1,3
c          read(ifi,*) symops(i,1:3,ig)
c        enddo
c      enddo
c      close(ifi)
      call readhamindex0()
      
      call getkeyvalue("GWinput", "n1n2n3", nnn,3)
      call getkeyvalue("GWinput", "QpGcut_psi",QpGx2)
      call getkeyvalue("GWinput", "QpGcut_cou",QpGcut_Cou)
      call getkeyvalue("GWinput", "unit_2pioa",unit2)
      if(unit2) then
        unit = 2d0*pi/alat
        QpGx2     = QpGx2      *unit
        QpGcut_cou= QpGcut_cou *unit
      endif
      QpGcut_psi = QpGx2
      if(iq0pin==4) then
         QpGcut_psi=0d0
         QpGcut_Cou=0d0
      endif
     
      qreduce0 = qreduce()
      newoffsetG=.true. 
c      call minv33tp(plat,qlat)
      ginv = transpose(plat)
      
!! NOTE: we use only mtet=[1,1,1]. If we like to recover general case, examine code again.
      call getkeyvalue("GWinput","multitet",mtet,3,default=[1,1,1])
      write(6,"('  SYMOPS ngrp=',i3)") ngrp
      write(6,"('  unit(a.u.) alat  =',f13.6 )") alat
      write(6,"('  ---  |k+G| < QpG(psi) QpG(Cou)=',2d13.6)") QpGcut_psi, QpGcut_Cou
      write(6,"('  --- k points for GW from GWinput =',3i3)") nnn(1:3) 
c debug      
c      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
c      imat33=0d0
c      imat33(1,1)=1d0
c      imat33(2,2)=1d0
c      imat33(3,3)=1d0
c      if(sum(abs(matmul(transpose(qlat),plat)-imat33))>tolq()) call rx('qlat*plat err')
c      if(sum(abs(ginv-transpose(plat)))       >tolq()) call rx('ginv=transpose(qlat) err')
      
      open(newunit=ifiqg ,file='QGpsi',form='unformatted')
      open(newunit=ifiqgc,file='QGcou',form='unformatted')

!     ! gammacellctrl=2 is a special mode. We only consider tetrahedron method within the Gammacell.
!! The Gammacell is a part of BZ made from three vectors following qlatbz=(qlat(:,1)/n1q,...)
!! Then the Gamma point is in the middle of micro_qlat = (qlat(:,1)/n1q,qlat(:,2)/n2q,...)
!! To get qbz which is in the Gamma cell, we use shift in the getbzdata1 for gammacellctrl=2.
!! Tetrahedron method is applied for such qbz.
      if(gammacellctrl==2) then
        do i=1,3
          qlatbz(:,i) = qlat(:,i)/nnn(i) !qlat for Gamma cell
        enddo
        call getkeyvalue("GWinput","GammaDivn1n2n3",nnng,3)  
        nnn = nnng          !division of Gamma cell
        dq_ = -matmul(qlatbz(1:3,1:3),(/.5d0,.5d0,.5d0/)) 
!This shift vector is to make the Gamma point centered in the Gamma cell.
        tetrai=.false.
        call minv33(qlatbz,ginv)
        write(6,*)'=== Gammacell qlatgz ==='
        write(6,"(3d23.15)") qlatbz
        write(6,*)'=== Gammacell ginv ==='
        write(6,"(3f9.4)") ginv
      else
        qlatbz(:,:) = qlat(:,:)
        tetrai = .true.         !used in heftet tetra_hsfp0()
        dq_ = 0d0
        if(.not.qbzreg()) dq_ = -matmul(qlat(1:3,1:3),(/.5d0/nnn(1),.5d0/nnn(2),.5d0/nnn(3)/)) 
                        ! dq_ is off-gamma mesh, used when qbzreg=F
      endif
      if(sum(abs(dq_))>tolq()) write(6,'(" Shift vector (skip Gamma) by dq_=",3f9.4)')dq_
      
!! Get BZ data by 'call Getbzdata1'. See m_get_bzdata1 
!! we usually use negative delta (tetrahedron).
      call getkeyvalue("GWinput","delta",aaa)
      if(aaa<0d0) then
         write(6,"(a)")'  GWinput delta<0: tetrahedron method for x0'
         tetraf=.true.
      else
         write(6,"(a)")'  GWinput delta>0: not use tetrahedron method for x0'
         tetraf=.false.
      endif
      call Getbzdata1(qlatbz,nnn,symops,ngrp,tetrai,tetraf,mtet,gammacellctrl) !all inputs
!! Write BZDATA
      print *,' Writing BZDATA...'
      open (newunit=ifbz, file='BZDATA',form='unformatted')
      write(ifbz) nqbz,nqibz, nqbzw, ntetf, nteti,ngrp,nnn ,qlat,ginv
      write(ifbz) qibz(1:3,1:nqibz),wibz(1:nqibz),nstar(1:nqibz),irk(1:nqibz,1:ngrp)
      write(ifbz) qbz(1:3,1:nqbz),wbz(1:nqbz),nstbz(1:nqbz)
      write(ifbz) idtetf(0:3,1:ntetf)
      write(ifbz) ib1bz(1:nqbzw), qbzw(1:3,1:nqbzw)
      write(ifbz) idteti(0:4,1:nteti)
      write(ifbz) dq_
      close(ifbz)
c!! Write QIBZ
c      write(6,*)' qibz are written in QIBZ file...'
c      open (newunit=ifiqibz, file='QIBZ',form='unformatted') !write q-points in IBZ.
c      write(ifiqibz) nqibz
c      write(ifiqibz) qibz(1:3,1:nqibz),wibz(1:nqibz)
c      close(ifiqibz)
!! Write QBZ
      open (newunit=ifiqbz, file='QBZ.chk') !write q-points in IBZ.
      write(ifiqbz,"(i10)") nqbz
      do iqbz = 1,nqbz
         write(ifiqbz,"(3d24.16,3x,d24.16)") qbz(1:3,iqbz)
      enddo
      close(ifiqbz)
!!  Write KPNTin1BZ.mkqg.chk (files *.chk is only for check.).
      open(newunit=ifkpt,file='KPTin1BZ.mkqg.chk')
      write(ifkpt,*)"  qbz --> shoten(qbz)"
      do i1 = 1,nqbz
         call shorbz(qbz(1,i1),qp,qlat,plat)
         write (ifkpt,"(1x,i7,4f10.5,'   ',3f10.5)") 
     &        i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
      enddo
      close (ifkpt)
      write(6,"('  --- TOTAL num of q nqbz and nqibz=)',2i6)") nqbz,nqibz
      write(6,'("  qibz = ",i6,3f12.5)')(i,qibz(1:3,i),i=1,min(10,nqibz))
      write(6,*)" ... QIBZ is written in QIBZ file ..."
      write(6,*)
      
!! Getallq0p; Q0P is offset Gamma or k point given in GWinput
!! See use m_q0p =>  q0i,wt,nq0i,nq0itrue are outputs
!! After 'call Getallq0p', we have q0i(:,nq0i+1,nq0i+nq0iadd). 
!!    q0i(:,1:nq0i+n0qiadd) contains all q0x(:,i)= qlat(:,i)/nnn(i)/2d0*deltaq_scale() for i=1,3.
! alpha is for auxially function of the offset Gamma method.
      call getkeyvalue("GWinput","alpha_OffG",alp,default=-1d60)
      alpv(:)=alp
      if(alp==-1d60) then
       call getkeyvalue("GWinput","alpha_OffG_vec",alpv,3,default=(/-1d50,0d0,0d0/))
       if(alpv(1)==-1d50) then
         call rx( ' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput')
       endif
      endif
      call Getallq0p(iq0pin,newoffsetG,alat,plat,qlat,nnn,alp,alpv,
     i     nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,lnq0iadd)
      do i=nq0i+1,nq0i+nq0iadd
        write(6,"('  q0iadd=  ', i3, 3f10.5)") i,q0i(:,i)
      enddo
      if (lmagnon) return !! exit because we get QforEPSL data (21 Feb, 2020)
      
!! Four kinds of mesh points setting. Q0P means offset Gamma (slightly different from Gamma).
!! Which we need?
!! 1. regular 
!! 2. offregular (not including Gamma)
!! 3. regular    + Q0P
!! 4. offregular + Q0P
      if(iq0pin==2) then        !this is just for dielectric case
         regmesh = qbzreg()
      else
         regmesh = .true.
      endif
      regmeshg = qbzreg()       !Gamma mesh based on regular mesh
      offmesh =  .not.qbzreg()  !we fix bzcase=1 now. apr2015.
      offmeshg = .not.qbzreg()  !Gamma mesh based on off-regular mesh
      write(6,*)
      write(6,"('  regmesh offmeshg=',2l)") regmesh,regmeshg !regular,     regular+shifted
      write(6,"('  offmesh offmeshg=',2l)") offmesh,offmeshg !offregmesh, offregular+shifted
!!  We check wether all q0i \in qbz or not. <--- Takao think this block is not necessary now.
      call minv33(qlat,ginv)
      nqnum = nqbz
      allocate( qq(1:3,nqnum),irr(nqnum) )
      qq(1:3,1:nqbz) = qbz(1:3,1:nqbz) 
      do iq0i=1,nq0i+nq0iadd
         do iq=1,nqbz
            if(sum(abs(q0i(:,iq0i)-qq(:,iq)))<tolq()) goto 2112
            call rangedq( matmul(ginv,q0i(:,iq0i)-qq(:,iq)), qx)
            if(sum(abs(qx))< tolq()) goto 2112
         enddo
         goto 2111
 2112    continue
         qq(:,iq) = q0i(:,iq0i) !replaced with equivalent q0i.
      enddo
      print *,' --- We find all q0i in qbz. Skip qreduce.'
      goto 2001
 2111 continue
      
!! Accumulate all required q points
      deallocate(qq,irr)
      nqnum = nqbz  + nqbz*(nq0i+nq0iadd)
      nqnum = nqnum + 1         !add Gamma
      nqnum = nqnum + nq0i + nq0iadd      !add Gamma + q0i
      allocate( qq(1:3,nqnum),irr(nqnum) )
      ix = 0
      if(regmesh) then
         qq(1:3,1:nqbz) = qbz(1:3,1:nqbz)
         ix = ix+ nqbz
      endif
!! Off Regular mesh.
      if(offmesh) then
         do iq = 1, nqbz
            ix = ix+1
            qq(1:3,ix) = qbz(1:3,iq) - dq_
         enddo
      endif
!! Shifted mesh
      if(regmeshg) then
         do iq00 = 1, nq0i+ nq0iadd 
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbz(1:3,iq) +  q0i(1:3,iq00)
            enddo
         enddo
      endif
      if(offmeshg) then
         do iq00 = 1, nq0i+ nq0iadd 
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbz(1:3,iq) - dq_ + q0i(1:3,iq00)
            enddo
         enddo
      endif
!!  - Add offset Gamma and Gamma point (these can be removed by qreduce and q0irre)
      do iq00 = 1, nq0i+ nq0iadd 
         ix = ix+1
         qq(1:3,ix) = q0i(1:3,iq00)
      enddo
      ix=ix+1
      qq(1:3,ix)=0d0
c$$$!! (this mtet block is not used now) Get qqm; q point for eigenvalues. 
c$$$!! Saved to Qmtet. Not so much used now... 
c$$$!! We need check again if we like to use this branch again (2016apr)
c$$$      if(sum(abs(mtet))/=3) then
c$$$         nqnumm= nqbzm * (nq0i+ nq0iadd +1)
c$$$         allocate( qqm(1:3,nqnumm) )
c$$$         ix=0
c$$$         do iq00 = 1, 1 + nq0i+ nq0iadd 
c$$$            do iq   = 1, nqbzm
c$$$               ix = ix+1
c$$$               if(iq00==1) then
c$$$                  qqm(1:3,ix) = qbzm(1:3,iq)
c$$$               else
c$$$                  qqm(1:3,ix) = q0i(1:3,iq00-1) + qbzm(1:3,iq)
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$c         ifiqmtet=ifile_handle()
c$$$         open(newunit=ifiqmtet, file='Qmtet')
c$$$         write(ifiqmtet,"(i10)") nqnumm
c$$$         do iq=1,nqnumm
c$$$            write(ifiqmtet,"(3d24.16)") qqm(1:3,iq)
c$$$         enddo
c$$$         close(ifiqmtet)
c$$$         deallocate(qqm)
c$$$  endif
      
!! Remove equivalent q point by the translational symmetry 
      if( qreduce0 ) then
         write(6,"(a,3i8)")'  nq0i nq0iadd nqnum=',nq0i,nq0iadd,nqnum
         call cputid (0)
         nmax= nq0i+nq0iadd+nqnum
         allocate(qsave(3,nmax)) !,qsavel(nmax))
         imx=0
         if(iq0pin /=1) then
            do iq=1,nq0i+ nq0iadd 
               call qqsave(q0i(1:3,iq),nmax,ginv,qsave,imx)
            enddo
         endif
         do iq=1,nqnum
            call qqsave(qq(1:3,iq),nmax,ginv,qsave,imx)
         enddo
         nqnum = imx
         qq(:,1:imx)=qsave(:,1:imx)
         deallocate(qsave)
      endif
!! ------------------------------------------
 2001 continue
!! ------------------------------------------

!! Here we get all requied q points. We do reduce them by space group symmetry.
      if(allocated(wt0)) deallocate(wt0)
      allocate(wt0(nqnum+nq0i+ nq0iadd ),qi(3,nqnum+nq0i+ nq0iadd ),wti(nqnum+nq0i+ nq0iadd ))
      wt0=1d0
!! Set irreducible k-point flag. irr=1 for (irredusible point) flag, otherwise =0. 
!! irr(iq)=1 for irreducile qq(:,iq), iq=1,nqnum
      call q0irre(qibz,nqibz,qq,wt0,nqnum,symops,ngrp, qi,nqi,wti,plat,.true.,0,irr)
!! nqnum is the finally obtained number of q points.
      allocate(ngpn(nqnum), ngcn(nqnum))
      if(debug) write(6,*) ' --- q vector in 1st BZ + Q0P shift. ngp ---'
      imx=0
      imxc=0
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q              
         if(iq0pin==1) then !use qxx on regular mesh points if q is on regular+Q0P(true).
            do iqbz=1,nqbz
            do i=1,nq0itrue+ nq0iadd  ! nq0itrue/=nq0i for anyq=F nov2015
               if(sum(abs(qbz(1:3,iqbz)-dq_+ q0i(:,i)-qxx))<tolq()) then
                  qxx=qbz(1:3,iqbz)
                  exit
               endif
            enddo
            enddo
         endif
         ngpn(iq)=1
!! get nqpn. # of G vector for |q+G| < QpGcut_psi
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi,1,ngpn(iq),imx11) !imx11 !nov2015
         imx0=imx11(1,1)
         if(imx0>imx) imx=imx0
         ngcn(iq)=1
!! get ngcn. # ofG vector for |q+G| < QpGcut_cou
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou,1,ngcn(iq),imx11) !imx11 to avoid warning.
         imx0c=imx11(1,1)
         if(imx0c>imxc) imxc=imx0c
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngpn(iq) !,ngcn(iq,iq00)
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngcn(iq) !,ngcn(iq,iq00)
      enddo
!! Get G vectors and Write q+G vectors -----------
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
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) )       !inverse mapping table for ngvecp (psi)
      allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc) ) !inverse mapping table for ngvecc (cou)
      ngvecprev=9999
      ngveccrev=9999
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q
         q0pf=''
         do iqbz=1,nqbz  !use qxx on regular mesh points if q is on regular+Q0P(true).
         do i=1,nq0itrue+ nq0iadd  !nq0itrue/=nq0i for anyq=F nov2015
            if(sum(abs(qbz(1:3,iqbz)-dq_+ q0i(:,i)-qxx))<tolq()) then
               if(sum(abs(q0i(:,i)-qxx))<tolq()) then
                 q0pf=' <--Q0P  '   ! offset Gamma points
               else
                 q0pf=' <--Q0P+R'   ! offset Gamma points-shifted nov2015
               endif
               if(iq0pin==1) then
                  qxx=qbz(1:3,iqbz)
               endif
               exit
            endif
         enddo
         enddo
         ngp = ngpn(iq)
         ngc = ngcn(iq)
         write(6,"(' iq=',i8,' q=',3f9.5,' ngp ngc= ',2i6,' irr.=',i2,a)")!irr=1 is irr. k points.
     &        iq, q, ngp, ngc, irr(iq),trim(q0pf)
         allocate( ngvecp(3,max(ngp,1)), ngvecc(3,max(ngc,1)) ) 
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi, 2, ngp,  ngvecp) ! for eigenfunctions (psi)
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou, 2, ngc,  ngvecc) ! for Coulomb        (cou)
         write (ifiqg) q, ngp, irr(iq)
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
      if(debug) print *,'--- end of mkqg2 ---'
      end
