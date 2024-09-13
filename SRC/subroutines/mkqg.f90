subroutine mkQG2(iq0pin,gammacellctrl,lnq0iadd,lmagnon)! Make required q and G to expand eigenfunctions for GW.
  use m_ftox
  use m_get_bzdata1,only: Getbzdata1, nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm
  use m_get_bzdata1,only: qbz,wbz,qibz,wibz, qbzw, idtetf, ib1bz, idteti, irk, nstar, nstbz 
  use m_q0p,only: Getallq0p, nq0i,nq0itrue,nq0iadd,  q0i,wt, epslgroup, lxklm, epinv,wklm, dmlx, epinvq0i,ixyz
  use m_keyvalue,only: getkeyvalue
  use m_hamindex0,only: Readhamindex0, symops,ngrp,alat,plat,qlat
  use m_lgunit,only: stdo
  implicit none
  intent(in)::   iq0pin, gammacellctrl,lnq0iadd,lmagnon
  !!     |q+G| < QpGcut_psi for eigenfunction psi.
  !!     |q+G| < QpGcut_Cou for coulomb interaction
  !!
  !! OUTPUT
  !!     file handle= ifiqg,  which contains q and G points for eigenfunction psi. --> QGpsi
  !!     file handle= ifiqgc, which contains q and G points for Coulomb            --> QGcou
  !!     QGpsi(ifiqg), QGcou(ifiqgc), Q0P are written.
  !!     BZDATA,QBZ.chk
  !!        The console output gives a list of q points for human check.
  !! ---------------------------------------------------
  integer ::nnn(3),ifiqg,ifiqgc,ngcxx,i,j,iq,iq00,ngp,ngpmx,ngc,ngcmx,nqnum,iq0pin,dummyia(1,1),iimx,irradd,nmax, &
       nline,nlinemax,ifsyml,iqq,is,nk,ix,nqnumx,i1,ifkpt, mxkp,ifiqibz,iqibz,ifigwin,mtet(3),nm1,nm2,nm3,&
       gammacellctrl,nnng(3),ifi0,itet,ifi00,ifidmlx,ierr,retval, nqi,ifix,ig,iq0i,lm,ifi, nqnumm,ifiqmtet,verbose,nn1,nn2,&
       ifiqbz,iqbz, ifbz, imxc,nnn3(3),imx0c,imx11(1,1), ifidml, llxxx,lm1,lm2, imx,ifinin,il,imx0
  integer,allocatable :: ngvecprev(:,:,:),ngveccrev(:,:,:),ngvecp(:,:), ngvecc(:,:), ngpn(:),ngcn(:), nqq(:),irr(:)
  real(8):: dq_(3),qlatbz(3,3),imat33(3,3),unit,QpGx2,aaa,q(3),dummy,qp(3), QpGcut_psi, QpGcut_Cou,QpGcut,alpv(3),q0smean,sumt,&
       alp, volum,voltot,q0(3),qlat0(3,3),tripl, xx,qqx(3),alpm,aaij,bbij,vol,ginv(3,3),dq(3),ddq(3)
  real(8):: qx(3),qxx(3),tolq, deltaq,delta5,delta8,deltaq_scale,volinv,wtrue00,qg(3),alpqg2,qg2,tpiba
  real(8),parameter:: pi=4d0* atan(1d0)
  real(8),allocatable :: qq(:,:),qq1(:,:),qq2(:,:),qqm(:,:),qsave(:,:), wt0(:),wti(:),qi(:,:)
  real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:)
  logical ::tetrai,tetraf,tetra_hsfp0,qbzreg, qreduce ,qreduce0
  logical :: offmesh=.false. ,offmeshg=.false.
  logical :: regmesh=.false. ,regmeshg=.false. ,  timereversal
  logical :: caca,debug=.false. !,newaniso
  logical :: newoffsetG !july2014
  logical :: lnq0iadd, lmagnon,unit2=.false. ,cmdopt0
  logical :: keepqg
  integer :: ifiqg2,ifiqgc2
  integer, allocatable :: ngvecp_tmp(:,:),ngvecc_tmp(:,:)
  character*99:: q0pf        !nov2012
  write(stdo,"('mkqg2: ')")
  call readhamindex0()
  call getkeyvalue("GWinput", "n1n2n3", nnn,3)
  if(lmagnon) call getkeyvalue("GWinput", "n1n2n3eps",nnn,3,default=nnn)
  call getkeyvalue("GWinput", "QpGcut_psi",QpGx2)
  call getkeyvalue("GWinput", "QpGcut_cou",QpGcut_Cou)
  call getkeyvalue("GWinput", "unit_2pioa",unit2)
  call getkeyvalue("GWinput","KeepQG",keepqg,default=.true.)
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
  ginv = transpose(plat)
  !! NOTE: we use only mtet=[1,1,1]. If we like to recover general case, examine code again.
  call getkeyvalue("GWinput","multitet",mtet,3,default=[1,1,1])
  write(stdo,"('  SYMOPS ngrp=',i3)") ngrp
  write(stdo,"('  unit(a.u.) alat  =',f13.6 )") alat
  write(stdo,"('  ---  |k+G| < QpG(psi) QpG(Cou)=',2d13.6)") QpGcut_psi, QpGcut_Cou
  write(stdo,"('  --- k points for GW from GWinput =',3i3)") nnn(1:3)
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
     ! his shift vector is to make the Gamma point centered in the Gamma cell.
     tetrai=.false.
     call minv33(qlatbz,ginv)
     write(stdo,*)'=== Gammacell qlatgz ==='
     write(stdo,"(3d23.15)") qlatbz
     write(stdo,*)'=== Gammacell ginv ==='
     write(stdo,"(3f9.4)") ginv
  else
     qlatbz(:,:) = qlat(:,:)
     tetrai = .true.         !used in heftet tetra_hsfp0()
     dq_ = 0d0
     if( .NOT. qbzreg()) dq_ = -matmul(qlat(1:3,1:3),(/.5d0/nnn(1),.5d0/nnn(2),.5d0/nnn(3)/))
     ! dq_ is off-gamma mesh, used when qbzreg=F
  endif
  if(sum(abs(dq_))>tolq()) write(stdo,'(" Shift vector (skip Gamma) by dq_=",3f9.4)')dq_
  !! Get BZ data by 'call Getbzdata1'. See m_get_bzdata1
  !! we usually use negative delta (tetrahedron).
  call getkeyvalue("GWinput","delta",aaa)
  if(aaa<0d0) then
     write(stdo,"(a)")'  GWinput delta<0: tetrahedron method for x0'
     tetraf=.true.
  else
     write(stdo,"(a)")'  GWinput delta>0: not use tetrahedron method for x0'
     tetraf=.false.
  endif
  call Getbzdata1(qlatbz,nnn,symops,ngrp,tetrai,tetraf,mtet,gammacellctrl) !all inputs
  open (newunit=ifiqbz, file='QBZ.chk') !write q-points in IBZ.
  write(ifiqbz,"(i10)") nqbz
  do iqbz = 1,nqbz
     write(ifiqbz,"(3d24.16,3x,d24.16)") qbz(1:3,iqbz)
  enddo
  close(ifiqbz)
  write(stdo,"('  --- TOTAL num of q nqbz and nqibz=)',2i6)") nqbz,nqibz
  write(stdo,'("  qibz = ",i6,3f12.5)')(i,qibz(1:3,i),i=1,min(10,nqibz))
  write(stdo,*)" ... QIBZ is written in QIBZ file ..."
  write(stdo,*)
  call getkeyvalue("GWinput","alpha_OffG",alp,default=-1d60)
  alpv(:)=alp
  if(alp==-1d60) then
     call getkeyvalue("GWinput","alpha_OffG_vec",alpv,3,default=(/-1d50,0d0,0d0/))
     if(alpv(1)==-1d50) then
        call rx( ' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput')
     endif
  endif
  call Getallq0p(iq0pin,newoffsetG,alat,plat,qlat,nnn,alp,alpv, &
       nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,lnq0iadd)
  do i=nq0i+1,nq0i+nq0iadd
     write(stdo,"('  q0iadd=  ', i3, 3f10.5)") i,q0i(:,i)
  enddo
  print *,' Writing BZDATA...'
  open(newunit=ifbz, file='BZDATA',form='unformatted')
  write(ifbz) nqbz,nqibz, nqbzw, ntetf, nteti,ngrp,nnn ,qlat,ginv
  write(ifbz) qibz(1:3,1:nqibz),wibz(1:nqibz),nstar(1:nqibz),irk(1:nqibz,1:ngrp)
  write(ifbz) qbz(1:3,1:nqbz),wbz(1:nqbz),nstbz(1:nqbz)
  write(ifbz) idtetf(0:3,1:ntetf)
  write(ifbz) ib1bz(1:nqbzw), qbzw(1:3,1:nqbzw)
  write(ifbz) idteti(0:4,1:nteti)
  write(ifbz) dq_
  write(ifbz) iq0pin,nq0i,iq0pin,nq0iadd
  write(ifbz) wt(1:nq0i+nq0iadd),q0i(1:3,1:nq0i+nq0iadd)
  if(iq0pin==1) then
     write(ifbz) ixyz(1:nq0i+nq0iadd)
     write(ifbz) lxklm,dmlx(1:nq0i,1:9),epinv(1:3,1:3,1:nq0i),epinvq0i(1:nq0i,1:nq0i)
     write(ifbz) wklm(1:(lxklm+1)**2)
  elseif(iq0pin==2.or.iq0pin==3) then
     write(ifbz) epslgroup(1:nq0i)
  endif
  close(ifbz)
  !! ---------
  if (lmagnon) return       !! exit because we get QforEPSL data only (21 Feb, 2020)
  !! Four kinds of mesh points setting. Q0P means offset Gamma (slightly different from Gamma).
  !! Which we need?
  !! 1. regular
  !! 2. offregular (not including Gamma)
  !! 3. regular    + Q0P
  !! 4. offregular + Q0P
  if(iq0pin==2.or.iq0pin==3) then        !this is just for dielectric case
     regmesh = qbzreg()
  else
     regmesh = .true.
  endif
  regmeshg = qbzreg()       !Gamma mesh based on regular mesh
  offmesh =  .not.qbzreg()  !we fix bzcase=1 now. apr2015.
  offmeshg = .not.qbzreg()  !Gamma mesh based on off-regular mesh
  write(stdo,*)
  write(stdo,"('  regmesh offmeshg=',2l)") regmesh,regmeshg !regular,     regular+shifted
  write(stdo,"('  offmesh offmeshg=',2l)") offmesh,offmeshg !offregmesh, offregular+shifted
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
2112 continue
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
  AddoffsetGammaandGammapoint:do iq00 = 1, nq0i+ nq0iadd !duplicated are removed by qreduce 
     ix = ix+1
     qq(1:3,ix) = q0i(1:3,iq00)
  enddo AddoffsetGammaandGammapoint
  ix=ix+1
  qq(1:3,ix)=0d0
  RemoveEquilvalentqBYTranslatioanlSymmetry: if( qreduce0 ) then
     call cputid (0)
     nmax= nq0i+nq0iadd+nqnum
     write(stdo,"(a,4i8)")'  nq0i nq0iadd nqnum=',nq0i,nq0iadd,nqnum,nmax
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
  endif RemoveEquilvalentqBYTranslatioanlSymmetry
2001 continue
  !! Here we get all requied q points. We do reduce them by space group symmetry.
  if(allocated(wt0)) deallocate(wt0)
  allocate(wt0(nqnum+nq0i+ nq0iadd ),qi(3,nqnum+nq0i+ nq0iadd ),wti(nqnum+nq0i+ nq0iadd ))
  wt0=1d0 
  !write(stdo,*)'ppppppppp',nqnum,nq0i,nq0iadd
  !! Set irreducible k-point flag. irr=1 for (irredusible point) flag, otherwise =0.
  !! irr(iq)=1 for irreducile qq(:,iq), iq=1,nqnum
  call q0irre(qibz,nqibz,qq,wt0,nqnum,symops,ngrp, qi,nqi,wti,plat,.true.,0,irr,nqbz)
  if(cmdopt0('--allqbz')) nqnum=nqbz
  !! nqnum is the finally obtained number of q points.
  allocate(ngpn(nqnum), ngcn(nqnum))
  if(debug) write(stdo,*) ' --- q vector in 1st BZ + Q0P shift. ngp ---'
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
     call getgv2(alat,plat,qlat, qxx, QpGcut_psi,1,ngpn(iq),imx11) ! get nqpn. # of G vector for |q+G| < QpGcut_psi
     imx0=imx11(1,1)
     if(imx0>imx) imx=imx0
     ngcn(iq)=1
     call getgv2(alat,plat,qlat, qxx, QpGcut_Cou,1,ngcn(iq),imx11) ! get ngcn. # ofG vector for |q+G| < QpGcut_cou
     imx0c=imx11(1,1) 
     if(imx0c>imxc) imxc=imx0c
     if(verbose()>150)write(stdo,'(3f12.5,3x,2i4)') q ,ngpn(iq) !,ngcn(iq,iq00)
     if(verbose()>150)write(stdo,'(3f12.5,3x,2i4)') q ,ngcn(iq) !,ngcn(iq,iq00)
  enddo
  ! Get G vectors and Write q+G vectors -----------
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
  write(stdo,*) ' --- Max number of G for psi =',ngpmx
  write(stdo,*) ' --- Max number of G for Cou =',ngcmx
  allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) )       !inverse mapping table for ngvecp (psi)
  allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc) ) !inverse mapping table for ngvecc (cou)
  if(.not.keepqg) allocate(ngvecp_tmp(3,ngpmx),ngvecc_tmp(3,ngcmx))
  if(.not.keepqg) then
    open(newunit=ifiqg2 ,file='QGpsi_rec',form='unformatted', access='direct', recl=4*(3*ngpmx+(2*imx+1)**3))
    open(newunit=ifiqgc2,file='QGcou_rec',form='unformatted', access='direct', recl=4*(3*ngcmx+(2*imxc+1)**3))
  endif
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
     write(stdo,"(' iq=',i8,' q=',3f9.5,' ngp ngc= ',2i6,' irr.=',i2,a)")iq, q, ngp, ngc, irr(iq),trim(q0pf) !irr=1 is irr. k points.
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
     if(.not.keepqg) then
        ngvecp_tmp(1:3,1:ngp) = ngvecp(1:3,1:ngp)
        ngvecc_tmp(1:3,1:ngc) = ngvecc(1:3,1:ngc)
        write(ifiqg2,rec=iq) ngvecp_tmp,ngvecprev
        write(ifiqgc2,rec=iq) ngvecc_tmp,ngveccrev
      endif
     deallocate(ngvecp,ngvecc)
  enddo
  if(.not.keepqg) then
    deallocate(ngvecp_tmp,ngvecc_tmp)
    close(ifiqg2)
    close(ifiqgc2)
  endif
  deallocate(ngpn,ngcn,ngvecprev,ngveccrev)
  if(debug) print *,'--- end of mkqg2 ---'
end subroutine mkQG2
