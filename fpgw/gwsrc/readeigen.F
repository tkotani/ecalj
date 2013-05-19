!!--------------------------------
      module m_readqgcou
!! this is somehow duplicated with reqdqg.F---> I hope they are unified... jun2012takao
!! ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,ikp) is added on jun2012takao
      implicit none
      integer:: nqnum,nqbz,imxc,ngcmx,nqbz_
      integer,allocatable:: ngc(:)
      real(8),allocatable:: qtt_(:,:)
      integer,allocatable:: ngvecc(:,:,:), ngveccrev(:,:,:,:),itimermap(:,:)
      contains
!!--- Readin QGcou ---
      subroutine readqgcou()
      implicit none
      integer:: ifiqg,iopen,nqi,ikp,nnnn,iclose,igc,ikpm,ikpp
      real(8):: QpGcut_cou,qmm(3)
      ifiqg  = iopen('QGcou',0,0,0)
      read(ifiqg) nqnum , ngcmx, QpGcut_cou, nqbz_, nqi,imxc
      write(6,*)'read(ifiqg)', nqnum , ngcmx, QpGcut_cou, nqbz_,nqi,imxc
      allocate( qtt_(3,nqnum),ngc(nqnum) )
      allocate( ngvecc(3,ngcmx,nqnum))
      allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,nqnum) )
      allocate( itimermap(ngcmx,nqnum))
      do ikp = 1,nqnum
c         print *,' ikp=',ikp,imxc
        read (ifiqg) qtt_(:,ikp), ngc(ikp)
c        print *,' ngc=',ngc(ikp)
        read (ifiqg) ngvecc(1:3, 1:ngc(ikp),ikp),ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,ikp)
!! time reversal mapping. only needed for time reversal case.
      enddo
c takao is trying to a mapping for time-reversal. But not yet. it can be complicated
c when we used current q-point mesh where -q do not exist in the list of q vector.
c      do ikp=1,nqnum
c        call iqindx2(-qtt_(:,ikp),ginv,qtt_,nqnum, ikpm,qmm) !qinv is true q
c        do igc=1,ngc(ikp)
c!!  ikp ---> ikpm
c!!  q+G (igc,ikp) is mapped to  qmm+G (itimermap,ikpm)
c!!
ccxxxxxxxxxxxxxxxxxxx qmm= -q + delta G xxxxxxxx 
ccxxxxxxxx this should be taken into account xxxxxxxx
ccxxxxxxxxx in addition, we have phase problem when there is delta G vector for inversion xxxxxx.
c          itimermap(igc,ikp)= ngveccrev(-ngvecc(1,igc,ikp),-ngvecc(2,igc,ikp),-ngvecc(3,igc,ikp),ikpm)
c          !ikp is mapped to itimermap
c        enddo  
c      enddo  
      ifiqg=iclose('QGcou')
      end subroutine readqgcou
      end module m_readqgcou

!! ===================================================
      module m_pbindex !generate index table for product basis.
      integer,parameter::NULLI=-9999999
      integer:: norbt=NULLI !number of PB block. 
      integer,allocatable:: 
     &    ibas_tbl(:),l_tbl(:),k_tbl(:),offset_tbl(:),offset_rev_tbl(:,:,:)
      integer::
     &    max_ibas_tbl=NULLI,max_l_tbl=NULLI,max_k_tbl=NULLI,max_offset_tbl=NULLI
      logical::init=.true.
      contains !-----------------------------------------------
      subroutine PBindex(natom,lx,l2nl,nx)
      implicit none
      integer::natom,l2nl,ibas,lb,nb,iorbt,kn,lx(natom),nx(0:l2nl,natom),mb
      if(.not.init) return
      init=.false.
      iorbt=0
      do ibas=1, natom
      do lb  = 0, lx (ibas)
      do nb  = 1, nx (lb,ibas)
        iorbt=iorbt+1
      enddo
      enddo
      enddo
      norbt = iorbt !number of product basis block
      print *,'npbindex: norbt=',norbt
      allocate( ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt) )
      iorbt=0
      offset_tbl(1)=0
      do ibas=1, natom
      do lb = 0, lx (ibas)
      do kn = 1, nx (lb,ibas)
          iorbt=iorbt+1
          ibas_tbl(iorbt)=ibas
          l_tbl(iorbt)=lb
          k_tbl(iorbt)=kn
          if(iorbt<norbt) offset_tbl(iorbt+1)=offset_tbl(iorbt)+2*lb+1
c        do mb  = -lb, lb
c          i = i+1  !The number of product basis is  =(i at the end of loop).
c          write(6,"(' === product basis index: iorbt ibas l n m',10i4)")iorbt,ic,lb,kn
c        enddo
      enddo
      enddo
      enddo
      max_ibas_tbl=natom
      max_l_tbl=maxval(l_tbl(1:norbt))
      max_k_tbl=maxval(k_tbl(1:norbt))
      allocate(offset_rev_tbl(max_ibas_tbl,0:max_l_tbl,max_k_tbl))
      do iorbt=1,norbt
        ibas=ibas_tbl(iorbt)
        lb= l_tbl(iorbt)
        kn= k_tbl(iorbt)
        offset_rev_tbl(ibas,lb,kn)=offset_tbl(iorbt)
c        write(6,"(' === product basis index: iorbt ibas l n offset',10i4)")iorbt,ibas,lb,kn,offset_rev_tbl(ibas,lb,kn)
      enddo
      end subroutine
      end module m_pbindex

!! ===================================================
      module m_rotMPB !WARN. rotMPB looks work OK. Used in eibz4x0()=T in switch.F  However,
!!---------------------------------------------------
!! == Mixed product basis rotator. ==
!! a curious strange segmention fault occurs in a Si444 case when we test a code qout= matmul(symops(:,:,igx),qin). 
!! I don't know why. 
!! This has developed for eibz mode jun2012, not tested completely. 
!! In future, this routine can be for some application.
!!takao2012jun
      use m_pbindex,only: norbt, ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
     &    max_ibas_tbl,max_l_tbl,max_k_tbl,max_offset_tbl
      use m_hamindex, only: qlat,plat,invgx,
     &  miat,tiat,shtvg,symops,nbas,ngrp
      use m_readqgcou,only: imxc,ngvecc,qtt_,nqnum,ngc,ngveccrev
      integer:: lxxa
      real(8),allocatable:: dlmm(:,:,:,:)
!! --- 
      contains
c$$$      subroutine rotMPB(zcousq,nbloch,ngbb,qin,igx,igxt,ginv
c$$$     o ,zcousqr) !zcousqr=Rotate_igx(zcousq) igxt=-1 means timereversal case.
c$$$      implicit none
c$$$      real(8):: qin(3),ginv(3,3),platt(3,3),qout(3),qu(3)
c$$$      integer:: igx,igxt,ngbb,iqin,iqout,ngcx,nbloch,iqindx,nl,i,j
c$$$      complex(8):: zcousq(ngbb,ngbb),zcousqr(ngbb,ngbb)
c$$$      integer,save:: init=1
c$$$      if(init==1) then
c$$$        lxxa= 2*max_l_tbl
c$$$        nl= lxxa +1 
c$$$        allocate(dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa,ngrp))
c$$$        call rotdlmm(symops,ngrp,nl,dlmm)
c$$$        init=0
c$$$      endif   
c$$$      call rotmto(qin,zcousq(1:nbloch,1:ngbb),nbloch,ngbb,
c$$$     i  norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl,
c$$$     i  ginv,symops(:,:,igx),shtvg(:,igx),dlmm(:,:,:,igx),lxxa,miat(:,igx),tiat(:,:,igx),igxt,nbas, 
c$$$     o  zcousqr(1:nbloch,1:ngbb))
c$$$c      print *,' end of rotmto '
c$$$      call iqindx2(qin,ginv,qtt_,nqnum, iqin,qu)
c$$$      if(sum(abs(qin-qu))>1d-6) stop 'rotMPB:qin is not included in QGcou'
c$$$      call iqindx2(matmul(symops(:,:,igx),qin),ginv,qtt_,nqnum, iqout,qout)
c$$$!!
c$$$      ngcx=ngc(iqin)
c$$$      if(ngcx/=ngc(iqout).or.ngcx/=ngbb-nbloch) then
c$$$         print *,' ngbb qin  =',ngbb,qin
c$$$         print *,' iqin iqout=',iqin,iqout
c$$$         print *,' qout=',qout
c$$$         print *,' ngc(iqin),ngc(iqout)=',ngc(iqin),ngc(iqout)
c$$$         print *,' ngcx nbloch ngbb=',ngcx,nbloch,ngbb
c$$$         print *,' ngc(iqout)=',ngc(iqout)
c$$$         if(ngcx/=ngbb-nbloch) stop 'rotMPB:ngc(iqin)/=ngbb-nbloch'
c$$$         stop 'rotMPB:ngc(iqin)/=ngc(iqout)'
c$$$      endif   
c$$$      platt=transpose(plat)
c$$$c      print *,' iqin iqout igx=',iqin,iqout,igx
c$$$c      print *,' goto  rotipw ',sum(zcousq(nbloch+1:nbloch+ngcx,1:ngbb))
c$$$c      print *,'qin  =',qin,iqin,igx
c$$$c      print *,'qout =',qout,iqout
c$$$c      print *,'ngcx ngbb=',ngcx,ngbb
c$$$c      print *,'platt=',platt
c$$$c      print *,' qlat=',qlat
c$$$      call rotipw(qin,qout,zcousq(nbloch+1:nbloch+ngcx,1:ngbb),ngcx,ngbb,
c$$$     &  platt,qlat,symops(1,1,igx),ngvecc(1,1,iqin),ngveccrev(:,:,:,iqout),shtvg(:,igx),igxt,imxc,
c$$$     o  zcousqr(nbloch+1:nbloch+ngcx,1:ngbb) )
c$$$c      if(igxt<0) then !time-reversal
c$$$c        zcousqr = dconjg(zcousqr)
c$$$c      endif  
c$$$      end subroutine rotMPB
!! ------------------------------------
      subroutine rotMPB2(nbloch,ngbb,qin,igx,igxt,ginv
     o ,zrotm) !zcousqr=Rotate_igx(zcousq) igxt=-1 means timereversal case.
!! --- zrotm(J,J') = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>. ---
!! See Eq.(51) around in PRB81 125102(2010). 
!! Eq. (51) can be written as
!!   P_IJ= \sum_i T_alpha_i [ zrotm_i^dagger (I,I') P'_I'J' zrom_i(J'J) ],
!! where P'_I'J' is ths sum not yet symmetrized. 

!!
!! Eexactrly speaking, we insert conversion matrix  between Enu basis and M_I basis.
!! It is zcousq(or zzr) in x0kf_v4h.F
!!
!! input qin = q
!! \hat{A}^k_i  is specified by symops(:,:,igx),and igxt (-1 for time-reversal).
!! Note that k= \hat{A}^k_i(k) (S_A^k)
      implicit none
      real(8):: qin(3),ginv(3,3),platt(3,3),qout(3),qu(3),sss
      integer:: igx,igxt,ngbb,iqin,iqout,ngcx,nbloch,nl,i,j
      complex(8):: zrotm(ngbb,ngbb)
      integer,save:: init=1
      logical ::debug=.false.
      if(debug) print *,' rotMPB2:'
      if(init==1) then
        lxxa= 2*max_l_tbl
        nl= lxxa +1 
        allocate(dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa,ngrp))
        call rotdlmm(symops,ngrp,nl,dlmm)
        init=0
      endif   
      zrotm=0d0
      call rotmto2(qin,nbloch,ngbb,
     i  norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl,
     i  ginv,symops(:,:,igx),shtvg(:,igx),dlmm(:,:,:,igx),lxxa,miat(:,igx),tiat(:,:,igx),igxt,nbas, 
     o  zrotm(1:nbloch,1:nbloch))
      if(debug) print *,' end of rotmto2 '
      call iqindx2(qin,ginv,qtt_,nqnum, iqin,qu) !iqindx2 is slow. If required, speed it up in a manner as iqindx2_.
      if(debug) print *,' end of iqindx2rotmto2 '
      if(sum(abs(qin-qu))>1d-6) stop 'rotMPB:qin is not included in QGcou'
      sss=1d0
      if(igxt==-1) sss=-1d0
      if(debug) print *,' iqindx2 111111111111 sss=',sss,igx,qin
      call iqindx2(sss*matmul(symops(:,:,igx),qin),ginv,qtt_,nqnum, iqout,qout)
      if(debug) print *,' iqindx2 111111111111222222222222'
!!
      ngcx=ngc(iqin)
      if(ngcx/=ngc(iqout).or.ngcx/=ngbb-nbloch) then
         print *,' ngbb qin  =',ngbb,qin
         print *,' iqin iqout=',iqin,iqout
         print *,' qout=',qout
         print *,' ngc(iqin),ngc(iqout)=',ngc(iqin),ngc(iqout)
         print *,' ngcx nbloch ngbb=',ngcx,nbloch,ngbb
         print *,' ngc(iqout)=',ngc(iqout)
         if(ngcx/=ngbb-nbloch) stop 'rotMPB2:ngc(iqin)/=ngbb-nbloch'
         stop 'rotMPB2:ngc(iqin)/=ngc(iqout)'
      endif   
      platt=transpose(plat)
c      print *,' iqin iqout igx=',iqin,iqout,igx
c      print *,' goto  rotipw ',sum(zcousq(nbloch+1:nbloch+ngcx,1:ngbb))
c      print *,'qin  =',qin,iqin,igx
c      print *,'qout =',qout,iqout
c      print *,'ngcx ngbb=',ngcx,ngbb
c      print *,'platt=',platt
c      print *,' qlat=',qlat
      call rotipw2(qin,qout,ngcx,ngbb,
     &  platt,qlat,symops(1,1,igx),ngvecc(1,1,iqin),ngveccrev(:,:,:,iqout),shtvg(:,igx),igxt,imxc,
     o  zrotm(nbloch+1:nbloch+ngcx,nbloch+1:nbloch+ngcx) )
c      if(igxt<0) then !time-reversal
c        zcousqr = dconjg(zcousqr)
c      endif  
      end subroutine rotMPB2
      end module m_rotMPB 



c$$$      module m_spacegrouprot
c$$$      integer(4),allocatable :: invgx(:), miat(:,:)
c$$$      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)
c$$$      contains !-------------------------------------------
c$$$      subroutine Spacegrouprot(symgg,ngrp,plat,natom,pos) !all input
c$$$      implicit none
c$$$      real(8)::symgg(*),plat(*),pos(*)
c$$$      integer::ngrp,natom,ificlass,iopen,ibas,ibasx,iclose
c$$$      integer,allocatable:: iclasst(:)
c$$$!!  Get data for space-group transformation. See header of mptaouof.
c$$$      ificlass   = iopen('CLASS',1,-1,0)
c$$$      allocate(iclasst(natom))
c$$$      print *,'  --- Spacegrouprot: Readingin true CLASS info ---'
c$$$      do ibas = 1,natom
c$$$        read(ificlass,*) ibasx, iclasst(ibas)
c$$$        write(6, "(2i10)") ibasx, iclasst(ibas)
c$$$      enddo
c$$$      ificlass= iclose('CLASS')
c$$$      allocate( invgx(ngrp),miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
c$$$      call mptauof(symgg,ngrp,plat,natom,pos,iclasst
c$$$     o     ,miat,tiat,invgx,shtvg )
c$$$      end subroutine
c$$$      end module 
c$$$!! ------------


      subroutine rotmto(qin,cphiin,nbloch,nband,
     i  norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl,
     i  ginv,sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,
     o  cphiout)
!!== Rotation of a function of ProductBasis or MT part of wave function.
!! (then MTO itself is rotated. a little diffent from here).
!! See also mto part of lm7K/subs/rotwv.F
!!
!! qin  = qtti(:,iqi)
!! cphiin cphi(1:ldim2,1:nband,iqi,isp) eigenfunction @ qtti(:,iqi)
!! dlmm   dlmm(-l:l,-l:l,l,igg)
!! sym   symops(:,:,igg)
!! tiat   tiat(:,ibas,igg)
!! shtvg  shtvg(:,igg)
!! miat   miat(ibas,igg)
!!
!!  cphin (iqq) ---> cphiout (iq)
      implicit none
      integer:: ibas,ibaso,nbloch,nband,norbt,iorb,nbas
      integer:: ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt)
      integer:: l,k,ini1,iend1,ini2,iend2,miat(nbas),lxxa
      real(8):: tiat(3,nbas)
      complex(8):: cphiin(nbloch,1:nband),cphiout(nbloch,1:nband),phase,ddd
      real(8):: sym(3,3), qin(3),qrot(3),shtvg(3) !shtvg
      real(8):: dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa),ginv(3,3)
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
      
      integer:: max_ibas_tbl,max_l_tbl,max_k_tbl,igxt
      integer::  offset_rev_tbl(max_ibas_tbl, 0:max_l_tbl, max_k_tbl)

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR

      qrot = matmul(sym,qin)
      if(igxt==-1) qrot=-qrot !july2012takao
      ibaso=-999
      do iorb=1,norbt !orbital-blocks
c         print *,'iorb=',iorb
        ibas = ibas_tbl(iorb)
        if(ibas/=ibaso) then
          phase = exp( -img2pi*sum(qrot*tiat(:,ibas))) 
        endif   
        ibaso=ibas
        l = l_tbl(iorb)
        k = k_tbl(iorb)
        ini1 = offset_tbl(iorb)+1
        iend1 = ini1+2*l
        ini2 = offset_rev_tbl(miat(ibas),l,k)+1
        iend2 = ini2+2*l
c        write(6,"(' === PBs index: iorb ibas l k ini1 ini2',10i4)")iorb,ibas,miat(ibas),l,k,ini1,ini2
        cphiout(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l),cphiin(ini1:iend1,:))*phase
      enddo  
      end subroutine rotmto
!! ------------------------------------------------------
c      subroutine rotmto2(qin,cphiin,nbloch,nband,
      subroutine rotmto2(qin,nbloch,nband,
     i  norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl,
     i  ginv,sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,
c     o  cphiout)
     o  zrotm)
!!== Rotation of a function of ProductBasis or MT part of wave function.
!! (then MTO itself is rotated. a little diffent from here).
!! See also mto part of lm7K/subs/rotwv.F
!!
!! qin  = qtti(:,iqi)
!! cphiin cphi(1:ldim2,1:nband,iqi,isp) eigenfunction @ qtti(:,iqi)
!! dlmm   dlmm(-l:l,-l:l,l,igg)
!! sym   symops(:,:,igg)
!! tiat   tiat(:,ibas,igg)
!! shtvg  shtvg(:,igg)
!! miat   miat(ibas,igg)
!!
!!  cphin (iqq) ---> cphiout (iq)
      implicit none
      integer:: ibas,ibaso,nbloch,nband,norbt,iorb,nbas
      integer:: ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt)
      integer:: l,k,ini1,iend1,ini2,iend2,miat(nbas),lxxa
      real(8):: tiat(3,nbas)
c      complex(8):: cphiin(nbloch,1:nband),cphiout(nbloch,1:nband),phaseq,phase,ddd
      complex(8):: phase,ddd
      real(8):: sym(3,3), qin(3),qrot(3),shtvg(3) !shtvg
      real(8):: dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa),ginv(3,3)
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
      integer:: max_ibas_tbl,max_l_tbl,max_k_tbl,igxt
      integer::  offset_rev_tbl(max_ibas_tbl, 0:max_l_tbl, max_k_tbl),m1,m2
!!
c      integer:: nrotm,irotm1(nrotmmx),irotm2(nrotmx)
      integer:: nrotm,irotm1,irotm2
      complex(8):: zrotm(nbloch,nbloch)

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR

      qrot = matmul(sym,qin)
      if(igxt==-1) qrot=-qrot !july2012takao
      ibaso=-999
      do iorb=1,norbt !orbital-blocks
        ibas = ibas_tbl(iorb)
        if(ibas/=ibaso) then
          phase = exp( -img2pi*sum(qrot*tiat(:,ibas)) ) 
        endif   
        ibaso=ibas
        l = l_tbl(iorb)
        k = k_tbl(iorb)
        ini1 = offset_tbl(iorb)+1
        iend1 = ini1+2*l
        ini2 = offset_rev_tbl(miat(ibas),l,k)+1
        iend2 = ini2+2*l
c        write(6,"(' === PBs index: iorb ibas l k ini1 ini2',10i4)")iorb,ibas,miat(ibas),l,k,ini1,ini2
!! this is the case ofr rotmto
c        cphiout(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l),cphiin(ini1:iend1,:))*phase*phaseq

!! We can get rotated rcxq as
!!      Rotated(rcxq)=  dconjg(chpiin(xx)) *  dconjg(zrotm)*rcxq(xx,irotm2) * zrotm(irotm2,irotm1)  * chpiin(irotm1)
        do m1=-l,l
        do m2=-l,l
          irotm1= ini1+m1+l
          irotm2= ini2+m2+l
          zrotm(irotm2,irotm1) = dlmm(m2,m1,l)*phase
        enddo
        enddo
      enddo
c      nrotm=ix
      end subroutine rotmto2
!! --------------------------------
c      subroutine rotipw2(qin,qtarget,geigenin,ngp,nband,
      subroutine rotipw2(qin,qtarget,ngp,nband,
     &  platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx,
c     o  geigenout)
     o  zrotm)
!! == Rotation of Plane wave part. by sym==
!!
!! Mapped from qtt(:,iqq) to qtt(:,iq)
!!   qtt(:,iq)= matmul(sym(igg),qtt(:,iqq))+some G vector
!!  geigenin (iqq) ---> geigenout (iq)
      implicit none
      real(8):: sym(3,3),qlat(3,3),platt(3,3)
     &  ,shtvg(3)
      integer:: ngp,imx,nband,
     &          ngvecp(3,ngp), !,iqq),
     &          ngvecprev(-imx:imx,-imx:imx,-imx:imx) !,iq)
      real(8):: qin(3),qpg(3),qpgr(3),qtarget(3)
      integer:: ig,ig2,nnn(3),igxt
c      complex(8):: geigenin(ngp,nband), geigenout(ngp,nband)
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
      complex(8):: zrotm(ngp,ngp)

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR

!!     qout = matmul(sym,qin)+some Gvector !qout \equiv qtt(:,iq) , qin=qtt(:,iqq)
      do ig = 1,ngp !ngp(iqq)
        qpg = qin + matmul( qlat(:,:),ngvecp(:,ig)) !iqq !q+G     qin=qtt(:,iqq)
        qpgr = matmul(sym,qpg)             !rotated q+G
        if(igxt==-1) qpgr=-qpgr !july2012
        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
c        if(maxval(abs(nnn))>imx) stop 'rotipw:maxval(abs(nnn))>imx'
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))  !,iq)
c        geigenout(ig2,:)= geigenin(ig,:) * exp( -img2pi*sum(qpgr*shtvg) )
        zrotm(ig2,ig)= exp( -img2pi*sum(qpgr*shtvg) )
      enddo
      end subroutine rotipw2
!! ------------------------------------------------------
!! --------------------------------
      subroutine rotipw(qin,qtarget,geigenin,ngp,nband,
     &  platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx,
     o  geigenout)
!! == Rotation of Plane wave part. by sym==
!!
!! Mapped from qtt(:,iqq) to qtt(:,iq)
!!   qtt(:,iq)= matmul(sym(igg),qtt(:,iqq))+some G vector
!!  geigenin (iqq) ---> geigenout (iq)
      implicit none
      real(8):: sym(3,3),qlat(3,3),platt(3,3)
     &  ,shtvg(3)
      integer:: ngp,imx,nband,
     &          ngvecp(3,ngp), !,iqq),
     &          ngvecprev(-imx:imx,-imx:imx,-imx:imx) !,iq)
      real(8):: qin(3),qout(3),qpg(3),qpgr(3),qtarget(3)
      integer:: ig,ig2,nnn(3),igxt
      complex(8):: geigenin(ngp,nband), geigenout(ngp,nband)
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR

!!     qout = matmul(sym,qin)+some Gvector !qout \equiv qtt(:,iq) , qin=qtt(:,iqq)
      do ig = 1,ngp !ngp(iqq)
        qpg = qin + matmul( qlat(:,:),ngvecp(:,ig)) !iqq !q+G     qin=qtt(:,iqq)
        qpgr = matmul(sym,qpg)             !rotated q+G
        if(igxt==-1) qpgr=-qpgr !july2012
        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
c        if(maxval(abs(nnn))>imx) stop 'rotipw:maxval(abs(nnn))>imx'
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))  !,iq)
        geigenout(ig2,:)= geigenin(ig,:) * exp( -img2pi*sum(qpgr*shtvg) )
      enddo
      end subroutine rotipw

!! --------------------------------------------------------
      module m_iqindx_qtt
      use m_hamindex, only: qtt,nqtt
      integer,allocatable,private:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      real(8),private:: epsd=1d-7
      integer,private:: nkey(3)
      real(8),private:: ginv(3,3)
      contains
!!---------------------------------------------------------
      subroutine iqindx2_(q,  iqindx,qu)
!!=== Find index as q=qq(:,iq) with modulo of premitive vector.===
!! ginv is the inverse of plat (premitive translation vector).
!! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
      implicit none
      integer:: i_out, iq,iqx,iqindx ,kkk3(3),ik1,ik2,ik3
      real(8):: q(3),qx(3),qu(3),qzz(3)
      logical::debug=.false.
      if(debug) write(*,"(' iqindx2_: q=',3f20.15)") q
      call rangedq(matmul(ginv,q), qzz)
      kkk3 = (qzz+0.5*epsd)/epsd
c      print *,'kkk3=',kkk3
c      print *,'nkey=',nkey
c      print *,'kk1=',kk1
      call tabkk(kkk3(1), kk1,nkey(1), ik1)
      call tabkk(kkk3(2), kk2,nkey(2), ik2)
      call tabkk(kkk3(3), kk3,nkey(3), ik3)
c      print *,' ik1,ik2,ik3=',ik1,ik2,ik3
      iqindx = iqkkk(ik1,ik2,ik3)
      qu =qtt(:,iqindx)
      end subroutine
!!----
      subroutine init_iqindx_qtt(ginv_)
!! === mapping of qtt ===
!! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
!!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
!!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
      real(8):: ginv_(3,3)
      real(8),allocatable:: qxx(:,:)
      integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3
      integer,allocatable:: ieord(:)
      ginv=ginv_
      allocate(ieord(nqtt))
      allocate(key(3,nqtt),qxx(3,nqtt))
      key=-99999
      do iq=1,nqtt
         call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
      enddo
!! get key and nkey for each ix.
      do ix =1,3
         call sortea(qxx(ix,:),ieord,nqtt,isig)
         ik=0
         do i=1,nqtt
            kkk=(qxx(ix,ieord(i))+0.5*epsd)/epsd  !kkk is digitized by 1/epsd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
               stop 'iqindx: bug not sorted well'
            endif
         enddo
         nkey(ix)=ik
      enddo
      deallocate(ieord)
!!  key is reallocated. inverse mattping, iqkkk
      allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
      kk1(:) = key(1,1:nkey(1))
      kk2(:) = key(2,1:nkey(2))
      kk3(:) = key(3,1:nkey(3))
      deallocate(key)
      allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
c      print *,' initqqq nqtt=',nqtt
      do i=1,nqtt
         kkk3= (qxx(:,i)+0.5*epsd)/epsd !kkk is digitized by 1/epsd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
c         print *, ' ik1,ik2,ik3,iqkkk=',ik1,ik2,ik3,i
      enddo
      deallocate(qxx)
      end subroutine

!! ---
      subroutine tabkk(kkin, kktable,n, nout)
      integer:: nout,n, kkin, kktable(n),i,mm,i1,i2
      i1=1
      i2=n
      if(kkin==kktable(1)) then
         nout=1
         return
      elseif(kkin==kktable(n)) then
         nout=n
         return
      endif
      do i=1,n
         mm=(i1+i2)/2
         if(kkin==kktable(mm)) then
            nout=mm
            return
         elseif(kkin>kktable(mm)) then
            i1=mm
         else
            i2=mm
         endif
      enddo
c$$$      do i=1,n                                                                                  
c$$$         if(kkin==kktable(i)) then                                                              
c$$$            nout=i                                                                              
c$$$            return                                                                              
c$$$         endif                                                                                  
c$$$      enddo                                                                                     
      print *, i1,i2,kkin
      print *, kktable(i1),kktable(i2)
      stop 'takk: error'
      end subroutine
      end module



!! ------------------------------------------------------
      module m_readeigen
!!-- Return eigenvalus and eigenfunctions for given q and isp.--
!! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
!! qtti(1:3,nqi) :eivenvalues, eigenvectors are calculated only for them.
!!                  See lmfgw (q-vector with irr flag in QGpsi).
!!
!!* we have to call init_foobar to call readeval, readcphi, readgeig.
!!----------------
      use m_iqindx_qtt,only: iqindx2_, init_iqindx_qtt
      implicit none
c      integer(4),private :: nqtt  
c      real(8),allocatable,private:: qtt(:,:) <---this is in m_hamindex
      real(8),allocatable,private:: evud(:,:,:), ginv(:,:)
      logical,private:: init=.true.,init2=.true.,keepeig
      integer(4),allocatable,private:: ngp(:)
      integer(4),private:: nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   nnnn, nsp,ldim2,mrecg,ifcphi,ifgeig
      complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
      real(8),private:: leval
ctm
      integer(4),private :: nwf,ifcphi_mlw,ifgeig_mlw, nqixx
      complex(8),allocatable,private:: geig_mlw(:,:,:,:),
     &           cphi_mlw(:,:,:,:)
      logical,private:: debug=.false.
      integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
      integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:)
     & ,offset_tbl(:),offset_rev_tbl(:,:,:)
      contains
!!--------------------------------------------------------------------------
      subroutine readeval(q,isp, ev)
!! == return ev(1:nband) for given q(1:3) and isp ==
!! *
      use m_hamindex, only: qtt,nqtt,iqimap
      implicit none
      integer(4):: isp,iq,iqindx,i
      real(8):: q(3),qu(3)
      real(8):: ev(nband)
      if(init) then
        stop 'readeigen: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
       call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
      ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
      if(debug) then
        print *,'iq iqimap(iq)=',iq,iqimap(iq)
        write(6,"('iq iqimap(iq)=',2i4,3f13.5)")iq,iqimap(iq),q
        write(6,"(9f9.4)")ev(1:9)
      endif
      end subroutine

      subroutine findemaxmin_r(emax,emin)
      implicit none
      real(8):: emax,emin,eee
      integer:: ib,iq,isp
      emax=-1d9
      do ib=1,nband
         do iq=1,nqixx
            do isp=1,nsp
               eee= evud(ib,iq,isp)
               if(eee>Emax.and.eee<1d9) Emax=eee
            enddo
         enddo
      enddo
      emin = minval(evud)
      end subroutine

!!-----------------------------------------------
      subroutine readgeig(q,ngp_in,isp, qu,geigen)
!!-- return ev(1:nband) for given q(1:3) and isp --
      use m_hamindex, only: iqmap,igmap,qtt,ngpmx,nqtt,shtvg,qlat,symops,plat,invgx,iqimap,imx
      implicit none
      integer:: isp,iq,iqindx,ngp_in,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi,igxt
      real(8)   :: q(3),qu(3),ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3),qout(3),qin(3)
      complex(8):: geigen(ngp_in,nband),geigenr(ngp_in,nband),img=(0d0,1d0),img2pi
      img2pi=2d0*4d0*datan(1d0)*img
      platt=transpose(plat) !this is inverse of qlat
      if(init2) then
        stop 'readgeig: modele is not initialized yet'
      endif
!! nqtt=nqnum. (aug2012. It was nqtt/2 before aug2012 when we use shorbz). See lm7k/fp/suham.F
c      call iqindx2(q,ginv,qtt,nqtt/2, iq, qu) !qu is used q. q-qu is a G vector.
      call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.

      if(debug) print *,' readgeig:xxx iq=',iq
c      do iq=1,nqtt
c         write(*,*)iq,qtt(:,iq)
c      enddo
      iqq=iqmap(iq)
      iqi=iqimap(iq)
      igg=igmap(iq)
      qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
      if(debug) then
        print *,' ddddd ginv=',ginv
        print *,' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        print *,' ddddd q      =',q
        print *,' ddddd qu      =',qu
        print *,' ddddd qtarget=',qtarget
        print *,' ddddd q (iqq)=',qtt(:,iqq)
        print *,' nqp(iq) ngp(iqq)=',ngp(iq),ngp(iqq)
      endif
!!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
      if(ngp(iq)==0) return
      if(ngp(iq)/=ngp(iqq)) stop 'readgeig:x ngp(iq)/=ngp(iqq)'
      if(debug) print *,' ngp=',ngp(iq)
      if(keepeig) then
        geigenr(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifgeig, rec=ikpisp) geigenr(1:ngpmx,1:nband)
      endif
      if(ngp_in < ngp(iq)) then
        print *,'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        stop 'readgeig: ngpmx<ngp(iq)'
      endif
!! new rotation of IPW. june2012
!!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
      igxt=1 !not timereversal
      call rotipw(qtt(:,iqq),qtt(:,iq),geigenr(1:ngp(iqq),1:nband),ngp(iqq),nband,
     &  platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx,
     &  geigen(1:ngp(iq),1:nband))
      return

!! Followings are not used now...  june2012
!! We can use the followings instead of rotipw (both should give the same result).
!! ----------
!! Rotation  See mto part of lm7K/subs/rotwv.F
!! NOTE that qrot = matmul(symops(:,:,igg),qtt(:,iqq))
!! rotation wave
c      ddd = shtvg(:,invgx(igg))
cccccccccccccccccccccccccccccccccc
c      print *,' xddddd1= ',-matmul(symops(:,:,igg),ddd)
c      print *,' xddddd2= ',shtvg(:,igg)
c      print *,' xddddd   '
ccccccccccccccccccccccccccccccccccc
c      if(debug) then
c        print *,' ddddd=',igg,ddd
c        write(*,"('ddddd iq =',2i5,3f13.5)") iq,  ngp(iq), qtt(:,iqq)
c        write(*,"('ddddd iqq=',2i5,3f13.5,2x,3f13.5)") iqq, ngp(iq), q,qu
c        print *,' ddddd '
c      endif
      do ig = 1,ngp(iqq)
        qpg = qtt(:,iqq) + matmul( qlat(:,:),ngvecp(:,ig,iqq)) !q+G
        qpgr = matmul(symops(:,:,igg),qpg)             !rotated q+G
        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3),iq)
ccccccccccccccccccccccccccccccccccc
cddd
c        geigen(ig2,:)= geigenr(ig,:) * exp( img2pi*sum(qpg*ddd) )
        geigen(ig2,:)= geigenr(ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
ccccccccccccccccccccccccccccccccccc
      enddo
      end subroutine readgeig

!!
      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
      use m_hamindex, only: iqmap,igmap,qtt,ngpmx,nqtt,qlat,plat,invgx,iqimap,
     &  miat,tiat,dlmm,shtvg,symops,lxxa,nbas ,  norbtx
!!-- return mto part of eigenfunction for given q(1:3) and isp
      implicit none
      integer:: isp,iq,iqindx,ldim2_in,ikpisp,iqq,iorb,ibaso,ibas,k,l,
     & ini1,ini2,iend1,iend2,igg,ig,iqi,i,igxt
      real(8)   :: q(3),qu(3) ,qrot(3) ,qout(3)
      complex(8):: phase,cphifr(ldim2_in,nband),cphif(ldim2_in,nband),phaseatom !takao 1->*->nband
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR


      if(init2) stop 'readcphi: modele is not initialized yet'
!! nqtt=nqnum. (aug2012. It was nqtt/2 before aug2012 when we use shorbz). See lm7k/fp/suham.F
c      call iqindx2(q,ginv,qtt,nqtt/2, iq, qu) !qu is used q. q-qu is a G vector.
      call iqindx2_(q, iq, qu) !qu is used q. q-qu= G vectors. qu=qtt(:,iq)
!!  For given iq, iqq = iqmap(iq)
!! qtt(:,iq)= matmul(sympos(  ,igg),qtt(:,iqq))  and  qtt(:,iqq) = qtti(:,iqi)
      igg=igmap(iq)
      iqq=iqmap(iq)  ! mapped from qtt(:,iqq) to qtt(:,iq);
      ! qtt(:,iq)=matmul(sym(igg),qtt(:,iqq))+some Gvector(see iqindx2 above)
      iqi=iqimap(iq) ! iqi is index for irr.=1 (cphi calculated. See lqg4gw and lm7K/subs/suham.F:L200).
      ! qtt(:,iqq) = qtti(:,iqi) is satisfied.
      ! qtti(:,iqi) are only when eigenfunctions are calculated).
      if(keepeig) then
        cphifr(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifcphi, rec=ikpisp) cphifr(1:ldim2,1:nband)
      endif
      if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ldim2_in,ldim2,norbtx,
     &                sum(cphifr(1:ldim2,1:nband)),sum(abs(cphifr(1:ldim2,1:nband)))
!!
      igxt=1 !not timereversal
      call rotmto(qtt(:,iqq),cphifr,ldim2,nband,
     i  norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
     i  maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl),
     i  ginv,symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lxxa,miat(:,igg),tiat(:,:,igg),igxt,nbas, 
     o  cphif)
      return

c$$$!! Followings are not used now...  june2012
c$$$!! We can use the followings instead of rotmto (both should give the same result).
c$$$!! ----------
c$$$!! Rotation  See mto part of lm7K/subs/rotwv.F
c$$$!! NOTE that qrot = matmul(symops(:,:,igg),qtt(:,iqq))
c$$$!! (then MTO itself is rotated. a little diffent from here).
c$$$!! qrot= symops(igg) qtt
c$$$      qrot=matmul(symops(:,:,igg),qtt(:,iqq)) !qrot-qtt(:,iq) is some G vector.
c$$$      phaseq   = exp( -img2pi*sum( qrot*shtvg(:,igg)) )
c$$$cddd
c$$$c      phaseq   = exp( img2pi*sum( qtt(:,iqq)*shtvg(:,invgx(igg)) ) )
c$$$cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c$$$cccccccccccccccccccccccccccccccccc
c$$$c      print *,' xddddd1= ', sum( qtt(:,iqq)*shtvg(:,invgx(igg)) )
c$$$c      print *,' xddddd2= ',-sum( qrot*shtvg(:,igg))
c$$$c      print *,' xddddd diff=   ',sum( qtt(:,iqq)*shtvg(:,invgx(igg)) )+sum( qrot*shtvg(:,igg))
c$$$ccccccccccccccccccccccccccccccccccc
c$$$      ibaso=-999
c$$$      do iorb=1,norbt !orbital-blocks are specified by ibas, l, and k.
c$$$        !Index of Hamiltonian is devided into these blocks.
c$$$        ibas = ibas_tbl(iorb)
c$$$        if(ibas/=ibaso) then
c$$$          phase = exp( -img2pi*sum(qrot*(tiat(:,ibas,igg)-shtvg(:,igg))) )
c$$$        endif
c$$$        ibaso=ibas
c$$$        l   = l_tbl(iorb)
c$$$        k   = k_tbl(iorb)
c$$$        ini1 = offset_tbl(iorb)+1
c$$$        iend1 = ini1+2*l
c$$$        ini2 = offset_rev_tbl(miat(ibas,igg),l,k)+1
c$$$        iend2 = ini2+2*l
c$$$        cphif(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),cphifr(ini1:iend1,:))
c$$$     &                       *phase*phaseq
c$$$        if(debug) then
c$$$          write(6,*)ibaso,ibas,phase
c$$$          write(6,*)ini1,iend1,ini2,iend2
c$$$          write(6,*)'sumd  dlmm =',igg,l,sum(abs(dlmm(-l:l,-l:l,l,igg)))
c$$$          write(6,*)'sumcr cphifr=',l,k,igg,sum(abs(cphifr(ini1:iend1,:)))
c$$$          write(6,*)'sumcf cphif=',l,k,igg,sum(abs(cphif(ini2:iend2,:)))
c$$$        endif
c$$$      enddo
c$$$ccccccccccccccccccccccccccccccccc
c$$$c       write(6,"('readcphi:: xxxx sum of cphifr=',2i4,4d23.16)")ldim2_in,ldim2,
c$$$c     &                sum(cphif(1:ldim2,1:nband)),sum(abs(cphif(1:ldim2,1:nband)))
c$$$c       stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxx takao test xxxxxxxxxxxxxxxxx'
c$$$cccccccccccccccccccccccccccccccccc
c$$$ccccccccccccccccccccccc
c$$$c      write(6,"('rrcc0:',i3,3f8.4)")iq,qu
c$$$c      do i=1,1
c$$$c      write(6,"('rrccb:',5(2x,2f8.4))")
c$$$c     &  cphif(14,i)/cphif(1,i) !*abs(cphif(1,i))
c$$$c      enddo
c$$$c      cphif(5:13, :)=0d0
c$$$c      cphif(14,:)=0d0
c$$$c      cphif(18:26,:)=0d0
c$$$
c$$$c      cphif(1:8,:)=0d0
c$$$c      cphif(19:26,:)=0d0
c$$$c      cphif(1,:)=0d0
c$$$c      cphif(19,:)=0d0
c$$$cccccccccccccccccccccccccccccccccccc
      end subroutine readcphi


c$$$c------------------------------
c$$$      subroutine readgeig(q,ngp_in,isp, qu,geigen)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: geigen(ngp_in,nband)
c$$$      if(init2) then
c$$$        stop 'readgeig: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(ngp_in < ngp(iq)) then
c$$$        print *,'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
c$$$        stop 'readgeig: ngpmx<ngp(iq)'
c$$$      endif
c$$$! for ngp=0 Aug2005
c$$$      if(ngp(iq)==0) return
c$$$      if(keepeig) then
c$$$        geigen(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifgeig, rec=ikpisp) geigen(1:ngpmx,1:nband)
c$$$      endif
c$$$      end subroutine readeig
c$$$c------------------------------
c$$$      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: cphif(ldim2_in,nband) !takao 1->*->nband
c$$$      if(init2) then
c$$$        stop 'readcphi: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(keepeig) then
c$$$        cphif(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifcphi, rec=ikpisp) cphif(1:ldim2,1:nband)
c$$$      endif
c$$$      end subroutine readcphi

c------------------------------
      subroutine readgeig_mlw(q,ngp_in,isp, qu,geigen)
      use m_hamindex
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) then
        stop 'readgeig_mlw: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgeig: q=',q
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.

      if(ngp_in < ngp(iq)) then
        print *,'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        stop 'readgeig_mlw: ngpmx<ngp(iq)'
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geig_mlw(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeig_mlw, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine

c------------------------------
      subroutine readcphi_mlw(q,ldim2_in,isp,  qu,cphif)
      use m_hamindex
!!-- return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2_in,nwf)
      if(init2) then
        stop 'readcphi_mlw: modele is not initialized yet'
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      print *,' readgcphi_mlw: q=',q
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.

      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphi(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine

ccccccccccccccccccccccccccc
      subroutine init_readeigen(ginv_in,nsp_in,nband_in,mrece_in)
      use m_hamindex,only: readhamindex,imx,ngpmx,nqi,nqnum  ,qtt,nqtt

!!-- initialization. Save QpGpsi EVU EVD to arrays.--
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,nsp_in, 
     & nband_in,mrece_in,iopen,iclose,isx,mrecb_in,ik,ib,verbose
      real(8):: ginv_in(3,3)
c      integer(4),allocatable:: ngvecp(:,:)
      real(8):: QpGcut_psi
      integer:: ifev,nband_ev, nqi_, nsp_ev !,nqnum
      real(8),allocatable:: qtt_(:,:),qtti_(:,:)

      write(6,*) 'init_readeigen:'
      allocate(ginv(3,3))
      ginv  = ginv_in
      nsp   = nsp_in
      if(nsp<0.or.nsp>2) stop 'init_reaeigen:nsp wrong'
      nband = nband_in
      mrece = mrece_in
c      write(6,*) 'init_readeigen:1'

!!--- store requied index
      call readhamindex()
c      write(*,*)ngrp,nbas,kxx,lxx
c      write(*,*)'nqtt=',nqtt,nqi,nqnum
c      write(*,*)'imx=',imx,ngpmx
c      do iq=1,nqtt
c        write(*,*) qtt(:,iq)
c      enddo
c      stop 'xxxxxxxxxxxxxxxxxx test xxxxx'

      call init_iqindx_qtt(ginv)

!!--- Readin QGpsi ---
      ifiqg  = 101
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum , ngpmx, QpGcut_psi, nnnn,nqi,imx
      write(6,*)'read(ifiqg)', nqnum , ngpmx, QpGcut_psi, nnnn
      allocate( qtt_(3,nqnum),ngp(nqnum) )
      allocate( ngvecp(3,ngpmx,nqnum))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
      do ikp = 1,nqnum
        read (ifiqg) qtt_(:,ikp), ngp(ikp)
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp)
      enddo
      close(ifiqg)
      deallocate(qtt_)

!!--- Read eigenvalue ---
      ifev= iopen('EValue',0,-1,0)    !feb2012takao
      read(ifev) nband_ev, nqi_, nsp_ev
      write(*,*)'read EValue: nband_ev,nqi,nsp_ev', nband_ev, nqi_, nsp_ev
      if(nband_ev/=nband) stop 'init_readeigen:nband_ev/=nband'
      if(nsp_ev  /=  nsp) stop 'init_readeigen:nsp_ev/=nsp'
      if(nqi  /=  nqi_) stop 'init_readeigen:nqi/=nqi_'
      nqixx=nqi
      allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
      read(ifev) qtti_(1:3,1:nqi_) !qtti
      read(ifev) evud(1:nband, 1:nqi, 1:nsp )
      close(ifev)
      if(debug) then
        do is= 1,nsp
          do ik= 1,nqi
            do ib= 1,nband
              write(6,"('ib ik e=',3i5,f13.5)") ib,ik,is,evud(ib,ik,is)
            enddo
          enddo
        enddo
        if(debug) print *,'init_readeigen:end'
      endif

c comment out mar2012takao
c      ifev(1) = iopen('EVU', 0,0,mrece)
c      if(nsp==2) ifev(2) = iopen('EVD', 0,0,mrece)
C$$$      do is= 1,nsp
C$$$        call rwdd (ifev(is),   !read eigenvalues
C$$$     d              nband,nqi,
C$$$     o              evud(1,1,is))
C$$$c        print *, 'evud is=',is,minval(evud(:,:,is))
C$$$        if(debug) then
C$$$          do ik=1,nqtt
C$$$            do ib=1,nband
C$$$              write(6,"('ib ik e=',2i5,f13.5)") ib,ik,evud(ib,ik,is)
C$$$            enddo
C$$$          enddo
C$$$        endif
C$$$      enddo
C$$$c      stop
C$$$      isx = iclose('EVU')
C$$$      if(nsp==2) isx = iclose('EVD')
      leval= minval(evud)
      init=.false.
      end subroutine

      real(8) function lowesteval()
      lowesteval=leval
      end function
c-------------------------------------------------------------
      subroutine init_readeigen2(mrecb_in,ldim2_in,mrecg_in)
      use m_hamindex,only: norbtx,nqi,ngpmx
c this should be called after init_readgeigen
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,ldim2_in, 
     & iopen,iclose,isx,mrecb_in,ikpisp,mrecg_in,verbose,ifoc,
     & i1,i2,i3,i4,i5,iorb,iorbold
      logical :: keepeigen

      mrecb = mrecb_in
      ldim2 = ldim2_in
      mrecg = mrecg_in
      keepeig = keepeigen()

      init2=.false.

!! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      read(ifoc,*)
      norbtx=0
c      write(*,*) ' goto readin @MNLA_CPHI: norbtx=',norbtx
      do
        read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
c        write(6,*) i1,i2,i3,i4,i5,iorb
        if(iorb>norbtx) norbtx=iorb
      enddo
 106  continue
      write(*,*) ' end of readin @MNLA_CPHI: norbtx=',norbtx
      rewind ifoc
      read(ifoc,*)
      allocate(l_tbl(norbtx),k_tbl(norbtx),ibas_tbl(norbtx),offset_tbl(norbtx))
      iorbold=0

      do
        read(ifoc,*,end=107)i1,i2,i3,i4,i5,iorb
        if(iorbold/=iorb) then
          k_tbl(iorb)=i2
          l_tbl(iorb)=i3
          ibas_tbl(iorb)=i4
          offset_tbl(iorb)=i5-1
          iorbold=iorb
        endif
      enddo
 107  continue

      allocate(offset_rev_tbl(maxval(ibas_tbl),0:maxval(l_tbl),maxval(k_tbl)))
      offset_rev_tbl=-99999
      do iorb=1,norbtx
        offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))= offset_tbl(iorb)
      enddo
      write(*,*)'init_readeigen2: table for rotation of cphi'
c      do iorb=1,norbtx
c        write(*,*)ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb),
c     &   offset_tbl(iorb),
c     &  offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))
c      enddo
      ifoc=iclose('@MNLA_CPHI')
c      stop 'init2 xxxxxxxxxxxxx'
c --- Readin CPHI and GEIG
      if(Keepeig     ) print *,' KeepEigen=T; readin geig and cphi into m_readeigen'
      if(.not.Keepeig) print *,' KeepEigen=F; not keep geig and cphi in m_readeigen'
      ifgeig = iopen('GEIG',0,0,mrecg)
      ifcphi = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction
      print *,' ifgeig ifcphi=',ifgeig,ifcphi
      if(.not.keepeig) return
c      ifgeig = iopen('GEIG',0,-1,0)
c      ifcphi = iopen('CPHI',0,-1,0) ! Coeeficients of eigenfunction
      allocate(geig(ngpmx,nband,nqi,nsp))
      allocate(cphi(ldim2,nband,nqi,nsp))
      do ikp= 1,nqi
        do is= 1,nsp
c        read(ifgeig) geig(1:ngpmx,1:nband,ikp,is) !geig(1:ngp(ikp),1:nband,ikp,1:is)
c        read(ifcphi) cphi(1:ldim2,1:nband,ikp,is)
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
ccccccccccccccccccccccc
c          if(verbose()>40)
c     &       write(6,'("Readin geig sumcheck ",i5,3f9.4,"  ",i5,2d13.6)')
c     &      ikp,qtti(:,ikp),ngp(ikp),sum(geig(1:ngp(ikp),1:nband,ikp,is))
cccccccccccccccccccccc
c        read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is)
          read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nband,ikp,is)
        enddo
      enddo
      isx= iclose('CPHI')
      isx= iclose('GEIG')

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw(nwf_o)
      use m_hamindex
c allocate cphi_mlw and geig_mlw
c for hx0fp0_mlw (ixc=112 mode)
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,ib,iwf,iwf2,iko_ix,iko_fx,in,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c allocate geig_mlw and cphi_mlw
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(dnk,uum)


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig_mlw(ngpmx,nwf,nqtt,nsp))
        allocate(cphi_mlw(ldim2,nwf,nqtt,nsp))
        geig_mlw = 0d0
        cphi_mlw = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig_mlw(:,iwf,ikp,is) = geig_mlw(:,iwf,ikp,is) + 
     &           geig(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi_mlw(:,iwf,ikp,is) = cphi_mlw(:,iwf,ikp,is) + 
     &           cphi(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

          enddo
        enddo


      else
        ifcphi_mlw = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_mlw = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_mlw, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_mlw, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
      endif


      deallocate(cbwf)


c replace dimensions
      nwf_o = nwf

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_noeval1D(
     &                 nwf_o,nband_o,mrecb_o,mrecg_o)
      use m_hamindex

c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2,
     & iset,nset,iset2,nset2,ini
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw_noeval1D'


c --- Readin MLWU/D and UU1dU/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifuu   = iopen('UU1dU',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifuu   = iopen('UU1dD',0,0,0)
        endif

c read UU1dU/D
        read(ifuu)nqbz,nq0i,nwf,nset

        iko_ix = 1
        iko_fx = nwf
        if (nqtt.ne.nset*nqbz) stop
     >      "init_readeigen_nowval1D: nqtt error"

c nqbz mesh-points

        if (is.eq.1) allocate(cbwf(nwf,nwf,nqtt,nsp))

        do iset = 1,nset
          ini = (iset-1)*nqbz
          read(ifmlw)nqbz2,nwf2,iset2,nset2
          if (nqbz.ne.nqbz2) stop
     >         "init_readeigen_nowval1D: nqbz error"
          if (nwf.ne.nwf2) stop
     >         "init_readeigen_nowval1D: nwf error"
          if (iset.ne.iset2) stop
     >         "init_readeigen_nowval1D: iset error"
          if (nset.ne.nset2) stop
     >         "init_readeigen_nowval1D: nset error"
          do iqbz = 1,nqbz
            read(ifmlw)iqbz2,q(1:3)
            if (iqbz2 .ne. iqbz)
     &           stop 'init_readeigen_mlw_noeval1D: iqbz error'
            read(ifmlw)cbwf(1:nwf,1:nwf,ini+iqbz,is)
          enddo
        enddo

c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifuu   = iclose('UU1dU')
        else
          ifmlw  = iclose('MLWD')
          ifuu   = iclose('UU1dD')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif


      deallocate(cbwf)


c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
      use m_hamindex

c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif

ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      deallocate(dnk,uum)


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif


      deallocate(cbwf)


c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_phi_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
      use m_hamindex

c replace cphi and geig
c for hwmat_phi
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbz2,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,nsp2,nlmto2,ngpmx2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,ifdim
      real(8):: q(3),rnorm,cnorm
      complex(8),allocatable :: geig2(:,:),cphi2(:,:)

      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_phi_noeval'

      ifdim = iopen('PHIG.d',1,0,0)
      read(ifdim,*) nsp2,nqbz2,nwf2,nlmto2,ngpmx2
      isx= iclose('PHIG.d')
      if(nsp2.ne.nsp) then
        write(*,*)'nsp,nsp2',nsp,nsp2
        stop 'init_readeigen_phi: ns'
      endif
      if(nlmto2.ne.ldim2) then
        write(*,*)'nlmto,nlmto2',ldim2,nlmto2
        stop 'init_readeigen_phi: nlmto'
      endif
      if(ngpmx2.ne.ngpmx) then
        write(*,*)'ngpmx,ngpmx2',ngpmx,ngpmx2
        stop 'init_readeigen_phi: ngpmx'
      endif
      nwf = nwf2
      nqbz=nqbz2

c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0

c replace geig and cphi

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband

      isx= iclose('CPHI')
      isx= iclose('GEIG')
      ifgeig = iopen('GEIGg',0,0,mrecg_o)
      ifcphi = iopen('CPHIg',0,0,mrecb_o)

      if(keepeig) then
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        do ikp= 1,nqtt
          iqbz=ikp
c         iqbz=mod(ikp,nqbz)
c         if (iqbz.eq.0) iqbz=nqbz
          do is= 1,nsp
            ikpisp = is + nsp*(iqbz-1)
            read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nwf,ikp,is)
            read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nwf,ikp,is)
          enddo
        enddo
      else
c         allocate(geig2(ngpmx,nwf),cphi2(ldim2,nwf))
c         do ikp= 1,nqtt
c           if (ikp.le.nqbz) cycle
c           iqbz=mod(ikp,nqbz)
c           if (iqbz.eq.0) iqbz=nqbz
c         do is= 1,nsp
c           ikpisp = is + nsp*(iqbz-1)
c           read(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           read(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c           ikpisp = is + nsp*(ikp-1)
c           write(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           write(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c         enddo
c         enddo
c         deallocate(geig2,cphi2)
      endif

c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o

c ---

      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_diag(nwf_o,nband_o,mrecb_o,mrecg_o)
      use m_hamindex

c replace cphi, geig and evud
c for hx0fp0_mlw, ixc=11 mode
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:),
     &                          dnk2(:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      print *,' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),
     >                        dnk2(iko_ix:iko_fx,nwf))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) stop "init_readeigen_mlw: nqbz2 error"
        if (iko_ix2.ne.iko_ix)
     >      stop "init_readeigen_mlw: iko_ix2 error"
        if (iko_fx2.ne.iko_fx)
     >      stop "init_readeigen_mlw: iko_fx2 error"
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c nqbze diag-points
        read(ifmlw)nqbze,nwf2
        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nqbze.ne.nqtt .and. nqbze.ne.nqtt-1)
     &      stop 'init_readeigen_mlw: nqbze error'
        read(ifmlwe)nqbze2,nwf2
        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nqbze2.ne.nqbze) stop 'init_readeigen_mlw: nqbze2 error'
        if (is.eq.1)allocate(evec(nwf,nwf,nqtt,nsp),eval(nwf,nqtt,nsp))
        evec = 0d0
        eval = 0d0
        do iqbz = 1,nqbze
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          read(ifmlw)evec(1:nwf,1:nwf,iqbz,is)

          read(ifmlwe)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) stop 'init_readeigen: iqbz error'
          read(ifmlwe)eval(1:nwf,iqbz,is)
        enddo

c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = eval
      deallocate(eval)

check write
c      do iq = 1,nqbze
c      write(6800,*)'***',iq
c      write(6800,*)evec(:,:,iq,1)
c      do iwf = 1,nwf
c         write(6800,*)evud(iwf,iq,1)
c      enddo
c      enddo

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqbze
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk(ib,in,iqbz,is)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
c
          else
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            dnk2 = 0
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                dnk2(ib,iwf) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
c
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk2(ib,in)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
          endif


ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      deallocate(dnk,evec,uum,dnk2)

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif

      deallocate(cbwf)

c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o
c ---
      end subroutine
      end module


c taken from lm7k/subs/rotwv.F
      subroutine conv2int(a,ii,n)
      implicit none
      integer::ii(n),i,n
      real(8):: a(n),eps=1d-6,aa !eps=1d-8 until aug2012
      do i=1,n
        if(a(i)>0) then
          aa = a(i)
          ii(i)=  int(aa+eps)
        else
          aa=-a(i)
          ii(i)= -int(aa+eps)
        endif
        if(abs(ii(i)-a(i))>eps)stop'conv2int: a is not close to integer'
      enddo
      end subroutine conv2int

c------------------------------------------------------
      subroutine rotwvigg(igg,q,qtarget,ndimh,napw_in,nband,evec,evecout,ierr)
      use m_hamindex,only: symops,invgx,miat,tiat,shtvg,qlat,plat,dlmm,ngrp,norbmto,
     &  ibastab,ltab,ktab,offl,offlrev,getikt,igv2,igv2rev,napwk
      implicit none
!! == Wave function rotator evec to evectout by symops(:,:,igg). ===
!! INPUT:
!!  need to supply qtarget=sym*q for check (you can modify this code to skip this.
!! OUTPUT evecout, ierr 
!! NOTE:
!! rotation of coefficients on PMT basis. 
!!  phi(r) = \sum_i evec(i,iband) |F_i> ==> Roteta[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
!!  Rotate[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
!!  miat tiat are atom-mapping table. mptauo.F
!! c.f. rotsig in lm7K/fp/rotwv.F also.
      integer   :: ig,ndimh,napw_in,nband,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k
      real(8)   :: q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3) !,plat(3,3),qlat(3,3)
      complex(8):: evec(ndimh,nband),evecout(ndimh,nband),phase,img=(0d0,1d0),img2pi
      img2pi=2*4d0*datan(1d0)*img
      platt = transpose(plat) !this is inverse of qlat
      ierr=1
!! check q is really rotated to qtarget by symops(:,:,igg)
      call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx)
      if(sum(abs(qx))>1d-8) then
        write(6,"(a,3f7.3,2x,3f7.3)")'  rotwvigg: qtarget is not a star of q',q,qtarget
        stop 'rotwvigg: qtarget is not symops(:,:,ig)*q'
      endif  
      evecout=0d0
      nlmto = ndimh-napw_in
!! mto part
c          write(6,"(3i3,2x,2i5,3x,a)")
c     &    ib,l,k, offl(iorb,ib)+1,offl(iorb,ib)+2*l+1,trim(spid)
c      print *,'rotwvigg :nlmto=',nlmto,norbmto
c      print *,' ibastab=',ibastab
      if(nlmto/=0) then
        ibaso=-999
        do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          !Index of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( -img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+1
          iend1 = offl(iorb)+2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+1
          iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
          evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase
        enddo
      endif
!! apw part
      if(napw_in/=0) then
        ikt  = getikt(q)    !index for q
        ikt2 = getikt(qtarget) !index for qtarget !you
        if(napw_in /= napwk(ikt) ) then
          call rx_('rotwv: napw_in /= napw(ikt)')
        endif
        do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)              !rotated q+G
          call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
        enddo
      endif
      ierr=0
      end

      subroutine rx_(string)
c  routine for error exit
      character*(*) string
      write(6,892) string
      write(6,890)
  890 format(' ---- Error exit')
  892 format(/' ---- ',a)
c ... write to log file also
      write(71,710) string
  710 format('++ ',a)
      call exit(-1)
      end
