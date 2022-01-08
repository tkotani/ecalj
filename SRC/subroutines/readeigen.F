!! -----------------------------------------------------------
!> Return eigenvalus and eigenfunctions for given q and isp.
!! We can get eigenfunctions for Wannier, as well. See hmagnon.F
!!
!! note: we have to call init_foobar to call readeval, readcphi, readgeig.
!!----------------
      module m_readeigen
      use m_iqindx_qtt,only: Iqindx2_, Init_iqindx_qtt
      use m_hamindex,only:   ngpmx, nqtt, nqi, nqnum, qtt,iqimap,
     &     iqmap,igmap,shtvg,qlat,symops,plat,invgx,
     &     miat,tiat,dlmm,shtvg,symops,lxxa,nbas 
      use m_read_bzdata,only: ginv
      use m_genallcf_v3,only: nsp =>nspin ,ldim2=>nlmto
      use m_readhbe,only: mrecb, mrece,nband, mrecg
!! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
!! qtti(1:3,nqi) :eivenvalues, eigenvectors are calculated only for them.
!!                  See lmfgw (q-vector with irr flag in QGpsi).
      implicit none
!-------------------------------
      public:: Readeval, Readgeigf, Readcphif, Readcphifq, Init_readeigen, Init_readeigen2, Lowesteval,
     &     Onoff_write_pkm4crpa,
     &     Init_readeigen_mlw_noeval, Init_readeigen_phi_noeval, Readcphiw, Readgeigw
      integer,public:: nwf
!------------------------------      
      private 
      integer:: ifgeigW,ifcphiW,norbtx,imx
      complex(8),allocatable:: geigW(:,:,:,:),cphiW(:,:,:,:)
      integer::nwf_ev
      complex(8),allocatable:: evud_w(:,:,:,:)
      real(8),allocatable:: eval_w(:,:,:)      
      logical:: evalwan=.true.
      real(8),allocatable,private:: evud(:,:,:)!, ginv(:,:)
      logical,private:: init=.true.,init2=.true.,keepeig
      integer,allocatable,private:: ngp(:)
      integer,private:: nprecb,
     &   nnnn, ifcphi,ifgeig !,mrecg,mrecb ,ldim2!,mrece,
      complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
      real(8),private:: leval
      integer,private :: ifcphi_mlw,ifgeig_mlw, nqixx
      complex(8),allocatable,private:: geig_mlw(:,:,:,:),
     &           cphi_mlw(:,:,:,:)
      logical,private:: debug=.false.
      integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
      integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:)
     & ,offset_tbl(:),offset_rev_tbl(:,:,:)
      logical,private:: Wpkm4crpa=.false.
      real(8),private :: quu(3)

      contains
      subroutine onoff_write_pkm4crpa(lll) 
      logical:: lll
      Wpkm4crpa=lll
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readgeigW(q,ngp_in,isp, qu,geigen)
      integer:: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) call rx( 'readgeig_mlw: modele is not initialized yet')
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        call rx( 'readgeig_mlw: ngpmx<ngp(iq)')
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geigW(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeigW, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readcphiW(q,ldim2_dummy,isp,  qu,cphif)
      integer:: isp,iq,iqindx,ldim2_dummy,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2,nwf)
      if(init2) call rx( 'readcphi_mlw: modele is not initialized yet')
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphiW(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
!> Return ev(1:nband) for given q(1:3) and isp 
      function readeval(q,isp) result(ev)
      intent(in)  ::    q,isp
      integer :: isp
      real(8) :: q(3), ev(nband)
      integer:: iq,iqindx,i
      real(8):: qu(3)
      if(init) call rx( 'readeigen: modele is not initialized yet')
c      iq = iqindx(q, ginv,qtt,nqtt)
       call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
      if(debug) then
        write(6,*)'iq iqimap(iq)=',iq,iqimap(iq)
        write(6,"('iq iqimap(iq) q=',2i8,3f13.5)")iq,iqimap(iq),q
      endif
      ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
      if(debug) then
        write(6,"(9f9.4)")ev(1:9)
      endif
      end function
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine findemaxmin_r(emax,emin)
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
!ssssssssssssssssssssssssssssssssssssssssssssss
!> Return ev(1:nband) for given q(1:3) and isp 
      function readgeigf(q,isp) result(geigen)
      real(8),intent(in):: q(3)
      integer,intent(in):: isp
      real(8):: qu(3)
      complex(8):: geigen(ngpmx,nband)
      call readgeig(q,ngpmx,isp,qu,geigen)
      end
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readgeig(q,ngp_in,isp, qu,geigen)
      use m_rotwave,only: Rotipw
      implicit none
      real(8),intent(in) :: q(3)
      integer,intent(in) :: isp,ngp_in
      real(8),intent(out) :: qu(3)
      complex(8), intent(out) :: geigen(ngp_in,nband)
      integer:: iq,iqindx,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi,igxt
      real(8)   :: ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3),qout(3),qin(3)
      complex(8):: geigenr(ngp_in,nband),img=(0d0,1d0),img2pi
      img2pi=2d0*4d0*datan(1d0)*img
      platt=transpose(plat) !this is inverse of qlat
      if(init2) call rx( 'readgeig: modele is not initialized yet')
!! nqtt=nqnum. (aug2012. It was nqtt/2 before aug2012 when we use shorbz). See lm7k/fp/suham.F
c      call iqindx2(q,ginv,qtt,nqtt/2, iq, qu) !qu is used q. q-qu is a G vector.
      call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
      if(debug) write(6,*)' readgeig:xxx iq=',iq
      iqq=iqmap(iq)
      iqi=iqimap(iq)
      igg=igmap(iq)
      qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
      if(debug) then
        write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        write(6,*)' ddddd q      =',q
        write(6,*)' ddddd qu      =',qu
        write(6,*)' ddddd qtarget=',qtarget
        write(6,*)' ddddd q (iqq)=',qtt(:,iqq)
        write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
      endif
!!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
      if(ngp(iq)==0) return
      if(ngp(iq)/=ngp(iqq)) then
        write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        write(6,*)' ddddd q      =',q
        write(6,*)' ddddd qu      =',qu
        write(6,*)' ddddd qtarget=',qtarget
        write(6,*)' ddddd q (iqq)=',qtt(:,iqq)
        write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
        call rx( 'readgeig:x2 ngp(iq)/=ngp(iqq)')
      endif  
      if(debug) write(6,*)' ngp=',ngp(iq)
      if(keepeig) then
        geigenr(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifgeig, rec=ikpisp) geigenr(1:ngpmx,1:nband)
      endif
      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        call rx( 'readgeig: ngpmx<ngp(iq)')
      endif
!! new rotation of IPW. june2012
!!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
      igxt=1 !not timereversal
      call rotipw(qtt(:,iqq),qtt(:,iq),geigenr(1:ngp(iqq),1:nband),ngp(iqq),nband,
     &  platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx,
     &  geigen(1:ngp(iq),1:nband))
      return
      do ig = 1,ngp(iqq)
        qpg = qtt(:,iqq) + matmul( qlat(:,:),ngvecp(:,ig,iqq)) !q+G
        qpgr = matmul(symops(:,:,igg),qpg)             !rotated q+G
c        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
        nnn=nint(matmul(platt,qpgr-qtarget))
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3),iq)
        geigen(ig2,:)= geigenr(ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
      enddo
      end subroutine readgeig
!ssssssssssssssssssssssssssssssssssssssssssssss
      function readcphif(q,isp) result(cphif)
      integer,intent(in):: isp
      real(8),intent(in):: q(3)
      real(8) :: qu(3)
      complex(8):: cphif(ldim2,nband)
      call readcphi(q,ldim2,isp, qu, cphif)
      quu=qu
      end
      function readcphifq() result(qu) ! I think this is not ncessary now
      real(8):: qu(3)                  ! I think qu=q now.
      qu=quu
      end function
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readcphi(q,ldim2,isp,  qu,cphif)
      use m_rotwave,only: Rotmto
      implicit none
!!-- return mto part of eigenfunction for given q(1:3) and isp
      real(8), intent(in) :: q(3)
      integer, intent(in)  :: ldim2
      integer, intent(in)  :: isp
      real(8), intent(out)  :: qu(3)
      complex(8), intent(out)  :: cphif(ldim2,nband)
      integer:: iq,iqindx,ikpisp,iqq,iorb,ibaso,ibas,k,l,
     & ini1,ini2,iend1,iend2,igg,ig,iqi,i,igxt
      real(8)   :: qrot(3) ,qout(3)
      complex(8):: phase,cphifr(ldim2,nband),phaseatom !takao 1->*->nband
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
      complex(8):: img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      if(init2) call rx( 'readcphi: modele is not initialized yet')
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
      if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ldim2,ldim2,norbtx,
     &                sum(cphifr(1:ldim2,1:nband)),sum(abs(cphifr(1:ldim2,1:nband)))
      igxt=1 !not timereversal
      call rotmto(qtt(:,iqq),ldim2,nband,
     i  norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
     i  maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl),
     i  symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lxxa,miat(:,igg),tiat(:,:,igg),igxt,nbas,
     i  cphifr,
     o  cphif)
      return
      end subroutine readcphi
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readgeig_mlw(q,ngp_in,isp, qu,geigen)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer:: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) then
        call rx( 'readgeig_mlw: modele is not initialized yet')
      endif
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        call rx( 'readgeig_mlw: ngpmx<ngp(iq)')
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geig_mlw(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeig_mlw, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine readcphi_mlw(q,ldim2_dummy,isp,  qu,cphif)
!!-- return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer:: isp,iq,iqindx,ldim2_dummy,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2,nwf)
      if(init2) then
        call rx( 'readcphi_mlw: modele is not initialized yet')
      endif
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphi(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine init_readeigen() !nband_in,mrece_in)
!!-- initialization. Save QpGpsi EVU EVD to arrays.--
      integer:: iq,is,ifiqg,nnnn,ikp,!,nsp_in, nband_in,mrece_in,
     & iopen,iclose,isx,mrecb_in,ik,ib,verbose,ifile_handle
      integer:: ifev,nband_ev, nqi_, nsp_ev ,ngpmx_ ,nqnum_
      real(8):: QpGcut_psi
      real(8),allocatable:: qtt_(:,:),qtti_(:,:)
      write(6,*) 'init_readeigen:'
      if(nsp<0.or.nsp>2) call rx( 'init_reaeigen:nsp wrong')
      write(*,*)'nqi=',nqi!,nqnum
      call init_iqindx_qtt() 
      ifiqg  = ifile_handle()
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum_ , ngpmx_, QpGcut_psi, nnnn,nqi_ ,imx 
      write(6,*)'read(ifiqg)', nqnum , ngpmx_, QpGcut_psi, nnnn,nqi
      if(nqi  /=  nqi_) call rx( 'init_readeigen:nqi/=nqi_ 11111')
      if(nqnum/=  nqnum_) call rx( 'init_readeigen:nqnum/=nqnum_ 11111')
      if(ngpmx_/=ngpmx) call rx('ngpmx error: 1111111 readeigen')
      allocate( qtt_(3,nqnum),ngp(nqnum) )
      allocate( ngvecp(3,ngpmx,nqnum))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
      do ikp = 1,nqnum
        read (ifiqg) qtt_(:,ikp), ngp(ikp)
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp)
      enddo
      close(ifiqg)
      deallocate(qtt_)
      ifev= iopen('EValue',0,-1,0)    !feb2012takao
      read(ifev) nband_ev, nqi_, nsp_ev
      write(6,*)'read EValue: nband_ev,nqi,nsp_ev', nband_ev, nqi_, nsp_ev
      write(6,*)'             nband,   nqi,nsp   ', nband, nqi, nsp
      if(nband_ev/=nband) call rx( 'init_readeigen:nband_ev/=nband')
      if(nsp_ev  /=  nsp) call rx( 'init_readeigen:nsp_ev/=nsp')
      if(nqi  /=  nqi_)   call rx( 'init_readeigen:nqi/=nqi_')
      nqixx=nqi
      allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
      read(ifev) qtti_(1:3,1:nqi_) !qtti
      read(ifev) evud(1:nband, 1:nqi, 1:nsp )
      close(ifev)
      if(debug) then
        do is= 1,nsp
          do ik= 1,nqi
            do ib= 1,nband
              write(6,"('ib ik e=',3i5,f13.5,2x,3f9.4)") ib,ik,is,evud(ib,ik,is), qtti_(1:3,ik)
            enddo
          enddo
        enddo
        if(debug) write(6,*)'init_readeigen:end'
      endif
      leval= minval(evud)
      init=.false.
      end subroutine init_readeigen
!ssssssssssssssssssssssssssssssssssssssssssssssssss
      real(8) function lowesteval()
      lowesteval=leval
      end function
!ssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine init_readeigen2() 
! this should be called after init_readgeigen
      implicit none
      integer:: iq,is,ifiqg,nnnn,ikp, !ldim2_in, 
     & iopen,iclose,isx,ikpisp,verbose,ifoc, !mrecg_in,mrecb_in,
     & i1,i2,i3,i4,i5,iorb,iorbold
      logical :: keepeigen
      keepeig = keepeigen()
      init2=.false.
!! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      read(ifoc,*)
      norbtx=0
      do
        read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
        if(iorb>norbtx) norbtx=iorb
      enddo
 106  continue
      write(6,*) ' end of readin @MNLA_CPHI: norbtx=',norbtx
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
      ifoc=iclose('@MNLA_CPHI')
      if(Keepeig     ) write(6,*)' KeepEigen=T; readin geig and cphi into m_readeigen'
      if(.not.Keepeig) write(6,*)' KeepEigen=F; not keep geig and cphi in m_readeigen'
      ifgeig = iopen('GEIG',0,0,mrecg)
      ifcphi = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction
      if(.not.keepeig) return
      allocate(geig(ngpmx,nband,nqi,nsp))
      allocate(cphi(ldim2,nband,nqi,nsp))
      write(6,*)' size geig=',ngpmx,nband,nqi,nsp,ldim2,size(geig),size(cphi)
      do ikp= 1,nqi
        do is= 1,nsp
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
          read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nband,ikp,is)
        enddo
      enddo
      write(6,*)'eeeee end of init_readeigen2 eeeee for keepeig=T'
      isx= iclose('CPHI')
      isx= iclose('GEIG')
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine init_readeigen_mlw()!nwf_o)
c allocate cphi_mlw and geig_mlw
c for hx0fp0_mlw (ixc=112 mode)
c this should be called after init_readgeigen2
      implicit none
      integer:: iq,is,ifiqg,ikp,
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
      write(6,*)' init_readeigen_mlw'
c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp
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
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo
        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix)  call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx)  call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo
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
      nwf_o = nwf
      end subroutine
!ssssssssssssssssssssssssssssssssssssssssssssss
      subroutine init_readeigen_mlw_noeval()!nwf_o,nband_o,mrecb_o,mrecg_o)
c      use m_hamindex,only: nqi,nqnum,qtti
c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer:: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm,qu(3),tolq=1d-8
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:),cphi2(:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)
      logical :: keepeigen
      integer:: ikpx,ifi,ifile_handle
      character*(8):: fname
      keepeig = keepeigen()
      write(6,*)' init_readeigen_mlw_noeval'
c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp
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
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1)  allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix) call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx) call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
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
      allocate(evud(nwf,nqi,nsp)) !nqtt
      evud = 0d0

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp)) 
      cbwf = 0d0
!!
      if(Wpkm4crpa) then
         fname='pkm4crpa'
         ifi=ifile_handle()
         open(ifi,file=fname,form='formatted',status='unknown')
         write(ifi,"('== p_km^alpha in PRB83,121101 ! weight in l-subspace ==')")
         write(ifi,"('( = c^sigma_km in book of 45th IFFK by Ersoy)')")
         write(ifi,"(8i8)") nqtt,nwf,nsp,iko_ix,iko_fx
         write(ifi,"('       |pkm|**2          ib      iq     is       q(1:3)')")
      endif

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
!! --- write pkm4crpa
          if(Wpkm4crpa) then
             do ib = iko_ix,iko_fx
                write(ifi,"(f19.15, 3i8, 3f13.6 )")
     &            sum(abs(cbwf(ib,1:nwf,ikp,is))**2),ib, ikp, is, qtt(1:3,ikp)
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
      if(Wpkm4crpa) then
         close(ifi)
      endif

      deallocate(dnk,uum)
      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband

      if(keepeig) then
        write(6,*)' xxx nband=',nband
        allocate(geig2(ngpmx,nband))  !nqtt -->nqi
        allocate(cphi2(ldim2,nband))
        allocate(geigW(ngpmx,nwf,nqtt,nsp))
        allocate(cphiW(ldim2,nwf,nqtt,nsp))
        geigW = 0d0
        cphiW = 0d0
        do ikp= 1,nqtt ! nqi
        do is= 1,nsp
           if(debug) write(6,"(' ikp=',i5,3f10.5)") ikp,qtt(:,ikp)
           call readgeig(qtt(:,ikp),ngpmx,is, qu,geig2)
           if(debug)print *,'qqqqqq1',qu
           if(debug)print *,'qqqqqq2',qtt(:,ikp)
           if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 1111')
           call readcphi(qtt(:,ikp),ldim2,is, qu,cphi2)
           if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 2222')
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geigW(:,iwf,ikp,is) = geigW(:,iwf,ikp,is) + geig2(:,ib)*cbwf(ib,iwf,ikp,is)
                cphiW(:,iwf,ikp,is) = cphiW(:,iwf,ikp,is) + cphi2(:,ib)*cbwf(ib,iwf,ikp,is)
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
        deallocate(geig2,cphi2,geig,cphi)
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
        ifgeigW = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphiW = iopen('CPHI.mlw',0,0,mrecb_o)
      endif
      deallocate(cbwf)
      end subroutine init_readeigen_mlw_noeval
c-------------------------------------------------------------
      subroutine init_readeigen_phi_noeval()!nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi and geig
c for hwmat_phi
c this should be called after init_readgeigen2
      implicit none
      integer:: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbz2,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,nsp2,nlmto2,ngpmx2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,ifdim
      real(8):: q(3),rnorm,cnorm
      complex(8),allocatable :: geig2(:,:),cphi2(:,:)
      logical :: keepeigen
      keepeig = keepeigen()
      write(6,*)' init_readeigen_phi_noeval'
      ifdim = iopen('PHIG.d',1,0,0)
      read(ifdim,*) nsp2,nqbz2,nwf2,nlmto2,ngpmx2
      isx= iclose('PHIG.d')
      if(nsp2.ne.nsp) then
        write(6,*)'nsp,nsp2',nsp,nsp2
        call rx( 'init_readeigen_phi: ns')
      endif
      if(nlmto2.ne.ldim2) then
        write(6,*)'nlmto,nlmto2',ldim2,nlmto2
        call rx( 'init_readeigen_phi: nlmto')
      endif
      if(ngpmx2.ne.ngpmx) then
        write(6,*)'ngpmx,ngpmx2',ngpmx,ngpmx2
        call rx( 'init_readeigen_phi: ngpmx')
      endif
      nwf = nwf2
      nqbz=nqbz2
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0
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
          do is= 1,nsp
            ikpisp = is + nsp*(iqbz-1)
            read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nwf,ikp,is)
            read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nwf,ikp,is)
          enddo
        enddo
      else
      endif
      end subroutine

c-------------------------------------------------------------
      subroutine init_readeigen_mlw_diag()!nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi, geig and evud
c for hx0fp0_mlw, ixc=11 mode
c this should be called after init_readgeigen2
      implicit none
      integer:: iq,is,ifiqg,ikp,
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
      write(6,*)' init_readeigen_mlw'
c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp
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
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp), dnk2(iko_ix:iko_fx,nwf))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix) call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx) call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo
c nqbze diag-points
        read(ifmlw)nqbze,nwf2
        if (nwf2 .ne. nwf) call rx( 'init_readeigen_mlw: nwf error')
        if (nqbze.ne.nqtt .and. nqbze.ne.nqtt-1)
     &      call rx( 'init_readeigen_mlw: nqbze error')
        read(ifmlwe)nqbze2,nwf2
        if (nwf2 .ne. nwf) call rx( 'init_readeigen_mlw: nwf error')
        if (nqbze2.ne.nqbze) call rx( 'init_readeigen_mlw: nqbze2 error')
        if (is.eq.1)allocate(evec(nwf,nwf,nqtt,nsp),eval(nwf,nqtt,nsp))
        evec = 0d0
        eval = 0d0
        do iqbz = 1,nqbze
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)evec(1:nwf,1:nwf,iqbz,is)
          read(ifmlwe)iqbz2,q(1:3)
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen: iqbz error')
          read(ifmlwe)eval(1:nwf,iqbz,is)
        enddo
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif
      enddo
c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = eval
      deallocate(eval)
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
      end subroutine
      end module m_readeigen

      
      
