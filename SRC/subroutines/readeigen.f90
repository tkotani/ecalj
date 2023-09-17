!! -----------------------------------------------------------
!> Return eigenvalus and eigenfunctions for given q and isp.
!! We can get eigenfunctions for Wannier, as well. See hmagnon.F
!! note: we have to call init_foobar to call readeval, readcphi, readgeig.
!!----------------
module m_readeigen
  use m_iqindx_qtt,only: Iqindx2_, Init_iqindx_qtt
  use m_hamindex,only:   ngpmx, nqtt, nqi, qtt,iqimap, iqmap,igmap,shtvg,qlat,symops
  use m_hamindex,only:   plat,invgx, miat,tiat,dlmm,shtvg,symops,lmxax,nbas
  use m_read_bzdata,only: ginv
  use m_genallcf_v3,only: nsp =>nspin ,ldim2=>nlmto
  use m_readhbe,only: mrecb, mrece,nband, mrecg
  !! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
  !! qtti(1:3,nqi)   :eivenvalues, eigenvectors are calculated only for them.
  !!                  See lmfgw (q-vector with irr flag in QGpsi (output of qg4gw).
  implicit none
  !-------------------------------
  public::Readeval, Readgeigf, Readcphif, Readcphifq, Init_readeigen, Init_readeigen2, Lowesteval
  public::Onoff_write_pkm4crpa, Init_readeigen_mlw_noeval, Init_readeigen_phi_noeval, Readcphiw, Readgeigw
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
  integer,private:: nprecb, nnnn, ifcphi,ifgeig
  complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
  real(8),private:: leval
  integer,private :: ifcphi_mlw,ifgeig_mlw, nqixx
  complex(8),allocatable,private:: geig_mlw(:,:,:,:), cphi_mlw(:,:,:,:)
  logical,private:: debug=.false.
  integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
  integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:),offset_tbl(:),offset_rev_tbl(:,:,:)
  logical,private:: Wpkm4crpa=.false.
  real(8),private :: quu(3)

contains
  subroutine onoff_write_pkm4crpa(lll)
    logical:: lll
    Wpkm4crpa=lll
  end subroutine onoff_write_pkm4crpa
  ! sssssssssssssssssssssssssssssssssssssssssssss
  !> Return ev(1:nband) for given q(1:3) and isp
  pure function readeval(q,isp) result(ev)
    intent(in)  ::       q,isp
    integer :: isp
    real(8) :: q(3), ev(nband)
    integer:: iq,iqindx,i
    real(8):: qu(3)
!    if(init) call rx( 'readeigen: modele is not initialized yet')
    !      iq = iqindx(q, ginv,qtt,nqtt)
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
!    if(debug) then
!       write(6,*)'iq iqimap(iq)=',iq,iqimap(iq)
!       write(6,"('iq iqimap(iq) q=',2i8,3f13.5)")iq,iqimap(iq),q
!    endif
    ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
!    if(debug) then
!       write(6,"(9f9.4)")ev(1:9)
!    endif
  end function readeval
  !> Return ev(1:nband) for given q(1:3) and isp
  function readgeigf(q,isp) result(geigen)
    real(8),intent(in):: q(3)
    integer,intent(in):: isp
    real(8):: qu(3)
    complex(8):: geigen(ngpmx,nband)
    call readgeig(q,ngpmx,isp,qu,geigen)
  end function readgeigf
  ! sssssssssssssssssssssssssssssssssssssssssssss
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
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
    if(debug) write(6,*)' readgeig:xxx iq=',iq
    iqq=iqmap(iq)
    iqi=iqimap(iq)
    igg=igmap(iq)
    qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
    !!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
    if(ngp(iq)==0) return
    if(ngp(iq)/=ngp(iqq)) then
       write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg,q,qu
       write(6,*)' ddddd qtarget=',qtarget,' ddddd q (iqq)=',qtt(:,iqq)
       write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
       call rx( 'readgeig:x2 ngp(iq)/=ngp(iqq)')
    endif
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
    !!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
    igxt=1 !not timereversal
    call rotipw(qtt(:,iqq),qtt(:,iq),ngp(iqq),nband, &
         platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx, &
         geigenr(1:ngp(iqq),1:nband), geigen(1:ngp(iq),1:nband))
  end subroutine readgeig
  ! sssssssssssssssssssssssssssssssssssssssssssss
  function readcphif(q,isp) result(cphif)
    integer,intent(in):: isp
    real(8),intent(in):: q(3)
    real(8) :: qu(3)
    complex(8):: cphif(ldim2,nband)
    call readcphi(q,ldim2,isp, qu, cphif)
    quu=qu
  end function readcphif
  ! sssssssssssssssssssssssssssssssssssssssssssss
  function readcphifq() result(qu)
    real(8):: qu(3)                  ! I think qu=q now.
    qu=quu
  end function readcphifq
  ! sssssssssssssssssssssssssssssssssssssssssssss
  subroutine readcphi(q,ldim2,isp,  qu,cphif)
    use m_rotwave,only: Rotmto
    implicit none
    !!-- return mto part of eigenfunction for given q(1:3) and isp
    real(8), intent(in) :: q(3)
    integer, intent(in)  :: ldim2
    integer, intent(in)  :: isp
    real(8), intent(out)  :: qu(3)
    complex(8), intent(out)  :: cphif(ldim2,nband)
    integer:: iq,iqindx,ikpisp,iqq,iorb,ibaso,ibas,k,l,ini1,ini2,iend1,iend2,igg,ig,iqi,i,igxt
    real(8)   :: qrot(3) ,qout(3)
    complex(8):: phase,cphifr(ldim2,nband),phaseatom !takao 1->*->nband
    complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
    complex(8):: img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR
    if(init2) call rx( 'readcphi: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !for given q, get iq. qu is used q. q-qu= G vectors. qu=qtt(:,iq)
    igg=igmap(iq)  ! qtt(:,iq)= matmul(sympos(  ,igg),qtt(:,iqq))
    iqq=iqmap(iq)  ! mapped from qtt(:,iqq) to qtt(:,iq);
    ! qtt(:,iq)=matmul(sym(igg),qtt(:,iqq))+some Gvector(see iqindx2 above)
    iqi=iqimap(iq) ! iqi is index for irr.=1 (cphi calculated. See qg4gw and sugw.F).
    ! qtt(:,iqq) = qtti(:,iqi) is satisfied.
    ! we have eigenfunctions calculated only for qtti(:,iqi).
    if(keepeig) then
       cphifr(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iqi,isp)
    else
       ikpisp= isp + nsp*(iqi-1)
       read(ifcphi, rec=ikpisp) cphifr(1:ldim2,1:nband)
    endif
    if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ldim2,ldim2,norbtx, &
         sum(cphifr(1:ldim2,1:nband)),sum(abs(cphifr(1:ldim2,1:nband)))
    igxt=1 !not timereversal (for future)
    call rotmto(qtt(:,iqq),ldim2,nband,norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl, &
         maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl), &
         symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lmxax,miat(:,igg),tiat(:,:,igg),igxt,nbas, &
         cphifr, cphif)
  end subroutine readcphi
  ! sssssssssssssssssssssssssssssssssssssssssssss
  subroutine init_readeigen() !nband_in,mrece_in)
    !!-- initialization. Save QpGpsi EVU EVD to arrays.--
    integer:: iq,is,ifiqg,nnnn,ikp,isx,mrecb_in,ik,ib,verbose
    integer:: ifev,nband_ev, nqi_, nsp_ev ,ngpmx_ ,nqtt_
    real(8):: QpGcut_psi
    real(8),allocatable:: qtt_(:,:),qtti_(:,:)
    write(6,*) 'init_readeigen:'
    if(nsp<0 .OR. nsp>2) call rx( 'init_reaeigen:nsp wrong')
    !write(*,*)'nqi=',nqi!,nqtt
    call init_iqindx_qtt()
    open(newunit=ifiqg ,file='QGpsi',form='unformatted')
    read(ifiqg) nqtt_ , ngpmx_, QpGcut_psi, nnnn,nqi_ ,imx
    write(6,*)'read(ifiqg)', nqtt , ngpmx_, QpGcut_psi, nnnn,nqi
    if(nqi  /=  nqi_) call rx( 'init_readeigen:nqi/=nqi_ 11111')
    if(nqtt/=  nqtt_) call rx( 'init_readeigen:nqtt/=nqtt_ 11111')
    if(ngpmx_/=ngpmx) call rx('ngpmx error: 1111111 readeigen')
    allocate( qtt_(3,nqtt),ngp(nqtt) )
    allocate( ngvecp(3,ngpmx,nqtt))
    allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqtt) )
    do ikp = 1,nqtt
       read (ifiqg) qtt_(:,ikp), ngp(ikp)
       read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp)
    enddo
    close(ifiqg)
    deallocate(qtt_)
    open(newunit=ifev,file='EValue',form='unformatted')
    read(ifev) nband_ev, nqi_, nsp_ev
    write(6,*)'read EValue: nband_ev,nqi,nsp_ev', nband_ev, nqi_, nsp_ev
    write(6,*)'             nband,   nqi,nsp   ', nband, nqi, nsp
    if(nband_ev/=nband) call rx( 'init_readeigen:nband_ev/=nband')
    if(nsp_ev  /=  nsp) call rx( 'init_readeigen:nsp_ev/=nsp')
    if(nqi  /=  nqi_)   call rx( 'init_readeigen:nqi/=nqi_')
    nqixx=nqi
    allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
    read(ifev) qtti_(1:3,1:nqi_)
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
  ! sssssssssssssssssssssssssssssssssssssssssssss
  real(8) function lowesteval()
    lowesteval=leval
  end function lowesteval
  ! sssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine init_readeigen2()    ! this should be called after init_readgeigen
    implicit none
    integer:: iq,is,ifiqg,nnnn,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold
    logical :: keepeigen
    call readmnla_cphi()
    keepeig = keepeigen()
    init2=.false.
    if(Keepeig     ) write(6,*)' KeepEigen=T; readin geig and cphi into m_readeigen'
    if( .NOT. Keepeig) write(6,*)' KeepEigen=F; not keep geig and cphi in m_readeigen'
    open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
    open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
    if( .NOT. keepeig) return
    allocate(geig(ngpmx,nband,nqi,nsp))
    allocate(cphi(ldim2,nband,nqi,nsp))
    !      write(6,*)' size geig=',ngpmx,nband,nqi,nsp,ldim2,size(geig),size(cphi)
    do ikp= 1,nqi
       do is= 1,nsp
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
          read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nband,ikp,is)
       enddo
    enddo
    close(ifgeig)
    close(ifcphi)
  end subroutine init_readeigen2
  subroutine readmnla_cphi()
    !! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
    implicit none
    integer:: iq,is,ifiqg,nnnn,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold
    open(newunit=ifoc,file='@MNLA_CPHI')
    read(ifoc,*)
    norbtx=0
    do
       read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
       if(iorb>norbtx) norbtx=iorb
    enddo
106 continue
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
107 continue
    allocate(offset_rev_tbl(maxval(ibas_tbl),0:maxval(l_tbl),maxval(k_tbl)))
    offset_rev_tbl=-99999
    do iorb=1,norbtx
       offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))= offset_tbl(iorb)
    enddo
    close(ifoc)
  end subroutine readmnla_cphi
  ! sssssssssssssssssssssssssssssssssssssssssssss
  subroutine init_readeigen_mlw_noeval()
    ! replace cphi and geig
    ! for hwmat
    ! this should be called after init_readgeigen2
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,mrecb_o,ikpisp,mrecg_o, &
         nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2, &
         ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o, &
         ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
    real(8):: q(3),rnorm,cnorm,qu(3),tolq=1d-8
    real(8),allocatable :: eval(:,:,:)
    complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:), &
         geig2(:,:),cphi2(:,:), &
         geig3(:,:),cphi3(:,:), &
         geig4(:,:),cphi4(:,:), &
         cbwf(:,:,:,:),uum(:,:,:,:,:)
    logical :: keepeigen
    integer:: ikpx,ifi
    character*(8):: fname
    keepeig = keepeigen()
    write(6,*)' init_readeigen_mlw_noeval'
    ! --- Readin MLWU/D, MLWEU/D, and UUq0U/D
    do is = 1,nsp
       if (is == 1) then
          open(newunit=ifmlw,file='MLWU',form='unformatted')
          open(newunit=ifmlwe,file='MLWEU', form='unformatted')
          open(newunit=ifuu,file='UUq0U', form='unformatted')
       else
          open(newunit=ifmlw  ,file='MLWD', form='unformatted')
          open(newunit=ifmlwe ,file='MLWED',form='unformatted')
          open(newunit=ifuu   ,file='UUq0D',form='unformatted')
       endif
       ! nqbz mesh-points
       read(ifmlw)nqbz,nwf,iko_ix,iko_fx
       if (is == 1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
       do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
          if (iqbz2 /= iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
       enddo
       read(ifuu)
       read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
       if (is == 1)  allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
       if (nqbz2 /= nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
       if (iko_ix2 /= iko_ix) call rx( "init_readeigen_mlw: iko_ix2 error")
       if (iko_fx2 /= iko_fx) call rx( "init_readeigen_mlw: iko_fx2 error")
       do iqbz = 1,nqbz
          do iq0i =1,nq0i
             read(ifuu)
             read(ifuu)iqbz2,iq0i2
             if (iqbz2 /= iqbz) call rx( 'init_readeigen_mlw: iqbz error')
             if (iq0i2 /= iq0i) call rx( 'init_readeigen_mlw: iq0i error')
             read(ifuu)((uum(j1,j2,iqbz,iq0i,is), j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
       enddo
       if (is == 1) then
          close(ifmlw)
          close(ifmlwe)
          close(ifuu)
       else
          close(ifmlw)
          close(ifmlwe)
          close(ifuu)
       endif
    enddo
    ! replace evud
    deallocate(evud)
    allocate(evud(nwf,nqi,nsp)) !nqtt
    evud = 0d0
    ! replace geig and cphi
    allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
    cbwf = 0d0
    if(Wpkm4crpa) then
       fname='pkm4crpa'
       open(newunit=ifi,file=fname,form='formatted',status='unknown')
       write(ifi,"('== p_km^alpha in PRB83,121101 ! weight in l-subspace ==')")
       write(ifi,"('( = c^sigma_km in book of 45th IFFK by Ersoy)')")
       write(ifi,"(8i8)") nqtt,nwf,nsp,iko_ix,iko_fx
       write(ifi,"('       |pkm|**2          ib      iq     is       q(1:3)')")
    endif
    do ikp = 1,nqtt
       iqbz = mod(ikp,nqbz)
       if (iqbz == 0) iqbz = nqbz
       iq0i = (ikp - iqbz)/nqbz
       do is= 1,nsp
          if (iq0i == 0) then
             do ib = iko_ix,iko_fx
                do iwf= 1,nwf
                   cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
                enddo
             enddo
          else
             !   <psi(k+q0,n) | psi(k+q0,m)^B>
             ! = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
             !      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
             ! ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>

             ! psi^B : bloch fn. corresponding to maxloc Wannier fn.
             do ib = iko_ix,iko_fx
                do iwf= 1,nwf
                   cbwf(ib,iwf,ikp,is) = sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is)) &
                        *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
                enddo
             enddo
          endif
          !! --- write pkm4crpa
          if(Wpkm4crpa) then
             do ib = iko_ix,iko_fx
                write(ifi,"(f19.15, 3i8, 3f13.6 )") &
                     sum(abs(cbwf(ib,1:nwf,ikp,is))**2),ib, ikp, is, qtt(1:3,ikp)
             enddo
          endif
          ! m norm check
          !         do iwf  = 1,nwf
          !         do iwf2 = 1,nwf
          !           rnorm = 0d0
          !           cnorm = 0d0
          !           do ib = iko_ix,iko_fx
          !              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
          !     &                                     *cbwf(ib,iwf2,ikp,is))
          !              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
          !     &                                     *cbwf(ib,iwf2,ikp,is))
          !              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
          !     &                                    *dnk(ib,iwf2,iqbz,is))
          !              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
          !     &                                    *dnk(ib,iwf2,iqbz,is))
          !           enddo
          !           do ib = 1,nwf
          !              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
          !     &                                    *evec(ib,iwf2,ikp,is))
          !              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
          !     &                                    *evec(ib,iwf2,ikp,is))
          !           enddo
          !           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
          !           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
          !         enddo
          !         enddo
          !         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
          !         write(7300,*)cbwf(:,:,ikp,is)
       enddo
    enddo
    if(Wpkm4crpa) close(ifi)
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
             ! eck write
             !            do iwf  = 1,nwf
             !            do iwf2 = 1,nwf
             !               rnorm = 0d0
             !               cnorm = 0d0
             !               do ib = 1,ldim2
             !                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
             !     &                                   cphi(ib,iwf2,ikp,is))
             !               enddo
             !               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
             !               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
             !            enddo
             !            enddo
             !            write(7500,*)ikp,ldim2,nwf
             !            write(7500,*)cphi(:,:,ikp,is)
          enddo
       enddo
       deallocate(geig2,cphi2,geig,cphi)
    else
       open(newunit=ifcphi_o,file='CPHI.mlw',form='unformatted')
       open(newunit=ifgeig_o,file='GEIG.mlw',form='unformatted')
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
       close(ifcphi)
       close(ifgeig)
       close(ifcphi_o)
       close(ifgeig_o)
       open(newunit=ifgeigW,file='GEIG.mlw',form='unformatted')
       open(newunit=ifcphiW,file='CPHI.mlw',form='unformatted')
    endif
    deallocate(cbwf)
  end subroutine init_readeigen_mlw_noeval
  !-------------------------------------------------------------
  subroutine init_readeigen_phi_noeval()
    ! replace cphi and geig
    ! for hwmat_phi
    ! this should be called after init_readgeigen2
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,mrecb_o,ikpisp,mrecg_o, &
         nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbz2,nqbze,nqbze2, &
         iqbz,iqbz2,nwf2,nsp2,nlmto2,ngpmx2, &
         ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,ifdim
    real(8):: q(3),rnorm,cnorm
    complex(8),allocatable :: geig2(:,:),cphi2(:,:)
    logical :: keepeigen
    keepeig = keepeigen()
    !      write(6,*)' init_readeigen_phi_noeval'
    open(newunit=ifdim,file='PHIG.d')
    read(ifdim,*) nsp2,nqbz2,nwf2,nlmto2,ngpmx2
    close(ifdim)
    if(nsp2 /= nsp) then
       write(6,*)'nsp,nsp2',nsp,nsp2
       call rx( 'init_readeigen_phi: ns')
    endif
    if(nlmto2 /= ldim2) then
       write(6,*)'nlmto,nlmto2',ldim2,nlmto2
       call rx( 'init_readeigen_phi: nlmto')
    endif
    if(ngpmx2 /= ngpmx) then
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
    close(ifcphi)
    close(ifgeig)
    open(newunit=ifgeig,file='GEIGg',form='unformatted',access='direct',recl=mrecg_o)
    open(newunit=ifcphi,file='CPHIg',form='unformatted',access='direct',recl=mrecb_o)
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
  end subroutine init_readeigen_phi_noeval
  ! sssssssssssssssssssssssssssssssssssssssssssss
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
  end subroutine readgeigW
  ! sssssssssssssssssssssssssssssssssssssssssssss
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
  end subroutine readcphiW
end module m_readeigen

