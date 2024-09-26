!> Return eigenvalus and eigenfunctions for given q and isp.
! -----------------------------------------------------------
! We can get eigenfunctions for Wannier, as well. See hmagnon.F
! note: we have to call init_foobar to call readeval, readcphi, readgeig.
! ----------------
module m_readeigen
  use m_mpiio,only: openm,readm,closem
  use m_ftox
  use m_lgunit,only:stdo
  use m_iqindx_qtt,only: Iqindx2_, Init_iqindx_qtt
  use m_hamindex,only:   ngpmx, nqtt, nqi, qtt,iqimap, iqmap,igmap,shtvg,qlat,symops
  use m_hamindex,only:   plat,invgx, miat,tiat,dlmm,shtvg,symops,lmxax,nbas
  use m_read_bzdata,only: ginv
  use m_genallcf_v3,only: nsp =>nspin ,ndima,ndimanspc, mrecb,mrece,mrecg,nband,nspc,nspx
  !! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
  !! qtti(1:3,nqi)   :eivenvalues, eigenvectors are calculated only for irr=1 in QGpsi (See lqg4gw).
  implicit none
  public:: Init_readeigen,Init_readeigen2, Lowesteval, Readeval,Readgeigf,Readcphif
  public:: Onoff_write_pkm4crpa,Readcphifq
  public:: Init_readeigen_mlw_noeval, Readcphiw, Readgeigw
  integer,public:: nwf
  private
  integer:: norbtx,imx,ifcphim,ifgeigm,nqixx  !ifgeigW,ifcphiW,
  real(8),private:: leval, quu(3)
  logical,private:: init=.true.,init2=.true.,keepeig, Wpkm4crpa=.false.
  logical,private:: debug=.false.
  character(8),external :: xt
  real(8),allocatable,private:: evud(:,:,:)
  complex(8),allocatable:: geigW(:,:,:,:),cphiW(:,:,:,:)
  complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
  integer,allocatable,private:: ngp(:),ngvecp(:,:,:), ngvecprev(:,:,:,:)
  integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:),offset_tbl(:),offset_rev_tbl(:,:,:)
  complex(8),allocatable,private:: geig_mlw(:,:,:,:), cphi_mlw(:,:,:,:)
contains
  subroutine onoff_write_pkm4crpa(lll)
    logical:: lll
    Wpkm4crpa=lll
  end subroutine onoff_write_pkm4crpa
  function readcphifq() result(qu)
    real(8):: qu(3)                  ! I think qu=q now.
    qu=quu
  end function readcphifq
  pure function readeval(q,isp) result(ev) ! Return ev(1:nband) for given q(1:3) and isp
    intent(in)  ::       q,isp
    integer :: isp
    real(8) :: q(3), ev(nband)
    integer:: iq,iqindx,i
    real(8):: qu(3)
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
    ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
    !if(debug) then
    !   write(6,*)'iq iqimap(iq)=',iq,iqimap(iq)
    !   write(6,"('iq iqimap(iq) q=',2i8,3f13.5)")iq,iqimap(iq),q
    !   write(6,"(9f9.4)")ev(1:9)
    !endif
  end function readeval
  !> Return ev(1:nband) for given q(1:3) and isp
  function readgeigf(q,isp) result(geigen)
    real(8),intent(in):: q(3)
    integer,intent(in):: isp
    real(8):: qu(3)
    complex(8):: geigen(ngpmx*nspc,nband)
    geigen=0d0 !2024-5-17 for ifort NaN initialization
    call readgeig(q,isp,qu,geigen)
  end function readgeigf
  subroutine readgeig(q,isp, qu,geigen)
    use m_ftox
    use m_rotwave,only: Rotipw
    use m_ftox
    implicit none
    !logical,optional,intent(in):: fpmt
    real(8),intent(in) :: q(3)
    integer,intent(in) :: isp
    real(8),intent(out) :: qu(3)
    complex(8), intent(out) :: geigen(ngpmx*nspc,nband)
    integer:: iq,iqindx,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi,igxt,i,ioff,ispc
    real(8)   :: ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3),qout(3),qin(3)
    complex(8):: geigenr(ngpmx*nspc,nband),img=(0d0,1d0),img2pi
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
       geigenr(1:ngpmx*nspc,1:nband) = geig(1:ngpmx*nspc,1:nband,iqi,isp)
    else
       ikpisp= isp + nsp*(iqi-1)
       i=readm(ifgeigm,rec=ikpisp, data=geigenr(1:ngpmx*nspc,1:nband) )
!       open(newunit=ifgeig, file='GEIG'//trim(xt(iqi))//trim(xt(isp)),form='unformatted')
!       read(ifgeig) geigenr(1:ngpmx*nspc,1:nband) 
!       close(ifgeig)
    endif
    !!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
    igxt=1 !not timereversal
    do ispc=1,nspc
       ioff=(ispc-1)*ngpmx
       call rotipw(qtt(:,iqq),qtt(:,iq),ngp(iqq),nband, &
            platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx, &
            geigenr(ioff+1:ioff+ngp(iqq),1:nband), geigen(ioff+1:ioff+ngp(iq),1:nband))
    enddo
  end subroutine readgeig
  function readcphif(q,isp) result(cphif)
    integer,intent(in):: isp
    real(8),intent(in):: q(3)
    real(8) :: qu(3)
    complex(8):: cphif(ndima*nspc,nband)
    call readcphi(q,isp, qu, cphif)
    quu=qu
  end function readcphif
  subroutine readcphi(q,isp,  qu,cphif)!, fpmt)
    use m_rotwave,only: Rotmto
    implicit none
    !logical,optional,intent(in):: fpmt
    !!-- return mto part of eigenfunction for given q(1:3) and isp
    real(8), intent(in) :: q(3)
    integer, intent(in)  :: isp
    real(8), intent(out)  :: qu(3)
    complex(8), intent(out)  :: cphif(ndima*nspc,nband)
    integer:: iq,iqindx,ikpisp,iqq,iorb,ibaso,ibas,k,l,ini1,ini2,iend1,iend2, igg,ig,iqi,i,igxt,ioff,ispc,ix
    real(8)   :: qrot(3) ,qout(3)
    complex(8):: phase,cphifr(ndima*nspc,nband),phaseatom !takao 1->*->nband
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
       cphifr(1:ndima*nspc,1:nband) = cphi(1:ndima*nspc,1:nband,iqi,isp)
    else 
       ikpisp= isp + nsp*(iqi-1)
       i=readm(ifcphim,rec=ikpisp, data=cphifr(1:ndima*nspc,1:nband)) ! , rec=ikpisp
!     open(newunit=ifcphi, file='CPHI'//trim(xt(iqi))//trim(xt(isp)),form='unformatted')
!       read(ifcphi) cphifr(1:ndima*nspc,1:nband) 
!       close(ifcphi)
    endif
    if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ndimanspc,ndimanspc,norbtx, &
         sum(cphifr(1:ndimanspc,1:nband)),sum(abs(cphifr(1:ndimanspc,1:nband)))
    igxt=1 !not timereversal (for future)
    do ispc=1,nspc
       ioff=ndima*(ispc-1)
       call rotmto(qtt(:,iqq),ndima,nband,norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl, &
            maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl), &
            symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lmxax,miat(:,igg),tiat(:,:,igg),igxt,nbas, &
            cphifr(ioff+1:ioff+ndima,:), cphif(ioff+1:ioff+ndima,:))
     enddo
  end subroutine readcphi
  subroutine init_readeigen() ! initialization. Save QpGpsi EVU EVD to arrays.--
    integer:: iq,is,ifiqg,nnnn,ikp,isx,ik,ib,verbose
    integer:: ifev,nband_ev, nqi_, nsp_ev ,ngpmx_ ,nqtt_,nspc_
    real(8):: QpGcut_psi
    real(8),allocatable:: qtt_(:,:),qtti_(:,:)
    write(stdo,ftox) 'init_readeigen:'
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
    read(ifev) nband_ev, nqi_, nsp_ev, nspc_
    write(stdo,ftox)'Read EValue: nband nqi nsp nspc nspx', nband, nqi, nsp,nspc,nspx
    if(nband_ev/=nband) call rx( 'init_readeigen:nband_ev/=nband')
    if(nsp_ev  /= nsp)  call rx( 'init_readeigen:nsp_ev/=nsp')
    if(nqi     /= nqi_) call rx( 'init_readeigen:nqi/=nqi_')
    if(nspc    /= nspc_)call rx( 'init_readeigen:nspc/=nspc_')
    nqixx=nqi
    allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
    read(ifev) qtti_(1:3,1:nqi_)
    read(ifev) evud(1:nband, 1:nqi, 1:nspx )
    close(ifev)
    if(debug) then
       do is= 1,nspx
          do ik= 1,nqi
             do ib= 1,nband
               if(evud(ib,ik,is)<1d10) & !Set huge number for padding in sugw.f90
                    write(6,"('ib ik e=',3i5,f13.5,2x,3f9.4)") ib,ik,is,evud(ib,ik,is), qtti_(1:3,ik)
             enddo
          enddo
       enddo
       if(debug) write(6,*)'init_readeigen:end'
       call rx0('xxxxxxxxxxxxxxxxxx')
    endif
    leval= minval(evud)
    init=.false.
  end subroutine init_readeigen
  real(8) function lowesteval()
    lowesteval=leval
  end function lowesteval
  subroutine init_readeigen2()    ! this should be called after init_readgeigen
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold,i
    logical :: keepeigen,cmdopt0
    character(8) :: xt
    call readmnla_cphi()
    keepeig = keepeigen()
    init2=.false.
    if(Keepeig       ) write(6,*)' KeepEigen=T; readin geig and cphi into m_readeigen'
    if( .NOT. Keepeig) write(6,*)' KeepEigen=F; not keep geig and cphi in m_readeigen'
    if( .NOT. keepeig) return
    allocate(geig(ngpmx*nspc,nband,nqi,nspx))
    allocate(cphi(ndima*nspc,nband,nqi,nspx))
    i=openm(newunit=ifcphim,file='CPHI',recl=mrecb)
    i=openm(newunit=ifgeigm,file='GEIG',recl=mrecg)
    do ikp= 1,nqi
       do is= 1,nspx
          ikpisp= is + nsp*(ikp-1)
          i=readm(ifcphim,rec=ikpisp, data=cphi(1:ndima*nspc,1:nband,ikp,is))
          if(ngpmx/=0) i=readm(ifgeigm,rec=ikpisp,data=geig(1:ngpmx*nspc,1:nband,ikp,is))
          !  open(newunit=ifcphi, file='CPHI'//trim(xt(ikp))//trim(xt(is)),form='unformatted')
          !  open(newunit=ifgeig, file='GEIG'//trim(xt(ikp))//trim(xt(is)),form='unformatted')
          !  if(ngpmx/=0) read(ifgeig) geig(1:ngpmx*nspc,1:nband,ikp,is) ! , rec=ikpisp)
          !  close(ifcphi)
          !  close(ifgeig)
       enddo
    enddo
    if(keepeig)i=closem(ifcphim)
    if(keepeig)i=closem(ifgeigm)
  end subroutine init_readeigen2
  subroutine readmnla_cphi()
    !! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold
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
  subroutine init_readeigen_mlw_noeval() ! replace cphi and geig for hwmat ! this should be called after init_readgeigen2
  !xxxxxxxxxxxxxx only for nspc=1. Need fixing for nspc=2  
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
    keepeig = .True. !keepeigen()
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
          !              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))*cbwf(ib,iwf2,ikp,is))
          !              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))*cbwf(ib,iwf2,ikp,is))
          !              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))*dnk(ib,iwf2,iqbz,is))
          !              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))*dnk(ib,iwf2,iqbz,is))
          !           enddo
          !           do ib = 1,nwf
          !              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))*evec(ib,iwf2,ikp,is))
          !              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))*evec(ib,iwf2,ikp,is))
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
       allocate(geig2(ngpmx*nspc,nband))  !nqtt -->nqi
       allocate(cphi2(ndima*nspc,nband))
       allocate(geigW(ngpmx*nspc,nwf,nqtt,nsp))
       allocate(cphiW(ndima*nspc,nwf,nqtt,nsp))
       geigW = 0d0
       cphiW = 0d0
       do ikp= 1,nqtt ! nqi
          do is= 1,nsp
             if(debug) write(6,"(' ikp=',i5,3f10.5)") ikp,qtt(:,ikp)
             call readgeig(qtt(:,ikp),is, qu,geig2)
             if(debug)print *,'qqqqqq1',qu
             if(debug)print *,'qqqqqq2',qtt(:,ikp)
             if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 1111')
             call readcphi(qtt(:,ikp),is, qu,cphi2)
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
             !               do ib = 1,ndimanspc
             !                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*cphi(ib,iwf2,ikp,is))
             !               enddo
             !               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
             !               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
             !            enddo
             !            enddo
             !            write(7500,*)ikp,ndimanspc,nwf
             !            write(7500,*)cphi(:,:,ikp,is)
          enddo
       enddo
       deallocate(geig2,cphi2,geig,cphi)
    else
       call rx('KeepEigen=F not implemented')
       ! open(newunit=ifcphi_o,file='CPHI.mlw',form='unformatted')
       ! open(newunit=ifgeig_o,file='GEIG.mlw',form='unformatted')
       ! allocate(geig3(ngpmx,nwf))
       ! allocate(cphi3(ndimanspc,nwf))
       ! allocate(geig4(ngpmx,nband))
       ! allocate(cphi4(ndimanspc,nband))
       ! do ikp= 1,nqtt
       !    do is= 1,nsp
       !       ikpisp= is + nsp*(ikp-1)
       !       read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
       !       read(ifcphi, rec=ikpisp) cphi4(1:ndimanspc,1:nband)
       !       geig3 = 0d0
       !       cphi3 = 0d0
       !       do iwf= 1,nwf
       !          do ib= iko_ix,iko_fx
       !             geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
       !             cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
       !          enddo
       !       enddo
       !       write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
       !       write(ifcphi_o, rec=ikpisp) cphi3(1:ndimanspc,1:nwf)
       !    enddo
       ! enddo
       ! deallocate(geig3,geig4,cphi3,cphi4)
       ! close(ifcphi)
       ! close(ifgeig)
       ! close(ifcphi_o)
       ! close(ifgeig_o)
       ! open(newunit=ifgeigW,file='GEIG.mlw',form='unformatted')
       ! open(newunit=ifcphiW,file='CPHI.mlw',form='unformatted')
    endif
    deallocate(cbwf)
  end subroutine init_readeigen_mlw_noeval
  subroutine readgeigW(q,ngp_in,isp, qu,geigen)
   integer:: isp,iq,iqindx,ngp_in,ikpisp
   real(8)   :: q(3),qu(3)
   complex(8):: geigen(ngp_in,nwf)
   if(init2) call rx( 'readgeig_mlw: modele is not initialized yet')
   call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
   if(ngp_in < ngp(iq)) then
      write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q,nspc
      call rx( 'readgeig_mlw: ngpmx<ngp(iq)')
   endif
!   if(keepeig) then
      geigen(1:ngp(iq),1:nwf) = geigW(1:ngp(iq),1:nwf,iq,isp)
!   else
!      ikpisp= isp + nsp*(iq-1)
!      read(ifgeigW) geigen(1:ngpmx,1:nwf)
!   endif
 end subroutine readgeigW
 subroutine readcphiW(q,ndimanspc_dummy,isp,  qu,cphif)
   integer:: isp,iq,iqindx,ndimanspc_dummy,ikpisp
   real(8)   :: q(3),qu(3)
   complex(8):: cphif(ndima*nspc,nwf)
   if(init2) call rx( 'readcphi_mlw: modele is not initialized yet')
   call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
!   if(keepeig) then
      cphif(1:ndima*nspc,1:nwf) = cphiW(1:ndima*nspc,1:nwf,iq,isp)
!   else
!      ikpisp= isp + nsp*(iq-1)
!      read(ifcphi_mlw) cphif(1:ndimanspc,1:nwf)
!   endif
 end subroutine readcphiW
end module m_readeigen
