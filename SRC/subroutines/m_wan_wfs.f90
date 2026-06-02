module m_wan_wfs
  use m_hamindex,only:   ngpmx, nqtt, nqi, qtt,iqimap, iqmap,igmap,shtvg,qlat,symops,ngrp
  use m_hamindex,only:   plat,invgx, miat,tiat,dlmm,shtvg,symops,lmxax,nbas
  use m_iqindx_qtt,  only: iqindx2_
  use m_genallcf_v3, only: ndima, nband, nspc
  use m_readeigen,   only: readgeigf => readgeigf_mpi, readcphif => readcphif_mpi, ngp
  use m_mpi, only: ipr, MPI__AllreduceAND, MPI__zBcast => MPI__zBcast_h, get_mpi_master
  use m_genallcf_v3,only: nsp =>nspin ,ndima,ndimanspc, mrecb,mrece,mrecg,nband,nspc,nspx
  use m_lgunit,only:stdo
  use m_ftox
  use m_cmdopt_registry, only: c0_mlo
  implicit none
  public ::  get_geig_wan, get_cphi_wan
  public :: Init_readeigen_mlw_noeval, Readcphiw, Readgeigw
  public:: Onoff_write_pkm4crpa
  integer,public:: nwf
  logical,private:: Wpkm4crpa=.false., keepeig_mlw = .true.
  logical,private:: debug=.false.
  private
  complex(8),allocatable:: geigW(:,:,:,:),cphiW(:,:,:,:)
  real(8),allocatable,private:: evud(:,:,:)
contains
function get_geig_wan(q, isp, mpi_mode, comm) result(geig_wan)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    complex(8) :: geig_wan(ngpmx*nspc,nwf)
    logical :: time_reversal_search
    integer :: iq, dummy
    real(8) :: quu(3)
    if(present(mpi_mode)) then
      if(mpi_mode) call rx('get_geig_wan: mpi_mode is not implemented')
    endif
    time_reversal_search = .false.
    call iqindx2_(q, iq)
    if(iq == 0) time_reversal_search = .true.
    if(time_reversal_search) then
      call readgeigW(-q, dummy, isp, quu, geig_wan)
      geig_wan = conjg(geig_wan)
    else
      call readgeigW( q, dummy, isp, quu, geig_wan)
    endif
  end function get_geig_wan

  function get_cphi_wan(q, isp, mpi_mode, comm) result(cphi_wan)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    complex(8) :: cphi_wan(ndima*nspc,nwf)
    logical :: time_reversal_search
    integer :: iq, dummy
    real(8) :: quu(3)
    if(present(mpi_mode)) then
      if(mpi_mode) call rx('get_cphi_wan: mpi_mode is not implemented')
    endif
    time_reversal_search = .false.
    call iqindx2_(q, iq)
    if(iq == 0) time_reversal_search = .true.
    if(time_reversal_search) then
      call readcphiW(-q, dummy, isp, quu, cphi_wan)
      cphi_wan = conjg(cphi_wan)
    else
      call readcphiW( q, dummy, isp, quu, cphi_wan)
    endif
  end function get_cphi_wan

  subroutine onoff_write_pkm4crpa(lll)
    logical:: lll
    Wpkm4crpa=lll
  end subroutine onoff_write_pkm4crpa
  subroutine init_readeigen_mlw_noeval() ! replace cphi and geig for hwmat ! this should be called after init_readgeigen2
  !xxxxxxxxxxxxxx only for nspc=1. Need fixing for nspc=2  
    use m_cmdopt_registry, only: c0_mlo
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
    logical :: keepeigen,mlocase
    integer:: ikpx,ifi,ndimMTO
    character*(8):: fname
    keepeig_mlw = .True. !keepeigen()
    if(ipr) write(6,*)' init_readeigen_mlw_noeval'
    ! --- Readin MLWU/D, MLWEU/D, and UUq0U/D
    mlocase=c0_mlo
    
    ! if(mlocase) then
    !   iko_ix=1
    !   iko_fx=nband
    !   open(newunit=ifi,file='__cmlo.info',form='unformatted') 
    !   read(ifi) ndimMTO,nqbz !,nqirr,nMTO,mrecb
    !   read(ifi) !ix(1:ndimMTO),qplistgw(1:3,nqirr)
    !   close(ifi)
    !   nwf=ndimMTO
    ! else
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
     ! deallocate(evud)
     ! allocate(evud(nwf,nqi,nsp)) !nqtt
     ! evud = 0d0
    ! endif  
    allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
    cbwf = 0d0
    ! allocate(ovlmW_inv(nwf,nwf,nqtt,nsp), source = (0d0,0d0))
    ! forall(iwf=1:nwf) ovlmW_inv(iwf,iwf,1:nqtt,1:nsp) = 1d0
    
    if(Wpkm4crpa) then
      fname='pkm4crpa'
      open(newunit=ifi,file=fname,form='formatted',status='unknown')
      write(ifi,"('== p_km^alpha in PRB83,121101 ! weight in l-subspace ==')")
      write(ifi,"('( = c^sigma_km in book of 45th IFFK by Ersoy)')")
      write(ifi,"(8i8)") nqtt,nwf,nsp,iko_ix,iko_fx
      write(ifi,"('       |pkm|**2          ib      iq     is       q(1:3)')")
    endif

    write(stdo,ftox) 'nqbz nqtt=',nqbz,nqtt
    do ikp = 1,nqtt
       iqbz = mod(ikp,nqbz)
       if (iqbz == 0) iqbz = nqbz
       iq0i = (ikp - iqbz)/nqbz
       do is= 1,nsp
         ! if(mlocase) then
         !   call readcmlo( qtt(:,ikp), is, nspx,nband, ndimMTO, cbwf(:,:,ikp,is), ovlmW_inv(:,:,ikp,is))
         !   !           write(*,*) 'readcmlo check ---', ikp,is,sum(abs(cbwf(1:nband,1:ndimMTO,ikp,is)))
         !   goto 889
         ! endif
         if (iq0i == 0) then !first block without adding Q0P
           do ib = iko_ix,iko_fx
             do iwf= 1,nwf
               cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
             enddo
           enddo
         else  !
           !2025-11-17. TK think this may be ok?
           !   <psi(k+q0,n) | psi(k+q0,m)^B>
           ! = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
           ! ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> * <psi(k,l) |psi(k,m)^B>
           ! psi^B : bloch fn. corresponding to maxloc Wannier fn.
           do ib = iko_ix,iko_fx
             do iwf= 1,nwf
               cbwf(ib,iwf,ikp,is) = sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is)) *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
             enddo
           enddo
         endif
889      continue
         !! --- write pkm4crpa
         if(Wpkm4crpa) then
           do ib = iko_ix,iko_fx
             write(ifi,"(f19.15, 3i8, 3f13.6 )") sum(abs(cbwf(ib,1:nwf,ikp,is))**2),ib, ikp, is, qtt(1:3,ikp)
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
    ! if(.not.mlocase) deallocate(dnk,uum)
    deallocate(dnk,uum)
    mrecb_o = mrecb * nwf / nband
    mrecg_o = mrecg * nwf / nband
    if(keepeig_mlw) then
       if(ipr) write(6,*)' xxx nband=',nband
       allocate(geig2(ngpmx*nspc,nband))  !nqtt -->nqi
       allocate(cphi2(ndima*nspc,nband))
       allocate(geigW(ngpmx*nspc,nwf,nqtt,nsp))
       allocate(cphiW(ndima*nspc,nwf,nqtt,nsp))
       geigW = 0d0
       cphiW = 0d0
       do ikp= 1,nqtt ! nqi
          do is= 1,nsp
             if(debug) write(6,"(' ikp=',i5,3f10.5)") ikp,qtt(:,ikp)
             ! call readgeig(qtt(:,ikp),is, qu,geig2)
             geig2 = readgeigf(qtt(:,ikp), is)
             ! if(debug)print *,'qqqqqq1',qu,'qqqqqq2',qtt(:,ikp)
             ! if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 1111')
             ! call readcphi(qtt(:,ikp),is, qu,cphi2)
             cphi2 = readcphif(qtt(:,ikp), is)
             ! if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 2222')
             do iwf= 1,nwf
               do ib= iko_ix,iko_fx !band index
                 geigW(:,iwf,ikp,is) = geigW(:,iwf,ikp,is) + geig2(:,ib)*cbwf(ib,iwf,ikp,is)
                 cphiW(:,iwf,ikp,is) = cphiW(:,iwf,ikp,is) + cphi2(:,ib)*cbwf(ib,iwf,ikp,is)
               enddo
             enddo
          enddo
       enddo
       deallocate(geig2,cphi2) !,geig,cphi)
       ! if(allocated(cphi)) deallocate(cphi)
       ! if(allocated(geig)) deallocate(geig)
    else
       call rx('KeepEigen_MLW=F not implemented')
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
    ! logical, intent(in), optional :: dual
    ! if(init2) call rx( 'readgeig_mlw: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
    if(ngp_in < ngp(iq)) then
      if(ipr) write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q,nspc
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
    ! if(init2) call rx( 'readcphi_mlw: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
    !   if(keepeig) then
    cphif(1:ndima*nspc,1:nwf) = cphiW(1:ndima*nspc,1:nwf,iq,isp)
    !   else
    !      ikpisp= isp + nsp*(iq-1)
    !      read(ifcphi_mlw) cphif(1:ndimanspc,1:nwf)
    !   endif
  end subroutine readcphiW
  end module m_wan_wfs
