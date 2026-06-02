module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
  use m_ftox
  use m_lgunit,only: stdo
  use m_struc_def,only: s_rv1
  use m_cmdopt_registry, only: c0_band, c0_cls, c0_fermisurface, c0_mkprocar, c0_pdos, c0_tdos, c0_zmel0
  type(s_rv1),allocatable,protected,public :: t_wtkb(:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
  integer,protected,public:: nevmx
  public :: m_subzi_init, m_subzi_bzintegration, m_subzi_bcast_wtkb, m_subzi_copy_wtkb
  private
  real(8),allocatable :: wtkb(:,:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
contains
  subroutine m_subzi_init() ! Set nevmx and allocate wtkb.
    use m_ext,only: sname
    use m_lmfinit, only: nspc,lmet=>bz_lmet, qbg=>zbak,nspx,lso
    use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
    !    use m_suham,only: ndhamx=>ham_ndhamx
    use m_igv2x,only: ndhamx=>nbandmx
    use m_mkpot,only:  qval
    use m_lmfinit,only: lso
    !   lmet/=0 : allocate tetrahedron weight wtkb
    !   ndhamx : leading dimension of wtkb
    !   nsp    : 2 for spin-polarized case, otherwise 1
    !   nkp    : number of irreducible k-points (bzmesh.f)
    !   nevmx  : maximum number of eigenvectors to find 
    implicit none
    integer :: nkp
    real(8) :: zval
    call tcn('m_subzi_init')
    if(lmet>0) then
      nkp  = bz_nkp
      if(allocated(t_wtkb)) deallocate(t_wtkb)       
      allocate(t_wtkb(nspx,nkp))
    endif
!    if(c0_pdos.or.c0_mkprocar.or.c0_zmel0.or.c0_cls) then
    if(c0_pdos.or.c0_mkprocar.or.c0_cls) then
      nevmx= ndhamx  !all bands
    elseif(c0_tdos.or. c0_band.or.c0_fermisurface) then !nevmx=0 implies eigenvalue-only mode
      nevmx = merge(ndhamx, 0, lso==1)
!    if(c0_tdos.or. c0_band.or.c0_fermisurface) then !nevmx=0 implies eigenvalue-only mode
!      nevmx = merge(ndhamx, 0, lso==1)
!    elseif(c0_pdos.or.c0_mkprocar.or.c0_zmel0.or.c0_cls) then
!      nevmx= ndhamx  !all bands
    else  !just above occipied bands. (tetrahedron method may require a little more than zval/2)
      zval = qval-qbg
      nevmx = ceiling(zval)/2
      if(lmet /= 0) nevmx = max(nevmx+nevmx/2,9) !probably safer setting nevmx for metal. Rough estimation.
      nevmx = min(nevmx*nspc, ndhamx) !nspc=2 for lso=1
      nevmx=nevmx+5 !+5 is for safer setting. At least +1 is required...
    endif
    call tcx('m_subzi_init')
  end subroutine m_subzi_init
  subroutine m_subzi_bzintegration(evlall,efermi,sev,qvalm,vmag)
    use m_bzintegration2,only: bzintegration2
    use m_igv2x,only: ndhamx=>nbandmx
    use m_lmfinit, only: nspx, lmet=>bz_lmet
    use m_mkqp,only: nkp=>bz_nkp
    implicit none
    intent(in)::                   evlall
    intent(out)::                         efermi,sev,qvalm,vmag
    real(8):: evlall(:,:,:),sev, qvalm(2),efermi,vmag
    integer::nx(3),nbmx
    nx=shape(evlall)
    nbmx=nx(1)
    if(lmet>0 .and. (.not.allocated(wtkb))) allocate(wtkb(ndhamx,nspx,nkp))
    call bzintegration2(nbmx,evlall, efermi,sev,wtkb,qvalm,vmag)
  end subroutine m_subzi_bzintegration
  subroutine m_subzi_bcast_wtkb()
    use m_lmfinit, only: nspx, lmet=>bz_lmet
    use m_mkqp,only: nkp=>bz_nkp
    use m_qplist, only: owner
    use m_igv2x,only: nbandmx
    use m_MPItk, only: master, master_mpi, procid, comm
    use mpi, only: mpi_double_precision, mpi_status_size
    implicit none
    integer:: status(MPI_Status_size), iq, isp, itag, ierr
    if(lmet==0) return
    do iq = 1, nkp
      do isp = 1, nspx
        itag = (iq-1)*nspx + isp
        if(owner(isp,iq) == procid .and. .not.allocated(t_wtkb(isp,iq)%v)) allocate(t_wtkb(isp,iq)%v(nbandmx))
        if(master_mpi) then
          if(owner(isp,iq) == procid) t_wtkb(isp,iq)%v = wtkb(:,isp,iq)
          if(owner(isp,iq) /= procid) call mpi_send(wtkb(:,isp,iq), nbandmx, mpi_double_precision, &
                                                  & owner(isp,iq), itag, comm, ierr)
        else
          if(owner(isp,iq) == procid) call mpi_recv(t_wtkb(isp,iq)%v, nbandmx, mpi_double_precision, &
                                                  & master, itag, comm, status, ierr)
        endif
      enddo
    enddo
  end subroutine m_subzi_bcast_wtkb
  subroutine m_subzi_copy_wtkb(isp1, ik1, isp2, ik2) !make copy of t_wtkb
    use m_igv2x,only: nbandmx
    integer, intent(in):: isp1, ik1, isp2, ik2
    if(.not.allocated(t_wtkb(isp1,ik1)%v)) call rx('m_subzi_copy_wtkb: t_wtkb not allocated')
    if(allocated(t_wtkb(isp2,ik2)%v)) deallocate(t_wtkb(isp2,ik2)%v)
    allocate(t_wtkb(isp2,ik2)%v(nbandmx))
    t_wtkb(isp2,ik2)%v = t_wtkb(isp1,ik1)%v
  end subroutine m_subzi_copy_wtkb
end module m_subzi
