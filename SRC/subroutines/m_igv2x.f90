!>Return G vectors (integer sets) for given q points specifiec by qplist(:,iq)
module m_igv2x 
  use m_struc_def,only: s_nv2
  public:: m_igv2xall_init, m_igv2x_setiq
  integer,protected,pointer,public :: igv2x(:,:)
  integer,protected,pointer,public :: napw,ndimh,ndimhx
  integer,allocatable,target,protected,public:: ndimhall(:)
  private
  integer,protected,private,target ::napw_z,ndimh_z,ndimhx_z
  type(s_nv2),allocatable,target,protected,private:: igv2xall(:)
  integer,allocatable,target,protected,private:: napwall(:),ndimhxall(:),igv2x_z(:,:)
  logical,private:: init=.true.
contains
  subroutine m_Igv2x_setiq(iq) ! Return G vectors for given qplist(:,iq)
    integer,intent(in):: iq
    napw  => napwall(iq)
    ndimh => ndimhall(iq)  !nlmto+napw
    ndimhx=> ndimhxall(iq) !nlmto+napw (but x2 when SO=1)
    igv2x => igv2xall(iq)%nv2
  end subroutine m_Igv2x_setiq
  subroutine m_igv2xall_init(iqini,iqend) !initialization for qplist(iqini:iqend)  
    use m_qplist,only: qplist 
    use m_MPItk,only: master_mpi,procid,master
    integer:: iqini,iqend,iq
    real(8):: qp(3)
    call tcn('m_igv2xall_init')
    allocate(igv2xall(iqini:iqend),napwall(iqini:iqend),ndimhall(iqini:iqend),ndimhxall(iqini:iqend))
    do iq = iqini, iqend
       qp = qplist(:,iq)
       call m_Igv2x_init(qp)   ! Get napw and so on for given qp
       allocate( igv2xall(iq)%nv2(3,napw_z) )
       igv2xall(iq)%nv2 = igv2x_z
       napwall(iq) = napw_z
       ndimhall(iq)= ndimh_z
       ndimhxall(iq)=ndimhx_z
    enddo
    call tcx('m_igv2xall_init')
  end subroutine m_igv2xall_init
  subroutine m_igv2x_init(qp) !Set napw and igv2x_z for given qp
    use m_lmfinit,only: pwmode=>ham_pwmode,pwemax,alat=>lat_alat,stdo,nspc
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
    use m_MPItk,only: master_mpi,procid,master
    use m_lmfinit,only: nlmto !    use m_shortn3_qlat,only: shortn3_qlat,nout,nlatout
    use m_ftox
    integer:: ifiese,imx11(1,1)
    integer,allocatable ::  kv_iv(:,:)
    real(8):: ppin(3),qp(3),qqq(3),pwgmax,dum,platt(3,3)
    logical:: debug,cmdopt0
    logical,save:: init=.true.
    integer:: iout,iapw,napwx,i !,nout,nlatout(3,noutmx)
    call tcn('m_igv2x_init')
    debug = cmdopt0('--debugbndfp')
    platt=transpose(plat)
    if(allocated(igv2x_z)) deallocate(igv2x_z)
    if (0<pwemax .and. mod(pwmode,10)>0) then !with APWs
       pwgmax = pwemax**.5d0
       qqq = 0d0
       if (mod(pwmode/10,10) == 1) qqq = qp !pwmode 1 in 10th digit means q-dependent nw
       call getgv2(alat,plat,qlat, qqq, pwgmax,1, napw_z, imx11)   ! get nqpn. # of G vector for |q+G| < pwgmax
       allocate(igv2x_z(3,napw_z))
       call getgv2(alat,plat,qlat, qqq, pwgmax,2, napw_z, igv2x_z) ! for eigenfunctions (psi)
!       write(stdo,ftox) ftof(qqq),' ',ftof(qp)
!       do i=1,napw_z
!          write(stdo,ftox) 'iiiiiiiggggggggg',i,'  ',igv2x_z(:,i)
!       enddo
    else !No APWs
       napw_z=0
       allocate(igv2x_z(1,1))  
    endif
    ndimh_z = nlmto+napw_z     ! ldim+napw
    if (mod(pwmode,10)==2) ndimh_z = napw_z !APW-only mode
    ndimhx_z = ndimh_z*nspc     !this is iq-dependent.
    call tcx('m_igv2x_init')
  end subroutine m_igv2x_init
end module m_igv2x
