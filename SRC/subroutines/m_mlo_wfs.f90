! MLO version of readeigenW and readcphiW in readeigen.f90
module m_mlo_wfs
  use m_lgunit,      only: stdo
  use m_iqindx_qtt,  only: iqindx2_
  use m_hamindex,    only: ngpmx, nqtt, qtt, symops, ngrp, plat
  use m_genallcf_v3, only: nsp => nspin, ndima, nband, nspc, nspx
  use m_readeigen,   only: readgeigf => readgeigf_mpi, readcphif => readcphif_mpi
  use m_keyvalue,    only: getkeyvalue
  use,intrinsic :: ieee_arithmetic
  use m_ftox
  implicit none
  public :: cmlo_init, get_geig_cmlo, get_cphi_cmlo
  integer, public, protected :: nmlo, nMTO
  private
  complex(8), allocatable :: cmlo(:,:,:,:)
  integer :: nqirr, ifile_cmlo
  logical :: keep_mlo, init = .true., debug = .false.
  integer, allocatable :: ix(:)
  real(8), allocatable :: qplistgw(:,:)
contains
  subroutine cmlo_init()
    integer :: ifihh, nqbz, mrecbb, istat
    if(.not.init) return
    call getkeyvalue("GWinput","KeepCMLO",keep_mlo,default=.true.)
    open(newunit=ifihh, file='__cmlo.info', form='unformatted')
    read(ifihh) nmlo, nqbz, nqirr, nMTO, mrecbb
    allocate(ix(nmlo), qplistgw(3,nqirr))
    read(ifihh) ix, qplistgw
    close(ifihh)
    open(newunit=ifile_cmlo,file='__cmlo.data',action='read',form='unformatted',access='direct',recl=mrecbb)
    init = .false.
  end subroutine cmlo_init

  function get_cmlo(q, isp) result(cmlo_q_isp)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    complex(8) :: cmlo_q_isp(nband,nmlo)
    integer :: iq, ikp, is
    if(keep_mlo .and. .not.allocated(cmlo)) then
      allocate(cmlo(nband,nmlo,nqtt,nsp))
      if(debug) write(stdo,ftox) 'xxx nqtt:', nqtt, nband
      do ikp = 1, nqtt
        do is = 1, nsp
          call read_cmlo(qtt(:,ikp), is, cmlo(:,:,ikp,is))
        enddo
      enddo
    endif
    if(keep_mlo) then
      call iqindx2_(q, iq)
      if(debug) write(stdo,ftox) 'xxx set cmlo', q, iq
      cmlo_q_isp = cmlo(:,:,iq,isp)
    else
      call read_cmlo(q, isp, cmlo_q_isp)
    endif
  end function get_cmlo

  function get_geig_cmlo(q, isp, mpi_mode, comm) result(geig_cmlo)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    logical :: time_reversal_search
    integer :: iq
    complex(8) :: geig_cmlo(ngpmx*nspc,nmlo), geig(ngpmx*nspc,nband), cmlo_ik_is(nband,nmlo)
    time_reversal_search = .false.
    call iqindx2_(q, iq)
    if(iq < 1) time_reversal_search = .true.
    if(time_reversal_search) then
      geig = conjg(readgeigf(-q, isp, mpi_mode, comm))
      cmlo_ik_is = conjg(get_cmlo(-q, isp))
    else
      geig = readgeigf(q, isp, mpi_mode, comm)
      cmlo_ik_is = get_cmlo(q, isp)
    endif
    geig_cmlo = matmul(geig, cmlo_ik_is)
  end function get_geig_cmlo

  function get_cphi_cmlo(q, isp, mpi_mode, comm) result(cphi_cmlo)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    logical :: time_reversal_search
    integer :: iq
    complex(8) :: cphi_cmlo(ndima*nspc,nmlo), cphi(ndima*nspc,nband), cmlo_ik_is(nband,nmlo)
    time_reversal_search = .false.
    call iqindx2_(q, iq)
    if(iq < 1) time_reversal_search = .true.
    if(time_reversal_search) then
      cphi = conjg(readcphif(-q, isp, mpi_mode, comm))
      cmlo_ik_is = conjg(get_cmlo(-q, isp))
    else
      cphi = readcphif(q, isp, mpi_mode, comm)
      cmlo_ik_is = get_cmlo(q, isp)
    endif
    cphi_cmlo = matmul(cphi, cmlo_ik_is)
  end function get_cphi_cmlo

  subroutine read_cmlo(qtarget, isp, cmlo_out, ovlm_inv)
    real(8), intent(in) :: qtarget(3)
    integer, intent(in) :: isp
    complex(8), intent(out) :: cmlo_out(nband,nmlo)
    complex(8), intent(out), optional :: ovlm_inv(nmlo,nmlo)
    integer :: i, igg, iqqisp, j, ig, iq, iqq, istat
    real(8) :: qp(3), qx(3), qxx(3)
    logical :: found
    real(8), external :: tolq !eps=1d-8
    ! find iq for given qtarget
    found = .false.
    FindIqIgg: do iq=1,nqirr
      qp = qplistgw(:,iq)
      do ig=1,ngrp
        qx = matmul(transpose(plat), qtarget-matmul(symops(:,:,ig),qp))
        qxx = qx-nint(qx) !qx-ndiff !translation of qx
        if(sum(abs(qxx))<tolq()) then
          iqq = iq
          igg = ig
          found = .true.
          exit FindIqIgg
        endif
      enddo
    enddo FindIqIgg
    if(.not.found) call rx('read_cmlo: can not find ig and iq')
    iqqisp = isp + nspx*(iqq-1)
    read(ifile_cmlo, rec=iqqisp) cmlo_out
    RotCMLO: block
      use m_rotwave, only: rotmatMTO
      use m_lapack, only: zminv => zminv_h
      complex(8) :: rotmatt(nmlo,nmlo), rotmat(nMTO,nMTO)
      complex(8) :: ovlm(nmlo,nmlo)
      call rotmatMTO(igg, qp,qtarget,nMTO, rotmat)
      forall(i=1:nmlo,j=1:nmlo) rotmatt(i,j)=rotmat(ix(i),ix(j))
      cmlo_out = matmul(cmlo_out,dconjg(transpose(rotmatt)))
      if(present(ovlm_inv)) then
        ovlm = matmul(dconjg(transpose(cmlo_out)), cmlo_out)
        istat = zminv(ovlm, n=nmlo)
        ovlm_inv = ovlm
      endif
    endblock RotCMLO
  end subroutine read_cmlo
end module m_mlo_wfs
