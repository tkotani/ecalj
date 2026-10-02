! MLO version of readeigenW and readcphiW in readeigen.f90
module m_mlo_wfs
  use m_lgunit,      only: stdo
  use m_iqindx_qtt,  only: iqindx2_
  use m_hamindex,    only: ngpmx, nqtt, qtt, symops, ngrp, plat
  use m_genallcf_v3, only: nsp => nspin, ndima, nband, nspc, nspx
  use m_readeigen,   only: readgeigf => readgeigf_mpi, readcphif => readcphif_mpi
  use m_GWinput,     only: gwinput_init, gwinput_loaded,  tg_KeepCMLO => KeepCMLO
  use m_ftox
  implicit none
  public :: cmlo_init, get_geig_cmlo, get_cphi_cmlo, get_cmlo, get_cmlo_qirr, write_pkm4crpa_mlo
  integer, public, protected :: nmlo, nMTO
  !> rnorm(i,isp): square integral of the real-space MLO i over the GW k mesh, eq. (1) below. A check only: the MLOs of
  !> __cmlo are normalized in real space already (Hreduction, with the norms of HamRsMLO), so it is close to 1. (2026-10-02)
  real(8), allocatable, public, protected :: rnorm(:,:)
  private
  complex(8), allocatable :: cmlo(:,:,:,:)
  integer :: nqirr, ifile_cmlo, nqbz_cmlo
  integer :: nbandmx_cmlo !rows per record in __cmlo.data; NOT necessarily nband
  logical :: keep_mlo, init = .true., debug = .false.
  integer, allocatable :: ix(:)
  real(8), allocatable :: qplistgw(:,:)
contains
  subroutine cmlo_init()
    integer :: ifihh, nqbz, mrecbb, istat
    if(.not.init) return
    call gwinput_init()
    if (gwinput_loaded) then
      keep_mlo = tg_KeepCMLO
    else
      call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
    endif
    open(newunit=ifihh, file='__cmlo.info', form='unformatted')
    read(ifihh) nmlo, nqbz, nqirr, nMTO, mrecbb
    nqbz_cmlo = nqbz
    allocate(ix(nmlo), qplistgw(3,nqirr))
    read(ifihh) ix, qplistgw
    close(ifihh)
    open(newunit=ifile_cmlo,file='__cmlo.data',action='read',form='unformatted',access='direct',recl=mrecbb)
    nbandmx_cmlo = mrecbb/(16*nmlo) !mrecbb = 2*nbandmx*nmlo*8
    write(stdo,ftox)'cmlo_init: nmlo nbandmx(record) nband(here) =',nmlo,nbandmx_cmlo,nband
    if(nbandmx_cmlo < nband) call rx('cmlo_init: __cmlo.data has fewer bands than nband')
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
    call set_rnorm()
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
    call set_rnorm()
    cphi_cmlo = matmul(cphi, cmlo_ik_is)
  end function get_cphi_cmlo

  !> The MLO at k is |F_i(k)> = sum_n |psi_kn> C_ni (C = cmlo); the MLO in real space is
  !> F_i0(r) = (1/N_k) sum_k F_i(k)(r) on the Born-von Karman supercell of the GW k mesh, whose square integral is
  !>    rnorm_i = int |F_i0(r)|^2 dr = (1/N_k) sum_k <F_i(k)|F_i(k)>_cell = (1/N_k) sum_k O_ii(k)          (1)
  !> The MLOs are normalized in real space where they are made (Hreduction with the norms of HamRsMLO, 2026-10-02), so this
  !> is a check: close to 1, apart from the difference of the k meshes (mlo_nkabc for HamRsMLO, the GW mesh here) and the
  !> bands above nband. (A normalization at each k, 1/sqrt(O_ii(k)), is not a constant: it changes the real-space orbital and
  !> the interpolated bands, C 0.74 eV off the mesh. It was removed 2026-10-01, user.)
  subroutine set_rnorm()
    integer :: ikp, is, j
    complex(8) :: c(nband,nmlo)
    if(allocated(rnorm)) return
    if(nqbz_cmlo < 1 .or. nqbz_cmlo > nqtt) call rx('m_mlo_wfs set_rnorm: nqbz of __cmlo.info is not in 1..nqtt (rerun lmf --jobgw=1 --mlo and mlo --mlo)')
    allocate(rnorm(nmlo,nspx), source=0d0)
    do ikp = 1, nqbz_cmlo        !qtt(:,1:nqbz) is the regular BZ mesh (the q0-shifted copies follow)
      do is = 1, nspx
        c = get_cmlo(qtt(:,ikp), is)
        forall(j=1:nmlo) rnorm(j,is) = rnorm(j,is) + sum(abs(c(:,j))**2)
      enddo
    enddo
    rnorm = rnorm/nqbz_cmlo
    do is = 1, nspx
      write(stdo,"(' m_mlo_wfs: square integral of the real-space MLOs on the GW mesh (check, ~1), isp=',i2,' nqbz=',i6,' :',100f8.4)") &
           is, nqbz_cmlo, rnorm(:,is)
    enddo
  end subroutine set_rnorm

  !> pkm4crpa for the cRPA of the MLO model (hwmatK_MPI --mlo, mode 10011; read by m_pkm4crpa in hx0fp0 10011).
  !> The MLOs at k are |phi_j> = sum_n |psi_kn> C_nj (C = cmlo, nband rows), not orthonormal: O = C^+ C.
  !> The weight of band n in the MLO subspace is the diagonal of the projector onto it,
  !>    p_kn = [ C O^-1 C^+ ]_nn ,  0 <= p_kn <= 1,  sum_n p_kn = nmlo when the MLOs lie within the nband bands,
  !> the MLO counterpart of p_kn = sum_m |<psi_kn|w_m>|^2 of the Wannier functions (PRB 83, 121101).
  !> hx0fp0 weights each transition kn -> k+q n' of chi0 by 1 - p_kn p_k+q,n'. (2026-10-01 23:20)
  !> Same file layout as m_wan_wfs wrote it; every q of qtt (also the q0-shifted ones hx0fp0 asks for).
  subroutine write_pkm4crpa_mlo()
    use m_lapack, only: zminv => zminv_h
    integer :: ifi, ikp, is, ib, istat
    complex(8) :: c(nband,nmlo), oinv(nmlo,nmlo), x(nband,nmlo)
    real(8) :: pkn(nband), psum, psummin, psummax
    if(nspc==2) call rx('write_pkm4crpa_mlo: cRPA of the MLO model with spin-orbit coupling is not implemented')
    open(newunit=ifi, file='pkm4crpa', form='formatted', status='replace')
    write(ifi,"('== p_kn: weight of band n in the MLO subspace, [C (C^+C)^-1 C^+]_nn (m_mlo_wfs) ==')")
    write(ifi,"('( PRB 83, 121101 with MLOs in place of Wannier functions )')")
    write(ifi,"(8i8)") nqtt, nmlo, nspx, 1, nband
    write(ifi,"('       |pkm|**2          ib      iq     is       q(1:3)')")
    psummin = 1d99; psummax = -1d99
    do ikp = 1, nqtt
      do is = 1, nspx
        c = get_cmlo(qtt(:,ikp), is)
        oinv = matmul(dconjg(transpose(c)), c)
        istat = zminv(oinv, n=nmlo)
        x = matmul(c, oinv)
        do ib = 1, nband
          pkn(ib) = dreal(sum(x(ib,:)*dconjg(c(ib,:))))
          write(ifi,"(f19.15, 3i8, 3f13.6 )") pkn(ib), ib, ikp, is, qtt(1:3,ikp)
        enddo
        psum = sum(pkn); psummin = min(psummin, psum); psummax = max(psummax, psum)
      enddo
    enddo
    close(ifi)
    write(stdo,ftox)'write_pkm4crpa_mlo: nqtt nmlo nband=', nqtt, nmlo, nband, &
         ' sum_n p_kn min max (= nmlo when the MLOs lie within nband)=', ftof(psummin,6), ftof(psummax,6)
  end subroutine write_pkm4crpa_mlo

  !> cmlo at a q that is literally in qplistgw (no symmetry rotation, no m_hamindex needed).
  !> ok=.false. if q is not on that list. Used by hqpe_sc, which only ever asks for those q.
  subroutine get_cmlo_qirr(qtarget, isp, cmlo_out, ok)
    real(8), intent(in) :: qtarget(3)
    integer, intent(in) :: isp
    complex(8), intent(out) :: cmlo_out(nband,nmlo)
    logical, intent(out) :: ok
    integer :: iq, iqqisp
    real(8), external :: tolq
    ok = .false.
    do iq = 1, nqirr
      if(sum(abs(qplistgw(:,iq)-qtarget)) < tolq()) then
        iqqisp = isp + nspx*(iq-1)
        call read_cmlo_rec(iqqisp, cmlo_out)
        ok = .true.
        return
      endif
    enddo
    cmlo_out = (0d0,0d0)
  end subroutine get_cmlo_qirr

  !> One record of __cmlo.data is (nbandmx_cmlo, nmlo); read it whole, then keep the
  !> first nband rows. Reading straight into an (nband,nmlo) array mis-aligns the
  !> columns whenever nband < nbandmx_cmlo.
  subroutine read_cmlo_rec(irec, cmlo_out)
    integer, intent(in) :: irec
    complex(8), intent(out) :: cmlo_out(nband,nmlo)
    complex(8), allocatable :: buf(:,:)
    allocate(buf(nbandmx_cmlo,nmlo))
    read(ifile_cmlo, rec=irec) buf
    cmlo_out(1:nband,1:nmlo) = buf(1:nband,1:nmlo)
    deallocate(buf)
  end subroutine read_cmlo_rec

  subroutine read_cmlo(qtarget, isp, cmlo_out, ovlm_inv)
    use m_rotwave, only: rotmatMTO
    use m_lapack,  only: zminv => zminv_h
    real(8), intent(in) :: qtarget(3)
    integer, intent(in) :: isp
    complex(8), intent(out) :: cmlo_out(nband,nmlo)
    complex(8), intent(out), optional :: ovlm_inv(nmlo,nmlo)
    integer :: i, igg, iqqisp, j, ig, iq, iqq, istat
    real(8) :: qp(3), qx(3), qxx(3)
    complex(8) :: rotmatt(nmlo,nmlo), rotmat(nMTO,nMTO), ovlm(nmlo,nmlo)
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
    call read_cmlo_rec(iqqisp, cmlo_out)
    call rotmatMTO(igg, qp, qtarget, nMTO, rotmat)
    forall(i=1:nmlo, j=1:nmlo) rotmatt(i,j) = rotmat(ix(i),ix(j))
    cmlo_out = matmul(cmlo_out, dconjg(transpose(rotmatt)))
    if(present(ovlm_inv)) then
      ovlm = matmul(dconjg(transpose(cmlo_out)), cmlo_out)
      istat = zminv(ovlm, n=nmlo)
      ovlm_inv = ovlm
    endif
  end subroutine read_cmlo
end module m_mlo_wfs
