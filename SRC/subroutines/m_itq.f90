!> m_itq — index map and band counts for self-energy / chi0 calculations.
!!
!! This module exposes three pieces of state that are used by build_zmel and the
!! self-energy / chi0 codepaths:
!!
!!   itq(:)      : index array. itq(i) = band index for self-energy output index i.
!!                  - QSGW (sc) and chi0/WV calculators (hrcxq, hx0fp0, hahc):
!!                    identity, size = nband. itq(i) = i for i = 1..nband.
!!                  - Single-pass non-sc hsfp0: size = nbmax-nbmin+1, values =
!!                    [nbmin..nbmax]. Possibly non-trivial mapping.
!!
!!   ntq         : number of self-energy output bands. ntq <= size(itq).
!!                  Determines the second dim of zsecall, omega, eqx, and the row
!!                  count in SEX/SEC output files.
!!                  - QSGW init / chi0 / WV: ntq = nband (full).
!!                  - QSGW correlation: ntq = maxval(nbandmx)+1, set by
!!                    set_nbandmx_for_sigma. itq is left at full nband size so
!!                    that build_zmel from x0kf_zxq can index up to nband; only
!!                    the first ntq entries of itq are read by the self-energy
!!                    output writer.
!!                  - Non-sc hsfp0: ntq = nbmax-nbmin+1.
!!
!!   nbandmx(ip,isp): per-(q-point ip, spin isp) actual band-summation cap.
!!                    Bounded by ebmx_sig (energy threshold) and nbmx_sig (count
!!                    threshold). Self-energy loops in m_sxcf_sc read this and
!!                    set ntqxx = nbandmx(ip,isp) as their effective inner cap
!!                    (ntqxx <= ntq <= size(itq)).
!!
!! Why the split? Historically setitq_hsfp0sc shrank itq to size ntq. This
!! conflicts with chi0/WV calculators that share build_zmel and pass band
!! indices up to nband: they would read itq() out of bounds. The refactor
!! keeps itq sized for nband (a harmless oversize for the QSGW output writer)
!! and only adjusts the output band count ntq.
!!
!! API:
!!   setitq()                      — itq=[1..nband], ntq=nband (full).
!!   set_itq_range(nbmin,nbmax)    — itq=[nbmin..nbmax], ntq=nbmax-nbmin+1.
!!   set_nbandmx_for_sigma(...)    — populate nbandmx; set ntq=maxval(nbandmx)+1;
!!                                   leave itq at its existing (nband-sized) state.
!!   setitq_hsfp0sc(...)           — backward-compat wrapper: setitq() (if needed)
!!                                   + set_nbandmx_for_sigma(...).
!!   setitq_hsfp0(...)             — backward-compat wrapper for set_itq_range.
module m_itq
  use m_struct_from_lmf, only: nband
  use m_lgunit, only: stdo
  use m_cmdopt_registry, only: c0_ntqxx
  implicit none
  integer, allocatable, protected, public :: itq(:)
  integer,             protected, public :: ntq
  integer, allocatable, protected, public :: nbandmx(:,:)
  public :: setitq, set_itq_range, set_nbandmx_for_sigma
  public :: setitq_hsfp0sc, setitq_hsfp0, setitq_mlo
  logical, save, private :: nbandmx_loaded_from_NTQXX = .false.
  private
contains
  subroutine setitq()
    !> Set itq=[1..nband] (identity, full band range), ntq=nband. Re-runnable.
    integer :: i
    if(allocated(itq)) deallocate(itq)
    ntq = nband
    allocate(itq, source=[(i, i=1, ntq)])
  end subroutine setitq
  subroutine setitq_mlo(nmlo)
    integer, intent(in) :: nmlo
    integer :: i
    if(allocated(itq)) deallocate(itq)
    ntq = nmlo
    allocate(itq,source=[(i,i=1,nmlo)])
  end subroutine setitq_mlo
  subroutine set_itq_range(nbmin, nbmax)
    !> Set itq=[nbmin..nbmax], ntq=nbmax-nbmin+1. Used by single-pass non-sc hsfp0.
    integer, intent(in) :: nbmin, nbmax
    integer :: i
    if(allocated(itq)) deallocate(itq)
    ntq = nbmax - nbmin + 1
    allocate(itq, source=[(i, i=nbmin, nbmax)])
  end subroutine set_itq_range
  subroutine set_nbandmx_for_sigma(nbmx_sig, ebmx_sig, eftrue, nspinmx)
    !> For QSGW correlation: compute nbandmx (per-q,per-spin self-energy band cap).
    !> Bounded by energy threshold ebmx_sig and count threshold nbmx_sig.
    !> Reads/writes NTQXX file when --ntqxx flag is active (caches across QSGW iterations).
    !> Reduces ntq to maxval(nbandmx)+1 so that self-energy output (zsecall, eqx,
    !> SEX/SEC files) is correctly sized. itq is NOT shrunk: callers (chi0/WV
    !> calculators) that index itq up to nband must still see full itq.
    use m_read_bzdata, only: qibz, nqibz
    use m_readeigen, only: readeval
    use m_nvfortran, only: findloc
    use m_ftox
    use m_mpi, only: mpi__root
    use m_cmdopt_registry, only: c0_ntqxx
    integer, intent(in) :: nbmx_sig, nspinmx
    real(8), intent(in) :: ebmx_sig, eftrue
    integer :: ifih, nspinmxin, is, ip, nqibzin, iqibz, ierr
    real(8), allocatable :: eqt(:)
    logical :: readntqxx
    if(nbandmx_loaded_from_NTQXX) return  ! NTQXX-mode idempotent: don't recompute
    if(allocated(nbandmx)) deallocate(nbandmx)
    allocate(nbandmx(nqibz, nspinmx), eqt(nband))
    readntqxx = .false.
    if(c0_ntqxx) then
      open(newunit=ifih, file='NTQXX', status='old', iostat=ierr)
      if(ierr == 0) then
         read(ifih,*) nqibzin, nspinmxin
         if(nqibzin == nqibz .or. nspinmxin == nspinmx) readntqxx = .true.
      endif
    endif
    if(readntqxx) then
      do iqibz = 1, nqibz
         read(ifih,*) nbandmx(iqibz,:)
      enddo
      close(ifih)
      nbandmx_loaded_from_NTQXX = .true.
    else
      do is = 1, nspinmx
         do ip = 1, nqibz
            eqt = readeval(qibz(:,ip), is)
            nbandmx(ip,is) = min(findloc(eqt-eftrue > ebmx_sig, value=.true., dim=1)-1, nbmx_sig)
         enddo
      enddo
      if(mpi__root .and. c0_ntqxx) then
         open(newunit=ifih, file='NTQXX')
         write(ifih,ftox) nqibz, nspinmx, ' !nqibz nspinmx. Note NTQXX is used when --ntqxx'
         do iqibz = 1, nqibz
            write(ifih,ftox) nbandmx(iqibz,:), ' !=nbandmx ', iqibz, ' ! =iqibz'
         enddo
         close(ifih)
      endif
    endif
    deallocate(eqt)
    ntq = maxval(nbandmx) + 1   ! reduce output count; itq stays at nband size
  end subroutine set_nbandmx_for_sigma
  subroutine setitq_hsfp0sc(nbmx_sig, ebmx_sig, eftrue, nspinmx)
    !> Backward-compatible wrapper for QSGW correlation phase setup.
    !> If itq is not yet allocated (some callers expect setitq_hsfp0sc to fully
    !> initialize), allocate it as full nband. Then compute nbandmx and reduce ntq.
    integer, intent(in) :: nbmx_sig, nspinmx
    real(8), intent(in) :: ebmx_sig, eftrue
    if(.not. allocated(itq)) call setitq()
    call set_nbandmx_for_sigma(nbmx_sig, ebmx_sig, eftrue, nspinmx)
  end subroutine setitq_hsfp0sc
  subroutine setitq_hsfp0(ngcmx_in, ngpmx_in, tote, nbmin, nbmax, noccxv)
    !> Backward-compatible wrapper for non-sc single-pass hsfp0.
    !> ngcmx_in, ngpmx_in are unused (historical artifact in the API).
    integer, intent(in) :: ngcmx_in, ngpmx_in, nbmin, nbmax, noccxv
    logical, intent(in) :: tote
    if(tote) then
       call set_itq_range(1, noccxv)
    else
       call set_itq_range(nbmin, nbmax)
    endif
  end subroutine setitq_hsfp0
end module m_itq
