!> Central registry of cmdline options whose name takes a value
!! (the `cmdopt2('--foo=', outs)` family). Two roles in one module:
!!
!!   1. Typo detection. `is_known_cmdopt2(flag)` checks a `--word`
!!      against CMDOPT2_REGISTRY before m_toml_override decides between
!!      "TOML override" and "user typo, abort".
!!
!!   2. Single-shot parsing into typed module variables. Callsites do
!!      `use m_cmdopt_registry, only: c2_jobgw` and read the cached
!!      value directly instead of the old `cmdopt2('--jobgw=', outs);
!!      read(outs,*) jobgw` pattern. load_cmdopt2_registry() is called
!!      once at startup (from m_args::m_setargs); subsequent calls are
!!      no-ops.
!!
!! Maintenance:
!!   - One entry per cmdopt2 callsite in SRC/subroutines/ and SRC/main/.
!!   - Drop the trailing `=` from CMDOPT2_REGISTRY; the matcher in
!!     m_toml_override compares `--<word>` after splitting on `=`.
!!   - When adding a new cmdopt2, add to CMDOPT2_REGISTRY AND declare
!!     a `c2_<name>` module variable AND parse it in
!!     load_cmdopt2_registry().
!!   - cmdopt0 (flag-only) options are NOT registered: typo detection
!!     applies only to value-bearing forms.
!!
!! Generated baseline (2026-06-02) by:
!!   grep -rEoh "cmdopt2\('-[^']+'" SRC/subroutines/ SRC/main/ | sort -u
module m_cmdopt_registry
  implicit none
  private
  public :: is_known_cmdopt2
  public :: load_cmdopt2_registry

  !========================================================================
  ! Cached cmdopt2 values. Populated once by load_cmdopt2_registry().
  ! `protected` so callers can read but not write.
  !
  ! Sentinel conventions:
  !   integer values use -1 for "not present" (none of the real cmdopt2
  !     ints can legitimately be -1 -- they are stage selectors / MPI
  !     group sizes / spin indices / DOS point counts).
  !   real values default to 0.0 and carry a separate *_set flag because
  !     0.0 is often a legitimate value (Ef shift, energy window edge).
  !   string values default to '' and carry a _set flag.
  !========================================================================

  integer, public, protected, save :: c2_jobgw         = -1
  integer, public, protected, save :: c2_job           = -1
  integer, public, protected, save :: c2_nb            = -1
  integer, public, protected, save :: c2_nk            = -1
  integer, public, protected, save :: c2_nww           = -1
  integer, public, protected, save :: c2_sp1           = -1
  integer, public, protected, save :: c2_sp2           = -1
  integer, public, protected, save :: c2_ndos          = -1

  real(8), public, protected, save :: c2_cutuu         = 0d0
  real(8), public, protected, save :: c2_emin_eV       = 0d0
  real(8), public, protected, save :: c2_emax_eV       = 0d0
  real(8), public, protected, save :: c2_EfermiShifteV = 0d0
  logical, public, protected, save :: c2_cutuu_set         = .false.
  logical, public, protected, save :: c2_emin_set          = .false.
  logical, public, protected, save :: c2_emax_set          = .false.
  logical, public, protected, save :: c2_EfermiShifteV_set = .false.

  character(len=32), public, protected, save :: c2_Wtype = ''
  logical,           public, protected, save :: c2_Wtype_set = .false.

  ! Edit this array when adding/removing cmdopt callsites that take
  ! `=value`. Two flavours live here:
  !
  !   (a) genuine cmdopt2 -- the binary calls cmdopt2('--foo=', outs)
  !       and reads <value> into outs.
  !   (b) cmdopt0 with a literal embedded `=` -- the binary calls
  !       cmdopt0('--foo=show') and the whole token is a flag.
  !       Both forms look like `--foo=...` at the cmdline, so for the
  !       purposes of classify_plain_dashdash they are equivalent.
  !
  ! Single source of truth: keep alphabetised by leading flag.
  character(len=20), parameter :: CMDOPT2_REGISTRY(*) = [character(len=20) :: &
       '--Wtype',           &  ! (a) main_hmagnon.f90: W matrix type selector
       '--cutuu',           &  ! (a) x0kf_ahc.f90: vc matrix cutoff
       '--diag',            &  ! (b) cmdopt0 --diag={default,tridiag,chefsi}
       '--dwnb',            &  ! (b) cmdopt0 --dwnb={mlo,wan}
       '--job',             &  ! (a) hbasfp0/hsfp0/hvccfp0/hahc/hx0fp0/qg4gw/...
       '--jobgw',           &  ! (a) main_lmf.f90: lmf-as-GW-driver stage
       '--nb',              &  ! (a) main_hahc.f90 / main_hrcxq: b-parallel group
       '--nk',              &  ! (a) main_hahc.f90 / main_hmagnon / hrcxq: k-parallel
       '--nww',             &  ! (a) x0kf_ahc.f90: frequency points
       '--quit',            &  ! (b) cmdopt0 --quit={show,ham,mkpot,dmat,band}
       '--sp1',             &  ! (a) hmlo_ovlppair / hmagnon / huumat /...: spin 1
       '--sp2',             &  ! (a) ditto, spin 2
       '-EfermiShifteV',    &  ! (a) m_tetwt / x0kf_ahc: AHC rigid Ef shift (eV)
       '-emax',             &  ! (a) m_writeband: DOS energy max
       '-emin',             &  ! (a) m_writeband: DOS energy min
       '-ndos'              ]  ! (a) m_writeband: DOS points

contains

  !> True if `flag` ("--foo" or "-foo", *without* trailing `=`) is a
  !! registered cmdopt2 name. Case-sensitive (preserves `-EfermiShifteV`).
  function is_known_cmdopt2(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    do i = 1, size(CMDOPT2_REGISTRY)
       if (trim(flag) == trim(CMDOPT2_REGISTRY(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt2

  !> Parse every cmdopt2 entry into the cached module variables above.
  !! Idempotent: subsequent calls are no-ops. Called once at startup
  !! from m_args::m_setargs (after arglist is populated).
  subroutine load_cmdopt2_registry()
    logical, save :: done = .false.
    character(len=120) :: outs
    logical, external :: cmdopt2
    if (done) return
    done = .true.

    ! Integer-valued cmdopt2
    if (cmdopt2('--jobgw=', outs)) then
       read(outs,*) c2_jobgw
       if (c2_jobgw /= 0 .and. c2_jobgw /= 1) &
            call rx('m_cmdopt_registry: --jobgw must be 0 or 1')
    endif
    if (cmdopt2('--job=',   outs)) read(outs,*) c2_job
    if (cmdopt2('--nb=',    outs)) read(outs,*) c2_nb
    if (cmdopt2('--nk=',    outs)) read(outs,*) c2_nk
    if (cmdopt2('--nww=',   outs)) read(outs,*) c2_nww
    if (cmdopt2('--sp1=',   outs)) read(outs,*) c2_sp1
    if (cmdopt2('--sp2=',   outs)) read(outs,*) c2_sp2
    if (cmdopt2('-ndos=',   outs)) read(outs,*) c2_ndos

    ! Real-valued cmdopt2 (raw eV / dimensionless units; conversion is
    ! the caller's responsibility because each callsite supplies its
    ! own default and Ry conversion).
    if (cmdopt2('--cutuu=', outs)) then
       read(outs,*) c2_cutuu;          c2_cutuu_set         = .true.
    endif
    if (cmdopt2('-emin=',   outs)) then
       read(outs,*) c2_emin_eV;        c2_emin_set          = .true.
    endif
    if (cmdopt2('-emax=',   outs)) then
       read(outs,*) c2_emax_eV;        c2_emax_set          = .true.
    endif
    if (cmdopt2('-EfermiShifteV=', outs)) then
       read(outs,*) c2_EfermiShifteV;  c2_EfermiShifteV_set = .true.
    endif

    ! String-valued cmdopt2
    if (cmdopt2('--Wtype=', outs)) then
       c2_Wtype     = trim(outs)
       c2_Wtype_set = .true.
    endif
  end subroutine load_cmdopt2_registry

end module m_cmdopt_registry
