!> Central registry of command-line option lookups. Three roles:
!!
!!   1. Single-shot parsing into typed module variables. Callsites do
!!      `use m_cmdopt_registry, only: c0_writeham` (bool) or
!!      `use m_cmdopt_registry, only: c2_jobgw` (int) and read the cached
!!      value directly instead of the old `cmdopt0(...) / cmdopt2(...);
!!      read(outs,*) jobgw` pattern. load_*_registry() are called once at
!!      startup (from m_args::m_setargs); subsequent calls are no-ops.
!!
!!   2. Strict typo detection. validate_arglist() walks every arglist
!!      token after both load_*_registry have populated the known-flag
!!      tables; any `--word` / `-word` not in either table aborts at
!!      startup with a hint. Runs from m_setargs so every program (both
!!      lmf-family and GW utilities) gets typo checking.
!!
!!   3. Retired-syntax diagnostics. validate_arglist also catches the
!!      legacy `-v...=`, `--[...]=`, `--toml.<path>=`, `--pr=`, etc.
!!      forms and aborts with a one-line migration hint pointing at the
!!      canonical `--ctrlg:<path>=<value>` syntax.
!!
!! Maintenance:
!!   - When adding a new cmdopt0 (pure flag), declare a `c0_<name>`
!!     logical above AND add a single `call set0('--<name>', c0_<name>)`
!!     line in load_cmdopt0_registry. The typo-detection table is built
!!     as a side effect of set0() -- no parallel list to keep in sync.
!!   - When adding a new cmdopt2 (=value), declare a `c2_<name>` of
!!     appropriate type AND add a `get2('--<name>', outs)` call in
!!     load_cmdopt2_registry; get2 registers the flag for typo detection.
!!   - `--quit=`, `--diag=`, `--dwnb=` are stored as cmdopt2 strings
!!     (c2_quit, c2_diag, c2_dwnb) even though callers historically used
!!     cmdopt0('--quit=show') etc. The use-site test becomes
!!     `if (c2_quit == 'show')` -- functionally identical, one fewer
!!     string-match per call.
!!
!! Sentinel conventions:
!!   integer:   -1   = "not present" (none of the real cmdopt2 ints can
!!                     legitimately be -1)
!!   real:       0.0 + a separate *_set flag (0.0 may be a real value)
!!   character:  ''  + a separate *_set flag
!!   logical:    .false. (c0_* are pure presence/absence flags)
module m_cmdopt_registry
  implicit none
  private
  public :: load_cmdopt2_registry
  public :: load_cmdopt0_registry
  public :: validate_arglist
  public :: list_cmdopts

  !========================================================================
  ! Cached cmdopt2 values (=value form). Populated by load_cmdopt2_registry.
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

  ! cmdopt0-with-embedded-`=`: --quit={show,ham,mkpot,dmat,band},
  ! --diag={default,tridiag,chefsi}, --dwnb={mlo,wan}. Parsed as cmdopt2.
  character(len=16), public, protected, save :: c2_quit = ''
  character(len=16), public, protected, save :: c2_diag = ''
  character(len=16), public, protected, save :: c2_dwnb = ''

  !========================================================================
  ! Cached cmdopt0 flags (pure presence/absence). Populated by
  ! load_cmdopt0_registry. Variable names mirror the flag name, with
  ! `-` -> `_` and `:` -> `_` for Fortran identifier rules.
  !========================================================================
  logical, public, protected, save :: c0_AHCMAT          = .false.
  logical, public, protected, save :: c0_UUMAT           = .false.
  logical, public, protected, save :: c0_afsym           = .false.
  logical, public, protected, save :: c0_ahc             = .false.
  logical, public, protected, save :: c0_allband         = .false.
  logical, public, protected, save :: c0_avoidgamma      = .false.
  logical, public, protected, save :: c0_band            = .false.
  logical, public, protected, save :: c0_boltztrap       = .false.
  logical, public, protected, save :: c0_cls             = .false.
  logical, public, protected, save :: c0_cmlo            = .false.
  logical, public, protected, save :: c0_corehole        = .false.
  logical, public, protected, save :: c0_cvK             = .false.  ! --cvK:
  logical, public, protected, save :: c0_debug           = .false.
  logical, public, protected, save :: c0_debugbndfp      = .false.
  logical, public, protected, save :: c0_debugpwmat      = .false.
  logical, public, protected, save :: c0_debugsugw       = .false.
  logical, public, protected, save :: c0_debugzmel       = .false.
  logical, public, protected, save :: c0_density         = .false.
  logical, public, protected, save :: c0_dos             = .false.
  logical, public, protected, save :: c0_eigen_at_k      = .false.  ! --eigen-at-k
  logical, public, protected, save :: c0_espot           = .false.
  logical, public, protected, save :: c0_estaticall      = .false.
  logical, public, protected, save :: c0_eszero          = .false.
  logical, public, protected, save :: c0_etot            = .false.
  logical, public, protected, save :: c0_fermisurface    = .false.
  logical, public, protected, save :: c0_fullmesh        = .false.
  logical, public, protected, save :: c0_fullstdo        = .false.
  logical, public, protected, save :: c0_geteta          = .false.
  logical, public, protected, save :: c0_getq            = .false.
  logical, public, protected, save :: c0_getwsr          = .false.
  logical, public, protected, save :: c0_gpu             = .false.
  logical, public, protected, save :: c0_gs              = .false.
  logical, public, protected, save :: c0_help            = .false.
  logical, public, protected, save :: c0_interbandonly   = .false.
  logical, public, protected, save :: c0_intrabandonly   = .false.
  logical, public, protected, save :: c0_jobgw           = .false.  ! bare --jobgw (not --jobgw=N)
  logical, public, protected, save :: c0_kchk            = .false.
  logical, public, protected, save :: c0_listcmdopt      = .false.  ! dump registered cmdopts and exit (for bash completion)
  logical, public, protected, save :: c0_mkprocar        = .false.
  logical, public, protected, save :: c0_mlo             = .false.
  logical, public, protected, save :: c0_mlo_diagnorm    = .false.
  logical, public, protected, save :: c0_mlo_feb4        = .false.
  logical, public, protected, save :: c0_mlo_ortho       = .false.
  logical, public, protected, save :: c0_mlo_orthonorm   = .false.
  logical, public, protected, save :: c0_mloahc          = .false.
  logical, public, protected, save :: c0_mlog            = .false.
  logical, public, protected, save :: c0_modifiedGS      = .false.
  logical, public, protected, save :: c0_n1n2n3eps       = .false.
  logical, public, protected, save :: c0_noinv           = .false.
  logical, public, protected, save :: c0_normcheck       = .false.
  logical, public, protected, save :: c0_nosym           = .false.
  logical, public, protected, save :: c0_nosymdm         = .false.
  logical, public, protected, save :: c0_novxc           = .false.
  logical, public, protected, save :: c0_nowritedw       = .false.
  logical, public, protected, save :: c0_ntqxx           = .false.
  logical, public, protected, save :: c0_onesp           = .false.
  logical, public, protected, save :: c0_pdos            = .false.
  logical, public, protected, save :: c0_phispinsym      = .false.
  logical, public, protected, save :: c0_q2q1test        = .false.
  logical, public, protected, save :: c0_qibzonly        = .false.
  logical, public, protected, save :: c0_quitecore       = .false.
  logical, public, protected, save :: c0_readQforGW      = .false.
  logical, public, protected, save :: c0_shorten         = .false.
  logical, public, protected, save :: c0_show_time       = .false.
  logical, public, protected, save :: c0_showdmat        = .false.
  logical, public, protected, save :: c0_skip1d          = .false.
  logical, public, protected, save :: c0_skip2nd         = .false.
  logical, public, protected, save :: c0_skip2ndd        = .false.
  logical, public, protected, save :: c0_skip2ndp        = .false.
  logical, public, protected, save :: c0_skip2nds        = .false.
  logical, public, protected, save :: c0_skipCPHI        = .false.
  logical, public, protected, save :: c0_skipGS          = .false.
  logical, public, protected, save :: c0_skip_qvalcheck  = .false.
  logical, public, protected, save :: c0_skipbstruxinit  = .false.
  logical, public, protected, save :: c0_skipd           = .false.
  logical, public, protected, save :: c0_skipf           = .false.
  logical, public, protected, save :: c0_skiphammsoc     = .false.
  logical, public, protected, save :: c0_skiplo          = .false.
  logical, public, protected, save :: c0_slat            = .false.
  logical, public, protected, save :: c0_socmatrix       = .false.
  logical, public, protected, save :: c0_tdos            = .false.
  logical, public, protected, save :: c0_tdostetf        = .false.
  logical, public, protected, save :: c0_testso          = .false.
  logical, public, protected, save :: c0_tetraw          = .false.
  logical, public, protected, save :: c0_tetwtk          = .false.
  logical, public, protected, save :: c0_use_fp32        = .false.
  logical, public, protected, save :: c0_use_gemmul8     = .false.
  logical, public, protected, save :: c0_use_sigm_fbz    = .false.
  logical, public, protected, save :: c0_v0fix           = .false.
  logical, public, protected, save :: c0_vbmonly         = .false.
  logical, public, protected, save :: c0_vesatom         = .false.
  logical, public, protected, save :: c0_vesdat          = .false.
  logical, public, protected, save :: c0_wanatom         = .false.
  logical, public, protected, save :: c0_wdsawada        = .false.
  logical, public, protected, save :: c0_wpotmt          = .false.
  logical, public, protected, save :: c0_wrhomt          = .false.
  logical, public, protected, save :: c0_writedw         = .false.
  logical, public, protected, save :: c0_writeeigen      = .false.
  logical, public, protected, save :: c0_writeham        = .false.
  logical, public, protected, save :: c0_writepdos       = .false.
  logical, public, protected, save :: c0_writesene       = .false.
  logical, public, protected, save :: c0_writev0         = .false.
  logical, public, protected, save :: c0_wsig_fbz        = .false.
  logical, public, protected, save :: c0_x0test          = .false.
  logical, public, protected, save :: c0_ylmc            = .false.
  logical, public, protected, save :: c0_zmel0           = .false.

  !========================================================================
  ! Runtime registry for typo detection. Populated as a side-effect of
  ! load_cmdopt0/2_registry; consulted by validate_arglist. Allocatable,
  ! grown one slot at a time via move_alloc -- no fixed capacity.
  !========================================================================
  character(len=24), private, save, allocatable :: known0(:)
  character(len=20), private, save, allocatable :: known2(:)

contains

  !> Append `flag` to the cmdopt0 typo-detection table. Idempotent.
  subroutine register0(flag)
    character(*), intent(in) :: flag
    character(len=24), allocatable :: tmp(:)
    integer :: n, i
    if (allocated(known0)) then
       do i = 1, size(known0)
          if (known0(i) == flag) return
       enddo
       n = size(known0)
       allocate(tmp(n + 1))
       tmp(1:n) = known0
       tmp(n + 1) = flag
       call move_alloc(from=tmp, to=known0)
    else
       allocate(known0(1))
       known0(1) = flag
    endif
  end subroutine register0

  !> Append `flag` (bare name, no trailing `=`) to the cmdopt2
  !! typo-detection table. Idempotent.
  subroutine register2(flag)
    character(*), intent(in) :: flag
    character(len=20), allocatable :: tmp(:)
    integer :: n, i
    if (allocated(known2)) then
       do i = 1, size(known2)
          if (known2(i) == flag) return
       enddo
       n = size(known2)
       allocate(tmp(n + 1))
       tmp(1:n) = known2
       tmp(n + 1) = flag
       call move_alloc(from=tmp, to=known2)
    else
       allocate(known2(1))
       known2(1) = flag
    endif
  end subroutine register2

  !> True if `flag` (full token) was registered by load_cmdopt0_registry.
  function is_known_cmdopt0(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    if (.not. allocated(known0)) return
    do i = 1, size(known0)
       if (trim(flag) == trim(known0(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt0

  !> True if `flag` (bare name, no trailing `=`) was registered by
  !! load_cmdopt2_registry.
  function is_known_cmdopt2(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    if (.not. allocated(known2)) return
    do i = 1, size(known2)
       if (trim(flag) == trim(known2(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt2

  !> Parse every cmdopt2 entry (--xxx=value form) into the cached
  !! c2_* module variables. Idempotent. Called once at startup from
  !! m_args::m_setargs after arglist is populated.
  !!
  !! narg/arglist are passed in (no `use m_args`) to keep this module
  !! self-contained and avoid a circular module-USE chain with m_args.
  !! They are forwarded explicitly to the module-level get2 helper.
  subroutine load_cmdopt2_registry(narg, arglist)
    integer,      intent(in) :: narg
    character(*), intent(in) :: arglist(narg)
    logical, save :: done = .false.
    character(len=120) :: outs
    if (done) return
    done = .true.

    ! Integer-valued
    if (get2('--jobgw', outs, narg, arglist)) then
       read(outs,*) c2_jobgw
       if (c2_jobgw /= 0 .and. c2_jobgw /= 1) &
            call rx('m_cmdopt_registry: --jobgw must be 0 or 1')
    endif
    if (get2('--job',   outs, narg, arglist)) read(outs,*) c2_job
    if (get2('--nb',    outs, narg, arglist)) read(outs,*) c2_nb
    if (get2('--nk',    outs, narg, arglist)) read(outs,*) c2_nk
    if (get2('--nww',   outs, narg, arglist)) read(outs,*) c2_nww
    if (get2('--sp1',   outs, narg, arglist)) read(outs,*) c2_sp1
    if (get2('--sp2',   outs, narg, arglist)) read(outs,*) c2_sp2
    if (get2('--ndos',  outs, narg, arglist)) read(outs,*) c2_ndos

    ! Real-valued
    if (get2('--cutuu', outs, narg, arglist)) then
       read(outs,*) c2_cutuu;          c2_cutuu_set         = .true.
    endif
    if (get2('--emin',  outs, narg, arglist)) then
       read(outs,*) c2_emin_eV;        c2_emin_set          = .true.
    endif
    if (get2('--emax',  outs, narg, arglist)) then
       read(outs,*) c2_emax_eV;        c2_emax_set          = .true.
    endif
    if (get2('--EfermiShifteV', outs, narg, arglist)) then
       read(outs,*) c2_EfermiShifteV;  c2_EfermiShifteV_set = .true.
    endif

    ! String-valued
    if (get2('--Wtype', outs, narg, arglist)) then
       c2_Wtype     = trim(outs)
       c2_Wtype_set = .true.
    endif
    ! cmdopt0-with-enum-value
    if (get2('--quit', outs, narg, arglist)) c2_quit = trim(outs)
    if (get2('--diag', outs, narg, arglist)) c2_diag = trim(outs)
    if (get2('--dwnb', outs, narg, arglist)) c2_dwnb = trim(outs)
  end subroutine load_cmdopt2_registry

  !> Look up `<flag>=<value>` in arglist; return .true. + value in
  !! outs when found. Also registers `flag` (bare name) in the cmdopt2
  !! typo-detection table.
  function get2(flag, outs, narg, arglist) result(found)
    character(*), intent(in)  :: flag
    character(*), intent(out) :: outs
    integer,      intent(in)  :: narg
    character(*), intent(in)  :: arglist(narg)
    logical :: found
    character(len=:), allocatable :: pref
    integer :: iarg, plen
    call register2(flag)
    pref = trim(flag) // '='
    plen = len(pref)
    found = .false.
    do iarg = 1, narg
       if (len_trim(arglist(iarg)) < plen) cycle
       if (arglist(iarg)(1:plen) == pref) then
          outs  = arglist(iarg)(plen + 1:)
          found = .true.
          return
       endif
    enddo
  end function get2

  !> Parse every cmdopt0 entry (pure flag form) into the cached c0_*
  !! module variables. Idempotent. Called once from m_args::m_setargs
  !! alongside load_cmdopt2_registry. Each `call set0()` also registers
  !! the flag in the typo-detection table, so adding a new cmdopt0
  !! only requires ONE new line here.
  subroutine load_cmdopt0_registry(narg, arglist)
    integer,      intent(in) :: narg
    character(*), intent(in) :: arglist(narg)
    logical, save :: done = .false.
    if (done) return
    done = .true.
    call set0('--AHCMAT',         c0_AHCMAT, narg, arglist)
    call set0('--UUMAT',          c0_UUMAT, narg, arglist)
    call set0('--afsym',          c0_afsym, narg, arglist)
    call set0('--ahc',            c0_ahc, narg, arglist)
    call set0('--allband',        c0_allband, narg, arglist)
    call set0('--avoidgamma',     c0_avoidgamma, narg, arglist)
    call set0('--band',           c0_band, narg, arglist)
    call set0('--boltztrap',      c0_boltztrap, narg, arglist)
    call set0('--cls',            c0_cls, narg, arglist)
    call set0('--cmlo',           c0_cmlo, narg, arglist)
    call set0('--corehole',       c0_corehole, narg, arglist)
    call set0('--cvK:',           c0_cvK, narg, arglist)
    call set0('--debug',          c0_debug, narg, arglist)
    call set0('--debugbndfp',     c0_debugbndfp, narg, arglist)
    call set0('--debugpwmat',     c0_debugpwmat, narg, arglist)
    call set0('--debugsugw',      c0_debugsugw, narg, arglist)
    call set0('--debugzmel',      c0_debugzmel, narg, arglist)
    call set0('--density',        c0_density, narg, arglist)
    call set0('--dos',            c0_dos, narg, arglist)
    call set0('--eigen-at-k',     c0_eigen_at_k, narg, arglist)
    call set0('--espot',          c0_espot, narg, arglist)
    call set0('--estaticall',     c0_estaticall, narg, arglist)
    call set0('--eszero',         c0_eszero, narg, arglist)
    call set0('--etot',           c0_etot, narg, arglist)
    call set0('--fermisurface',   c0_fermisurface, narg, arglist)
    call set0('--fullmesh',       c0_fullmesh, narg, arglist)
    call set0('--fullstdo',       c0_fullstdo, narg, arglist)
    call set0('--geteta',         c0_geteta, narg, arglist)
    call set0('--getq',           c0_getq, narg, arglist)
    call set0('--getwsr',         c0_getwsr, narg, arglist)
    call set0('--gpu',            c0_gpu, narg, arglist)
    call set0('--gs',             c0_gs, narg, arglist)
    call set0('--help',           c0_help, narg, arglist)
    call set0('--interbandonly',  c0_interbandonly, narg, arglist)
    call set0('--intrabandonly',  c0_intrabandonly, narg, arglist)
    call set0('--jobgw',          c0_jobgw, narg, arglist)
    call set0('--kchk',           c0_kchk, narg, arglist)
    call set0('--listcmdopt',     c0_listcmdopt, narg, arglist)
    call set0('--mkprocar',       c0_mkprocar, narg, arglist)
    call set0('--mlo',            c0_mlo, narg, arglist)
    call set0('--mlo_diagnorm',   c0_mlo_diagnorm, narg, arglist)
    call set0('--mlo_feb4',       c0_mlo_feb4, narg, arglist)
    call set0('--mlo_ortho',      c0_mlo_ortho, narg, arglist)
    call set0('--mlo_orthonorm',  c0_mlo_orthonorm, narg, arglist)
    call set0('--mloahc',         c0_mloahc, narg, arglist)
    call set0('--mlog',           c0_mlog, narg, arglist)
    call set0('--modifiedGS',     c0_modifiedGS, narg, arglist)
    call set0('--n1n2n3eps',      c0_n1n2n3eps, narg, arglist)
    call set0('--noinv',          c0_noinv, narg, arglist)
    call set0('--normcheck',      c0_normcheck, narg, arglist)
    call set0('--nosym',          c0_nosym, narg, arglist)
    call set0('--nosymdm',        c0_nosymdm, narg, arglist)
    call set0('--novxc',          c0_novxc, narg, arglist)
    call set0('--nowritedw',      c0_nowritedw, narg, arglist)
    call set0('--ntqxx',          c0_ntqxx, narg, arglist)
    call set0('--onesp',          c0_onesp, narg, arglist)
    call set0('--pdos',           c0_pdos, narg, arglist)
    call set0('--phispinsym',     c0_phispinsym, narg, arglist)
    call set0('--q2q1test',       c0_q2q1test, narg, arglist)
    call set0('--qibzonly',       c0_qibzonly, narg, arglist)
    call set0('--quitecore',      c0_quitecore, narg, arglist)
    call set0('--readQforGW',     c0_readQforGW, narg, arglist)
    call set0('--shorten',        c0_shorten, narg, arglist)
    call set0('--show_time',      c0_show_time, narg, arglist)
    call set0('--showdmat',       c0_showdmat, narg, arglist)
    call set0('--skip1d',         c0_skip1d, narg, arglist)
    call set0('--skip2nd',        c0_skip2nd, narg, arglist)
    call set0('--skip2ndd',       c0_skip2ndd, narg, arglist)
    call set0('--skip2ndp',       c0_skip2ndp, narg, arglist)
    call set0('--skip2nds',       c0_skip2nds, narg, arglist)
    call set0('--skipCPHI',       c0_skipCPHI, narg, arglist)
    call set0('--skipGS',         c0_skipGS, narg, arglist)
    call set0('--skip_qvalcheck', c0_skip_qvalcheck, narg, arglist)
    call set0('--skipbstruxinit', c0_skipbstruxinit, narg, arglist)
    call set0('--skipd',          c0_skipd, narg, arglist)
    call set0('--skipf',          c0_skipf, narg, arglist)
    call set0('--skiphammsoc',    c0_skiphammsoc, narg, arglist)
    call set0('--skiplo',         c0_skiplo, narg, arglist)
    call set0('--slat',           c0_slat, narg, arglist)
    call set0('--socmatrix',      c0_socmatrix, narg, arglist)
    call set0('--tdos',           c0_tdos, narg, arglist)
    call set0('--tdostetf',       c0_tdostetf, narg, arglist)
    call set0('--testso',         c0_testso, narg, arglist)
    call set0('--tetraw',         c0_tetraw, narg, arglist)
    call set0('--tetwtk',         c0_tetwtk, narg, arglist)
    call set0('--use_fp32',       c0_use_fp32, narg, arglist)
    call set0('--use_gemmul8',    c0_use_gemmul8, narg, arglist)
    call set0('--use_sigm_fbz',   c0_use_sigm_fbz, narg, arglist)
    call set0('--v0fix',          c0_v0fix, narg, arglist)
    call set0('--vbmonly',        c0_vbmonly, narg, arglist)
    call set0('--vesatom',        c0_vesatom, narg, arglist)
    call set0('--vesdat',         c0_vesdat, narg, arglist)
    call set0('--wanatom',        c0_wanatom, narg, arglist)
    call set0('--wdsawada',       c0_wdsawada, narg, arglist)
    call set0('--wpotmt',         c0_wpotmt, narg, arglist)
    call set0('--wrhomt',         c0_wrhomt, narg, arglist)
    call set0('--writedw',        c0_writedw, narg, arglist)
    call set0('--writeeigen',     c0_writeeigen, narg, arglist)
    call set0('--writeham',       c0_writeham, narg, arglist)
    call set0('--writepdos',      c0_writepdos, narg, arglist)
    call set0('--writesene',      c0_writesene, narg, arglist)
    call set0('--writev0',        c0_writev0, narg, arglist)
    call set0('--wsig_fbz',       c0_wsig_fbz, narg, arglist)
    call set0('--x0test',         c0_x0test, narg, arglist)
    call set0('--ylmc',           c0_ylmc, narg, arglist)
    call set0('--zmel0',          c0_zmel0, narg, arglist)
  end subroutine load_cmdopt0_registry

  !> Look up `flag` verbatim in arglist; set `var` to the result AND
  !! register `flag` in the cmdopt0 typo-detection table.
  subroutine set0(flag, var, narg, arglist)
    character(*), intent(in)  :: flag
    logical,      intent(out) :: var
    integer,      intent(in)  :: narg
    character(*), intent(in)  :: arglist(narg)
    integer :: iarg
    call register0(flag)
    var = .false.
    do iarg = 1, narg
       if (trim(arglist(iarg)) == flag) then
          var = .true.
          return
       endif
    enddo
  end subroutine set0

  !> Print every registered cmdopt0 / cmdopt2 name, one per line, to
  !! stdout, then exit cleanly. cmdopt2 names are emitted with a
  !! trailing `=` so a bash completion script can distinguish them
  !! from bare flags. Used by `<binary> --listcmdopt` to feed shell
  !! completion (`compgen -W "$(lmf --listcmdopt)"`).
  !!
  !! Triggered from main_lmf / main_lmfa / main_lmchk after m_setargs
  !! populates the runtime tables and before MPI is initialized; we
  !! call `exit` directly so MPI need not have been spun up.
  subroutine list_cmdopts()
    use mpi
    integer :: i, ierr
    logical :: mpi_was_init
    ! Some MPI launchers (e.g. HPCX on kt1) flag an "abnormal
    ! termination" if a process exits without going through
    ! MPI_Init/MPI_Finalize, even when its rank's exit code is 0.
    ! Bracket the dump with MPI init/finalize so `mpirun -np 1 lmf
    ! --listcmdopt` (the path InstallAll.py uses to generate the
    ! completion data) returns 0 cleanly.
    call MPI_Initialized(mpi_was_init, ierr)
    if (.not. mpi_was_init) call MPI_Init(ierr)
    if (allocated(known0)) then
       do i = 1, size(known0)
          write(*,'(a)') trim(known0(i))
       enddo
    endif
    if (allocated(known2)) then
       do i = 1, size(known2)
          write(*,'(a)') trim(known2(i)) // '='
       enddo
    endif
    flush(6)
    if (.not. mpi_was_init) call MPI_Finalize(ierr)
    call exit(0)
  end subroutine list_cmdopts

  !> Final-pass strict typo / retired-syntax check. Walks every token
  !! and aborts on:
  !!   - retired forms (-v...=, --[...]=, --toml.<path>=, --pr=, etc.)
  !!     with a one-line migration hint
  !!   - any --word / -word (no `=`) not in the cmdopt0 table
  !!   - any --word=... whose bare prefix isn't in the cmdopt2 table
  !! Skips:
  !!   - bare positional args (no leading `-`)
  !!   - --ctrlg:<path>=<value> (canonical TOML override syntax)
  !!
  !! narg/arglist are passed in (rather than `use m_args, only: ...`) to
  !! avoid a circular module dependency: m_args::m_setargs already uses
  !! m_cmdopt_registry to call this routine.
  subroutine validate_arglist(narg, arglist)
    integer,      intent(in) :: narg
    character(*), intent(in) :: arglist(narg)
    logical, save :: done = .false.
    integer :: i, alen, eq
    character(len=:), allocatable :: arg
    if (done) return
    done = .true.
    do i = 1, narg
       arg = trim(arglist(i))
       alen = len(arg)
       if (alen == 0) cycle
       if (arg(1:1) /= '-') cycle                         ! positional
       if (alen >= 8 .and. arg(1:8) == '--ctrlg:') cycle  ! TOML override
       call check_retired(arg)                            ! aborts on legacy forms
       eq = index(arg, '=')
       if (eq == 0) then
          if (is_known_cmdopt0(arg)) cycle
       else
          if (is_known_cmdopt2(arg(1:eq-1))) cycle
       endif
       call die_unknown(arg)
    enddo
  end subroutine validate_arglist

  subroutine check_retired(arg)
    character(*), intent(in) :: arg
    integer :: alen
    alen = len_trim(arg)
    if (alen < 2) return
    if (alen >= 3 .and. arg(1:2) == '-v' .and. index(arg, '=') > 0) &
         call die_retired(arg, &
              '-v...=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
    if (alen >= 5 .and. arg(1:3) == '--[' .and. index(arg, ']=') > 0) &
         call die_retired(arg, &
              '--[<path>]=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
    if (alen >= 9 .and. arg(1:7) == '--toml.' .and. index(arg, '=') > 0) &
         call die_retired(arg, &
              '--toml.<path>=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
    if (alen >= 5 .and. arg(1:5) == '--pr=') &
         call die_retired(arg, '--pr=N is retired. Use --ctrlg:verbose=N.')
    if (alen >= 7 .and. arg(1:7) == '--time=') &
         call die_retired(arg, '--time=<...> is retired. Use --ctrlg:time=[N,M].')
    if (arg(1:alen) == '--phispinsym') &
         call die_retired(arg, &
              '--phispinsym is retired. Use --ctrlg:ham.phispinsym=true, '// &
              'or set [ham] phispinsym = true in ctrlg.<sname>.toml.')
    ! Single-dash legacy variants of cmdopt2 entries -- now double-dash only.
    if (alen >= 6 .and. arg(1:6) == '-emin=') &
         call die_retired(arg, '-emin=<value> is retired. Use --emin=<value> (double dash).')
    if (alen >= 6 .and. arg(1:6) == '-emax=') &
         call die_retired(arg, '-emax=<value> is retired. Use --emax=<value> (double dash).')
    if (alen >= 6 .and. arg(1:6) == '-ndos=') &
         call die_retired(arg, '-ndos=<value> is retired. Use --ndos=<value> (double dash).')
    if (alen >= 15 .and. arg(1:15) == '-EfermiShifteV=') &
         call die_retired(arg, '-EfermiShifteV=<value> is retired. Use --EfermiShifteV=<value> (double dash).')
  end subroutine check_retired

  subroutine die_unknown(arg)
    character(*), intent(in) :: arg
    write(*,'(a)') ' '
    write(*,'(a)') 'ERROR: unknown option `'//trim(arg)//'`.'
    write(*,'(a)') '       To override a TOML key, use --ctrlg:<dotted.path>=<value>.'
    write(*,'(a)') '       Registered options: see m_cmdopt_registry.f90 (c0_*/c2_*)'
    write(*,'(a)') '       or https://ecalj.github.io/ecaljdoc/manual/cmdopts'
    call rx('unknown option: '//trim(arg))
  end subroutine die_unknown

  subroutine die_retired(arg, hint)
    character(*), intent(in) :: arg, hint
    write(*,'(a)') ' '
    write(*,'(a)') 'ERROR: retired syntax `'//trim(arg)//'`'
    write(*,'(a)') '       '//trim(hint)
    write(*,'(a)') '       https://ecalj.github.io/ecaljdoc/manual/toml_migration'
    call rx('retired syntax: '//trim(arg))
  end subroutine die_retired

end module m_cmdopt_registry
