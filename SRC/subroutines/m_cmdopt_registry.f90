!> Central registry of command-line option lookups. Two roles in one module:
!!
!!   1. Single-shot parsing into typed module variables. Callsites do
!!      `use m_cmdopt_registry, only: c0_writeham` (bool) or
!!      `use m_cmdopt_registry, only: c2_jobgw` (int) and read the cached
!!      value directly instead of the old `cmdopt0(...) / cmdopt2(...);
!!      read(outs,*) jobgw` pattern. load_*_registry() are called once at
!!      startup (from m_args::m_setargs); subsequent calls are no-ops.
!!
!!   2. Strict typo detection. classify_non_toml_arg() in m_toml_override
!!      consults is_known_cmdopt0()/is_known_cmdopt2() before deciding
!!      between "TOML override", "registered cmdopt, pass through", and
!!      "user typo, abort with hint".
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
  public :: is_known_cmdopt0
  public :: is_known_cmdopt2
  public :: load_cmdopt2_registry
  public :: load_cmdopt0_registry

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
  logical, public, protected, save :: c0_terse           = .false.  ! --terse
  logical, public, protected, save :: c0_terse_short     = .false.  ! -terse (single dash alias)
  logical, public, protected, save :: c0_testso          = .false.
  logical, public, protected, save :: c0_tetraw          = .false.
  logical, public, protected, save :: c0_tetwtk          = .false.
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
  ! Runtime registry for typo detection.
  !
  ! Populated as a side-effect of load_cmdopt0/2_registry: each flag we
  ! look up is also appended to one of these tables. is_known_cmdopt0/2
  ! linearly scan until the first empty slot. No separate counter or
  ! parameter array to keep in sync -- the load functions ARE the
  ! source of truth, and the array itself records how many entries
  ! exist (the populated prefix ends at the first '' slot).
  !
  ! Capacity is a fixed upper bound (~9 KB total) deliberately well
  ! above the present 105+16 entries.
  !========================================================================
  integer, parameter :: REG_CAP = 200
  character(len=24), private, save :: known0(REG_CAP) = ''
  character(len=20), private, save :: known2(REG_CAP) = ''

contains

  !> Append `flag` to the cmdopt0 typo-detection table. Idempotent.
  !! Scans for the first empty slot; if we hit REG_CAP without finding
  !! one, abort -- means someone added > REG_CAP flags and the cap
  !! needs bumping.
  subroutine register0(flag)
    character(*), intent(in) :: flag
    integer :: i
    do i = 1, REG_CAP
       if (known0(i) == flag) return         ! already registered
       if (known0(i) == '') then             ! first empty slot
          known0(i) = flag
          return
       endif
    enddo
    call rx('m_cmdopt_registry: bump REG_CAP')
  end subroutine register0

  !> Append `flag` (bare name, no trailing `=`) to the cmdopt2
  !! typo-detection table. Idempotent.
  subroutine register2(flag)
    character(*), intent(in) :: flag
    integer :: i
    do i = 1, REG_CAP
       if (known2(i) == flag) return
       if (known2(i) == '') then
          known2(i) = flag
          return
       endif
    enddo
    call rx('m_cmdopt_registry: bump REG_CAP')
  end subroutine register2

  !> True if `flag` (full token, e.g. "--writeham" or "-terse") was
  !! registered by load_cmdopt0_registry. Case-sensitive.
  function is_known_cmdopt0(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    do i = 1, REG_CAP
       if (known0(i) == '') return           ! end of populated entries
       if (trim(flag) == trim(known0(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt0

  !> True if `flag` ("--foo" or "-foo", *without* trailing `=`) was
  !! registered by load_cmdopt2_registry. Case-sensitive.
  function is_known_cmdopt2(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    do i = 1, REG_CAP
       if (known2(i) == '') return
       if (trim(flag) == trim(known2(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt2

  !> Look up a pure cmdopt0 flag in arglist, cache the boolean, AND
  !! register the flag in the runtime typo-detection table.
  subroutine set0(flag, var)
    character(*), intent(in)  :: flag
    logical,      intent(out) :: var
    logical, external :: cmdopt0
    var = cmdopt0(flag)
    call register0(flag)
  end subroutine set0

  !> Wrapper around cmdopt2 that also registers `flag` (bare name; the
  !! trailing `=` is appended internally before the cmdopt2 call) in
  !! the cmdopt2 typo-detection table.
  function get2(flag, outs) result(found)
    character(*), intent(in)  :: flag
    character(*), intent(out) :: outs
    logical :: found
    logical, external :: cmdopt2
    call register2(flag)
    found = cmdopt2(trim(flag) // '=', outs)
  end function get2

  !> Parse every cmdopt2 entry into the cached module variables above.
  !! Idempotent: subsequent calls are no-ops. Called once at startup
  !! from m_args::m_setargs (after arglist is populated).
  subroutine load_cmdopt2_registry()
    logical, save :: done = .false.
    character(len=120) :: outs
    if (done) return
    done = .true.

    ! Integer-valued
    if (get2('--jobgw', outs)) then
       read(outs,*) c2_jobgw
       if (c2_jobgw /= 0 .and. c2_jobgw /= 1) &
            call rx('m_cmdopt_registry: --jobgw must be 0 or 1')
    endif
    if (get2('--job',   outs)) read(outs,*) c2_job
    if (get2('--nb',    outs)) read(outs,*) c2_nb
    if (get2('--nk',    outs)) read(outs,*) c2_nk
    if (get2('--nww',   outs)) read(outs,*) c2_nww
    if (get2('--sp1',   outs)) read(outs,*) c2_sp1
    if (get2('--sp2',   outs)) read(outs,*) c2_sp2
    if (get2('-ndos',   outs)) read(outs,*) c2_ndos

    ! Real-valued
    if (get2('--cutuu', outs)) then
       read(outs,*) c2_cutuu;          c2_cutuu_set         = .true.
    endif
    if (get2('-emin',   outs)) then
       read(outs,*) c2_emin_eV;        c2_emin_set          = .true.
    endif
    if (get2('-emax',   outs)) then
       read(outs,*) c2_emax_eV;        c2_emax_set          = .true.
    endif
    if (get2('-EfermiShifteV', outs)) then
       read(outs,*) c2_EfermiShifteV;  c2_EfermiShifteV_set = .true.
    endif

    ! String-valued
    if (get2('--Wtype', outs)) then
       c2_Wtype     = trim(outs)
       c2_Wtype_set = .true.
    endif
    ! cmdopt0-with-enum-value
    if (get2('--quit', outs)) c2_quit = trim(outs)
    if (get2('--diag', outs)) c2_diag = trim(outs)
    if (get2('--dwnb', outs)) c2_dwnb = trim(outs)
  end subroutine load_cmdopt2_registry

  !> Parse every cmdopt0 entry (pure flags) into the cached c0_* logicals.
  !! Idempotent. Called once from m_args::m_setargs alongside
  !! load_cmdopt2_registry. Each call to set0() also registers the flag
  !! in the typo-detection table, so adding a new cmdopt0 only requires
  !! ONE new line here.
  subroutine load_cmdopt0_registry()
    logical, save :: done = .false.
    if (done) return
    done = .true.
    call set0('--AHCMAT',         c0_AHCMAT)
    call set0('--UUMAT',          c0_UUMAT)
    call set0('--afsym',          c0_afsym)
    call set0('--ahc',            c0_ahc)
    call set0('--allband',        c0_allband)
    call set0('--avoidgamma',     c0_avoidgamma)
    call set0('--band',           c0_band)
    call set0('--boltztrap',      c0_boltztrap)
    call set0('--cls',            c0_cls)
    call set0('--cmlo',           c0_cmlo)
    call set0('--corehole',       c0_corehole)
    call set0('--cvK:',           c0_cvK)
    call set0('--debug',          c0_debug)
    call set0('--debugbndfp',     c0_debugbndfp)
    call set0('--debugpwmat',     c0_debugpwmat)
    call set0('--debugsugw',      c0_debugsugw)
    call set0('--debugzmel',      c0_debugzmel)
    call set0('--density',        c0_density)
    call set0('--dos',            c0_dos)
    call set0('--eigen-at-k',     c0_eigen_at_k)
    call set0('--espot',          c0_espot)
    call set0('--estaticall',     c0_estaticall)
    call set0('--eszero',         c0_eszero)
    call set0('--etot',           c0_etot)
    call set0('--fermisurface',   c0_fermisurface)
    call set0('--fullmesh',       c0_fullmesh)
    call set0('--fullstdo',       c0_fullstdo)
    call set0('--geteta',         c0_geteta)
    call set0('--getq',           c0_getq)
    call set0('--getwsr',         c0_getwsr)
    call set0('--gpu',            c0_gpu)
    call set0('--gs',             c0_gs)
    call set0('--help',           c0_help)
    call set0('--interbandonly',  c0_interbandonly)
    call set0('--intrabandonly',  c0_intrabandonly)
    call set0('--jobgw',          c0_jobgw)
    call set0('--kchk',           c0_kchk)
    call set0('--mkprocar',       c0_mkprocar)
    call set0('--mlo',            c0_mlo)
    call set0('--mlo_diagnorm',   c0_mlo_diagnorm)
    call set0('--mlo_feb4',       c0_mlo_feb4)
    call set0('--mlo_ortho',      c0_mlo_ortho)
    call set0('--mlo_orthonorm',  c0_mlo_orthonorm)
    call set0('--mloahc',         c0_mloahc)
    call set0('--mlog',           c0_mlog)
    call set0('--modifiedGS',     c0_modifiedGS)
    call set0('--n1n2n3eps',      c0_n1n2n3eps)
    call set0('--noinv',          c0_noinv)
    call set0('--normcheck',      c0_normcheck)
    call set0('--nosym',          c0_nosym)
    call set0('--nosymdm',        c0_nosymdm)
    call set0('--novxc',          c0_novxc)
    call set0('--nowritedw',      c0_nowritedw)
    call set0('--ntqxx',          c0_ntqxx)
    call set0('--onesp',          c0_onesp)
    call set0('--pdos',           c0_pdos)
    call set0('--phispinsym',     c0_phispinsym)
    call set0('--q2q1test',       c0_q2q1test)
    call set0('--qibzonly',       c0_qibzonly)
    call set0('--quitecore',      c0_quitecore)
    call set0('--readQforGW',     c0_readQforGW)
    call set0('--shorten',        c0_shorten)
    call set0('--show_time',      c0_show_time)
    call set0('--showdmat',       c0_showdmat)
    call set0('--skip1d',         c0_skip1d)
    call set0('--skip2nd',        c0_skip2nd)
    call set0('--skip2ndd',       c0_skip2ndd)
    call set0('--skip2ndp',       c0_skip2ndp)
    call set0('--skip2nds',       c0_skip2nds)
    call set0('--skipCPHI',       c0_skipCPHI)
    call set0('--skipGS',         c0_skipGS)
    call set0('--skip_qvalcheck', c0_skip_qvalcheck)
    call set0('--skipbstruxinit', c0_skipbstruxinit)
    call set0('--skipd',          c0_skipd)
    call set0('--skipf',          c0_skipf)
    call set0('--skiphammsoc',    c0_skiphammsoc)
    call set0('--skiplo',         c0_skiplo)
    call set0('--slat',           c0_slat)
    call set0('--socmatrix',      c0_socmatrix)
    call set0('--tdos',           c0_tdos)
    call set0('--tdostetf',       c0_tdostetf)
    call set0('--terse',          c0_terse)
    call set0('-terse',           c0_terse_short)  ! single-dash alias for lmaux
    call set0('--testso',         c0_testso)
    call set0('--tetraw',         c0_tetraw)
    call set0('--tetwtk',         c0_tetwtk)
    call set0('--use_gemmul8',    c0_use_gemmul8)
    call set0('--use_sigm_fbz',   c0_use_sigm_fbz)
    call set0('--v0fix',          c0_v0fix)
    call set0('--vbmonly',        c0_vbmonly)
    call set0('--vesatom',        c0_vesatom)
    call set0('--vesdat',         c0_vesdat)
    call set0('--wanatom',        c0_wanatom)
    call set0('--wdsawada',       c0_wdsawada)
    call set0('--wpotmt',         c0_wpotmt)
    call set0('--wrhomt',         c0_wrhomt)
    call set0('--writedw',        c0_writedw)
    call set0('--writeeigen',     c0_writeeigen)
    call set0('--writeham',       c0_writeham)
    call set0('--writepdos',      c0_writepdos)
    call set0('--writesene',      c0_writesene)
    call set0('--writev0',        c0_writev0)
    call set0('--wsig_fbz',       c0_wsig_fbz)
    call set0('--x0test',         c0_x0test)
    call set0('--ylmc',           c0_ylmc)
    call set0('--zmel0',          c0_zmel0)
  end subroutine load_cmdopt0_registry

end module m_cmdopt_registry

!=========================================================================
! cmdopt0 / cmdopt2: top-level argument-lookup functions.
!
! Kept OUTSIDE m_cmdopt_registry so that m_args::m_setargs (which calls
! load_cmdopt0/2_registry on first invocation) does not introduce a
! circular module USE with m_cmdopt_registry. Module-level USE chains
! must be acyclic in Fortran, so cmdopt0/cmdopt2 live as standalone
! external functions instead -- m_cmdopt_registry's load_*_registry
! routines declare them `logical, external` and call them directly.
!
! Migrated 2026-06-02 from m_ext.f90 alongside the cmdopt0/cmdopt2
! single-shot caching refactor.
!=========================================================================

!> True if `argstr` exists verbatim in arglist. Triggers lazy m_setargs
!! init on first call; subsequent calls reuse arglist.
logical function cmdopt0(argstr)
  use m_args, only: m_setargs, arglist, narg
  implicit none
  character(*) :: argstr
  integer :: iarg
  cmdopt0 = .false.
  call m_setargs()
  do iarg = 1, narg
     if (trim(arglist(iarg)) == trim(argstr)) then
        cmdopt0 = .true.
        return
     endif
  enddo
end function cmdopt0

!> True if some arglist token starts with `argstr` (typically
!! `'--foo='`); the remainder of the matching token is returned in
!! `outstr`. Triggers lazy m_setargs init on first call.
logical function cmdopt2(argstr, outstr)
  use m_args, only: m_setargs, arglist, narg
  implicit none
  character(*) :: argstr, outstr
  integer :: iarg, strlnx
  cmdopt2 = .false.
  call m_setargs()
  do iarg = 1, narg
     strlnx = len_trim(argstr)
     if (arglist(iarg)(1:strlnx) == trim(argstr)) then
        cmdopt2 = .true.
        outstr = arglist(iarg)(strlnx + 1:)
        return
     endif
  enddo
end function cmdopt2
