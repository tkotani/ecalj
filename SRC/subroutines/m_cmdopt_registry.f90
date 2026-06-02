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
!!   - When adding a new cmdopt0 (pure flag), append to CMDOPT0_REGISTRY
!!     AND declare a `c0_<name>` logical AND set it in load_cmdopt0_registry.
!!   - When adding a new cmdopt2 (=value), append to CMDOPT2_REGISTRY AND
!!     declare a `c2_<name>` of appropriate type AND parse it in
!!     load_cmdopt2_registry.
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
  ! Registry tables (single source of truth for typo detection).
  !========================================================================

  ! Pure cmdopt0 flag names. Strict typo check: if arglist has a `--word`
  ! or `-word` (no `=`), and it isn't in this list, abort.
  character(len=24), parameter :: CMDOPT0_REGISTRY(*) = [character(len=24) :: &
       '--AHCMAT',          &
       '--UUMAT',           &
       '--afsym',           &
       '--ahc',             &
       '--allband',         &
       '--avoidgamma',      &
       '--band',            &
       '--boltztrap',       &
       '--cls',             &
       '--cmlo',            &
       '--corehole',        &
       '--cvK:',            &
       '--debug',           &
       '--debugbndfp',      &
       '--debugpwmat',      &
       '--debugsugw',       &
       '--debugzmel',       &
       '--density',         &
       '--dos',             &
       '--eigen-at-k',      &
       '--espot',           &
       '--estaticall',      &
       '--eszero',          &
       '--etot',            &
       '--fermisurface',    &
       '--fullmesh',        &
       '--fullstdo',        &
       '--geteta',          &
       '--getq',            &
       '--getwsr',          &
       '--gpu',             &
       '--gs',              &
       '--help',            &
       '--interbandonly',   &
       '--intrabandonly',   &
       '--jobgw',           &
       '--kchk',            &
       '--mkprocar',        &
       '--mlo',             &
       '--mlo_diagnorm',    &
       '--mlo_feb4',        &
       '--mlo_ortho',       &
       '--mlo_orthonorm',   &
       '--mloahc',          &
       '--mlog',            &
       '--modifiedGS',      &
       '--n1n2n3eps',       &
       '--noinv',           &
       '--normcheck',       &
       '--nosym',           &
       '--nosymdm',         &
       '--novxc',           &
       '--nowritedw',       &
       '--ntqxx',           &
       '--onesp',           &
       '--pdos',            &
       '--phispinsym',      &
       '--q2q1test',        &
       '--qibzonly',        &
       '--quitecore',       &
       '--readQforGW',      &
       '--shorten',         &
       '--show_time',       &
       '--showdmat',        &
       '--skip1d',          &
       '--skip2nd',         &
       '--skip2ndd',        &
       '--skip2ndp',        &
       '--skip2nds',        &
       '--skipCPHI',        &
       '--skipGS',          &
       '--skip_qvalcheck',  &
       '--skipbstruxinit',  &
       '--skipd',           &
       '--skipf',           &
       '--skiphammsoc',     &
       '--skiplo',          &
       '--slat',            &
       '--socmatrix',       &
       '--tdos',            &
       '--tdostetf',        &
       '--terse',           &
       '-terse',            &  ! single-dash alias retained for lmaux
       '--testso',          &
       '--tetraw',          &
       '--tetwtk',          &
       '--use_gemmul8',     &
       '--use_sigm_fbz',    &
       '--v0fix',           &
       '--vbmonly',         &
       '--vesatom',         &
       '--vesdat',          &
       '--wanatom',         &
       '--wdsawada',        &
       '--wpotmt',          &
       '--wrhomt',          &
       '--writedw',         &
       '--writeeigen',      &
       '--writeham',        &
       '--writepdos',       &
       '--writesene',       &
       '--writev0',         &
       '--wsig_fbz',        &
       '--x0test',          &
       '--ylmc',            &
       '--zmel0'            ]

  ! cmdopt2 names (=value form). Drop the trailing `=`; the matcher in
  ! m_toml_override compares `--<word>` after splitting on `=`.
  !
  ! Three of these (`--quit`, `--diag`, `--dwnb`) are conceptually
  ! cmdopt0-with-enum-value; we cache them as cmdopt2 strings (c2_quit
  ! etc.) and rewrite use sites as `if (c2_quit == 'show')` etc.
  character(len=20), parameter :: CMDOPT2_REGISTRY(*) = [character(len=20) :: &
       '--Wtype',           &
       '--cutuu',           &
       '--diag',            &
       '--dwnb',            &
       '--job',             &
       '--jobgw',           &
       '--nb',              &
       '--nk',              &
       '--nww',             &
       '--quit',            &
       '--sp1',             &
       '--sp2',             &
       '-EfermiShifteV',    &
       '-emax',             &
       '-emin',             &
       '-ndos'              ]

contains

  !> True if `flag` ("--foo" or "-foo", *without* trailing `=`) is a
  !! registered cmdopt2 name. Case-sensitive.
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

  !> True if `flag` (full token, e.g. "--writeham" or "-terse") is a
  !! registered cmdopt0 flag. Case-sensitive.
  function is_known_cmdopt0(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    do i = 1, size(CMDOPT0_REGISTRY)
       if (trim(flag) == trim(CMDOPT0_REGISTRY(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt0

  !> Parse every cmdopt2 entry into the cached module variables above.
  !! Idempotent: subsequent calls are no-ops. Called once at startup
  !! from m_args::m_setargs (after arglist is populated).
  subroutine load_cmdopt2_registry()
    logical, save :: done = .false.
    character(len=120) :: outs
    logical, external :: cmdopt2
    if (done) return
    done = .true.

    ! Integer-valued
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

    ! Real-valued
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

    ! String-valued
    if (cmdopt2('--Wtype=', outs)) then
       c2_Wtype     = trim(outs)
       c2_Wtype_set = .true.
    endif
    ! cmdopt0-with-enum-value
    if (cmdopt2('--quit=', outs)) c2_quit = trim(outs)
    if (cmdopt2('--diag=', outs)) c2_diag = trim(outs)
    if (cmdopt2('--dwnb=', outs)) c2_dwnb = trim(outs)
  end subroutine load_cmdopt2_registry

  !> Parse every cmdopt0 entry (pure flags) into the cached c0_* logicals.
  !! Idempotent. Called once from m_args::m_setargs alongside
  !! load_cmdopt2_registry.
  subroutine load_cmdopt0_registry()
    logical, save :: done = .false.
    logical, external :: cmdopt0
    if (done) return
    done = .true.
    c0_AHCMAT          = cmdopt0('--AHCMAT')
    c0_UUMAT           = cmdopt0('--UUMAT')
    c0_afsym           = cmdopt0('--afsym')
    c0_ahc             = cmdopt0('--ahc')
    c0_allband         = cmdopt0('--allband')
    c0_avoidgamma      = cmdopt0('--avoidgamma')
    c0_band            = cmdopt0('--band')
    c0_boltztrap       = cmdopt0('--boltztrap')
    c0_cls             = cmdopt0('--cls')
    c0_cmlo            = cmdopt0('--cmlo')
    c0_corehole        = cmdopt0('--corehole')
    c0_cvK             = cmdopt0('--cvK:')
    c0_debug           = cmdopt0('--debug')
    c0_debugbndfp      = cmdopt0('--debugbndfp')
    c0_debugpwmat      = cmdopt0('--debugpwmat')
    c0_debugsugw       = cmdopt0('--debugsugw')
    c0_debugzmel       = cmdopt0('--debugzmel')
    c0_density         = cmdopt0('--density')
    c0_dos             = cmdopt0('--dos')
    c0_eigen_at_k      = cmdopt0('--eigen-at-k')
    c0_espot           = cmdopt0('--espot')
    c0_estaticall      = cmdopt0('--estaticall')
    c0_eszero          = cmdopt0('--eszero')
    c0_etot            = cmdopt0('--etot')
    c0_fermisurface    = cmdopt0('--fermisurface')
    c0_fullmesh        = cmdopt0('--fullmesh')
    c0_fullstdo        = cmdopt0('--fullstdo')
    c0_geteta          = cmdopt0('--geteta')
    c0_getq            = cmdopt0('--getq')
    c0_getwsr          = cmdopt0('--getwsr')
    c0_gpu             = cmdopt0('--gpu')
    c0_gs              = cmdopt0('--gs')
    c0_help            = cmdopt0('--help')
    c0_interbandonly   = cmdopt0('--interbandonly')
    c0_intrabandonly   = cmdopt0('--intrabandonly')
    c0_jobgw           = cmdopt0('--jobgw')
    c0_kchk            = cmdopt0('--kchk')
    c0_mkprocar        = cmdopt0('--mkprocar')
    c0_mlo             = cmdopt0('--mlo')
    c0_mlo_diagnorm    = cmdopt0('--mlo_diagnorm')
    c0_mlo_feb4        = cmdopt0('--mlo_feb4')
    c0_mlo_ortho       = cmdopt0('--mlo_ortho')
    c0_mlo_orthonorm   = cmdopt0('--mlo_orthonorm')
    c0_mloahc          = cmdopt0('--mloahc')
    c0_mlog            = cmdopt0('--mlog')
    c0_modifiedGS      = cmdopt0('--modifiedGS')
    c0_n1n2n3eps       = cmdopt0('--n1n2n3eps')
    c0_noinv           = cmdopt0('--noinv')
    c0_normcheck       = cmdopt0('--normcheck')
    c0_nosym           = cmdopt0('--nosym')
    c0_nosymdm         = cmdopt0('--nosymdm')
    c0_novxc           = cmdopt0('--novxc')
    c0_nowritedw       = cmdopt0('--nowritedw')
    c0_ntqxx           = cmdopt0('--ntqxx')
    c0_onesp           = cmdopt0('--onesp')
    c0_pdos            = cmdopt0('--pdos')
    c0_phispinsym      = cmdopt0('--phispinsym')
    c0_q2q1test        = cmdopt0('--q2q1test')
    c0_qibzonly        = cmdopt0('--qibzonly')
    c0_quitecore       = cmdopt0('--quitecore')
    c0_readQforGW      = cmdopt0('--readQforGW')
    c0_shorten         = cmdopt0('--shorten')
    c0_show_time       = cmdopt0('--show_time')
    c0_showdmat        = cmdopt0('--showdmat')
    c0_skip1d          = cmdopt0('--skip1d')
    c0_skip2nd         = cmdopt0('--skip2nd')
    c0_skip2ndd        = cmdopt0('--skip2ndd')
    c0_skip2ndp        = cmdopt0('--skip2ndp')
    c0_skip2nds        = cmdopt0('--skip2nds')
    c0_skipCPHI        = cmdopt0('--skipCPHI')
    c0_skipGS          = cmdopt0('--skipGS')
    c0_skip_qvalcheck  = cmdopt0('--skip_qvalcheck')
    c0_skipbstruxinit  = cmdopt0('--skipbstruxinit')
    c0_skipd           = cmdopt0('--skipd')
    c0_skipf           = cmdopt0('--skipf')
    c0_skiphammsoc     = cmdopt0('--skiphammsoc')
    c0_skiplo          = cmdopt0('--skiplo')
    c0_slat            = cmdopt0('--slat')
    c0_socmatrix       = cmdopt0('--socmatrix')
    c0_tdos            = cmdopt0('--tdos')
    c0_tdostetf        = cmdopt0('--tdostetf')
    c0_terse           = cmdopt0('--terse')
    c0_terse_short     = cmdopt0('-terse')
    c0_testso          = cmdopt0('--testso')
    c0_tetraw          = cmdopt0('--tetraw')
    c0_tetwtk          = cmdopt0('--tetwtk')
    c0_use_gemmul8     = cmdopt0('--use_gemmul8')
    c0_use_sigm_fbz    = cmdopt0('--use_sigm_fbz')
    c0_v0fix           = cmdopt0('--v0fix')
    c0_vbmonly         = cmdopt0('--vbmonly')
    c0_vesatom         = cmdopt0('--vesatom')
    c0_vesdat          = cmdopt0('--vesdat')
    c0_wanatom         = cmdopt0('--wanatom')
    c0_wdsawada        = cmdopt0('--wdsawada')
    c0_wpotmt          = cmdopt0('--wpotmt')
    c0_wrhomt          = cmdopt0('--wrhomt')
    c0_writedw         = cmdopt0('--writedw')
    c0_writeeigen      = cmdopt0('--writeeigen')
    c0_writeham        = cmdopt0('--writeham')
    c0_writepdos       = cmdopt0('--writepdos')
    c0_writesene       = cmdopt0('--writesene')
    c0_writev0         = cmdopt0('--writev0')
    c0_wsig_fbz        = cmdopt0('--wsig_fbz')
    c0_x0test          = cmdopt0('--x0test')
    c0_ylmc            = cmdopt0('--ylmc')
    c0_zmel0           = cmdopt0('--zmel0')
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
