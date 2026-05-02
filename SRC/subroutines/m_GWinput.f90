!> m_GWinput: single source of truth for GWinput.toml contents.
!  Loads once via toml-f, exposes all values as protected module variables.
!  Callers `use m_GWinput, only: niw, deltaw, ...` to access values.
!
!  This replaces scattered `getkeyvalue("GWinput", key, var)` calls — all
!  GWinput keys are loaded eagerly at startup, then read-only thereafter.
!
!  Schema based on survey of 66 Samples GWinput files.
!
module m_GWinput
  use tomlf, only: toml_table, toml_array, toml_load, toml_error, &
                   get_value, len
  implicit none
  private

  !-----------------------------------------------------------------
  ! [gw] section: scalar / vector keys (defaults match legacy ecalj)
  !-----------------------------------------------------------------

  ! BZ mesh and Q vector cutoffs
  integer, protected, public :: n1n2n3(3)    = [0, 0, 0]
  integer, protected, public :: n1n2n3eps(3) = [0, 0, 0]
  integer, protected, public :: BZmesh       = 1
  real(8), protected, public :: QpGcut_psi   = 4.0d0
  real(8), protected, public :: QpGcut_cou   = 3.0d0
  real(8), protected, public :: alpha_OffG   = 1.0d0
  logical, protected, public :: unit_2pioa   = .false.

  ! Sigma / chi0 mode
  integer, protected, public :: iSigMode     = 3
  integer, protected, public :: niw          = 10
  integer, protected, public :: nband_chi0   = 999
  integer, protected, public :: EMINforGW    = -9999
  integer, protected, public :: EMAXforGW    = 9999
  real(8), protected, public :: emax_sigm    = 3.0d0
  real(8), protected, public :: emax_chi0    = 999.0d0
  real(8), protected, public :: HistBin_ratio = 1.03d0
  real(8), protected, public :: HistBin_dw   = 1.0d-6
  real(8), protected, public :: deltaw       = 0.02d0
  real(8), protected, public :: esmr         = 0.003d0
  real(8), protected, public :: delta        = -1.0d-6
  real(8), protected, public :: dw           = 0.005d0
  real(8), protected, public :: omg_c        = 0.04d0
  real(8), protected, public :: WgtQ0P       = 0.01d0
  real(8), protected, public :: SmearX0      = 0.0d0
  real(8), protected, public :: GaussianFilterX0 = 0.0d0
  logical, protected, public :: GaussSmear   = .false.

  ! Optional flags
  logical, protected, public :: KeepEigen    = .false.
  logical, protected, public :: KeepPPOVL    = .false.
  logical, protected, public :: NormChk      = .false.
  logical, protected, public :: AnyQ         = .false.
  logical, protected, public :: QforEPSau    = .false.
  logical, protected, public :: QforEPSunita = .false.
  logical, protected, public :: QforEPSLIncLeft = .false.
  logical, protected, public :: tetrakbt     = .false.
  integer, protected, public :: t_tetrakbt   = 0
  ! MagAtom: variable-length integer array of magnetic-atom site indices.
  ! Allocated to size(>=1) on load; consumers use size(MagAtom) for count.
  integer, protected, public, allocatable :: MagAtom(:)
  ! nband_sigm: legacy reads as a single integer (first token if vector).
  integer, protected, public :: nband_sigm = 9999999

  ! Additional GW config (added 2026-05-02 for m_readgwinput migration)
  real(8), protected, public :: ecut_p          = 1.0d10
  real(8), protected, public :: ecuts_p         = 1.0d10
  integer, protected, public :: multitet(3)     = [1, 1, 1]
  real(8), protected, public :: gauss_img       = 1.0d0
  logical, protected, public :: KeepPositiveCou = .true.

  ! Wannier-related
  integer, protected, public :: mlo_emax     = 0
  integer, protected, public :: mlo_method   = 0
  integer, protected, public :: wan_maxit_1st = 0
  integer, protected, public :: wan_maxit_2nd = 0
  integer, protected, public :: wan_tb_cut   = 0
  real(8), protected, public :: wan_conv_1st = 1.0d-7
  real(8), protected, public :: wan_conv_end = 1.0d-8
  real(8), protected, public :: wan_max_1st  = 0.1d0
  real(8), protected, public :: wan_max_2nd  = 0.3d0
  real(8), protected, public :: wan_in_emin  = -10.0d0
  real(8), protected, public :: wan_in_emax  = 4.0d0
  real(8), protected, public :: wan_out_emin = -1.05d0
  real(8), protected, public :: wan_out_emax = 2.4d0
  logical, protected, public :: wan_in_ewin  = .false.

  !-----------------------------------------------------------------
  ! [product_basis]: structured data
  !-----------------------------------------------------------------
  real(8), protected, public, allocatable :: pb_tolerance(:)
  integer, protected, public, allocatable :: pb_lcutmx(:)
  integer, protected, public, allocatable :: pb_nlx(:,:)      ! (4, n_nlx) [iatom,l,nnvv,nnc]
  integer, protected, public, allocatable :: pb_valence(:,:)  ! (5, n_val) [iatom,l,n,occ,unocc]
  integer, protected, public, allocatable :: pb_core(:,:)     ! (7, n_core) [iatom,l,n,occ,unocc,forX0,forSxc]
  integer, protected, public :: pb_n_nlx = 0, pb_n_val = 0, pb_n_core = 0

  !-----------------------------------------------------------------
  ! [blocks]: raw text (future: parse to structured form)
  !-----------------------------------------------------------------
  character(len=:), protected, public, allocatable :: block_QPNT
  character(len=:), protected, public, allocatable :: block_QforEPS
  character(len=:), protected, public, allocatable :: block_QforEPSL
  character(len=:), protected, public, allocatable :: block_QforGW
  character(len=:), protected, public, allocatable :: block_Worb
  character(len=:), protected, public, allocatable :: block_hrotr

  !-----------------------------------------------------------------
  ! State
  !-----------------------------------------------------------------
  logical, protected, public :: gwinput_loaded = .false.

  public :: gwinput_load

contains

  subroutine gwinput_load(filename, error)
    !> Load GWinput.toml. Idempotent: returns immediately if already loaded.
    !  filename defaults to 'GWinput.toml' if absent.
    !  On parse failure, sets allocatable error message and returns
    !  without setting gwinput_loaded=.true.
    character(*),               intent(in),  optional :: filename
    character(len=:), allocatable, intent(out), optional :: error

    type(toml_table), allocatable, target :: root
    type(toml_table), pointer :: gw, pb, blocks
    type(toml_error), allocatable :: terr
    character(len=:), allocatable :: fname

    if (gwinput_loaded) return
    if (present(error)) error = ""

    if (present(filename)) then
       fname = trim(filename)
    else
       fname = 'GWinput.toml'
    endif

    call toml_load(root, fname, error=terr)
    if (allocated(terr)) then
       if (present(error)) error = "m_GWinput: parse error: " // terr%message
       return
    endif

    !---- [gw] ----
    call get_value(root, 'gw', gw)
    if (associated(gw)) call load_gw_section(gw)

    !---- [product_basis] ----
    call get_value(root, 'product_basis', pb)
    if (associated(pb)) call load_pb_section(pb)

    !---- [blocks] ----
    call get_value(root, 'blocks', blocks)
    if (associated(blocks)) call load_blocks_section(blocks)

    gwinput_loaded = .true.
  end subroutine gwinput_load


  subroutine load_gw_section(gw)
    type(toml_table), pointer, intent(in) :: gw
    ! Integer scalars
    call gv_i(gw, 'iSigMode',      iSigMode)
    call gv_i(gw, 'niw',           niw)
    call gv_i(gw, 'nband_chi0',    nband_chi0)
    call gv_i(gw, 'EMINforGW',     EMINforGW)
    call gv_i(gw, 'EMAXforGW',     EMAXforGW)
    call gv_i(gw, 'BZmesh',        BZmesh)
    call gv_i(gw, 't_tetrakbt',    t_tetrakbt)
    call gv_i(gw, 'mlo_emax',      mlo_emax)
    call gv_i(gw, 'mlo_method',    mlo_method)
    call gv_i(gw, 'wan_maxit_1st', wan_maxit_1st)
    call gv_i(gw, 'wan_maxit_2nd', wan_maxit_2nd)
    call gv_i(gw, 'wan_tb_cut',    wan_tb_cut)
    ! nband_sigm: legacy reads as integer; TOML may have list[float] -- take first as int
    call gv_iv_first(gw, 'nband_sigm', nband_sigm)
    ! MagAtom: VLA -- accept scalar or vector
    call gv_iv_alloc(gw, 'MagAtom', MagAtom)

    ! Real scalars
    call gv_r(gw, 'QpGcut_psi',    QpGcut_psi)
    call gv_r(gw, 'QpGcut_cou',    QpGcut_cou)
    call gv_r(gw, 'alpha_OffG',    alpha_OffG)
    call gv_r(gw, 'emax_sigm',     emax_sigm)
    call gv_r(gw, 'emax_chi0',     emax_chi0)
    call gv_r(gw, 'HistBin_ratio', HistBin_ratio)
    call gv_r(gw, 'HistBin_dw',    HistBin_dw)
    call gv_r(gw, 'deltaw',        deltaw)
    call gv_r(gw, 'esmr',          esmr)
    call gv_r(gw, 'delta',         delta)
    call gv_r(gw, 'dw',            dw)
    call gv_r(gw, 'omg_c',         omg_c)
    call gv_r(gw, 'WgtQ0P',        WgtQ0P)
    call gv_r(gw, 'SmearX0',       SmearX0)
    call gv_r(gw, 'GaussianFilterX0', GaussianFilterX0)
    call gv_r(gw, 'wan_conv_1st',  wan_conv_1st)
    call gv_r(gw, 'wan_conv_end',  wan_conv_end)
    call gv_r(gw, 'wan_max_1st',   wan_max_1st)
    call gv_r(gw, 'wan_max_2nd',   wan_max_2nd)
    call gv_r(gw, 'wan_in_emin',   wan_in_emin)
    call gv_r(gw, 'wan_in_emax',   wan_in_emax)
    call gv_r(gw, 'wan_out_emin',  wan_out_emin)
    call gv_r(gw, 'wan_out_emax',  wan_out_emax)

    ! Newly added (m_readgwinput migration)
    call gv_r(gw, 'ecut_p',          ecut_p)
    call gv_r(gw, 'ecuts_p',         ecuts_p)
    call gv_r(gw, 'gauss_img',       gauss_img)

    ! Boolean flags
    call gv_l(gw, 'GaussSmear',      GaussSmear)
    call gv_l(gw, 'KeepEigen',       KeepEigen)
    call gv_l(gw, 'KeepPPOVL',       KeepPPOVL)
    call gv_l(gw, 'NormChk',         NormChk)
    call gv_l(gw, 'unit_2pioa',      unit_2pioa)
    call gv_l(gw, 'AnyQ',            AnyQ)
    call gv_l(gw, 'QforEPSau',       QforEPSau)
    call gv_l(gw, 'QforEPSunita',    QforEPSunita)
    call gv_l(gw, 'QforEPSLIncLeft', QforEPSLIncLeft)
    call gv_l(gw, 'tetrakbt',        tetrakbt)
    call gv_l(gw, 'wan_in_ewin',     wan_in_ewin)
    call gv_l(gw, 'KeepPositiveCou', KeepPositiveCou)

    ! Integer vectors
    call gv_iv3(gw, 'n1n2n3',    n1n2n3)
    call gv_iv3(gw, 'n1n2n3eps', n1n2n3eps)
    call gv_iv3(gw, 'multitet',  multitet)
  end subroutine load_gw_section


  subroutine load_pb_section(pb)
    type(toml_table), pointer, intent(in) :: pb
    type(toml_array), pointer :: arr
    integer :: i, n

    ! tolerance
    call get_value(pb, 'tolerance', arr, requested=.false.)
    if (associated(arr)) then
       n = len(arr)
       allocate(pb_tolerance(n))
       do i = 1, n
          call get_value(arr, i, pb_tolerance(i))
       enddo
    endif

    ! lcutmx
    call get_value(pb, 'lcutmx', arr, requested=.false.)
    if (associated(arr)) then
       n = len(arr)
       allocate(pb_lcutmx(n))
       do i = 1, n
          call get_value(arr, i, pb_lcutmx(i))
       enddo
    endif

    ! 2-D integer arrays
    call load_int_2darray(pb, 'nlx',     4, pb_nlx,     pb_n_nlx)
    call load_int_2darray(pb, 'valence', 5, pb_valence, pb_n_val)
    call load_int_2darray(pb, 'core',    7, pb_core,    pb_n_core)
  end subroutine load_pb_section


  subroutine load_int_2darray(tbl, key, ncol, mat, nrow_out)
    type(toml_table), pointer, intent(in)  :: tbl
    character(*),              intent(in)  :: key
    integer,                   intent(in)  :: ncol
    integer, allocatable,      intent(out) :: mat(:,:)
    integer,                   intent(out) :: nrow_out
    type(toml_array), pointer :: outer, row
    integer :: i, j, n

    call get_value(tbl, key, outer, requested=.false.)
    if (.not. associated(outer)) then
       nrow_out = 0
       return
    endif
    n = len(outer)
    allocate(mat(ncol, n))
    do i = 1, n
       call get_value(outer, i, row)
       if (.not. associated(row)) cycle
       do j = 1, min(ncol, len(row))
          call get_value(row, j, mat(j, i))
       enddo
    enddo
    nrow_out = n
  end subroutine load_int_2darray


  subroutine load_blocks_section(blocks)
    type(toml_table), pointer, intent(in) :: blocks
    call gv_c_alloc(blocks, 'QPNT',     block_QPNT)
    call gv_c_alloc(blocks, 'QforEPS',  block_QforEPS)
    call gv_c_alloc(blocks, 'QforEPSL', block_QforEPSL)
    call gv_c_alloc(blocks, 'QforGW',   block_QforGW)
    call gv_c_alloc(blocks, 'Worb',     block_Worb)
    call gv_c_alloc(blocks, 'hrotr',    block_hrotr)
  end subroutine load_blocks_section


  ! ===== Helper getters (silently keep default on miss) =====

  subroutine gv_i(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    integer,                   intent(inout) :: var
    integer :: stat
    call get_value(tbl, key, var, stat=stat)
  end subroutine

  subroutine gv_r(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    real(8),                   intent(inout) :: var
    integer :: stat
    call get_value(tbl, key, var, stat=stat)
  end subroutine

  subroutine gv_l(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    logical,                   intent(inout) :: var
    integer :: stat
    call get_value(tbl, key, var, stat=stat)
  end subroutine

  subroutine gv_iv3(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    integer,                   intent(inout) :: var(3)
    type(toml_array), pointer :: arr
    integer :: i
    call get_value(tbl, key, arr, requested=.false.)
    if (.not. associated(arr)) return
    do i = 1, min(3, len(arr))
       call get_value(arr, i, var(i))
    enddo
  end subroutine

  subroutine gv_rv_alloc(tbl, key, var)
    type(toml_table), pointer, intent(in)  :: tbl
    character(*),              intent(in)  :: key
    real(8), allocatable,      intent(out) :: var(:)
    type(toml_array), pointer :: arr
    integer :: i, n
    call get_value(tbl, key, arr, requested=.false.)
    if (.not. associated(arr)) return
    n = len(arr)
    allocate(var(n))
    do i = 1, n
       call get_value(arr, i, var(i))
    enddo
  end subroutine

  !> Read an int that may be stored as scalar or as the first element of an array.
  !  Used for keys like 'nband_sigm' where TOML may have list[float] but legacy
  !  consumes a single int.
  subroutine gv_iv_first(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    integer,                   intent(inout) :: var
    type(toml_array), pointer :: arr
    integer :: stat
    real(8) :: rval
    ! Try scalar integer first
    call get_value(tbl, key, var, stat=stat)
    if (stat == 0) return
    ! Try scalar real (truncate to int)
    call get_value(tbl, key, rval, stat=stat)
    if (stat == 0) then
       var = int(rval)
       return
    endif
    ! Try array, take first element
    call get_value(tbl, key, arr, requested=.false.)
    if (.not. associated(arr)) return
    if (len(arr) < 1) return
    call get_value(arr, 1, rval, stat=stat)
    if (stat == 0) var = int(rval)
  end subroutine

  !> Allocate and fill integer array. If TOML key absent, leaves var unallocated.
  subroutine gv_iv_alloc(tbl, key, var)
    type(toml_table), pointer, intent(in)  :: tbl
    character(*),              intent(in)  :: key
    integer, allocatable,      intent(out) :: var(:)
    type(toml_array), pointer :: arr
    integer :: i, n, stat, scal
    ! Scalar form 'MagAtom 1' is common; try scalar first
    call get_value(tbl, key, scal, stat=stat)
    if (stat == 0) then
       allocate(var(1)); var(1) = scal
       return
    endif
    call get_value(tbl, key, arr, requested=.false.)
    if (.not. associated(arr)) return
    n = len(arr)
    allocate(var(n))
    do i = 1, n
       call get_value(arr, i, var(i))
    enddo
  end subroutine

  subroutine gv_c_alloc(tbl, key, var)
    type(toml_table), pointer,        intent(in)  :: tbl
    character(*),                     intent(in)  :: key
    character(len=:), allocatable,    intent(out) :: var
    integer :: stat
    call get_value(tbl, key, var, stat=stat)
  end subroutine

end module m_GWinput
