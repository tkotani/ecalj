!> Which implementation (backend) runs each kind of GPU matrix product and the epstilde inverse.
!>
!> One table per precision level (tf32 / fp32 / fp64; the level is fixed by the build (MP or not) and --use_fp32).
!> A row is  <level>.<op>.<small|large> = <backend>[:<moduli>]  plus the size thresholds that separate small from large.
!>   op      : cgemm (complex single), zgemm (complex double), dgemm (real double), epsinv (inverse of epstilde)
!>   backend : cublas | realsgemm | gemmul8[:moduli]      (epsinv: lu64 | mixed1 | mixed2)
!> The table is filled in this order, later ones win:
!>   1. built-in defaults: cuBLAS everywhere, lu64 (what the code did before the table existed)
!>   2. the policy file written at installation by linalgtune, <bindir>/ecalj_linalg_policy.toml,
!>      used only when its gpu line matches this GPU (or the file sets gpu = "any")
!>   3. --use_gemmul8 (old switch: GEMMul8 for the large products)
!>   4. --linalg=<row>,<row>,...   e.g. --linalg=fp32.cgemm.large=realsgemm,fp64.zgemm.large=gemmul8:14
!> The table in use is printed once.  Backends that cannot do a given call (e.g. realsgemm with opB /= N) fall back
!> to cuBLAS inside m_blas.  The level also fixes the arithmetic of cuBLAS single precision (tf32: TF32, else FP32).
!> linalgtune (SRC/main/linalgtune.f90) sets the level and the rows itself (la_set_level, la_apply).  2026-09-27.
!> --sigma_tf32 (gwsc --prec=tf32): the products of Sigma_c (callers pass sigma=.true.) take the rows of level tf32 and
!> TF32 arithmetic, everything else stays at the level of the run (fp32).  In Sigma_c the TF32 error enters linearly;
!> in chi0 -> W it is amplified by (1 - v chi0)^-1 (LiTi2O4 6^3: Re Sigma_c within 0.8 meV of FP64 near E_F, against
!> 5 meV with every product in TF32).
module m_linalg_policy
  implicit none
  private
  public :: la_init, la_backend, la_moduli, la_level, la_print, this_gpu, la_set_level, la_apply, la_policy_path, la_sigma_tf32
  integer, parameter, public :: BK_CUBLAS = 0, BK_REALSGEMM = 1, BK_GEMMUL8 = 2
  integer, parameter, public :: BK_LU64 = 10, BK_MIXED1 = 11, BK_MIXED2 = 12
  integer, parameter, public :: OP_CGEMM = 1, OP_ZGEMM = 2, OP_DGEMM = 3, OP_EPSINV = 4
  integer, parameter :: nop = 4
  character(6), parameter :: opname(nop) = ['cgemm ', 'zgemm ', 'dgemm ', 'epsinv']
  type rule_t
    integer :: small = BK_CUBLAS, large = BK_CUBLAS
    integer :: msmall = 0, mlarge = 0      ! GEMMul8 moduli (0 = the module default)
    integer :: minm = 256, minn = 512, mink = 256
    real(8) :: minmnk = 1d9                ! large: m>=minm, n>=minn, k>=mink and m*n*k>=minmnk
  end type
  type(rule_t), save :: rule(nop), rule_s(nop)   ! rule_s: rows of level tf32 for the Sigma_c products (--sigma_tf32)
  character(4), save :: level = 'fp64'
  logical, save :: sigma_tf32 = .false.
  character(256), save :: source = 'built-in defaults'
  logical, save :: inited = .false.
contains
  subroutine la_init()
    use m_cmdopt_registry, only: c0_use_fp32, c0_use_gemmul8, c2_linalg, c0_sigma_tf32
    if (inited) return
    inited = .true.
#ifdef __MP
    level = merge('fp32', 'tf32', c0_use_fp32)
    sigma_tf32 = c0_sigma_tf32 .and. level /= 'tf32'
#else
    level = 'fp64'
#endif
    rule(OP_EPSINV)%small = BK_LU64
    rule(OP_EPSINV)%large = BK_LU64
    call read_policy_file()
    if (c0_use_gemmul8) then                        ! the old switch: GEMMul8 for the large products
      rule(OP_CGEMM)%large = BK_GEMMUL8
      rule(OP_CGEMM)%minm = 1000; rule(OP_CGEMM)%minn = 1000; rule(OP_CGEMM)%mink = 1000; rule(OP_CGEMM)%minmnk = 1d10
      rule(OP_ZGEMM)%large = BK_GEMMUL8; rule(OP_DGEMM)%large = BK_GEMMUL8
      rule(OP_ZGEMM)%minm = 64; rule(OP_ZGEMM)%minn = 64; rule(OP_ZGEMM)%mink = 64; rule(OP_ZGEMM)%minmnk = 1d8
      rule(OP_DGEMM)%minm = 64; rule(OP_DGEMM)%minn = 64; rule(OP_DGEMM)%mink = 64; rule(OP_DGEMM)%minmnk = 1d8
      source = trim(source)//' + --use_gemmul8'
    endif
    if (len_trim(c2_linalg) > 0) then
      call apply_rows(c2_linalg, ',')
      source = trim(source)//' + --linalg'
    endif
    call print_on_rank0()
  end subroutine la_init

  subroutine print_on_rank0()
    !> The table in use, once, from MPI rank 0 (m_mpi is not used here: m_mpi -> m_gpu -> m_blas -> this module).
    use mpi
    use m_lgunit, only: stdo
    logical :: ini
    integer :: rank, ierr
    rank = 0
    call MPI_Initialized(ini, ierr)
    if (ini) call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) call la_print(stdo)
  end subroutine print_on_rank0

  integer function la_backend(op, m, n, k, sigma) result(bk)
    !> Backend for one call: the small or the large row of this level's table (of level tf32 for a Sigma_c product,
    !> sigma=.true., under --sigma_tf32).
    integer, intent(in) :: op, m, n, k
    logical, intent(in), optional :: sigma
    type(rule_t) :: r
    call la_init()
    r = pick(op, sigma)
    if (islarge(r, m, n, k)) then
      bk = r%large
    else
      bk = r%small
    endif
  end function la_backend

  integer function la_moduli(op, m, n, k, sigma) result(nm)
    !> GEMMul8 moduli of the row la_backend used (0 = the default of m_gemmul8).
    integer, intent(in) :: op, m, n, k
    logical, intent(in), optional :: sigma
    type(rule_t) :: r
    call la_init()
    r = pick(op, sigma)
    nm = merge(r%mlarge, r%msmall, islarge(r, m, n, k))
  end function la_moduli

  logical function la_sigma_tf32()
    !> Sigma_c products run with TF32 arithmetic (--sigma_tf32 at level fp32).
    call la_init()
    la_sigma_tf32 = sigma_tf32
  end function la_sigma_tf32

  type(rule_t) function pick(op, sigma) result(r)
    integer, intent(in) :: op
    logical, intent(in), optional :: sigma
    r = rule(op)
    if (present(sigma)) then
      if (sigma .and. sigma_tf32) r = rule_s(op)
    endif
  end function pick

  character(4) function la_level()
    call la_init()
    la_level = level
  end function la_level

  subroutine la_set_level(lv)
    !> For linalgtune: work at level lv (tf32, fp32 or fp64) whatever the build.
    character(*), intent(in) :: lv
    call la_init()
    level = lv
  end subroutine la_set_level

  subroutine la_apply(rows)
    !> For linalgtune: apply rows '<level>.<op>.<field>=<value>,...' (as --linalg) to the table in use.
    character(*), intent(in) :: rows
    call la_init()
    call apply_rows(rows, ',')
  end subroutine la_apply

  logical function islarge(r, m, n, k)
    type(rule_t), intent(in) :: r
    integer, intent(in) :: m, n, k
    islarge = m >= r%minm .and. n >= r%minn .and. k >= r%mink .and. real(m,8)*real(n,8)*real(k,8) >= r%minmnk
  end function islarge

  subroutine la_print(iout)
    !> The table in use, one line per op.
    integer, intent(in) :: iout
    integer :: op
    write(iout,'(a)') ' linalg policy ('//level//', from '//trim(source)//'):'
    do op = 1, nop
      write(iout,'(3x,a,a,a,a,a,i5,i6,i5,es9.1)') opname(op), ' small=', trim(bkname(rule(op)%small, rule(op)%msmall)), &
           ' large=', trim(bkname(rule(op)%large, rule(op)%mlarge)), rule(op)%minm, rule(op)%minn, rule(op)%mink, rule(op)%minmnk
    enddo
    if (sigma_tf32) write(iout,'(3x,a,a,a,a,a)') 'Sigma_c products in TF32 (--sigma_tf32, rows of tf32): cgemm small=', &
         trim(bkname(rule_s(OP_CGEMM)%small, rule_s(OP_CGEMM)%msmall)), ' large=', &
         trim(bkname(rule_s(OP_CGEMM)%large, rule_s(OP_CGEMM)%mlarge)), ''
  end subroutine la_print

  function bkname(bk, nm) result(s)
    integer, intent(in) :: bk, nm
    character(20) :: s
    select case (bk)
    case (BK_CUBLAS);    s = 'cublas'
    case (BK_REALSGEMM); s = 'realsgemm'
    case (BK_GEMMUL8);   s = 'gemmul8'
      if (nm > 0) write(s,'(a,i0)') 'gemmul8:', nm
    case (BK_LU64);      s = 'lu64'
    case (BK_MIXED1);    s = 'mixed1'
    case (BK_MIXED2);    s = 'mixed2'
    case default;        s = '?'
    end select
  end function bkname

  subroutine read_policy_file()
    !> <bindir>/ecalj_linalg_policy.toml: dotted keys, one per line, e.g.
    !>   gpu = "NVIDIA GeForce RTX 5090"
    !>   fp32.cgemm.large = "realsgemm"
    !>   fp32.cgemm.minn = 512
    !> Rows of other levels are read but only this level's are used.  The file is written once at installation;
    !> it is only read here.  ECALJ_LINALG_POLICY overrides the path.
    character(1024) :: path, line
    character(256) :: gpu
    integer :: ifi, ios, lb
    path = la_policy_path()
    open(newunit=ifi, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    gpu = ''
    do
      read(ifi,'(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      if (line(1:3) == 'gpu') then
        lb = index(line, '=')
        gpu = unquote(line(lb+1:))
        if (trim(gpu) /= 'any' .and. trim(gpu) /= trim(this_gpu())) then
          close(ifi)
          source = 'built-in defaults (policy file is for "'//trim(gpu)//'")'
          return
        endif
      endif
    enddo
    rewind(ifi)
    do
      read(ifi,'(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line(1:1) == '#' .or. len_trim(line) == 0 .or. line(1:3) == 'gpu') cycle
      call apply_row(line)
    enddo
    close(ifi)
    source = trim(path)
  end subroutine read_policy_file

  function this_gpu() result(nm)
    !> Name of the current CUDA device ('none' without GPU), as nvidia-smi prints it.
#ifdef __GPU
    use cudafor
#endif
    character(256) :: nm
#ifdef __GPU
    type(cudaDeviceProp) :: prop
    integer :: dev, istat
    istat = cudaGetDevice(dev)
    istat = cudaGetDeviceProperties(prop, dev)
    nm = trim(prop%name)
#else
    nm = 'none'
#endif
  end function this_gpu

  function la_policy_path() result(path)
    !> ECALJ_LINALG_POLICY, else ecalj_linalg_policy.toml next to the running executable (where InstallAll puts it).
    character(1024) :: path
    integer :: st
    call get_environment_variable('ECALJ_LINALG_POLICY', path, status=st)
    if (st == 0 .and. len_trim(path) > 0) return
    call exe_dir(path)
    path = trim(path)//'/ecalj_linalg_policy.toml'
  end function la_policy_path

  subroutine exe_dir(dir)
    !> Directory of the running executable as it was called (the bin directory, not the target of its symlink);
    !> a bare name is looked up in PATH.
    character(*), intent(out) :: dir
    character(1024) :: exe
    character(4096) :: pth
    integer :: i, i0, i1, st
    logical :: ex
    call get_command_argument(0, exe)
    i = index(exe, '/', back=.true.)
    dir = '.'
    if (i > 1) then
      dir = exe(1:i-1)
      return
    endif
    call get_environment_variable('PATH', pth, status=st)
    if (st /= 0) return
    i0 = 1
    do while (i0 <= len_trim(pth))
      i1 = index(pth(i0:), ':')
      if (i1 == 0) then
        i1 = len_trim(pth) + 1
      else
        i1 = i0 + i1 - 1
      endif
      if (i1 > i0) then
        inquire(file=pth(i0:i1-1)//'/'//trim(exe), exist=ex)
        if (ex) then
          dir = pth(i0:i1-1)
          return
        endif
      endif
      i0 = i1 + 1
    enddo
  end subroutine exe_dir

  subroutine apply_rows(rows, sep)
    character(*), intent(in) :: rows
    character(1), intent(in) :: sep
    integer :: i0, i1
    i0 = 1
    do while (i0 <= len_trim(rows))
      i1 = index(rows(i0:), sep)
      if (i1 == 0) then
        call apply_row(rows(i0:len_trim(rows)))
        exit
      endif
      call apply_row(rows(i0:i0+i1-2))
      i0 = i0 + i1
    enddo
  end subroutine apply_rows

  subroutine apply_row(row)
    !> <level>.<op>.<field> = <value>; rows of other levels are ignored.
    character(*), intent(in) :: row
    character(64) :: key, lv, opn, field
    character(64) :: val
    integer :: ie, i1, i2, op, bk, nm, ios, i
    ie = index(row, '=')
    if (ie == 0) return
    key = adjustl(row(1:ie-1))
    val = unquote(row(ie+1:))
    i1 = index(key, '.')
    if (i1 == 0) return
    i2 = index(key(i1+1:), '.')
    if (i2 == 0) return
    lv = key(1:i1-1)
    opn = key(i1+1:i1+i2-1)
    field = trim(key(i1+i2+1:))
    if (trim(lv) /= level .and. .not. (sigma_tf32 .and. trim(lv) == 'tf32')) return
    op = 0
    do i = 1, nop                                ! not findloc: nvfortran compares strings of unequal length unpadded
      if (trim(opname(i)) == trim(opn)) op = i
    enddo
    if (op == 0) call rx('linalg policy: unknown op in "'//trim(row)//'"')
    if (trim(lv) == level) call set_field(rule(op))
    if (sigma_tf32 .and. trim(lv) == 'tf32') call set_field(rule_s(op))
  contains
    subroutine set_field(r)
      type(rule_t), intent(inout) :: r
      select case (trim(field))
      case ('small', 'large')
        call parse_backend(val, bk, nm)
        if (trim(field) == 'small') then
          r%small = bk; r%msmall = nm
        else
          r%large = bk; r%mlarge = nm
        endif
      case ('minm');   read(val,*,iostat=ios) r%minm
      case ('minn');   read(val,*,iostat=ios) r%minn
      case ('mink');   read(val,*,iostat=ios) r%mink
      case ('minmnk'); read(val,*,iostat=ios) r%minmnk
      case default
        call rx('linalg policy: unknown field in "'//trim(row)//'"')
      end select
    end subroutine set_field
  end subroutine apply_row

  subroutine parse_backend(val, bk, nm)
    character(*), intent(in) :: val
    integer, intent(out) :: bk, nm
    integer :: ic, ios
    nm = 0
    ic = index(val, ':')
    if (ic > 0) then
      read(val(ic+1:),*,iostat=ios) nm
      if (ios /= 0) nm = 0
    endif
    select case (val(1:merge(ic-1, len_trim(val), ic > 0)))
    case ('cublas');    bk = BK_CUBLAS
    case ('realsgemm'); bk = BK_REALSGEMM
    case ('gemmul8');   bk = BK_GEMMUL8
    case ('lu64');      bk = BK_LU64
    case ('mixed1');    bk = BK_MIXED1
    case ('mixed2');    bk = BK_MIXED2
    case default
      call rx('linalg policy: unknown backend "'//trim(val)//'"')
    end select
  end subroutine parse_backend

  function unquote(s) result(u)
    character(*), intent(in) :: s
    character(256) :: u
    integer :: ic
    u = adjustl(s)
    ic = index(u, '#')                          ! trailing comment
    if (ic > 0) u = u(1:ic-1)
    u = trim(adjustl(u))
    if (len_trim(u) >= 2 .and. u(1:1) == '"') u = u(2:len_trim(u)-1)
  end function unquote
end module m_linalg_policy
