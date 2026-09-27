!> linalgtune: choose the GPU backend of each matrix operation (the table of m_linalg_policy) by measuring them here.
!>
!>   linalgtune_gpu [--out=<file>] [--ecalj=<revision>]
!>
!> Every backend runs on hgw-like shapes of each size class (small / large, split as in m_linalg_policy).  Its error is
!> checked against FP64 on the same inputs and it is timed (best of 3 batches).  For each precision level and
!> operation the fastest backend within the level's error bound is chosen when it is at least 10% faster than cuBLAS
!> over the class (geometric mean of the time ratios, each shape weighted by its rough share of hgw time) and not more
!> than 10% slower on any shape that carries weight (>= 0.1); otherwise cuBLAS (lu64 for the inverse) stays.
!> Error bounds (relative, Frobenius): cgemm at most 4 times the error of cuBLAS at the same level on the same shape
!> (and below 5e-3 at tf32, 1e-5 at fp32): the promise is FP32 (TF32) accuracy, which GEMMul8 with few moduli can
!> miss on long sums (it quantizes each row against its largest element; 2026-09-27: 7 moduli were 5 times worse
!> than cuBLAS and moved Im Sigma_c of hgw by 1e-4 eV); zgemm and dgemm 1e-12; the epstilde
!> inverse 1e-6 at tf32/fp32 and 1e-13 at fp64.
!> realhgemm (FP16 inputs, FP32 accumulation) has the error of TF32, so it can only pass at level tf32.
!> The sizes are deliberately not multiples of 64, as in real GW runs: aligned sizes favour some backends (GEMMul8)
!> and would bias the choice (2026-09-27: a table measured on 1024 chose GEMMul8 where hgw with 1053 runs slower).
!> The table goes to --out, else to ecalj_linalg_policy.toml next to this executable (ECALJ_LINALG_POLICY overrides),
!> through a temporary file and rename.  Run it on an idle GPU; InstallAll.py does so after a GPU build.  2026-09-27.
program linalgtune
  use cudafor
  use cublas_v2
  use iso_c_binding
  use m_linalg_policy, only: la_set_level, la_apply, la_policy_path, this_gpu
  use m_blas, only: cmm_d, zmm_d, dmm_d, la_cache_reset, cublas_init, cublas_handle
  use m_lapack, only: zminv_eps_d, zminv_d
  implicit none
  type shape_t
    character(10) :: name
    integer :: m, n, k
    character :: opa
    logical :: keyed
    real(8) :: weight           ! rough share of hgw time within its class
  end type
  ! The products of hgw (opB = N): Sigma_c on the imaginary / real axis W^H zmel with W kept under a key, chi0 bins
  ! zmel^H (w zmel), the change of basis in m_llw, the final zsec (small m = n, long k).
  integer, parameter :: nshape = 6, minm = 256, minn = 512, mink = 256
  real(8), parameter :: minmnk = 1d9
  ! Weights: rough shares of hgw time within each class (the Sigma_c products dominate, then the chi0 bins).
  type(shape_t), parameter :: shp(nshape) = [ &
       shape_t('sigma-imag', 1037, 8191,  1037, 'C', .true.,  0.80d0), &
       shape_t('x0-bin',     1037, 1037,  4099, 'C', .false., 0.15d0), &
       shape_t('square',     1037, 1037,  1037, 'N', .false., 0.05d0), &
       shape_t('sigma-real', 1037,  389,  1037, 'C', .true.,  0.80d0), &
       shape_t('x0-bin-s',   1037, 1037,   263, 'C', .false., 0.15d0), &
       shape_t('zsec',        131,  131, 65539, 'C', .false., 0.05d0)]
  integer, parameter :: ninv = 2, ninvn(ninv) = [389, 1037]
  character(4), parameter :: levels(3) = ['tf32', 'fp32', 'fp64']
#ifdef __GEMMUL8
  character(10), parameter :: cand_c(4) = [character(10) :: 'cublas', 'realsgemm', 'realhgemm', 'gemmul8:7']
  character(10), parameter :: cand_z(2) = [character(10) :: 'cublas', 'gemmul8:14']
#else
  character(10), parameter :: cand_c(3) = [character(10) :: 'cublas', 'realsgemm', 'realhgemm']
  character(10), parameter :: cand_z(1) = [character(10) :: 'cublas']
#endif
  character(10), parameter :: cand_e(3) = [character(10) :: 'lu64', 'mixed1', 'mixed2']
  complex(4), device, allocatable :: a4(:), b4(:), c4(:)
  complex(8), device, allocatable :: a8(:), b8(:), c8(:), cref(:), e8(:), w8(:), eref(:)
  real(8), device, allocatable :: ad(:), bd(:), cd(:), crefd(:)
  real(8) :: tc(size(cand_c), nshape, 2), ec(size(cand_c), nshape, 2)   ! cgemm at tf32 (1) and fp32 (2)
  real(8) :: tz(size(cand_z), nshape), ez(size(cand_z), nshape)
  real(8) :: td(size(cand_z), nshape), ed(size(cand_z), nshape)
  real(8) :: te(size(cand_e), ninv), ee(size(cand_e), ninv)
  character(10) :: pick_c(2, 2), pick_z(2), pick_d(2), pick_e(2, 3)  ! (class small/large, level)
  character(len=:), allocatable :: out, rev
  character(1024) :: arg
  character(16) :: sd
  integer :: istat, i, is, ic, il, maxa, maxb, maxc, nmax, ifi, ios, cudaver
  logical :: large(nshape)
  real(8) :: wsig(nshape)     ! level tf32 serves only the products of Sigma_c (--sigma_tf32): the other shapes weigh 0
  interface
    integer(c_int) function c_rename(old, new) bind(C, name='rename')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: old(*), new(*)
    end function c_rename
  end interface

  out = trim(la_policy_path()); rev = ''
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if (arg(1:6) == '--out=') out = trim(arg(7:))
    if (arg(1:8) == '--ecalj=') rev = trim(arg(9:))
  enddo
  istat = cublas_init()
  istat = cudaRuntimeGetVersion(cudaver)
  write(6,'(a)') ' linalgtune on '//trim(this_gpu())
  do is = 1, nshape
    large(is) = shp(is)%m >= minm .and. shp(is)%n >= minn .and. shp(is)%k >= mink .and. &
                real(shp(is)%m,8)*shp(is)%n*shp(is)%k >= minmnk
    wsig(is) = merge(shp(is)%weight, 0d0, index(shp(is)%name, 'sigma') == 1 .or. shp(is)%name == 'zsec')
  enddo
  maxa = maxval(shp%m*shp%k); maxb = maxval(shp%k*shp%n); maxc = maxval(shp%m*shp%n)
  allocate(a4(maxa), b4(maxb), c4(maxc), a8(maxa), b8(maxb), c8(maxc), cref(maxc))
  allocate(ad(maxa), bd(maxb), cd(maxc), crefd(maxc))
  write(6,'(a)') ' op     level shape         m     n     k opA key  backend        ms/call   error'

  ! zgemm and dgemm (FP64 promise, same choice at every level), inputs in full double precision
  call fill(a8, maxa, 1d0, 3d0); call fill(b8, maxb, 2d0, 2d0)
  call copy_re(a8, ad, maxa); call copy_re(b8, bd, maxb)
  do is = 1, nshape
    call ref_z(shp(is)); call ref_d(shp(is))
    do ic = 1, size(cand_z)
      call measure('zgemm', 'fp64', cand_z(ic), shp(is), tz(ic,is), ez(ic,is))
      call measure('dgemm', 'fp64', cand_z(ic), shp(is), td(ic,is), ed(ic,is))
    enddo
  enddo
  do i = 1, 2
    pick_z(i) = choose(cand_z, tz, ez, 1d-12, merge(large, .not.large, i == 2), shp%weight)
    pick_d(i) = choose(cand_z, td, ed, 1d-12, merge(large, .not.large, i == 2), shp%weight)
  enddo

  ! cgemm at tf32 and fp32: inputs rounded to single precision, reference FP64 on the rounded inputs
  call round4(a8, a4, maxa); call round4(b8, b4, maxb)
  do is = 1, nshape
    call ref_z(shp(is))
    do il = 1, 2
      do ic = 1, size(cand_c)
        call measure('cgemm', levels(il), cand_c(ic), shp(is), tc(ic,is,il), ec(ic,is,il))
      enddo
    enddo
  enddo
  do il = 1, 2
    do i = 1, 2
      pick_c(i,il) = choose(cand_c, tc(:,:,il), ec(:,:,il), merge(5d-3, 1d-5, il == 1), merge(large, .not.large, i == 2), &
                          merge(wsig, shp%weight, il == 1), relcublas=.true.)
    enddo
  enddo
  deallocate(a4, b4, c4, a8, b8, c8, cref, ad, bd, cd, crefd)

  ! Inverse of epstilde; mixed1/mixed2 use zgemm, so the zgemm rows chosen above are in place.  Timed once (at fp64);
  ! the choice differs by level only through the error bound.
  call la_set_level('fp64')
  call la_apply('fp64.zgemm.small='//trim(pick_z(1))//',fp64.zgemm.large='//trim(pick_z(2)))
  nmax = maxval(ninvn)
  allocate(e8(nmax*nmax), w8(nmax*nmax), eref(nmax*nmax))
  do i = 1, ninv
    call make_eps(ninvn(i))
    do ic = 1, size(cand_e)
      call measure_inv(cand_e(ic), ninvn(i), te(ic,i), ee(ic,i))
    enddo
  enddo
  do il = 1, 3
    pick_e(1,il) = choose_inv(1, merge(1d-13, 1d-6, il == 3))
    pick_e(2,il) = choose_inv(2, merge(1d-13, 1d-6, il == 3))
  enddo
  deallocate(e8, w8, eref)

  ! the policy file
  open(newunit=ifi, file=out//'.tmp', status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    write(6,'(a)') ' linalgtune: cannot write '//out//'.tmp; nothing written.'
    stop 1
  endif
  call date_and_time(date=sd)
  write(ifi,'(a)') '# ecalj linalg policy: which GPU backend runs each matrix operation (read by m_linalg_policy).'
  write(ifi,'(a)') '# Written by linalgtune on '//sd(1:8)//' from measurements on this GPU.  Rerun linalgtune after'
  write(ifi,'(a)') '# changing the GPU, the driver or CUDA.  Rows: <level>.<op>.<small|large|minm|minn|mink|minmnk>.'
  write(ifi,'(a)') 'gpu = "'//trim(this_gpu())//'"'
  write(ifi,'(a,i0)') 'cuda = ', cudaver
  if (len(rev) > 0) write(ifi,'(a)') 'ecalj = "'//rev//'"'
  do il = 1, 3
    write(ifi,'(a)') ''
    call row(ifi, levels(il), 'cgemm', pick_c(:, min(il,2)))   ! fp64 builds do no complex(4) products; fp32 choice
    call row(ifi, levels(il), 'zgemm', pick_z)
    call row(ifi, levels(il), 'dgemm', pick_d)
    call row(ifi, levels(il), 'epsinv', pick_e(:,il))
  enddo
  write(ifi,'(a)') ''
  write(ifi,'(a)') '# Measurements (ms per call, best of 3 batches; error relative to FP64, Frobenius):'
  do is = 1, nshape
    do ic = 1, size(cand_z)
      call line(ifi, 'zgemm', 'fp64', shp(is), cand_z(ic), tz(ic,is), ez(ic,is))
      call line(ifi, 'dgemm', 'fp64', shp(is), cand_z(ic), td(ic,is), ed(ic,is))
    enddo
    do il = 1, 2
      do ic = 1, size(cand_c)
        call line(ifi, 'cgemm', levels(il), shp(is), cand_c(ic), tc(ic,is,il), ec(ic,is,il))
      enddo
    enddo
  enddo
  do i = 1, ninv
    do ic = 1, size(cand_e)
      call line(ifi, 'epsinv', 'fp64', shape_t('eps', ninvn(i), ninvn(i), ninvn(i), 'N', .false., 1d0), &
                cand_e(ic), te(ic,i), ee(ic,i))
    enddo
  enddo
  close(ifi)
  if (c_rename(out//'.tmp'//c_null_char, out//c_null_char) /= 0) then
    write(6,'(a)') ' linalgtune: rename to '//out//' failed; the table is in '//out//'.tmp'
    stop 1
  endif
  write(6,'(a)') ' linalgtune: wrote '//out
  do il = 1, 3
    write(6,'(3x,a,": cgemm ",a,"/",a,"  zgemm ",a,"/",a,"  dgemm ",a,"/",a,"  epsinv ",a,"/",a," (small/large)")') &
         levels(il), (trim(pick_c(i,min(il,2))), i=1,2), (trim(pick_z(i)), i=1,2), (trim(pick_d(i)), i=1,2), &
         (trim(pick_e(i,il)), i=1,2)
  enddo

contains

  subroutine measure(op, lv, cand, s, tms, err)
    !> Error and time of one backend (cand) for op at level lv on shape s, called through m_blas as hgw calls it.
    character(*), intent(in) :: op, lv, cand
    type(shape_t), intent(in) :: s
    real(8), intent(out) :: tms, err
    integer :: key
    key = merge(7, -1, s%keyed)
    call la_set_level(lv)
    call la_apply(lv//'.'//op//'.small='//trim(cand)//','//lv//'.'//op//'.large='//trim(cand))
    call la_cache_reset()
    call run(op, s, key)                                  ! warm-up; with a key it also keeps A
    istat = cudaDeviceSynchronize()
    select case (op)
    case ('cgemm'); err = relerr_c(c4, cref, s%m*s%n)
    case ('zgemm'); err = relerr_z(c8, cref, s%m*s%n)
    case ('dgemm'); err = relerr_d(cd, crefd, s%m*s%n)
    end select
    tms = timed(op, s, key)
    call la_cache_reset()
    write(6,'(1x,a6,1x,a4,1x,a10,3i6,2x,a,2x,l1,4x,a12,f10.4,es10.2)') op, lv, s%name, s%m, s%n, s%k, s%opa, &
         s%keyed, cand, tms, err
  end subroutine measure

  real(8) function timed(op, s, key) result(tms)
    !> ms per call: best of 3 batches, each about 50 ms long.
    character(*), intent(in) :: op
    type(shape_t), intent(in) :: s
    integer, intent(in) :: key
    real(8) :: t0, t1
    integer :: nrep, ib, ir
    t0 = now()
    call run(op, s, key)
    istat = cudaDeviceSynchronize()
    t1 = now() - t0
    nrep = max(2, min(200, nint(0.05d0/max(t1, 1d-6))))
    tms = huge(1d0)
    do ib = 1, 3
      t0 = now()
      do ir = 1, nrep
        call run(op, s, key)
      enddo
      istat = cudaDeviceSynchronize()
      tms = min(tms, (now() - t0)/nrep*1d3)
    enddo
  end function timed

  subroutine run(op, s, key)
    character(*), intent(in) :: op
    type(shape_t), intent(in) :: s
    integer, intent(in) :: key
    select case (op)
    case ('cgemm'); istat = cmm_d(a4, b4, c4, s%m, s%n, s%k, opa=s%opa, key=key)
    case ('zgemm'); istat = zmm_d(a8, b8, c8, s%m, s%n, s%k, opa=s%opa, key=key)
    case ('dgemm'); istat = dmm_d(ad, bd, cd, s%m, s%n, s%k, opa=merge('T', s%opa, s%opa == 'C'), key=key)
    end select
  end subroutine run

  subroutine ref_z(s)
    !> cref = op(A) B in FP64 (cuBLAS zgemm, not the 3M variant) from a8, b8.
    type(shape_t), intent(in) :: s
    istat = cublasZgemm(cublas_handle, merge(CUBLAS_OP_C, CUBLAS_OP_N, s%opa == 'C'), CUBLAS_OP_N, s%m, s%n, s%k, &
                        (1d0,0d0), a8, merge(s%k, s%m, s%opa == 'C'), b8, s%k, (0d0,0d0), cref, s%m)
  end subroutine ref_z

  subroutine ref_d(s)
    type(shape_t), intent(in) :: s
    istat = cublasDgemm(cublas_handle, merge(CUBLAS_OP_T, CUBLAS_OP_N, s%opa == 'C'), CUBLAS_OP_N, s%m, s%n, s%k, &
                        1d0, ad, merge(s%k, s%m, s%opa == 'C'), bd, s%k, 0d0, crefd, s%m)
  end subroutine ref_d

  subroutine measure_inv(cand, n, tms, err)
    !> The epstilde inverse through zminv_eps_d, with the row epsinv set to cand; w8 = E^-1, compared with eref.
    character(*), intent(in) :: cand
    integer, intent(in) :: n
    real(8), intent(out) :: tms, err
    real(8) :: t0, t1
    integer :: nrep, ib, ir
    call la_set_level('fp64')
    call la_apply('fp64.epsinv.small='//trim(cand)//',fp64.epsinv.large='//trim(cand))
    call inv_once(n)
    istat = cudaDeviceSynchronize()
    err = relerr_z(w8, eref, n*n)
    t0 = now()
    call inv_once(n)
    istat = cudaDeviceSynchronize()
    t1 = now() - t0
    nrep = max(2, min(50, nint(0.05d0/max(t1, 1d-6))))
    tms = huge(1d0)
    do ib = 1, 3
      t0 = now()
      do ir = 1, nrep
        call inv_once(n)
      enddo
      istat = cudaDeviceSynchronize()
      tms = min(tms, (now() - t0)/nrep*1d3)
    enddo
    write(6,'(1x,a6,1x,a4,1x,a10,3i6,2x,a,2x,l1,4x,a12,f10.4,es10.2)') 'epsinv', 'fp64', 'eps', n, n, n, 'N', &
         .false., cand, tms, err
  end subroutine measure_inv

  subroutine inv_once(n)
    integer, intent(in) :: n
    integer :: i
    !$cuf kernel do <<<*,*>>>
    do i = 1, n*n
      w8(i) = e8(i)
    enddo
    istat = zminv_eps_d(w8, n, n)
  end subroutine inv_once

  subroutine make_eps(n)
    !> A test epstilde: diagonal 1..10, off-diagonal of norm about 0.3, so that cond is about 10-30; eref = its
    !> inverse by FP64 LU.
    integer, intent(in) :: n
    integer :: i, j
    real(8) :: u1, u2, x
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, n
      do i = 1, n
        x = real(i + (j-1)*n, 8)
        u1 = sin(x*12.9898d0 + 5d0*78.233d0)*43758.5453d0; u1 = u1 - floor(u1)
        u2 = sin(x*12.9898d0 + 6d0*78.233d0)*43758.5453d0; u2 = u2 - floor(u2)
        e8(i + (j-1)*n) = cmplx(u1 - 0.5d0, u2 - 0.5d0, kind=8)*0.6d0/sqrt(real(n,8))
        if (i == j) e8(i + (j-1)*n) = e8(i + (j-1)*n) + 1d0 + 9d0*real(i-1,8)/real(n,8)
      enddo
    enddo
    !$cuf kernel do <<<*,*>>>
    do i = 1, n*n
      eref(i) = e8(i)
    enddo
    istat = zminv_d(eref, n, n)
  end subroutine make_eps

  character(10) function choose(cand, t, e, bound, inclass, w, relcublas) result(pick)
    !> cand(1) (cuBLAS) unless another backend is within the error bound on every shape of the class, at least 10%
    !> faster over the class (geometric mean of the time ratios, weights w) and at most 10% slower on each shape of
    !> weight >= 0.1.  Earlier candidates win ties within 5%.  relcublas: the bound on each shape is also 4 times the
    !> error of cand(1) there.
    character(10), intent(in) :: cand(:)
    real(8), intent(in) :: t(:,:), e(:,:), bound, w(:)
    logical, intent(in) :: inclass(:)
    logical, intent(in), optional :: relcublas
    real(8) :: g, best, r, ws, bnd(size(inclass))
    integer :: ic, is
    bnd = bound
    if (present(relcublas)) then
      if (relcublas) bnd = min(bound, 4d0*e(1,:))
    endif
    pick = cand(1); best = 0.9d0
    do ic = 2, size(cand)
      if (any(inclass .and. e(ic,:) > bnd)) cycle
      g = 0d0; ws = 0d0; r = 0d0
      do is = 1, size(inclass)
        if (.not. inclass(is)) cycle
        g = g + w(is)*log(t(ic,is)/t(1,is)); ws = ws + w(is)
        if (w(is) >= 0.1d0) r = max(r, t(ic,is)/t(1,is))
      enddo
      if (ws == 0d0 .or. r > 1.1d0) cycle
      g = exp(g/ws)
      if (g < best*0.95d0 .or. (pick == cand(1) .and. g < best)) then
        pick = cand(ic); best = g
      endif
    enddo
  end function choose

  character(10) function choose_inv(i, bound) result(pick)
    !> The same rule for the inverse; one size per class (1: small, 2: large).
    integer, intent(in) :: i
    real(8), intent(in) :: bound
    logical :: one(1)
    real(8) :: w1(1)
    one = .true.
    w1 = 1d0
    pick = choose(cand_e, te(:,i:i), ee(:,i:i), bound, one, w1)
  end function choose_inv

  subroutine row(ifi, lv, op, pick)
    integer, intent(in) :: ifi
    character(*), intent(in) :: lv, op
    character(10), intent(in) :: pick(2)
    write(ifi,'(a)') lv//'.'//op//'.small = "'//trim(pick(1))//'"'
    write(ifi,'(a)') lv//'.'//op//'.large = "'//trim(pick(2))//'"'
    write(ifi,'(a,i0,a)') lv//'.'//op//'.minm = ', minm, '   # large: m >= minm, n >= minn, k >= mink, m n k >= minmnk'
    write(ifi,'(a,i0)') lv//'.'//op//'.minn = ', minn
    write(ifi,'(a,i0)') lv//'.'//op//'.mink = ', mink
    write(ifi,'(a,es8.1)') lv//'.'//op//'.minmnk = ', minmnk
  end subroutine row

  subroutine line(ifi, op, lv, s, cand, tms, err)
    integer, intent(in) :: ifi
    character(*), intent(in) :: op, lv, cand
    type(shape_t), intent(in) :: s
    real(8), intent(in) :: tms, err
    write(ifi,'(a,a6,1x,a4,1x,a10,3i6,1x,a,1x,a5,1x,a10,f10.4," ms",es10.2)') '# ', op, lv, s%name, s%m, s%n, s%k, &
         s%opa, merge('key  ', 'nokey', s%keyed), cand, tms, err
  end subroutine line

  real(8) function now()
    integer(8) :: c, r
    call system_clock(c, r)
    now = real(c,8)/real(r,8)
  end function now

  subroutine fill(a, n, seed, decades)
    !> a(i) = (u1 - 1/2, u2 - 1/2) 10^(-decades u3): entries over several decades, as W and zmel have.
    complex(8), device :: a(*)
    integer, intent(in) :: n
    real(8), intent(in) :: seed, decades
    integer :: i
    real(8) :: x, h1, h2, h3
    !$cuf kernel do <<<*,*>>>
    do i = 1, n                                     ! h = frac(43758.5453 sin(12.9898 x + 78.233 seed)), a hash in [0,1)
      x = real(i,8)
      h1 = sin(x*12.9898d0 + seed*78.233d0)*43758.5453d0;          h1 = h1 - floor(h1)
      h2 = sin(x*12.9898d0 + (seed+0.5d0)*78.233d0)*43758.5453d0;  h2 = h2 - floor(h2)
      h3 = sin(x*12.9898d0 + (seed+0.25d0)*78.233d0)*43758.5453d0; h3 = h3 - floor(h3)
      a(i) = cmplx(h1 - 0.5d0, h2 - 0.5d0, kind=8)*10d0**(-decades*h3)
    enddo
  end subroutine fill

  subroutine copy_re(a, b, n)
    complex(8), device :: a(*)
    real(8), device :: b(*)
    integer, intent(in) :: n
    integer :: i
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      b(i) = real(a(i))
    enddo
  end subroutine copy_re

  subroutine round4(a, b, n)
    !> b = a rounded to single precision, and a = b again (the reference then sees the rounded inputs).
    complex(8), device :: a(*)
    complex(4), device :: b(*)
    integer, intent(in) :: n
    integer :: i
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      b(i) = cmplx(a(i), kind=4)
      a(i) = cmplx(b(i), kind=8)
    enddo
  end subroutine round4

  real(8) function relerr_c(c, r, n)
    complex(4), device :: c(*)
    complex(8), device :: r(*)
    integer, intent(in) :: n
    real(8) :: num, den
    integer :: i
    num = 0d0; den = 0d0
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      num = num + abs(cmplx(c(i), kind=8) - r(i))**2
    enddo
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      den = den + abs(r(i))**2
    enddo
    relerr_c = sqrt(num/den)
  end function relerr_c

  real(8) function relerr_z(c, r, n)
    complex(8), device :: c(*), r(*)
    integer, intent(in) :: n
    real(8) :: num, den
    integer :: i
    num = 0d0; den = 0d0
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      num = num + abs(c(i) - r(i))**2
    enddo
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      den = den + abs(r(i))**2
    enddo
    relerr_z = sqrt(num/den)
  end function relerr_z

  real(8) function relerr_d(c, r, n)
    real(8), device :: c(*), r(*)
    integer, intent(in) :: n
    real(8) :: num, den
    integer :: i
    num = 0d0; den = 0d0
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      num = num + (c(i) - r(i))**2
    enddo
    !$cuf kernel do <<<*,*>>>
    do i = 1, n
      den = den + r(i)**2
    enddo
    relerr_d = sqrt(num/den)
  end function relerr_d
end program linalgtune
