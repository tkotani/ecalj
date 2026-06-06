!> Generate the GW-related part of ctrlg.<sname>.toml + PB.<sname>.toml.
!  Input: ctrlg.<sname>.toml (must already exist with ctrl sections in
!         place -- typically created by ctrlgenToml.py or Legacy2toml.py).
!  Output:
!    ctrlg.<sname>.toml   (existing ctrl sections preserved; any old
!                          [gw]/[product_basis]/[blocks] stripped and
!                          replaced with freshly-generated ones)
!    PB.<sname>.toml              (per-atom product basis tables: nlx/valence/core)
!
!  Atomic write: each output is first written to <name>.partial, then
!  renamed in place. Caller (ctrlgenToml.py / mkGWinput) is responsible
!  for taking a .bakup snapshot of any pre-existing files before invoking.
module m_gwinit
  private
  public :: gwinit_v2
contains

  subroutine gwinit_v2() bind(C)
    use m_hamindex0,only: readhamindex0,nbas,lmxax,konft, &
         lmxaa => lmxa, nindx, lindx, mnla => ndima, iat => ibasindx, &
         caption, spid, npqn, zz, pqn
    use m_ext, only: sname
    implicit none
    integer :: ifi, ifpb
    integer :: ibas, lk, izz, nx, isp
    integer :: iatbk, kkk
    integer :: nocc, nunocc, noccc, nunoccc, ncinc, ncinc2
    integer, parameter :: n1q = 4, n2q = 4, n3q = 4   ! default BZ mesh
    integer, allocatable :: konf(:,:), nncx(:,:), nnvv(:,:), lcutmx(:)
    character(len=1) :: seg2
    character(len=:), allocatable :: ctrlg, ctrlg_partial, pb_partial

    write(6,"(a,3i5)") ' Default n1q n2q n3q = ', n1q, n2q, n3q
    call readhamindex0()

    allocate(konf(0:lmxax,nbas), nncx(0:lmxax,nbas))
    isp  = 1                  ! konf is not spin-dependent
    konf = konft(:,:,isp)
    do ibas = 1, nbas
       write(6,"(i4,f9.4,100i4)") ibas, zz(ibas), lmxaa(ibas), konf(0:lmxaa(ibas),ibas)
       do lk = 0, lmxaa(ibas)
          nncx(lk, ibas) = konf(lk, ibas) - 1 - lk   ! number of cores per (l, ibas)
       enddo
    enddo

    allocate(nnvv(0:lmxax, nbas), source = 0)
    do izz = 1, mnla
       if (nnvv(lindx(izz), iat(izz)) < nindx(izz)) &
            nnvv(lindx(izz), iat(izz)) = nindx(izz)
    enddo

    ctrlg          = 'ctrlg.'//trim(sname)//'.toml'
    ctrlg_partial  = ctrlg//'.partial'
    pb_partial     = 'PB.'//trim(sname)//'.toml.partial'

    !! ===== Build ctrlg.<sname>.toml.partial =====
    !!   step 1: copy ctrl sections from existing ctrlg.<sname>.toml,
    !!           dropping any prior [gw]/[product_basis]/[blocks] sections.
    !!   step 2: append fresh [gw]/[product_basis]/[blocks] sections.
    call copy_ctrl_sections(ctrlg, ctrlg_partial, ifi)

    !! ----- [gw] -----
    write(ifi,'(a)') '[gw]'
    write(ifi,'(a,3(i0,a),a)') 'n1n2n3 = [', n1q, ', ', n2q, ', ', n3q, ']', '   # BZ mesh'
    write(ifi,'(a)') '# n1n2n3eps = [8, 8, 8]   # Optional finer mesh for eps/chipm/magnon (qg4gw --job=2/4/40001).'
    write(ifi,'(a)') '#                          # Silently ignored by GW (job=1) if not set.'
    write(ifi,'(a)') 'QpGcut_psi = 4.0    # |q+G| cutoff for eigenfunctions  (a.u. unless unit_2pioa=true)'
    write(ifi,'(a)') 'QpGcut_cou = 3.0    # |q+G| cutoff for Coulomb / W'
    write(ifi,'(a)') 'unit_2pioa = false  # false: a.u.; true: 2*pi/alat'
    write(ifi,'(a)') 'alpha_OffG = 1.0    # offset-Gamma auxiliary function'
    write(ifi,'(a)') '# emax_chi0 = 999.0   # (Ry) optional emax cutoff for chi0'
    write(ifi,'(a)') 'emax_sigm  = 3.0    # (Ry) emax cutoff for Sigma'
    write(ifi,'(a)')
    write(ifi,'(a)') '# ----- Frequencies -----'
    write(ifi,'(a)') 'niw           = 10      # # of frequencies along Im axis (try 6 or 12 for tests)'
    write(ifi,'(a)') 'HistBin_dw    = 1e-4    # bin width along real axis at omega=0'
    write(ifi,'(a)') 'HistBin_ratio = 1.05    # bin: frhis(iw) = (dw/(r-1))*(exp((r-1)*(iw-1)) - 1)'
    write(ifi,'(a)') '# SmearX0 = 0.01        # (Ha) Gaussian smear for X0/chi0 along freq (metals); 0/absent = off'
    write(ifi,'(a)') '# SmearX0q0 = 0.01      # (Ha) SmearX0 override applied at offset-Gamma q0 ONLY (interband-nesting safeguard)'
    write(ifi,'(a)') 'GaussSmear = true       # Gaussian smearing for poles of G^LDA in hsfp0'
    write(ifi,'(a)') 'deltaw = 0.02           # (a.u.) numerical-derivative mesh for Z factor'
    write(ifi,'(a)') 'esmr   = 0.003          # (Ry) hsfp0 smearing; keep < band gap for insulators'
    write(ifi,'(a)')
    write(ifi,'(a)') '# ----- Q for diagonal Sigma=GW (gw_lmfh) -----'
    write(ifi,'(a)') '# QforGWIBZ = true   # use IBZ instead of <QforGW>'
    write(ifi,'(a)') 'EMAXforGW = 15      # eV (above EFermi) bands cutoff for Sigma=GW'
    write(ifi,'(a)')
    write(ifi,'(a)') '# ----- Q for dielectric eps -----'
    write(ifi,'(a)') 'QforEPSau = true    # interpret <QforEPS> as a.u.'
    write(ifi,'(a)')
    write(ifi,'(a)') '# ----- Wannier (uncomment to use) -----'
    write(ifi,'(a)') '# wan_out_emin  = -1.05   # eV relative to EFermi'
    write(ifi,'(a)') '# wan_out_emax  =  2.4'
    write(ifi,'(a)') '# wan_maxit_1st = 300'
    write(ifi,'(a)') '# wan_conv_1st  = 1e-7'
    write(ifi,'(a)') '# wan_max_1st   = 0.1'
    write(ifi,'(a)') '# wan_maxit_2nd = 1500'
    write(ifi,'(a)') '# wan_max_2nd   = 0.3'
    write(ifi,'(a)') '# wan_conv_end  = 1e-8'
    write(ifi,'(a)')

    !! ----- [product_basis] (slim: only the user-tuned scalars) -----
    write(ifi,'(a)') '[product_basis]'
    write(ifi,'(a)') '# Tolerance to drop linearly-dependent products. Larger gives smaller PB.'
    write(ifi,'(a)') 'pb_tolerance = [1e-3]'
    write(ifi,'(a)')
    write(ifi,'(a)') '# lcutmx(atom): max l-cutoff for the product basis per atom.'
    write(ifi,'(a)') '# Use 4 for atoms with valence d (Ni, Ga, etc.), 6 for f.'
    allocate(lcutmx(nbas), source = 4)
    do ibas = 1, nbas
       if (zz(ibas) <  10.5d0)                          lcutmx(ibas) = 2
       if (57.001d0 < zz(ibas) .and. zz(ibas) < 71.001d0) lcutmx(ibas) = 6
       if (89.001d0 < zz(ibas))                         lcutmx(ibas) = 6
    enddo
    call write_int_vec(ifi, 'pb_lcutmx', lcutmx)
    write(ifi,'(a)')
    write(ifi,'(a)') '# Per-atom product-basis tables (nlx / valence / core) live in PB.<sname>.toml'
    write(ifi,'(a)')

    !! ----- [blocks] -----
    write(ifi,'(a)') '[blocks]'
    write(ifi,'(a)') 'QforEPS = """'
    write(ifi,'(a)') ' 0 0 0.00050'
    write(ifi,'(a)') ' 0 0 0.00100'
    write(ifi,'(a)') ' 0 0 0.00200'
    write(ifi,'(a)') '"""'
    write(ifi,'(a)')
    write(ifi,'(a)') '# QforEPSL = """'
    write(ifi,'(a)') '#  0 0 0   1 0 0  20'
    write(ifi,'(a)') '# """'
    write(ifi,'(a)')

    write(ifi,'(a)') 'QforGW = """'
    write(ifi,'(a)') ' 0.0 0.0 0.0'
    write(ifi,'(a)') ' 0.1 0.0 0.0'
    write(ifi,'(a)') ' 0.2 0.0 0.0'
    write(ifi,'(a)') ' 0.3 0.0 0.0'
    write(ifi,'(a)') '"""'
    write(ifi,'(a)')

    write(ifi,'(a)') '# Worb: atomic orbitals for MLWF / MLO modelling.'
    write(ifi,'(a)') '# Each row: <iatom> <label> <lm1> <lm2> ...'
    write(ifi,'(a)') '# lm index: 1=s, 2=py, 3=pz, 4=px, 5=xy, 6=yz, 7=3z^2-1, 8=xz, 9=x^2-y^2, ... (real harmonics)'
    write(ifi,'(a)') 'Worb = """'
    do ibas = 1, nbas
       write(ifi,'(a,i0,1x,a,a)') '! ', ibas, trim(spid(ibas)), '   1 2 3 4 5 6 7 8 9'
    enddo
    write(ifi,'(a)') '"""'

    close(ifi)

    !! atomically replace ctrlg.<sname>.toml with the new file
    call execute_command_line('mv '//ctrlg_partial//' '//ctrlg, wait=.true.)

    !! ===== Write PB.<sname>.toml (per-atom product-basis tables) =====
    open(newunit=ifpb, file=pb_partial)
    write(ifpb,'(a)') '# PB.<sname>.toml -- per-atom product basis tables.'
    write(ifpb,'(a)') '# Auto-generated by gwinit. Loaded by m_GWinput together with ctrlg.<sname>.toml.'
    write(ifpb,'(a)') '# (pb_tolerance and pb_lcutmx live in ctrlg.<sname>.toml; do not duplicate here.)'
    write(ifpb,'(a)')

    write(ifpb,'(a)') '[product_basis]'
    write(ifpb,'(a)')
    write(ifpb,'(a)') '# nlx: per (iatom, l) row [iatom, l, nnvv, nnc]'
    write(ifpb,'(a)') 'nlx = ['
    do ibas = 1, nbas
       do lk = 0, lmxaa(ibas)
          write(ifpb,'(a,4(i0,a))') '  [', ibas, ', ', lk, ', ', nnvv(lk, ibas), ', ', nncx(lk, ibas), '],'
       enddo
    enddo
    write(ifpb,'(a)') ']'
    write(ifpb,'(a)')

    !! valence rows
    write(ifpb,'(a)') '# valence: [iatom, l, n, occ, unocc]'
    write(ifpb,'(a)') 'valence = ['
    iatbk = 0
    do ibas = 1, nbas
       do lk = 0, lmxaa(ibas)
          do nx = 1, npqn
             do izz = 1, mnla
                if (iat(izz) == ibas .and. lk == lindx(izz) .and. nx == nindx(izz)) then
                   call valence_occ(zz(ibas), pqn(izz), lindx(izz), nindx(izz), nocc, nunocc)
                   write(ifpb,'(a,5(i0,a),a)') '  [', iat(izz), ', ', lindx(izz), ', ', &
                        nindx(izz), ', ', nocc, ', ', nunocc, '],', '   # '//trim(caption(izz))
                   exit
                endif
             enddo
          enddo
       enddo
    enddo
    write(ifpb,'(a)') ']'
    write(ifpb,'(a)')

    !! core rows
    write(ifpb,'(a)') '# core: [iatom, l, n, occ, unocc, forX0, forSxc]'
    write(ifpb,'(a)') 'core = ['
    do ibas = 1, nbas
       do lk = 0, lmxaa(ibas)
          if (lk == 0) seg2 = 'S'
          if (lk == 1) seg2 = 'P'
          if (lk == 2) seg2 = 'D'
          if (lk == 3) seg2 = 'F'
          if (lk == 4) seg2 = 'G'
          if (lk == 5) seg2 = 'H'
          if (lk == 6) seg2 = 'I'
          do kkk = lk + 1, konf(lk, ibas) - 1
             noccc   = 0
             nunoccc = 0
             ncinc   = 0
             ncinc2  = 0
             write(ifpb,'(a,7(i0,a),a)') '  [', ibas, ', ', lk, ', ', kkk - lk, ', ', &
                  noccc, ', ', nunoccc, ', ', ncinc, ', ', ncinc2, '],', &
                  '   # '//char(48 + kkk)//seg2
          enddo
       enddo
    enddo
    write(ifpb,'(a)') ']'

    close(ifpb)
    call execute_command_line('mv '//pb_partial//' PB.'//trim(sname)//'.toml', wait=.true.)

    stop ' OK! gwinit upserted [gw]/[product_basis]/[blocks] into '//ctrlg// &
         ' and wrote PB.'//trim(sname)//'.toml.'
  end subroutine gwinit_v2


  !> Copy the ctrl sections of a ctrlg.<sname>.toml into a fresh file,
  !  dropping any pre-existing [gw], [product_basis] and [blocks]
  !  sections so we can rewrite them. The fresh file is left open at
  !  the unit number `ifi_out`, positioned at end-of-file ready for
  !  appending the new GW sections.
  subroutine copy_ctrl_sections(src, dst, ifi_out)
    character(*), intent(in)  :: src, dst
    integer,      intent(out) :: ifi_out
    integer :: u_in, ios
    logical :: skip, src_exists
    character(len=4096) :: line
    inquire(file=src, exist=src_exists)
    if (.not. src_exists) call rx('gwinit: '//src//' not found '// &
         '(generate ctrl sections via ctrlgenToml.py first).')
    open(newunit=u_in, file=src, status='old', action='read')
    open(newunit=ifi_out, file=dst, status='replace', action='write')
    skip = .false.
    do
       read(u_in,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (line(1:5) == '[gw]'              .or. &
           line(1:16) == '[product_basis]' .or. &
           line(1:8) == '[blocks]') then
          skip = .true.
          cycle
       endif
       if (skip .and. len_trim(line) > 0 .and. line(1:1) == '[') then
          skip = .false.   ! a fresh non-target section starts: keep it
       endif
       if (.not. skip) write(ifi_out,'(a)') trim(line)
    enddo
    close(u_in)
  end subroutine copy_ctrl_sections


  subroutine write_int_vec(ifi, key, v)
    integer, intent(in) :: ifi, v(:)
    character(*), intent(in) :: key
    integer :: i
    write(ifi, '(a)', advance='no') key//' = ['
    do i = 1, size(v)
       if (i < size(v)) then
          write(ifi, '(i0,a)', advance='no') v(i), ', '
       else
          write(ifi, '(i0)', advance='no') v(i)
       endif
    enddo
    write(ifi, '(a)') ']'
  end subroutine write_int_vec


  !> Compute occ / unocc flags for valence rows from atomic Z and pqn / l.
  pure subroutine valence_occ(zc, pq, lx, nxx, nocc, nunocc)
    real(8),    intent(in)  :: zc
    integer,    intent(in)  :: pq, lx, nxx
    integer,    intent(out) :: nocc, nunocc
    nocc = 0
    nunocc = 1
    ! s
    if (lx == 0 .and. pq == 2 .and. zc > 1.5d0)  nocc = 1
    if (lx == 0 .and. pq == 3 .and. zc > 4.5d0)  nocc = 1
    if (lx == 0 .and. pq == 4 .and. zc > 12.5d0) nocc = 1
    if (lx == 0 .and. pq == 5 .and. zc > 30.5d0) nocc = 1
    if (lx == 0 .and. pq == 6 .and. zc > 48.5d0) nocc = 1
    if (lx == 0 .and. pq == 7 .and. zc > 80.5d0) nocc = 1
    if (lx == 0 .and. pq == 2 .and. zc > 10.5d0) nunocc = 0
    if (lx == 0 .and. pq == 3 .and. zc > 18.5d0) nunocc = 0
    if (lx == 0 .and. pq == 4 .and. zc > 36.5d0) nunocc = 0
    if (lx == 0 .and. pq == 5 .and. zc > 54.5d0) nunocc = 0
    if (lx == 0 .and. pq == 6 .and. zc > 86.5d0) nunocc = 0
    ! p
    if (lx == 1 .and. pq == 2 .and. zc > 1.5d0)  nocc = 1
    if (lx == 1 .and. pq == 3 .and. zc > 4.5d0)  nocc = 1
    if (lx == 1 .and. pq == 4 .and. zc > 12.5d0) nocc = 1
    if (lx == 1 .and. pq == 5 .and. zc > 30.5d0) nocc = 1
    if (lx == 1 .and. pq == 6 .and. zc > 48.5d0) nocc = 1
    if (lx == 1 .and. pq == 7 .and. zc > 80.5d0) nocc = 1
    if (lx == 1 .and. pq == 2 .and. zc > 10.5d0) nunocc = 0
    if (lx == 1 .and. pq == 3 .and. zc > 18.5d0) nunocc = 0
    if (lx == 1 .and. pq == 4 .and. zc > 36.5d0) nunocc = 0
    if (lx == 1 .and. pq == 5 .and. zc > 54.5d0) nunocc = 0
    if (lx == 1 .and. pq == 6 .and. zc > 86.5d0) nunocc = 0
    ! d
    if (lx == 2 .and. pq == 3 .and. zc > 20.5d0) nocc = 1
    if (lx == 2 .and. pq == 4 .and. zc > 38.5d0) nocc = 1
    if (lx == 2 .and. pq == 5 .and. zc > 56.5d0) nocc = 1
    if (lx == 2 .and. pq == 3 .and. zc > 30.5d0) nunocc = 0
    if (lx == 2 .and. pq == 4 .and. zc > 48.5d0) nunocc = 0
    ! f
    if (lx == 3 .and. pq == 4 .and. zc > 57.5d0) nocc = 1
    if (lx == 3 .and. pq == 5 .and. zc > 89.5d0) nocc = 1
    if (lx == 2 .and. pq == 4 .and. zc > 71.5d0) nunocc = 0
    ! g+
    if (lx > 3) then
       nocc   = 0
       nunocc = 0
    endif
    ! phidot skip
    if (nxx == 2) then
       nocc   = 0
       nunocc = 0
    endif
  end subroutine valence_occ

end module m_gwinit
