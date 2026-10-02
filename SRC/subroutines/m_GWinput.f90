!> m_GWinput: single source of truth for the GW-side input (ctrlg.<sname>.toml).
!  Loads once via toml-f, exposes all values as protected module variables.
!  Callers `use m_GWinput, only: niw, deltaw, ...` to access values.
!
!  This replaces scattered `getkeyvalue("GWinput", key, var)` calls — all
!  GWinput keys are loaded eagerly at startup, then read-only thereafter.
!
!  Schema based on survey of 66 Samples GWinput files.
!
!  ---------------------------------------------------------------------------
!  THE INPUT IS ctrlg.<sname>.toml. NOTHING ELSE.
!  ---------------------------------------------------------------------------
!  There is no GWinput.toml. That name was an intermediate of the 2026-05
!  migration and nothing in ecalj reads it; a file by that name in a run
!  directory is dead and can be moved to trash. The legacy plain-text GWinput
!  is not read either (the reader is disabled and aborts with a message).
!
!  Converting a legacy GWinput by hand (what Legacy2toml.py does, written out
!  here so it can be done without the script if it ever bit-rots):
!
!    GWinput line              -> ctrlg.<sname>.toml
!    ---------------------------------------------------------------------
!    scalar/vector keys        -> [gw]      key = value   (same names, TOML
!                                 syntax: 1d-7 -> 1e-7, "on"/"off" -> true/false,
!                                 "n1n2n3 4 4 4" -> n1n2n3 = [4, 4, 4])
!    mlo_*  (any)              -> [mlo]     key = value   (MLO model)
!    wan_*  (any)              -> [gw]      only needed by hmaxloc (cRPA, magnon);
!                                 comment them out otherwise
!    <PRODUCT_BASIS> tolerance -> [product_basis] pb_tolerance = [..]
!    <PRODUCT_BASIS> lcutmx    -> [product_basis] pb_lcutmx    = [..]
!    <PRODUCT_BASIS> nlx / valence / core tables
!                              -> [product_basis] nlx = [[iatom,l,nnvv,nnc],...],
!                                 valence = [[iatom,l,n,occ,unocc],...],
!                                 core = [[iatom,l,n,occ,unocc,forX0,forSxc],...]
!                                 (one inner list per row; the section is the
!                                 last one in the file. Until 2026-09 these three
!                                 were a separate PB.<sname>.toml; that file is
!                                 not read any more -- ctrlg_absorb.py folds it in.)
!    <QforEPS> <QforGW>        -> [gw]      NAME = """<raw lines verbatim>"""
!    <Worb>                    -> [mlo]     mlo_lm = """<raw lines verbatim>"""
!    <QPNT> <QforEPSL> <hrotr> -> [blocks]  NAME = """<raw lines verbatim>"""
!
!    renamed / retired on the way:
!      zmel_max_size     -> zmel_batch_gb
!      MEMnmbatch        -> dropped (different meaning)
!      GaussSmear, dw, omg_c, delta, WgtQ0P, tetrakbt -> removed 2026-09-19 (dead keys; the logical tetrakbt
!                           was replaced by t_tetrakbt > 0). Left in a file they are ignored, as is any key
!                           that this module does not ask for (there is no check for unknown keys).
!      SmearX0, SmearX0q0 -> t_tetrakbt = -T (2026-09-28). These two are NOT ignored: load_gw_section stops
!                           when either is present (ctrlg_update.py converts), and it stops when
!                           t_tetrakbt is missing.
!      esmr, t_sigmakbt  -> t_sigmaw (2026-09-20); both are still read when t_sigmaw is absent.
!                           (Comment corrected 2026-09-30 01:20: it listed SmearX0q0 among the ignored keys.)
!      mlo_emax          -> retired (mlo_method=4 does not use it)
!      esm_input.dat     -> [esm] section (see m_lmfinit; migrated in place)
!
!  Two files that look like GW input but are not: ctrls.<sname> is the seed
!  for ctrlgenToml.py, and GWinput in Samples/Legacy/ is the converter's input.
!  ---------------------------------------------------------------------------
!
module m_GWinput
  use tomlf, only: toml_table, toml_array, toml_error, &
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
  ! alpha_OffG: legacy default -1d60 (sentinel "not given" -> fall through to alpha_OffG_vec)
  real(8), protected, public :: alpha_OffG   = -1.0d60
  logical, protected, public :: unit_2pioa   = .false.

  ! Sigma / chi0 mode
  integer, protected, public :: iSigMode     = 3
  integer, protected, public :: niw          = 10
  integer, protected, public :: nband_chi0   = 99999  ! legacy extension.f90 default; effectively "all bands"
  ! EMINforGW/EMAXforGW: legacy reads as REAL with default ±99999d0
  real(8), protected, public :: EMINforGW    = -99999.0d0
  real(8), protected, public :: EMAXforGW    = 99999.0d0
  real(8), protected, public :: emax_sigm    = 3.0d0
  real(8), protected, public :: emax_chi0    = 999.0d0
  real(8), protected, public :: HistBin_ratio = 1.03d0
  real(8), protected, public :: HistBin_dw   = 1.0d-5  ! legacy m_freq default
  real(8), protected, public :: deltaw       = 0.02d0
  real(8), protected, public :: SmearX0      = 0.0d0   ! (Ha) Gaussian smearing of Im chi0 along omega (dpsion5); 0=off. Not an input key
                                                      ! since 2026-09-28: set from t_tetrakbt < 0 (see there).

  ! Optional flags
  logical, protected, public :: KeepEigen    = .true.
  logical, protected, public :: KeepCMLO     = .true.
  logical, protected, public :: KeepPPOVL    = .false.
  logical, protected, public :: NormChk      = .false.
  logical, protected, public :: AnyQ         = .false.
  logical, protected, public :: QforEPSau    = .false.
  logical, protected, public :: QforEPSunita = .false.
  logical, protected, public :: QforEPSLIncLeft = .false.
  ! t_tetrakbt (K, required in [gw] since 2026-09-28): how chi0 is smeared, by its sign.
  !   0  = none (T=0 tetrahedron lindtet6; insulators and semiconductors).
  !   >0 = finite-T tetrahedron (method B'): Fermi-Dirac occupations at this T (lindtet6_kbt), and heftet writes
  !        EFERMI_kbt, the Fermi level of chi0 and Sigma.
  !   <0 = T=0 tetrahedron, Im chi0 smeared along omega by a Gaussian as wide as the Fermi-Dirac distribution at |T|
  !        (std = pi kB|T|/sqrt3; -1000 K = 0.00574 Ha), i.e. the former [gw] SmearX0, which is no longer read.
  !   (The old logical 'tetrakbt' is gone.)
  real(8), protected, public :: t_tetrakbt   = 0.0d0
  ! t_sigmaw: width (K) of the Fermi-Dirac kernel that smears the intermediate levels of the
  !   self-energy (Sigma_x = Gv and the pole term of Sigma_c = G(W-v)).  A numerical width, not a
  !   full finite-T self-energy.  Default 1000 K ~ the former Gaussian esmr = 0.01 Ry.  Metals may
  !   want 300-1000 K together with a finer k mesh.  The Fermi level is EFERMI_kbt when
  !   t_tetrakbt > 0 (chi0 also at finite T), EFERMI otherwise.
  !   (2026-09-20: replaces esmr [Ry, Gaussian] and t_sigmakbt; both are read for compatibility.)
  real(8), protected, public :: t_sigmaw     = 1000.0d0
  real(8), parameter :: kb_ev = 8.6171d-5, rydberg_ev = 13.605693d0   ! for the esmr -> t_sigmaw conversion
  ! wcsmear (default true): in the real-axis pole term of Sigma_c, apply the Fermi-Dirac smearing of
  !   the intermediate level (t_sigmaw) to W_c(omega) itself instead of evaluating W_c at the mean
  !   energy (false = the old 3-point interpolation at the mean energy).  Removes the knife-edge sensitivity to sharp plasmon poles of W
  !   (LiTi2O4 at <=1000 K, MD/research_kotani_log.md 2026-09-19).  Unchanged where W_c is smooth.
  logical, protected, public :: wcsmear      = .true.
  ! chi0_skip_window = [emin, emax] (eV, relative to E_F): band pairs whose occupied AND unoccupied
  !   states both lie inside the window are left out of chi0 (cRPA-like removal of the intraband /
  !   low-energy transitions of a partially filled band).  Empty (default) = off.
  real(8), protected, public :: chi0_skip_window(2) = 0d0
  logical, protected, public :: chi0_skip_window_set = .false.
  ! chi0_filterw = [wc, dw] (eV): the tetrahedron weight of every band pair is multiplied, bin by bin,
  !   by the Fermi-like step 1/(1+exp((wc-|omega|)/dw)) when chi0 is accumulated (x0kf_v4h), i.e. the
  !   transitions below wc are removed from chi0 (omega-space cRPA).  Empty (default) = off.
  ! chi0_filterw_drude (default true): the intraband pairs (same band, occupied part -> unoccupied part
  !   across E_F; the Drude term, same test as --intrabandonly of the eps code) are exempt from the
  !   filter, so the static metallic screening and the Drude weight stay.  false: filtered as well.
  !   Meaningless without chi0_filterw.  Experimental (2026-09-20, MD/research_kotani_log.md 23:26).
  real(8), protected, public :: chi0_filterw(2) = 0d0
  logical, protected, public :: chi0_filterw_set = .false.
  logical, protected, public :: chi0_filterw_drude = .true.
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
  ! mlo_method: how the weight theta in A = S*theta is built (see m_hreduction).
  !   4 (default, recommended) ecut_j = max(ecbot + mlo_delta, eps^MTO_j), one sigmoid.
  !                            No per-material input; mlo_emax is ignored.
  !   0/1/2  legacy, one sigmoid with a hand-set mlo_emax as the floor (0), the only
  !          cut (1), or no floor at all (2, the orbital's own energy).
  !   3      two-stage sigmoid; superseded by 4 (kept for comparison).
  integer, protected, public :: mlo_method   = 4
  ! mlo_emax: used by mlo_method 0 and 1 only. Absent -> evl(ndimMTO+nskip), which is
  !   an accident of band counting rather than a window; that is why 0/1 needed a
  !   hand-set value per material and why 4 exists.
  real(8), protected, public :: mlo_emax     = huge(0d0)  ! sentinel: key absent ⇒ runtime default
! integer, protected, public :: wan_maxit_1st = 100 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_maxit_2nd = 100 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_tb_cut   = 1.01d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_conv_1st = 1.0d-5 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_conv_end = 1.0d-8 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_max_1st  = 0.1d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_max_2nd  = 0.3d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
  ! wan_*_emin/emax sentinels: legacy uses 999/-999 to force user spec via leout/lein gate.
! real(8), protected, public :: wan_in_emin  = 999.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_in_emax  = -999.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_out_emin = 999.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_out_emax = -999.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_in_ewin  = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)

  !-----------------------------------------------------------------
  ! Additional keys (batch 2 -- 2026-05-02 migration of remaining callers)
  ! Defaults reflect the most common legacy `default=` value at the
  ! callsite. Where the legacy default is non-constant (e.g., default=nnn),
  ! the value is set to a sentinel and the caller emulates the fallback.
  !-----------------------------------------------------------------

  ! Reals
  real(8), protected, public :: BZadiv             = 1.0d0
  real(8), protected, public :: ene_sppola         = 0.0d0
  ! All mlo_* energies are in eV, like mlo_emax.
  ! mlo_delta: how far ABOVE the band edge (the global CBM; EF in metals) the model
  !   is required to be accurate. mlo_method=4 puts its floor at ecbot+mlo_delta.
  !   This is a STATEMENT OF WHAT YOU WANT, not a fitting parameter: set it to the
  !   top of the energy window you care about, and evaluate in that same window.
  ! mlo_w: width of the sigmoid that falls off above that floor. THIS is the knob
  !   to turn when the residual is too large. Measured optima (eV): semiconductors
  !   ~1.8-2.0, single-atom transition metals Fe/Cu ~11, RuO2 ~1.8.
  ! Defaults 2.0/2.0 eV beat the earlier 2.45/2.72 on the sample set
  !   (window rms 16.9 -> 15.0 meV, occupied side 11.6 -> 7.4 meV, |gap| 16 -> 13 meV).
  ! mlo_w default is METHOD-DEPENDENT and resolved in m_hreduction:
  !   mlo_method=4 -> 2.0 eV (the value optimized for it)
  !   mlo_method=0/1/2/3 -> 2.7211 eV (= 0.2 Ry, the historical value, so that the
  !     shipped samples and their reference band files are reproduced exactly).
  real(8), protected, public :: mlo_w            = huge(0d0) ! sentinel: key absent
  real(8), protected, public :: mlo_delta           = 2.0d0    ! eV
  real(8), protected, public :: mlo_wfrz           = 1.36d0   ! eV: freeze-edge width, mlo_method=3 only
  real(8), protected, public :: mlo_down           = 0.0d0    ! eV: own-energy floor lift, mlo_method=3 only
  ! mlo_low / mlo_wlow (hidden, any mlo_method): also cut the model OFF BELOW an
  ! energy. theta_j gets an extra factor sigma((mlo_low - eps)/mlo_wlow), so PMT
  ! states below EF + mlo_low (eV) stop pulling on the fit. Off unless mlo_low is
  ! set; mlo_wlow defaults to mlo_w. Meant for e.g. a d-only model of Cu, where
  ! the s band bottom far below the d band otherwise enters with weight 1.
  ! (2026-09-17: tried on Cu's d model, kept as a record, not active)
  !real(8), protected, public :: mlo_low            = huge(0d0) ! eV rel. EF; huge = off
  !real(8), protected, public :: mlo_wlow           = huge(0d0) ! eV; huge = same as mlo_w
  ! mlo_pcut / mlo_pw (hidden): cut by CHARACTER instead of energy. PMT state i
  ! enters the fit with the extra factor sigma((p_i - mlo_pcut)/mlo_pw), where
  ! p_i = sum_j |<Psi^PMT_i|Psi^MTO_j>|^2 is its weight in the model's MTO
  ! subspace (0..1). States that are not model-like (an s band crossing a d-only
  ! model) drop out wherever they are in energy. Off unless mlo_pcut is set.
  ! (not tried yet; kept as a record, not active)
  !real(8), protected, public :: mlo_pcut           = huge(0d0) ! 0..1; huge = off
  !real(8), protected, public :: mlo_pw             = 0.1d0     ! width of the sigmoid in p
  real(8), protected, public :: mixbeta            = 1.0d0
  real(8), protected, public :: mixtj              = 0.0d0
  real(8), protected, public :: TFscreen           = 1.0d-5**0.5d0
  real(8), protected, public :: removed_r0c        = 1.0d60
  real(8), protected, public :: q0scale            = 0.8d0
  real(8), protected, public :: shift_majority     = 0.0d0
  real(8), protected, public :: output_ddmat_atom  = 1.0d0
  real(8), protected, public :: dRdIatRmax         = 0.003d0
  real(8), parameter,  public :: zmel_batch_gb_min     = 0.4d0  ! minimum zmel batch, CPU (GB/rank)
  real(8), parameter,  public :: zmel_batch_gb_min_gpu = 2.0d0  ! minimum zmel batch, GPU (GB/rank)
  real(8), protected, public :: zmel_batch_gb         = zmel_batch_gb_min
  real(8), protected, public :: magnon_delta       = 0.0d0
  real(8), protected, public :: magnon_delta_dos   = 1.0d-6
  real(8), protected, public :: magnon_HistBin_ratio = 1.03d0
  real(8), protected, public :: magnon_HistBin_dw  = 1.0d-5
  ! mlo
! real(8), protected, public :: mlo_conv           = 1.0d-6 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_mix            = 0.5d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_EUinner        = 1.0d8 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_CUouter        = 0.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_CUinner        = 0.9d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_WTinner        = 2048.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_WTband         = 64.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_WTseed         = 32.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_ELinner        = -1.0d8 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_ewid           = 1.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_WTouter        = 32768.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_CLhard         = 0.33d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: mlo_ELhard         = -1.0d8 !removed 2026-10-02 (the Wannier functions / lmfham2)
  ! wmat
  real(8), protected, public :: wmat_rcut1         = 0.01d0
  real(8), protected, public :: wmat_rcut2         = 0.01d0
  ! wan
! real(8), protected, public :: wan_mix_1st        = 0.1d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_mix_2nd        = 0.1d0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_conv_2nd       = 1.0d-5 !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_tbcut_rcut     = -1.0d50  ! sentinel; default=rcut !removed 2026-10-02 (the Wannier functions / lmfham2)
! real(8), protected, public :: wan_tbcut_heps     = 0.0d0 !removed 2026-10-02 (the Wannier functions / lmfham2)

  ! Integers
  integer, protected, public :: ngcell             = 1
  integer, protected, public :: nkeep_wfs          = 2
  integer, protected, public :: mlo_nskip          = -huge(0)  ! sentinel
  integer, protected, public :: nbcutlow_sig       = 0
! integer, protected, public :: mlo_maxit          = 100 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_nb_below       = 0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_nb_above       = 0 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_out_bmin       = 999 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_out_bmax       = -999 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_in_bmin        = 999 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_in_bmax        = -999 !removed 2026-10-02 (the Wannier functions / lmfham2)
  integer, protected, public :: mixpriorit         = 3
  integer, protected, public :: Q0Pchoice          = 1
  real(8), protected, public :: DeltaQscale        = -1.0d0  ! >0: directly set offset-Gamma deltaq_scale; <0: use Q0Pchoice (1->0.1, 2->1/sqrt3)
  integer, protected, public :: Verbose            = 0
  integer, protected, public :: Q0P_Choice         = 0
  ! NormChk_int: switch.f90 reads NormChk as integer (default=1) — keep
  ! both forms (logical NormChk above for default=.false. callsite,
  ! plus integer NormChk_int derived from the same TOML key).
  integer, protected, public :: NormChk_int        = 1

  ! Logicals
  logical, protected, public :: EIBZmode           = .true.
  logical, protected, public :: QforEPSIBZ         = .false.
  logical, protected, public :: QforGWIBZ          = .false.
  logical, protected, public :: TestOnlyQ0P        = .false.
  logical, protected, public :: TestNoQ0P          = .false.
  logical, protected, public :: NoQ0P              = .false.
  logical, protected, public :: KeepPpb            = .false.
#ifdef __GPU
  logical, protected, public :: KeepWV             = .true.
#else
  logical, protected, public :: KeepWV             = .false.
#endif
  logical, protected, public :: KeepQG             = .true.
  logical, protected, public :: TimeReversal       = .true.
  logical, protected, public :: rmeshrefine        = .true.
  logical, protected, public :: chi_RegQbz         = .true.
  logical, protected, public :: tetrahedron_matrix_linear = .false.
  logical, protected, public :: magnon_w_onsite_dddd = .true.
  logical, protected, public :: magnon_negative_cut = .false.
  logical, protected, public :: allq0i             = .false.
! logical, protected, public :: wan_out_ewin       = .true. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_in_bwin        = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
  logical, protected, public :: wmat_static        = .false.
  logical, protected, public :: wmat_all           = .false.
  logical, protected, public :: wmat_WSsuper       = .true.
! logical, protected, public :: wan_gauss_head     = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_truncate       = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: mlo_EUinnerAUTOsp  = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_out_emax_auto  = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_in_emax_auto   = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! logical, protected, public :: wan_small_ham      = .false. !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_nsh1           = 1 !removed 2026-10-02 (the Wannier functions / lmfham2)
! integer, protected, public :: wan_nsh2           = 2 !removed 2026-10-02 (the Wannier functions / lmfham2)
  ! MPI layout overrides: 0 = auto (MPI__AutoSetup decides)
  integer, protected, public :: mpi_worker_exch    = 0
  integer, protected, public :: mpi_worker_corr    = 0

  ! Vectors of 3
  integer, protected, public :: n1n2n3dos(3)       = [0, 0, 0]
  integer, protected, public :: GammaDivn1n2n3(3)  = [0, 0, 0]
  real(8), protected, public :: alpha_OffG_vec(3)  = [-1.0d50, 0.0d0, 0.0d0]
  real(8), protected, public :: wmat_rsite(3)      = [0.0d0, 0.0d0, 0.0d0]

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
  character(len=:), protected, public, allocatable :: block_Worb2
  character(len=:), protected, public, allocatable :: block_Worb3
  character(len=:), protected, public, allocatable :: block_hrotr

  !-----------------------------------------------------------------
  ! [blocks] structured form -- parsed from raw text in load_blocks_section.
  ! Callers read these directly instead of opening a unit on raw text.
  !-----------------------------------------------------------------

  ! QforEPS: list of q-vectors, 3 reals per line.
  real(8), protected, public, allocatable :: q_eps(:,:)         ! (3, n_eps)
  integer, protected, public               :: n_eps = 0

  ! QforEPSL: 7 fields per line -- q(3), qend(3), idx(int).
  real(8), protected, public, allocatable :: q_epsl(:,:)        ! (3, n_epsl)
  real(8), protected, public, allocatable :: qend_epsl(:,:)     ! (3, n_epsl)
  integer, protected, public, allocatable :: idx_epsl(:)        ! (n_epsl)
  integer, protected, public               :: n_epsl = 0

  ! QforGW: list of q-vectors, 3 reals per line.
  real(8), protected, public, allocatable :: q_qgw(:,:)         ! (3, n_qgw)
  integer, protected, public               :: n_qgw = 0

  ! QPNT: structured.
  integer, protected, public               :: qpnt_allq      = 0
  integer, protected, public               :: qpnt_spinonly  = 0
  integer, protected, public               :: qpnt_nstates   = 0
  integer, protected, public, allocatable  :: qpnt_bands(:)
  integer, protected, public               :: qpnt_nq        = 0
  real(8), protected, public, allocatable  :: qpnt_q(:,:)        ! (3, qpnt_nq)

  ! mlo_lm (formerly Worb): list of records: iatom, label(8), lm(1..nlm).
  integer, protected, public               :: n_worb         = 0
  integer, protected, public, allocatable  :: worb_iatom(:)
  character(8), protected, public, allocatable :: worb_label(:)
  integer, protected, public, allocatable  :: worb_lm(:,:)       ! (16, n_worb), -999 for unused slots
  integer, protected, public, allocatable  :: worb_nlm(:)        ! actual count per row
  ! mlo_lm2: the same, but for the SECOND radial set (EH2) of the PMT basis.
  ! Optional; absent means "no EH2 channel in the MLO model", which is the
  ! historical behaviour (m_HamPMT skipped k_table==2 unconditionally).
  integer, protected, public               :: n_worb2        = 0
  integer, protected, public, allocatable  :: worb2_iatom(:)
  character(8), protected, public, allocatable :: worb2_label(:)
  integer, protected, public, allocatable  :: worb2_lm(:,:)      ! (16, n_worb2)
  integer, protected, public, allocatable  :: worb2_nlm(:)
  ! mlo_lm3: the same, for the local-orbital (PZ, k=3) channels.  Listing a
  ! channel here takes the LO IN ADDITION to the EH function; without it the
  ! historical either/or applies (LO replaces EH when use_lo, else LO dropped).
  integer, protected, public               :: n_worb3        = 0
  integer, protected, public, allocatable  :: worb3_iatom(:)
  character(8), protected, public, allocatable :: worb3_label(:)
  integer, protected, public, allocatable  :: worb3_lm(:,:)
  integer, protected, public, allocatable  :: worb3_nlm(:)

  !-----------------------------------------------------------------
  ! State
  !-----------------------------------------------------------------
  logical, protected, public :: gwinput_loaded = .false.

  public :: gwinput_load, gwinput_init, gwinput_available

contains

  !> Idempotent helper for callers: ensure GW input has been loaded.
  !  Reads ctrlg.<sname>.toml: [gw], [mlo], [blocks] and [product_basis]
  !  (cut-offs plus the per-atom tables nlx/valence/core). Legacy
  !  ctrl/GWinput must be pre-converted via Legacy2toml.py first.
  !> True when the GW-side input (ctrlg.<sname>.toml) is present.
  !  For callers that must work with and without GW input (verbose(), BZadiv, ...):
  !  gwinput_init() aborts when the file is missing, this only reports.
  logical function gwinput_available()
    use m_ext, only: sname
    inquire(file='ctrlg.'//trim(sname)//'.toml', exist=gwinput_available)
  end function gwinput_available

  subroutine gwinput_init()
    use m_ext, only: sname
    logical :: have_ctrlg
    character(len=:), allocatable :: errmsg
    if (gwinput_loaded) return
    inquire(file='ctrlg.'//trim(sname)//'.toml', exist=have_ctrlg)
    if (.not. have_ctrlg) call rx('m_GWinput: ctrlg.'//trim(sname)// &
         '.toml not found in cwd (run Legacy2toml.py to convert legacy inputs).')
    call gwinput_load(error=errmsg)
    if (.not. gwinput_loaded) call rx('m_GWinput: failed to load ctrlg.'//trim(sname)//'.toml')
  end subroutine gwinput_init

  subroutine gwinput_load(filename, error)
    !> Load ctrlg.<sname>.toml. Idempotent.
    !  filename defaults to 'ctrlg.<sname>.toml' if absent.
    !
    !  Sections consumed:
    !    [gw]              -- run-level scalars (n1n2n3, QpGcut_*, etc.) + QforEPS/QforGW
    !    [mlo]             -- mlo_method/mlo_delta/mlo_w + mlo_lm
    !    [blocks]          -- raw text blocks (QPNT, QforEPSL, hrotr)
    !    [product_basis]   -- pb_tolerance, pb_lcutmx and the per-atom tables
    !                         nlx / valence / core (before 2026-09 those three
    !                         sat in a separate PB.<sname>.toml; such a file
    !                         is not read: abort -> ctrlg_absorb.py).
    use m_ext, only: sname
    character(*),               intent(in),  optional :: filename
    character(len=:), allocatable, intent(out), optional :: error

    type(toml_table), allocatable, target :: root
    type(toml_table), pointer :: gw, pb, blocks, mlo
    type(toml_error), allocatable :: terr
    character(len=:), allocatable :: fname

    if (gwinput_loaded) return
    if (present(error)) error = ""

    if (present(filename)) then
       fname = trim(filename)
    else
       fname = 'ctrlg.'//trim(sname)//'.toml'
    endif

    block
      use m_toml_override, only: load_toml_with_overrides
      use tomlf, only: toml_loads
      character(len=:), allocatable :: text
      call load_toml_with_overrides(fname, text)
      call toml_loads(root, text, error=terr)
    end block
    if (allocated(terr)) then
       if (present(error)) error = "m_GWinput: parse error in "//fname//": " // terr%message
       return
    endif

    !---- [gw] ----
    call get_value(root, 'gw', gw)
    if (associated(gw)) call load_gw_section(gw)

    !---- [product_basis]: cut-offs + per-atom tables nlx / valence / core ----
    call get_value(root, 'product_basis', pb)
    if (associated(pb)) then
       call load_pb_section(pb)
       call load_pb_tables(pb)
    endif
    if (pb_n_nlx == 0) call refuse_pb_file('PB.'//trim(sname)//'.toml')   ! pre-2026-09 layout: not read

    !---- [blocks] ----
    call get_value(root, 'blocks', blocks)
    if (associated(blocks)) call load_blocks_section(blocks)

    !---- [mlo] (mlo_* keys + mlo_lm; overrides anything left in [gw]/[blocks]) ----
    call get_value(root, 'mlo', mlo, requested=.false.)   ! do not create it when absent
    if (associated(mlo)) call load_mlo_section(mlo)
    if (allocated(block_QforEPS)) call parse_qvec_list(block_QforEPS, q_eps, n_eps)
    if (allocated(block_QforGW))  call parse_qvec_list(block_QforGW,  q_qgw, n_qgw)
    if (allocated(block_Worb))    call parse_Worb(block_Worb)
    if (allocated(block_Worb2))   call parse_Worb2(block_Worb2)
    if (allocated(block_Worb3))   call parse_Worb3(block_Worb3)

    gwinput_loaded = .true.
  end subroutine gwinput_load


  subroutine load_gw_section(gw)
    use m_mpi, only: ipr
    use m_lgunit, only: stdo
    type(toml_table), pointer, intent(in) :: gw
    logical :: ok
    ! Integer scalars
    call gv_i(gw, 'iSigMode',      iSigMode)
    call gv_i(gw, 'niw',           niw)
    call gv_i(gw, 'nband_chi0',    nband_chi0)
    call gv_r(gw, 'EMINforGW',     EMINforGW)
    call gv_r(gw, 'EMAXforGW',     EMAXforGW)
    call gv_i(gw, 'BZmesh',        BZmesh)
    if (.not. gw%has_key('t_tetrakbt')) call rx('m_GWinput: [gw] t_tetrakbt is required (2026-09-28). '// &
         '0 = no smearing of chi0 (insulators, semiconductors); T > 0 = finite-T tetrahedron at T (K); '// &
         '-T = Im chi0 smeared by a Gaussian as wide as the Fermi-Dirac distribution at T (the former SmearX0; '// &
         'metals, e.g. -1000). ctrlg_update.py adds it (0, or -T from SmearX0).')
    call gv_r(gw, 't_tetrakbt',    t_tetrakbt)
    block   ! t_sigmaw, with the retired esmr (Ry, Gaussian) and t_sigmakbt (K) read for compatibility
      real(8) :: tsk, esm
      integer :: st1, st2, st3
      call get_value(gw, 't_sigmaw', t_sigmaw, stat=st1)
      tsk = -1d0; call get_value(gw, 't_sigmakbt', tsk, stat=st2)
      esm = -1d0; call get_value(gw, 'esmr', esm, stat=st3)
      if (st2 == 0 .and. tsk > 0d0) then
        if (st1 /= 0) t_sigmaw = tsk
        if (ipr) write(stdo,"(a)") ' m_GWinput: t_sigmakbt is retired -> use t_sigmaw (K); value copied'
      endif
      if (st3 == 0 .and. st1 /= 0 .and. (st2 /= 0 .or. tsk <= 0d0) .and. esm > 0d0) then
        t_sigmaw = esm*rydberg_ev/(1.81d0*kb_ev)   ! Gaussian sigma (Ry) -> kBT with the same std (FD std = pi kBT/sqrt3 = 1.81 kBT)
        if (ipr) write(stdo,"(a,f8.1,a)") ' m_GWinput: esmr is retired (Gaussian smearing); t_sigmaw =', &
             t_sigmaw, ' K set from it. Please write t_sigmaw in [gw].'
      elseif (st3 == 0) then
        if (ipr) write(stdo,"(a)") ' m_GWinput: esmr is retired and ignored (the Sigma kernel is Fermi-Dirac of width t_sigmaw)'
      endif
      if (t_sigmaw < 0d0) call rx('m_GWinput: t_sigmaw must be >= 0 (K)')
    end block
    call gv_l(gw, 'wcsmear',       wcsmear)
    block
      type(toml_array), pointer :: arr
      call get_value(gw, 'chi0_skip_window', arr, requested=.false.)
      if (associated(arr)) then
        if (len(arr) /= 2) call rx('m_GWinput: chi0_skip_window must be [emin, emax] (eV)')
        call get_value(arr, 1, chi0_skip_window(1)); call get_value(arr, 2, chi0_skip_window(2))
        chi0_skip_window_set = .true.
      endif
      call get_value(gw, 'chi0_filterw', arr, requested=.false.)
      if (associated(arr)) then
        if (len(arr) /= 2) call rx('m_GWinput: chi0_filterw must be [wc, dw] (eV)')
        call get_value(arr, 1, chi0_filterw(1)); call get_value(arr, 2, chi0_filterw(2))
        if (chi0_filterw(2) <= 0d0) call rx('m_GWinput: chi0_filterw: dw must be > 0')
        chi0_filterw_set = .true.
      endif
      call gv_l(gw, 'chi0_filterw_drude', chi0_filterw_drude)
    end block
    ! QforEPS / QforGW: q-point lists for eps and for the one-shot GW driver.
    ! They live in [gw] since 2026-09-17; a copy left under [blocks] is still read.
    ok = take_block(gw, 'QforEPS', block_QforEPS)
    ok = take_block(gw, 'QforGW',  block_QforGW)
    ! mlo_* used to live in [gw]. They now belong to [mlo] (see load_mlo_section);
    ! a [gw] that still carries them is read the old way, with a one-line notice,
    ! and any [mlo] value loaded afterwards overrides it.
    call load_mlo_keys(gw, legacy=.true.)
!   call gv_i(gw, 'wan_maxit_1st', wan_maxit_1st) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_maxit_2nd', wan_maxit_2nd) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_tb_cut',    wan_tb_cut)      ! legacy reads as REAL (default 1.01) !removed 2026-10-02 (the Wannier functions / lmfham2)
    ! nband_sigm: legacy reads as integer; TOML may have list[float] -- take first as int
    call gv_iv_first(gw, 'nband_sigm', nband_sigm)
    ! MagAtom: VLA -- accept scalar or vector
    call gv_iv_alloc(gw, 'MagAtom', MagAtom)

    ! Real scalars
    call gv_r(gw, 'QpGcut_psi',    QpGcut_psi)
    call gv_r(gw, 'QpGcut_cou',    QpGcut_cou)
    call gv_r(gw, 'alpha_OffG',    alpha_OffG)
    ! emax_sigm/emax_chi0: legacy reads as scalar; TOML may have list -- take first.
    call gv_rv_first(gw, 'emax_sigm',  emax_sigm)
    call gv_rv_first(gw, 'emax_chi0',  emax_chi0)
    call gv_r(gw, 'HistBin_ratio', HistBin_ratio)
    call gv_r(gw, 'HistBin_dw',    HistBin_dw)
    call gv_r(gw, 'deltaw',        deltaw)
    if (gw%has_key('SmearX0') .or. gw%has_key('SmearX0q0')) call rx('m_GWinput: [gw] SmearX0 is no longer read '// &
         '(2026-09-28): write t_tetrakbt = -T instead, T such that pi kB T/sqrt3 is the Gaussian width '// &
         '(SmearX0 = 0.0057 Ha <-> t_tetrakbt = -992; -1000 K = 0.00574 Ha). ctrlg_update.py converts it.')
    ! t_tetrakbt < 0: the Gaussian smearing of Im chi0 (dpsion5), std = that of the Fermi-Dirac distribution at |T|
    if (t_tetrakbt < 0d0) SmearX0 = 4d0*atan(1d0)/sqrt(3d0)*kb_ev*abs(t_tetrakbt)/(2d0*rydberg_ev)   ! Ha
!   call gv_r(gw, 'wan_conv_1st',  wan_conv_1st) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_conv_end',  wan_conv_end) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_max_1st',   wan_max_1st) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_max_2nd',   wan_max_2nd) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_in_emin',   wan_in_emin) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_in_emax',   wan_in_emax) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_out_emin',  wan_out_emin) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_out_emax',  wan_out_emax) !removed 2026-10-02 (the Wannier functions / lmfham2)

    ! Newly added (m_readgwinput migration)
    call gv_r(gw, 'ecut_p',          ecut_p)
    call gv_r(gw, 'ecuts_p',         ecuts_p)
    call gv_r(gw, 'gauss_img',       gauss_img)

    ! Batch 2 reals
    call gv_r(gw, 'BZadiv',             BZadiv)
    call gv_r(gw, 'ene_sppola',         ene_sppola)
    call gv_r(gw, 'mixbeta',            mixbeta)
    call gv_r(gw, 'mixtj',              mixtj)
    call gv_r(gw, 'TFscreen',           TFscreen)
    call gv_r(gw, 'removed_r0c',        removed_r0c)
    call gv_r(gw, 'q0scale',            q0scale)
    call gv_r(gw, 'shift_majority',     shift_majority)
    call gv_r(gw, 'output_ddmat_atom',  output_ddmat_atom)
    call gv_r(gw, 'dRdIatRmax',         dRdIatRmax)
    call gv_r(gw, 'zmel_batch_gb',      zmel_batch_gb)
    call gv_r(gw, 'zmel_max_size',      zmel_batch_gb)  ! legacy alias
    clamp_zmel: block
      use m_gpu, only: use_gpu
      real(8) :: zmin
      zmin = merge(zmel_batch_gb_min_gpu, zmel_batch_gb_min, use_gpu)
      zmel_batch_gb = max(zmel_batch_gb, zmin)
    endblock clamp_zmel
    call gv_r(gw, 'magnon_delta',       magnon_delta)
    call gv_r(gw, 'magnon_delta_dos',   magnon_delta_dos)
    call gv_r(gw, 'magnon_HistBin_ratio', magnon_HistBin_ratio)
    call gv_r(gw, 'magnon_HistBin_dw',  magnon_HistBin_dw)
    call gv_r(gw, 'wmat_rcut1',         wmat_rcut1)
    call gv_r(gw, 'wmat_rcut2',         wmat_rcut2)
!   call gv_r(gw, 'wan_mix_1st',        wan_mix_1st) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_mix_2nd',        wan_mix_2nd) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_conv_2nd',       wan_conv_2nd) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_tbcut_rcut',     wan_tbcut_rcut) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(gw, 'wan_tbcut_heps',     wan_tbcut_heps) !removed 2026-10-02 (the Wannier functions / lmfham2)

    ! Batch 2 integers
    call gv_i(gw, 'ngcell',             ngcell)
    call gv_i(gw, 'nkeep_wfs',          nkeep_wfs)
    call gv_i(gw, 'nbcutlow_sig',       nbcutlow_sig)
!   call gv_i(gw, 'wan_nb_below',       wan_nb_below) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_nb_above',       wan_nb_above) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_out_bmin',       wan_out_bmin) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_out_bmax',       wan_out_bmax) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_in_bmin',        wan_in_bmin) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_in_bmax',        wan_in_bmax) !removed 2026-10-02 (the Wannier functions / lmfham2)
    call gv_i(gw, 'mixpriorit',         mixpriorit)
    call gv_i(gw, 'Q0Pchoice',          Q0Pchoice)
    call gv_r(gw, 'deltaq_scale',       DeltaQscale)
    call gv_i(gw, 'Verbose',            Verbose)
    call gv_i(gw, 'Q0P_Choice',         Q0P_Choice)
    call gv_i(gw, 'NormChk',            NormChk_int)  ! integer form (switch.f90)

    ! Boolean flags
    call gv_l(gw, 'KeepEigen',       KeepEigen)
    call gv_l(gw, 'KeepCMLO',        KeepCMLO)
    call gv_l(gw, 'KeepPPOVL',       KeepPPOVL)
    call gv_l(gw, 'NormChk',         NormChk)
    call gv_l(gw, 'unit_2pioa',      unit_2pioa)
    call gv_l(gw, 'AnyQ',            AnyQ)
    call gv_l(gw, 'QforEPSau',       QforEPSau)
    call gv_l(gw, 'QforEPSunita',    QforEPSunita)
    call gv_l(gw, 'QforEPSLIncLeft', QforEPSLIncLeft)
!   call gv_l(gw, 'wan_in_ewin',     wan_in_ewin) !removed 2026-10-02 (the Wannier functions / lmfham2)
    call gv_l(gw, 'KeepPositiveCou', KeepPositiveCou)

    ! Batch 2 logicals
    call gv_l(gw, 'EIBZmode',                  EIBZmode)
    call gv_l(gw, 'QforEPSIBZ',                QforEPSIBZ)
    call gv_l(gw, 'QforGWIBZ',                 QforGWIBZ)
    call gv_l(gw, 'TestOnlyQ0P',               TestOnlyQ0P)
    call gv_l(gw, 'TestNoQ0P',                 TestNoQ0P)
    call gv_l(gw, 'NoQ0P',                     NoQ0P)
    call gv_l(gw, 'KeepPpb',                   KeepPpb)
    call gv_l(gw, 'KeepWV',                    KeepWV)
    call gv_l(gw, 'KeepQG',                    KeepQG)
    call gv_l(gw, 'TimeReversal',              TimeReversal)
    call gv_l(gw, 'rmeshrefine',               rmeshrefine)
    call gv_l(gw, 'chi_RegQbz',                chi_RegQbz)
    call gv_l(gw, 'tetrahedron_matrix_linear', tetrahedron_matrix_linear)
    call gv_l(gw, 'magnon_w_onsite_dddd',      magnon_w_onsite_dddd)
    call gv_l(gw, 'magnon_negative_cut',       magnon_negative_cut)
    call gv_l(gw, 'allq0i',                    allq0i)
!   call gv_l(gw, 'wan_out_ewin',              wan_out_ewin) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(gw, 'wan_in_bwin',               wan_in_bwin) !removed 2026-10-02 (the Wannier functions / lmfham2)
    call gv_l(gw, 'wmat_static',               wmat_static)
    call gv_l(gw, 'wmat_all',                  wmat_all)
    call gv_l(gw, 'wmat_WSsuper',              wmat_WSsuper)
!   call gv_l(gw, 'wan_gauss_head',            wan_gauss_head) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(gw, 'wan_truncate',              wan_truncate) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(gw, 'wan_out_emax_auto',         wan_out_emax_auto) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(gw, 'wan_in_emax_auto',          wan_in_emax_auto) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(gw, 'wan_small_ham',             wan_small_ham) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_nsh1',                  wan_nsh1) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_i(gw, 'wan_nsh2',                  wan_nsh2) !removed 2026-10-02 (the Wannier functions / lmfham2)
    call gv_i(gw, 'mpi_worker_exch',           mpi_worker_exch)
    call gv_i(gw, 'mpi_worker_corr',           mpi_worker_corr)

    ! Integer vectors
    call gv_iv3(gw, 'n1n2n3',         n1n2n3)
    call gv_iv3(gw, 'n1n2n3eps',      n1n2n3eps)
    call gv_iv3(gw, 'multitet',       multitet)
    call gv_iv3(gw, 'n1n2n3dos',      n1n2n3dos)
    call gv_iv3(gw, 'GammaDivn1n2n3', GammaDivn1n2n3)

    ! Real 3-vectors
    call gv_rv3(gw, 'alpha_OffG_vec', alpha_OffG_vec)
    call gv_rv3(gw, 'wmat_rsite',     wmat_rsite)
  end subroutine load_gw_section

  !> All mlo_* scalar keys. Shared by the new [mlo] section and, for backward
  !  compatibility, by a [gw] that still carries them.
  subroutine load_mlo_keys(tbl, legacy)
    use m_mpi, only: master_mpi
    use m_lgunit, only: stdo
    type(toml_table), pointer, intent(in) :: tbl
    logical, intent(in) :: legacy
    integer :: st, idum
    if (legacy) then
       call get_value(tbl, 'mlo_method', idum, stat=st)   ! only asks whether the key exists
       if (st == 0 .and. master_mpi) write(stdo,'(a)') &
            ' m_GWinput: NOTE mlo_* keys found in [gw]; they belong in a [mlo] section now.'// &
            ' [gw] is still read, but please move them.'
    endif
    call gv_r(tbl, 'mlo_emax',      mlo_emax)        ! legacy reads as REAL
    call gv_i(tbl, 'mlo_method',    mlo_method)
    call gv_r(tbl, 'mlo_w',            mlo_w)
    call gv_r(tbl, 'mlo_delta',           mlo_delta)
    call gv_r(tbl, 'mlo_wfrz',           mlo_wfrz)
    call gv_r(tbl, 'mlo_down',           mlo_down)
    !call gv_r(tbl, 'mlo_low',            mlo_low)    ! hidden lower cut, see the comment at the declaration
    !call gv_r(tbl, 'mlo_wlow',           mlo_wlow)
    !call gv_r(tbl, 'mlo_pcut',           mlo_pcut)   ! hidden character cut, see the comment at the declaration
    !call gv_r(tbl, 'mlo_pw',             mlo_pw)
!   call gv_r(tbl, 'mlo_conv',           mlo_conv) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_mix',            mlo_mix) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_EUinner',        mlo_EUinner) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_CUouter',        mlo_CUouter) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_CUinner',        mlo_CUinner) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_WTinner',        mlo_WTinner) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_WTband',         mlo_WTband) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_WTseed',         mlo_WTseed) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_ELinner',        mlo_ELinner) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_ewid',           mlo_ewid) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_WTouter',        mlo_WTouter) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_CLhard',         mlo_CLhard) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_r(tbl, 'mlo_ELhard',         mlo_ELhard) !removed 2026-10-02 (the Wannier functions / lmfham2)
    !call gv_i(tbl, 'mlo_nskip',          mlo_nskip)   ! retired 2026-09-18: nskip is automatic (min over k of the non-model count)
!   call gv_i(tbl, 'mlo_maxit',          mlo_maxit) !removed 2026-10-02 (the Wannier functions / lmfham2)
!   call gv_l(tbl, 'mlo_EUinnerAUTOsp',         mlo_EUinnerAUTOsp) !removed 2026-10-02 (the Wannier functions / lmfham2)
  end subroutine load_mlo_keys

  !> [mlo] section: everything that defines the MLO model -- the mlo_* keys and
  !  the mlo_lm block (formerly Worb). It takes precedence over a Worb left in [blocks].
  subroutine load_mlo_section(mlo)
    type(toml_table), pointer, intent(in) :: mlo
    logical :: got
    call load_mlo_keys(mlo, legacy=.false.)
    ! mlo_lm: which lm channels (per atom) make the MLO model. This is the block
    ! that was called Worb under [blocks] until 2026-09-17; the old name is still
    ! accepted here and under [blocks].
    got = take_block(mlo, 'mlo_lm', block_Worb)
    if (.not. got) got = take_block(mlo, 'Worb', block_Worb)
    ! mlo_lm2: optional, same syntax, for the EH2 (second energy) radial set.
    got = take_block(mlo, 'mlo_lm2', block_Worb2)
    ! mlo_lm3: optional, same syntax, for the local-orbital (PZ) channels.
    got = take_block(mlo, 'mlo_lm3', block_Worb3)
  end subroutine load_mlo_section


  subroutine load_pb_section(pb)
    !> Read the cut-offs (pb_tolerance, pb_lcutmx) of [product_basis];
    !  the per-atom tables of the same section go through load_pb_tables.
    type(toml_table), pointer, intent(in) :: pb
    type(toml_array), pointer :: arr
    integer :: i, n

    call get_value(pb, 'pb_tolerance', arr, requested=.false.)
    if (associated(arr)) then
       n = len(arr)
       if (allocated(pb_tolerance)) deallocate(pb_tolerance)
       allocate(pb_tolerance(n))
       do i = 1, n
          call get_value(arr, i, pb_tolerance(i))
       enddo
    endif

    call get_value(pb, 'pb_lcutmx', arr, requested=.false.)
    if (associated(arr)) then
       n = len(arr)
       if (allocated(pb_lcutmx)) deallocate(pb_lcutmx)
       allocate(pb_lcutmx(n))
       do i = 1, n
          call get_value(arr, i, pb_lcutmx(i))
       enddo
    endif
  end subroutine load_pb_section


  !> Per-atom product-basis tables nlx / valence / core of a [product_basis] table.
  subroutine load_pb_tables(pb)
    type(toml_table), pointer, intent(in) :: pb
    call load_int_2darray(pb, 'nlx',     4, pb_nlx,     pb_n_nlx)
    call load_int_2darray(pb, 'valence', 5, pb_valence, pb_n_val)
    call load_int_2darray(pb, 'core',    7, pb_core,    pb_n_core)
  end subroutine load_pb_tables

  !> Pre-2026-09 layout: the tables in a separate PB.<sname>.toml. Not read;
  !  abort and point at ctrlg_absorb.py, which appends them to ctrlg.
  subroutine refuse_pb_file(filename)
    use m_ext, only: sname
    character(*), intent(in) :: filename
    logical :: exists
    inquire(file=filename, exist=exists)
    if (.not. exists) return
    call rx('m_GWinput: '//trim(filename)//' is retired and not read. Run  ctrlg_absorb.py '// &
         trim(sname)//'  to move its nlx/valence/core into [product_basis] of ctrlg.'// &
         trim(sname)//'.toml, then remove the file.')
  end subroutine refuse_pb_file


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


  !> Read a multi-line-string key into var only if the key is present
  !  (gv_c_alloc's intent(out) would otherwise deallocate a value loaded
  !  from another section). Returns .true. when something was read.
  logical function take_block(tbl, key, var)
    type(toml_table), pointer,     intent(in)    :: tbl
    character(*),                  intent(in)    :: key
    character(len=:), allocatable, intent(inout) :: var
    character(len=:), allocatable :: w
    call gv_c_alloc(tbl, key, w)
    take_block = allocated(w)
    if (take_block) call move_alloc(w, var)
  end function take_block

  subroutine load_blocks_section(blocks)
    use m_mpi, only: master_mpi
    use m_lgunit, only: stdo
    type(toml_table), pointer, intent(in) :: blocks
    logical :: old_eps, old_gw, old_worb
    call gv_c_alloc(blocks, 'QPNT',     block_QPNT)
    call gv_c_alloc(blocks, 'QforEPSL', block_QforEPSL)
    call gv_c_alloc(blocks, 'hrotr',    block_hrotr)
    ! Old homes of QforEPS / QforGW (now [gw]) and Worb (now [mlo] mlo_lm).
    ! Read them if the new home has not supplied them, and say so.
    old_eps  = .false.; old_gw = .false.; old_worb = .false.
    if (.not. allocated(block_QforEPS)) old_eps  = take_block(blocks, 'QforEPS', block_QforEPS)
    if (.not. allocated(block_QforGW))  old_gw   = take_block(blocks, 'QforGW',  block_QforGW)
    if (.not. allocated(block_Worb))    old_worb = take_block(blocks, 'Worb',    block_Worb)
    if ((old_eps .or. old_gw .or. old_worb) .and. master_mpi) write(stdo,'(a)') &
         ' m_GWinput: NOTE [blocks] still carries QforEPS/QforGW (now [gw]) and/or Worb'// &
         ' (now [mlo] mlo_lm). Still read, but please move them.'

    ! Parse raw text into structured arrays (caller-friendly).
    if (allocated(block_QforEPSL)) call parse_QforEPSL(block_QforEPSL, q_epsl, qend_epsl, idx_epsl, n_epsl)
    if (allocated(block_QPNT))     call parse_QPNT(block_QPNT)
  end subroutine load_blocks_section


  !> Open a scratch unit pre-loaded with the raw block text split into records.
  !  The unit is positioned at the start. Caller is responsible for closing.
  subroutine open_block_unit(text, unit)
    character(len=:), allocatable, intent(in)  :: text
    integer,                       intent(out) :: unit
    integer :: i, j, n
    open(newunit=unit, status='scratch', form='formatted')
    i = 1
    n = len(text)
    do while (i <= n)
       j = index(text(i:n), char(10))
       if (j == 0) then
          write(unit,'(a)') text(i:n)
          exit
       endif
       if (j == 1) then
          write(unit,'(a)') ''
       else
          write(unit,'(a)') text(i:i+j-2)
       endif
       i = i + j
    enddo
    rewind(unit)
  end subroutine open_block_unit


  !> Parse list of q-vectors (3 reals per line) from block text.
  !  Skips blank/comment lines (starting with '!' or '#').
  subroutine parse_qvec_list(text, qvec, nq)
    character(len=:), allocatable, intent(in)    :: text
    real(8), allocatable,          intent(inout) :: qvec(:,:)
    integer,                       intent(out)   :: nq
    integer :: u, ios, n, i
    real(8) :: q(3)
    character(256) :: line
    nq = 0
    if (.not. allocated(text)) return
    if (len(text) == 0) return
    call open_block_unit(text, u)
    n = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '!' .or. line(1:1) == '#') cycle
       read(line,*,iostat=ios) q
       if (ios == 0) n = n + 1
    enddo
    if (n == 0) then
       close(u); return
    endif
    if (allocated(qvec)) deallocate(qvec)
    allocate(qvec(3, n))
    rewind(u)
    i = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '!' .or. line(1:1) == '#') cycle
       read(line,*,iostat=ios) q
       if (ios /= 0) cycle
       i = i + 1
       qvec(:,i) = q
    enddo
    close(u)
    nq = n
  end subroutine parse_qvec_list


  !> Parse QforEPSL: q(3), qend(3), idx(int) per line.
  subroutine parse_QforEPSL(text, qv, qe, idx, nq)
    character(len=:), allocatable, intent(in)    :: text
    real(8), allocatable,          intent(inout) :: qv(:,:), qe(:,:)
    integer, allocatable,          intent(inout) :: idx(:)
    integer,                       intent(out)   :: nq
    integer :: u, ios, n, i, ii
    real(8) :: q(3), qq(3)
    character(256) :: line
    nq = 0
    if (.not. allocated(text)) return
    if (len(text) == 0) return
    call open_block_unit(text, u)
    n = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '!' .or. line(1:1) == '#') cycle
       read(line,*,iostat=ios) q, qq, ii
       if (ios == 0) n = n + 1
    enddo
    if (n == 0) then
       close(u); return
    endif
    if (allocated(qv))  deallocate(qv)
    if (allocated(qe))  deallocate(qe)
    if (allocated(idx)) deallocate(idx)
    allocate(qv(3,n), qe(3,n), idx(n))
    rewind(u)
    i = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '!' .or. line(1:1) == '#') cycle
       read(line,*,iostat=ios) q, qq, ii
       if (ios /= 0) cycle
       i = i + 1
       qv(:,i)  = q
       qe(:,i)  = qq
       idx(i)   = ii
    enddo
    close(u)
    nq = n
  end subroutine parse_QforEPSL


  !> Parse QPNT: multi-section format mirroring legacy block.
  !    line: allq spinonly       (2 ints)
  !    --- comment lines (start with *** or ! or other non-numeric) skipped ---
  !    line: nstates             (1 int)
  !    line: bands               (nstates ints)
  !    line: nq                  (1 int)
  !    nq lines: id qx qy qz    (1 int + 3 reals)
  subroutine parse_QPNT(text)
    character(len=:), allocatable, intent(in) :: text
    integer :: u, ios, ii, jj, idummy
    real(8) :: qq(3)
    character(256) :: line
    integer, allocatable :: bands(:)
    real(8), allocatable :: qs(:,:)
    if (.not. allocated(text)) return
    if (len(text) == 0) return
    call open_block_unit(text, u)
    ! 1) Find first non-comment line: allq spinonly
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) goto 999
       if (skip_line(line)) cycle
       read(line,*,iostat=ios) ii, jj
       if (ios == 0) then
          qpnt_allq = ii; qpnt_spinonly = jj
          exit
       endif
    enddo
    ! 2) Next data line: nstates
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) goto 999
       if (skip_line(line)) cycle
       read(line,*,iostat=ios) ii
       if (ios == 0) then
          qpnt_nstates = ii
          exit
       endif
    enddo
    ! 3) Next data line: bands(1:nstates)
    if (qpnt_nstates > 0) then
       allocate(bands(qpnt_nstates))
       do
          read(u,'(a)',iostat=ios) line
          if (ios /= 0) goto 999
          if (skip_line(line)) cycle
          read(line,*,iostat=ios) bands(1:qpnt_nstates)
          if (ios == 0) exit
       enddo
       if (allocated(qpnt_bands)) deallocate(qpnt_bands)
       call move_alloc(bands, qpnt_bands)
    endif
    ! 4) Next data line: nq
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) goto 999
       if (skip_line(line)) cycle
       read(line,*,iostat=ios) ii
       if (ios == 0) then
          qpnt_nq = ii
          exit
       endif
    enddo
    ! 5) qpnt_nq lines: idummy qx qy qz
    if (qpnt_nq > 0) then
       allocate(qs(3, qpnt_nq))
       ii = 0
       do
          if (ii >= qpnt_nq) exit
          read(u,'(a)',iostat=ios) line
          if (ios /= 0) exit
          if (skip_line(line)) cycle
          read(line,*,iostat=ios) idummy, qq
          if (ios /= 0) cycle
          ii = ii + 1
          qs(:, ii) = qq
       enddo
       if (allocated(qpnt_q)) deallocate(qpnt_q)
       call move_alloc(qs, qpnt_q)
       qpnt_nq = ii
    endif
999 close(u)
  end subroutine parse_QPNT


  !> Worb: each non-comment line is "iatom label lm1 lm2 ... lmN".
  subroutine parse_Worb(text)
    character(len=:), allocatable, intent(in) :: text
    call parse_worb_generic(text, n_worb, worb_iatom, worb_label, worb_lm, worb_nlm)
  end subroutine parse_Worb

  !> Same syntax as mlo_lm, for the EH2 radial set (mlo_lm2).
  subroutine parse_Worb2(text)
    character(len=:), allocatable, intent(in) :: text
    call parse_worb_generic(text, n_worb2, worb2_iatom, worb2_label, worb2_lm, worb2_nlm)
  end subroutine parse_Worb2

  !> Same syntax as mlo_lm, for the local-orbital (PZ) channels (mlo_lm3).
  subroutine parse_Worb3(text)
    character(len=:), allocatable, intent(in) :: text
    call parse_worb_generic(text, n_worb3, worb3_iatom, worb3_label, worb3_lm, worb3_nlm)
  end subroutine parse_Worb3

  subroutine parse_worb_generic(text, nrec, iat, lab_out, lm_out, nlm_out)
    character(len=:), allocatable, intent(in) :: text
    integer, intent(out) :: nrec
    integer, allocatable, intent(inout) :: iat(:), lm_out(:,:), nlm_out(:)
    character(8), allocatable, intent(inout) :: lab_out(:)
    integer :: u, ios, n, i, ib, lmtmp(16), nlm, k
    character(256) :: line
    character(8)   :: lab
    nrec = 0
    if (.not. allocated(text)) return
    if (len(text) == 0) return
    call open_block_unit(text, u)
    n = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (skip_line(line)) cycle
       n = n + 1
    enddo
    if (n == 0) then
       close(u); return
    endif
    if (allocated(iat))     deallocate(iat)
    if (allocated(lab_out)) deallocate(lab_out)
    if (allocated(lm_out))  deallocate(lm_out)
    if (allocated(nlm_out)) deallocate(nlm_out)
    allocate(iat(n), lab_out(n), lm_out(16,n), nlm_out(n))
    lm_out = -999
    rewind(u)
    i = 0
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (skip_line(line)) cycle
       lmtmp = -999
       read(line,*,iostat=ios) ib, lab, lmtmp(1:16)
       ! ios may be nonzero when fewer than 16 lm tokens present; that's OK
       if (ios /= 0 .and. ios > 0) then
          ! Try minimal: ib + lab + at least one lm
          read(line,*,iostat=ios) ib, lab, lmtmp(1)
          if (ios /= 0) cycle
       endif
       i = i + 1
       iat(i)     = ib
       lab_out(i) = lab
       lm_out(:,i)= lmtmp
       nlm = 0
       do k = 1, 16
          if (lmtmp(k) /= -999) nlm = k
       enddo
       nlm_out(i) = nlm
    enddo
    close(u)
    nrec = i
  end subroutine parse_worb_generic


  pure logical function skip_line(line)
    character(*), intent(in) :: line
    character(:), allocatable :: t
    skip_line = .false.
    t = adjustl(line)
    if (len_trim(t) == 0) then
       skip_line = .true.; return
    endif
    if (t(1:1) == '!' .or. t(1:1) == '#') then
       skip_line = .true.; return
    endif
    if (len(t) >= 3) then
       if (t(1:3) == '***') then
          skip_line = .true.; return
       endif
       if (t(1:3) == '---') then
          skip_line = .true.; return
       endif
    endif
  end function skip_line


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

  subroutine gv_rv3(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    real(8),                   intent(inout) :: var(3)
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

  !> Read a real that may be stored as scalar or as the first element of an array.
  !  Used for keys like 'emax_sigm' where TOML may have list[float] but legacy
  !  reads it as a scalar (first element).
  subroutine gv_rv_first(tbl, key, var)
    type(toml_table), pointer, intent(in)    :: tbl
    character(*),              intent(in)    :: key
    real(8),                   intent(inout) :: var
    type(toml_array), pointer :: arr
    integer :: stat
    call get_value(tbl, key, var, stat=stat)
    if (stat == 0) return
    call get_value(tbl, key, arr, requested=.false.)
    if (.not. associated(arr)) return
    if (len(arr) < 1) return
    call get_value(arr, 1, var, stat=stat)
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
