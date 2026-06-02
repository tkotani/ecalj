"""Inline help text for ctrlG.<sname>.toml and PB.toml.

Used by ctrl2ctrltoml.py (Legacy2toml.py) and ctrlgenToml.py to embed
self-contained explanations next to each TOML section/key, mirroring the
documentation that ctrlgenM1.py used to emit into legacy ctrl files.

Two dictionaries:
  SECTION_HEADER[sec]  -> multi-line block printed BEFORE the [sec] line
  KEY_INLINE[sec][k]   -> short trailing '# ...' appended after the value

For top-level scalars (no section), use SECTION_HEADER['top'] / KEY_INLINE['top'].
Sections without entries simply emit no comments.
"""

SECTION_HEADER = {
    # ----------------------------------------------------------------- top-level
    'top': (
        "# === Top-level keys ===",
        "# symgrp: 'find' lets lmf detect maximum space-group symmetry from",
        "#   lattice + species. Set explicitly (e.g. 'r4z r3d r2x') to lower",
        "#   the symmetry by hand. 'lmchk foobar --pr60' lists generators.",
        "#   See https://ecalj.github.io/ecaljdoc/manual/lmf#symgrp",
    ),

    'io': (
        "# === IO ===",
        "# verbose: console verbosity (default 31). Larger -> more diagnostic",
        "#   output. Use 35-41 when investigating failures.",
    ),

    'struc': (
        "# === STRUC: lattice ===",
        "# alat in a.u.; plat is the primitive cell in units of alat (rows = a1,a2,a3).",
        "# nbas/nspec are filled automatically from the [[site]]/[[spec]] arrays.",
    ),

    'site': (
        "# === SITE: atomic positions ===",
        "# Each [[site]] table is one atom. xpos = fractional coords in plat;",
        "# pos = Cartesian coords in units of alat (use one or the other).",
    ),

    'spec': (
        "# === SPEC: per-species parameters ===",
        "# Touch mainly mmom (initial magnetic moment, NSPIN=2),",
        "# idu/uh/jh (LDA+U), and rsmh/eh (smoothed-Hankel basis).",
        "# kmxa = radial expansion order at tail sites (default 5).",
        "# Rule of thumb: kmxa > pwemax (in Ry). r = MT radius (a.u.).",
    ),

    # ----------------------------------------------------------------- BZ
    'bz': (
        "# === BZ: Brillouin-zone integration ===",
        "# nkabc = k-mesh divisions (1 entry -> isotropic, 3 entries -> per axis).",
        "# metal=3 + tetra=true is the safe default for metals.",
        "# For molecules in a cell use tetra=false, n=-1, w=0.001 (~157K).",
        "# fsmom (fixed-spin moment) constrains total magnetic moment when set.",
        "#   FSMOMMETHOD=0 solids, =1 discrete (LUMO-HOMO bias). Cannot mix with so=1.",
    ),

    # ----------------------------------------------------------------- ITER
    'iter': (
        "# === ITER: self-consistency ===",
        "# mix:  'A#' Anderson (stable), 'B#'/'B3' Broyden (faster, default).",
        "#       The trailing digit is history length kept (e.g. B3 = keep 3).",
        "# b:    mixing ratio (smaller -> stable but slower convergence).",
        "# conv: max |dE| between iterations.  convc: max d(rho_out - rho_in).",
        "# umix/tolu: LDA+U density-matrix mixing parameters.",
    ),

    # ----------------------------------------------------------------- HAM
    'ham': (
        "# === HAM: Hamiltonian / basis ===",
        "# nspin: 1 nonmagnetic, 2 spin-polarized (set spec.mmom for guess).",
        "# rel:   scalar-relativistic Schroedinger equation (default true).",
        "# xcfun: 1=VWN, 2=Barth-Hedin, 103=PBE-GGA.",
        "# gmax / ftmesh: real-space mesh for charge density. Use one.",
        "#   If gmax is set, ftmesh is auto-derived. Watch sugcut output: it",
        "#   reports the gmax actually required for the given tol.",
        "# pwmode: 0=MTO only (LMTO), 1=APW+MTO (PMT, |G|cut), 2=APW only,",
        "#         11=APW+MTO with |q+G|cut (recommended), 12=LAPW.",
        "#   GW driver (--jobgw) forces 11 internally.",
        "# pwemax (Ry): APW cutoff; large values may need larger spec.kmxa.",
        "# oveps: drop near-linear-dependent basis in (H-eO)z=0 diagonalisation.",
        "#   Default 1d-8.  Use 1d-6 if pwemax is large and basis is over-complete.",
        "# so: 0 no SO, 1 add L.S to H, 2 add Lz.Sz only (collinear).",
        "# scaledsigma: QSGW mixing factor.  1.0 = full QSGW, 0.8 = QSGW80.",
        "# readp / readpskipf / pnufix:",
        "#   readp=true reads P,PZ from atom-calculation results;",
        "#   pnufix=true freezes log-derivative B.C. of radial functions.",
        "# forces: 0 no force calc, 1 compute Hellmann-Feynman forces.",
    ),

    # ----------------------------------------------------------------- OPTIONS
    'options': (
        "# === OPTIONS ===",
        "# hf=true: Hartree-Fock-style non-self-consistent diag (debugging).",
    ),

    # ----------------------------------------------------------------- DYN
    'dyn': (
        "# === DYN: structural relaxation (LDA/GGA only, no cell relaxation) ===",
        "# mode: 0 skip, 4 conjugate-gradient, 5 Fletcher-Powell, 6 Broyden.",
        "# hess: true reads hessian from disk; false starts from H=I.",
        "# xtol/gtol: stop when displacement / force is below these (set one >0).",
        "# step:  initial step in units of alat.  nit: max relax steps.",
    ),

    # ----------------------------------------------------------------- STR
    'str': (
        "# === STR: real-space structure (Ewald, neighbour tables) ===",
        "# rmax: cutoff radius for structure constants.",
    ),

    'ewald': (
        "# === EWALD: Ewald summation ===",
        "# tol controls accuracy; default 1e-12 is fine.",
    ),

    # ----------------------------------------------------------------- GW
    'gw': (
        "# === GW: driver (used by lmf --jobgw=N and gwsc) ===",
        "# nband_chi0/emax_chi0: bands/energy window for chi0.",
        "# Touch only if you know what you are doing -- defaults are sensible.",
    ),

    'product_basis': (
        "# === PRODUCT_BASIS ===",
        "# pb_tolerance / pb_lcutmx : ctrlG.<sname>.toml (cut-off tunables).",
        "# nlx / valence / core     : PB.toml (sname-free, shared per-spec tables).",
    ),

    'blocks': (
        "# === BLOCKS: GW additional blocks (epsilon, etc.) ===",
        "# QforEPS / QforEPSL : q-point lists for dielectric function calc.",
    ),
}

KEY_INLINE = {
    'top': {
        'symgrp': "# 'find' = auto-detect from lattice (see header)",
    },
    'io': {
        'verbose': "# 31 default, 35 verbose, 41+ debug",
    },
    'bz': {
        'metal':   "# 0 insulator, 3 tetra+broadening (safe metal default)",
        'tetra':   "# true = tetrahedron BZ integration",
        'w':       "# energy broadening (Ry); 0.002 ~ 0.01",
        'npts':    "# division of DOS plot (larger = finer)",
        'savdos':  "# write DOS files",
        'nkabc':   "# k-mesh divisions",
    },
    'iter': {
        'mix':     "# B3 = Broyden hist=3 (default); A3 = Anderson hist=3",
        'b':       "# mixing ratio (smaller -> more stable)",
        'conv':    "# max dE between iterations (Ry)",
        'convc':   "# max d(rho_out - rho_in) (Ry)",
        'nit':     "# max self-consistency iterations",
        'umix':    "# LDA+U density-matrix mixing parameter",
    },
    'ham': {
        'nspin':       "# 1 nonmag, 2 spin-polarized",
        'rel':         "# scalar-relativistic",
        'xcfun':       "# 1=VWN, 2=Barth-Hedin, 103=PBE-GGA",
        'gmax':        "# real-space mesh G-cutoff (Ry); alt: ftmesh",
        'ftmesh':      "# real-space mesh divisions",
        'tol':         "# sugcut tolerance",
        'pwmode':      "# 0=MTO 1=PMT 2=APW 11=PMT(|q+G|, default) 12=LAPW",
        'pwemax':      "# APW cutoff (Ry)",
        'oveps':       "# overlap-eps for (H-eO)z=0; 1d-8 default",
        'so':          "# 0 none, 1 L.S, 2 Lz.Sz",
        'scaledsigma': "# QSGW mixing: 1.0 full, 0.8 = QSGW80",
        'readp':       "# read P,PZ from atom calc results",
        'pnufix':      "# fix B.C. of radial functions",
        'frzwf':       "# freeze wavefunctions",
        'forces':      "# 0 none, 1 compute forces",
    },
    'struc': {
        'alat': "# lattice constant (a.u.)",
        'plat': "# primitive vectors in units of alat",
    },
    'product_basis': {
        'pb_tolerance': "# drop linearly-dep products below this (default 1e-3)",
        'pb_lcutmx':    "# max l-cutoff per atom (4 for d-valence, 6 for f)",
        'tolerance':    "# (legacy alias, prefer pb_tolerance)",
        'lcutmx':       "# (legacy alias, prefer pb_lcutmx)",
    },
    'gw': {
        'n1n2n3':        "# BZ mesh for GW",
        'QpGcut_psi':    "# |q+G| cutoff for eigenfunctions (a.u. unless unit_2pioa=true)",
        'QpGcut_cou':    "# |q+G| cutoff for Coulomb / W",
        'unit_2pioa':    "# false: a.u.;  true: 2*pi/alat",
        'alpha_OffG':    "# offset-Gamma auxiliary function (dimensionless)",
        'emax_chi0':     "# (Ry) emax cutoff for chi0",
        'emax_sigm':     "# (Ry) Sigma exact below; extrapolated above",
        'nband_chi0':    "# explicit # of bands for chi0 (overrides emax_chi0)",
        'nband_sigm':    "# explicit # of bands for Sigma (overrides emax_sigm)",
        'iSigMode':      "# self-energy mode flag (3 = standard QSGW)",
        'niw':           "# # of imag-axis frequencies (try 6/10/12)",
        'HistBin_dw':    "# bin width on real-axis at omega=0",
        'HistBin_ratio': "# frhis(iw)=(dw/(r-1))*(exp((r-1)*(iw-1))-1)",
        'SmearX0':       "# (Ha) Gaussian smear for X0 (metals)",
        'GaussSmear':    "# Gaussian smearing for poles of G^LDA in hsfp0",
        'deltaw':        "# (a.u.) numerical-derivative mesh for Z factor",
        'delta':         "# (Ry) small imaginary shift (~1e-6)",
        'esmr':          "# (Ry) hsfp0 smearing; keep < band gap for insulators",
        'EMINforGW':     "# (eV, rel. EFermi) lower band cutoff for GW Sigma",
        'EMAXforGW':     "# (eV, rel. EFermi) upper band cutoff for GW Sigma",
        'QforGWIBZ':     "# true = use IBZ instead of <QforGW> list",
        'QforEPSau':     "# true = interpret <QforEPS> as a.u.",
        'wan_out_emin':  "# (eV, rel. EFermi) Wannier outer window low",
        'wan_out_emax':  "# (eV, rel. EFermi) Wannier outer window high",
        'wan_maxit_1st': "# 1st-stage Wannier max iterations",
        'wan_conv_1st':  "# 1st-stage Wannier convergence threshold",
        'wan_max_1st':   "# 1st-stage Wannier step size",
        'wan_maxit_2nd': "# 2nd-stage Wannier max iterations",
        'wan_max_2nd':   "# 2nd-stage Wannier step size",
        'wan_conv_end':  "# Wannier final convergence threshold",
    },
    'blocks': {
        'QforEPS':     "# q-list for eps(omega) head",
        'QforEPSL':    "# q-list for off-diagonal eps",
        'QforGW':      "# q-list for diagonal GW",
        'Worb':        "# orbital list for MLWF / MLO",
    },
}


def fmt_section_header(sec):
    """Return a multi-line string (with trailing '\n' on each line) of the
    section header comment block, or '' if no header is defined."""
    lines = SECTION_HEADER.get(sec)
    if not lines:
        return ''
    return '\n'.join(lines) + '\n'


def fmt_key_inline(sec, key):
    """Return inline trailing comment ('  # ...') or empty string."""
    return KEY_INLINE.get(sec, {}).get(key, '')


__all__ = ['SECTION_HEADER', 'KEY_INLINE',
           'fmt_section_header', 'fmt_key_inline']


import re as _re_app

def apply_toml_annotations(text):
    """Insert SECTION_HEADER blocks before each [section] / [[array]] line,
    and KEY_INLINE comments after each `key = value` line.
    Idempotent: section headers already containing "# === SEC" markers
    are skipped, and key lines that already have a trailing "#" comment
    are left alone.
    """
    lines = text.split("\n")
    have_header = set()
    for ln in lines:
        m = _re_app.match(r"^#\s*===\s*([A-Za-z_][A-Za-z0-9_]*)", ln)
        if m:
            have_header.add(m.group(1).lower())

    out = []
    seen_sec_in_pass = set()
    last_sec = "top"
    for ln in lines:
        m = _re_app.match(r"\s*\[\[?([A-Za-z_][A-Za-z0-9_]*)\]\]?", ln)
        if m:
            sec = m.group(1)
            if sec not in seen_sec_in_pass and sec not in have_header:
                h = fmt_section_header(sec)
                if h:
                    out.extend(h.rstrip("\n").split("\n"))
            seen_sec_in_pass.add(sec)
            last_sec = sec
            out.append(ln)
            continue
        m2 = _re_app.match(r"^(\s*)([A-Za-z0-9_/\"]+)(\s*=\s*)(.+)$", ln)
        if m2 and last_sec:
            indent, key, sep, val = m2.groups()
            if ("#" not in val
                    and "\"\"\"" not in val
                    and "'''" not in val):
                k = key.strip().strip("\"")
                c = fmt_key_inline(last_sec, k)
                if c and "\n" not in val:
                    ln = f"{indent}{key}{sep}{val}  {c}"
        out.append(ln)

    has_top_scalars = False
    for _ln in lines:
        if _re_app.match(r"\s*\[", _ln):
            break
        if _re_app.match(r"^\s*[A-Za-z_][A-Za-z0-9_]*\s*=", _ln):
            has_top_scalars = True
            break
    if "top" not in have_header and has_top_scalars:
        h_top = fmt_section_header("top")
        if h_top:
            idx = 0
            while idx < len(out):
                ln = out[idx]
                if _re_app.match(r"^#\s*===", ln):
                    break
                if not (ln.startswith("#") or ln.strip() == ""):
                    break
                idx += 1
            for ln in reversed(h_top.rstrip("\n").split("\n")):
                out.insert(idx, ln)
    return "\n".join(out)


__all__ = list(__all__) + ["apply_toml_annotations"]
