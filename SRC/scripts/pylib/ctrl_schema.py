#!/usr/bin/env python3
"""
ctrl.<sname>.toml のスキーマ定義 (single source of truth).

各エントリ: legacy_cattok -> (toml_section, toml_key, toml_type, legacy_subkey)

  toml_section : None for top-level scalars, 'site'/'spec' for [[site]]/[[spec]]
                 array-of-table, otherwise [section] table name (lowercase)
  toml_key     : key under that section (lowercase, snake_case)
  toml_type    : 'int','real','bool','str','int_vec','real_vec',
                 'int_vec3','real_vec3','int_mat3x3','real_mat3x3'
  legacy_subkey: the case-exact sub-token used by m_lmfinit's rval2 calls.
                 The loader emits "<SECTION>_<legacy_subkey>" so rval2 matches.

Coverage: live keys actually consumed by m_lmfinit::rval2. Anything not in
this dict is "dead emission" from old ctrl files (CHARGE/DOS/GW/START/OPTICS,
plus STR_DELTR/SIGMA/etc.) and is silently dropped by the converter with a
warning.
"""

# (section, key, type, legacy_subkey)
SCHEMA = {
    # -------- top-level (no section) --------
    'SYMGRP':       (None, 'symgrp',   'str', 'SYMGRP'),
    'SYMGRPAF':     (None, 'symgrpaf', 'str', 'SYMGRPAF'),

    # -------- [io] --------
    'IO_VERBOS':    ('io', 'verbos', 'int',     'VERBOS'),
    'IO_TIM':       ('io', 'tim',    'int_vec', 'TIM'),

    # -------- [struc] --------
    'STRUC_NSPEC':  ('struc', 'nspec', 'int',         'NSPEC'),
    'STRUC_NBAS':   ('struc', 'nbas',  'int',         'NBAS'),
    'STRUC_ALAT':   ('struc', 'alat',  'real',        'ALAT'),
    'STRUC_PLAT':   ('struc', 'plat',  'real_mat3x3', 'PLAT'),

    # -------- [[site]] (array of tables) --------
    'SITE_ATOM':    ('site', 'atom',  'str',       'ATOM'),
    'SITE_POS':     ('site', 'pos',   'real_vec3', 'POS'),
    'SITE_XPOS':    ('site', 'xpos',  'real_vec3', 'XPOS'),
    'SITE_RELAX':   ('site', 'relax', 'int_vec3',  'RELAX'),
    'SITE_AF':      ('site', 'af',    'int',       'AF'),

    # -------- [[spec]] (array of tables) --------
    'SPEC_ATOM':    ('spec', 'atom',   'str',      'ATOM'),
    'SPEC_Z':       ('spec', 'z',      'int',      'Z'),
    'SPEC_R':       ('spec', 'r',      'real',     'R'),
    'SPEC_R/W':     ('spec', 'r/w',    'real',     'R/W'),  # TOML quoted key
    'SPEC_R/A':     ('spec', 'r/a',    'real',     'R/A'),
    'SPEC_A':       ('spec', 'a',      'real',     'A'),
    'SPEC_NR':      ('spec', 'nr',     'int',      'NR'),
    'SPEC_RSMH':    ('spec', 'rsmh',   'real_vec', 'RSMH'),
    'SPEC_EH':      ('spec', 'eh',     'real_vec', 'EH'),
    'SPEC_RSMH2':   ('spec', 'rsmh2',  'real_vec', 'RSMH2'),
    'SPEC_EH2':     ('spec', 'eh2',    'real_vec', 'EH2'),
    'SPEC_LMX':     ('spec', 'lmx',    'int',      'LMX'),
    'SPEC_LMXA':    ('spec', 'lmxa',   'int',      'LMXA'),
    'SPEC_LMXL':    ('spec', 'lmxl',   'int',      'LMXL'),
    'SPEC_P':       ('spec', 'p',      'real_vec', 'P'),
    'SPEC_Q':       ('spec', 'q',      'real_vec', 'Q'),
    'SPEC_MMOM':    ('spec', 'mmom',   'real_vec', 'MMOM'),
    'SPEC_NMCORE':  ('spec', 'nmcore', 'int',      'NMCORE'),
    'SPEC_PZ':      ('spec', 'pz',     'real_vec', 'PZ'),
    'SPEC_LFOCA':   ('spec', 'lfoca',  'int',      'LFOCA'),
    'SPEC_KMXA':    ('spec', 'kmxa',   'int',      'KMXA'),
    'SPEC_RSMA':    ('spec', 'rsma',   'real',     'RSMA'),
    'SPEC_IDMOD':   ('spec', 'idmod',  'int_vec',  'IDMOD'),
    'SPEC_FRZWF':   ('spec', 'frzwf',  'bool',     'FRZWF'),
    'SPEC_IDU':     ('spec', 'idu',    'int_vec',  'IDU'),
    'SPEC_UH':      ('spec', 'uh',     'real_vec', 'UH'),
    'SPEC_JH':      ('spec', 'jh',     'real_vec', 'JH'),
    'SPEC_C-HOLE':  ('spec', 'c-hole', 'str',      'C-HOLE'),
    'SPEC_C-HQ':    ('spec', 'c-hq',   'real_vec', 'C-HQ'),
    'SPEC_EREF':    ('spec', 'eref',   'real',     'EREF'),

    # -------- [bz] --------
    'BZ_NKABC':       ('bz', 'nkabc',       'int_vec', 'NKABC'),
    'BZ_BZJOB':       ('bz', 'bzjob',       'int_vec', 'BZJOB'),
    'BZ_METAL':       ('bz', 'metal',       'int',     'METAL'),
    'BZ_TETRA':       ('bz', 'tetra',       'bool',    'TETRA'),
    'BZ_N':           ('bz', 'n',           'int',     'N'),
    'BZ_W':           ('bz', 'w',           'real',    'W'),
    'BZ_ZBAK':        ('bz', 'zbak',        'real',    'ZBAK'),
    'BZ_SAVDOS':      ('bz', 'savdos',      'bool',    'SAVDOS'),
    'BZ_NPTS':        ('bz', 'npts',        'int',     'NPTS'),
    'BZ_DOSMAX':      ('bz', 'dosmax',      'real',    'DOSMAX'),
    'BZ_EFMAX':       ('bz', 'efmax',       'real',    'EFMAX'),
    'BZ_FSMOM':       ('bz', 'fsmom',       'real',    'FSMOM'),
    'BZ_FSMOMMETHOD': ('bz', 'fsmommethod', 'int',     'FSMOMMETHOD'),

    # -------- [options] --------
    'OPTIONS_HF':   ('options', 'hf', 'bool', 'HF'),

    # -------- [ham] --------
    'HAM_NSPIN':       ('ham', 'nspin',        'int',       'NSPIN'),
    'HAM_REL':         ('ham', 'rel',          'bool',      'REL'),
    'HAM_SO':          ('ham', 'so',           'int',       'SO'),
    'HAM_SOCAXIS':     ('ham', 'socaxis',      'real_vec3', 'SOCAXIS'),
    'HAM_GMAX':        ('ham', 'gmax',         'real',      'GMAX'),
    'HAM_FTMESH':      ('ham', 'ftmesh',       'int_vec3',  'FTMESH'),
    'HAM_TOL':         ('ham', 'tol',          'real',      'TOL'),
    'HAM_FRZWF':       ('ham', 'frzwf',        'bool',      'FRZWF'),
    'HAM_XCFUN':       ('ham', 'xcfun',        'int',       'XCFUN'),
    'HAM_FORCES':      ('ham', 'forces',       'int',       'FORCES'),
    'HAM_RDSIG':       ('ham', 'rdsig',        'int',       'RDSIG'),
    'HAM_ScaledSigma': ('ham', 'scaledsigma',  'real',      'SCALEDSIGMA'),
    'HAM_EWALD':       ('ham', 'ewald',        'bool',      'EWALD'),
    'HAM_OVEPS':       ('ham', 'oveps',        'real',      'OVEPS'),
    'HAM_PWMODE':      ('ham', 'pwmode',       'int',       'PWMODE'),
    'HAM_PWEMAX':      ('ham', 'pwemax',       'real',      'PWEMAX'),
    'HAM_READP':       ('ham', 'readp',        'bool',      'READP'),
    'HAM_PHISPINSYM':  ('ham', 'phispinsym',   'bool',      'PHISPINSYM'),
    'HAM_PNUFIX':      ('ham', 'pnufix',       'bool',      'PNUFIX'),

    # -------- [iter] --------
    # NOTE: ITER_b/wc/w/k were lower-case in legacy. We standardize the
    # legacy_subkey to upper-case here; m_lmfinit.f90 must be updated to call
    # rval2 with the upper-case form (5 lines).
    'ITER_NIT':    ('iter', 'nit',   'int',      'NIT'),
    'ITER_NRMIX':  ('iter', 'nrmix', 'int',      'NRMIX'),
    'ITER_MIX':    ('iter', 'mix',   'str',      'MIX'),
    'ITER_CONV':   ('iter', 'conv',  'real',     'CONV'),
    'ITER_CONVC':  ('iter', 'convc', 'real',     'CONVC'),
    'ITER_UMIX':   ('iter', 'umix',  'real',     'UMIX'),
    'ITER_TOLU':   ('iter', 'tolu',  'real',     'TOLU'),
    'ITER_b':      ('iter', 'b',     'real',     'B'),
    'ITER_wc':     ('iter', 'wc',    'real',     'WC'),
    'ITER_w':      ('iter', 'w',     'real_vec', 'W'),
    'ITER_k':      ('iter', 'k',     'int',      'K'),

    # -------- [dyn] --------
    'DYN_MODE':   ('dyn', 'mode',  'int',  'MODE'),
    'DYN_NIT':    ('dyn', 'nit',   'int',  'NIT'),
    'DYN_HESS':   ('dyn', 'hess',  'bool', 'HESS'),
    'DYN_XTOL':   ('dyn', 'xtol',  'real', 'XTOL'),
    'DYN_GTOL':   ('dyn', 'gtol',  'real', 'GTOL'),
    'DYN_STEP':   ('dyn', 'step',  'real', 'STEP'),
    'DYN_NKILL':  ('dyn', 'nkill', 'int',  'NKILL'),

    # -------- [str] --------
    'STR_RMAX':   ('str', 'rmax',  'real', 'RMAX'),
    'STR_RMAXS':  ('str', 'rmaxs', 'real', 'RMAXS'),
    'STR_MXNBR':  ('str', 'mxnbr', 'real', 'MXNBR'),

    # -------- [ewald] --------
    'EWALD_TOL':  ('ewald', 'tol', 'real', 'TOL'),
}

# Section enumeration order (for emit ordering)
SECTION_ORDER = [
    None,        # top-level scalars first (symgrp, symgrpaf)
    'io',
    'struc',
    'site',      # [[site]]
    'spec',      # [[spec]]
    'bz',
    'iter',
    'ham',
    'options',
    'dyn',
    'str',
    'ewald',
]

# Sections that are arrays-of-tables ([[name]])
ARRAY_SECTIONS = {'site', 'spec'}

# Inverse lookup for the loader: by section, what legacy keys to emit
def reverse_lookup_by_section(section):
    """Return dict toml_key -> (legacy_subkey, type) for the section."""
    out = {}
    for legacy_full, (sec, k, typ, legacy_sub) in SCHEMA.items():
        if sec == section:
            out[k] = (legacy_sub, typ)
    return out
