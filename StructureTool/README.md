# StructureTool: converters between ecalj structure files and POSCAR, and a VESTA launcher

Small helpers for the crystal-structure part of the input. They handle the structure only (lattice and sites), and are simple:
check the result (`lmchk`, or `viewvesta`) before using it.

| command | from | to |
|---|---|---|
| `vasp2ctrl.py POSCAR_foo [--alat=10.66]` | POSCAR (vasp5, Cartesian) | `ctrls.POSCAR_foo.vasp2ctrl` (structure-only ctrls; `ctrlgenToml.py` makes the `ctrlg.<sname>.toml` from a `ctrls.<sname>`) |
| `ctrl2vasp.py ctrlg.foo.toml` (also `ctrls.foo`, `ctrl.foo`) | ecalj input | `POSCAR_foo.vasp` (Cartesian) |
| `viewvesta.py <file>` | `ctrlg.foo.toml`, `ctrls.foo`, `ctrl.foo` or a POSCAR | opens VESTA (a ctrl file is converted with `ctrl2vasp` first) |
| `refineposcar.py --poscar IN --refine OUT [--primitive]` | POSCAR | the structure symmetrized by pymatgen's SpacegroupAnalyzer |
| `superlattice/generatepos.py` | two zinc-blende compounds and layer numbers | POSCAR of the (001)/(110) superlattice, strained by the elastic data of Van de Walle, PRB 39, 1871 (1989) |

- `--alat=` of `vasp2ctrl.py` sets ALAT in a.u.; give the lattice constant to get a clean ctrls.
- `viewvesta.py` finds VESTA by `which VESTA` (or `vesta`). VESTA: <http://jp-minerals.org/vesta/>. Bonds: Edit → Bonds in VESTA.
- `InstallAll.py` links `vasp2ctrl`, `ctrl2vasp` and `viewvesta` (without `.py`) into the bindir. The scripts import `convctrl.py`
  from this directory, so a copy of a script alone does not run.
- `superlattice/POSCAR.*` are the superlattices made by `generatepos.py` (A3B3, (001) and (110)).

POSCARs to try: `ecalj_auto/INPUT/gw1500/POSCARALL/`.
