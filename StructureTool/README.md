# StructureTool — converters between ctrl (or ctrls) and POSCAR, and a VESTA launcher

(2026-10-01: converted from `README.txt` to Markdown; the text below is as it was, with the spelling fixed. The `sample/` directory used in the examples was moved to
`trash/` — see `MD/past_log.md` §12; use your own POSCAR, or one from `ecalj_auto/INPUT/gw1500/POSCARALL/`. Errors in this text are listed in
`MD/TODOandQuestion.md`.)

NOTE: these are still very primitive, and not well written... (one of my students made them and t.kotani modified them a little).
So, do not believe them too much.

These commands below just treat the crystal structure part of the ctrl file and POSCAR (for VASP).
If you see something wrong or inconvenient, let t.kotani know it.

## Install

As a viewer, we can use VESTA (<http://jp-minerals.org/vesta/jp/download.html>).
This is required to invoke the viewer by the command `viewvesta.py` below. At the beginning of `viewvesta.py`, set the VESTA path.

## Usage

We now have three commands below. These commands use functions in `convert/`, so the relative path is important.

### `viewvesta.py`

A simple utility to invoke VESTA with `POSCAR_foo`, `ctrl.*` and `ctrls.*`.

You MUST set the VESTA path at the beginning of `viewvesta.py` as `VESTA='~/VESTA-x86_64/VESTA'`. You can replace VESTA with any
command of a viewer (VESTA is a little inconvenient because of its GUI).

```bash
cp -r sample sample.backup
cd sample
../viewvesta.py ctrl.nio
../viewvesta.py ctrl.cu2gase2     # this calls ctrl2vasp internally; thus POSCAR_cu2gase2 is generated
```

### `ctrl2vasp.py`

Convert POSCAR to ctrl file. For help, type this command without arguments.

```bash
mkdir TEST
cd TEST
cp ../sample/ctrl.nio .
../vasp2ctrl.py ctrl.nio
```

### `vasp2ctrl.py`

POSCAR file to ctrl (the current version is for Cartesian, but not so difficult if you like Direct).

```bash
mkdir TEST
cd TEST
cp ../sample/10-Opal.cif.vasp POSCAR_opal
../vasp2ctrl.py POSCAR_opal
../viewvesta.py POSCAR_cugase2
```

or

```bash
mkdir TEST
cd TEST
cp ../sample/ctrl.cu2gase2 .
../vasp2ctrl.py ctrl.cu2gase2
../viewvesta.py POSCAR_cugase2
```

To show bonds, use Edit-Bond in the VESTA console.

You can test it in your new directory with some files in the `sample/` directory.
Crystal structure files `sample/STRUC-CIF/*.cif` are taken from AIST (Osaka), <http://staff.aist.go.jp/nomura-k/japanese/itscgallary.htm>.
