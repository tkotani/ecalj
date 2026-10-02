# `syml.<sname>`

Symmetry line for the band plot.
We can generate `syml.<sname>` by `getsyml`.


## Get symmetry lines for band plot and Brillouwin zone plot.
`syml.<sname>` is generated from `ctrlg.<sname>.toml` (legacy: `ctrl.<sname>`).
`syml.<sname>` is needed for band plot.
After generated, you can easily edit `syml.<sname>` for `job_band`.

At ecalj/GetSyml, we have getsyml.py, which is based on the
seekpath at https://github.com/giovannipizzi/seekpath/
and spglib at https://anaconda.org/conda-forge/spglib

### Usage:
We make a softlink `getsyml` to `ecalj/GetSyml/getsyml.py` during the install by `InstallAll.py`.
Run 
```
getsyml nio
(or)
getsyml ctrls.nio
```
. This show 3D Brillouin zone together with symmetry lines for band plot.
See [BZsamples](https://ecalj.sakura.ne.jp/BZgetsyml/) here.
The symmetry lines are written into the `syml.<sname>` file for ecalj.
The number of divisions for `syml` is give by a crude algorism, so edit it if necessary.

### Needed citations for getsyml
  In addition to usual ecalj acknowledgement, following citations are required when you make a publication.

   1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
     Band structure diagram paths based on crystallography,
     Comp. Mat. Sci. 128, 140 (2017)
     
   2.You should also cite spglib that is an essential library used in the implementation.
     https://github.com/atztogo/spglib.git

* See Lincence.txt for spglib and seekpath.

### (memo for developer)
   a.Modify lmchk to write required information to supply reasonable.
     For example, ndiv (mesh size along lines).
   b.Numerical accuracy of calculations. 
     np.set_printoptions(precision=16) is not meaningful since we read output of lmchk

 symmetry-line file : input for plotting energy bands along
    selected symmetry lines or for generating constant-energy
    contours such as a Fermi surface.  This file (whose name is
    specified as a modifier with the command-line argument
    --band, described in the "Command-line switches" section) can
    take on of several forms.

### format of syml
    generate bands along specific symmetry lines.
    The following sample input illustrates input for lines
    X->Gamma and Gamma->M for the simple cubic lattice.
       21  .5 0  0    0  0  0         X     Gamma
       21  0  0  0    .5 .5 0         Gamma  M
       0   0  0  0    0  0  0
    The first number designates how many points along each line.
    The next six label the starting and ending q-points,
    respectively.  Note that the last line must contain zeros.

