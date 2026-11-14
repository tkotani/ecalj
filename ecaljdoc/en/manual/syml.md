# `syml.<sname>`

> English version of the page [manual/syml](../../manual/syml.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/syml.md of commit 694da29ec (source blob 661c6db1c957).

This is the file of the symmetry lines for band plots.
`syml.<sname>` can be made with `getsyml`.

## Getting the symmetry lines for band plots and a figure of the Brillouin zone
`syml.<sname>` is made from `ctrlg.<sname>.toml` (`ctrl.<sname>` in the old format).
Band plots need `syml.<sname>`.
After making it, it is also easy to edit `syml.<sname>` by hand for `job_band`.

There is getsyml.py in ecalj/GetSyml. It is
based on [seekpath](https://github.com/giovannipizzi/seekpath/) and [spglib](https://anaconda.org/conda-forge/spglib).

### Usage
When installing with `InstallAll.py`, a symlink `getsyml` to `ecalj/GetSyml/getsyml.py` is made.
Run it as follows.
```
getsyml nio
(or)
getsyml ctrls.nio
```
Together with the symmetry lines for the band plot, it shows the three-dimensional Brillouin zone.
Examples are in [BZsamples](https://ecalj.sakura.ne.jp/BZgetsyml/).
The symmetry lines are written to the file `syml.<sname>` for ecalj.
The numbers of divisions in `syml` are determined by a rough method, so correct them if necessary.

### Citation when you use getsyml
  When you publish a paper, in addition to citing ecalj, please cite the following (see [License](./license)).

   1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
     Band structure diagram paths based on crystallography,
     Comp. Mat. Sci. 128, 140 (2017)
     
   2.Please also cite spglib, a library essential to the implementation.
     A. Togo, K. Shinohara, I. Tanaka, Sci. Technol. Adv. Mater., Meth. 4, 2384822 (2024).

* Licenses of spglib and seekpath: [License](./license).

### (Notes for developers)
   a. Fix lmchk so that it writes the information needed to give reasonable values.
     For example, ndiv (the number of mesh points along a line).
   b. Numerical precision of the calculation.
     Since we read the output of lmchk, np.set_printoptions(precision=16) is meaningless.

 The symmetry-line file is the input for drawing energy bands along selected symmetry lines,
    or for making constant-energy surfaces such as Fermi surfaces.
    This file (its name is specified by a modifier of the command-line argument --band, described
    in the section "Command-line switches") can take several forms.

### Format of syml
    Bands are made along given symmetry lines.
    The following example is the input for the lines X→Γ and Γ→M of the simple cubic lattice.
       21  .5 0  0    0  0  0         X     Gamma
       21  0  0  0    .5 .5 0         Gamma  M
       0   0  0  0    0  0  0
    The first number is the number of points along each line.
    The next six are the q points of the start and the end.
    The last line is all zeros.
