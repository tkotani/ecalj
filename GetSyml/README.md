### Get symmetry lines, along which we make band plot,  and Brillouwin zone plot.  syml.* is generated from ctrl.*

In this directory, we have getsyml.py, which is based on the
seekpath at https://github.com/giovannipizzi/seekpath/
and spglib at https://anaconda.org/conda-forge/spglib

===========================
Requirement and Install:
With python3,
```
pip install --user q
spglib seekpath  plotly
```

===========================
Usage: 
Make softlink getsyml.py as getsyml. Then
```getsyml nio```
or
```getsyml ctrls.nio```
This show 3D Brillouin zone together with symmetry lines for band plot.
See [BZsamples](https://ecalj.sakura.ne.jp/BZgetsyml/) here.
The symmetry lines are written into the syml.* file for ecalj.
You can edit syml.* for bandplot by the job_band command.
The number of divisions for syml is give by simple algorism, so edit it if necessary.

===========================
Needed citations when we use.
  In addition to usual ecalj acknowledgement,
  following citations are required when you make a publication.

   1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
     Band structure diagram paths based on crystallography,
     Comp. Mat. Sci. 128, 140 (2017) 
   2.You should also cite spglib that is an essential library used in the implementation.
     https://github.com/atztogo/spglib.git

============
See Lincence.txt for spglib and seekpath.

============
TODO:
   a.Modify lmchk to write required information to supply reasonable.
     For example, ndiv (mesh size along lines).
   b.Numerical accuracy of calculations. 
     np.set_printoptions(precision=16) is not meaningful since we read output of lmchk
