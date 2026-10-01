# GetSyml: symmetry lines for band plots

`getsyml <sname>` writes `syml.<sname>`, the k-point path that `job_band` draws the bands along,
and shows the Brillouin zone with the path (plotly, in the browser).

- The cell is read from `ctrlg.<sname>.toml`: `getsyml` runs `lmchk <sname>` and reads `PlatQlat.chk` and `SiteInfo.lmchk`.
- The path is the one of seekpath (Hinuma et al., below) on the crystal found by spglib.
- The number of divisions of each line (the first column of `syml.<sname>`) is set by a simple rule; edit it if needed.

```
getsyml nio               # writes syml.nio and opens the Brillouin zone view
getsyml nio --nobzview    # writes syml.nio only (batch runs)
job_band nio -np 4        # bands along syml.nio
```

Requirements (python3): `pip install spglib seekpath plotly`.

## Citations

When you publish bands along these paths, cite in addition to ecalj:

1. Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography,
   Comp. Mat. Sci. 128, 140 (2017).
2. spglib, https://github.com/spglib/spglib

Licences: seekpath (MIT) `GetSyml/LICENSE.txt`, spglib (BSD-3) `SRC/external/spglib/COPYING`; the list is in the top `LICENCE`.

## Known limits

- `getsyml` takes the cell from the text output of `lmchk`, so its precision is that of the printed digits.
- `lmchk` could write what `getsyml` needs (e.g. the divisions) directly.
