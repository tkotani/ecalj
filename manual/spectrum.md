# Spectrum function of G.

> ⚠️ Binaries read `ctrlg.<sname>.toml` only. Use `Legacy2toml.py <sname>` to convert legacy `ctrl.<sname>` / `GWinput`. See [TOML migration](./toml_migration).

How to calculate $\langle {\bf q} n|\Sigma(\omega)|{\bf q} n\rangle$

We calculate the diagonal elements
$\langle \psi({\bf q},n)|\Sigma_{\rm c}(\omega) |\psi({\bf q},n)\rangle$
by the program `hsfp0` in the mode `--job=4`, after a one-shot GW calculation by `gw_lmfh`
(`gw_lmfh` leaves $W$ in the files `__WVR.*` and `__WVI.*`, which `hsfp0 --job=4` reads).
If you have `sigm.*`, it is read as in the case of `gwsc`.

We need to set $\bf q$ and band index $n$ for which we calcualte.
These are keys of `[gw]` in `ctrlg.<sname>.toml`, the same ones as for `gw_lmfh`.

(A) q points

The q points are those of `QforGW`, one q per line in the unit of `2pi/alat`:
```toml
[gw]
QforGW = """
 0.0 0.0 0.0
 0.5 0.0 0.0
 1.0 0.0 0.0
"""
```
Without `QforGW`, or with `QforGWIBZ = true`, you will have self-energy for all irreducible k points
of the mesh `n1n2n3`. This may be needed for A(omega).
Because of the shifted-mesh method, q points which are not on the mesh can be given in `QforGW`.

(B) Band index

The bands are chosen by their energies. `EMINforGW` and `EMAXforGW` (eV, relative to the Fermi energy)
give the range of the bands for which the self-energy is calculated:
```toml
[gw]
EMINforGW = -2.0
EMAXforGW =  3.0
```

(C) energy mesh

The mesh of the self-energy is $\omega - E_{\rm F} = 0.01\,{\rm Ry} \times i$ ($i=0,\pm1,\pm2,...$), up to twice the
highest energy of the mesh of $W$ along the real axis.

Note that imaginary part of Sigma is given as the comvolution of ImW(omega) and the pole of Green's function
(`[gw].t_sigmaw` (K) gives the width of the Fermi-Dirac smearing of the pole). Resolution for Im W (near omega=0) is set by `[gw].HistBin_dw`.

Run
```bash
gw_lmfh si -np 24
mpirun -np 24 hsfp0 si --job=4 > lsc_spec
```
Then we have SEComg.UP (DN) files. A line of the file is written in `main_hsfp0.f90` by
```
    write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,2f16.8,x,3f16.8)")
         iw,itq(i),ip,is, q(1:3,ip), wibz(ip), eqx(i,ip,is),
         (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip)
```
This means we use energy in eV.

**Table 1**. Columns of `SEComg.UP` (`SEComg.DN`)

| column | variable | meaning |
|---|---|---|
| 1 | `iw` | omega index |
| 2 | `itq(i)` | band index |
| 3 | `ip` | q point index (in the order of `QforGW`) |
| 4 | `is` | spin index |
| 5-7 | `q` | q vector (cartesian in 2pi/alat) |
| 8 | `wibz` | weight of the q point in the irreducible BZ |
| 9 | `eqx` | eigenvalue in eV, relative to the Fermi energy |
| 10 | `(omega(i,iw)-ef)*rydberg()` | omega relative to the Fermi energy |
| 11, 12 | `hartree*zsec(iw,i,ip)` | Self energy. real and imaginary part |

The data of a band at a q point make a block; the blocks are divided by two blank lines.
When you change `QforGW`, run `gw_lmfh` again (it makes the eigenfunctions at the q points of `QforGW`).

Cautions for q points of `QforGW` that are not on the mesh:

- Set `EMAXforGW`. At such q the eigenvalue solver returns fewer states than `nband` and pads the rest with 1d20;
  without an upper limit these padded states enter $\Sigma$ and `hsfp0` stops with `sxcf 222: |w-e| out of range`.
- When you change `EMINforGW` or `EMAXforGW`, redo the exchange part as well (`gw_lmfh` from the start). The exchange files
  (`SEXU`, `XCU`) and the correlation file (`SECU`) must have the same bands; otherwise `hqpe` stops at the end of `SEXU`.
- The standard output of `hgw` and `hx0fp0` has lines `epsWVR: iq iw omg eps(wFC) eps(woLFC)` for the small q around $\Gamma$:
  $\varepsilon(\omega)$ there shows the plasmon (Re $\varepsilon = 0$) without an extra run. A level that falls on a sharp plasmon structure
  is sensitive to the real-axis mesh; check it with a finer `HistBin_ratio` (e.g. 1.03 → 1.015).

* Example.
```
  plot 'SEComg.UP' u ($10):($11) w l,'' u ($10):($12) w l
```
  can give a plot for Re (Sigma_c(omega)) and Im(Sigma_c) (columns of Table 1).

## 4
 To get integrated spectrum function (DOS), we need to superpose all the spectrum function
 (All q points and all band index). Be careful about the degeneracy (multiplicity) for each q points.
  You have to build it from SEComg file.
  The weight of a q point is in the 8th column (Table 1); see also the lines with the keyword `Multiplicity`
  in the console output of qg4gw (lqg4gw).

  Anyway, consider about ``is it worth to do?''
  To confirm your result, use sum rule (sum of spectrum weight). And pay attention to the relation
  between real and imag parts (Hilbert transformation).
