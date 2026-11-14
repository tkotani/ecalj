# How to make doping?

> English version of the page [manual/Memo_bgcharge](../../manual/Memo_bgcharge.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/Memo_bgcharge.md of the tree of c9f9a3003 (source blob 3d8e6a246f78).

<div class="clfix-old">

> ⚠️ **TOML migration (2026-05)** — The legacy `Category_Token`s in this
> page map onto `ctrlg.<sname>.toml` as: `SPEC_ATOM_Z` → `[[spec]].z`,
> `SPEC_ATOM_Q` → `[[spec]].q`, `BZ_ZBAK` → `[bz].zbak`. The conceptual
> recipe below is identical; only the file/key path changes. See
> [Input file](./toml_migration).

</div>

<div class="clfix">

The old keys below are, in `ctrlg.<sname>.toml`: `SPEC_ATOM_Z` → `[[spec]].z`, `SPEC_ATOM_Q` → `[[spec]].q`, `BZ_ZBAK` → `[bz].zbak`.

</div>

The 1st principle calculations are performed under the charge neutrality.

We have two method (or combined ) to treat fractional number of electrons in the primitive cell.

 1. You can use fractional numbers for SPEC_ATOM_Z, 
 2. Set valence charge by BZ_ZBAK.

-----
  Following sample is doping for Si, Z=14 and two atoms per cell.


##  Fractional Z
, such as Z=14.2 
Together with Z, you have to add SPEC_ATOM_Q= for atoms. 
For example, for Z=14.2, we may need to add Q=2,2.2 (Note Q= is the initial valence electron density for {s} {p} {d} {f}, successively.) 
Thus Q=2,2.2 implies adding +0.2 electron for p channel to keep charge neutrality of the 
atom when we calculate spherical atomic densities.
Run lmfa with 'grep conf'. 
It shows electron distribution of initial condition.
Then Q=0 (check of charge neutrality) is shown in the output of lmfa.
Pay attention that lmfa finishs normally. 

Run lmf to keep the console output to llmf. Then
```bash
 grep 'z=',llmf 
```
shows Z for atoms.


```bash
 grep -A5 'Charges:'
```
shows 
```
 Charges:  valence     8.40000   cores    20.00000   nucleii   -28.40000
    hom background     0.00000   deviation from neutrality:     -0.00000
```
Here we see   nuclei = - 2* 14.2  = -28.4  # given by Z

CAUTION: at the first iteration, Charges: shows such as
```
  Charges:  valence     8.00000   cores    20.00000   nucleii   -28.00000
  hom background     0.12300   deviation from neutrality: 0.12300
```
because of the initial condition is given by superposition of atomic densities.
Deviation is nonzero ==> But 'deviation from neutrality: 0' from the next iteration.


##  Back ground charge by ZBAK

Set BZ_ZBAK. If ZBAK=0.4, the back ground charge is negative as -0.4|e|, 
thus the number of electron becomes 0.4 smaller.
(cores + valence electron + ZBAK = |total nucleus charge|).
```bash
 grep -A5 'Charges:'
```
shows 
```
Charges:  valence     7.60000   cores    20.00000   nucleii   -28.00000
    hom background     0.40000   deviation from neutrality:     -0.00000
```
That is, 7.6+0.4+20.0  =  -28.0


## Combination of fractional Z and BZ_ZBAK
Let us try Z=14.2, Q=2,2.2, together with ZBAK=0.3
In this case, we have 
```
   Charges:  valence     8.10000   cores    20.00000   nucleii   -28.40000
    hom background     0.30000   deviation from neutrality:     -0.00000
```
This means that (2 atoms per cell).
  nuclei     = - 2* 14.2  = -28.4  # given by Z
  cores      =   2* 10    =  20    # given by the atom of integer part of Z.
  background =  0.3                # given by ZBAK
Thus we have 
  valence electron = 28.4 - 20- 0.3 = 8.1 

## Check in the console output {#background-charge-and-fractional-z}
You see console output, e.g, 
```
Charges: valence 7.80000 cores 20.00000 nucleii -28.00000
hom background .20000 deviation from neutrality: 0.00000
```
<span class="clfix">(The valence and core charges were swapped here (19.8 and 8.0); for Si with `zbak = 0.2` they are 7.8 and 20.0.)</span>
This is a case with `zbak = 0.2`.

NOTE: at the first iteration, Charges: shows such as
```
Charges: valence 8.00000 cores 20.00000 nucleii -28.00000
hom background 0.12300 deviation from neutrality: 0.12300
```
because of the initial condition by superposition of atoms. It shows
deviation seems nonzero. But charge should be conserved from the next iteration.
