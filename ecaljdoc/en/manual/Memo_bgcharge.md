# How to dope

> English version of the page [manual/Memo_bgcharge](../../manual/Memo_bgcharge.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/Memo_bgcharge.md of commit e1dd775c3 (source blob af34ce36c9d3).

First-principles calculations are done under charge neutrality.
There are two ways (or their combination) to have a fractional number of electrons in the cell:

 1. Make the atomic number `[[spec]].z` fractional.
 2. Set the valence charge with the background charge `[bz].zbak`.

The examples below dope Si (Z=14, two atoms per cell).

## Fractional Z

Take Z=14.2. Together with Z, write `[[spec]].q` for the atom.
For Z=14.2 use e.g. `q = [2, 2.2]` (`q` is the initial number of valence electrons in s, p, d, f).
So `q = [2, 2.2]` adds 0.2 electron to p, keeping the atom neutral when the spherical atomic densities are computed.

Run `lmfa` with `grep conf` to see the initial configuration.
The output of `lmfa` shows Q=0 (the check of neutrality). Make sure that `lmfa` finished normally.

Run `lmf` and keep the console output in `llmf`.
```bash
 grep 'z=' llmf
```
shows Z of the atoms, and
```bash
 grep -A5 'Charges:' llmf
```
shows
```
 Charges:  valence     8.40000   cores    20.00000   nucleii   -28.40000
    hom background     0.00000   deviation from neutrality:     -0.00000
```
The nuclei are $-2 \times 14.2 = -28.4$ (from Z).

## Background charge zbak

Set `[bz].zbak`. With `zbak = 0.4` the background charge is −0.4|e|, and the number of electrons becomes 0.4 smaller
(core + valence + zbak = |total nuclear charge|).
```
Charges:  valence     7.60000   cores    20.00000   nucleii   -28.00000
    hom background     0.40000   deviation from neutrality:     -0.00000
```
That is, 7.6 + 0.4 + 20.0 = 28.0. With `zbak = 0.2` the valence is 7.8 and the cores 20.0.

## Fractional Z combined with zbak

With Z=14.2, `q = [2, 2.2]` and `zbak = 0.3`:
```
   Charges:  valence     8.10000   cores    20.00000   nucleii   -28.40000
    hom background     0.30000   deviation from neutrality:     -0.00000
```
With two atoms per cell,
- nuclei = $-2 \times 14.2 = -28.4$ (from Z)
- cores = $2 \times 10 = 20$ (from the atom of the integer part of Z)
- background = 0.3 (from zbak)

so the valence = 28.4 − 20 − 0.3 = 8.1.

## Caution at the first iteration {#background-charge-and-fractional-z}

At the first iteration `Charges:` shows, for example,
```
  Charges:  valence     8.00000   cores    20.00000   nucleii   -28.00000
  hom background     0.12300   deviation from neutrality: 0.12300
```
because the initial density is the superposition of atomic densities, so the deviation from neutrality looks nonzero.
From the next iteration it is `deviation from neutrality: 0` (the charge is conserved).
