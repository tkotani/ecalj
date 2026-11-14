# Forces and atomic position relaxation {#forces-and-atomic-position-relaxation}
See the sample of the relaxation (LaGaO3).
We have to set the `[dyn]` section of `ctrlg.<sname>.toml` (only relaxation). 
We can set directions for relaxation with `relax = [0, 0, 1]` or so in `[[site]]`.
No cell optimizations.  

* CAUTION:
    For DYN mode, we treat the atomic position file `AtomPos.lagao3`. If it exists, we use atomic positions written in the file. 
    See `ecalj/SRC/subroutines/main_lmf.f90`. You find `call ReadAtomPos(irpos)` which override positions specified by `ctrlg.<sname>.toml`. 
    I have checked relaxation of the Perovskite LaGaO3 worked well. Since we do not use so much about the relaxation mode, be careful
    (I think some possible developments since atomic forces calculated well).
