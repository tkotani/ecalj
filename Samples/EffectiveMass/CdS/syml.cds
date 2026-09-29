### Lines from Gamma for the effective mass (mass mode): the columns after the labels are
### ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry): points ninit2..nend2 of a fine mesh of ndiv2 points on the line
### etolv = 0.3 Ry here (0.1 in ../GaAs): with SOC the top of the valence band lies slightly off Gamma, and lmf
### then counts the window of the bands to write from there, VBM - etolv to VBM + etolv; 0.3 Ry = 4.08 eV reaches
### the conduction band at Gamma (2026-09-30)
#ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry)
51    0 0 0      .5 .5  .5   Gamma  L       513     1    81    0.3       0.01
51    0 0 0      1.  0  0    Gamma  X       513     1    81    0.3       0.01
51    0 0 0      .75 .75 0   Gamma  K       513     1    81    0.3       0.01
