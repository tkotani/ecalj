HEADER  sc C atom
IO      SHOW=F HELP=F VERBOS=30 40 50 60 WKP=F
%const ef0=0 nk=2 lmxa=3  nit=10 hf=f zbak=0
OPTIONS NSPIN=2 REL=t FRZ=0 NRMIX=2 TPAN=0 HF={hf} ESP=F XCN=0 LMH=0
        XCFUN=2 FORCES=12
SYMGRP  find
ITER    MIX=A2,b=.5 CONV=1e-5 CONVC=.0005 NIT={nit}
MIX     MODE=A2,b=.5 CONV=1e-5 CONVC=.0005
        NMIX=2
HAM     NSPIN=2  REL=t XCFUN=2
        FTMESH=25 25 25  TOL=1E-6 FRZ=f
        FORCES=12 CONV=1e-5 CONVC=.0005
STRUC   NBAS=1 NSPEC=1 
        ALAT=7.9370052598409968 PLAT=   1 1 0  1 0 1  0 1 1
FIT     WVS=1 1  NFIT=2 EFIT=-.5 -2
BZ      NKABC={nk} {nk} {nk}  BZJOB=f GETQP=f
        METAL=2 TETRA=t  SAVDOS=f DOS=ef0-1.5 ef0+0.5
        EF0={ef0} DELEF=.2 W=.004 NPTS=200 NEVMX=5 EFMAX=5
        ZBAK={zbak}
SITE    ATOM=C POS=  0   0   0
SPEC   ATOM=C Z=6 R=3 EREF=-74.9949 A=0.02  LMXA={lmxa} RMXA=0.6
       P=2.9,2.85,3.18,4.12 IDMOD=0,1
#      for the ionized case
#      Q=2,1  MMOM=0,1
       RSMH= 1.3 1.1 -1 -1  EH= -0.7,-0.2,0 0 MMOM=0,2
       RSMH2=0.8 0.8 EH2= -1.5,-1
