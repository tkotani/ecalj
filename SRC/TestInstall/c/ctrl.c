HEADER  sc C atom
VERS    LMASA-6 LMF-6 LM:7 FP:7
IO      SHOW=F HELP=F VERBOS=30 40 60 60 WKP=F
TESTLMF lmfa --no-iactiv c -vzbak=0
        lmf  --no-iactiv c -vzbak=0
TESTION lmfa --no-iactiv c -vzbak=1
        lmf  --no-iactiv c -vzbak=1
CLEAN   rm -f ctrl.c moms.c atm.c mixm.c rst.c save.c log.c hssn.c wkp.c bsmv.c syml.c bnds.c
OPTIONS NSPIN=2 REL=t FRZ=0 NRMIX=2 TPAN=0 HF=hf ESP=F XCN=0 LMH=0
        XCFUN=2 FORCES=12
% ifdef minx
        MSTAT: t,f,t,.0005,.002 NMOVE=10
% endif
SYMGRP  find
# for Version 7:
ITER    MIX=A2,b=.5 CONV=1e-5 CONVC=.0005 NIT=nit
MIX     MODE=A2,b=.5 CONV=1e-5 CONVC=.0005
        ELIND=-1 NMIX=2
# for Version 7:
HAM     NSPIN=2  REL=t XCFUN=2
        FTMESH=50 50 50  TOL=1E-6 FRZ=f
        FORCES=12 ELIND=-1 CONV=1e-5 CONVC=.0005
STRUC   NBAS=1 NSPEC=1 NL=lmxa+1
#       ALAT=10 PLAT=   1 0 0 0 1 0 0 0 1
        ALAT=10/2^(1/3) PLAT=   1 1 0  1 0 1  0 1 1
FIT     WVS=1 1  NFIT=2 EFIT=-.5 -2
BZ      NKABC=nk nk nk  BZJOB=f GETQP=f
        METAL=2 TETRA=t  SAVDOS=f DOS=ef0-1.5 ef0+0.5
       EF0=ef0 DELEF=.2 W=.004 NPTS=200 NEVMX=5 EFMAX=5
        ZBAK=zbak
SITE    ATOM=C POS=  0   0   0
SPEC   ATOM=C Z=6 R=3 EREF=-74.9949 A=0.02  LMXA=lmxa RMXA=0.6
       P=2.9,2.85,3.18,4.12 IDMOD=0,1
#      for the ionized case
#      Q=2,1  MMOM=0,1
       RSMH= 1.3 1.1 -1 -1  EH= -0.7,-0.2,0 0 MMOM=0,2
       RSMH2=0.8 0.8 EH2= -1.5,-1
CONST   ef0=0 nk=4
        rsma=.6 lmxa=3  nit=10
        hf=f zbak=0
START   NIT=nit
