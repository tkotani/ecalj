IO      SHOW=f HELP=F VERBOS=31,30 WKP=F IACTIV=0
% const bigbas=0 cupval=0
% const a=6.798 da=0 nit=12 vol=a^3/4 avw=(3/4/pi*vol)^(1/3) rmt=.87*avw nk=12 gmax=9 
STRUC   NBAS=1 NSPEC=1 NL=5
        ALAT={a} PLAT=  .0 .5 .5  .5  .0 .5  .5 .5  .0
        DALAT={da}
SITE    ATOM=A POS= 0 0 0
% const pwmode=0 pwemin=1 pwemax=3 oveps=0
HAM     NSPIN=1 REL=t XCFUN=2 
        FORCES=0 TOL=1e-6
        GMAX={gmax}
        FTMESH=10 10 10
        PWMODE={pwmode} PWEMIN={pwemin} PWEMAX={pwemax} OVEPS={oveps}
OPTIONS NSPIN=1 REL=t XCFUN=2 HF=f
BZ      NKABC={nk} BZJOB=1 W=.002 NPTS=1001 SAVDOS=t
%const metal=3
        EF0=0 DELEF=.1 TETRA=t DOS=0-1 0+.5 METAL={metal}
EWALD   AS=2.0 TOL=1D-12 ALAT0={a} NKRMX=600 NKDMX=600
ITER    MIX=A3 CONV=1e-5 CONVC=1e-5 NIT={nit}

SPEC    ATOM=A Z=29 R={rmt} IDMOD=0,0,0,1,1
        P=4.65,4.34,3.87,4.11
        PZ=5.5,5.5,4.5
%const lmx=3
        KMXA=4
        A=.025
        EREF=-3304.4345
        RSMG=.25*{rmt} RFOCA=0.4*{rmt} LFOCA=1 LMX=3 LMXA={lmx}
%const rsm1=2.5 rsmd1=1 ed1=-.01
%const rsm2=0 rsmd1x=0
        RSMH={rsm1},{rsm1},{rsmd1} EH=-.01,-.01,{ed1},-.01,-.01
        RSMH2={rsm2},0,{rsmd1x},{rsm2},0 EH2=-1,-1,-1,-.01,-.01
