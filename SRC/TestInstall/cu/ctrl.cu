# This is an input file for Cu, with comments documenting the input.
# Lines beginning with '#' are comment lines.

# ... Version control
VERS    LMF-6 LMASA-6 LM:7 FP:7
IO      SHOW=f HELP=F VERBOS=31,30 WKP=F IACTIV=0
TESTLMF lmfa --no-iactiv cu -vnk=8 -vbigbas=f
        lmf  --no-iactiv cu -vnk=8 -vbigbas=f
        rm mixm.cu
        lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7
        lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7 --band:fn=syml
TEST3P  lmfa --no-iactiv cu -vnk=8 -vbigbas=t -vcupval=1
        lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vcupval=1 --rs=0
TESTGW  echo -1 | lmfgw cu -vnk=8 -vbigbas=t 
        echo  0 | lmfgw cu -vnk=8 -vbigbas=t 
        echo  1 | lmfgw cu -vnk=8 -vbigbas=t 
        echo cu | lmf2gw 
CLEAN   rm -f atm.cu dos.cu mixm.cu out.lmf.cu save.cu ctrl.cu log.cu moms.cu rst.cu wkp.cu
# ... Preprocessor variable declarations:
#   bigbas   bigbas=1 uses a larger basis
#   cupval   Puts cu p state in the valence as a local orbital
#            uses foca=0
% const bigbas=0 cupval=0

# ... Variable declarations used as input in other categories
#   a     the lattice constant, in a.u.
#   da    change in lattice constant: a+da is lattice constant
#   nit   maximum number of band passes for LDA self-consistency
#   nk    number k-points divisions for BZ integrations
#   dist  parameter defining lattice shear; see SHEAR below
#   bzj   0=> k-points include gamma; 1=> k-points offset from gamma; see BZJOB
#   vol   cell volume
#   avw   average WS radius
#   rmt   MT radius.  lmf can handle up to about
#                     10% overlaps with negligible loss in accuracy.
#   gmax  energy cutoff for specifying FT mesh
% const a=6.798
CONST   a={a} da=0 nit=12
        nk=12 dist=0 bzj=1
        vol=a^3/4 avw=(3/4/pi*vol)^(1/3) rmt=.87*avw
        gmax={cupval?12:9} nkgw=8
STRUC   NBAS=1 NSPEC=1 NL=5
        ALAT=a PLAT=  .0 .5 .5  .5  .0 .5  .5 .5  .0
        DALAT=da
# Use one the two following line with dist<>0 for a volume-conserving shear
# SHEAR=0 0 1 ... => tetragonal  SHEAR=1 1 1 ... => trigonal
        SHEAR=0 0 1 1+dist
#       SHEAR=1 1 1 1+dist
SITE    ATOM=A POS= 0 0 0
% const pwmode=0 pwemin=1 pwemax=3 oveps=0
HAM     NSPIN=1 REL=t XCFUN=2 
        FORCES=0 TOL=1e-6
        GMAX=gmax
        FTMESH=10 10 10
        PWMODE={pwmode} PWEMIN={pwemin} PWEMAX={pwemax} OVEPS={oveps}
GW      NKABC=nkgw GCUTB=2.7 GCUTX=2.2
% const hf=f
OPTIONS NSPIN=1 REL=t XCFUN=2 HF={hf}
BZ      NKABC=nk BZJOB=bzj W=.002 NPTS=1001 SAVDOS=t
# Because bigbas=t is really large, use a more cautious metal treatment
% ifdef bigbas
        METAL=3
% endif
        EF0=0 DELEF=.1 TETRA=t DOS=0-1 0+.5 METAL=2
% ifdef hf
        NEVMX=-1
% endif
EWALD   AS=2.0 TOL=1D-12 ALAT0=a NKRMX=600 NKDMX=600
# Because bigbas=t is really large, use a smaller mixing to help convergence
%ifdef bigbas&f
# for version 7
ITER    MIX=A1,b=.5,n=1;A0,b=.5,n=2 CONV=1e-5 CONVC=1e-5 NIT=nit
MIX     MODE=A1,b=.5,n=1;A0,b=.5,n=2 CONV=1e-5 CONVC=1e-5
%endif
# for version 7
ITER    MIX=A3 CONV=1e-5 CONVC=1e-5 NIT=nit
MIX     MODE=A3 CONV=1e-5 CONVC=1e-5
START   NIT=nit

#  ... Tokens for SPEC category
#  KMXA  defines the cutoff in the polynomial expansion of augmented basis functions.
#        (If not specified, a default is chosen.) KMXA=4 is a rather strict cutoff.
#  A=    parameter defining radial mesh for tabulation of augmented w.f. and density
#  EREF= reference energy, subtracted from the total energy.
#  RSMG= smoothing radius used in electrostatics
#  RFOCA=smoothing radius used in fitting core tails
#  LFOCA=1=> frozen core with tails expanded into the interstitial
#        2=> frozen core with xc pot from tails treated in perturbation theory
#  LMXA= l-cutoff for the basis function in augmentation spheres.
#
#  RSMH,EH, RSMH2, EH2 below define the basis set.
#  NB: the %const construct defines variables in a manner similar to the
#  CONST category above.  But variables defined with the %const or %var
#  are defined for the preprocessor stages, and are cleared once the
#  preprocessor is complete.  See doc/input-file-style.txt
SPEC    ATOM=A Z=29 R=rmt IDMOD=0,0,0,1,1
        P=4.65,4.34,3.87,4.11

%ifdef cupval&bigbas
        PZ=5.5,3.9,4.5
        LMXA=4 LMX=4
        LFOCA=0
%elseifd cupval
%stop 1 cupval must be used in conjunction with bigbas
%elseifd bigbas
        PZ=5.5,5.5,4.5
        LMXA=4 LMX=4
%endif
        KMXA=4

        A=.025
        EREF=-3304.4345
# The following line is not needed since these are the defaults.
        RSMG=.25*rmt RFOCA=0.4*rmt LFOCA=1 LMX=3 LMXA=3

% const rsm1=2.5 rsmd1=1 ed1=-.01
        RSMH={rsm1},{rsm1},{rsmd1} EH=-.01,-.01,{ed1},-.01,-.01
        PZ=5.5,5.5,4.5

%ifdef cupval&bigbas
% const rsmp=0.7 ep=-5 rsm2=1.3
        RSMH2={rsm2},{rsmp},{rsmd1},{rsm2},{rsm2*0} EH2=-1,{ep},-1,-.01,-.01
%elseifd bigbas
% const rsm2=1.3
        RSMH2={rsm2},0,{rsmd1},{rsm2},{rsm2*0} EH2=-1,-1,-1,-.01,-.01
%endif
