% macro xp(x1,x2,x3,x4) x1+2*x2+3*x3+4*x4
{xp(1,2,3,4)}
% trace 2
.SYMGRP R4Z MX I
.START   NIT=0 FREE=F BEGMOM=T CNTROL=t CNVG=1D-6 RDVES=t
% const yesvec=t
% const na=4 nb=12 fm=t dble=fm?((na+nb)%2+1):2 twoc=f dbxy=f fe=t v=t
% const ss=0 qss=0 nc=f theta=0 beta=90 sdyn=0 sdmod=1
% const nca=(na+1)/2 ncb=(nb+1)/2 nk1=60 nk2=8 vfe=0 dvfe=0
% const sc=f bulkm=~sc nspin=2 nit=1 xclas=0 ordalloy=f
% ifdef yesvec
% vec silly[12] 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
% vec silly(10:12) 110 111 112 113 114 115 116
% endif
% var dble=ordalloy|nc?2:dble
% var nk2=dble==2?nk2+nk2%4:nk2
% stop ~fm&(dble==1) : case fm={fm} and dble={dble} not allowed
% stop (na+nb)%2==1&dble==1 : na+nb={na+nb} and dble={dble} not allowed
% stop dbxy&na%2==0 dbxy=t requires that na be odd
% const asa=t lmxf=4 lmxb=221 ef0=0
% const ksdw=1-1/(nb+1)*(((fm&nb%2)|(~fm&(1-nb%2)))?0:1)
% const xv=1 ecnst=0

% getenv abc HOME
% echo HOME environment variable is: {abc}

% ifdef yesvec
try vectors ... {silly(9)} {silly(9+1)} {silly(9)*silly(10)}
% endif

  show that indirection works ...
% const essi=10 esge=20
% char nam=si
% const es=es{nam}
% echo es is {ES}, es{nam}
% char nam=ge
% var es=es{nam}
% echo is is {es}, es{nam}

  show that cchar works ...
% char eltc Ga elta N eltb In
% var i=4
% while i-=1 i>=0
% cchar elt i==1 {eltc} i==2 {elta} i==3 {eltb}
% echo for i={i}, elt is {elt}
% end

% cchar elt elt=='In' {eltc} i==2 {elta} i==3 {eltb}
% echo now elt is {elt}
% var i=2
% cchar elt elt=='In' {eltc} i==2 {elta} i==3 {eltb}
% echo now elt is {elt}

# --- Magnetic element use sbbcm=1 for m bcc, sqrt(2) for m fcc ---
% cconst fe Za=26 abccm=5.42 erefa=-2540.5681 sbccm=1
% ifdef fe
% char me fe
% endif
% ifdef co
% var Za=27 abccm=6.71/sqrt(2) erefa=-2781.6549 sbccm=sqrt(2)
% char me co
% endif
# ... Spacer element
% ifdef al
% var Zb=13 abccs=7.40/sqrt(2) erefb=lmh?-484.294:-483.546 sbccs=sqrt(2)
% char se al
% endif
% ifdef v
% var Zb=23 abccs=5.72 erefb=-1894.2951 sbccs=1
% char se v
% endif
% ifdef cr
% var Zb=24 abccs=5.42 erefb=-2097.3034 sbccs=1
% char se cr
% endif
% ifdef mn
% var Zb=25 abccs=5.42 erefb=-2312.6092 sbccs=1
% char se mn
% endif
% ifdef cu
% var Zb=29 abccs=6.82/sqrt(2) erefb=-3304.4346 sbccs=sqrt(2)
% char se cu
% endif
% ifdef ag
% var Zb=47 abccs=7.73/sqrt(2) erefb=-10621.9681 sbccs=sqrt(2)
% char se ag
% endif
% char strn "{me}/{se}"
% ifdef xv<>1
% char strn2 "{strn}-cr xv={xv}"  strn {strn2}
% endif
% ifdef ss<>0|qss<>0
% char strn2 "{strn} qss={qss} b={beta}"  strn {strn2}
% endif
% ifdef nc<>0
% char strn2 "{strn} theta={theta}" strn {strn2}
% endif
% echo n={na},{nb} fm={fm}({dble}) nk={nk1},{nk2} {strn}
% var Zs=(1-xv)*24+xv*Zb erefb=(1-xv)*-2097.3034+xv*erefb
% var abccs=(1-xv)*5.42+xv*abccs
% ifdef ordalloy
% var Zs=24
% endif

VERS    LMASA-4.1 LMFP-4
IO      SHOW=f HELP=F VERBOS=31 20 WKP=f IACTIV=f
% ifdef asa
MIX     MODE=B4,w=2,1,b=.2,n=4;A2,w=1,2,n=4 XIPMX={sc} BETV=.2
        AMODE=A,w=0,0,wa=1,fn=ma,n=2,b=1
MIX     MODE=B4,w=2,1,b=.2,n=4;A2,w=1,2,n=4 XIPMX={sc} BETV=.05
        AMODE=A,w=0,0,wa=1,fn=ma,n=2,b=1
# use this one for fine convergence of rho
MIX     MODE=B6,w=2,1,b=.2,k=4 XIPMX={sc} BETV=.05
#MIX     MODE=B6,w=2,1,b=.2,k=3 XIPMX={sc} BETV=.05
# use this one for simultaneous convergence of rho, EA
MIX     MODE=B6,w=2,1,wa=1,b=.2,k=4 XIPMX={sc} BETV=.05
        AMODE=A,w=0,0,wa=1,fn=ma,n=2,b=1
% else
MIX     AMIX=t BETA=0.1 BETSW=F CONV=.000001 CONVC=.001*5
        XIPMX=t BETAX=.003 20.3 .1 NIPMX=0 2 8
MIX     AMIX=t BETA=.1 BETSW=F CONV=.0001 CONVC=.001
        XIPMX=t BETAX=.01 .1 NIPMX=1 10 1 0 0
        XIPMX=t BETAX=.01 .01 NIPMX=1 7 3
        XIPMX=t BETAX=.01 .01 NIPMX=1 5 5
% endif
SOLVE   SMODE=0 TOL=.01 H=.01 EEPS=1d-5 VAR=hi
MASTER  a={abccs} nb={nb} fm={fm}
        JOB  1
CONST   pi4b3=4*pi/3
        nit=200 beta=.3 betv=.015 nv1=1 nv2=-6

        hb={sbccs/2} ha={sbccm/2}*{abccm/abccs}^3 hi=(ha+hb)/2
        hc=2*hi+({na}-1)*ha+hb*({nb}-1)
        nplane={na+nb} hi=(ha+hb)/2

        vola=a^3*ha                   volb=a^3*hb
        Ra=(vola/pi4b3)^(1/3)         Rb=(volb/pi4b3)^(1/3)
        delV=(2*hi-ha-hb)*a^3/4/pi4b3 hend=hc-hi delV*=
        Rai=(delV+Ra^3)^(1/3)         Rbi=(delV+Rb^3)^(1/3)

    bzj=11 rwa={asa?1:.86} rwb={asa?1:.88} nsp={nspin} mom=2.2*(nsp-1)
    zb0=ha*{na-1}+hi-hb

OPTIONS NSPIN=nsp REL=t TWOC={twoc} F ELIN=-.10 INVIT=t XCN=11 TPAN=t
% ifdef nc<>0 | qss<>0|ss<>0
        NONCOL=t {sdyn} SDMOD={sdmod} SDPRM=100
% endif
% ifdef qss<>0|ss<>0
        SS=0 0 {qss} pi/2*0
% endif
% if ~sc
        Q=BAND
% endif
EWALD   TOL=1d-12 NKDMX=1500 NKRMX=1500
% ifdef dir==110
SYMGRP  M(-1,1,0) R2Z I
% elseifd dbxy
SYMGRP  M(-1,1,0) R2Z I
% elseif beta<>90|sdyn
SYMGRP  MX MY R4Z
% else
SYMGRP  MX MY R4Z I
% endif
BZ      NKABC={nk1} {nk1} {nk2/dble} BZJOB=bzj NKABC2=3 BZJOB2=0
        TETRA=t METAL=t DOS=-.8 .8 EF0={ef0} DELEF=.001
% if ~sc
        NEVMX=-1
% endif
STR     RMAX={twoc?4.1:3.5} MODE={twoc?10:0}

# ---------------- Maps --------
% ifdef ordalloy
MAP  F  a*:  1: {fm?1:11}  a.{me}
       xa*:  1: {fm?1:01}  a.{me}
 *b*: 2: zc==24?1:(zc=={zb}?0:.5),Zc=={Zb}?1:(Zc==24?0:.5) a.cr a.{se}
% endif
%ifdef map==0 ... bulk moments
%  echo  ... map=0: bulk moments
%  ifdef xv==1
MAP  F *b*:  1: 1          a.{se}
%   else
MAP  F *b*:  2: {1-xv},{xv} a.cr a.{se}
%   endif
        a*:  1: {fm?1:11}  a.{me}
       xa*:  1: {fm?1:01}  a.{me}
%endif
% ifdef map==-1
% echo  ... map=-1: copy moments to files
MAP   T  *:  1: 1 .
% endif
% ifdef map==1 ... Spin-flip classes [ab]* when fm=0 (for eg bulkm=12)
% echo  ... map=1: Spin-flip classes [ab]* when fm=0
MAP   T  [ab]*:  1: {fm?1:11}  x\h\t
% endif
% ifdef map==2 ... contract cell by 1 nonmag side  (odd->even)
% stop t fix for alloy
MAP   T   b{flor(ncb)}:  3: 1,.5,-.5  \h\t b{1+flor(ncb)}\t a.{se}
         xb{flor(ncb)}:  3: 1,.5,-.5  \h\t xb{1+flor(ncb)}\t a.{se}
% endif
% ifdef map==3&nb%2 ... enlarge cell by 1 nonmag side (even->odd)
% stop t fix for alloy
MAP   T  b{flor(ncb)}:  1: 1 a.{se}
        xb{flor(ncb)}:  1: 1 a.{se}
% endif
% ifdef map==3 ... enlarge cell by 1 nonmag side (odd->even)
% stop t fix for alloy
MAP   T  b{flor(ncb)}:  2: .5,.5  \h\t a.{se}
        xb{flor(ncb)}:  2: .5,.5  \h\t a.{se}
% endif
% ifdef map==4
% echo  ... map=4: average :f with :a
MAP  F  [ab]*:  2: {fm?.5:10.5},{fm?10.5:.5} \h\t:f \h\t:a
       x[ab]*:  2: .5,.5                     \h\t:f \h\t:a
% endif
% ifdef map==4
% echo  ... map=5: sdw
MAP  F *b*: 2: 1,-xs*cos(pi/hb*{ksdw}*(z-zb0)) a.cr mm.cr
% endif
... enlarge cell by 2 nonmagnetic side (use 1st map only if nb even)
MAP     F xb{flor(ncb)}: a.cr
          xb{flor(ncb)-1}:  2: .5 .5  xb{flor(ncb)-1}\t a.cr
... enlarge cell by 1 mag side (odd->even)
MAP     F  a{flor(nca)}:  2: .5 {fm?0:10}.5  \h\t a.fe
          xa{flor(nca)}:  2: .5 .5    \h\t a.fe

# ---------------- CLASS ----------------
% var nclass=(flor(nca)+flor(ncb))*dble-(xclas?0:nb%dble)
# ... extra classes when lowered symmetry in random Euler angles
#% var nclass=sdyn&dble==1?flor(nca)*dble+nb:nclass
% ifdef dbxy
% var nclass=nclass+(flor(ncb)*dble-nb%dble)+(dble==2?2:0)
% endif dbxy
STRUC   NBAS=nplane*{dble*(dbxy?2:1)} NCLASS={nclass} NL=3
% ifdef dir==110
        ALAT=a PLAT= -.5 .5 .5   .5 -.5 .5  hc*{dble} hc*{dble} 0 TET=1
% elseifd dbxy
        ALAT=a PLAT= 1 1 0   1 -1 0   0 0 hc*{dble} TET=1
% endif
        ALAT=a PLAT= 1 0 0   0 1 0   0 0 hc*{dble} TET=1
CLASS
% ifdef ordalloy
% var Zs=(Zs==Zb)?24:Zb
% endif
% repeat id= 1:dble
% cchar x  id==1       X
% cchar g  ~fm&id==2   - t " "
% char i i
% repeat k= 1:nca
 ATOM={x}A{k} Z={Za} IDMOD=0 0 0 EREF={erefa} GROUP={k}
  R=Ra{i}*rwa LMXB={lmxb} LMXF={lmxf} NR=601 A=.02 GRP2={g}{k}
% if    dbxy&k==1
 ATOM={x}Z{k} Z={Zb} IDMOD=0 0 0 EREF={erefb} GROUP=-{k}
  R=Rb{i}*rwb LMXB={lmxb} LMXF={lmxf} NR=601 A=.02 GRP2={g}{nca+k}
% endif dbxy&k==1
% char i
% end
% char i i
% repeat k= 1:ncb
% ifdef ordalloy
%   var Zs=(Zs==Zb)?24:Zb
 ATOM={x}B{k} EREF={erefb} GROUP=-{k} Z={k==ncb&nb%2==1?(Zb+24)/2:Zs}
% else
 ATOM={x}B{k} Z={Zs} IDMOD=0 0 0 EREF={erefb} GROUP=-{k}
% endif
  R=Rb{i}*rwb LMXB={lmxb} LMXF={lmxf} NR=601 A=.02 GRP2={g}{nca+k}
% if    dbxy
 ATOM={x}B{k}2 Z={Zs} IDMOD=0 0 0 EREF=erefb GROUP=-{k}
  R=Rb{i}*rwb LMXB={lmxb} LMXF={lmxf} NR=601 A=.02 GRP2={g}{nca+k}
% endif dbxy
% char i
% end  loop over k
% ifdef ordalloy
%   var Zs=(Zs==Zb)?24:Zb
% endif
% end  loop over dble
# ... handle case low symmetry dynamics ...
% char i i
% repeat k= 1:ncb
 ATOM=B{k} Z={24} IDMOD=0 0 0 EREF={erefb} GROUP=-{100+k}
  R=Rb{i}*rwb LMXB={lmxb} LMXF={lmxf} NR=601 A=.02 GRP2={g}{100+nca+k}
% char i
% end  loop over k

# ---------------- SITE -----------------
SITE
# ... Significance of variables:   n: toggles between 0,1 as the planes
# ... are built up.  k2 marks the evolution of the classes which ascend
# ... to the center of the magnetic layers and descend away from the
# ... center.  x is string X for first na contigous planes of magnet;
# ... x is blank for the second (only present if cell doubled).
% var n=1
# ... Loop over cell doubling
% repeat id=1:dble
% cchar x  id==1  X
% repeat k= 1:na
% var n=(n+1)%2 k2= (k<nca)?k:2*nca-k
% ifdef dir==110
   ATOM={x}A{k2} POS=hc*{id-1}+ha*{k-1} hc*{id-1}+ha*{k-1} {n}/2
   ROT=z:{theta*(id-1)},y:pi/180*{beta} RELAX=0
% else
%   if ~dbxy
        ATOM={x}A{k2}  POS= {n}/2 {n}/2 hc*{id-1}+ha*{k-1}
     ROT=z:{theta*(id-1)},y:pi/180*{beta} RELAX={k2==flor(nca)?ecnst:1}
%   else
%     if k2<>1
        ATOM={x}A{k2} POS={n}/2   {n}/2 hc*{id-1}+ha*{k-1}
        ATOM={x}A{k2} POS={n}/2-1 {n}/2 hc*{id-1}+ha*{k-1}
%     endif k2<>1
%     if k==1&k2==1
        ATOM={x}A{k2} POS={n}/2   {n}/2 hc*{id-1}+ha*{k-1}
        ATOM={x}Z{k2} POS={n}/2-1 {n}/2 hc*{id-1}+ha*{k-1}
%     endif
%     if k<>1&k2==1
        ATOM={x}Z{k2} POS={n}/2   {n}/2 hc*{id-1}+ha*{k-1}
        ATOM={x}A{k2} POS={n}/2-1 {n}/2 hc*{id-1}+ha*{k-1}
%     endif k<>1&k2==1
%   endif dbxy
% endif dir
% end   loop over k=1:na

# ... Loop over planes of the spacer: vars n and k2 as in the magnetic
# ... planes.  Now char variable x is X for the ascending classes,
# ... and blank for the descending when dble is 2 or X when dble is 1.
% repeat k= 1:nb
% var n=(n+1)%2 k2=(k<ncb)?k:2*ncb-k
% cchar xx  dble==1|id==1&k2==k|id==2&(k2<>k|(k==ncb)&~xclas) X
% char x {xx}
# this line undoes site equivalence for spin dynamics
#% cchar x  dble==1&sdyn&k>=ncb "" t {xx}
% ifdef dir==110
   ATOM={x}B{k2}  POS=hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
                      hc*{id-1}+hb*{k-1}+ha*{na-1}+hi {n}/2
        ROT=z:pi*{k},y:pi/180*{beta}
% elseifd ~dbxy
       ATOM={x}B{k2}  POS= {n}/2 {n}/2 hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
#           ROT=z:pi*{k2},y:pi/180*{beta}
            ROT=z:pi*{k},y:pi/180*{beta}
% else
# ... Next line handles case dble is 1 & nb even: B{nca} equiv. B{nca}2
% cchar y  ~(k==ncb&nb%2==1&dble==1) 2
# ... Classes ascending
% if    k==k2&(id<>2|k<>ncb)
      ATOM={x}B{k2}    POS={n}/2 {n}/2 hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
           ROT=z:pi*{k},y:pi/180*{beta}
      ATOM={x}B{k2}{y} POS={n}/2-1 {n}/2 hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
           ROT=z:pi*{k},y:pi/180*{beta}
% endif k==k2
# ... Classes descending
% if    k>k2|(id==2&k==ncb)
      ATOM={x}B{k2}{y} POS={n}/2 {n}/2 hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
      ROT=z:pi*{k-1},y:pi/180*{beta}
      ATOM={x}B{k2}    POS={n}/2-1 {n}/2 hc*{id-1}+hb*{k-1}+ha*{na-1}+hi
      ROT=z:pi*{k-1},y:pi/180*{beta}
% endif k2<>1&k>k2
% endif dir==110,100,or dbxy
% end   loop over k
% end   loop over id

# ---------------- START (default) ------------
% ifdef nit<>0&bulkm<2
START   NIT=nit FREE=F BEGMOM=T CNTROL=f CNVG=1D-5
% endif
# ---------------- START (self-consistent solutions) ------------
% ifdef xv==1
%   include q.{me}{se}
% elseifd xv==0
%   include q.cr
% else
%   include q.{se}-cr
% endif
# --------------- START (bulk moments) --------------
%ifdef bulkm==1&xv==1
% echo reading START bulk moments (bulkm==1) ...
. nit<0 skips reading ctrl file for moments, nit=0 reads, no iter
START NIT={nit==0?0:1} BEGMOM=t CNTROL={nit<0?0:1} CNVG=1D-5 RDVES=T
% repeat id= 1:dble
% cchar x   id==1  X
% repeat k= 1:nca
 ATOM={x}A{k}
% if nspin==1
% include q0.{me}
% endif nspin==1
% if nspin==2&(id==1|fm)
% include q1.{me}
% endif
% if nspin==2&(id==2&~fm)
% include q2.{me}
% endif nspin==2
                   V=  vfe
% end  loop over k

% repeat k= 1:ncb
 ATOM={x}B{k}
% if nspin==1
% include q0.{se}
% endif nspin==1
% if nspin==2
% include q12.{se}
% endif nspin==2
                   V= 0
% end  loop over k
% end  loop over id
% endif bulkm==1
# --------------- START (potentials only bulkm==2) --------------
% ifdef bulkm==2
% echo reading START to shift V (bulkm==2) ...  dvfe={dvfe}
. nit<0 skips reading ctrl file for moments, nit=0 reads, no iter
START NIT={nit} BEGMOM=t CNTROL={nit<0?0:1} CNVG=1D-5 RDVES=T
% repeat id= 1:dble
% cchar x   id==1  X

#% repeat k= 1:nca
#   ATOM={x}A{k} V= vfe   dV = dvfe
#% end  loop over k

#% repeat k= 1:ncb
#   ATOM={x}B{k} V=0
#% end  loop over k

% repeat k= 1:nca
   ATOM={x}A{k} DV={dvfe}
% end  loop over k

% end  loop over id
% endif bulkm==2
# ---------------- START (in the absence of any other) ------------
START   NIT=0 FREE=F BEGMOM=T CNTROL=f CNVG=1D-5

  show that literal quoting of characters works:
we should see that \{this doesn't get expanded}


  show that multiple nesting works
% const es{1{2+{3+4}1}} = 2

  illustrate the syntax \{?~expr~strn1~strn2}

    MODE={?~k~B~C}{k+3}

% ifdef yesvec
  show that multiple nesting and indirection works
% vec  v_{es{nam}}_{2}[3]  -11 -22 -33
{v_{es{nam}}_{2}}
% endif

% char bassi " " nam si
% ifdef bas{nam}
%   echo test nesting of ifdef
% endif

  illustrate subexpressions in character variables
% char beta  4567654534
% char delta1 ->{beta(2,4)}<-
% char delta2 ->{beta('56',4)}<-
  delta1 should look like: '->567<-'
  delta1 is :              '{delta1}'
  delta2 should look like: '->6<-'
  delta2 is :              '{delta2}'
  show string substitution: |{beta(/45/' -a b c- '/)}| and |{beta(/45/' -d e f-'/,2,2)}| and |{beta(/45/' -hij-'/,1,2)}|
% show vars

