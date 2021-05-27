#NOTE: when you add your foobar.F, choose one of *_list.

ELIB_list=\
m_bzintegration.o\
m_bandcal.o\
m_procar.o\
m_igv2x.o\
m_qplist.o\
m_clsmode.o\
m_mkpot.o\
m_lmfgw.o\
writefs.o\
writeham.o\
vbmmode.o\
getqmode.o\
phispinsym.o\
m_MPItk.o \
writedossawada.o\
m_clsprm.o\
xcfunc.o\
m_struc_func.o\
pqmix.o\
m_lmfinit.o\
latsym.o\
toksw_init.o\
writeband.o\
writepdos.o\
lmhelp.o\
show_programinfo.o\
util.o\
m_toksw.o \
m_psigd.o \
m_hamindex.o \
m_sstrnsname.o \
m_struc_def.o \
m_gtv.o \
addbas.o \
addrwf.o \
aiocor.o \
aiogen.o \
aiomom.o \
aiomp.o \
aiopar.o \
aiopot.o \
aiorme.o \
aiosop.o \
aiova.o \
asymop.o \
atomsr.o \
atwf.o \
avwsr.o \
baslst.o \
besslr.o \
bloch.o \
bzints.o \
bzmesh.o \
bzmio.o \
bzwts.o \
bzwtsf.o \
ccutup.o \
chgmsh.o \
chkdmu.o \
clabel.o \
clsset.o \
cmplat.o \
config.o \
contet.o \
cshear.o \
defpq.o \
defspc.o \
defwsr.o \
delstd.o \
delstp.o \
dfphi.o \
dlmn.o \
dosio.o \
dostet.o \
drr2.o \
dstrbp.o \
dvdos.o \
efrng2.o \
f100.o \
f102.o \
fermi.o \
fixef0.o \
fixpos.o \
fovlp.o \
freeat.o \
fswgts.o \
gengrp.o \
gensym.o \
getbzp.o \
getq.o \
getqp.o \
groupg.o \
grpfnd.o \
grpgen.o \
grpop.o \
gtbsl1.o \
gtbvec.o \
gtpcor.o \
gvctof.o \
gvlist.o \
gvlst2.o \
gvmtch.o \
hamfb3.o \
hanr.o \
hansmd.o \
hansmr.o \
hansr4.o \
hansr5.o \
hansr.o \
hansrz.o \
hft2rs.o \
hnsmft.o \
hsmq.o \
iax2fd.o \
iclbas.o \
intnos.o \
iobsmv.o \
ioeula.o \
iofa.o \
ioorbp.o \
iorbtm.o \
iors.o \
iosave.o \
latlim.o \
lattc.o \
lattdf.o \
lattic.o \
latvec.o \
ldau.o \
lx2vec.o \
madmat.o \
makdla.o \
makdos.o \
makidx.o \
maknos.o \
maknrs.o \
makrm0.o \
makrot.o \
makrwf.o \
mchan.o \
mkcond.o \
mkplat.o \
mkqp.o \
mkrtab.o \
mksym.o \
mshsiz.o \
mtchae.o \
mxmymz.o \
nghbor.o \
nnrl.o \
nwit.o \
orbl.o \
ovlchk.o \
ovmin.o \
pairc.o \
parmxp.o \
phidx.o \
phvsfp.o \
pkli.o \
pptrns.o \
prdmts.o \
prjpos.o \
prmx.o \
projql.o \
prrmsh.o \
psymop.o \
qdist.o \
radgkl.o \
radmsh.o \
radpkl.o \
rdeq.o \
relax.o \
rlxstp.o \
rmesh.o \
rnatm.o \
ropbes.o \
ropcsm.o \
roplat.o \
ropyln1.o \
ropyln.o \
rotevs.o \
roth.o \
rothrm.o \
rotpnt.o \
rotspu.o \
rotycs.o \
rseq.o \
rsmsym.o \
rsort.o \
rstr0.o \
s2sph.o \
scg.o \
sclwsr.o \
shear.o \
shftpp.o \
shoist.o \
shorbz.o \
shorps.o \
shortn.o \
shoshl.o \
siteid.o \
slinz.o \
slist.o \
soprm.o \
spcgrp.o \
splcls.o \
splwts.o \
stoner.o \
strx00.o \
subzi.o \
suclst.o \
sudmtu.o \
sugcut.o \
suham2.o \
suham.o \
suidx.o \
suldau.o \
sumlst.o \
sumsro.o \
suqlst.o \
sylm.o \
sylmnc.o \
symcry.o \
symdmu.o \
symfor.o \
symiax.o \
symlat.o \
symprj.o \
symq.o \
symstr.o \
symtbl.o \
tailsm.o \
tetirr.o \
trysop.o \
uspecb.o \
vecpkl.o \
veecomp.o \
vldau.o \
volsph.o \
vxc0sp.o \
xlgen.o \
xxxdif.o \
ylmrtg.o \
writebasis.o \
spackv_array_copy.o \
struc_strn.o \
switch.o \
rdsigm.o\
rdsigm2.o\
rotwv.o \
dpdbyl.o \
dpdftr.o \
pvsms2.o \
getcor.o \
radgra.o \
xcpbe.o \
vxcnls.o \
vxc0gc.o \
m_ext.o\
fsubs.o \
a2bin.o \
a2bina.o \
a2rotm.o \
a2vec.o \
amix.o \
arrprt.o \
awrite.o \
bin2a.o \
bitops.o \
brentm.o \
brmin.o \
broyj.o \
cdasin.o \
cdiv.o \
chcase.o \
chkhss.o \
cmdopt.o \
cmdstr.o \
cplxdm.o \
cpstr.o \
cpudel.o \
cross.o \
d1mach.o \
dbesj.o \
dcabsr.o \
derfc.o \
dfdump.o \
dgedi.o \
dgefa.o \
dgesl.o \
dinv33.o \
discop.o \
dmcpy.o \
dmpy.o \
dmpy22.o \
dmpyt.o \
dmsadd.o \
dmscop.o \
dot3.o \
dpcopy.o \
dpdump.o \
dpsadd.o \
dpscop.o \
dqinv.o \
dqinvb.o \
dsbak.o \
dschd.o \
dsev1.o \
dshell.o \
dsidi.o \
dsifa.o \
dsisl.o \
dsmpy.o \
dspbmm.o \
dsred.o \
dsredx.o \
dsum.o \
dsumdf.o \
dtribx.o \
dtridx.o \
dvheap.o \
dvshel.o \
epslon.o \
fftz3.o \
fgtcat.o \
finits.o \
fopna.o \
fpiint.o \
garg.o \
getcat.o \
gettime.o \
gradfl.o \
gradzr.o \
hrmint.o \
htribk.o \
htridi.o \
hunti.o \
huntx.o \
i1mach.o \
icopy.o \
idmxmn.o \
iinear.o \
imtql2.o \
imtqlv.o \
imxmn.o \
info.o \
ipdump.o \
iprint.o \
isanrg.o \
ishell.o \
isum.o \
isw.o \
ivheap.o \
ivshel.o \
iyamax.o \
ll.o \
macro.o \
matm.o \
minfit.o \
mkdlst.o \
mkilst.o \
mklagw.o \
mklegw.o \
mpibc.o \
mxint.o \
nintsp.o \
nxtarg.o \
nxtgrp.o \
parchv.o \
pars1v.o \
parsvc.o \
pertev.o \
platsl.o \
pmblk.o \
poldvm.o \
polint.o \
poseof.o \
ppfac.o \
praxis.o \
pretty.o \
pythag.o \
ran1.o \
rdfiln.o \
rdm.o \
rdtok.o \
rfalsi.o \
ropylg.o \
s8tor8.o \
setfac.o \
strings.o \
strip.o \
symvar.o \
symvec.o \
syscall.o \
tinvit.o \
tql2.o \
tqlrat.o \
tred2.o \
wordg.o \
words.o \
ymsadd.o \
ymscop.o \
ymtrns.o \
yscal.o \
ywrm.o \
yyaxpy.o \
yydot.o \
yydotc.o \
yygefa.o \
yygemm.o \
yygesl.o \
yyhifa.o \
yyhmpy.o \
yyqinv.o \
yyqnvb.o \
yysbmm.o \
yysbnv.o \
yysp2a.o \
yyxcpy.o \
z2herm.o \
zampy.o \
zerfc.o \
zhev.o \
zmcpy.o \
zmsadd.o \
zmscop.o \
zpmblk.o \
zqinv.o \
zqinvb.o \
ztoy.o \
ztoyy.o \
tc.o \
dpzero.o \
scglp1.o \
addbkg.o \
addrbl.o \
atqval.o \
augmat.o \
augmbl.o \
bndfp.o \
bstrux.o \
corprm.o \
dfaugm.o \
dfqkkl.o \
dfrce.o \
efldos.o \
elocp.o \
esmsmves.o \
fklbl.o \
flocbl.o \
fpchk.o \
fradhd.o \
fradpk.o \
fsmbl.o \
fsmbpw.o \
ftlxp.o \
gaugm.o \
gfigbl.o \
ggugbl.o \
ghibl.o \
ghigbl.o \
ghios.o \
gklbl.o \
gklft.o \
gklq.o \
gtbsl2.o \
gvsym.o \
hambl.o \
hgugbl.o \
hhibl.o \
hhigbl.o \
hhugbl.o \
hklbl.o \
hklft.o \
hklgbl.o \
hklos.o \
hsibl.o \
hsmbl.o \
hxpbl.o \
hxpgbl.o \
hxpos.o \
ioden.o \
jxpos.o \
lgstar.o \
lindsc.o \
lmfp.o \
locpot.o \
loctsh.o \
makusp.o \
makusq.o \
mixrho.o \
mkdmtu.o \
mkehkf.o \
mkekin.o \
mkewgt.o \
mkorbm.o \
mkpdos.o \
mkpot.o \
mkrout.o \
momusl.o \
mshdot.o \
mshint.o \
mshn3p.o \
mshvmt.o \
mullmf.o \
ovlocr.o \
ovlpfa.o \
paug1.o \
pkl2ro.o \
pnunew.o \
poinsp.o \
potpus.o \
praugm.o \
prrhat.o \
prtev.o \
rdovfa.o \
rhgcmp.o \
rhogkl.o \
rhomom.o \
rlocbl.o \
rsibl.o \
rsibl_ev.o \
setofl.o \
sgvsym.o \
smcorm.o \
smhsbl.o \
smves.o \
sopert.o \
splrho.o \
suphas.o \
supot.o \
surho.o \
suylg.o \
symrho.o \
tbhsi.o \
totfrc.o \
ugcomp.o \
vcdmel.o \
vesft.o \
vesgcm.o \
vlm2us.o \
vxtrap.o \
wrhomt.o \
ncutcorrect.o \
stonerpb.o \
stonerrsa.o \
polinta.o \
lmaux.o \
bndconn_v2.o \
gwcphi.o \
pwmat.o \
sugw.o \
smshft.o \
m_events.o \
grfmsh.o \
smvxcm.o \
vxcnsp.o \
vxcnlm.o \
w_psir.o \
gwinput_v2.o \
gintxx.o 


MPI_list= MPI_fpgw2.o


MLOC_list= \
maxloc0.o \
maxloc1.o \
maxloc2.o \
maxloc3.o \
maxloc4.o \
m_read_Worb.o\
m_DATA4GW.o\
m_QG.o\
m_FFT3D.o\
mymath.o\
myYlm.o\
wfn2dx_2.o\
m_MLWF.o\
wmatK.o \
wansub.o \
cubeformat.o\
xsfformat2.o\
expand_mesh.o\
wfn2dx_abc.o\
wfnrho_abc.o\
write_mesh2d.o\
RSMPI_mod.o \
efsimplef_RSMPI.o\
RSMPI_qkgroup_mod.o \
RSMPI_rotkindex_mod.o \
RSMPI_utils.o\
ent_MPI.o\
wmatK_mpi.o

GWLIB_list=\
hilbertmat.o\
dpsion5.o\
nwordr.o\
m_homoelectron.o \
m_readwan.o \
diagcv3.o \
m_w0w0i.o \
getwemax.o \
genallcf_dump.o \
wse.o \
bzints2.o \
wintzsg.o \
pplmat.o \
rs.o \
getbzdata1.o \
wcf.o \
tetwt5.o \
m_tetrakbt.o \
m_tetwt.o \
x0kf_v4h.o \
cinvrx.o \
zsvd.o \
m_zmel.o \
m_freq.o \
readpomat.o \
keyvalue.o \
rppovl.o \
nocctotg.o \
ppbafp.fal.o \
psi2b_v2.o \
psi2b_v3.o \
wfacx.o \
rydberg.o \
efsimplef.o \
extension.o \
rdpp.o \
genallcf_mod.o \
rgwinf_mod.o \
iopen.o \
cputid.o \
rw.o \
ext2.o \
mate.o \
mate1.o \
bib1.o \
index.o \
idxk.o \
maxnn.o \
reindx.o \
bz.o \
genqbz.o \
switches.o \
rwbzdata.o \
llnew.o  \
readeigen.o \
readqg.o \
iqindx.o  \
alloclist.o \
m_pkm4crpa.o \
m_anf.o \
qpe1.sc.o \
icompvv2.o \
iopenxx.o \
qpe1.o \
mopen.o \
checksymlon.o \
mkqg.o \
m_q0p.o \
q0irre.o \
basnfp.o \
excore.o \
mkjp.o \
strxq.o \
sxcf_fal2.sc.o \
sxcf_fal2.o \
diagcv2.o \
wronkj.o \
bessl.o \
lgen.o \
rpaq.o \
m_readeps.o \
m_readq0p.o \
m_readQGcou.o\
m_pbindex.o \
m_rotMPB.o \
m_iqindx_qtt.o \
rotwave.o \
rotcg.o \
rotgvec.o\
m_readVcoud.o\
m_itq.o\
nan.o\
m_readgwinput.o\
m_qbze.o\
m_eibz.o\
m_gaussianfilter.o\
m_llw.o\
m_readfreq_r.o\
m_readepswklm.o\
qpe1.qsgw.o\
mptauof.o
