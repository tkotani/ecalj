#NOTE: when you add your foobar.F, choose one of *_list.

ELIB_list=\
m_shortn3.o\
m_shortn3_plat.o\
m_smhankel.o\
derfc.o\
m_orbl.o\
m_ftox.o\
getgv2.o\
m_gennlat.o\
rydberg.o \
mkppovl2.o\
m_density.o\
m_ldau.o\
boltztrap.o\
m_bandcal.o\
m_procar.o\
m_igv2x.o\
m_qplist.o\
m_clsmode.o\
m_mkpot.o\
writefs.o\
writeham.o\
vbmmode.o\
getqmode.o\
phispinsym.o\
m_MPItk.o \
writedossawada.o\
xcfunc.o\
m_struc_func.o\
m_lmfinit.o\
toksw_init.o\
writeband.o\
writepdos.o\
lmhelp.o\
util.o\
m_toksw.o \
m_psigd.o \
m_hamindex0.o \
m_hamindex.o \
m_struc_def.o \
m_gtv.o \
addrwf.o \
atomsr.o \
atwf.o \
avwsr.o \
besslr.o \
bloch.o \
bzints.o \
bzmesh.o \
bzwts.o \
bzwtsf.o \
ccutup.o \
chgmsh.o \
clabel.o \
config.o \
cshear.o \
defpq.o \
defwsr.o \
delstd.o \
delstp.o \
dfphi.o \
dostet.o \
dstrbp.o \
dvdos.o \
fermi.o \
fixef0.o \
freeat.o \
fswgts.o \
getq.o \
gtbsl1.o \
gtbvec.o \
gtpcor.o \
gvctof.o \
gvlst2.o \
gvmtch.o \
hansmr.o \
hansr.o \
hansrz.o \
hnsmft.o \
hsmq.o \
iclbas.o \
intnos.o \
iofa.o \
iorbtm.o \
iors.o \
iors_old.o \
iosave.o \
latlim.o \
lattc.o \
lattic.o \
ldau.o \
madmat.o \
makdos.o \
maknos.o \
makrot.o \
makrwf.o \
mkqp.o \
mkrtab.o \
m_mksym.o \
m_mksym_util.o \
mshsiz.o \
mtchae.o \
mxmymz.o \
nwit.o \
phidx.o \
phvsfp.o \
pkli.o \
pptrns.o \
prdmts.o \
prjpos.o \
projql.o \
qdist.o \
radgkl.o \
radmsh.o \
radpkl.o \
relax.o \
rmesh.o \
rnatm.o \
ropbes.o \
roplat.o \
ropyln.o \
rseq.o \
rsort.o \
rstr0.o \
s2sph.o \
scg.o \
shear.o \
shftpp.o \
shoist.o \
siteid.o \
slinz.o \
soprm.o \
splcls.o \
splwts.o \
strx00.o \
subzi.o \
sugcut.o \
suham.o \
suldau.o \
sylm.o \
sylmnc.o \
symdmu.o \
symfor.o \
symprj.o \
symq.o \
symstr.o \
tailsm.o \
tetirr.o \
uspecb.o \
vecpkl.o \
xlgen.o \
ylmrtg.o \
switch.o \
rdsigm2.o\
pvsms2.o \
getcor.o \
radgra.o \
xcpbe.o \
m_ext.o\
fsubs.o \
a2bin.o \
a2bina.o \
a2vec.o \
amix.o \
awrite.o \
bitops.o \
brmin.o \
broyj.o \
cmdopt.o \
cpstr.o \
cpudel.o \
cross.o \
dinv33.o \
dmpy.o \
dmpy22.o \
dot3.o \
dpcopy.o \
dpdump.o \
dpscop.o \
dsmpy.o \
dspbmm.o \
dsred.o \
dsum.o \
dsumdf.o \
dvheap.o \
dvshel.o \
epslon.o \
fftz3.o \
finits.o \
fpiint.o \
garg.o \
gettime.o \
gradfl.o \
gradzr.o \
hunti.o \
huntx.o \
i1mach.o \
icopy.o \
imtql2.o \
imtqlv.o \
info.o \
iprint.o \
isanrg.o \
isw.o \
ivheap.o \
llnew.o  \
matm.o \
minfit.o \
mkdlst.o \
mkilst.o \
mklagw.o \
mklegw.o \
mpibc.o \
mxint.o \
nintsp.o \
platsl.o \
poldvm.o \
polint.o \
poseof.o \
ppfac.o \
pretty.o \
pythag.o \
ran1.o \
rdfiln.o \
rfalsi.o \
ropylg.o \
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
zerfc.o \
zhev.o \
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
dfqkkl.o \
dfrce.o \
efldos.o \
elocp.o \
esmsmves.o \
flocbl.o \
fradhd.o \
fradpk.o \
fsmbl.o \
fsmbpw.o \
ftlxp.o \
gaugm.o \
gklft.o \
gvsym.o \
hambl.o \
hklft.o \
hsibl.o \
lindsc.o \
lmfp.o \
locpot.o \
makusp.o \
makusq.o \
mixrho.o \
mkdmtu.o \
mkehkf.o \
mkekin.o \
mkewgt.o \
mkpdos.o \
mkrout.o \
momusl.o \
mshdot.o \
mshint.o \
mshn3p.o \
mshvmt.o \
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
setofl.o \
sgvsym.o \
smhsbl.o \
smves.o \
sopert.o \
splrho.o \
suphas.o \
supot.o \
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
gwcphi.o \
pwmat.o \
sugw.o \
smshft.o \
m_events.o \
grfmsh.o \
smvxcm.o \
vxcnsp.o \
w_psir.o \
gintxx.o 


MPI_list= MPI_fpgw2.o


MLOC_list= \
iopen.o \
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
lmf2gw.o\
zsecsym.o \
hilbertmat.o\
dpsion5.o\
nwordr.o\
m_homoelectron.o \
m_readwan.o \
diagcv3.o \
m_w0w0i.o \
getwemax.o \
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
zsvd.o \
m_zmel.o \
m_freq.o \
readpomat.o \
keyvalue.o \
rppovl.o \
nocctotg.o \
wfacx.o \
efsimplef.o \
extension.o \
rdpp.o \
genallcf_mod.o \
rgwinf_mod.o \
cputid.o \
rw.o \
mate.o \
bib1.o \
index.o \
maxnn.o \
reindx.o \
bz.o \
genqbz.o \
switches.o \
rwbzdata.o \
readeigen.o \
readqg.o \
iqindx.o  \
alloclist.o \
m_pkm4crpa.o \
m_anf.o \
qpe1.sc.o \
qpe1.o \
mopen.o \
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
rpaq.o \
m_readeps.o \
m_readq0p.o \
m_readQGcou.o\
m_pbindex.o \
m_rotMPB.o \
m_iqindx_qtt.o \
m_iqindx_wan.o \
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
eibzgen.o \
getw0.o

MLO_list=\
m_hamMTO.o