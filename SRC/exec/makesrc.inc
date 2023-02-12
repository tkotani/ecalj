#NOTE: when you add your foobar.F, choose one of *_list.

ELIB_list=\
m_cmdpath.o\
diagcv2.o \
m_shortn3.o\
m_shortn3_plat.o\
m_smhankel.o\
m_orbl.o\
m_ftox.o\
getgv2.o\
m_gennlat.o\
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
writeband.o\
writepdos.o\
lmhelp.o\
util.o\
m_hamindex0.o \
m_hamindex.o \
m_struc_def.o \
m_gtv.o \
addrwf.o \
atomsr.o \
m_rhocor.o\
atwf.o \
besslr.o \
bzints.o \
bzmesh.o \
bzwts.o \
chgmsh.o \
config.o \
defpq.o \
defwsr.o \
delstp.o \
dostet.o \
dstrbp.o \
dvdos.o \
freeat.o \
fswgts.o \
getq.o \
gtbsl1.o \
gtpcor.o \
gvctof.o \
gvlst2.o \
m_hansr.o \
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
lattic.o \
ldau.o \
madmat.o \
maknos.o \
makrwf.o \
mkqp.o \
m_mksym.o \
m_mksym_util.o \
mtchae.o \
nwit.o \
phidx.o \
prdmts.o \
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
s2sph.o \
scg.o \
slinz.o \
splcls.o \
splwts.o \
strx00.o \
subzi.o \
sugcut.o \
suham.o \
suldau.o \
sylm.o \
symdmu.o \
symfor.o \
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
xcpbe.o \
m_ext.o\
amix.o \
broyj.o \
cmdopt.o \
cpudel.o \
sort.o \
fftz3.o \
finits.o \
fpiint.o \
gradfl.o \
gradzr.o \
hunti.o \
huntx.o \
iprint.o \
llnew.o  \
matm.o \
mklagw.o \
mklegw.o \
mpibc.o \
nintsp.o \
platsl.o \
poldvm.o \
polint.o \
rfalsi.o \
strings.o \
words.o \
zhev.o \
tc.o \
scglp1.o \
addrbl.o \
atqval.o \
augmat.o \
augmbl.o \
bndfp.o \
bstrux.o \
dfqkkl.o \
dfrce.o \
elocp.o \
esmsmves.o \
fradhd.o \
fradpk.o \
gaugm.o \
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
mkrout.o \
ovlpfa.o \
pnunew.o \
poinsp.o \
potpus.o \
prtev.o \
rdovfa.o \
rhogkl.o \
rhomom.o \
rlocbl.o \
rsibl.o \
sgvsym.o \
smhsbl.o \
smves.o \
sopert.o \
supot.o \
suylg.o \
symrho.o \
tbhsi.o \
totfrc.o \
vcdmel.o \
vesgcm.o \
vxtrap.o \
wrhomt.o \
polinta.o \
lmaux.o \
pwmat.o \
sugw.o \
grfmsh.o \
smvxcm.o \
m_vxcatom.o \
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
m_homoelectron.o \
m_readwan.o \
m_w0w0i.o \
getwemax.o \
wse.o \
bzints2.o \
wintzsg.o \
pplmat.o \
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
wfacx.o \
efsimplef.o \
extension.o \
rdpp.o \
genallcf_mod.o \
rgwinf_mod.o \
cputid.o \
mate.o \
reindx.o \
bz.o \
genqbz.o \
switches.o \
rwbzdata.o \
readeigen.o \
readqg.o \
iqindx.o  \
m_pkm4crpa.o \
m_anf.o \
qpe1.sc.o \
qpe1.o \
mkqg.o \
m_q0p.o \
q0irre.o \
basnfp.o \
excore.o \
mkjp.o \
strxq.o \
sxcf_fal2.sc.o \
sxcf_fal2.o \
wronkj.o \
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