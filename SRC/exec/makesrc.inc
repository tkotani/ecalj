#NOTE: when you add your foobar.F, choose one of *_list.

ELIB_list=\
main_hrcxq.o \
main_lmv7.o \
main_qg4gw.o \
main_lmfa.o \
main_lmfham1.o \
main_lmfham2.o \
main_hbasfp0.o \
main_gwinit.o \
main_hmagnon.o \
main_heftet.o \
main_hmlo.o \
main_hqpe.o \
main_hqpe.sc.o \
main_hqpe_qsgw.o \
main_hsfp0.o \
main_hsfp0.sc.o \
main_hsigmconv.o \
main_hvccfp0.o \
main_hx0fp0.o \
main_rdata4gw_v2.o \
main_lmchk.o \
main_hx0init.o \
tbhsi.o\
m_fatom.o\
m_HamPMT.o\
m_cmdpath.o\
diagcv2.o \
m_shortn3.o\
m_shortn3_plat.o\
m_smhankel.o\
m_orbl.o\
m_ftox.o\
getgv2.o\
m_gennlat.o\
m_density.o\
m_ldau.o\
m_ldau_util.o \
boltztrap.o\
m_bandcal.o\
m_procar.o\
m_igv2x.o\
m_qplist.o\
m_clsmode.o\
m_mkpot.o\
writeham.o\
vbmmode.o\
getqmode.o\
m_MPItk.o \
xcfunc.o\
m_lmfinit.o\
util.o\
m_hamindex0.o \
m_hamindex.o \
m_struc_def.o \
m_gtv.o \
m_rhocor.o\
m_atwf.o \
besslr.o \
m_bzints.o \
bzmesh.o \
m_bzwts.o \
chgmsh.o \
defpq.o \
delstp.o \
dstrbp.o \
freeat.o \
getq.o \
gtbsl1.o \
gtpcor.o \
m_hansr.o \
hsmq.o \
iofa.o \
iors.o \
iors_old.o \
lattic.o \
madmat.o \
mkqp.o \
m_mksym.o \
m_mksym_util.o \
mtchae.o \
radgkl.o \
radmsh.o \
radpkl.o \
relax.o \
rmesh.o \
ropbes.o \
roplat.o \
ropyln.o \
rseq.o \
scg.o \
slinz.o \
subzi.o \
sugcut.o \
suham.o \
sylm.o \
symq.o \
symstr.o \
tetirr.o \
uspecb.o \
xlgen.o \
ylmrtg.o \
switch.o \
rdsigm2.o\
xcpbe.o \
m_ext.o\
m_amix.o \
cmdopt.o \
cputid.o \
sort.o \
fftz3.o \
finits.o \
fpiint.o \
huntx.o \
iprint.o \
llnew.o  \
matm.o \
mpibc.o \
poldvm.o \
polint.o \
rfalsi.o \
zhev.o \
tc.o \
scglp1.o \
m_addrbl.o \
atqval.o \
m_augmat.o \
m_augmbl.o \
m_bndfp.o \
m_bndfp_util.o\
m_writeband.o\
m_bstrux.o \
dfrce.o \
elocp.o \
esmsmves.o \
fradhd.o \
fradpk.o \
gaugm.o \
hambl.o \
hklft.o \
hsibl.o \
m_lmfp.o \
locpot.o \
makusp.o \
makusq.o \
mixrho.o \
mkehkf.o \
mkrout.o \
pnunew.o \
poinsp.o \
potpus.o \
rdovfa.o \
rhomom.o \
rlocbl.o \
rsibl.o \
smves.o \
supot.o \
symrho.o \
totfrc.o \
vcdmel.o \
vesgcm.o \
lmaux.o \
pwmat.o \
sugw.o \
smvxcm.o \
m_vxcatom.o \
w_psir.o \
gintxx.o \
rotwave.o\
\
\
MPI_fpgw2.o\
\
\
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
wmatK_mpi.o \
\
\
m_readwan.o \
lmf2gw.o\
zsecsym.o \
hilbertmat.o\
dpsion5.o\
m_homoelectron.o \
m_w0w0i.o \
wse.o \
wintzsg.o \
pplmat.o \
getbzdata1.o \
tetwt5.o \
m_tetrakbt.o \
m_tetwt.o \
x0kf_v4h.o \
m_zmel.o \
m_freq.o \
keyvalue.o \
rppovl.o \
wfacx.o \
efsimplef.o \
extension.o \
rdpp.o \
genallcf_mod.o \
rgwinf_mod.o \
mate.o \
bz.o \
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
m_hvccfp0_util.o \
sxcf_fal2.sc.o \
sxcf_fal2_count.sc.o \
sxcf_fal2.o \
wronkj.o \
rpaq.o \
m_readeps.o \
m_readQGcou.o\
m_pbindex.o \
m_rotMPB.o \
m_iqindx_qtt.o \
m_iqindx_wan.o \
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
qpe1.qsgw.o\
eibzgen.o \
\
m_hamMTO.o