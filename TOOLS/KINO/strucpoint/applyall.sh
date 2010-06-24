#/bin/bash
dir=/home/kino/kit/GW/7K/ecalj_2010_0618/lm7K

#gawk -f struc_delunused.awk m_struc_def.F |  python struc_makeinit.py > m_struc_def.F.changed
#\cp -f m_struc_def.F.changed $dir/subs/m_struc_def.F
#exit


gawk -f struc_change2.awk m_struc_def.F | python struc_makeinit.py > m_struc_def.F.changed
\cp -f m_struc_def.F.changed $dir/subs/m_struc_def.F

exit

#oipq
for name in fp/bndfp.F subs/asados.F subs/mkqp.F subs/rdsigm.F subs/rdsigm2.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# oics
for name in fp/dfrce.F fp/lmaux.F subs/asados.F subs/clsprm.F subs/mksym.F subs/rdctrl2.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# ohave
for name in fp/lmaux.F subs/clsprm.F subs/rdctrl2.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# ov0 ov1
for name in fp/locpot.F fp/mkrout.F fp/rdovfa.F fp/rsedit.F subs/iors.F fp/elocp.F fp/pnunew.F fp/sugw.F fp/vcdmel.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit




#osmpot
for name in fp/bndfp.F fp/lmfp.F fp/rsedit.F fp/supot.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#orhrmx
for name in fp/lmaux.F subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#oqt
for name in fp/lmaux.F subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#oqpp
for name in subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#oqnu
for name in fp/lmaux.F subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#oqc
for name in fp/lmaux.F subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#opp
for name in fp/lmaux.F subs/clsprm.F subs/suham.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#opnu
for name in fp/lmaux.F subs/clsprm.F subs/rdctrl2.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# opmpol
for name in subs/clsprm.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# oorhat
for name in fp/bndfp.F fp/lmfp.F fp/rsedit.F fp/supot.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# omad
for name in fp/supot.F 
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# osymgr
for name in fp/bndfp.F fp/chimedit.F fp/lmaux.F fp/lmfp.F fp/mshvmt.F fp/sugw.F fp/supot.F \
fp/symrho.F fp/totfrc.F subs/mkqp.F subs/mksym.F subs/rdsigm.F subs/rdsigm2.F subs/suham.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#oqlv
for name in fp/fklbl.F fp/gklbl.F fp/hklbl.F fp/hsmbl.F fp/lmaux.F fp/supot.F subs/lattic.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#okv
for name in fp/dfrce.F fp/hsibl.F fp/ioden.F fp/mixrho.F fp/rdovfa.F fp/rhgcmp.F fp/rsibl.F \
fp/smshft.F fp/smves.F fp/smvxcm.F fp/supot.F fp/symrho.F fp/vxcnlm.F subs/chgmsh.F \
subs/gvmtch.F subs/hft2rs.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#ojcg
for name in  fp/bndfp.F fp/fsmbl.F fp/ggugbl.F fp/hgugbl.F fp/hhugbl.F fp/locpot.F fp/makusq.F \
fp/mkrout.F fp/ovlocr.F fp/rlocbl.F fp/smhsbl.F fp/symrho.F fp/vcdmel.F subs/setcg.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit




# oistab
for name in fp/bndfp.F fp/chimedit.F fp/lmfp.F fp/sugw.F fp/symrho.F fp/totfrc.F subs/mksym.F \
subs/rdsigm.F subs/rdsigm2.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# oips0
for name in fp/mixrho.F fp/supot.F fp/symrho.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#oidxcg
for name in fp/augmbl.F fp/bndfp.F fp/fsmbl.F fp/ggugbl.F fp/hgugbl.F fp/hhugbl.F fp/locpot.F \
fp/makusq.F fp/mkrout.F fp/ovlocr.F fp/rlocbl.F fp/smhsbl.F fp/symrho.F fp/vcdmel.F subs/setcg.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#ogv
for name in fp/dfrce.F fp/hsibl.F fp/mixrho.F fp/rdovfa.F fp/rhgcmp.F fp/rsibl.F fp/smshft.F \
fp/smves.F fp/smvxcm.F fp/supot.F fp/symrho.F fp/vxcnlm.F slatsm/alloc.F slatsm/rdtok.F \
subs/chgmsh.F subs/gvlist.F subs/gvmtch.F subs/hft2rs.F subs/suham.F subs/suham2.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


#odlv
for name in fp/fklbl.F fp/gklbl.F fp/hklbl.F fp/hsmbl.F fp/lmaux.F fp/supot.F subs/lattic.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



#ocy
for name in fp/augmbl.F fp/fsmbl.F fp/ggugbl.F fp/hgugbl.F fp/hhugbl.F fp/makusq.F fp/mkrout.F \
fp/ovlocr.F fp/rlocbl.F fp/smcorm.F fp/smhsbl.F fp/smves.F fp/smvxcm.F fp/symrho.F fp/vcdmel.F subs/setcg.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# ocg
for name in fp/augmbl.F fp/fsmbl.F fp/ggugbl.F fp/hgugbl.F fp/hhugbl.F fp/locpot.F \
fp/makusq.F fp/mkrout.F fp/ovlocr.F fp/rlocbl.F fp/smhsbl.F fp/symrho.F fp/vcdmel.F subs/setcg.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# obgv
for name in fp/supot.F fp/symrho.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# oag
for name in fp/chimedit.F fp/lmaux.F fp/lmfp.F fp/mshvmt.F fp/supot.F fp/symrho.F subs/mksym.F subs/rdsigm.F subs/rdsigm2.F subs/suham.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit



# oindxo
for name in fp/addrbl.F fp/bndfp.F fp/hambl.F fp/hambls.F fp/makusq.F fp/mkekin.F fp/mkrout.F \
fp/sugw.F subs/rdsigm.F subs/rdsigm2.F subs/suham.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# oidtet
for name in fp/bndfp.F subs/asados.F subs/mkqp.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# owtkp
for name in fp/bndfp.F subs/asados.F subs/mkqp.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# ostar
for name in subs/asados.F subs/mkqp.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# oqp
for name in fp/bndfp.F subs/asados.F subs/clsprm.F subs/mkqp.F subs/rdsigm.F subs/rdsigm2.F lmv7util.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# oipq
for name in subs/asados.F subs/mkqp.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# ormax
for name in fp/lmaux.F subs/asados.F subs/clsprm.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit


# onrc onrcp
for name in fp/lmaux.F subs/mksym.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit

## oips
for name in fp/lmaux.F fp/supot.F fp/symrho.F subs/mcasim.F subs/mksym.F subs/rdctrl2.F subs/suham.F subs/asados.F
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit

# oclabl
for name in  fp/lmaux.F subs/mksym.F subs/asados.F 
do python delw.py < $dir/$name > x; mv -f x $dir/$name; done
exit

