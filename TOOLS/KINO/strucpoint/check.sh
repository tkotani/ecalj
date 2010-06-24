dir=/home/kino/kit/GW/7K/ecalj_2010_0618/lm7K
for name in  \
fp/elocp.F fp/locpot.F fp/mkrout.F fp/pnunew.F fp/rdovfa.F fp/rsedit.F fp/sugw.F \
fp/vcdmel.F subs/iors.F
do echo $name; python chp2.py "real(8)" rv_p_ ov0 < $dir/$name >x; mv -f x $dir/$name
done

for name in \
fp/locpot.F fp/mkrout.F fp/rdovfa.F fp/rsedit.F subs/iors.F
do echo $name; python chp2.py "real(8)" rv_p_ ov1 < $dir/$name >x; mv -f x $dir/$name
done


