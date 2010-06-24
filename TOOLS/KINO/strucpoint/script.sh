dir=../ecalj_2010_0510.2nd/lm7K
for name in `cat list.subs`
do echo $name
python delw.py   < $dir/$name >x
mv -f x $dir/$name
done
