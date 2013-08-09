for dir in main gwsrc nfpsrc tote slatsmlib
do 
for n in ${dir}/*.F; 
do echo $n
gawk -f stop2rx.awk $n > x
mv x $n
done

done
