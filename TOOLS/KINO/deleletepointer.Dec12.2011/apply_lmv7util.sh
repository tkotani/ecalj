
gawk -f tmp/replace1change.awk tmp/replace_lmv7util   > tmp/change_lmv7util.awk
gawk -f tmp/replace1change.awk tmp/replace_lmv7util_2 > tmp/change_lmv7util_2.awk
gawk -f tmp/replace1change.awk tmp/replace_lmv7util_4 > tmp/change_lmv7util_4.awk
gawk -f tmp/replace1change.awk tmp/replace_lmv7util_5 > tmp/change_lmv7util_5.awk
gawk -f tmp/replace1change.awk tmp/replace_lmv7util_6 > tmp/change_lmv7util_6.awk
gawk -f tmp/replace1change.awk tmp/replace_lmv7util_7 > tmp/change_lmv7util_7.awk

cat lmv7util.F | gawk -f tmp/change_lmv7util.awk | gawk -f tmp/change_lmv7util_2.awk | gawk -f tmp/change_lmv7util_3.awk| gawk -f tmp/change_lmv7util_4.awk | gawk -f tmp/change_lmv7util_5.awk | gawk -f tmp/change_lmv7util_6.awk | gawk -f tmp/change_lmv7util_7.awk | gawk -f tmp/change_lmv7util_8.awk> x
echo x is made.
