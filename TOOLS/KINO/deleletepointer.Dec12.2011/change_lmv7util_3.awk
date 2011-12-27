#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*integer,pointer :: iv_p_oinitc/{
 print comment,$0
 next
}
/^ .*allocate\(iv_p_oinitc/{
 print comment,$0
 gsub("iv_p_oinitc","sarray%iv_p_ohave")
}
/^ .*if \(-nclasp<0\) iv_p_oinitc/{
 print comment,$0
 gsub("iv_p_oinitc","sarray%iv_p_ohave")
}
/^ .*sarray%iv_p_ohave => iv_p_oinitc/{
print comment,$0
next
}
{print}
