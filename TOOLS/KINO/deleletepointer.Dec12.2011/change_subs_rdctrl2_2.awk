#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_rdctrl2_2.awk
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ *real\(8\),pointer :: rv_p_opos\(:\) =>NULL\(\)/{
 print comment,$0
 next
}
/^ *allocate\(rv_p_opos\(abs\(3\*nsite\)\)\)/{
 print comment,$0
 sub("rv_p_opos","v_slat%rv_p_opos")
}
/^ *if \(3\*nsite<0\) rv_p_opos\(:\)=0.0d0/{
 print comment,$0
 sub("rv_p_opos","v_slat%rv_p_opos")
}
/^ *\.   , i_copy_size , i_spackv \+ 1 \- 1 , rv_p_opos \)/ {
 print comment,$0
 sub("rv_p_opos","v_slat%rv_p_opos")
}
/^ *v_slat%rv_p_opos => rv_p_opos/{
 print comment,$0
 next
}
/^ *v_sarry%rv_p_opos => rv_p_opos/{
 print comment,$0
 print "        v_sarry%rv_p_opos => v_slat%rv_p_opos"
 next
}
{print}
