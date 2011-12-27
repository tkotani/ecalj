#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util_4
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_opnu/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_opnu *=> *spot%rv_p_opnu")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_opnu *=> *rv_p_opnu")) { print comment,$0;next;}
  gsub("rv_p_opnu","spot%rv_p_opnu")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*iv_p_ohave/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"iv_p_ohave *=> *sarray%iv_p_ohave")) { print comment,$0;next;}
  if (match($0,"sarray%iv_p_ohave *=> *iv_p_ohave")) { print comment,$0;next;}
  gsub("iv_p_ohave","sarray%iv_p_ohave")
  gsub("sarray%sarray%","sarray%")
  if (old != $0 ) {print comment,$0}
}
{print}
