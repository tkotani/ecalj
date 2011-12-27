#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util_6
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ovintr/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ovintr *=> *spot%rv_p_ovintr")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_ovintr *=> *rv_p_ovintr")) { print comment,$0;next;}
  gsub("rv_p_ovintr","spot%rv_p_ovintr")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
{print}
