#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util_7
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_opmpol/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_opmpol *=> *spot%rv_p_opmpol")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_opmpol *=> *rv_p_opmpol")) { print comment,$0;next;}
  gsub("rv_p_opmpol","spot%rv_p_opmpol")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
{print}
