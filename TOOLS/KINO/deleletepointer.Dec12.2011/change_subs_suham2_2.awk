#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_suham2_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ogv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ogv *=> *slat%rv_p_ogv")) { next;}
  if (match($0,"slat%rv_p_ogv *=> *rv_p_ogv")) { next;}
  gsub("rv_p_ogv","slat%rv_p_ogv")
  if (match($0,"slat%slat%")) {next;}
}
{print}
