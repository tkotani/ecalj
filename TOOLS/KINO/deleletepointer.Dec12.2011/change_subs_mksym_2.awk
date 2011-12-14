#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_mksym_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_osymgr/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_osymgr *=> *slat%rv_p_osymgr")) { next;}
  if (match($0,"slat%rv_p_osymgr *=> *rv_p_osymgr")) { next;}
  gsub("rv_p_osymgr","slat%rv_p_osymgr")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  if (match($0,"slat%rv_p_oag *=> *rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oistab/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oistab *=> *slat%iv_p_oistab")) { next;}
  if (match($0,"slat%iv_p_oistab *=> *iv_p_oistab")) { next;}
  gsub("iv_p_oistab","slat%iv_p_oistab")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oclabl/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oclabl *=> *sarray%rv_p_oclabl")) { next;}
  if (match($0,"sarray%rv_p_oclabl *=> *rv_p_oclabl")) { next;}
  gsub("rv_p_oclabl","sarray%rv_p_oclabl")
  if (match($0,"sarray%sarray%")) {next;}
}
{print}
