#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_rdsigm2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  if (match($0,"slat%rv_p_opos *=> *rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_og/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_og *=> *slat%rv_p_osymgr")) { next;}
  if (match($0,"slat%rv_p_osymgr *=> *rv_p_og")) { next;}
  gsub("rv_p_og","slat%rv_p_osymgr")
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
/^ .*iv_p_oiprmb/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oiprmb *=> *sham%iv_p_oindxo")) { next;}
  if (match($0,"sham%iv_p_oindxo *=> *iv_p_oiprmb")) { next;}
  gsub("iv_p_oiprmb","sham%iv_p_oindxo")
  if (match($0,"sham%sham%")) {next;}
}
/^ .*iv_p_ooffh/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_ooffh *=> *sham%iv_p_ooffh")) { next;}
  if (match($0,"sham%iv_p_ooffh *=> *iv_p_ooffh")) { next;}
  gsub("iv_p_ooffh","sham%iv_p_ooffh")
  if (match($0,"sham%sham%")) {next;}
}
{print}
