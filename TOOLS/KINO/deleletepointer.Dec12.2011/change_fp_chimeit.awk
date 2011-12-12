BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_og/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_og *=> *slat%rv_p_osymgr")) { next;}
  gsub("rv_p_og","slat%rv_p_osymgr")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oistab/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oistab *=> *slat%iv_p_oistab")) { next;}
  gsub("iv_p_oistab","slat%iv_p_oistab")
  if (match($0,"slat%slat%")) {next;}
}
{print}
