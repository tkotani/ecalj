BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
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
{print}
