BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_ov0/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ov0 *=> *ssite\\(ib\\)%rv_p_ov0")) { next;}
  if (match($0,"ssite\\(ib\\)%rv_p_ov0 *=> *rv_p_ov0")) { next;}
  gsub("rv_p_ov0","ssite(ib)%rv_p_ov0")
  if (match($0,"ssite\\(ib\\)%ssite\\(ib\\)%")) {next;}
}
/^ .*rv_p_og/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_og *=> *slat%rv_p_osymgr")) { next;}
  if (match($0,"slat%rv_p_osymgr *=> *rv_p_og")) { next;}
  gsub("rv_p_og","slat%rv_p_osymgr")
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
/^ .*rv_p_oqsig/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqsig *=> *sham%rv_p_oqsig")) { next;}
  if (match($0,"sham%rv_p_oqsig *=> *rv_p_oqsig")) { next;}
  gsub("rv_p_oqsig","sham%rv_p_oqsig")
  if (match($0,"sham%sham%")) {next;}
}
{print}
