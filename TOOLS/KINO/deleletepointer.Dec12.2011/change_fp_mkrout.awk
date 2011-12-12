BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_orhoca/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_orhoca *=> *sspec\\(is\\)%rv_p_orhoc")) { next;}
  if (match($0,"sspec\\(is\\)%rv_p_orhoc *=> *rv_p_orhoca")) { next;}
  gsub("rv_p_orhoca","sspec(is)%rv_p_orhoc")
  if (match($0,"sspec\\(is\\)%sspec\\(is\\)%")) {next;}
}
/^ .*rv_p_ov01/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ov01 *=> *ssite\\(ib\\)%rv_p_ov1")) { next;}
  if (match($0,"ssite\\(ib\\)%rv_p_ov1 *=> *rv_p_ov01")) { next;}
  gsub("rv_p_ov01","ssite(ib)%rv_p_ov1")
  if (match($0,"ssite\\(ib\\)%ssite\\(ib\\)%")) {next;}
}
/^ .*iv_p_ojcg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_ojcg *=> *slat%iv_p_ojcg")) { next;}
  if (match($0,"slat%iv_p_ojcg *=> *iv_p_ojcg")) { next;}
  gsub("iv_p_ojcg","slat%iv_p_ojcg")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oidxcg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oidxcg *=> *slat%iv_p_oidxcg")) { next;}
  if (match($0,"slat%iv_p_oidxcg *=> *iv_p_oidxcg")) { next;}
  gsub("iv_p_oidxcg","slat%iv_p_oidxcg")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ocy/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ocy *=> *slat%rv_p_ocy")) { next;}
  if (match($0,"slat%rv_p_ocy *=> *rv_p_ocy")) { next;}
  gsub("rv_p_ocy","slat%rv_p_ocy")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ocg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ocg *=> *slat%rv_p_ocg")) { next;}
  if (match($0,"slat%rv_p_ocg *=> *rv_p_ocg")) { next;}
  gsub("rv_p_ocg","slat%rv_p_ocg")
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
{print}
