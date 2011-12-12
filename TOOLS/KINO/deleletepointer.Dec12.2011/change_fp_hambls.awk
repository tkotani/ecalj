BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_ontabs/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_ontabs *=> *sham%iv_p_ontabs")) { next;}
  if (match($0,"sham%iv_p_ontabs *=> *iv_p_ontabs")) { next;}
  gsub("iv_p_ontabs","sham%iv_p_ontabs")
  if (match($0,"sham%sham%")) {next;}
}
/^ .*iv_p_oiaxs/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oiaxs *=> *sham%iv_p_oiaxs")) { next;}
  if (match($0,"sham%iv_p_oiaxs *=> *iv_p_oiaxs")) { next;}
  gsub("iv_p_oiaxs","sham%iv_p_oiaxs")
  if (match($0,"sham%sham%")) {next;}
}
/^ .*rv_p_ohrs/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ohrs *=> *sham%rv_p_ohrs")) { next;}
  if (match($0,"sham%rv_p_ohrs *=> *rv_p_ohrs")) { next;}
  gsub("rv_p_ohrs","sham%rv_p_ohrs")
  if (match($0,"sham%sham%")) {next;}
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
