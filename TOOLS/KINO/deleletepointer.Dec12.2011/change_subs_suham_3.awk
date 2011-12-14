#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_suham_3
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*iv_p_oidxsh/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oidxsh *=> *sham%iv_p_oindxo")) { next;}
  if (match($0,"sham%iv_p_oindxo *=> *iv_p_oidxsh")) { next;}
  gsub("iv_p_oidxsh","sham%iv_p_oindxo")
  if (match($0,"sham%sham%")) {next;}
}
/^ .*iv_p_offh/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_offh *=> *sham%iv_p_ooffh")) { next;}
  if (match($0,"sham%iv_p_ooffh *=> *iv_p_offh")) { next;}
  gsub("iv_p_offh","sham%iv_p_ooffh")
  if (match($0,"sham%sham%")) {next;}
}
{print}
