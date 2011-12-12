# please set LANG=C before invoking this script to show date in English
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*iv_p_offh/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_offh *=> *sham%iv_p_ooffh")) { next;}
  if (match($0,"sham%iv_p_ooffh *=> *iv_p_offh")) { next;}
  gsub("iv_p_offh[^a-zA-Z0-9]","sham%iv_p_ooffh")
  if (match($0,"sham%sham%")) {next;}
}
{print}
