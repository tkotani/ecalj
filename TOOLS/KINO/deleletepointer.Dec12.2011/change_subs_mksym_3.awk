#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_mksym_3
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*iv_p_oics/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oics *=> *sarray%iv_p_oics")) { next;}
  if (match($0,"sarray%iv_p_oics *=> *iv_p_oics")) { next;}
  gsub("iv_p_oics","sarray%iv_p_oics")
  if (match($0,"sarray%sarray%")) {next;}
}
{print}
