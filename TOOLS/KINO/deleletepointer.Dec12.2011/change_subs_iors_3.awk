#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_iors_3
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*zv_p_osmrho/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"zv_p_osmrho *=> *spot%zv_p_osmrho")) { next;}
  if (match($0,"spot%zv_p_osmrho *=> *zv_p_osmrho")) { next;}
  gsub("zv_p_osmrho","spot%zv_p_osmrho")
  if (match($0,"spot%spot%")) {next;}
}
{print}
