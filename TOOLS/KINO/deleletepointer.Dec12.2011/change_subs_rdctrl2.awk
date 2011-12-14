#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_rdctrl2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*iv_p_oips/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oips *=> *v_sarry%iv_p_oips")) { next;}
  if (match($0,"v_sarry%iv_p_oips *=> *iv_p_oips")) { next;}
  gsub("iv_p_oips","v_sarry%iv_p_oips")
  if (match($0,"v_sarry%v_sarry%")) {next;}
}
{print}
