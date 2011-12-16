#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_asados_3
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_oqp/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqp *=> *sbz%rv_p_oqp")) { next;}
  if (match($0,"sbz%rv_p_oqp *=> *rv_p_oqp")) { next;}
  gsub("rv_p_oqp","sbz%rv_p_oqp")
  if (match($0,"sbz%sbz%")) {next;}
}
{print}
