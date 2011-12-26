#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_rdsigm2_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ohrs/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ohrs *=> *sham%rv_p_ohrs")) { next;}
  if (match($0,"sham%rv_p_ohrs *=> *rv_p_ohrs")) { next;}
  gsub("rv_p_ohrs","sham%rv_p_ohrs")
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
