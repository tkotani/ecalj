#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_fp_lmfp_4
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  if (match($0,"slat%rv_p_opos *=> *rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
{print}
