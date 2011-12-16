#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_iors_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_orhoca/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_orhoca *=> *sspec\\(is\\)%rv_p_orhoc")) { next;}
  if (match($0,"sspec\\(is\\)%rv_p_orhoc *=> *rv_p_orhoca")) { next;}
  gsub("rv_p_orhoca","sspec(is)%rv_p_orhoc")
  if (match($0,"sspec\\(is\\)%sspec\\(is\\)%")) {next;}
}
{print}
