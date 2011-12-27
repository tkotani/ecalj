#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_fp_rsedit
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*sv_p_oorhat/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"sv_p_oorhat *=> *spot%sv_p_oorhat")) { next;}
  if (match($0,"spot%sv_p_oorhat *=> *sv_p_oorhat")) { next;}
  gsub("sv_p_oorhat","spot%sv_p_oorhat")
  if (match($0,"spot%spot%")) {next;}
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
