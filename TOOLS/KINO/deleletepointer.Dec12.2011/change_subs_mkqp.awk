# please set LANG=C before invoking this script to show date in English
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_osymgr/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_osymgr *=> *slat%rv_p_osymgr")) { next;}
  if (match($0,"slat%rv_p_osymgr *=> *rv_p_osymgr")) { next;}
  gsub("rv_p_osymgr[^a-zA-Z0-9]","slat%rv_p_osymgr")
  if (match($0,"slat%slat%")) {next;}
}
{print}
