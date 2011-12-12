# please set LANG=C before invoking this script to show date in English
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*iv_p_oipc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oipc *=> *sarray%iv_p_oipc")) { next;}
  if (match($0,"sarray%iv_p_oipc *=> *iv_p_oipc")) { next;}
  gsub("iv_p_oipc[^a-zA-Z0-9]","sarray%iv_p_oipc")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_ogstar/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_ogstar *=> *sbz%iv_p_ostar")) { next;}
  if (match($0,"sbz%iv_p_ostar *=> *iv_p_ogstar")) { next;}
  gsub("iv_p_ogstar[^a-zA-Z0-9]","sbz%iv_p_ostar")
  if (match($0,"sbz%sbz%")) {next;}
}
/^ .*iv_p_oidtet/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oidtet *=> *sbz%iv_p_oidtet")) { next;}
  if (match($0,"sbz%iv_p_oidtet *=> *iv_p_oidtet")) { next;}
  gsub("iv_p_oidtet[^a-zA-Z0-9]","sbz%iv_p_oidtet")
  if (match($0,"sbz%sbz%")) {next;}
}
/^ .*rv_p_ormax/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ormax *=> *sarray%rv_p_ormax")) { next;}
  if (match($0,"sarray%rv_p_ormax *=> *rv_p_ormax")) { next;}
  gsub("rv_p_ormax[^a-zA-Z0-9]","sarray%rv_p_ormax")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_oipq/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oipq *=> *sbz%iv_p_oipq")) { next;}
  if (match($0,"sbz%iv_p_oipq *=> *iv_p_oipq")) { next;}
  gsub("iv_p_oipq[^a-zA-Z0-9]","sbz%iv_p_oipq")
  if (match($0,"sbz%sbz%")) {next;}
}
/^ .*rv_p_owtkp/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_owtkp *=> *sbz%rv_p_owtkp")) { next;}
  if (match($0,"sbz%rv_p_owtkp *=> *rv_p_owtkp")) { next;}
  gsub("rv_p_owtkp[^a-zA-Z0-9]","sbz%rv_p_owtkp")
  if (match($0,"sbz%sbz%")) {next;}
}
/^ .*iv_p_oics/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oics *=> *sarray%iv_p_oics")) { next;}
  if (match($0,"sarray%iv_p_oics *=> *iv_p_oics")) { next;}
  gsub("iv_p_oics[^a-zA-Z0-9]","sarray%iv_p_oics")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_oclabl/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oclabl *=> *sarray%rv_p_oclabl")) { next;}
  if (match($0,"sarray%rv_p_oclabl *=> *rv_p_oclabl")) { next;}
  gsub("rv_p_oclabl[^a-zA-Z0-9]","sarray%rv_p_oclabl")
  if (match($0,"sarray%sarray%")) {next;}
}
{print}
